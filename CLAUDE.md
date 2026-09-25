# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

A standalone Julia port of quokka's gas-dust-radiation energy exchange solver
(`RadSystem<problem_t>::SolveGasDustRadiationEnergyExchange`, in
`src/radiation/radiation_dust_system.hpp` of the main quokka repo). It exists to debug one
specific bug outside the C++ build: a `Newton-Raphson iteration for matter-radiation coupling
failed to converge!` abort (`radiation_dust_system.hpp:671`) that the `DTypeFront3D` problem
hits on Coarse STEP 1 of `DTypeFront3D_crash.toml`.

This is not a general-purpose radiation solver reimplementation. It is a scratch/debug harness:
port the exact code path, replay the exact failing input, and compare bit-for-bit against the
C++ result to prove or disprove a fix.

## File map

- `main.jl` — generic driver with hand-picked inputs (`rho`, `T_gas`, `T_rad`, `dt`, ...). Include
  order at the top of the file doubles as the dependency graph.
- `hyperparameters.jl` — line-for-line mirror of `radiation_system.hpp` lines 37-70.
- `types.jl` — `OpacityModel` enum (integer values match the C++ enum exactly) and the mutable
  structs (`OpacityTerms`, `JacobianResult`, `NewtonIterationResult`, `SolverIterationState`).
- `support.jl` — C++/AMReX shims: `amrex_assert`, `std_sqrt/max/min/pow/exp/log10` (calls libm via
  `ccall` where Julia's own implementation can differ from libm in the last bit — see the comment
  there), and `BackwardEulerOneVariable` (port of the bracketed Newton/bisection root finder).
- `planck_integral.jl` — port of `planck_integral.hpp`, including its 1000-entry interpolation
  table extracted mechanically from the C++ header.
- `radiation_defaults.jl` — problem-independent quokka defaults (Planck emission and its
  temperature derivative, zero cooling rates). These ARE faithful ports.
- `problem.jl` — **stand-in** EOS and opacity, used only by `main.jl`. Read its header before
  trusting any number derived from it: it replaces per-problem code each real quokka problem
  defines for itself (`hydro/EOS.hpp`, `DefineOpacityExponentsAndLowerValues`), and is not a port
  of anything.
- `opacity.jl` — `ComputeModelDependentKappa*` (group-mean opacity assembly) and
  `ComputeDustTemperatureBateKeto`, ported from `source_terms_multi_group.hpp` /
  `radiation_system.hpp`.
- `jacobian.jl` — `ComputeJacobianForGasAndDust[Decoupled]`, `RebaseThinGroupsOntoErad!`,
  `SolveLinearEqs`, ported from `radiation_dust_system.hpp` / `radiation_system.hpp`. Currently has
  debug `println`s left in `SolveLinearEqs` and `RebaseThinGroupsOntoErad!`.
- `solver.jl` — `SolveGasDustRadiationEnergyExchange` itself, aligned 1:1 with
  `radiation_dust_system.hpp` lines 239-698 (a stale commented-out C++ debug block near the end of
  the Newton loop marks exactly where the two sources are kept in step).
- `DTypeFront3D/problem_DTypeFront3D.jl` — the real `DTypeFront3D` problem definitions (EOS,
  opacity), but read entirely from a C++ debug dump's `SGTRAITS` lines rather than transcribed by
  hand, so there is nothing to get wrong by copying. Requires `SGDUMP_PATH` to be set by the
  includer.
- `tests/dump_io.jl` — parses the `SGIN`/`SGOUT`/`SGTRAITS` debug dump format, computes ULP
  distance between Julia and C++ doubles, and drives `SolveGasDustRadiationEnergyExchange` with a
  dump record's exact inputs.
- `tests/repro_DTypeFront3D.jl` — replays the one aborted (crashing) call from a dump and checks
  whether the Julia port reproduces the non-convergence failure. Usage:
  `julia repro_DTypeFront3D.jl [dump.txt]`.
- `tests/crosscheck_DTypeFront3D.jl` — replays every *completed* call in a dump and diffs the
  Julia port against the C++ result bit-for-bit (exact IEEE-754 equality, not "close enough").
  Usage: `julia crosscheck_DTypeFront3D.jl dump.txt [julia_out.txt]`.
- `tests/dump.txt`, `tests/julia_out.txt` — a captured dump and the crosscheck's output; regenerate
  per the instructions at the top of `crosscheck_DTypeFront3D.jl`, don't hand-edit.
- `crash/crash.jl` — self-contained reproduction of the same crash with the crashing call's
  inputs (`SGIN 15855`, cell `(15,15,15)`, outer iter 0) hardcoded verbatim, no dump-parsing
  machinery. For quick edits/experiments against this one failing case; also writes a JSON dump of
  the full per-iteration Newton history to `crash/result.json` for inspection outside Julia (e.g.
  from `crash/temp.ipynb`).

## How to generate/refresh a dump

1. Build quokka with the debug dump enabled in `src/radiation/source_terms_multi_group.hpp` (and
   the `SGTRAITS` block in the problem file).
2. Run the instrumented binary and grep its stderr/stdout for the dump lines:
   ```
   ./src/problems/DTypeFront3D/DTypeFront3D ./DTypeFront3D_crash.toml 2>&1 \
       | grep -E '^SGTRAITS|^SGIN |^SGOUT ' > dump.txt
   ```
3. Feed that file to `repro_DTypeFront3D.jl` or `crosscheck_DTypeFront3D.jl`.

## Rules for porting C++ into this repo

- **libm, not Julia's builtins, for `pow`/`exp`/`log10`.** Julia implements `^`, `exp`, `log10` in
  Julia rather than calling the system libm, and the two disagree in the last bit often enough to
  matter (measured: `x^3.0` differs from `pow(x, 3.0)` for about a quarter of arguments). Anywhere
  the C++ writes `std::pow`/`std::exp`/`std::log10`, use `std_pow`/`std_exp`/`std_log10` from
  `support.jl`. `sqrt` is exempt — IEEE-754 requires correct rounding, so Julia's `sqrt` and
  libm's always agree.
- **`std_max`/`std_min`, not Julia's `max`/`min`, when a NaN operand can occur.** `std::max(a, b)`
  is `a < b ? b : a`, so a NaN operand returns the *other* operand; Julia's `max`/`min` propagate
  NaN. The solver is written to survive a negative/NaN temperature by continuing past assertions
  and recording a failure counter — reproducing that requires the C++ semantics exactly.
- **`pi^4` is not `PI*PI*PI*PI`.** Julia's `pi` is an `Irrational`, so `pi^4` is evaluated exactly
  and rounded once, landing 1 ULP above `Float64(pi) * Float64(pi) * Float64(pi) * Float64(pi)`,
  which is what the C++ actually computes. See `gInf` in `planck_integral.jl` and the `PlanckFunction`
  comment in `radiation_defaults.jl` for the pattern to copy.
- **`amrex_debug_ = false` reproduces Release-build behaviour** (assertions compile away and the
  code continues, incrementing failure counters). Only flip it to `true` for a Debug-style run.
- Before trusting any comparison against quokka, re-check the file header of whatever you're
  reading — `problem.jl` in particular is explicitly *not* a reference implementation, whereas
  `radiation_defaults.jl`, `solver.jl`, `opacity.jl`, `jacobian.jl` are meant to be faithful,
  checkable ports (each carries the exact source file/line range it mirrors).
- When adding a new problem's stand-ins (a second `<Problem>/problem_<Problem>.jl`), prefer
  deriving constants from a dump's `SGTRAITS` lines over hand-transcribing them, following the
  `DTypeFront3D` example — it eliminates an entire class of transcription bugs.

## Running

Plain Julia, no package environment/`Project.toml` in this directory — just `julia <file>.jl`
from this directory (or pass a dump path as `ARGS[1]` where the driver supports it).

## Known repo state

- `crash/temp.ipynb` is gitignored (see `.gitignore`); it's a scratch notebook for inspecting
  `crash/result.json` outside Julia.
- `jacobian.jl` has debug `println`s in `SolveLinearEqs` and `RebaseThinGroupsOntoErad!` that are
  not gated behind the `debug` flag `solver.jl` otherwise uses — noisy but harmless; don't be
  surprised by the console spam when running anything that hits those functions.
