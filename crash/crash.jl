# A minimal, self-contained reproduction of the dust-temperature abort
#   "Newton-Raphson iteration for dust temperature failed to converge or dust temperature is negative!"
# (QuokkaSimulation.hpp, raised when p_iteration_failure_counter[1] > 0) that DTypeFront3D hits on
# Coarse STEP 70 (t = 9.530016024e+11 s, 3.02% of stop_time) of the 256^3, 8-GPU run of inputs/DTypeFront3D.toml.
#
# Every number below is hardcoded, copied verbatim (%.17g, round-trips exactly) from the DUSTFAIL debug dump
# that the instrumented SolveGasDustRadiationEnergyExchange (radiation_dust_system.hpp) prints on failure.
# The full dump -- 64 failing calls, 8 per MPI rank -- is in crash/dustfail_dump.txt; lines are prefixed with
# the srun rank label ("r: DUSTFAIL id kind ..."), and a record is identified by (rank, id). This file
# replays record (rank 0, id 0). All 64 records fail the same way: dust_model == 2 (coeff_n = 0), and the
# first Newton step drives T_d negative (where = 1, n = 1), with T_gas0 ~ 1e15 K.
#
# Q_dust is not an argument of the C++ solver in this quokka tree; passing 0.0 reduces every Julia
# expression that uses it to the C++ one.
#
# Run with:  julia crash.jl

include(joinpath(@__DIR__, "..", "types.jl")) # defines OpacityModel

# --- minimal JSON writer, so the result can be read from a notebook without a Julia kernel. Every value
# --- here is a Float64, Int or a flat vector of one of those, so this needs no escaping/nesting logic
# --- beyond what is written below -- not a general-purpose JSON encoder. -------------------------------
json_num(x::Float64) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Infinity" : "-Infinity") : repr(x)
json_num(x::Int) = repr(x)
json_arr(v) = "[" * join(json_num.(v), ", ") * "]"

function json_iteration_state(s::SolverIterationState)
	return "{\"n\": $(json_num(s.n)), \"T_gas\": $(json_num(s.T_gas)), \"T_d\": $(json_num(s.T_d)), " *
	       "\"Egas_guess\": $(json_num(s.Egas_guess)), \"EradVec_guess\": $(json_arr(s.EradVec_guess)), " *
	       "\"Rvec\": $(json_arr(s.Rvec)), \"tau\": $(json_arr(s.tau)), \"delta_x\": $(json_num(s.delta_x)), " *
	       "\"delta_R\": $(json_arr(s.delta_R)), \"F0\": $(json_num(s.F0)), \"Fg\": $(json_arr(s.Fg)), " *
	       "\"Fg_abs_sum\": $(json_num(s.Fg_abs_sum)), \"Fg_roundoff\": $(json_num(s.Fg_roundoff)), " *
	       "\"Etot0\": $(json_num(s.Etot0)), \"F0_resid_ratio\": $(json_num(s.F0_resid_ratio)), " *
	       "\"Fg_resid_ratio\": $(json_num(s.Fg_resid_ratio)), " *
	       "\"Fg_roundoff_ratio\": $(json_num(s.Fg_roundoff_ratio)), \"relax\": $(json_num(s.relax))}"
end

function json_opacity_terms(o::OpacityTerms)
	return "{\"kappaE\": $(json_arr(o.kappaE)), \"kappaP\": $(json_arr(o.kappaP)), " *
	       "\"kappaF\": $(json_arr(o.kappaF)), \"kappaPoverE\": $(json_arr(o.kappaPoverE)), " *
	       "\"delta_nu_kappa_B_at_edge\": $(json_arr(o.delta_nu_kappa_B_at_edge)), " *
	       "\"alpha_P\": $(json_arr(o.alpha_P)), \"alpha_E\": $(json_arr(o.alpha_E))}"
end

function json_result(r::NewtonIterationResult)
	iterations_json = "[" * join(json_iteration_state.(r.iterations), ", ") * "]"
	return "{\"Egas\": $(json_num(r.Egas)), \"T_gas\": $(json_num(r.T_gas)), \"T_d\": $(json_num(r.T_d)), " *
	       "\"EradVec\": $(json_arr(r.EradVec)), \"work\": $(json_arr(r.work)), " *
	       "\"opacity_terms\": $(json_opacity_terms(r.opacity_terms)), \"iterations\": $iterations_json}"
end

# =====================================================================================================
# Traits -- the DUSTFAIL "traits" line (identical in all 64 records)
# =====================================================================================================

const nGroups_ = 3
const nGroupsThermal_ = 2
const beta_order_ = 1
const opacity_model_ = OpacityModel(1) # piecewise_constant_opacity
const gamma_ = 1.6666666666666667
const c_light_ = 29979245800.0
const c_hat_ = 29979245.800000001	# c / 1000
const radiation_constant_ = 7.5657332502800087e-15	# a_rad, copied verbatim -- NOT 4 sigma_SB / c
const boltzmann_constant_ = 1.3806490000000002e-16
const energy_unit_ = 1.6021766339999999e-12	# eV in erg -- group boundaries below are in eV
const Erad_floor_ = 7.2632007407999998e-22	# already divided by nGroups_, as radiation_system.hpp does
const mean_molecular_mass_ = 1.0
const gas_dust_coupling_threshold = 9.9999999999999995e-07

# DUSTFAIL "edge" lines: kappa_lower per group (all exponents 0) = [kappa_ir, kappa_optical, 0]
const kappa_ir_ = 0.01
const kappa_optical_ = 1000.0

# EOS species table (spmasses / gammas), filled by actual_eos_init from extern/Microphysics/EOS/photoionization/_parameters,
# in the network's species order: H, H+, e- (build/3d/src/problems/DTypeFront3D/network_properties.H). Not in the dump;
# the T_gas0 / c_v0 cross-check below verifies them (they reproduce the dumped values bit-for-bit).
const spmasses_ = [1.673773e-24, 1.6728620616289998e-24, 9.10938371e-28]
const gammas_ = [1.6666666666666667, 1.6666666666666667, 1.6666666666666667]

include(joinpath(@__DIR__, "..", "hyperparameters.jl"))
include(joinpath(@__DIR__, "..", "support.jl"))
include(joinpath(@__DIR__, "..", "planck_integral.jl"))
include(joinpath(@__DIR__, "..", "radiation_defaults.jl"))

# --- DTypeFront3D's own definitions, same shapes as DTypeFront3D/problem_DTypeFront3D.jl -------------

# RadSystem<DTypeFront3D>::DefineOpacityExponentsAndLowerValues. Each thermal group carries its own
# constant gray opacity; the ionizing (chemistry) band is transparent to it.
function DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, Tgas)
	kappa_g = vcat([kappa_ir_, kappa_optical_], zeros(nGroups_ - 2))
	exponents_and_values = [zeros(nGroups_ + 1), zeros(nGroups_ + 1)]
	for i in 1:(nGroups_ + 1)
		exponents_and_values[1][i] = 0.0
		exponents_and_values[2][i] = (i <= nGroups_) ? kappa_g[i] : 0.0
	end
	return exponents_and_values
end

# Port of quokka::EOSMicrophysics -- see DTypeFront3D/problem_DTypeFront3D.jl for the full derivation.
function eos_sums_(massScalars)
	rhotot = 0.0
	sum_ni_fi_over_2 = 0.0
	for n in 1:length(massScalars)
		xn = massScalars[n] / spmasses_[n]
		rhotot += xn * spmasses_[n]
		sum_ni_fi_over_2 += xn / (gammas_[n] - 1.0)
	end
	return rhotot, sum_ni_fi_over_2
end

function ComputeTgasFromEint(rho, Eint, massScalars)
	rhotot, s = eos_sums_(massScalars)
	return (Eint / rho) * rhotot / (s * boltzmann_constant_)
end

function ComputeEintFromTgas(rho, Tgas, massScalars)
	rhotot, s = eos_sums_(massScalars)
	return (s * boltzmann_constant_ * Tgas / rhotot) * rho
end

function ComputeEintTempDerivative(rho, Tgas, massScalars)
	rhotot, s = eos_sums_(massScalars)
	return (s * boltzmann_constant_ / rhotot) * rho
end

include(joinpath(@__DIR__, "..", "opacity.jl"))
include(joinpath(@__DIR__, "..", "jacobian.jl"))
include(joinpath(@__DIR__, "..", "solver.jl"))

# Distance in units in the last place between two Float64s (0 = bit-identical).
function ulp_distance(a::Float64, b::Float64)
	a == b && return 0
	(isnan(a) || isnan(b)) && return typemax(Int)
	ia = reinterpret(Int64, a)
	ib = reinterpret(Int64, b)
	ia = ia < 0 ? typemin(Int64) - ia : ia
	ib = ib < 0 ? typemin(Int64) - ib : ib
	return abs(ia - ib)
end

function check(name, julia_value, cpp_value)
	d = ulp_distance(julia_value, cpp_value)
	println(rpad(name, 24), " julia = ", repr(julia_value), "  C++ = ", repr(cpp_value), "  ULP = ", d)
	return d
end

# =====================================================================================================
# The crashing call: DUSTFAIL record (rank 0, id 0), Coarse STEP 70
# =====================================================================================================

function main()
	# DUSTFAIL "scalars" line
	Egas0 = 2.7931396745623172
	rho = 5.7256058727513376e-24
	coeff_n = 0.0
	dt = 10.446078398865433
	n_outer_iter = 0
	tol = 1.0e-10
	tol_rel = -1.0
	tempFloor = 10.0
	Q_dust = 0.0

	# DUSTFAIL "group" lines, g = 0, 1, 2
	Erad0Vec = [3.0000434388425533e-10, 3.2864713504122893e-09, 2.8873297972796427e-08]
	work = [0.0, 0.0, 0.0]
	vel_times_F = [6.6417262252495831e-14, 6.8923430825582663e-11, 6.0995425657430077e-10]
	Src = [0.0, 0.0, 0.0]

	# DUSTFAIL "edge" lines (eV) and "massScalar" lines (H, H+, e-)
	rad_boundaries = [9.9999999999999995e-07, 0.41356700000000002, 13.6, 26.0]
	massScalars = [1.0094063843583376e-28, 5.7223888691007661e-24, 3.1160630121356934e-27]

	# What the C++ saw when it flagged the failure (DUSTFAIL "state" line)
	cpp_T_gas0 = 1971365050104450.2
	cpp_T_d0 = 1971365050104450.2
	cpp_c_v0 = 1.4168556323012454e-15
	cpp_T_d_at_failure = -6382108661.5 # Newton iterate n = 1, dust_model = 2
	cpp_fourPiBoverC_Tgas0 = [84626729463.325348, 3009362029319383.0, 0.0]

	println("--- cross-checks against the C++ dump (EOS, Planck integrals, initial dust temperature) ---")
	T_gas0 = ComputeTgasFromEint(rho, Egas0, massScalars)
	check("T_gas0", T_gas0, cpp_T_gas0)
	check("c_v0", ComputeEintTempDerivative(rho, T_gas0, massScalars), cpp_c_v0)
	fourPiBoverC = ComputeThermalRadiationMultiGroup(T_gas0, rad_boundaries)
	for g in 1:nGroups_
		check("fourPiBoverC_Tgas0[$g]", fourPiBoverC[g], cpp_fourPiBoverC_Tgas0[g])
	end
	check("T_d0 (Bate-Keto)", ComputeDustTemperatureBateKeto(T_gas0, T_gas0, rho, Erad0Vec, coeff_n, dt, NaN, 0, Q_dust, rad_boundaries),
	      cpp_T_d0)
	println()

	p_iteration_counter = zeros(Int, 4)
	p_iteration_failure_counter = zeros(Int, 3)

	result = SolveGasDustRadiationEnergyExchange(
		Egas0, Erad0Vec, rho, coeff_n, dt, massScalars, n_outer_iter, work, vel_times_F, Src, Q_dust,
		rad_boundaries, tol, tol_rel, tempFloor, p_iteration_counter, p_iteration_failure_counter, true)

	n = p_iteration_counter[2] - 1
	println("Newton iterations n = $n  (maxIter = 100)")
	println("failure counters [non-convergence, negative T_d, outer] = $p_iteration_failure_counter")
	println("decoupled dust branch (dust_model == 2) taken: $(p_iteration_counter[4] == 1)")
	println("T_gas = $(result.T_gas) K,  T_d = $(result.T_d) K")
	println("Egas  = $(result.Egas)")
	println("Erad  = $(result.EradVec)")
	if length(result.iterations) >= 2
		check("T_d at n = 1", result.iterations[2].T_d, cpp_T_d_at_failure)
	end
	println()
	if p_iteration_failure_counter[2] > 0
		println("REPRODUCED: the dust temperature went negative (the C++ dust-temperature abort)")
	else
		println("NOT reproduced: the dust temperature stayed non-negative")
	end

	dump_path = joinpath(@__DIR__, "result.json")
	write(dump_path, json_result(result))
	println("\nwrote result (with per-iteration history) to $dump_path")
	println("reload from Python with: json.load(open(\"$dump_path\"))")

	return result
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
