# A minimal, self-contained reproduction of the Newton-Raphson non-convergence abort
# (radiation_dust_system.hpp:671) that DTypeFront3D hits on Coarse STEP 1 of DTypeFront3D_crash.toml.
#
# Every number below is hardcoded, not read from a dump -- unlike tests/repro_DTypeFront3D.jl. They are
# copied verbatim (full Float64 precision, via Julia's round-tripping repr) from the crashing call's
# SGIN 15855 record and the SGTRAITS lines in tests/dump.txt, which tests/repro_DTypeFront3D.jl already
# proved reproduces the abort. This file exists to call the solver with no dump-parsing machinery in the
# way, for quick edits/experiments against this one failing case.
#
# Run with:  julia crash.jl

include(joinpath(@__DIR__, "..", "types.jl")) # defines OpacityModel

# --- minimal JSON writer, so the result can be read from temp.ipynb without a Julia kernel. Every value
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
# Traits (SGTRAITS line, DTypeFront3D run) -- Physics_Traits / RadSystem_Traits / ISM_Traits / EOS_Traits
# =====================================================================================================

const nGroups_ = 3
const nGroupsThermal_ = 2
const beta_order_ = 1
const opacity_model_ = OpacityModel(1) # piecewise_constant_opacity
const gamma_ = 1.6666666666666667
const c_light_ = 2.99792458e10
const c_hat_ = 2.99792458e7
const radiation_constant_ = 7.565733250280009e-15	# a_rad, copied verbatim -- NOT 4 sigma_SB / c (see solver.jl history)
const boltzmann_constant_ = 1.3806490000000002e-16
const energy_unit_ = 6.62607015e-27	# h -- group boundaries below are in Hz
const Erad_floor_ = 1.246805533225e-21	# already divided by nGroups_, as radiation_system.hpp does
const mean_molecular_mass_ = 1.0
const gas_dust_coupling_threshold = 1.0e-6
const kappa_ir_ = 0.01
const kappa_optical_ = 1000.0

# EOS species table (spmasses / eos_gammas), filled by actual_eos_init from the network's parameters
const spmasses_ = [9.10938291e-28, 1.673532715291e-24, 1.672621777e-24]
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

# =====================================================================================================
# The crashing call: SGIN 15855, cell (15,15,15), outer iter 0
# =====================================================================================================

function main()
	Egas0 = 8.689688612484974e-13
	rho = 1.6735327152926736e-22
	coeff_n = 0.0
	dt = 3.0878236770052437e10
	Q_dust = 2.6813418069295e-20
	tol = 1.0e-10
	tol_rel = -1.0
	tempFloor = 10.0
	n_outer_iter = 0

	Erad0Vec = [6.260021242303187e-17, 3.740416599676655e-21, 3.740416599675e-41]
	work = [0.0, 0.0, 0.0]
	vel_times_F = [0.0, 0.0, 0.0]
	Src = [0.0, 4.259979144172045e-10, 0.0]
	rad_boundaries = [1.0e8, 1.0e14, 3.29e15, 8.0e15]
	massScalars = [6.7639341145484424e-37, 1.673532715280247e-22, 1.2419615697312844e-33]

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
	println()
	if p_iteration_failure_counter[1] > 0
		println("REPRODUCED: the Newton-Raphson iteration hit maxIter without converging,")
	else
		println("NOT reproduced: the Julia solve converged in $n iterations")
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
