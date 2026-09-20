# DustSolver -- a standalone Julia driver for quokka's gas-dust-radiation energy exchange solver.
#
#   solver.jl          SolveGasDustRadiationEnergyExchange, aligned 1:1 with
#                      quokka src/radiation/radiation_dust_system.hpp lines 239-698
#   hyperparameters.jl mirror of radiation_system.hpp lines 37-70
#   types.jl           OpacityTerms, JacobianResult, NewtonIterationResult, OpacityModel
#   jacobian.jl        ComputeJacobianForGasAndDust[Decoupled], RebaseThinGroupsOntoErad, SolveLinearEqs
#   opacity.jl         ComputeModelDependentKappa*, ComputeDustTemperatureBateKeto
#   support.jl         amrex_assert, std_sqrt/max/min, BackwardEulerOneVariable
#   problem.jl         EOS, Planck, opacity and cooling -- read its header, half of it is stand-ins
#
# Run with:  julia main.jl

include("types.jl") # defines OpacityModel, needed by opacity_model_ below

# =====================================================================================================
# Problem constants. In quokka these come from Physics_Traits / RadSystem_Traits / ISM_Traits.
# =====================================================================================================

const nGroups_ = 1		# problem.jl only ports the nGroups_ == 1 Planck path
const nGroupsThermal_ = 1	# nGroups_ minus the number of chemical (ionizing) bands
const opacity_model_ = piecewise_constant_opacity
const gamma_ = 5.0 / 3.0
const beta_order_ = 1

const c_light_ = 2.99792458e10		# cgs
const c_hat_ = 0.1 * c_light_		# reduced speed of light
const sigma_SB_ = 5.670374419184432e-05
const radiation_constant_ = 4.0 * sigma_SB_ / c_light_ # a_rad, formed as quokka does (4 sigma_SB / c)
const boltzmann_constant_cgs_ = 1.380649e-16
# energy_unit_ and boltzmann_constant_ set the units of the group boundaries: PlanckFunction and
# ComputePlanckEnergyFractions both form x = nu * energy_unit_ / (boltzmann_constant_ * T). With
# energy_unit_ = 1 and k_B in cgs, boundaries are in ergs. Quokka problems usually set energy_unit to
# an eV (or h, for boundaries in Hz) and scale boltzmann_constant_ to match -- if you switch to
# quokka-style boundaries and leave these alone, x is wrong by ~12 orders of magnitude and
# PlanckFunction silently returns 0 at every group edge.
const boltzmann_constant_ = boltzmann_constant_cgs_
const energy_unit_ = 1.0		# RadSystem_Traits::energy_unit
const Erad_floor_ = 1.0e-20
const mean_molecular_mass_ = 0.6 * 1.6726231e-24 # cgs

# ISM_Traits::gas_dust_coupling_threshold -- below this the gas and dust are treated as decoupled
const gas_dust_coupling_threshold = 1.0e-5

# Stand-in opacity, used by problem.jl (see its PART 2 header)
const kappa0_ = 1.0e1 # cm^2 g^-1

# =====================================================================================================

include("hyperparameters.jl")
include("support.jl")
include("planck_integral.jl")
include("radiation_defaults.jl")
include("problem.jl")
include("opacity.jl")
include("jacobian.jl")
include("solver.jl")

# =====================================================================================================
# Driver: one call to the solver for one set of inputs.
# =====================================================================================================

function main()
	# gas and radiation state
	rho = 3.0e-8			 # g cm^-3
	T_gas = 1.0e4			 # K
	T_rad = 2.0e3			 # K -- radiation colder than the gas, so the gas heats the radiation
	dt = 1.0e-3			 # s
	dustGasCoeff = 1.0e-30		 # ISM_Traits::dust_gas_interaction_coeff

	Egas0 = ComputeEintFromTgas(rho, T_gas)
	Erad0Vec = fill(radiation_constant_ * T_rad^4 / nGroups_, nGroups_)
	rad_boundaries = [1.0e-3, 1.0e3] # nGroups_ + 1 group boundaries, in units of energy_unit_

	# coeff_n as assembled in AddSourceTermsMultiGroup (source_terms_multi_group.hpp)
	H_num_den = ComputeNumberDensityH(rho, Float64[])
	cscale = c_light_ / c_hat_
	coeff_n = dt * dustGasCoeff * H_num_den * H_num_den / cscale

	# terms the full code carries in from the rest of the timestep; all inactive here
	massScalars = Float64[]
	n_outer_iter = 0
	work = zeros(nGroups_)
	vel_times_F = zeros(nGroups_)
	Src = zeros(nGroups_)
	Q_dust = 0.0

	resid_tol = 1.0e-11
	rel_change_tol = 0.0
	tempFloor = 0.0

	p_iteration_counter = zeros(Int, 4)
	p_iteration_failure_counter = zeros(Int, 3)

	result = SolveGasDustRadiationEnergyExchange(
		Egas0, Erad0Vec, rho, coeff_n, dt, massScalars, n_outer_iter, work, vel_times_F, Src, Q_dust,
		rad_boundaries, resid_tol, rel_change_tol, tempFloor, p_iteration_counter, p_iteration_failure_counter)

	tau = dt * rho * kappa0_ * c_hat_
	E0 = Egas0 + cscale * sum(Erad0Vec)
	E1 = result.Egas + cscale * sum(result.EradVec)

	println("input   rho = $rho g/cm^3, T_gas = $T_gas K, T_rad = $T_rad K, dt = $dt s")
	println("        coeff_n = $coeff_n, tau = $(round(tau, sigdigits = 4)), nGroups = $nGroups_")
	println()
	println("T_gas   $T_gas -> $(round(result.T_gas, sigdigits = 8)) K")
	println("T_d     $(round(result.T_d, sigdigits = 8)) K")
	println("Egas    $(round(Egas0, sigdigits = 8)) -> $(round(result.Egas, sigdigits = 8)) erg/cm^3")
	println("Erad    $(round(sum(Erad0Vec), sigdigits = 8)) -> $(round(sum(result.EradVec), sigdigits = 8)) erg/cm^3")
	println("kappaP  $(result.opacity_terms.kappaP)")
	println()
	println("total energy Egas + (c/chat) Erad: relative change = $((E1 - E0) / abs(E0))")
	println("Newton iterations = $(p_iteration_counter[2]), decoupled branch = $(p_iteration_counter[4] == 1)")
	println("failures [non-convergence, negative T_d, outer] = $p_iteration_failure_counter")

	return result
end

# run the driver only when this file is executed directly, so a test can include("main.jl") to pick up
# the constants and the solver without also running it
if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
