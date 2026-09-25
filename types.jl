# Types shared by the radiation solver. Ports of the structs in src/radiation/radiation_system.hpp.
# The C++ structs are templated on problem_t and sized by Physics_Traits<problem_t>::nGroups; here
# they are plain mutable structs of length-nGroups_ vectors. They must be mutable: the solver assigns
# fields one at a time and RebaseThinGroupsOntoErad updates a Jacobian in place.

# Integer values match the C++ enum in radiation_system.hpp exactly, so a dumped opacity_model can be
# mapped straight through with OpacityModel(i).
@enum OpacityModel begin
	single_group = 0
	piecewise_constant_opacity = 1
	PPL_opacity_fixed_slope_spectrum = 2
	PPL_opacity_full_spectrum = 3
end

# A struct to hold the opacity terms for the radiation-matter energy exchange, containing the following
# elements: kappaE, kappaP, kappaF, kappaPoverE, delta_nu_kappa_B_at_edge, alpha_P, alpha_E
mutable struct OpacityTerms
	kappaE::Vector{Float64}
	kappaP::Vector{Float64}
	kappaF::Vector{Float64}
	kappaPoverE::Vector{Float64}
	delta_nu_kappa_B_at_edge::Vector{Float64} # Delta (nu * kappa * B)
	alpha_P::Vector{Float64}
	alpha_E::Vector{Float64}
end

OpacityTerms() = OpacityTerms(zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_))

# A struct to hold the results of the ComputeJacobian functions, containing the following elements:
# J00, F0, Fg_abs_sum, J0g, Jg0, Jgg, Fg
mutable struct JacobianResult
	J00::Float64		# (0, 0) component of the Jacobian matrix
	F0::Float64		# (0) component of the residual
	Fg_abs_sum::Float64	# sum of the absolute values of the (g) components of the residual, with tau(g) > 0
	J0g::Vector{Float64}	# (0, g) components of the Jacobian matrix
	Jg0::Vector{Float64}	# (g, 0) components of the Jacobian matrix
	Jgg::Vector{Float64}	# (g, g) components of the Jacobian matrix
	Jg1::Vector{Float64}	# (g, 1) components of the Jacobian matrix
	Fg::Vector{Float64}	# (g) components of the residual
end

JacobianResult() = JacobianResult(0.0, 0.0, 0.0, zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_), zeros(nGroups_))

# Snapshot of one pass through the Newton-Raphson loop in SolveGasDustRadiationEnergyExchange. Not a
# port of any C++ struct -- the loop there has no per-iteration record, it just overwrites its locals
# each pass. This exists purely for inspecting/debugging the Julia solver's convergence history.
mutable struct SolverIterationState
	n::Int				# n = 0 is the initial condition; n = k > 0 is the state after Newton step k
	T_gas::Float64
	T_d::Float64
	Egas_guess::Float64
	EradVec_guess::Vector{Float64}
	Rvec::Vector{Float64}		# R_g, the radiation exchange term (or Newton unknown for thin groups)
	tau::Vector{Float64}		# optical depth across chat * dt at the new state
	delta_x::Float64		# Newton step in the gas-energy (or dust-temperature) unknown
	delta_R::Vector{Float64}	# Newton step in the R_g (or Erad_g, for rebased groups) unknowns
	F0::Float64			# gas/dust-temperature residual
	Fg::Vector{Float64}		# per-group residual (g) components, before taking abs and summing
	Fg_abs_sum::Float64		# sum of |Fg| over groups with tau > 0
	Fg_roundoff::Float64		# round-off floor on Fg_abs_sum, used in the convergence check (solver.jl)
	Etot0::Float64			# total energy scale used to non-dimensionalize the convergence check (solver.jl)
	F0_resid_ratio::Float64	# abs(F0 / Etot0), compared against resid_tol
	Fg_resid_ratio::Float64	# cscale * Fg_abs_sum / Etot0, compared against resid_tol
	Fg_roundoff_ratio::Float64	# Fg_abs_sum / Fg_roundoff, compared against newton_resid_roundoff_factor
	relax::Float64			# damping factor applied to the step (decoupled dust branch only)
end

SolverIterationState() = SolverIterationState(0, NaN, NaN, NaN, zeros(nGroups_), zeros(nGroups_), zeros(nGroups_),
					      NaN, zeros(nGroups_), NaN, zeros(nGroups_), NaN, NaN, NaN, NaN, NaN, NaN)

# A struct to hold the results of the Newton-Raphson iteration for energy update, containing the
# following elements: Egas, T_gas, T_d, EradVec, work, opacity_terms
mutable struct NewtonIterationResult
	Egas::Float64			# gas internal energy
	T_gas::Float64			# gas temperature
	T_d::Float64			# dust temperature
	EradVec::Vector{Float64}	# radiation energy density
	work::Vector{Float64}		# work term
	opacity_terms::OpacityTerms
	iterations::Vector{SolverIterationState} # per-iteration history, debug builds only (see solver.jl)
end

NewtonIterationResult() = NewtonIterationResult(NaN, NaN, NaN, zeros(nGroups_), zeros(nGroups_), OpacityTerms(), SolverIterationState[])
