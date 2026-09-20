# =====================================================================================================
#                          READ THIS BEFORE COMPARING NUMBERS AGAINST QUOKKA
#
# Everything in this file is a STAND-IN. It replaces code that each quokka problem defines for itself:
# the EOS comes from hydro/EOS.hpp, the opacity from the problem's own
# DefineOpacityExponentsAndLowerValues. None of it is a port of anything and none of it is a reference
# implementation -- replace it with your problem's real versions before trusting any number that
# depends on it. repro_DTypeFront3D.jl shows what a real problem's definitions look like.
#
# The problem-independent quokka defaults (Planck, cooling) ARE ports, and live in radiation_defaults.jl.
# =====================================================================================================

# Ideal-gas EOS, standing in for quokka::EOS<problem_t> (hydro/EOS.hpp), which may carry a full
# Microphysics network. Egas = rho * c_v * T with c_v = k_B / (mu * m_H * (gamma - 1)).
cv_specific_() = boltzmann_constant_cgs_ / (mean_molecular_mass_ * (gamma_ - 1.0))
ComputeTgasFromEint(rho, Eint, massScalars) = Eint / (rho * cv_specific_())
ComputeEintFromTgas(rho, Tgas) = rho * cv_specific_() * Tgas
ComputeEintTempDerivative(rho, Tgas, massScalars) = rho * cv_specific_()

# Constant opacity, standing in for the problem-defined DefineOpacityExponentsAndLowerValues (the
# quokka default returns NaN, forcing each problem to supply its own). Row 1 holds the power-law
# exponent across each group, row 2 the opacity at the lower group boundary; a flat spectrum means
# exponent 0 and a single constant value.
function DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, Tgas)
	return [zeros(nGroups_ + 1), fill(kappa0_, nGroups_ + 1)]
end

# Single-group opacities, standing in for the problem-defined ComputePlanckOpacity (the quokka default
# returns NaN). ComputeEnergyMeanOpacity defaults to ComputePlanckOpacity.
ComputePlanckOpacity(rho, Tgas) = kappa0_
ComputeEnergyMeanOpacity(rho, Tgas) = ComputePlanckOpacity(rho, Tgas)
