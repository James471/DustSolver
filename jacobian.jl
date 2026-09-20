# Jacobians of the energy-update equations, and the linear solve for their sparse structure.
# Ports of src/radiation/radiation_dust_system.hpp (ComputeJacobianForGasAndDust*) and
# src/radiation/radiation_system.hpp (RebaseThinGroupsOntoErad, SolveLinearEqs).

# Compute the Jacobian of energy update equations for the gas-dust-radiation system. The result is a struct containing the following elements:
# J00: (0, 0) component of the Jacobian matrix. = d F0 / d Egas
# F0: (0) component of the residual. = Egas residual
# Fg_abs_sum: sum of the absolute values of the each component of Fg that has tau(g) > 0
# J0g: (0, g) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d F0 / d R_g
# Jg0: (g, 0) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d Fg / d Egas
# Jgg: (g, g) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d Fg / d R_g
# Fg: (g) components of the residual, g = 1, 2, ..., nGroups. = Erad residual
function ComputeJacobianForGasAndDust(T_gas, T_d, Egas_diff, Erad_diff, Rvec, Src, Q_dust, coeff_n, tau, c_v,
				      lambda_gd_time_dt, kappaPoverE, d_fourpiboverc_d_t, num_den, dt)::JacobianResult
	result = JacobianResult()

	cscale = c_light_ / c_hat_

	# compute cooling/heating terms
	cooling = DefineNetCoolingRate(T_gas, num_den) * dt
	cooling_derivative = DefineNetCoolingRateTempDerivative(T_gas, num_den) * dt
	CR_heating = DefineCosmicRayHeatingRate(num_den) * dt

	# Q_dust is already c/c_hat scaled.
	result.F0 = Egas_diff + cscale * sum(Rvec) + sum(cooling) - CR_heating - Q_dust
	result.Fg = Erad_diff .- (Rvec .+ Src)
	if add_line_cooling_to_radiation_in_jac
		result.Fg .-= (1.0 / cscale) * cooling
	end
	result.Fg_abs_sum = 0.0
	for g in 1:nGroups_
		if tau[g] > 0.0
			result.Fg_abs_sum += abs(result.Fg[g])
		else
			result.Fg_abs_sum += abs(result.Fg[g] + Rvec[g])
		end
	end

	# compute Jacobian elements
	# I assume (kappaPVec / kappaEVec) is constant here. This is usually a reasonable assumption. Note that this assumption
	# only affects the convergence rate of the Newton-Raphson iteration and does not affect the converged solution at all.

	dEg_dT = kappaPoverE .* d_fourpiboverc_d_t

	result.J00 = 1.0 + sum(cooling_derivative) / c_v
	fill!(result.J0g, cscale)
	d_Td_d_T = 3.0 / 2.0 - T_d / (2.0 * T_gas)
	dEg_dT = dEg_dT * d_Td_d_T
	dTd_dRg = -1.0 / (coeff_n * std_sqrt(T_gas))
	rg = kappaPoverE .* d_fourpiboverc_d_t * dTd_dRg
	result.Jg0 = 1.0 / c_v * dEg_dT .- (1 / cscale) * cooling_derivative .- 1.0 / cscale * rg * result.J00
	# Note that Fg is modified here, but it does not change Fg_abs_sum, which is used to check the convergence.
	result.Fg = result.Fg .- 1.0 / cscale * rg * result.F0
	for g in 1:nGroups_
		if tau[g] <= 0.0
			result.Jgg[g] = -Inf
		else
			result.Jgg[g] = -1.0 * kappaPoverE[g] / tau[g] - 1.0
		end
	end

	return result
end

# Compute the Jacobian of energy update equations for the gas-dust-radiation system with gas and dust decoupled. The result is a struct containing the
# following elements: J00: (0, 0) component of the Jacobian matrix. = d F0 / d T_d F0: (0) component of the residual. = sum_g R_g - lambda_gd_time_dt
# Fg_abs_sum: sum of the absolute values of the each component of Fg that has tau(g) > 0
# J0g: (0, g) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d F0 / d R_g
# Jg0: (g, 0) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d Fg / d T_d
# Jgg: (g, g) components of the Jacobian matrix, g = 1, 2, ..., nGroups. = d Fg / d R_g
# Fg: (g) components of the residual, g = 1, 2, ..., nGroups. = Erad residual
function ComputeJacobianForGasAndDustDecoupled(T_gas, T_d, Egas_diff, Erad_diff, Rvec, Src, Q_dust, coeff_n, tau, c_v,
					       lambda_gd_time_dt, kappaPoverE, d_fourpiboverc_d_t)::JacobianResult
	c_hat_over_c = c_hat_ / c_light_
	result = JacobianResult()

	# lambda_gd_time_dt and Rvec are c_hat/c scaled. So, we are solving sum(R_g) - lambda_gd_time_dt - c_hat / c * Q_dust = 0, because Q_dust is already c
	# / c_hat scaled. We could have equivalently moved the c_hat / c factor to the other terms.
	result.F0 = -lambda_gd_time_dt + sum(Rvec) - c_hat_over_c * Q_dust
	result.Fg = Erad_diff .- (Rvec .+ Src)
	result.Fg_abs_sum = 0.0
	for g in 1:nGroups_
		if tau[g] > 0.0
			result.Fg_abs_sum += abs(result.Fg[g])
		end
	end

	# compute Jacobian elements
	# I assume (kappaPVec / kappaEVec) is constant here. This is usually a reasonable assumption. Note that this assumption
	# only affects the convergence rate of the Newton-Raphson iteration and does not affect the converged solution at all.

	dEg_dT = kappaPoverE .* d_fourpiboverc_d_t

	result.J00 = 0.0
	fill!(result.J0g, 1.0)
	result.Jg0 = dEg_dT
	for g in 1:nGroups_
		if tau[g] <= 0.0
			result.Jgg[g] = -Inf
		else
			result.Jgg[g] = -1.0 * kappaPoverE[g] / tau[g] - 1.0
		end
	end

	return result
end

# Rebase the optically thin groups of a Jacobian built in the R_g unknowns onto the Erad_g unknowns.
# See newton_erad_base_tau_threshold for why. The map is R_g = (4 pi B_g / c - Erad_g / kappaPoverE_g) *
# tau_g + w_g, so dR_g/dErad_g = -tau_g / kappaPoverE_g: the column of group g is scaled by that factor.
# The (g, 0) entry is scaled too, because dF_g/dx at fixed Erad_g is not the same partial derivative as
# dF_g/dx at fixed R_g. The (0, 0) entry instead gains a term, since R_g now varies with x: F0 contains
# J0g[g] * R_g, contributing J0g[g] * dR_g/dx = -J0g[g] * scale * Jg0[g]. The residuals are unchanged --
# they are the same numbers regardless of which variable is held independent.
function RebaseThinGroupsOntoErad!(jacobian::JacobianResult, tau::Vector{Float64}, kappaPoverE::Vector{Float64})
	for g in 1:nGroups_
		if (tau[g] > 0.0) && (tau[g] < newton_erad_base_tau_threshold)
			scale = -tau[g] / kappaPoverE[g] # dR_g / dErad_g
			jacobian.J00 -= jacobian.J0g[g] * scale * jacobian.Jg0[g]
			jacobian.Jgg[g] *= scale
			jacobian.J0g[g] *= scale
			jacobian.Jg0[g] *= scale
			println("scale: ", scale)
		end
	end
end

# Linear equation solver for matrix with non-zeros at the first row, first column, and diagonal only.
# solve the linear system
#   [a00 a0i] [x0] = [y0]
#   [ai0 aii] [xi]   [yi]
# for x0 and xi, where a0i = (a01, a02, ...); ai0 = (a10, a20, ...); aii = (a11, a22, ...)
# The C++ version writes x0 and xi through out-parameters; Julia returns them instead.
function SolveLinearEqs(jacobian::JacobianResult)
	ratios = jacobian.J0g ./ jacobian.Jgg
	x0 = (sum(ratios .* jacobian.Fg) - jacobian.F0) / (-sum(ratios .* jacobian.Jg0) + jacobian.J00)
	xi = (-1.0 * jacobian.Fg .- jacobian.Jg0 * x0) ./ jacobian.Jgg
	println("jacobian.J00: ", jacobian.J00)
	println("jacobian.J0g: ", jacobian.J0g)
	println("jacobian.Jg0: ", jacobian.Jg0)
	println("jacobian.Jgg: ", jacobian.Jgg)
	println("ratios: ", ratios)
	println("jacobian.Fg: ", jacobian.Fg)
	println("jacobian.F0: ", jacobian.F0)
	println("radios .* jacobian.Fg: ", ratios .* jacobian.Fg)
	println("sum(ratios .* jacobian.Fg): ", sum(ratios .* jacobian.Fg))
	println("sum(ratios .* jacobian.Jg0): ", sum(ratios .* jacobian.Jg0))
	println("num = sum(ratios .* jacobian.Fg) - jacobian.F0: ", sum(ratios .* jacobian.Fg) - jacobian.F0)
	println("den = -sum(ratios .* jacobian.Jg0) + jacobian.J00: ", -sum(ratios .* jacobian.Jg0) + jacobian.J00)
	println("x0: ", x0)
	println("xi: ", xi)
	return x0, xi
end
