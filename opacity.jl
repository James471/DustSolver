# Group-mean opacity assembly and the dust-temperature root solve.
# Ports of src/radiation/source_terms_multi_group.hpp (ComputeModelDependentKappa*) and
# src/radiation/radiation_system.hpp (ComputeDustTemperatureBateKeto).

# Compute kappaE and kappaP based on the opacity model. The result is stored in the last five arguments: alpha_P, alpha_E, kappaP, kappaE, and kappaPoverE.
function ComputeModelDependentKappaEAndKappaP(T, rho, rad_boundaries, rad_boundary_ratios, fourPiBoverC, Erad, n_iter,
					      alpha_E = zeros(nGroups_), alpha_P = zeros(nGroups_))::OpacityTerms
	result = OpacityTerms()

	kappa_expo_and_lower_value = DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, T)

	if opacity_model_ == piecewise_constant_opacity
		for g in 1:nGroups_
			result.kappaP[g] = kappa_expo_and_lower_value[2][g]
			result.kappaE[g] = kappa_expo_and_lower_value[2][g]
		end
	elseif opacity_model_ == PPL_opacity_fixed_slope_spectrum
		alpha_quant_minus_one = zeros(nGroups_)
		if !special_edge_bin_slopes
			for g in 1:nGroups_
				alpha_quant_minus_one[g] = -1.0
			end
		else
			alpha_quant_minus_one[1] = 2.0
			alpha_quant_minus_one[nGroups_] = -4.0
			for g in 2:(nGroups_ - 1)
				alpha_quant_minus_one[g] = -1.0
			end
		end
		result.kappaP = ComputeGroupMeanOpacity(kappa_expo_and_lower_value, rad_boundary_ratios, alpha_quant_minus_one)
		result.kappaE = copy(result.kappaP)
	elseif opacity_model_ == PPL_opacity_full_spectrum
		if n_iter < max_iter_to_update_alpha_E
			result.alpha_E = ComputeRadQuantityExponents(Erad, rad_boundaries)
			result.alpha_P = ComputeRadQuantityExponents(fourPiBoverC, rad_boundaries)
		else
			result.alpha_E = copy(alpha_E)
			result.alpha_P = copy(alpha_P)
		end
		result.kappaE = ComputeGroupMeanOpacity(kappa_expo_and_lower_value, rad_boundary_ratios, result.alpha_E)
		result.kappaP = ComputeGroupMeanOpacity(kappa_expo_and_lower_value, rad_boundary_ratios, result.alpha_P)
	end
	amrex_assert(!any(isnan, result.kappaP))
	amrex_assert(!any(isnan, result.kappaE))
	for g in 1:nGroups_
		if result.kappaE[g] > 0.0
			result.kappaPoverE[g] = result.kappaP[g] / result.kappaE[g]
		else
			result.kappaPoverE[g] = 1.0
		end
	end

	return result
end

# Compute kappaF and the delta_nu_kappa_B_at_edge term. kappaF is used to compute the work term and the delta_nu_kappa_B_at_edge term is used to compute the
# transport between groups in the momentum function. Only opacity_terms is modified in this function.
function ComputeModelDependentKappaFAndDeltaTerms(T, rho, rad_boundaries, fourPiBoverC, opacity_terms::OpacityTerms)
	delta_nu_B_at_edge = zeros(nGroups_)
	kappa_expo_and_lower_value = DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, T)
	for g in 1:nGroups_
		nu_L = rad_boundaries[g]
		nu_R = rad_boundaries[g + 1]
		B_L = PlanckFunction(nu_L, T) # 4 pi B(nu) / c
		B_R = PlanckFunction(nu_R, T) # 4 pi B(nu) / c
		kappa_L = kappa_expo_and_lower_value[2][g]
		kappa_R = kappa_L * std_pow(nu_R / nu_L, kappa_expo_and_lower_value[1][g]) # C++ std::pow
		opacity_terms.delta_nu_kappa_B_at_edge[g] = nu_R * kappa_R * B_R - nu_L * kappa_L * B_L
		delta_nu_B_at_edge[g] = nu_R * B_R - nu_L * B_L
	end
	if opacity_model_ == piecewise_constant_opacity
		opacity_terms.kappaF = copy(opacity_terms.kappaP)
	else
		if use_diffuse_flux_mean_opacity
			opacity_terms.kappaF =
			    ComputeDiffusionFluxMeanOpacity(opacity_terms.kappaP, opacity_terms.kappaE, fourPiBoverC, opacity_terms.delta_nu_kappa_B_at_edge,
							    delta_nu_B_at_edge, kappa_expo_and_lower_value[1])
		else
			# for simplicity, I assume kappaF = kappaE when opacity_model_ ==
			# OpacityModel::PPL_opacity_full_spectrum, if !use_diffuse_flux_mean_opacity. We won't use this
			# option anyway.
			opacity_terms.kappaF = copy(opacity_terms.kappaE)
		end
	end
end

function ComputeDustTemperatureBateKeto(T_gas, T_d_init, rho, Erad, N_d, dt, R_sum, n_step, Q_dust, rad_boundaries)::Float64
	c_hat_over_c = c_hat_ / c_light_
	if n_step > 0
		T_d = T_gas - (R_sum - c_hat_over_c * Q_dust) / (N_d * std_sqrt(T_gas))
		amrex_assert(T_d >= 0.0, "Dust temperature is negative!")
		return T_d
	end

	rad_boundary_ratios = zeros(nGroups_)

	if nGroups_ > 1 && opacity_model_ != piecewise_constant_opacity
		for g in 1:nGroups_
			rad_boundary_ratios[g] = rad_boundaries[g + 1] / rad_boundaries[g]
		end
	end

	# the RHS of the equation 0 = c_hat_ dt rho (kappa_E * E_g - kappa_P * B_g) + N_d sqrt(T_gas) (T_gas - T_d) + Q_dust
	rhs = function (T_d)
		LHS = NaN

		if nGroups_ == 1
			fourPiBoverC = ComputeThermalRadiationSingleGroup(T_d)
			kappaE = ComputeEnergyMeanOpacity(rho, T_d)
			kappaP = ComputePlanckOpacity(rho, T_d)
			LHS = c_hat_ * dt * rho * (kappaE * Erad[1] - kappaP * fourPiBoverC) + N_d * std_sqrt(T_gas) * (T_gas - T_d) + Q_dust * c_hat_over_c
		else
			fourPiBoverC = ComputeThermalRadiationMultiGroup(T_d, rad_boundaries)
			opacity_terms = ComputeModelDependentKappaEAndKappaP(T_d, rho, rad_boundaries, rad_boundary_ratios, fourPiBoverC, Erad, 0)
			LHS = c_hat_ * dt * rho * sum(opacity_terms.kappaE .* Erad .- opacity_terms.kappaP .* fourPiBoverC) +
			      N_d * std_sqrt(T_gas) * (T_gas - T_d) + Q_dust * c_hat_over_c
		end

		return LHS
	end

	# the Jacobian of the RHS of the equation 0 = c_hat_ dt rho (kappa_E * E_g - kappa_P * B_g) + N_d sqrt(T_gas) (T_gas - T_d) + Q_dust
	jac = function (T_d)
		dLHS_dTd = NaN

		if nGroups_ == 1
			kappaP = ComputePlanckOpacity(rho, T_d)
			d_fourpib_over_c_d_t = ComputeThermalRadiationTempDerivativeSingleGroup(T_d)
			dLHS_dTd = -c_hat_ * dt * rho * (kappaP * d_fourpib_over_c_d_t) - N_d * std_sqrt(T_gas)
		else
			fourPiBoverC = ComputeThermalRadiationMultiGroup(T_d, rad_boundaries)
			opacity_terms = ComputeModelDependentKappaEAndKappaP(T_d, rho, rad_boundaries, rad_boundary_ratios, fourPiBoverC, Erad, 0)
			d_fourpib_over_c_d_t = ComputeThermalRadiationTempDerivativeMultiGroup(T_d, rad_boundaries)
			dLHS_dTd = -c_hat_ * dt * rho * sum(opacity_terms.kappaP .* d_fourpib_over_c_d_t) - N_d * std_sqrt(T_gas)
		end

		return dLHS_dTd
	end

	# Scale for the convergence test. The residual balances the radiative term against the gas-dust collisional
	# term, so the scale must contain both: the collisional term vanishes identically when the gas-dust coupling
	# coefficient N_d is zero, and a scale of zero would make the convergence test unsatisfiable.
	Lambda_compare = N_d * std_sqrt(T_gas) * T_gas
	if nGroups_ == 1
		fourPiBoverC = ComputeThermalRadiationSingleGroup(T_d_init)
		kappaE = ComputeEnergyMeanOpacity(rho, T_d_init)
		kappaP = ComputePlanckOpacity(rho, T_d_init)
		Lambda_compare += c_hat_ * dt * rho * (kappaE * Erad[1] + kappaP * fourPiBoverC)
	else
		fourPiBoverC = ComputeThermalRadiationMultiGroup(T_d_init, rad_boundaries)
		opacity_terms = ComputeModelDependentKappaEAndKappaP(T_d_init, rho, rad_boundaries, rad_boundary_ratios, fourPiBoverC, Erad, 0)
		Lambda_compare += c_hat_ * dt * rho * sum(opacity_terms.kappaE .* Erad .+ opacity_terms.kappaP .* fourPiBoverC)
	end

	T_d = BackwardEulerOneVariable(rhs, jac, T_d_init, Lambda_compare)
	amrex_assert(T_d >= 0.0, "Dust temperature is negative!")

	return T_d
end
