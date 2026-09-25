# translated from RadSystem<problem_t>::SolveGasDustRadiationEnergyExchange (src/radiation/radiation_dust_system.hpp)
function SolveGasDustRadiationEnergyExchange(
    Egas0::Float64, Erad0Vec::Vector{Float64}, rho::Float64, coeff_n::Float64, dt::Float64,
    massScalars::Vector{Float64}, n_outer_iter::Int, work::Vector{Float64},
    vel_times_F::Vector{Float64}, Src::Vector{Float64}, Q_dust::Float64,
    rad_boundaries::Vector{Float64}, resid_tol::Float64, rel_change_tol::Float64, tempFloor::Float64,
    p_iteration_counter::Vector{Int}, p_iteration_failure_counter::Vector{Int}, debug::Bool = false, # NOSONAR: High cognitive complexity is expected for this numerical solver
)::NewtonIterationResult
	# 1. Compute energy exchange

	# BEGIN NEWTON-RAPHSON LOOP
	# Define the source term: S = dt chat gamma rho (kappa_P B - kappa_E E) + dt chat c^-2 gamma rho kappa_F v * F_i, where gamma =
	# 1 / sqrt(1 - v^2 / c^2) is the Lorentz factor. Solve for the new radiation energy and gas internal energy using a
	# Newton-Raphson method using the base variables (Egas, D_0, D_1,
	# ...), where D_i = R_i / tau_i^(t) and tau_i^(t) = dt * chat * gamma * rho * kappa_{P,i}^(t) is the optical depth across chat
	# * dt for group i at time t. Compared with the old base (Egas, Erad_0, Erad_1, ...), this new base is more stable and
	# converges faster. Furthermore, the PlanckOpacityTempDerivative term is not needed anymore since we assume d/dT (kappa_P /
	# kappa_E) = 0 in the calculation of the Jacobian. Note that this assumption only affects the convergence rate of the
	# Newton-Raphson iteration and does not affect the result at all once the iteration is converged.
	#
	# The Jacobian of F(E_g, D_i) is
	#
	# dF_G / dE_g = 1
	# dF_G / dD_i = c / chat * tau0_i
	# dF_{D,i} / dE_g = 1 / (chat * C_v) * (kappa_{P,i} / kappa_{E,i}) * d/dT (4 \pi B_i)
	# dF_{D,i} / dD_i = - (1 / (chat * dt * rho * kappa_{E,i}) + 1) * tau0_i = - ((1 / tau_i)(kappa_Pi / kappa_Ei) + 1) * tau0_i

	c = c_light_ # make a copy of c_light_ to avoid compiler error "undefined in device code"
	chat = c_hat_
	cscale = c / chat

	dust_model = 1
	T_d0 = NaN
	lambda_gd_times_dt = NaN
	T_gas0 = ComputeTgasFromEint(rho, Egas0, massScalars)
	amrex_assert(T_gas0 >= 0.0)
	T_d0 = ComputeDustTemperatureBateKeto(T_gas0, T_gas0, rho, Erad0Vec, coeff_n, dt, NaN, 0, Q_dust, rad_boundaries)
	amrex_assert(T_d0 >= 0.0, "Dust temperature is negative!")
	if T_d0 < 0.0
		p_iteration_failure_counter[2] += 1
	end

	max_Gamma_gd = coeff_n * std_max(std_sqrt(T_gas0) * T_gas0, std_sqrt(T_d0) * T_d0)
	# A zero collisional coupling coefficient means the gas and the dust exchange no energy at all, which is
	# exactly what the decoupled model describes, so select it on that ground alone. The comparison below
	# cannot be relied on to do it: it flips sense once Egas0 is negative, and the coupled branch it then
	# selects divides by coeff_n, filling the state with NaNs that surface later as a spurious
	# "Newton-Raphson iteration failed to converge".
	if !(coeff_n > 0.0) || (cscale * max_Gamma_gd < gas_dust_coupling_threshold * Egas0)
		dust_model = 2
		lambda_gd_times_dt = coeff_n * std_sqrt(T_gas0) * (T_gas0 - T_d0)
	end

	# Etot0 = Egas0 + cscale * (sum(Erad0Vec) + sum(Src)) + Q_dust
	Etot0 = NaN
	if dust_model == 1
		Etot0 = Egas0 + cscale * (sum(Erad0Vec) + sum(Src)) + Q_dust
	else
		# for dust_model == 2 (decoupled gas and dust), Egas0 is not involved in the iteration
		fourPiBoverC = sum(ComputeThermalRadiationMultiGroup(T_d0, rad_boundaries))
		Etot0 = abs(lambda_gd_times_dt) + fourPiBoverC + (sum(Erad0Vec) + sum(Src)) + Q_dust / cscale
	end

	T_gas = NaN
	T_d = NaN
	delta_x = NaN
	delta_R = zeros(nGroups_)
	Rvec = zeros(nGroups_)
	tau0 = zeros(nGroups_)		# optical depth across c * dt at old state
	tau = zeros(nGroups_)		# optical depth across c * dt at new state
	work_local = zeros(nGroups_)	# work term used in the Newton-Raphson iteration of the current outer iteration
	fourPiBoverC = zeros(nGroups_)
	rad_boundary_ratios = zeros(nGroups_)
	kappa_expo_and_lower_value = [zeros(nGroups_ + 1) for _ in 1:2]
	opacity_terms = OpacityTerms()

	# fill kappa_expo_and_lower_value with NaN to get warned when there are uninitialized values
	for i in 1:2
		for j in 1:(nGroups_ + 1)
			kappa_expo_and_lower_value[i][j] = NaN
		end
	end

	if !(opacity_model_ == piecewise_constant_opacity)
		for g in 1:nGroups_
			rad_boundary_ratios[g] = rad_boundaries[g + 1] / rad_boundaries[g]
		end
	end

	# define a list of alpha_quant for the model PPL_opacity_fixed_slope_spectrum
	alpha_quant_minus_one = zeros(nGroups_)
	if (opacity_model_ == PPL_opacity_fixed_slope_spectrum) ||
	   (gamma_ == 1.0 && opacity_model_ == PPL_opacity_full_spectrum)
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
	end

	Egas_guess = Egas0
	EradVec_guess = copy(Erad0Vec)

	Egas_guess_prev = Egas_guess
	EradVec_guess_prev = copy(EradVec_guess)

	H_num_den = ComputeNumberDensityH(rho, massScalars)

	T_gas = ComputeTgasFromEint(rho, Egas_guess, massScalars)
	amrex_assert(T_gas >= 0.0)

	maxIter = 100
	# Rebasing an optically thin group onto Erad_g (see newton_erad_base_tau_threshold) is applied only in
	# the decoupled branch. With dust_model == 1 the dust temperature is itself a function of sum(Rvec), and
	# ComputeJacobianForGasAndDust has already eliminated that coupling in the R_g unknowns, so the columns
	# are no longer a plain change of variable away from the Erad_g ones.
	rebase_thin = (dust_model == 2) && !use_D_as_base
	# Adaptive damping state. In the decoupled-dust branch the Newton iteration can oscillate with growing
	# amplitude rather than converge, so the step is shortened whenever the radiation residual fails to
	# decrease, and allowed to grow back towards a full step when it does. See newton_damping_* below.
	relax = 1.0
	Fg_abs_sum_prev = floatmax(Float64)
	delta_x_prev = 0.0
	delta_R_prev = zeros(nGroups_)
	iterations = SolverIterationState[]
	if debug
		# n = 0: the initial condition the loop starts from, before any Newton step. Rvec/tau are not yet
		# defined at this point (the first loop pass computes them from T_d0), and no step or residual
		# exists yet either, so those fields are left at their SolverIterationState() defaults (0 or NaN).
		# tau specifically is filled with NaN rather than 0 so it isn't mistaken for an actual zero optical
		# depth at this point (it has no defined value yet).
		push!(iterations, SolverIterationState(0, T_gas, T_d0, Egas_guess, copy(EradVec_guess), zeros(nGroups_),
							fill(NaN, nGroups_), NaN, zeros(nGroups_), NaN, zeros(nGroups_), NaN, NaN, Etot0,
							NaN, NaN, NaN, relax))
	end
	n = 0
	while n < maxIter # NOSONAR
		# if relative change is within tol, break
		if rel_change_tol > 0.0 && n > 0
			Erad_tot_guess_prev = sum(EradVec_guess_prev)
			Erad_rel_diff = abs.(EradVec_guess .- EradVec_guess_prev)
			Egas_rel_diff = abs(Egas_guess - Egas_guess_prev)

			if (sum(Erad_rel_diff) <= rel_change_tol * Erad_tot_guess_prev) && (Egas_rel_diff <= rel_change_tol * Egas_guess_prev)
				break
			end
		end

		Egas_guess_prev = Egas_guess
		EradVec_guess_prev = copy(EradVec_guess)

		# 1. Compute dust temperature
		# If the dust model is turned off, ComputeDustTemperature should be a function that returns T_gas.

		if n > 0
			T_gas = ComputeTgasFromEint(rho, Egas_guess, massScalars)
			amrex_assert(T_gas >= 0.0)
		end

		if dust_model == 1
			if n == 0
				T_d = T_d0
			else
				T_d = T_gas - (sum(Rvec) - Q_dust / cscale) / (coeff_n * std_sqrt(T_gas))
				amrex_assert(T_d >= 0.0,
					     "Dust temperature is negative! Consider increasing ISM_Traits::gas_dust_coupling_threshold")
			end
		else
			if n == 0
				T_d = T_d0
			end
		end
		if T_d < 0.0
			p_iteration_failure_counter[2] += 1
		end

		# 2. Compute kappaP and kappaE at dust temperature

		fourPiBoverC = ComputeThermalRadiationMultiGroup(T_d, rad_boundaries)

		opacity_terms = ComputeModelDependentKappaEAndKappaP(T_d, rho, rad_boundaries, rad_boundary_ratios, fourPiBoverC, EradVec_guess, n,
								    opacity_terms.alpha_E, opacity_terms.alpha_P)

		if n == 0
			# Compute kappaF and the delta_nu_kappa_B term. kappaF is used to compute the work term.
			# Will update opacity_terms in place
			ComputeModelDependentKappaFAndDeltaTerms(T_d, rho, rad_boundaries, fourPiBoverC, opacity_terms) # update opacity_terms in place
		end

		# 3. In the first loop, calculate kappaF, work, tau0, R

		if n == 0

			if (beta_order_ == 1) && (include_work_term_in_source)
				# compute the work term at the old state
				# gamma = 1.0 / sqrt(1.0 - vsqr / (c * c))
				if n_outer_iter == 0
					for g in 1:nGroups_
						if opacity_model_ == piecewise_constant_opacity
							work_local[g] = vel_times_F[g] * opacity_terms.kappaF[g] * chat / (c * c) * dt
						else
							kappa_expo_and_lower_value = DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, T_d)
							work_local[g] = vel_times_F[g] * opacity_terms.kappaF[g] * chat / (c * c) * dt *
									(1.0 + kappa_expo_and_lower_value[1][g])
						end
					end
				else
					# If n_outer_iter > 0, use the work term from the previous outer iteration, which is passed as the parameter 'work'
					work_local = copy(work)
				end
			else
				fill!(work_local, 0.0)
			end

			tau0 = dt * rho * opacity_terms.kappaP * chat
			tau = copy(tau0)
			Rvec = (fourPiBoverC .- EradVec_guess ./ opacity_terms.kappaPoverE) .* tau0 .+ work_local
			if use_D_as_base
				# tau0 is used as a scaling factor for Rvec
				for g in 1:nGroups_
					if tau0[g] <= 1.0
						tau0[g] = 1.0
					end
				end
			end
		else # in the second and later loops, calculate tau, then recover whichever of E and R is not
		     # the unknown for that group (see newton_erad_base_tau_threshold)
			tau = dt * rho * opacity_terms.kappaP * chat
			for g in 1:nGroups_
				# If tau = 0.0, Erad_guess shouldn't change
				if tau[g] > 0.0
					if rebase_thin && (tau[g] < newton_erad_base_tau_threshold)
						# Erad_g is the unknown for this thin group; R_g follows from it without
						# cancellation
						Rvec[g] = (fourPiBoverC[g] - EradVec_guess[g] / opacity_terms.kappaPoverE[g]) * tau[g] + work_local[g]
					else
						EradVec_guess[g] = opacity_terms.kappaPoverE[g] * (fourPiBoverC[g] - (Rvec[g] - work_local[g]) / tau[g])
					end
					if force_rad_floor_in_iteration
						if EradVec_guess[g] < 0.0
							Egas_guess -= cscale * (Erad_floor_ - EradVec_guess[g])
							EradVec_guess[g] = Erad_floor_
						end
					end
				end
			end
		end

		d_fourpiboverc_d_t = ComputeThermalRadiationTempDerivativeMultiGroup(T_d, rad_boundaries)
		amrex_assert(!any(isnan, d_fourpiboverc_d_t))
		c_v = ComputeEintTempDerivative(rho, T_gas, massScalars) # Egas = c_v * T

		Egas_diff = Egas_guess - Egas0
		Erad_diff = EradVec_guess .- Erad0Vec

		local jacobian

		if dust_model == 1
			jacobian = ComputeJacobianForGasAndDust(T_gas, T_d, Egas_diff, Erad_diff, Rvec, Src, Q_dust, coeff_n, tau, c_v, lambda_gd_times_dt,
								opacity_terms.kappaPoverE, d_fourpiboverc_d_t, H_num_den, dt)
		else
			jacobian = ComputeJacobianForGasAndDustDecoupled(T_gas, T_d, Egas_diff, Erad_diff, Rvec, Src, Q_dust, coeff_n, tau, c_v,
									 lambda_gd_times_dt, opacity_terms.kappaPoverE, d_fourpiboverc_d_t)
		end

		if use_D_as_base
			jacobian.J0g = jacobian.J0g .* tau0
			jacobian.Jgg = jacobian.Jgg .* tau0
		elseif rebase_thin
			RebaseThinGroupsOntoErad!(jacobian, tau, opacity_terms.kappaPoverE)
		end

		# Round-off floor on the radiation residual, as in SolveGasRadiationEnergyExchange: a group whose
		# unknown is R_g has its energy recovered from a cancelling difference, so |Fg| cannot fall below
		# the double-precision round-off of the larger operand and a purely relative test is unreachable.
		Fg_roundoff = 0.0
		for g in 1:nGroups_
			if tau[g] > 0.0
				# A group rebased onto Erad_g reaches its residual without that cancellation, so its floor
				# is just the round-off of the terms of the residual itself.
				operand =
				    (rebase_thin && (tau[g] < newton_erad_base_tau_threshold)) ?
					std_max(std_max(abs(EradVec_guess[g]), abs(Rvec[g])), std_max(abs(Erad0Vec[g]), abs(Src[g]))) :
					std_max(fourPiBoverC[g], EradVec_guess[g])
				Fg_roundoff += eps(Float64) * operand
			end
		end

		# check relative convergence of the residuals, or that the radiation residual has bottomed out at
		# the round-off floor and cannot be reduced any further
		if (abs(jacobian.F0 / Etot0) < resid_tol) &&
		   ((cscale * jacobian.Fg_abs_sum / Etot0 < resid_tol) || (jacobian.Fg_abs_sum < newton_resid_roundoff_factor * Fg_roundoff))
			break
		end

#=
		// For debugging: print (Egas0, Erad0Vec, tau0), which defines the initial condition for a Newton-Raphson iteration
		if (n == 0) {
			std::cout << "Egas0 = " << Egas0 << ", Erad0Vec = [";
			for (int g = 0; g < nGroups_; ++g) {
				std::cout << Erad0Vec[g] << ", ";
			}
			std::cout << "], tau0 = [";
			for (int g = 0; g < nGroups_; ++g) {
				std::cout << tau0[g] << ", ";
			}
			std::cout << "]";
			std::cout << "; C_V = " << c_v << ", a_rad = " << radiation_constant_ << ", coeff_n = " << coeff_n << "\n";
		} else if (n >= 0) {
			std::cout << "n = " << n << ", Egas_guess = " << Egas_guess << ", EradVec_guess = [";
			for (int g = 0; g < nGroups_; ++g) {
				std::cout << EradVec_guess[g] << ", ";
			}
			std::cout << "], tau = [";
			for (int g = 0; g < nGroups_; ++g) {
				std::cout << tau[g] << ", ";
			}
			std::cout << "]";
			std::cout << ", F_G = " << jacobian.F0 << ", F_D_abs_sum = " << jacobian.Fg_abs_sum << ", Etot0 = " << Etot0 << "\n";
		}
=#

		# update variables
		delta_x, delta_R = SolveLinearEqs(jacobian) # This is modify delta_x and delta_R in place
		amrex_assert(!isnan(delta_x))
		amrex_assert(!any(isnan, delta_R))

		# Update independent variables (Egas_guess, Rvec)
		# enable_dE_constrain is used to prevent the gas temperature from dropping/increasing below/above the radiation
		# temperature
		if dust_model == 2
			if n > 0
				# The bound is copied into a local first: std::max binds its argument by reference, and
				# that would ODR-use the namespace-scope constexpr constant, which nvcc does not make
				# available in device code.
				damping_min = newton_damping_min
				relax *= (jacobian.Fg_abs_sum > Fg_abs_sum_prev) ? newton_damping_down : newton_damping_up
				relax = std_min(std_max(relax, damping_min), 1.0)
			end
			Fg_abs_sum_prev = jacobian.Fg_abs_sum
			# Oscillation catch. When the iteration cycles about the root rather than approaching it, the
			# steps alternate in sign and each one overshoots past the root; the mean of two consecutive
			# steps is what actually points at it (for a clean period-two cycle the mean lands on it).
			# Advance by that mean instead of the raw Newton step, and let the usual convergence test
			# below decide -- this damps the cycle without bypassing the criterion.
			step_x = delta_x
			step_R = copy(delta_R)
			if n > 0 && delta_x * delta_x_prev < 0.0
				step_x = 0.5 * (delta_x + delta_x_prev)
				step_R = 0.5 * (delta_R .+ delta_R_prev)
			end
			delta_x_prev = delta_x
			delta_R_prev = copy(delta_R)
			T_d += relax * step_x
			for g in 1:nGroups_
				if rebase_thin && (tau[g] > 0.0) && (tau[g] < newton_erad_base_tau_threshold)
					# step_R holds the Erad_g step for a rebased group. Rvec is advanced to first order
					# as well, so that it stays usable if the group leaves the thin regime; when it does
					# not, Rvec is recovered exactly at the top of the next iteration.
					EradVec_guess[g] += relax * step_R[g]
					Rvec[g] += -tau[g] / opacity_terms.kappaPoverE[g] * relax * step_R[g]
				else
					Rvec[g] += relax * step_R[g]
				end
			end
		else
			T_rad = std_sqrt(std_sqrt(sum(EradVec_guess) / radiation_constant_))
			if enable_dE_constrain && delta_x / c_v > std_max(T_gas, T_rad)
				Egas_guess = ComputeEintFromTgas(rho, T_rad)
				# fill!(Rvec, 0.0)
			else
				Egas_guess += delta_x
				if use_D_as_base
					Rvec .+= tau0 .* delta_R
				else
					Rvec .+= delta_R
				end
			end
		end
		if debug
			push!(iterations, SolverIterationState(n + 1, T_gas, T_d, Egas_guess, copy(EradVec_guess), copy(Rvec),
								copy(tau), delta_x, copy(delta_R), jacobian.F0, copy(jacobian.Fg),
								jacobian.Fg_abs_sum, Fg_roundoff, Etot0, abs(jacobian.F0 / Etot0),
								cscale * jacobian.Fg_abs_sum / Etot0, jacobian.Fg_abs_sum / Fg_roundoff, relax))
		end
		n += 1
		# check relative and absolute convergence of E_r
		# if (std::abs(deltaEgas / Egas_guess) < 1e-7) {
		# 	break;
		# }
	end # END NEWTON-RAPHSON LOOP

	# Inject the source directly into transparent groups (tau ~ 0). The Newton solve above excludes such
	# groups from its residual and Jacobian (Fg_abs_sum and Jgg skip tau <= 0) and leaves their radiation
	# energy at Erad0, so an injected source in a transparent group would otherwise be silently dropped
	# while UpdateFlux still applies the matching flux source, leaving |F| > c E. This mirrors the loop in
	# the gas-only solver (source_terms_multi_group.hpp) and the single-group negligible-optical-depth
	# branch; Src is already counted in Etot0, so it is energy-consistent. Groups with tau > 0 (the usual
	# case) and groups without a source are unaffected.
	for g in 1:nGroups_
		if !(tau[g] > 0.0)
			EradVec_guess[g] = Erad0Vec[g] + Src[g]
		end
	end

	cooling_tend = DefineNetCoolingRate(T_gas, H_num_den) * dt
	if dust_model == 2
		# include line cooling/heating, cosmic ray heating terms; implicitly update Egas_guess

		CR_heating = DefineCosmicRayHeatingRate(H_num_den) * dt

		# Sum the magnitudes of the terms in the residual: lambda_gd_times_dt is signed, so the raw sum could
		# cancel to zero and leave the convergence test unsatisfiable. Egas_guess > 0 keeps the scale positive.
		compare = Egas_guess + abs(cscale * lambda_gd_times_dt) + sum(abs.(cooling_tend)) + abs(CR_heating)

		# RHS of the equation 0 = Egas - Egas0 + cscale * lambda_gd_times_dt + sum(cooling)
		rhs = function (Egas_::Float64)
			T_gas_ = ComputeTgasFromEint(rho, Egas_, massScalars)
			cooling_ = DefineNetCoolingRate(T_gas_, H_num_den) * dt
			return Egas_ - Egas0 + cscale * lambda_gd_times_dt + sum(cooling_) - CR_heating
		end

		# Jacobian of the RHS of the equation 0 = Egas - Egas0 + cscale * lambda_gd_times_dt + sum(cooling)
		jac = function (Egas_::Float64)
			T_gas_ = ComputeTgasFromEint(rho, Egas_, massScalars)
			d_cooling_d_Tgas_ = DefineNetCoolingRateTempDerivative(T_gas_, H_num_den) * dt
			c_v_ = ComputeEintTempDerivative(rho, T_gas_, massScalars) # Egas = c_v * T
			# The residual is a function of Egas_, so convert dCooling/dT to dCooling/dEgas via the chain rule (dT/dEgas = 1 / c_v).
			return 1.0 + sum(d_cooling_d_Tgas_) / c_v_
		end

		Egas_guess = BackwardEulerOneVariable(rhs, jac, Egas0, compare)
	end

	if !add_line_cooling_to_radiation_in_jac
		amrex_assert(minimum(cooling_tend) >= 0.0, "add_line_cooling_to_radiation has to be enabled when there is negative cooling rate!")
		# TODO(CCH): potential GPU-related issue here.
		EradVec_guess .+= (1 / cscale) * cooling_tend
	end

	amrex_assert(Egas_guess > 0.0)
	amrex_assert(minimum(EradVec_guess) >= 0.0)

	amrex_assert(n < maxIter, "Newton-Raphson iteration for matter-radiation coupling failed to converge!")
	if n >= maxIter
		p_iteration_failure_counter[1] += 1
	end

	p_iteration_counter[1] += 1			       # total number of radiation updates.
	p_iteration_counter[2] += n + 1			       # total number of Newton-Raphson iterations.
	p_iteration_counter[3] = max(p_iteration_counter[3], n + 1) # maximum number of Newton-Raphson iterations.
	if dust_model == 2
		p_iteration_counter[4] += 1 # total number of decoupled gas-dust iterations.
	end

	result = NewtonIterationResult()

	if n > 0
		# calculate kappaF since the temperature has changed
		# Will update opacity_terms in place
		ComputeModelDependentKappaFAndDeltaTerms(T_d, rho, rad_boundaries, fourPiBoverC, opacity_terms) # update opacity_terms in place
	end

	result.Egas = Egas_guess
	result.EradVec = copy(EradVec_guess)
	result.work = copy(work_local)
	result.T_gas = T_gas
	result.T_d = T_d
	result.opacity_terms = opacity_terms
	result.iterations = iterations
	return result
end
