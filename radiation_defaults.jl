# Ports of the problem-independent quokka defaults in src/radiation/radiation_system.hpp. Faithful
# translations; compare them line-for-line. Anything a quokka problem is expected to define for
# itself (EOS, opacity) lives in the driver instead -- see problem.jl.

# rho / mean_molecular_mass_
ComputeNumberDensityH(rho, massScalars) = rho / mean_molecular_mass_

# define ComputeThermalRadiation for single-group, returns the thermal radiation power = a_r * T^4
function ComputeThermalRadiationSingleGroup(temperature)::Float64
	power = radiation_constant_ * std_pow(temperature, 4.0) # C++ std::pow(temperature, 4)
	# set floor
	if power < Erad_floor_
		power = Erad_floor_
	end
	return power
end

# by default, d emission/dT = 4 emission / T
ComputeThermalRadiationTempDerivativeSingleGroup(temperature) = 4.0 * radiation_constant_ * temperature^3

# Compute radiation energy fractions for each photon group from a Planck function, given nGroups,
# radBoundaries, and temperature. Uses the ported lookup table in planck_integral.jl.
function ComputePlanckEnergyFractions(boundaries, temperature)::Vector{Float64}
	radEnergyFractions = zeros(nGroups_)
	if nGroups_ == 1
		radEnergyFractions[1] = 1.0
		return radEnergyFractions
	else
		energy_unit_over_kT = energy_unit_ / (boltzmann_constant_ * temperature)
		# (P, Q) = lower/upper Planck integrals at the group's lower edge; the first group's lower edge is
		# taken as x = 0.
		p_prev = 0.0
		q_prev = 1.0
		prev_in_tail = false
		# Only the thermal groups (the leading nGroupsThermal_ groups) receive blackbody emission. When
		# chemical bands are present the thermal fractions are NOT renormalized: the blackbody radiation
		# above the first chemical-band boundary is simply dropped, so the fractions sum to < 1.
		for g in 1:nGroupsThermal_
			in_tail = true
			if g == nGroups_
				# no chemical bands: the last group carries all remaining blackbody, total fraction = 1.0
				p, q = 1.0, 0.0
			else
				x = boundaries[g + 1] * energy_unit_over_kT
				in_tail = x >= X_TAIL
				p, q = integrate_planck_below_and_above_x(x)
			end
			# See X_TAIL in planck_integral.jl (ported to quokka-james471): a group whose lower edge is in
			# the Wien tail takes its fraction as a difference of upper integrals, Q(x_lo) - Q(x_hi), not of
			# lower integrals that are both ~1. Otherwise P(x_hi) - P(x_lo), as before.
			radEnergyFractions[g] = prev_in_tail ? q_prev - q : p - p_prev
			p_prev, q_prev, prev_in_tail = p, q, in_tail
		end
		# chemical bands (g > nGroupsThermal_) emit no blackbody radiation; left at 0.
		amrex_assert(sum(radEnergyFractions) < 1.0 + 1.0e-10)

		return radEnergyFractions
	end
end

# ComputeThermalRadiationMultiGroup, returns the thermal radiation power for each photon group.
# = a_r * T^4 * radEnergyFractions
function ComputeThermalRadiationMultiGroup(temperature, boundaries)::Vector{Float64}
	power = radiation_constant_ * std_pow(temperature, 4.0) # C++ std::pow(temperature, 4)
	radEnergyFractions = ComputePlanckEnergyFractions(boundaries, temperature)
	Erad_g = power * radEnergyFractions
	# set floor on the thermal groups only; chemical bands emit no blackbody radiation and are left at 0.
	for g in 1:nGroupsThermal_
		if Erad_g[g] < Erad_floor_
			Erad_g[g] = Erad_floor_
		end
	end
	return Erad_g
end

function ComputeThermalRadiationTempDerivativeMultiGroup(temperature, boundaries)::Vector{Float64}
	d_fourpiboverc_d_t = zeros(nGroups_)
	a_T3 = radiation_constant_ * temperature * temperature * temperature
	if nGroups_ == 1
		d_fourpiboverc_d_t[1] = 4.0 * a_T3
		return d_fourpiboverc_d_t
	else
		# The group emission is a T^4 f_g(T), so its temperature derivative is 4 a T^3 f_g + a T^4 df_g/dT.
		# Integrating the exact kernel by parts gives the cumulative form
		#     D(x) = (15/pi^4) \int_0^x s^4 e^s / (e^s - 1)^2 ds = 4 P(x) - (15/pi^4) x^4 / (e^x - 1),
		# where P is the same normalized Planck integral used for the energy fractions. D(inf) = 4
		# recovers d(a T^4)/dT.
		# See X_TAIL in planck_integral.jl (ported to quokka-james471): alongside D(x) this also tracks its complement
		#     U(x) = 4 - D(x) = 4 Q(x) + (15/pi^4) x^4 / (e^x - 1),
		# and a group whose lower edge is in the Wien tail takes U(x_lo) - U(x_hi) instead of
		# D(x_hi) - D(x_lo), mirroring ComputePlanckEnergyFractions. There Q comes from the exact tail
		# series, so the derivative is that of the B actually emitted, not of the true Planck function
		# that the table only approximates.
		energy_unit_over_kT = energy_unit_ / (boltzmann_constant_ * temperature)
		d_prev = 0.0
		u_prev = 4.0
		prev_in_tail = false
		# Only the thermal groups emit; the chemical bands are left at 0, as in ComputePlanckEnergyFractions.
		for g in 1:nGroupsThermal_
			in_tail = true
			if g == nGroups_
				# no chemical bands: the last group carries all remaining blackbody, so D = D(inf) = 4
				d, u = 4.0, 0.0
			else
				x = boundaries[g + 1] * energy_unit_over_kT
				if x >= X_TAIL
					q = integrate_planck_from_x_to_inf_series(x)
					kernel = (x * x * x * x / (std_exp(x) - 1.0)) / gInf
					d, u = 4.0 * (1.0 - q) - kernel, 4.0 * q + kernel
				else
					in_tail = false
					d = 4.0 * integrate_planck_from_0_to_x(x) - (x * x * x * x / (std_exp(x) - 1.0)) / gInf
					u = 4.0 - d
				end
			end
			d_fourpiboverc_d_t[g] = a_T3 * (prev_in_tail ? u_prev - u : d - d_prev)
			d_prev, u_prev, prev_in_tail = d, u, in_tail
		end

		return d_fourpiboverc_d_t
	end
end

# returns 4 pi B(nu) / c
function PlanckFunction(nu, T)::Float64
	coeff = energy_unit_ / (boltzmann_constant_ * T)
	x = coeff * nu
	if x > 100.0
		return 0.0
	end
	planck_integral = NaN
	if x <= 1.0e-10
		# Taylor series
		planck_integral = x * x - x * x * x / 2.0
	else
		planck_integral = std_pow(x, 3.0) / (std_exp(x) - 1.0) # C++ std::pow(x, 3)
	end
	# C++: coeff / (std::pow(PI, 4) / 15.0) * (radiation_constant_ * std::pow(T, 4)). Julia's pi^4
	# would be the exactly-rounded value, 1 ULP above what std::pow(PI, 4) returns.
	return coeff / (std_pow(Float64(pi), 4.0) / 15.0) * (radiation_constant_ * std_pow(T, 4.0)) * planck_integral
end

# The quokka defaults are zero; a problem overrides them to add line cooling and cosmic-ray heating.
DefineNetCoolingRate(temperature, num_density) = zeros(nGroups_)
DefineNetCoolingRateTempDerivative(temperature, num_density) = zeros(nGroups_)
DefineCosmicRayHeatingRate(num_density) = 0.0
