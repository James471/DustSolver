# The DTypeFront3D problem, defined entirely from the C++ dump rather than transcribed by hand.
#
# Every trait, physical constant and EOS table entry below is read from the SGTRAITS lines the
# instrumented binary writes, so there is nothing here for me to get wrong by copying: if the problem,
# the network or a unit system changes, the next dump changes with it. The only things written out as
# code are the *shapes* of quokka functions (the EOS relation, the opacity layout), and each of those
# is checked against a dumped value by the drivers.
#
# The driver must set SGDUMP_PATH before including this file.

isdefined(Main, :SGDUMP_PATH) ||
	error("define SGDUMP_PATH (path to the C++ SGIN/SGOUT/SGTRAITS dump) before including problem_DTypeFront3D.jl")

include(joinpath(@__DIR__, "..", "types.jl")) # defines OpacityModel

# Parse "SGTRAITS key v [v ...] key v ..." into key => Vector{Float64}
function _read_traits(path)
	d = Dict{String, Vector{Float64}}()
	for line in eachline(path)
		startswith(line, "SGTRAITS") || continue
		t = split(line)[2:end]
		key = ""
		for tok in t
			v = tryparse(Float64, tok)
			if v === nothing
				key = tok
				d[key] = Float64[]
			else
				push!(d[key], v)
			end
		end
	end
	isempty(d) && error("no SGTRAITS lines in $path")
	return d
end

const TRAITS = _read_traits(SGDUMP_PATH)
_t(k) = (haskey(TRAITS, k) || error("SGTRAITS has no key $k"); TRAITS[k][1])
_tv(k) = (haskey(TRAITS, k) || error("SGTRAITS has no key $k"); TRAITS[k])

# Physics_Traits / RadSystem_Traits / ISM_Traits / EOS_Traits, straight from the dump
const nGroups_ = Int(_t("nGroups"))
const nGroupsThermal_ = Int(_t("nGroupsThermal"))
const beta_order_ = Int(_t("beta_order"))
const opacity_model_ = OpacityModel(Int(_t("opacity_model")))
const gamma_ = _t("gamma")
const c_light_ = _t("c_light")
const c_hat_ = _t("c_hat")
const radiation_constant_ = _t("a_rad")
const boltzmann_constant_ = _t("k_B")
const energy_unit_ = _t("energy_unit")
const Erad_floor_ = _t("Erad_floor")
const mean_molecular_mass_ = _t("mean_molecular_mass")
const gas_dust_coupling_threshold = _t("gas_dust_coupling_threshold")
const radBoundaries_ = _tv("radBoundaries")

# the problem's two gray dust opacities (photoionize.kappa_ir / kappa_optical)
const kappa_ir_ = _t("kappa_ir")
const kappa_optical_ = _t("kappa_optical")

# EOS species table, filled by actual_eos_init from the network's runtime parameters
const spmasses_ = _tv("spmasses")
const gammas_ = _tv("eos_gammas")

include(joinpath(@__DIR__, "..", "hyperparameters.jl"))
include(joinpath(@__DIR__, "..", "support.jl"))
include(joinpath(@__DIR__, "..", "planck_integral.jl"))
include(joinpath(@__DIR__, "..", "radiation_defaults.jl")) # Planck fractions, emission and its T-derivative, zero cooling rates

# --- DTypeFront3D's own definitions. problem.jl is deliberately NOT included: its constant-kappa
# --- opacity and ideal-gas EOS are not this problem's. ------------------------------------------------

# RadSystem<DTypeFront3D>::DefineOpacityExponentsAndLowerValues. Each thermal group carries its own
# constant gray opacity; the ionizing (chemistry) band is transparent to it. Checked against the dumped
# kappaP_Td0.
function DefineOpacityExponentsAndLowerValues(rad_boundaries, rho, Tgas)
	kappa_g = vcat([kappa_ir_, kappa_optical_], zeros(nGroups_ - 2))
	exponents_and_values = [zeros(nGroups_ + 1), zeros(nGroups_ + 1)]
	for i in 1:(nGroups_ + 1)
		exponents_and_values[1][i] = 0.0
		exponents_and_values[2][i] = (i <= nGroups_) ? kappa_g[i] : 0.0
	end
	return exponents_and_values
end

# Port of quokka::EOSMicrophysics (hydro/EOS.hpp), which is what DTypeFront3D resolves to: its
# EOS_Traits names no EOSBackend and the build defines PHOTOCHEMISTRY, so DefaultEOSBackend picks
# EOSMicrophysics over EOSIdeal. It forwards to the "multigamma" photoionization EOS in
# extern/Microphysics/EOS/photoionization/actual_eos.H, where each species is an ideal gas with its own
# gamma:
#
#     xn[n]            = massScalars[n] / spmasses[n]          (a number density, not a mass fraction)
#     rhotot           = sum_n xn[n] * spmasses[n]
#     sum_ni_fi_over_2 = sum_n xn[n] / (gammas[n] - 1)
#     eos_input_re:  T = (Eint / rho) * rhotot / (sum_ni_fi_over_2 * k_B)
#     eos_input_rt:  e = sum_ni_fi_over_2 * k_B * T / rhotot,  dedT = sum_ni_fi_over_2 * k_B / rhotot
#
# Checked against the dumped Tgas0 and c_v.
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
