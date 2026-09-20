# Reading the C++ dump. The instrumented binary writes, per call to
# SolveGasDustRadiationEnergyExchange:
#
#   SGIN  <idx> <every argument> <every quantity the solver derives from them>
#   SGOUT <idx> <every field of the returned NewtonIterationResult, plus the iteration count n>
#
# SGIN is written before the call and SGOUT after, so a call that aborts inside the solver leaves an
# SGIN with no matching SGOUT -- that unmatched record is the crashing call, and it is the only way to
# reproduce it without copying numbers by hand.

struct DumpRecord
	idx::Int
	ijk::NTuple{3, Int}
	iter::Int
	scal::Dict{String, Float64}   # Egas0, rho, coeff_n, dt, Q_dust, tol, tol_rel, tempFloor, Tgas0, ...
	vecs::Dict{String, Vector{Float64}} # Erad0Vec, work, vel_times_F, Src, radBoundaries, massScalars, ...
	out::Union{Nothing, Dict{String, Vector{Float64}}}
end

# how many numbers follow each key; anything not listed is a scalar
const _VEC_ARITY = Dict(
	"Erad0Vec" => :ng, "work" => :ng, "vel_times_F" => :ng, "Src" => :ng,
	"radBoundaries" => :ngp1, "massScalars" => :nspec,
	"fourPiBoverC_Td0" => :ng, "dfourPiBoverCdT_Td0" => :ng, "kappaP_Td0" => :ng,
	"kappaE_Td0" => :ng, "kappaPoverE_Td0" => :ng, "tau0_Td0" => :ng,
	"EradVec" => :ng, "workOut" => :ng, "kappaP" => :ng, "kappaE" => :ng,
	"kappaF" => :ng, "kappaPoverE" => :ng, "delta_nu_kappa_B_at_edge" => :ng,
	"alpha_P" => :ng, "alpha_E" => :ng,
)

# Generic "key value..." walk: a token that does not parse as a number starts a new key.
function _kv(tokens)
	scal = Dict{String, Float64}()
	vecs = Dict{String, Vector{Float64}}()
	key = ""
	for tok in tokens
		v = tryparse(Float64, tok)
		if v === nothing
			key = tok
			vecs[key] = Float64[]
		else
			push!(vecs[key], v)
		end
	end
	for (k, v) in vecs
		length(v) == 1 && !haskey(_VEC_ARITY, k) && (scal[k] = v[1])
	end
	return scal, vecs
end

function read_dump(path)
	ins = Dict{Int, DumpRecord}()
	outs = Dict{Int, Dict{String, Vector{Float64}}}()
	order = Int[]
	for line in eachline(path)
		if startswith(line, "SGIN ")
			t = split(line)
			idx = parse(Int, t[2])
			@assert t[3] == "ijk"
			ijk = (parse(Int, t[4]), parse(Int, t[5]), parse(Int, t[6]))
			@assert t[7] == "iter"
			iter = parse(Int, t[8])
			scal, vecs = _kv(t[9:end])
			ins[idx] = DumpRecord(idx, ijk, iter, scal, vecs, nothing)
			push!(order, idx)
		elseif startswith(line, "SGOUT ")
			t = split(line)
			idx = parse(Int, t[2])
			_, vecs = _kv(t[3:end])
			outs[idx] = vecs
		end
	end
	isempty(ins) && error("no SGIN records in $path")
	completed = [DumpRecord(r.idx, r.ijk, r.iter, r.scal, r.vecs, outs[r.idx]) for r in (ins[i] for i in order) if haskey(outs, r.idx)]
	aborted = [ins[i] for i in order if !haskey(outs, i)]
	return completed, aborted
end

# ULP distance, so a near-miss is quantified rather than merely flagged
function ulps(a::Float64, b::Float64)
	(isnan(a) && isnan(b)) && return 0
	a == b && return 0
	(isnan(a) || isnan(b) || isinf(a) || isinf(b)) && return typemax(Int)
	ia = reinterpret(Int64, a); ib = reinterpret(Int64, b)
	ia = ia < 0 ? typemin(Int64) - ia : ia
	ib = ib < 0 ? typemin(Int64) - ib : ib
	return abs(ia - ib)
end

# call the ported solver with a record's inputs, exactly as quokka called it
function call_solver(r::DumpRecord)
	cnt = zeros(Int, 4); fail = zeros(Int, 3)
	res = SolveGasDustRadiationEnergyExchange(
		r.scal["Egas0"], r.vecs["Erad0Vec"], r.scal["rho"], r.scal["coeff_n"], r.scal["dt"],
		r.vecs["massScalars"], r.iter, r.vecs["work"], r.vecs["vel_times_F"], r.vecs["Src"],
		r.scal["Q_dust"], r.vecs["radBoundaries"], r.scal["tol"], r.scal["tol_rel"],
		r.scal["tempFloor"], cnt, fail)
	return res, cnt, fail
end

# recompute, in Julia, every derived quantity the C++ dumped in SGIN
function derived(r::DumpRecord)
	rho = r.scal["rho"]; ms = r.vecs["massScalars"]; rb = r.vecs["radBoundaries"]
	Tg0 = ComputeTgasFromEint(rho, r.scal["Egas0"], ms)
	Td0 = ComputeDustTemperatureBateKeto(Tg0, Tg0, rho, r.vecs["Erad0Vec"], r.scal["coeff_n"],
					     r.scal["dt"], NaN, 0, r.scal["Q_dust"], rb)
	fpb = ComputeThermalRadiationMultiGroup(Td0, rb)
	dfpb = ComputeThermalRadiationTempDerivativeMultiGroup(Td0, rb)
	op = ComputeModelDependentKappaEAndKappaP(Td0, rho, rb, zeros(nGroups_), fpb, r.vecs["Erad0Vec"], 0)
	return Dict{String, Vector{Float64}}(
		"Tgas0" => [Tg0], "Td0" => [Td0],
		"c_v" => [ComputeEintTempDerivative(rho, Tg0, ms)],
		"H_num_den" => [ComputeNumberDensityH(rho, ms)],
		"fourPiBoverC_Td0" => fpb, "dfourPiBoverCdT_Td0" => dfpb,
		"kappaP_Td0" => op.kappaP, "kappaE_Td0" => op.kappaE, "kappaPoverE_Td0" => op.kappaPoverE,
		"tau0_Td0" => r.scal["dt"] * rho * op.kappaP * c_hat_,
	)
end

const DERIVED_FIELDS = ("Tgas0", "Td0", "c_v", "H_num_den", "fourPiBoverC_Td0", "dfourPiBoverCdT_Td0",
			"kappaP_Td0", "kappaE_Td0", "kappaPoverE_Td0", "tau0_Td0")
# Every field of the returned NewtonIterationResult, including all seven members of its OpacityTerms.
const OUTPUT_FIELDS = ("n", "Egas", "T_gas", "T_d", "EradVec", "workOut", "kappaP", "kappaE", "kappaF",
		       "kappaPoverE", "delta_nu_kappa_B_at_edge", "alpha_P", "alpha_E")

# alpha_P and alpha_E are only ever written by the PPL_opacity_full_spectrum branch of
# ComputeModelDependentKappaEAndKappaP. Under any other opacity model that function returns a
# default-initialised OpacityTerms, whose amrex::GpuArray members are left indeterminate, so the C++
# hands back whatever was on the stack -- and never reads it again. A port cannot and should not
# reproduce that, so these two are compared but reported apart from the rest.
const UNINITIALISED_UNLESS_PPL = ("alpha_P", "alpha_E")

function out_dict(res, cnt)
	return Dict{String, Vector{Float64}}(
		"n" => [Float64(cnt[2] - 1)], "Egas" => [res.Egas], "T_gas" => [res.T_gas], "T_d" => [res.T_d],
		"EradVec" => res.EradVec, "workOut" => res.work,
		"kappaP" => res.opacity_terms.kappaP, "kappaE" => res.opacity_terms.kappaE,
		"kappaF" => res.opacity_terms.kappaF, "kappaPoverE" => res.opacity_terms.kappaPoverE,
		"delta_nu_kappa_B_at_edge" => res.opacity_terms.delta_nu_kappa_B_at_edge,
		"alpha_P" => res.opacity_terms.alpha_P, "alpha_E" => res.opacity_terms.alpha_E)
end
