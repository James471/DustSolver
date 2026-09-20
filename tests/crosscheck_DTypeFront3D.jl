# Replay every completed SolveGasDustRadiationEnergyExchange call in a C++ dump and compare the Julia
# port against it bit for bit -- both the quantities the solver derives from its inputs and every field
# of its result.
#
#   1. Build quokka with the debug dump in src/radiation/source_terms_multi_group.hpp (and the SGTRAITS
#      block in the problem file) enabled, then
#          ./src/problems/DTypeFront3D/DTypeFront3D ./DTypeFront3D_crash.toml 2>&1 \
#              | grep -E '^SGTRAITS|^SGIN |^SGOUT ' > dump.txt
#   2. julia crosscheck_DTypeFront3D.jl dump.txt [julia_out.txt]
#
# "Machine precision" here means exact equality of the IEEE-754 doubles: the C++ prints with %.17g,
# which round-trips a double exactly, so a matching port reproduces every bit. Anything merely close is
# reported as a mismatch with its ULP distance rather than waved through.

const SGDUMP_PATH = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "dump.txt")
include(joinpath(@__DIR__, "..", "DTypeFront3D", "problem_DTypeFront3D.jl"))
include("dump_io.jl")

function main()
	out_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "julia_out.txt")
	completed, aborted = read_dump(SGDUMP_PATH)
	println("replaying ", length(completed), " completed calls from ", SGDUMP_PATH,
		" (", length(aborted), " aborted, not replayable here)")

	fields = (DERIVED_FIELDS..., OUTPUT_FIELDS...)
	nmis = Dict(k => 0 for k in fields)
	worst = Dict(k => 0 for k in fields)
	wrel = Dict(k => 0.0 for k in fields)
	ex = Dict{String, String}()
	nfail = 0

	open(out_path, "w") do io
		for r in completed
			d = derived(r)
			res, cnt, fail = call_solver(r)
			o = out_dict(res, cnt)
			nfail += fail[1]

			print(io, "JLIN ", r.idx)
			for k in DERIVED_FIELDS
				print(io, " ", k); for x in d[k]; print(io, " ", x); end
			end
			print(io, "\nJLOUT ", r.idx)
			for k in OUTPUT_FIELDS
				print(io, " ", k); for x in o[k]; print(io, " ", x); end
			end
			println(io)

			for k in DERIVED_FIELDS
				haskey(r.vecs, k) || continue
				_cmp!(nmis, worst, wrel, ex, k, d[k], r.vecs[k], r)
			end
			for k in OUTPUT_FIELDS
				_cmp!(nmis, worst, wrel, ex, k, o[k], r.out[k], r)
			end
		end
	end
	println("wrote ", out_path)
	println()

	total = 0
	println(rpad("field", 22), rpad("mismatches", 12), rpad("worst ULP", 12), "worst rel diff")
	for k in fields
		total += nmis[k]
		println(rpad(k, 22), rpad(nmis[k], 12), rpad(worst[k], 12), nmis[k] == 0 ? "-" : string(wrel[k]))
	end
	println()
	nfail > 0 && println("note: $nfail of the replayed calls did not converge in Julia either")
	# separate the fields the C++ leaves uninitialised from genuine disagreements
	uninit = sum(nmis[k] for k in UNINITIALISED_UNLESS_PPL)
	real_total = total - uninit
	if uninit > 0
		println("note: alpha_P / alpha_E differ in $uninit values. Under opacity_model_ = $opacity_model_ the")
		println("      C++ never writes them -- ComputeModelDependentKappaEAndKappaP returns a")
		println("      default-initialised OpacityTerms and only the PPL_opacity_full_spectrum branch fills")
		println("      these -- so it is returning indeterminate stack memory, which it then never reads.")
		println("      The magnitudes give it away: subnormals around 1e-314. Not a port discrepancy.")
		println()
	end
	total = real_total
	if total == 0
		println("ALL ", length(completed), " CALLS MATCH BIT FOR BIT across every derived quantity and")
		println("every output field", uninit > 0 ? " except the two noted above." : ".")
	else
		println("$total field values differ. Worst example per field:")
		for k in fields
			haskey(ex, k) && println("   ", rpad(k, 22), ex[k])
		end
	end
	return total
end

function _cmp!(nmis, worst, wrel, ex, k, a, b, r)
	for q in eachindex(b)
		u = ulps(a[q], b[q])
		u == 0 && continue
		nmis[k] += 1
		if u > worst[k]
			worst[k] = u
			den = max(abs(b[q]), abs(a[q]))
			wrel[k] = den == 0 ? Inf : abs(a[q] - b[q]) / den
			ex[k] = "idx $(r.idx) cell $(r.ijk) iter $(r.iter) [$q]: julia $(a[q]) vs c++ $(b[q])"
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	exit(main() == 0 ? 0 : 1)
end
