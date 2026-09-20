# Reproduce the "Newton-Raphson iteration for matter-radiation coupling failed to converge!" abort
# (radiation_dust_system.hpp:671) that DTypeFront3D hits on Coarse STEP 1 of DTypeFront3D_crash.toml.
#
# Nothing in this file is transcribed by hand. The crashing call is found in the C++ dump as the SGIN
# record with no matching SGOUT -- the solver aborted before returning, so no output record was written
# -- and its arguments are read straight from that line. The problem's traits, constants and EOS table
# come from the dump's SGTRAITS lines (see problem_DTypeFront3D.jl), and every quantity the solver
# derives from the inputs is checked against the C++ value the same record carries.
#
#   julia repro_DTypeFront3D.jl [dump.txt]

const SGDUMP_PATH = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "dump.txt")
include(joinpath(@__DIR__, "..", "DTypeFront3D", "problem_DTypeFront3D.jl"))
include("dump_io.jl")

function main()
	completed, aborted = read_dump(SGDUMP_PATH)
	println("dump: ", length(completed), " completed calls, ", length(aborted), " aborted")
	if isempty(aborted)
		println("no aborted call in this dump -- nothing to reproduce")
		return 0
	end
	length(aborted) == 1 || println("note: more than one aborted call; using the last")
	r = aborted[end]

	println("crashing call: index ", r.idx, ", cell ", r.ijk, ", outer iter ", r.iter)
	println()

	# Check every quantity the C++ derived from these inputs before trusting the replay. A mismatch here
	# means the port disagrees about the setup, so a matching failure downstream would prove nothing.
	println("derived quantities vs the C++ values in the same record (exact equality required):")
	d = derived(r)
	nbad = 0
	for k in DERIVED_FIELDS
		haskey(r.vecs, k) || continue
		a = d[k]; b = r.vecs[k]
		for q in eachindex(b)
			u = ulps(a[q], b[q])
			u != 0 && (nbad += 1)
			lbl = length(b) == 1 ? k : "$k[$q]"
			println("   ", rpad(lbl, 22), rpad(a[q], 26), rpad(b[q], 26), u == 0 ? "MATCH" : "DIFFER ($u ULP)")
		end
	end
	println()
	nbad == 0 || println("WARNING: $nbad derived values differ; the reproduction below is not trustworthy\n")

	res, cnt, fail = call_solver(r)
	n = cnt[2] - 1
	println("Newton iterations n = $n  (maxIter = 100)")
	println("failure counters [non-convergence, negative T_d, outer] = $fail")
	println("decoupled dust branch (dust_model == 2) taken: $(cnt[4] == 1)")
	println("T_gas = $(res.T_gas) K,  T_d = $(res.T_d) K")
	println("Egas  = $(res.Egas)")
	println("Erad  = $(res.EradVec)")
	println()
	if fail[1] > 0
		println("REPRODUCED: the Newton-Raphson iteration hit maxIter without converging,")
		println("which is the C++ abort at radiation_dust_system.hpp:671")
	else
		println("NOT reproduced: the Julia solve converged in $n iterations")
	end
	return fail[1] > 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
	exit(main())
end
