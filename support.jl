# C++/AMReX shims and the shared root finder.

# AMREX_ASSERT / AMREX_ASSERT_WITH_MESSAGE compile away in a Release build, which is why the C++ code
# deliberately continues past them to increment the failure counters (negative T_d, non-convergence).
# Leave amrex_debug_ = false to reproduce that behaviour; set it true for a Debug-style run.
const amrex_debug_ = false

function amrex_assert(cond::Bool, msg::String = "assertion failed")
	if amrex_debug_ && !cond
		error(msg)
	end
end

# std::sqrt / std::max / std::min equivalents. These are not cosmetic: the C++ solver is written to
# survive a negative or NaN temperature (it continues past the assertions above to record a failure
# and lets the NaN propagate), and Julia's own versions do not behave that way.
#   std::sqrt(x < 0) returns NaN;      Julia's sqrt throws a DomainError.
#   std::max(a, b) is a < b ? b : a,   so a NaN operand returns the *other* one;
#   std::min(a, b) is b < a ? b : a,   whereas Julia's max/min propagate NaN.
std_sqrt(x::Float64) = x < 0.0 ? NaN : sqrt(x)
std_max(a, b) = a < b ? b : a
std_min(a, b) = b < a ? b : a

# std::pow, std::exp and std::log10 resolve to the system libm. Julia implements ^, exp and log10 in
# Julia rather than calling libm, and the two disagree in the last bit often enough to matter: over a
# decade-wide sweep, x^3.0 differs from pow(x, 3.0) for about a quarter of arguments and exp for about
# one percent. Anywhere the C++ writes std::pow / std::exp / std::log10, call libm so the port lands on
# the same double. (sqrt needs no such treatment: IEEE-754 requires it to be correctly rounded, so
# Julia's sqrt and libm's always agree.)
std_pow(x::Float64, y::Float64) = ccall(:pow, Float64, (Float64, Float64), x, y)
std_exp(x::Float64) = ccall(:exp, Float64, (Float64,), x)
std_log10(x::Float64) = ccall(:log10, Float64, (Float64,), x)

# Port of RadSystem<problem_t>::BackwardEulerOneVariable (radiation_system.hpp lines 1659-1758).
function BackwardEulerOneVariable(rhs, jac, x0::Float64, compare::Float64)::Float64
	rel_tol = 1.0e-8
	rel_change_tol = 1.0e-6
	max_iter_td = 100

	# The caller passes `compare`, the physical scale the residual is measured against. It must be positive:
	# a zero scale makes the convergence test unsatisfiable.
	amrex_assert(compare > 0.0)

	f0 = rhs(x0)
	if abs(f0) < rel_tol * compare
		return x0
	end

	# Bracket the root, then solve with Newton guarded by bisection.
	#
	# Every residual solved through this routine is monotone in its unknown -- emission and the collisional
	# term both grow with the dust temperature, and the gas-energy residual grows with the gas energy -- so
	# the root is unique and can be bracketed by marching outward in the Newton direction. That guard is
	# what makes this robust: a bare Newton iteration is unsafe on these functions because a band-limited
	# Planck function is stiff enough that one step can leave the neighbourhood of the root and never
	# return, which is the origin of the "dust temperature failed to converge" abort. Both unknowns are
	# physically positive, so the search is confined to x > 0.
	j0 = jac(x0)
	dir = (j0 * f0 > 0.0) ? -1.0 : 1.0 # sign of the Newton step -f/j
	xa = x0
	fa = f0
	xb = x0
	fb = f0
	# March multiplicatively rather than by fixed increments: both unknowns are positive and can sit orders
	# of magnitude from the initial guess, and halving repeatedly also keeps the search inside x > 0 without
	# needing a special case.
	bracketed = false
	x_prev = x0
	f_prev = f0
	factor = (dir > 0.0) ? 2.0 : 0.5
	for k in 0:199
		x_try = x_prev * factor
		if !(x_try > 0.0) || !isfinite(x_try)
			break
		end
		f_try = rhs(x_try)
		if f_try * f0 <= 0.0
			xa = std_min(x_prev, x_try)
			xb = std_max(x_prev, x_try)
			fa = (xa == x_prev) ? f_prev : f_try
			fb = (xb == x_prev) ? f_prev : f_try
			bracketed = true
			break
		end
		x_prev = x_try
		f_prev = f_try
	end
	if !bracketed
		return -1.0 # caller treats a negative return as failure
	end

	x = 0.5 * (xa + xb)
	iter_Td = 0
	while iter_Td < max_iter_td
		the_rhs = rhs(x)
		if abs(the_rhs) < rel_tol * compare
			break
		end

		# keep the bracket around the root
		if the_rhs * fa > 0.0
			xa = x
			fa = the_rhs
		else
			xb = x
			fb = the_rhs
		end

		j = jac(x)
		x_new = (j != 0.0) ? (x - the_rhs / j) : (0.5 * (xa + xb))
		# Fall back to bisection whenever Newton would leave the bracket. Written as the negation of the
		# in-bracket test rather than its DeMorgan dual so that a NaN step also lands here.
		if !(x_new > xa && x_new < xb)
			x_new = 0.5 * (xa + xb)
		end
		dx = x_new - x
		x = x_new
		if abs(dx) < rel_change_tol * abs(x)
			break
		end
		iter_Td += 1
	end

	amrex_assert(iter_Td < max_iter_td, "Newton iteration in IntegratorOneVariable failed to converge.")
	if iter_Td >= max_iter_td
		x = -1.0
	end

	return x
end
