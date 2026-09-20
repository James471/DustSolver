# Hyper parameters for the radiation solver.
# Line-for-line mirror of src/radiation/radiation_system.hpp lines 37-70, comments trimmed to the
# one-liners; the long block comments there explain the reasoning and are worth reading alongside.

const add_line_cooling_to_radiation_in_jac = false
const include_delta_B = true
const use_diffuse_flux_mean_opacity = true
const special_edge_bin_slopes = false	   # Use 2 and -4 as the slopes for the first and last bins, respectively
const force_rad_floor_in_iteration = false # force radiation energy density to be positive (and above the floor value) in the Newton iteration
const include_work_term_in_source = true
const max_iter_to_update_alpha_E = 5	   # Apply to the PPL_opacity_full_spectrum only
const enable_dE_constrain = false
# Multiple of the estimated double-precision round-off floor at which the radiation residual is
# accepted as converged.
const newton_resid_roundoff_factor = 10.0
# Adaptive damping for the decoupled-dust branch, where the iteration can oscillate with growing
# amplitude instead of converging.
const newton_damping_down = 0.5
const newton_damping_up = 1.5
const newton_damping_min = 0.05
const use_D_as_base = false
# Optical depth below which a group's Newton unknown is its radiation energy Erad_g rather than its
# exchange term R_g.
const newton_erad_base_tau_threshold = 1.0
const PPL_free_slope_st_total = false
