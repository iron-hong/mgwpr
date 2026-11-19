#' Simulate a spatial panel with known multiscale coefficient surfaces
#'
#' Generates a regular spatial grid, known location-specific coefficients, a
#' time-invariant spatial context effect, four spatially confounded predictors,
#' and a continuous outcome in long panel format.
#'
#' @param n_grid Number of rows and columns in the square spatial grid.
#' @param n_time Number of periods per spatial unit.
#' @param sigma_x Standard deviation of the independent component of each
#'   predictor.
#' @param spatial_coupling Strength of the time-invariant context effect in
#'   each predictor.
#' @param sigma_epsilon Standard deviation of the outcome disturbance.
#' @param context_scale Multiplier controlling the scale of the true spatial
#'   context effect.
#' @param time_values Optional vector of distinct period labels. Defaults to
#'   `0, ..., n_time - 1`.
#' @param seed Optional non-negative integer random seed. When supplied, the
#'   caller's existing random-number state is restored after simulation.
#' @return A list of class `mgwpr_simulation` with `units` (one row per
#'   location, including true coefficient surfaces), `panel` (long-format
#'   observations), `diagnostics` (two generated-data correlations), and
#'   `settings`.
#' @examples
#' sim <- simulate_mgwpr_panel(seed = 42)
#' head(sim$panel)
#' head(sim$units)
#' @export
simulate_mgwpr_panel <- function(n_grid = 15L, n_time = 3L,
                                 sigma_x = 0.5,
                                 spatial_coupling = 0.05,
                                 sigma_epsilon = 0.5,
                                 context_scale = 30,
                                 time_values = NULL,
                                 seed = NULL) {
  validate_integer <- function(x, name, minimum) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x != as.integer(x) || x < minimum) {
      stop(name, " must be an integer greater than or equal to ", minimum, ".")
    }
    as.integer(x)
  }
  n_grid <- validate_integer(n_grid, "n_grid", 2L)
  n_time <- validate_integer(n_time, "n_time", 2L)
  positive_settings <- list(
    sigma_x = sigma_x,
    sigma_epsilon = sigma_epsilon,
    context_scale = context_scale
  )
  for (nm in names(positive_settings)) {
    value <- positive_settings[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(nm, " must be one positive finite number.")
    }
  }
  if (!is.numeric(spatial_coupling) || length(spatial_coupling) != 1L ||
      !is.finite(spatial_coupling) || spatial_coupling < 0) {
    stop("spatial_coupling must be one non-negative finite number.")
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) ||
                         !is.finite(seed) || seed < 0 ||
                         seed > .Machine$integer.max ||
                         seed != as.integer(seed))) {
    stop("seed must be NULL or one non-negative integer.")
  }

  if (is.null(time_values)) {
    time_values <- 0:(n_time - 1L)
  } else if (length(time_values) != n_time || anyNA(time_values) ||
             anyDuplicated(time_values)) {
    stop("time_values must contain n_time distinct, non-missing values.")
  }

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  n_units <- n_grid * n_grid
  n_obs <- n_units * n_time
  grid_i <- rep(seq_len(n_grid), each = n_grid)
  grid_j <- rep(seq_len(n_grid), times = n_grid)
  units <- data.frame(
    unit_id = seq_len(n_units),
    grid_i = grid_i,
    grid_j = grid_j
  )

  q <- ceiling(n_grid / 4)
  units$beta1_true <- 1 +
    ((q^2 - (q - units$grid_i / 2)^2) *
       (q^2 - (q - units$grid_j / 2)^2)) / q^4
  units$beta2_true <- 1 +
    (units$grid_i + units$grid_j) / (2 * n_grid)
  units$beta3_true <- 1.5
  units$beta4_true <- 0
  units$alpha_true <- context_scale *
    (exp(units$grid_j / n_grid) - 1)

  unit_index <- rep(seq_len(n_units), each = n_time)
  panel <- data.frame(
    unit_id = units$unit_id[unit_index],
    time_id = rep(time_values, times = n_units)
  )
  panel <- cbind(panel, units[unit_index, setdiff(names(units), "unit_id"), drop = FALSE])
  rownames(panel) <- NULL

  spatial_context <- panel$alpha_true
  panel$x1 <- sigma_x * stats::rnorm(n_obs) + spatial_coupling * spatial_context
  panel$x2 <- sigma_x * stats::rnorm(n_obs) + spatial_coupling * spatial_context
  panel$x3 <- sigma_x * stats::rnorm(n_obs) + spatial_coupling * spatial_context
  panel$x4 <- sigma_x * stats::rnorm(n_obs) + spatial_coupling * spatial_context
  panel$epsilon <- sigma_epsilon * stats::rnorm(n_obs)
  panel$y <- panel$alpha_true +
    panel$beta1_true * panel$x1 +
    panel$beta2_true * panel$x2 +
    panel$beta3_true * panel$x3 +
    panel$epsilon

  diagnostics <- c(
    cor_x1_spatial_context = stats::cor(panel$x1, spatial_context),
    cor_x4_y = stats::cor(panel$x4, panel$y)
  )
  settings <- list(
    n_grid = n_grid, n_time = n_time, n_units = n_units, n_obs = n_obs,
    sigma_x = sigma_x, spatial_coupling = spatial_coupling,
    sigma_epsilon = sigma_epsilon, context_scale = context_scale,
    time_values = time_values, seed = seed
  )
  structure(list(units = units, panel = panel,
                 diagnostics = diagnostics, settings = settings),
            class = "mgwpr_simulation")
}
