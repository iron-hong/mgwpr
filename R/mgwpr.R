#' Multiscale geographically weighted panel regression
#'
#' Fits pooled or individual fixed-effects geographically weighted panel
#' regression by iterative backfitting. Gaussian responses support pooled and
#' within fixed-effects fits. Binomial responses support pooled logistic fits.
#' A separate spatial bandwidth is selected for every model-matrix term.
#'
#' @param formula A model formula.
#' @param data A data frame in long panel format.
#' @param coords Names of two numeric coordinate columns in data. Supply either
#'   `coords` or `spatial`.
#' @param spatial Optional `sf` object with one geometry per panel unit and an
#'   identifier column named by `id`. Polygon centroids/representative points
#'   are used to calculate distances; original geometries are retained for maps.
#' @param id Name of the panel-unit identifier. Required for method="fe";
#'   may also be supplied for pooled models to identify locations.
#' @param time Optional name of the time variable, used for validation only.
#' @param method "fe" for individual fixed effects (MGWFER) or "pooled".
#' @param family "gaussian" for continuous outcomes or "binomial" for binary
#'   outcomes coded 0/1 or as a two-level factor. Binomial fits currently use
#'   `method = "pooled"`.
#' @param bw Optional numeric vector of bandwidths, one per model-matrix
#'   coefficient. For adaptive bandwidths, values are neighbor counts.
#' @param adaptive If TRUE, search/select numbers of nearest locations;
#'   otherwise use coordinate-distance bandwidths.
#' @param kernel Spatial kernel: "bisquare", "gaussian", "exponential",
#'   "tricube", or "boxcar".
#' @param selector Bandwidth criterion; currently "cv" (leave-one-location-
#'   out cross-validation).
#' @param bw_min,bw_max Optional lower and upper search bounds. Adaptive bounds
#'   are counts of locations; fixed bounds are distances.
#' @param max_iter Maximum backfitting iterations.
#' @param tol Convergence tolerance on the largest coefficient change.
#' @param scale If true, center and scale model-matrix columns and the outcome
#'   before fitting.
#' @param verbose Print iteration progress.
#' @return An object of class mgwpr with location-level coefficients,
#'   observation-level fitted values and residuals, selected bandwidths, and
#'   convergence information. Fixed-effects fits also include recovered
#'   location-specific context effects.
#' @export
mgwpr <- function(formula, data, coords = NULL, spatial = NULL,
                  id = NULL, time = NULL,
                  method = c("fe", "pooled"),
                  family = c("gaussian", "binomial"), bw = NULL,
                  adaptive = TRUE,
                  kernel = c("bisquare", "gaussian", "exponential",
                             "tricube", "boxcar"),
                  selector = "cv", bw_min = NULL, bw_max = NULL,
                  max_iter = 30L, tol = 1e-5, scale = FALSE,
                  verbose = FALSE) {
  method <- match.arg(method)
  family <- match.arg(family)
  kernel <- match.arg(kernel)
  if (!identical(selector, "cv")) stop("Only selector='cv' is supported.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (is.null(coords) == is.null(spatial)) stop("Supply exactly one of coords or spatial.")
  if (!is.null(coords)) {
    if (length(coords) != 2L || !all(coords %in% names(data))) stop("coords must name two columns in data.")
    if (!all(vapply(data[coords], is.numeric, logical(1)))) stop("Coordinate columns must be numeric.")
  }
  if (!is.null(spatial) && !inherits(spatial, "sf")) stop("spatial must be an sf object.")
  id_name <- if (is.character(id) && length(id) == 1L && id %in% names(data)) id else NULL
  time_name <- if (is.character(time) && length(time) == 1L && time %in% names(data)) time else NULL
  if (is.character(id) && length(id) == 1L && id %in% names(data)) id <- data[[id]]
  if (is.character(time) && length(time) == 1L && time %in% names(data)) time <- data[[time]]
  if (is.character(id) && length(id) == 1L) stop("id does not match a column in data.")
  if (is.character(time) && length(time) == 1L) stop("time does not match a column in data.")
  if (method == "fe" && is.null(id)) stop("id is required when method='fe'.")
  if (!is.null(spatial) && is.null(id_name)) {
    stop("When spatial is supplied, id must be the name of an ID column in data and spatial.")
  }
  if (!is.null(spatial) && is.null(id)) {
    stop("id is required when spatial is supplied.")
  }
  if (family == "binomial" && method != "pooled") {
    stop("family='binomial' currently requires method='pooled'.")
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 || max_iter != as.integer(max_iter)) stop("max_iter must be a positive integer.")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive.")

  mf <- model.frame(formula, data = data, na.action = na.fail)
  response <- model.response(mf)
  if (family == "binomial") {
    if (is.factor(response) || is.character(response)) {
      response <- factor(response)
      if (nlevels(response) != 2L) stop("A binomial response factor must have exactly two observed levels.")
      response_levels <- levels(response)
      y_raw <- as.numeric(response) - 1
    } else {
      if (!is.numeric(response) || !all(na.omit(unique(response)) %in% c(0, 1))) {
        stop("A binomial response must be coded 0/1 or be a two-level factor.")
      }
      response_levels <- c("0", "1")
      y_raw <- as.numeric(response)
    }
    if (length(unique(y_raw)) != 2L) stop("Both outcome classes must be present.")
  } else {
    if (!is.numeric(response)) stop("A gaussian response must be numeric.")
    y_raw <- as.numeric(response)
    response_levels <- NULL
  }
  X0 <- model.matrix(attr(mf, "terms"), mf)
  if (method == "fe") X0 <- X0[, colnames(X0) != "(Intercept)", drop = FALSE]
  if (ncol(X0) == 0L) stop("No slope terms remain after model-matrix construction.")
  Xraw <- X0
  xy_obs <- if (!is.null(coords)) as.matrix(data[coords]) else NULL
  if (nrow(mf) != nrow(data)) stop("Missing model variables are not supported; remove incomplete rows first.")
  if ((!is.null(xy_obs) && any(!is.finite(xy_obs))) ||
      any(!is.finite(y_raw)) || any(!is.finite(X0))) stop("Model variables and coordinates must be complete and finite.")

  if (is.null(id)) {
    loc_key <- paste(format(xy_obs[, 1], digits = 15), format(xy_obs[, 2], digits = 15), sep = "\r")
    unit <- match(loc_key, unique(loc_key))
    id_values <- unique(loc_key)
  } else {
    if (length(id) != nrow(data) || anyNA(id)) stop("id must identify every row.")
    unit <- match(id, unique(id))
    id_values <- unique(as.character(id))
  }
  n <- nrow(data); m <- max(unit)
  if (!is.null(spatial)) {
    if (is.null(id_name) || !id_name %in% names(spatial)) {
      stop("When spatial is supplied, id must name its unit-ID column.")
    }
    spatial_ids <- as.character(spatial[[id_name]])
    if (anyNA(spatial_ids) || anyDuplicated(spatial_ids)) {
      stop("spatial must contain one non-missing geometry per unit ID.")
    }
    spatial_index <- match(id_values, spatial_ids)
    if (anyNA(spatial_index)) stop("spatial is missing one or more panel units.")
    spatial_aligned <- spatial[spatial_index, , drop = FALSE]
    spatial_aligned$location_id <- id_values
    representative <- sf::st_point_on_surface(sf::st_geometry(spatial_aligned))
    xy <- sf::st_coordinates(representative)[, 1:2, drop = FALSE]
    if (nrow(xy) != m) stop("Each spatial geometry must yield one representative point.")
    D <- as.matrix(sf::st_distance(representative))
    D <- matrix(as.numeric(D), nrow = m, ncol = m)
    spatial_aligned$coord_x <- xy[, 1]
    spatial_aligned$coord_y <- xy[, 2]
  } else {
    xy <- matrix(NA_real_, m, 2L)
    for (g in seq_len(m)) {
      rows <- which(unit == g)
      if (any(apply(xy_obs[rows, , drop = FALSE], 2, function(z) diff(range(z)) > 1e-10)))
        stop("Coordinates must be constant within each panel unit.")
      xy[g, ] <- xy_obs[rows[1L], ]
    }
    D <- as.matrix(stats::dist(xy))
    spatial_aligned <- sf::st_as_sf(
      data.frame(location_id = id_values, coord_x = xy[, 1], coord_y = xy[, 2]),
      coords = c("coord_x", "coord_y"), remove = FALSE, crs = NA
    )
  }
  if (!is.null(time)) {
    if (length(time) != n) stop("time must identify each observation.")
    if (!is.null(id) && anyDuplicated(data.frame(unit = unit, time = time))) stop("Each id-time combination must be unique.")
  }

  if (method == "fe") {
    y <- y_raw - ave(y_raw, unit, FUN = mean)
    X0 <- apply(X0, 2, function(z) z - ave(z, unit, FUN = mean))
    X0 <- as.matrix(X0)
    keep <- apply(X0, 2, function(z) stats::sd(z) > sqrt(.Machine$double.eps))
    if (!all(keep)) warning("Time-invariant or zero-variance terms were removed by the within transformation.")
    X0 <- X0[, keep, drop = FALSE]
    Xraw <- Xraw[, keep, drop = FALSE]
    if (!ncol(X0)) stop("No time-varying regressors remain after within transformation.")
  } else y <- y_raw

  # Standardization is applied after any within transformation.
  y_center <- if (scale && family == "gaussian") mean(y) else 0
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd <= 0) stop("Outcome has zero variance after any within transformation.")
  y_scale <- if (scale && family == "gaussian") y_sd else 1
  x_center <- if (scale) colMeans(X0) else rep(0, ncol(X0))
  x_scale <- if (scale) apply(X0, 2, stats::sd) else rep(1, ncol(X0))
  fixed_scale <- !is.finite(x_scale) | x_scale <= 0
  x_scale[fixed_scale] <- 1
  x_center[fixed_scale] <- 0
  yy <- if (family == "gaussian") (y - y_center) / y_scale else y
  XX <- sweep(sweep(X0, 2, x_center, "-"), 2, x_scale, "/")

  if (m < 4L) stop("At least four distinct panel locations are required.")
  if (is.null(bw_min)) bw_min <- if (adaptive) min(max(5L, ncol(XX) + 2L), m - 1L) else as.numeric(stats::quantile(D[D > 0], .15))
  if (is.null(bw_max)) bw_max <- if (adaptive) m else max(D)
  if (adaptive) {
    bw_min <- max(2L, as.integer(bw_min)); bw_max <- min(m, as.integer(bw_max))
    candidates <- unique(round(seq(bw_min, bw_max, length.out = min(12L, bw_max - bw_min + 1L))))
  } else {
    if (bw_min <= 0 || bw_max <= bw_min) stop("For fixed bandwidths require 0 < bw_min < bw_max.")
    candidates <- unique(as.numeric(stats::quantile(D[D > 0], probs = seq(.15, 1, length.out = 12), names = FALSE)))
    candidates <- candidates[candidates >= bw_min & candidates <= bw_max]
    candidates <- unique(c(bw_min, candidates, bw_max))
  }
  auto_bw <- is.null(bw)
  if (auto_bw) bw <- rep(NA_real_, ncol(XX)) else {
    if (length(bw) == 1L) bw <- rep(bw, ncol(XX))
    if (length(bw) != ncol(XX) || any(!is.finite(bw)) || any(bw <= 0)) stop("bw must be positive and have length 1 or one value per coefficient.")
  }

  B <- matrix(0, m, ncol(XX), dimnames = list(NULL, colnames(XX)))
  bws <- bw
  converged <- FALSE
  for (iter in seq_len(as.integer(max_iter))) {
    old <- B
    if (family == "binomial") {
      eta0 <- as.vector(rowSums(XX * B[unit, , drop = FALSE]))
      mu0 <- pmin(pmax(stats::plogis(eta0), 1e-6), 1 - 1e-6)
      var0 <- pmax(mu0 * (1 - mu0), 1e-6)
      working_y <- eta0 + (yy - mu0) / var0
    }
    for (j in seq_len(ncol(XX))) {
      if (family == "gaussian") {
        partial <- yy - rowSums(XX * B[unit, , drop = FALSE]) + XX[, j] * B[unit, j]
        working_weights <- NULL
      } else {
        partial <- working_y - rowSums(XX * B[unit, , drop = FALSE]) + XX[, j] * B[unit, j]
        working_weights <- var0
      }
      if (auto_bw) {
        score <- vapply(candidates, function(h) {
          .mgwpr_cv(partial, XX[, j], unit, D, h, adaptive, kernel,
                    observation_weights = working_weights)
        }, numeric(1))
        bws[j] <- candidates[which.min(score)]
      }
      B[, j] <- .mgwpr_smooth(
        partial, XX[, j], unit, D, bws[j], adaptive, kernel,
        observation_weights = working_weights
      )
    }
    delta <- max(abs(B - old))
    if (verbose) message(sprintf("Iteration %d: max coefficient change = %.6g", iter, delta))
    if (is.finite(delta) && delta < tol) { converged <- TRUE; break }
  }
  coef_scale <- y_scale / x_scale
  names(coef_scale) <- colnames(XX)
  B <- sweep(B, 2, x_scale, "/") * y_scale
  B_fit <- sweep(sweep(B, 2, x_scale, "*"), 2, y_scale, "/")
  fitted_fit <- rowSums(XX * B_fit[unit, , drop = FALSE])
  residuals_fit <- yy - fitted_fit
  context_effects <- NULL
  if (method == "fe") {
    Xbar <- t(vapply(seq_len(m), function(g) colMeans(Xraw[unit == g, , drop = FALSE]), numeric(ncol(Xraw))))
    ybar <- vapply(seq_len(m), function(g) mean(y_raw[unit == g]), numeric(1))
    context_effects <- ybar - rowSums(Xbar * B)
    fit <- rowSums(Xraw * B[unit, , drop = FALSE]) + context_effects[unit]
  } else fit <- rowSums(X0 * B[unit, , drop = FALSE])
  if (scale && method == "pooled") fit <- fit + y_center - rowSums(sweep(B, 2, x_center, "*"))[unit]
  if (family == "binomial") {
    eta <- as.vector(rowSums(XX * B_fit[unit, , drop = FALSE]))
    fit <- stats::plogis(eta)
    fitted_fit <- fit
    residuals_fit <- yy - fit
  }
  resid <- y_raw - fit
  B_output <- B
  intercept_term <- which(colnames(B_output) == "(Intercept)")
  if (scale && method == "pooled" && length(intercept_term) == 1L) {
    B_output[, intercept_term] <- B_output[, intercept_term] + y_center -
      rowSums(sweep(B, 2, x_center, "*"))
  }
  names(bws) <- colnames(XX)
  rownames(B_output) <- rownames(B) <- rownames(xy) <- id_values
  panel_data <- data
  panel_data$.mgwpr_unit <- unit
  panel_data$.mgwpr_time <- if (is.null(time)) {
    ave(seq_len(n), unit, FUN = seq_along)
  } else time
  structure(list(coefficients = B_output, bandwidths = bws, fitted.values = fit,
                 residuals = resid, coords = xy, unit = unit, formula = formula,
                 method = method, family = family, response_levels = response_levels,
                 kernel = kernel, adaptive = adaptive,
                 context_effects = context_effects,
                 nobs = n, nloc = m, iterations = iter, converged = converged,
                 response = y_raw,
                 x = X0, x_raw = Xraw, x_fit = XX, y_fit = yy,
                 coefficients_fit = B_fit, residuals_fit = residuals_fit,
                 coef_scale = coef_scale,
                 distance = D,
                 spatial = spatial_aligned, time = time,
                 spatial_input_sf = !is.null(spatial),
                 panel_data = panel_data, id_name = id_name, time_name = time_name,
                 call = match.call()), class = "mgwpr")
}

.mgwpr_weights <- function(D, focal, bw, adaptive, kernel) {
  d <- D[focal, ]
  if (adaptive) {
    ord <- order(d)
    h <- d[ord[min(length(ord), max(2L, as.integer(bw)))]]
  } else h <- bw
  if (!is.finite(h) || h <= 0) h <- max(d[d > 0], 1)
  z <- d / h
  switch(kernel,
    bisquare = pmax(0, (1 - z^2)^2),
    gaussian = exp(-.5 * z^2),
    exponential = exp(-abs(z)),
    tricube = pmax(0, (1 - abs(z)^3)^3),
    boxcar = as.numeric(abs(z) <= 1),
    stop("Unknown kernel.")
  )
}

.mgwpr_local <- function(y, x, unit, weights, observation_weights = NULL) {
  w <- weights[unit]
  if (!is.null(observation_weights)) w <- w * observation_weights
  sxx <- sum(w * x * x); sxy <- sum(w * x * y)
  if (!is.finite(sxx) || sxx < 1e-12) return(0)
  sxy / sxx
}

.mgwpr_smooth <- function(y, x, unit, D, bw, adaptive, kernel,
                          observation_weights = NULL) {
  m <- nrow(D); out <- numeric(m)
  for (g in seq_len(m)) out[g] <- .mgwpr_local(
    y, x, unit, .mgwpr_weights(D, g, bw, adaptive, kernel), observation_weights
  )
  out
}

.mgwpr_cv <- function(y, x, unit, D, bw, adaptive, kernel,
                      observation_weights = NULL) {
  m <- nrow(D); pred <- numeric(length(y))
  for (g in seq_len(m)) {
    w <- .mgwpr_weights(D, g, bw, adaptive, kernel)
    w[g] <- 0
    b <- .mgwpr_local(y, x, unit, w, observation_weights)
    rows <- which(unit == g)
    pred[rows] <- b * x[rows]
  }
  loss <- (y - pred)^2
  if (!is.null(observation_weights)) loss <- loss * observation_weights
  mean(loss)
}

#' @export
coef.mgwpr <- function(object, ...) object$coefficients
#' @export
fitted.mgwpr <- function(object, ...) object$fitted.values
#' @export
residuals.mgwpr <- function(object, ...) object$residuals

#' @export
summary.mgwpr <- function(object, ...) {
  y <- object$response
  r <- object$residuals
  structure(list(call = object$call, method = object$method,
                 family = object$family,
                 nobs = object$nobs, nloc = object$nloc,
                 bandwidths = object$bandwidths, coefficients = object$coefficients,
                 context_effects = object$context_effects,
                 sigma = sqrt(mean(r^2)), rmse = sqrt(mean(r^2)),
                 r.squared = if (object$family == "gaussian") {
                   1 - sum(r^2) / sum((y - mean(y))^2)
                 } else NA_real_,
                 iterations = object$iterations, converged = object$converged),
            class = "summary.mgwpr")
}

#' @export
print.mgwpr <- function(x, ...) {
  cat("Multiscale geographically weighted panel regression\n")
  cat("Family:", x$family, "\n")
  cat("Method:", if (x$method == "fe") "individual fixed effects (within)" else "pooled", "\n")
  cat("Locations:", x$nloc, " Observations:", x$nobs, "\n")
  cat("Backfitting iterations:", x$iterations, " Converged:", x$converged, "\n")
  cat("Bandwidths:\n"); print(x$bandwidths)
  invisible(x)
}

#' @export
print.summary.mgwpr <- function(x, ...) {
  cat("Multiscale geographically weighted panel regression\n")
  cat("Method:", x$method, " | Locations:", x$nloc, " | Observations:", x$nobs, "\n")
  cat("RMSE:", signif(x$rmse, 5), " R-squared:", signif(x$r.squared, 5), "\n")
  cat("Converged:", x$converged, "after", x$iterations, "iterations\n\n")
  cat("Bandwidths:\n"); print(x$bandwidths)
  cat("\nLocation-specific coefficients (first six locations):\n")
  print(utils::head(x$coefficients))
  if (!is.null(x$context_effects)) {
    cat("\nRecovered location-specific context effects (first six):\n")
    print(utils::head(x$context_effects))
  }
  invisible(x)
}
