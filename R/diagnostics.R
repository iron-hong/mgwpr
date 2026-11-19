#' Moran's I diagnostic for panel residuals
#'
#' Computes Moran's I for residuals separately by panel period and a
#' permutation p-value for the mean statistic. If no spatial weights are
#' supplied, inverse-distance weights are constructed from the fitted
#' locations using the nearest `k` neighbours.
#'
#' @param object An object returned by [mgwpr()].
#' @param spatial_weights Optional non-negative `nloc` by `nloc` spatial weight
#'   matrix. The diagonal is set to zero and rows are standardized.
#' @param k Number of nearest neighbours used to construct default weights.
#' @param permutations Number of random permutations; set to zero to omit the
#'   permutation p-value.
#' @param seed Optional random seed. The caller's RNG state is restored.
#' @return A list with period-specific results, mean Moran's I, permutation
#'   p-value, and the row-standardized weights.
#' @export
diagnose_moran <- function(object, spatial_weights = NULL, k = 4L,
                           permutations = 999L, seed = NULL) {
  .check_mgwpr_diag_object(object)
  if (length(k) != 1L || is.na(k) || k < 1L || k != as.integer(k)) {
    stop("k must be a positive integer.")
  }
  if (length(permutations) != 1L || is.na(permutations) || permutations < 0 ||
      permutations != as.integer(permutations)) stop("permutations must be a non-negative integer.")
  m <- object$nloc
  if (is.null(spatial_weights)) {
    k <- min(as.integer(k), m - 1L)
    D <- object$distance
    W <- matrix(0, m, m)
    for (i in seq_len(m)) {
      ord <- order(D[i, ], seq_len(m))
      nbr <- ord[ord != i][seq_len(k)]
      W[i, nbr] <- 1 / pmax(D[i, nbr], .Machine$double.eps)
    }
  } else {
    W <- as.matrix(spatial_weights)
    if (!is.numeric(W) || !identical(dim(W), c(m, m)) ||
        any(!is.finite(W)) || any(W < 0)) {
      stop("spatial_weights must be a finite, non-negative nloc by nloc matrix.")
    }
  }
  diag(W) <- 0
  rs <- rowSums(W)
  W[rs > 0, ] <- W[rs > 0, , drop = FALSE] / rs[rs > 0]
  if (sum(W) <= 0) stop("spatial_weights has no positive off-diagonal weights.")

  unit <- object$unit
  time <- object$panel_data$.mgwpr_time
  residual <- object$residuals
  residual_type <- "model residuals"
  if (identical(object$family, "binomial")) {
    fitted_prob <- pmin(pmax(object$fitted.values, 1e-8), 1 - 1e-8)
    residual <- (object$response - fitted_prob) / sqrt(fitted_prob * (1 - fitted_prob))
    residual_type <- "Pearson residuals (exploratory for a binomial model)"
  }
  has_time <- !is.null(object$time)
  periods <- if (!has_time) "all" else unique(time)
  residual_by_period <- lapply(periods, function(tt) {
    out <- rep(NA_real_, m)
    rows <- if (!has_time) seq_along(time) else which(time == tt)
    for (g in unique(unit[rows])) {
      ix <- rows[unit[rows] == g]
      out[g] <- mean(residual[ix])
    }
    out
  })
  names(residual_by_period) <- as.character(periods)

  moran_value <- function(z, Wsub) {
    keep <- is.finite(z) & rowSums(Wsub) > 0
    if (sum(keep) < 3L) return(NA_real_)
    z <- z[keep] - mean(z[keep])
    Wsub <- Wsub[keep, keep, drop = FALSE]
    s0 <- sum(Wsub)
    den <- sum(z^2)
    if (s0 <= 0 || den <= .Machine$double.eps) return(NA_real_)
    length(z) / s0 * sum(z * as.vector(Wsub %*% z)) / den
  }
  observed <- vapply(residual_by_period, moran_value, numeric(1), Wsub = W)
  per_period <- data.frame(
    time = names(residual_by_period),
    moran_i = as.numeric(observed),
    n = vapply(residual_by_period, function(z) sum(is.finite(z)), integer(1)),
    stringsAsFactors = FALSE
  )
  statistic <- if (all(!is.finite(observed))) NA_real_ else mean(observed, na.rm = TRUE)
  p_value <- NA_real_
  if (permutations > 0L && is.finite(statistic)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    sims <- replicate(as.integer(permutations), {
      permutation <- sample.int(m)
      vals <- vapply(residual_by_period, function(z) {
        moran_value(z[permutation], Wsub = W)
      }, numeric(1))
      if (all(!is.finite(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    })
    sims <- sims[is.finite(sims)]
    if (length(sims)) p_value <- (1 + sum(abs(sims) >= abs(statistic))) / (length(sims) + 1)
  }
  list(by_period = per_period, statistic = statistic,
       p_value = p_value, permutations = as.integer(permutations),
       weights = W, residual_type = residual_type)
}

#' Local F test comparing pooled and unit fixed-effects panel fits
#'
#' At each focal location, compares spatially weighted local least-squares
#' models with and without unit indicators. Tests are unavailable for binomial
#' fits.
#'
#' @param object An object returned by [mgwpr()].
#' @param bw Optional diagnostic bandwidth. Defaults to the median fitted
#'   bandwidth; adaptive values are neighbour counts.
#' @return A data frame with location-specific F statistics, degrees of
#'   freedom, p-values, and status.
#' @export
diagnose_local_f <- function(object, bw = NULL) {
  .check_mgwpr_diag_object(object)
  if (object$family != "gaussian") return(.mgwpr_unavailable_diag(object, "local_f", "Gaussian response required."))
  .mgwpr_local_panel_test(object, test = "f", bw = bw)
}

#' Local Hausman test comparing within and random-effects slopes
#'
#' Computes a local weighted Hausman statistic at each spatial unit. The
#' variance-component estimate uses a Swamy-Arora-style within/between
#' decomposition. Tests are unavailable for binomial fits.
#'
#' @param object An object returned by [mgwpr()].
#' @param bw Optional diagnostic bandwidth. Defaults to the median fitted
#'   bandwidth; adaptive values are neighbour counts.
#' @return A data frame with location-specific Hausman statistics, degrees of
#'   freedom, p-values, and status.
#' @export
diagnose_hausman <- function(object, bw = NULL) {
  .check_mgwpr_diag_object(object)
  if (object$family != "gaussian") return(.mgwpr_unavailable_diag(object, "hausman", "Gaussian response required."))
  .mgwpr_local_panel_test(object, test = "hausman", bw = bw)
}

#' Local Breusch-Pagan LM test for random panel effects
#'
#' Computes a local Breusch-Pagan LM statistic for random intercepts versus
#' pooling. The balanced-panel formula is used; locations with unbalanced local
#' neighborhoods are marked as not testable.
#'
#' @param object An object returned by [mgwpr()].
#' @param bw Optional diagnostic bandwidth. Defaults to the median fitted
#'   bandwidth; adaptive values are neighbour counts.
#' @return A data frame with location-specific LM statistics, degrees of
#'   freedom, p-values, and status.
#' @export
diagnose_lm <- function(object, bw = NULL) {
  .check_mgwpr_diag_object(object)
  if (object$family != "gaussian") return(.mgwpr_unavailable_diag(object, "lm", "Gaussian response required."))
  .mgwpr_local_panel_test(object, test = "lm", bw = bw)
}

#' Run selected MGWPR diagnostics
#'
#' @param object An object returned by [mgwpr()].
#' @param tests Diagnostics to run: "moran", "f_test", "hausman", and/or
#'   "lm_test".
#' @param spatial_weights Optional weights matrix passed to [diagnose_moran()].
#' @param k Neighbours for automatic Moran weights.
#' @param permutations Permutations for the Moran test.
#' @param seed Optional seed for the Moran permutation test.
#' @param bw Optional local diagnostic bandwidth for F, Hausman, and LM tests.
#' @return A named list of diagnostic results.
#' @export
diagnose_mgwpr <- function(object,
                           tests = c("moran", "f_test", "hausman", "lm_test"),
                           spatial_weights = NULL, k = 4L,
                           permutations = 999L, seed = NULL, bw = NULL) {
  tests <- match.arg(tests, c("moran", "f_test", "hausman", "lm_test"), several.ok = TRUE)
  out <- list()
  if ("moran" %in% tests) out$moran <- diagnose_moran(object, spatial_weights, k, permutations, seed)
  if ("f_test" %in% tests) out$f_test <- diagnose_local_f(object, bw)
  if ("hausman" %in% tests) out$hausman <- diagnose_hausman(object, bw)
  if ("lm_test" %in% tests) out$lm_test <- diagnose_lm(object, bw)
  structure(out, class = "mgwpr_diagnostics")
}

.check_mgwpr_diag_object <- function(object) {
  if (!inherits(object, "mgwpr")) stop("object must be an mgwpr fit.")
  if (is.null(object$residuals) || is.null(object$distance) || is.null(object$unit)) {
    stop("The fit is missing residual, distance, or panel-unit data.")
  }
  invisible(TRUE)
}

.mgwpr_unavailable_diag <- function(object, test, message) {
  data.frame(location_id = rownames(object$coefficients), test = test,
             statistic = NA_real_, df = NA_real_, p_value = NA_real_,
             status = message, stringsAsFactors = FALSE)
}

.mgwpr_local_panel_test <- function(object, test, bw = NULL) {
  m <- object$nloc
  if (is.null(object$panel_data$.mgwpr_unit)) stop("The fit is missing unit IDs.")
  if (is.null(bw)) bw <- stats::median(object$bandwidths, na.rm = TRUE)
  if (length(bw) != 1L || !is.finite(bw) || bw <= 0) stop("bw must be one positive finite value.")
  unit <- object$unit
  y <- object$response
  X <- object$x_raw
  ids <- rownames(object$coefficients)
  result <- lapply(seq_len(m), function(g) {
    wloc <- .mgwpr_weights(object$distance, g, bw, object$adaptive, object$kernel)
    wobs_all <- wloc[unit]
    keep <- which(wobs_all > 1e-10)
    yg <- y[keep]
    Xg <- X[keep, , drop = FALSE]
    ug <- unit[keep]
    wg <- wobs_all[keep]
    n <- length(keep)
    ng <- length(unique(ug))
    fail <- function(msg) data.frame(location_id = ids[g], test = test,
                                     statistic = NA_real_, df = NA_real_,
                                     p_value = NA_real_, status = msg)
    if (n <= ncol(Xg) + ng + 1L || ng < 2L) return(fail("insufficient_local_sample"))

    if (test == "f") {
      # FE fits omit the intercept after within transformation. Rebuild the
      # conventional pooled intercept for the nested pooled-versus-unit-FE test.
      Xs <- Xg[, colnames(Xg) != "(Intercept)", drop = FALSE]
      Xpool <- cbind(`(Intercept)` = 1, Xs)
      U <- stats::model.matrix(~ factor(ug) - 1L)
      Xfe <- cbind(Xs, U)
      fit0 <- stats::lm.wfit(Xpool, yg, w = wg)
      fit1 <- stats::lm.wfit(Xfe, yg, w = wg)
      df1 <- fit1$rank - fit0$rank
      df2 <- n - fit1$rank
      if (df1 < 1L || df2 < 1L || fit0$rank < ncol(Xpool) ||
          fit1$rank < ncol(Xfe)) return(fail("rank_deficient"))
      ss0 <- sum(wg * fit0$residuals^2)
      ss1 <- sum(wg * fit1$residuals^2)
      stat <- max(0, (ss0 - ss1) / df1) / max(ss1 / df2, .Machine$double.eps)
      return(data.frame(location_id = ids[g], test = test, statistic = stat,
                        df = df1, df2 = df2, p_value = stats::pf(stat, df1, df2, lower.tail = FALSE),
                        status = "ok"))
    }

    Xs <- Xg[, colnames(Xg) != "(Intercept)", drop = FALSE]
    if (ncol(Xs) < 1L) return(fail("no_slope_terms"))
    local_units <- unique(ug)
    group_index <- match(ug, local_units)
    T_i <- tabulate(group_index, nbins = length(local_units))
    means_x <- t(vapply(seq_along(local_units), function(i) colMeans(Xs[group_index == i, , drop = FALSE]), numeric(ncol(Xs))))
    means_y <- vapply(seq_along(local_units), function(i) mean(yg[group_index == i]), numeric(1))
    xwithin <- Xs - means_x[group_index, , drop = FALSE]
    ywithin <- yg - means_y[group_index]
    fe <- stats::lm.wfit(xwithin, ywithin, w = wg)
    if (fe$rank < ncol(Xs)) return(fail("within_rank_deficient"))
    sigma_e <- sum(wg * fe$residuals^2) / max(n - length(local_units) - fe$rank, 1)
    if (!is.finite(sigma_e) || sigma_e <= 0) return(fail("invalid_within_variance"))
    xtwx_fe <- crossprod(xwithin, xwithin * wg)
    Vfe <- tryCatch(sigma_e * solve(xtwx_fe), error = function(e) NULL)
    if (is.null(Vfe)) return(fail("within_covariance_singular"))

    if (test == "hausman") {
      between_weights <- vapply(seq_along(local_units), function(i) {
        mean(wg[group_index == i]) * T_i[i]
      }, numeric(1))
      between <- stats::lm.wfit(cbind(`(Intercept)` = 1, means_x), means_y,
                                w = between_weights)
      Tbar <- mean(T_i)
      sigma_between <- sum(between_weights * between$residuals^2) / max(length(local_units) - between$rank, 1)
      sigma_u <- max(0, sigma_between - sigma_e / Tbar)
      if (sigma_u <= .Machine$double.eps) return(fail("random_effect_variance_zero"))
      theta <- 1 - sqrt(sigma_e / (sigma_e + T_i * sigma_u))
      xq <- Xs - theta[group_index] * means_x[group_index, , drop = FALSE]
      yq <- yg - theta[group_index] * means_y[group_index]
      re <- stats::lm.wfit(cbind(`(Intercept)` = 1 - theta[group_index], xq), yq, w = wg)
      if (re$rank < ncol(Xs) + 1L) return(fail("random_rank_deficient"))
      sigma_re <- sum(wg * re$residuals^2) / max(n - re$rank, 1)
      Vre_full <- tryCatch(
        sigma_re * solve(crossprod(cbind(1 - theta[group_index], xq),
                                   cbind(1 - theta[group_index], xq) * wg)),
        error = function(e) NULL
      )
      Vre <- if (is.null(Vre_full)) NULL else Vre_full[-1L, -1L, drop = FALSE]
      if (is.null(Vre)) return(fail("random_covariance_singular"))
      d <- fe$coefficients - re$coefficients[-1L]
      Vd <- (Vfe - Vre + t(Vfe - Vre)) / 2
      ee <- tryCatch(eigen(Vd, symmetric = TRUE), error = function(e) NULL)
      if (is.null(ee) || any(!is.finite(ee$values))) return(fail("hausman_covariance_invalid"))
      keep_e <- ee$values > max(abs(ee$values), 1) * 1e-8
      df <- sum(keep_e)
      if (!df) return(fail("hausman_covariance_not_positive"))
      projection <- crossprod(ee$vectors[, keep_e, drop = FALSE], d)
      stat <- sum(as.numeric(projection)^2 / ee$values[keep_e])
      return(data.frame(location_id = ids[g], test = test, statistic = stat,
                        df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE),
                        status = "ok"))
    }

    if (test == "lm") {
      Xpool <- cbind(`(Intercept)` = 1,
                     Xg[, colnames(Xg) != "(Intercept)", drop = FALSE])
      pooled <- stats::lm.wfit(Xpool, yg, w = wg)
      counts <- tabulate(group_index, nbins = length(local_units))
      if (length(unique(counts)) != 1L || min(counts) < 2L) return(fail("unbalanced_local_panel"))
      Tn <- counts[1L]
      rss <- sum(wg * pooled$residuals^2)
      score <- sum(vapply(seq_along(local_units), function(i) {
        rows <- which(group_index == i)
        mean(wg[rows]) * sum(pooled$residuals[rows])^2
      }, numeric(1)))
      ratio <- score / max(rss, .Machine$double.eps)
      stat <- length(local_units) * Tn / (2 * (Tn - 1)) * (ratio - 1)^2
      return(data.frame(location_id = ids[g], test = test, statistic = stat,
                        df = 1, p_value = stats::pchisq(stat, 1, lower.tail = FALSE),
                        status = "ok"))
    }
    fail("unknown_test")
  })
  do.call(rbind, result)
}
