#' Calculate conditional local inference for an MGWPR fit
#'
#' Builds the coefficient-specific linear smoother operators from the fitted
#' spatial kernels and bandwidths. It returns local t-values using a
#' panel-unit cluster sandwich variance and coefficient-specific effective
#' parameter counts, ENP_j = tr(H_j). Inference is conditional
#' on the selected bandwidths and does not account for bandwidth-selection
#' uncertainty.
#'
#' @param object An object returned by [mgwpr()].
#' @param tol Convergence tolerance for the linear smoother backfitting.
#' @param max_iter Maximum linear smoother backfitting iterations.
#' @param finite_sample Apply a unit-cluster finite-sample correction to the
#'   sandwich variance.
#' @return A list with `t_values`, `standard_errors`, `enp`, `df`, and
#'   `iterations`, of class `mgwpr_inference`.
#' @export
coefficient_inference <- function(object, tol = 1e-8, max_iter = 500L,
                                  finite_sample = TRUE) {
  if (!inherits(object, "mgwpr")) stop("object must be an mgwpr fit.")
  if (!is.null(object$family) && object$family != "gaussian") {
    stop("coefficient_inference() currently supports family='gaussian' only.")
  }
  required <- c("x_fit", "y_fit", "coefficients_fit", "residuals_fit",
                "unit", "coords", "bandwidths", "coef_scale")
  if (!all(required %in% names(object))) {
    stop("This fit does not contain the design data required for inference; refit with MGWPR >= 0.5.0.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one positive finite number.")
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop("max_iter must be a positive integer.")
  }
  if (!is.logical(finite_sample) || length(finite_sample) != 1L || is.na(finite_sample)) {
    stop("finite_sample must be TRUE or FALSE.")
  }

  X <- as.matrix(object$x_fit)
  y <- as.numeric(object$y_fit)
  resid <- as.numeric(object$residuals_fit)
  unit <- as.integer(object$unit)
  coef_fit <- as.matrix(object$coefficients_fit)
  n <- nrow(X); m <- nrow(coef_fit); p <- ncol(X)
  if (length(y) != n || length(resid) != n || length(unit) != n ||
      ncol(coef_fit) != p || max(unit) != m) {
    stop("The fitted object has inconsistent stored model dimensions.")
  }
  terms <- colnames(X)
  if (is.null(terms)) terms <- paste0("V", seq_len(p))
  ids <- rownames(coef_fit)
  if (is.null(ids)) ids <- as.character(seq_len(m))
  if (length(object$bandwidths) != p) stop("Stored bandwidths do not match the model terms.")

  distance <- if (!is.null(object$distance)) object$distance else as.matrix(stats::dist(object$coords))
  A <- vector("list", p)
  D <- vector("list", p)
  for (j in seq_len(p)) {
    xj <- X[, j]
    Dj <- matrix(0, n, m)
    Dj[cbind(seq_len(n), unit)] <- xj
    D[[j]] <- Dj

    Aj <- matrix(0, m, n)
    for (g in seq_len(m)) {
      wloc <- .mgwpr_weights(distance, g, object$bandwidths[j],
                             object$adaptive, object$kernel)
      wobs <- wloc[unit]
      denom <- sum(wobs * xj^2)
      if (is.finite(denom) && denom > .Machine$double.eps) {
        Aj[g, ] <- wobs * xj / denom
      }
    }
    A[[j]] <- Aj
  }

  # C[[j]][[k]] maps coefficient influence for term k into the next
  # backfitting update for term j.
  C <- lapply(seq_len(p), function(j) {
    lapply(seq_len(p), function(k) if (j == k) NULL else A[[j]] %*% D[[k]])
  })
  L <- lapply(seq_len(p), function(j) matrix(0, m, n))
  converged <- FALSE
  for (iter in seq_len(as.integer(max_iter))) {
    old <- L
    for (j in seq_len(p)) {
      updated <- A[[j]]
      if (p > 1L) for (k in seq_len(p)) {
        if (k != j) updated <- updated - C[[j]][[k]] %*% L[[k]]
      }
      L[[j]] <- updated
    }
    delta <- max(vapply(seq_len(p), function(j) max(abs(L[[j]] - old[[j]])), numeric(1)))
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) warning("Linear smoother backfitting did not converge within max_iter.")

  enp <- stats::setNames(numeric(p), terms)
  se <- matrix(NA_real_, m, p, dimnames = list(ids, terms))
  correction <- if (finite_sample && m > 1L) m / (m - 1L) else 1
  for (j in seq_len(p)) {
    Hj <- D[[j]] %*% L[[j]]
    # MGWR's term-specific effective parameter count is the trace of the
    # covariate-specific component smoother, not the whole-model residual-DF
    # expression 2 tr(H) - tr(H' H).
    enp[j] <- sum(diag(Hj))

    for (g in seq_len(m)) {
      influence <- L[[j]][g, ]
      cluster_score <- vapply(seq_len(m), function(unit_id) {
        rows <- which(unit == unit_id)
        sum(influence[rows] * resid[rows])
      }, numeric(1))
      variance <- correction * sum(cluster_score^2)
      if (is.finite(variance) && variance >= 0) se[g, j] <- sqrt(variance)
    }
  }
  estimates <- object$coefficients
  se <- sweep(se, 2, object$coef_scale[colnames(se)], "*")
  t_values <- estimates / se
  t_values[!is.finite(t_values)] <- NA_real_

  structure(
    list(t_values = t_values, standard_errors = se, enp = enp,
         df = max(m - 1L, 1L), iterations = iter, converged = converged,
         variance = "panel-unit cluster sandwich",
         bandwidth_conditioning = TRUE),
    class = "mgwpr_inference"
  )
}

#' Extract local t-values from MGWPR inference
#'
#' @param object An object returned by [mgwpr()].
#' @param ... Arguments passed to [coefficient_inference()].
#' @return A location-by-coefficient matrix of local t-values.
#' @export
mgwpr_tvalues <- function(object, ...) {
  coefficient_inference(object, ...)$t_values
}

#' Extract coefficient-specific effective parameter counts
#'
#' @param object An object returned by [mgwpr()].
#' @param ... Arguments passed to [coefficient_inference()].
#' @return A named vector of coefficient-specific ENP_j values.
#' @export
mgwpr_enp <- function(object, ...) {
  coefficient_inference(object, ...)$enp
}
