#' Align estimated coefficient surfaces and calculate accuracy metrics
#'
#' Joins pooled and fixed-effects local coefficients to known coefficient
#' surfaces by location ID, then calculates RMSE and correlation by term.
#'
#' @param pooled An MGWPR fit with method "pooled".
#' @param fixed An MGWPR fit with method "fe".
#' @param truth A data frame containing a location ID and true coefficient
#'   columns.
#' @param id Name of the location ID column in truth.
#' @param variables Coefficient names to compare. Defaults to common slope
#'   terms in both fits, excluding the intercept.
#' @param truth_vars Optional named character vector mapping model terms to
#'   truth columns, for example c(x1 = "beta1_true").
#' @return A list with coefficients (one row per location) and metrics (one
#'   row per coefficient), of class mgwpr_comparison.
#' @export
compare_coefficient_surfaces <- function(pooled, fixed, truth, id = "unit_id",
                                         variables = NULL,
                                         truth_vars = NULL) {
  .check_mgwpr_plot_object(pooled)
  .check_mgwpr_plot_object(fixed)
  if (pooled$method != "pooled") stop("pooled must be fit with method='pooled'.")
  if (fixed$method != "fe") stop("fixed must be fit with method='fe'.")
  if (!is.data.frame(truth)) stop("truth must be a data frame.")
  if (length(id) != 1L || !id %in% names(truth)) {
    stop("id must name the location ID column in truth.")
  }
  if (anyNA(truth[[id]]) || anyDuplicated(truth[[id]])) {
    stop("truth must contain one non-missing row per location ID.")
  }

  pooled_terms <- colnames(pooled$coefficients)
  fixed_terms <- colnames(fixed$coefficients)
  available <- intersect(pooled_terms, fixed_terms)
  available <- setdiff(available, c("(Intercept)", "Intercept"))
  if (is.null(variables)) variables <- available
  if (!length(variables) || any(!variables %in% available)) {
    stop("variables must be common non-intercept terms in both fits.")
  }

  if (is.null(truth_vars)) {
    truth_cols <- vapply(variables, function(v) {
      if (v %in% names(truth)) return(v)
      match <- regmatches(v, regexec("^x([0-9]+)$", v))[[1L]]
      candidate <- if (length(match) == 2L) {
        paste0("beta", match[2L], "_true")
      } else {
        NA_character_
      }
      if (!is.na(candidate) && candidate %in% names(truth)) return(candidate)
      stop(
        "No truth column found for coefficient '", v,
        "'. Supply truth_vars, e.g. c(x1 = 'beta1_true').",
        call. = FALSE
      )
    }, character(1))
    names(truth_cols) <- variables
  } else {
    if (is.null(names(truth_vars)) ||
        !all(variables %in% names(truth_vars))) {
      stop("truth_vars must be named for every coefficient in variables.")
    }
    truth_cols <- truth_vars[variables]
    if (anyNA(truth_cols) || any(!truth_cols %in% names(truth))) {
      stop("Every truth_vars entry must name a column in truth.")
    }
  }
  if (anyDuplicated(unname(truth_cols))) {
    stop("Each coefficient must map to a distinct truth column.")
  }
  if (!all(vapply(truth[unname(truth_cols)], is.numeric, logical(1)))) {
    stop("True coefficient columns must be numeric.")
  }

  pooled_ids <- rownames(pooled$coefficients)
  fixed_ids <- rownames(fixed$coefficients)
  if (is.null(pooled_ids) || is.null(fixed_ids) ||
      anyDuplicated(pooled_ids) || anyDuplicated(fixed_ids)) {
    stop("Both fits must have unique location IDs as coefficient row names.")
  }
  fixed_index <- match(pooled_ids, fixed_ids)
  truth_index <- match(pooled_ids, as.character(truth[[id]]))
  if (anyNA(fixed_index)) {
    stop("The pooled and fixed-effects fits do not contain the same locations.")
  }
  if (anyNA(truth_index)) {
    stop("truth is missing one or more locations included in the fitted models.")
  }
  if (!isTRUE(all.equal(
    unname(pooled$coords),
    unname(fixed$coords[fixed_index, , drop = FALSE]),
    tolerance = 1e-8
  ))) {
    stop("The pooled and fixed-effects coordinates do not match by location.")
  }

  result <- data.frame(
    location_id = pooled_ids,
    coord_x = pooled$coords[, 1],
    coord_y = pooled$coords[, 2],
    stringsAsFactors = FALSE
  )
  for (v in variables) {
    result[[paste0(v, "_true")]] <- truth[[truth_cols[[v]]]][truth_index]
    result[[paste0(v, "_pmgwr")]] <- pooled$coefficients[, v]
    result[[paste0(v, "_mgwfer")]] <- fixed$coefficients[fixed_index, v]
  }

  metrics <- do.call(rbind, lapply(variables, function(v) {
    actual <- result[[paste0(v, "_true")]]
    pooled_est <- result[[paste0(v, "_pmgwr")]]
    fixed_est <- result[[paste0(v, "_mgwfer")]]
    data.frame(
      variable = v,
      truth_column = unname(truth_cols[[v]]),
      PMGWR_RMSE = .mgwpr_rmse(pooled_est, actual),
      MGWFER_RMSE = .mgwpr_rmse(fixed_est, actual),
      PMGWR_correlation = .mgwpr_cor(pooled_est, actual),
      MGWFER_correlation = .mgwpr_cor(fixed_est, actual),
      stringsAsFactors = FALSE
    )
  }))
  rownames(metrics) <- NULL

  structure(
    list(coefficients = result, metrics = metrics, variables = variables),
    class = "mgwpr_comparison"
  )
}

#' Extract coefficient-surface accuracy metrics
#'
#' Convenience wrapper returning the metrics table from
#' compare_coefficient_surfaces().
#'
#' @inheritParams compare_coefficient_surfaces
#' @return A data frame of RMSE and correlation by coefficient.
#' @export
coef_metrics <- function(pooled, fixed, truth, id = "unit_id",
                         variables = NULL, truth_vars = NULL) {
  compare_coefficient_surfaces(
    pooled = pooled,
    fixed = fixed,
    truth = truth,
    id = id,
    variables = variables,
    truth_vars = truth_vars
  )$metrics
}

.mgwpr_rmse <- function(estimate, truth) {
  keep <- is.finite(estimate) & is.finite(truth)
  if (!any(keep)) return(NA_real_)
  sqrt(mean((estimate[keep] - truth[keep])^2))
}

.mgwpr_cor <- function(estimate, truth) {
  keep <- is.finite(estimate) & is.finite(truth)
  estimate <- estimate[keep]
  truth <- truth[keep]
  if (length(estimate) < 2L ||
      !is.finite(stats::sd(estimate)) ||
      !is.finite(stats::sd(truth)) ||
      stats::sd(estimate) == 0 ||
      stats::sd(truth) == 0) {
    return(NA_real_)
  }
  stats::cor(estimate, truth)
}
