#' Fit multiscale geographically weighted regression
#'
#' Fits a cross-sectional MGWR model with \pkg{GWmodel}. By default,
#' \code{GWmodel::gwr.multiscale()} first selects an initial bandwidth for
#' every model-matrix term and then optimizes the term-specific bandwidths
#' during back-fitting. This step may take some time.
#'
#' @param formula A regression formula with a numeric response.
#' @param data A data frame with one row per spatial observation. Repeated
#'   panel rows are not supported by this cross-sectional MGWR function.
#' @param coords Names of two numeric coordinate columns in \code{data}.
#'   Supply either \code{coords} or \code{spatial}. If \code{data} is already
#'   an sf object, its geometry is used automatically.
#' @param spatial An sf object with one geometry per row of \code{data}. Rows
#'   must be in the same order as the model data.
#' @param id Optional location-ID column name or vector. If supplied, IDs must
#'   be unique.
#' @param kernel Spatial kernel supported by GWmodel: \code{"bisquare"},
#'   \code{"gaussian"}, \code{"exponential"}, \code{"tricube"}, or
#'   \code{"boxcar"}.
#' @param adaptive If TRUE, bandwidths are nearest-neighbour counts; otherwise
#'   they are distances in the coordinate reference units.
#' @param approach Criterion used by GWmodel to select bandwidths, such as
#'   \code{"AICc"}, \code{"AIC"}, \code{"BIC"}, or \code{"CV"}.
#' @param criterion Back-fitting convergence criterion, \code{"dCVR"} or
#'   \code{"CVR"}.
#' @param max_iter Maximum number of back-fitting iterations.
#' @param tol Convergence threshold passed to GWmodel.
#' @param bws0 Optional initial bandwidth for each model-matrix term. If NULL,
#'   GWmodel selects initial bandwidths automatically.
#' @param bws_thresholds Optional bandwidth-change thresholds, one per model
#'   term or a single value.
#' @param bws_reopts Number of consecutive converged selections required
#'   before a term bandwidth is fixed.
#' @param hatmatrix If TRUE, calculate the hat matrix, standard errors, and
#'   t-values. This can increase memory and computation time.
#' @param predictor_centered Logical vector indicating which predictors are
#'   centered internally. Defaults to centering all non-intercept terms.
#' @param nlower Lower adaptive-neighbour bound used for initial bandwidth
#'   selection by GWmodel.
#' @param longlat If NULL, infer geographic coordinates from the sf CRS;
#'   otherwise pass the value to GWmodel.
#' @param crs CRS assigned when \code{coords} are converted to sf points.
#' @param parallel_method Parallel approach passed to GWmodel; FALSE, "omp",
#'   "cluster", or "cuda" depending on the GWmodel installation.
#' @param parallel_arg Optional parallel setting passed to GWmodel.
#' @param verbose If TRUE, print GWmodel bandwidth-selection progress.
#' @param ... Additional named arguments passed to
#'   \code{GWmodel::gwr.multiscale()}.
#' @return An object of class \code{mgwr}. It contains local coefficients,
#'   selected bandwidths, fitted values, residuals, an sf result layer, model
#'   diagnostics, and the original GWmodel fit in \code{gwmodel_fit}.
#' @export
mgwr <- function(formula, data, coords = NULL, spatial = NULL, id = NULL,
                 kernel = "bisquare", adaptive = FALSE,
                 approach = "AICc", criterion = "dCVR",
                 max_iter = 2000L, tol = 1e-5, bws0 = NULL,
                 bws_thresholds = NULL, bws_reopts = 5L,
                 hatmatrix = TRUE, predictor_centered = NULL,
                 nlower = 10L, longlat = NULL, crs = NA,
                 parallel_method = FALSE, parallel_arg = NULL,
                 verbose = FALSE, ...) {
  if (!requireNamespace("GWmodel", quietly = TRUE)) {
    stop("Package 'GWmodel' is required for mgwr(). Install it with install.packages('GWmodel').")
  }
  if (!is.data.frame(data)) stop("data must be a data frame or sf object.")
  if (inherits(data, "sf") && is.null(spatial) && is.null(coords)) {
    spatial <- data
  }
  if (!is.null(spatial) && !is.null(coords)) {
    stop("Supply either coords or spatial, not both.")
  }
  if (is.null(spatial) && is.null(coords)) {
    stop("Supply coords, spatial, or an sf data object.")
  }
  if (!kernel %in% c("bisquare", "gaussian", "exponential", "tricube", "boxcar")) {
    stop("kernel must be one of: bisquare, gaussian, exponential, tricube, boxcar.")
  }
  if (!approach %in% c("AICc", "AIC", "BIC", "CV", "aic", "bic", "cv", "aicc")) {
    stop("approach must be AICc, AIC, BIC, or CV.")
  }
  approach <- switch(toupper(approach), AICC = "AICc", AIC = "AIC",
                     BIC = "BIC", CV = "CV")
  if (!criterion %in% c("dCVR", "CVR")) stop("criterion must be 'dCVR' or 'CVR'.")
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) stop("max_iter must be a positive integer.")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive.")
  if (length(bws_reopts) != 1L || !is.finite(bws_reopts) || bws_reopts < 1 ||
      bws_reopts != as.integer(bws_reopts)) stop("bws_reopts must be a positive integer.")

  map_as_sf <- is.null(coords) && (inherits(data, "sf") || inherits(spatial, "sf"))
  model_data <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  if (nrow(model_data) < 4L) stop("At least four spatial observations are required.")
  if (!is.null(spatial) && !inherits(spatial, "sf")) {
    stop("spatial must be an sf object.")
  }
  if (inherits(spatial, "sf") && nrow(spatial) != nrow(model_data)) {
    stop("spatial must contain exactly one geometry per row of data, in the same order.")
  }

  if (!is.null(coords)) {
    if (length(coords) != 2L || !all(coords %in% names(model_data))) {
      stop("coords must name two columns in data.")
    }
    if (!all(vapply(model_data[coords], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(model_data[coords])))) {
      stop("Coordinate columns must be numeric, finite, and complete.")
    }
    spatial_data <- sf::st_as_sf(model_data, coords = coords, remove = FALSE, crs = crs)
  } else {
    if (!inherits(spatial, "sf")) stop("spatial must be an sf object.")
    spatial_data <- sf::st_sf(model_data, geometry = sf::st_geometry(spatial))
  }

  mf <- stats::model.frame(formula, data = model_data, na.action = stats::na.fail)
  if (nrow(mf) != nrow(model_data)) {
    stop("Missing model values are not supported; remove incomplete rows first.")
  }
  response <- stats::model.response(mf)
  if (!is.numeric(response) || any(!is.finite(response))) {
    stop("MGWR requires a complete, finite, numeric response.")
  }
  design <- stats::model.matrix(attr(mf, "terms"), mf)
  if (any(!is.finite(design))) stop("Predictors must be complete and finite.")
  terms <- colnames(design)

  if (!is.null(id)) {
    id_values <- if (is.character(id) && length(id) == 1L && id %in% names(model_data)) {
      model_data[[id]]
    } else id
    if (length(id_values) != nrow(model_data) || anyNA(id_values) || anyDuplicated(id_values)) {
      stop("id must contain one unique, non-missing ID per spatial observation.")
    }
  } else {
    id_values <- seq_len(nrow(model_data))
  }
  location_id <- as.character(id_values)

  if (is.null(longlat)) longlat <- isTRUE(sf::st_is_longlat(spatial_data))
  if (!is.logical(longlat) || length(longlat) != 1L || is.na(longlat)) {
    stop("longlat must be TRUE, FALSE, or NULL.")
  }
  if (is.null(predictor_centered)) {
    predictor_centered <- rep(TRUE, ncol(design) -
      as.integer("(Intercept)" %in% colnames(design)))
  }
  if (!is.null(bws_thresholds) &&
      !(length(bws_thresholds) %in% c(1L, ncol(design)))) {
    stop("bws_thresholds must be NULL, a scalar, or one value per model-matrix term.")
  }
  if (!is.null(bws_thresholds) && any(!is.finite(bws_thresholds) | bws_thresholds < 0)) {
    stop("bws_thresholds must contain non-negative finite values.")
  }
  if (!is.null(bws0) && (!is.numeric(bws0) ||
                         !(length(bws0) %in% c(1L, ncol(design))) ||
                         any(!is.finite(bws0) | bws0 <= 0))) {
    stop("bws0 must be NULL, a positive scalar, or one positive value per model-matrix term.")
  }
  if (length(bws0) == 1L && ncol(design) > 1L) {
    bws0 <- rep(bws0, ncol(design))
  }
  if (length(bws_thresholds) == 1L && ncol(design) > 1L) {
    bws_thresholds <- rep(bws_thresholds, ncol(design))
  }
  n_predictors <- ncol(design) - as.integer("(Intercept)" %in% colnames(design))
  if (length(predictor_centered) == 1L) {
    predictor_centered <- rep(predictor_centered, n_predictors)
  }
  if (!is.logical(predictor_centered) || anyNA(predictor_centered) ||
      length(predictor_centered) != n_predictors) {
    stop("predictor_centered must be a logical scalar or a vector for all non-intercept terms.")
  }

  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("Additional GWmodel arguments in ... must be named.")
  }
  if (any(names(dots) %in% c("formula", "data"))) {
    stop("Pass formula and data through mgwr(), not through ... .")
  }

  gw_args <- list(
    formula = formula,
    data = spatial_data,
    kernel = kernel,
    adaptive = adaptive,
    criterion = criterion,
    max.iterations = as.integer(max_iter),
    threshold = tol,
    dMats = NULL,
    var.dMat.indx = NULL,
    p.vals = NULL,
    theta.vals = NULL,
    longlat = longlat,
    bws0 = bws0,
    bw.seled = NULL,
    approach = approach,
    bws.thresholds = bws_thresholds,
    bws.reOpts = as.integer(bws_reopts),
    verbose = verbose,
    hatmatrix = hatmatrix,
    predictor.centered = predictor_centered,
    nlower = nlower,
    parallel.method = parallel_method,
    parallel.arg = parallel_arg
  )
  for (nm in names(dots)) gw_args[[nm]] <- dots[[nm]]

  fit <- do.call(GWmodel::gwr.multiscale, gw_args)
  sdf <- fit$SDF
  if (!inherits(sdf, "sf")) sdf <- sf::st_as_sf(sdf)
  sdf_names <- names(sf::st_drop_geometry(sdf))
  if (ncol(design) > length(sdf_names)) {
    stop("GWmodel returned an unexpected coefficient layout in SDF.")
  }
  names(sdf)[seq_len(ncol(design))] <- terms
  sdf$location_id <- location_id
  coef_mat <- as.matrix(sf::st_drop_geometry(sdf)[, terms, drop = FALSE])
  rownames(coef_mat) <- location_id
  bw_values <- fit$GW.arguments$bws
  names(bw_values) <- terms

  fitted_values <- if ("yhat" %in% names(sdf)) sdf$yhat else rep(NA_real_, nrow(sdf))
  residual_values <- if ("residual" %in% names(sdf)) sdf$residual else response - fitted_values
  geometry <- sf::st_geometry(sdf)
  map_geometry <- if (any(sf::st_is(geometry, c("POLYGON", "MULTIPOLYGON")))) {
    sf::st_point_on_surface(geometry)
  } else geometry
  coords_out <- sf::st_coordinates(map_geometry)[, 1:2, drop = FALSE]

  structure(list(
    coefficients = coef_mat,
    bandwidths = bw_values,
    fitted.values = as.numeric(fitted_values),
    residuals = as.numeric(residual_values),
    response = as.numeric(response),
    coords = coords_out,
    unit = seq_len(nrow(model_data)),
    formula = formula,
    method = "mgwr",
    family = "gaussian",
    adaptive = adaptive,
    kernel = kernel,
    context_effects = NULL,
    nobs = nrow(model_data),
    nloc = nrow(model_data),
    spatial = sdf,
    spatial_input_sf = map_as_sf,
    gwmodel_fit = fit,
    GW.diagnostic = fit$GW.diagnostic,
    call = match.call()
  ), class = "mgwr")
}

#' @export
coef.mgwr <- function(object, ...) object$coefficients

#' @export
fitted.mgwr <- function(object, ...) object$fitted.values

#' @export
residuals.mgwr <- function(object, ...) object$residuals

#' @export
print.mgwr <- function(x, ...) {
  cat("Multiscale geographically weighted regression (GWmodel)\n")
  cat("Observations:", x$nobs, "\n")
  cat("Kernel:", x$kernel, " | Adaptive:", x$adaptive, "\n")
  cat("Selected bandwidths:\n")
  print(x$bandwidths)
  if (!is.null(x$GW.diagnostic)) {
    cat("\nGWmodel diagnostics:\n")
    print(x$GW.diagnostic)
  }
  invisible(x)
}

#' @export
summary.mgwr <- function(object, ...) {
  structure(list(
    call = object$call,
    nobs = object$nobs,
    kernel = object$kernel,
    adaptive = object$adaptive,
    bandwidths = object$bandwidths,
    coefficients = object$coefficients,
    diagnostics = object$GW.diagnostic
  ), class = "summary.mgwr")
}

#' @export
print.summary.mgwr <- function(x, ...) {
  cat("Multiscale geographically weighted regression (GWmodel)\n")
  cat("Observations:", x$nobs, "\n")
  cat("Kernel:", x$kernel, " | Adaptive:", x$adaptive, "\n\n")
  cat("Selected bandwidths:\n")
  print(x$bandwidths)
  cat("\nLocation-specific coefficients (first six locations):\n")
  print(utils::head(x$coefficients))
  if (!is.null(x$diagnostics)) {
    cat("\nGWmodel diagnostics:\n")
    print(x$diagnostics)
  }
  invisible(x)
}

#' @export
plot.mgwr <- function(x, variables = NULL, ncol = 2,
                      tile_width = NULL, tile_height = NULL, ...) {
  plot_coefficients(x, variables = variables, ncol = ncol,
                    tile_width = tile_width, tile_height = tile_height, ...)
}
