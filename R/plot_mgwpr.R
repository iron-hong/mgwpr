#' Plot location-specific coefficient surfaces
#'
#' Creates faceted tile maps for one or more coefficient surfaces in an
#' MGWPR fit. Coordinates are taken from the fitted object.
#'
#' @param object An object returned by mgwpr.
#' @param variables Coefficient names to plot. Defaults to all coefficients.
#' @param ncol Number of facet columns.
#' @param tile_width,tile_height Optional tile dimensions in coordinate units.
#'   If omitted, ggplot2 infers them from the coordinate spacing.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
#' @export
plot_coefficients <- function(object, variables = NULL, ncol = 2,
                              tile_width = NULL, tile_height = NULL,
                              xlab = "X coordinate",
                              ylab = "Y coordinate") {
  .check_mgwpr_plot_object(object)
  available <- colnames(object$coefficients)
  if (is.null(variables)) variables <- available
  if (!length(variables) || any(!variables %in% available)) {
    stop("variables must be coefficient names in the fitted model.")
  }
  if (length(ncol) != 1L || is.na(ncol) || ncol < 1) {
    stop("ncol must be a positive integer.")
  }

  long <- .mgwpr_long_map_data(
    object, variables,
    lapply(variables, function(v) object$coefficients[, v]),
    "term", "estimate"
  )
  long$term <- factor(long$term, levels = variables)

  p <- ggplot2::ggplot(long, ggplot2::aes(fill = estimate))
  if (isTRUE(object$spatial_input_sf)) {
    p <- p + ggplot2::geom_sf(color = NA)
    coord_layer <- ggplot2::coord_sf(datum = NA)
  } else {
    p <- p + ggplot2::aes(x = x, y = y) +
      .mgwpr_tile_layer(tile_width, tile_height)
    coord_layer <- ggplot2::coord_equal()
  }
  p + ggplot2::facet_wrap(~term, ncol = as.integer(ncol)) +
    ggplot2::scale_fill_viridis_c(option = "C", name = "Estimate") +
    coord_layer +
    ggplot2::labs(x = xlab, y = ylab, title = "Spatial coefficient surfaces") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Plot recovered unit-level context effects
#'
#' For a fixed-effects fit, maps the recovered location-specific context
#' effects. An optional truth vector adds a reference panel.
#'
#' @param object An MGWPR fit with method "fe".
#' @param truth Optional numeric vector of true context effects, in the same
#'   order as the fitted locations. A named vector is aligned by location ID.
#' @param tile_width,tile_height Optional tile dimensions in coordinate units.
#'   If omitted, ggplot2 infers them from the coordinate spacing.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
#' @export
plot_context_effects <- function(object, truth = NULL,
                                 tile_width = NULL, tile_height = NULL,
                                 xlab = "X coordinate",
                                 ylab = "Y coordinate") {
  .check_mgwpr_plot_object(object)
  if (object$method != "fe" || is.null(object$context_effects)) {
    stop("Context effects are available only for method='fe'.")
  }
  ids <- rownames(object$coefficients)
  estimate <- as.numeric(object$context_effects)
  long <- .mgwpr_long_map_data(object, "Recovered MGWFER", list(estimate),
                               "model", "estimate")
  if (!is.null(truth)) {
    if (!is.numeric(truth)) stop("truth must be a numeric vector.")
    if (!is.null(names(truth))) {
      if (!all(ids %in% names(truth))) {
        stop("Named truth values must include every fitted location ID.")
      }
      truth <- truth[ids]
    }
    if (length(truth) != length(estimate)) {
      stop("truth must have one value per fitted location.")
    }
    reference <- .mgwpr_long_map_data(
      object, "True context effect", list(as.numeric(truth)), "model", "estimate"
    )
    long <- rbind(reference, long)
  }
  long$model <- factor(
    long$model,
    levels = c("True context effect", "Recovered MGWFER")
  )
  limits <- .mgwpr_plot_limits(long$estimate)

  p <- ggplot2::ggplot(long, ggplot2::aes(fill = estimate))
  if (isTRUE(object$spatial_input_sf)) {
    p <- p + ggplot2::geom_sf(color = NA)
    coord_layer <- ggplot2::coord_sf(datum = NA)
  } else {
    p <- p + ggplot2::aes(x = x, y = y) +
      .mgwpr_tile_layer(tile_width, tile_height)
    coord_layer <- ggplot2::coord_equal()
  }
  p + ggplot2::facet_wrap(~model, nrow = 1) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = limits, name = "Context effect"
    ) +
    coord_layer +
    ggplot2::labs(
      x = xlab, y = ylab, title = "Spatial context effects"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Compare pooled and fixed-effects coefficient maps
#'
#' Maps one coefficient from pooled and fixed-effects MGWPR fits. Optionally
#' includes the true coefficient surface as a reference panel.
#'
#' @param pooled A fit returned by mgwpr with method "pooled".
#' @param fixed A fit returned by mgwpr with method "fe".
#' @param variable Name of the coefficient to map.
#' @param truth Optional numeric vector of true values, in the same order as
#'   the pooled fit locations. A named vector is aligned by location ID.
#' @param tile_width,tile_height Optional tile dimensions in coordinate units.
#'   If omitted, ggplot2 infers them from the coordinate spacing.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
#' @export
plot_mgwpr_comparison <- function(pooled, fixed, variable, truth = NULL,
                                  tile_width = NULL, tile_height = NULL,
                                  xlab = "X coordinate",
                                  ylab = "Y coordinate") {
  .check_mgwpr_plot_object(pooled)
  .check_mgwpr_plot_object(fixed)
  if (pooled$method != "pooled") stop("pooled must be fit with method='pooled'.")
  if (fixed$method != "fe") stop("fixed must be fit with method='fe'.")
  if (length(variable) != 1L || !variable %in% colnames(pooled$coefficients) ||
      !variable %in% colnames(fixed$coefficients)) {
    stop("variable must be a coefficient present in both fitted models.")
  }

  ids <- rownames(pooled$coefficients)
  fixed_index <- match(ids, rownames(fixed$coefficients))
  if (anyNA(fixed_index) || nrow(pooled$coefficients) != nrow(fixed$coefficients)) {
    stop("The pooled and fixed-effects fits must contain the same locations.")
  }
  if (!isTRUE(all.equal(
    unname(pooled$coords),
    unname(fixed$coords[fixed_index, , drop = FALSE]),
    tolerance = 1e-8
  ))) {
    stop("The pooled and fixed-effects coordinates do not match by location.")
  }

  long <- .mgwpr_long_map_data(pooled, "PMGWR", list(pooled$coefficients[, variable]),
                               "model", "estimate")
  if (!is.null(truth)) {
    if (!is.numeric(truth)) stop("truth must be a numeric vector.")
    if (!is.null(names(truth))) {
      if (!all(ids %in% names(truth))) {
        stop("Named truth values must include every fitted location ID.")
      }
      truth <- truth[ids]
    }
    if (length(truth) != length(ids)) {
      stop("truth must have one value per fitted location.")
    }
    reference <- .mgwpr_long_map_data(pooled, "True", list(as.numeric(truth)),
                                      "model", "estimate")
    long <- rbind(reference, long)
  }
  fe_panel <- .mgwpr_long_map_data(
    pooled, "MGWFER", list(fixed$coefficients[fixed_index, variable]),
    "model", "estimate"
  )
  long <- rbind(long, fe_panel)
  long$model <- factor(long$model, levels = c("True", "PMGWR", "MGWFER"))
  limits <- .mgwpr_plot_limits(long$estimate)

  p <- ggplot2::ggplot(long, ggplot2::aes(fill = estimate))
  if (isTRUE(pooled$spatial_input_sf)) {
    p <- p + ggplot2::geom_sf(color = NA)
    coord_layer <- ggplot2::coord_sf(datum = NA)
  } else {
    p <- p + ggplot2::aes(x = x, y = y) +
      .mgwpr_tile_layer(tile_width, tile_height)
    coord_layer <- ggplot2::coord_equal()
  }
  p + ggplot2::facet_wrap(~model, nrow = 1) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = limits, name = variable
    ) +
    coord_layer +
    ggplot2::labs(
      x = xlab, y = ylab,
      title = paste("Spatial coefficient comparison:", variable)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' @export
plot.mgwpr <- function(x, type = c("coefficients", "context"),
                       variables = NULL, ncol = 2,
                       tile_width = NULL, tile_height = NULL, ...) {
  type <- match.arg(type)
  if (type == "context") {
    return(plot_context_effects(
      x, tile_width = tile_width, tile_height = tile_height, ...
    ))
  }
  plot_coefficients(
    x, variables = variables, ncol = ncol,
    tile_width = tile_width, tile_height = tile_height, ...
  )
}

.check_mgwpr_plot_object <- function(object) {
  if (!inherits(object, "mgwpr") && !inherits(object, "mgwr")) {
    stop("object must be an mgwpr or mgwr fit.")
  }
  if (is.null(object$coefficients) || is.null(object$coords)) {
    stop("The fitted object is missing coefficients or coordinates.")
  }
  invisible(TRUE)
}

.mgwpr_plot_limits <- function(x) {
  limits <- range(x, finite = TRUE)
  if (!all(is.finite(limits))) stop("No finite values are available to plot.")
  if (diff(limits) == 0) {
    padding <- max(abs(limits[1]) * 0.01, 0.01)
    limits <- limits + c(-padding, padding)
  }
  limits
}

.mgwpr_tile_layer <- function(tile_width = NULL, tile_height = NULL) {
  args <- list()
  if (!is.null(tile_width)) {
    if (!is.numeric(tile_width) || length(tile_width) != 1L ||
        !is.finite(tile_width) || tile_width <= 0) {
      stop("tile_width must be NULL or one positive finite number.")
    }
    args$width <- tile_width
  }
  if (!is.null(tile_height)) {
    if (!is.numeric(tile_height) || length(tile_height) != 1L ||
        !is.finite(tile_height) || tile_height <= 0) {
      stop("tile_height must be NULL or one positive finite number.")
    }
    args$height <- tile_height
  }
  do.call(ggplot2::geom_tile, args)
}

.mgwpr_long_map_data <- function(object, groups, values, group_name, value_name) {
  if (length(groups) != length(values)) stop("Internal map inputs do not align.")
  ids <- rownames(object$coefficients)
  index <- match(ids, as.character(object$spatial$location_id))
  if (anyNA(index)) stop("Stored sf geometries do not match fitted locations.")
  pieces <- lapply(seq_along(groups), function(i) {
    if (isTRUE(object$spatial_input_sf)) {
      piece <- object$spatial[index, , drop = FALSE]
    } else {
      piece <- data.frame(x = object$coords[, 1], y = object$coords[, 2])
    }
    piece[[group_name]] <- groups[i]
    piece[[value_name]] <- as.numeric(values[[i]])
    piece
  })
  if (length(pieces) == 1L) return(pieces[[1L]])
  do.call(rbind, pieces)
}

#' Return fitted results joined to spatial geometries
#'
#' @param object An object returned by [mgwpr()].
#' @param include_context Include recovered fixed effects when available.
#' @return An `sf` object with location-specific coefficients and geometry.
#' @export
mgwpr_sf <- function(object, include_context = TRUE) {
  .check_mgwpr_plot_object(object)
  if (!inherits(object$spatial, "sf")) stop("The fitted object has no sf geometries.")
  out <- object$spatial
  idx <- match(rownames(object$coefficients), as.character(out$location_id))
  if (anyNA(idx)) stop("Stored geometries do not match fitted locations.")
  out <- out[idx, , drop = FALSE]
  for (term in colnames(object$coefficients)) {
    out[[paste0("coef_", make.names(term))]] <- object$coefficients[, term]
  }
  if (include_context && !is.null(object$context_effects)) {
    out$context_effect <- as.numeric(object$context_effects)
  }
  out
}
