#' Apply the MGWR effective-number-of-parameters correction to local t-values
#'
#' The function filters supplied local t-values using a per-coefficient
#' significance level alpha / ENP_j, following the multiscale adjustment of
#' da Silva and Fotheringham (2016). It does not estimate standard errors or
#' t-values; those must come from an inferential procedure appropriate to the
#' fitted model.
#'
#' @param object An object returned by [mgwpr()].
#' @param t_values Optional matrix or data frame of local t-values, with one
#'   row per fitted location and one column per coefficient. Supply either
#'   `t_values` or `standard_errors`.
#' @param standard_errors Optional matrix or data frame of coefficient standard
#'   errors with the same dimensions and dimnames as `coef(object)`. T-values
#'   are calculated as coefficient / standard error.
#' @param enp Numeric vector of effective numbers of parameters, one per
#'   coefficient. Names must match coefficient names; an unnamed vector is
#'   interpreted in coefficient order.
#' @param alpha Family-wise target significance level, between 0 and 1.
#' @param df Degrees of freedom used for the two-sided critical t-value. May be
#'   one positive number or a named vector by coefficient. The default `Inf`
#'   uses the standard normal critical value.
#' @return A long data frame with estimates, t-values, adjusted alpha, critical
#'   values, significance flags and filtered t-values; class
#'   `mgwpr_significance`.
#' @references da Silva, A. R. and Fotheringham, A. S. (2016). The multiple
#'   testing issue in geographically weighted regression. *Geographical
#'   Analysis*, 48(3), 233-247. <doi:10.1111/gean.12084>.
#' @export
coefficient_significance <- function(object, t_values = NULL,
                                     standard_errors = NULL, enp,
                                     alpha = 0.05, df = Inf) {
  .check_mgwpr_plot_object(object)
  if (is.null(t_values) == is.null(standard_errors)) {
    stop("Supply exactly one of t_values or standard_errors.")
  }
  if (missing(enp) || !is.numeric(enp) || !length(enp)) {
    stop("enp must contain effective parameter counts.")
  }
  invalid_enp <- !is.finite(enp) | enp <= 0
  if (any(invalid_enp)) {
    bad_terms <- names(enp)[invalid_enp]
    if (is.null(bad_terms)) bad_terms <- which(invalid_enp)
    stop("enp must be positive and finite. Invalid coefficient(s): ",
         paste(bad_terms, collapse = ", "), ".")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number strictly between 0 and 1.")
  }

  estimates <- object$coefficients
  ids <- rownames(estimates)
  terms <- colnames(estimates)
  if (is.null(ids) || is.null(terms)) stop("The fit must have named locations and coefficients.")

  if (is.null(names(enp))) {
    if (length(enp) != length(terms)) stop("Unnamed enp must have one value per coefficient.")
    names(enp) <- terms
  }
  if (anyDuplicated(names(enp)) || !all(terms %in% names(enp))) {
    stop("enp must have one named value for every fitted coefficient.")
  }
  enp <- enp[terms]

  if (is.null(t_values)) {
    standard_errors <- .mgwpr_align_matrix(standard_errors, ids, terms, "standard_errors")
    if (any(standard_errors < 0, na.rm = TRUE)) stop("standard_errors cannot be negative.")
    t_values <- estimates / standard_errors
    zero_se <- !is.na(standard_errors) & standard_errors == 0
    t_values[zero_se] <- NA_real_
  } else {
    t_values <- .mgwpr_align_matrix(t_values, ids, terms, "t_values")
  }

  if (length(df) == 1L) {
    df <- stats::setNames(rep(df, length(terms)), terms)
  } else {
    if (is.null(names(df)) || anyDuplicated(names(df)) || !all(terms %in% names(df))) {
      stop("df must be one positive value or a named value for every coefficient.")
    }
    df <- df[terms]
  }
  if (!is.numeric(df) || anyNA(df) || any(df <= 0)) stop("df must be positive.")

  rows <- lapply(terms, function(term) {
    adjusted_alpha <- min(alpha / enp[[term]], 1)
    critical_t <- stats::qt(1 - adjusted_alpha / 2, df = df[[term]])
    tv <- t_values[, term]
    significant <- is.finite(tv) & abs(tv) > critical_t
    result_term <- data.frame(
      location_id = ids,
      coord_x = object$coords[, 1],
      coord_y = object$coords[, 2],
      term = term,
      estimate = as.numeric(estimates[, term]),
      t_value = as.numeric(tv),
      enp = unname(enp[[term]]),
      alpha_adjusted = adjusted_alpha,
      df = unname(df[[term]]),
      critical_t = critical_t,
      significant = significant,
      direction = ifelse(!significant, "not significant",
                         ifelse(tv > 0, "positive", "negative")),
      filtered_t = ifelse(significant, tv, 0),
      stringsAsFactors = FALSE
    )
    if (isTRUE(object$spatial_input_sf)) {
      spatial_index <- match(ids, as.character(object$spatial$location_id))
      if (anyNA(spatial_index)) stop("Stored sf geometries do not match fitted locations.")
      result_term <- sf::st_sf(
        result_term,
        geometry = sf::st_geometry(object$spatial[spatial_index, , drop = FALSE])
      )
    }
    result_term
  })
  result <- do.call(rbind, rows)
  result$term <- factor(result$term, levels = terms)
  rownames(result) <- NULL
  result_class <- if (inherits(result, "sf")) {
    c("mgwpr_significance", class(result))
  } else {
    c("mgwpr_significance", "data.frame")
  }
  structure(result, class = result_class,
            alpha = alpha, correction = "da Silva-Fotheringham ENP_j")
}

#' Plot coefficient significance by location
#'
#' Maps the output of [coefficient_significance()] as filled polygons or tiles. Significant
#' positive and negative coefficients receive separate colors; non-significant
#' locations are dark blue.
#'
#' @param x An object returned by [coefficient_significance()].
#' @param ncol Number of facet columns.
#' @param tile_width,tile_height Optional tile dimensions in coordinate units.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
#' @export
plot_coefficient_significance <- function(x, ncol = 2,
                                          tile_width = NULL,
                                          tile_height = NULL,
                                          xlab = "X coordinate",
                                          ylab = "Y coordinate") {
  if (!inherits(x, "mgwpr_significance")) {
    stop("x must be returned by coefficient_significance().")
  }
  if (length(ncol) != 1L || is.na(ncol) || ncol < 1) {
    stop("ncol must be a positive integer.")
  }
  plot_data <- x
  plot_data$direction <- factor(
    plot_data$direction,
    levels = c("positive", "negative", "not significant")
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(fill = direction))
  if (inherits(plot_data, "sf")) {
    p <- p + ggplot2::geom_sf(color = NA) + ggplot2::coord_sf(datum = NA)
  } else {
    p <- p + ggplot2::aes(x = coord_x, y = coord_y) +
      .mgwpr_tile_layer(tile_width, tile_height) + ggplot2::coord_equal()
  }
  p +
    ggplot2::facet_wrap(~term, ncol = as.integer(ncol)) +
    ggplot2::scale_fill_manual(
      values = c("positive" = "#E69F00", "negative" = "#B30000",
                 "not significant" = "#17365D"),
      drop = FALSE, name = "Coefficient"
    ) +
    ggplot2::labs(
      x = xlab, y = ylab,
      title = "Local coefficient significance",
      subtitle = paste0("da Silva–Fotheringham correction; target alpha = ",
                        signif(attr(x, "alpha"), 3))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

.mgwpr_align_matrix <- function(x, ids, terms, arg) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) stop(arg, " must be a numeric matrix or data frame.")
  if (!is.null(rownames(x))) {
    if (!all(ids %in% rownames(x))) stop(arg, " is missing one or more fitted locations.")
    x <- x[ids, , drop = FALSE]
  } else if (nrow(x) != length(ids)) {
    stop(arg, " must have one row per fitted location.")
  }
  if (!is.null(colnames(x))) {
    if (!all(terms %in% colnames(x))) stop(arg, " is missing one or more fitted coefficients.")
    x <- x[, terms, drop = FALSE]
  } else {
    if (ncol(x) != length(terms)) stop(arg, " must have one column per fitted coefficient.")
    colnames(x) <- terms
  }
  storage.mode(x) <- "double"
  x
}
