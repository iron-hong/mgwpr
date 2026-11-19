\name{plot_coefficients}
\alias{plot_coefficients}
\alias{plot_context_effects}
\alias{plot_mgwpr_comparison}
\alias{plot.mgwpr}
\title{Spatial Maps for MGWPR Fits}
\description{
Creates coefficient-surface maps, fixed-effects context maps, and pooled
versus fixed-effects coefficient comparisons from MGWPR model objects.
}
\usage{
plot_coefficients(object, variables = NULL, ncol = 2,
  tile_width = NULL, tile_height = NULL,
  xlab = "X coordinate", ylab = "Y coordinate")

plot_context_effects(object, truth = NULL,
  tile_width = NULL, tile_height = NULL,
  xlab = "X coordinate", ylab = "Y coordinate")

plot_mgwpr_comparison(pooled, fixed, variable, truth = NULL,
  tile_width = NULL, tile_height = NULL,
  xlab = "X coordinate", ylab = "Y coordinate")

\method{plot}{mgwpr}(x, type = c("coefficients", "context"),
  variables = NULL, ncol = 2, tile_width = NULL, tile_height = NULL, ...)
}
\arguments{
  \item{object,x}{An MGWPR fitted object.}
  \item{pooled}{An MGWPR fit with method "pooled".}
  \item{fixed}{An MGWPR fit with method "fe".}
  \item{variables}{Coefficient names to map. NULL maps all coefficients.}
  \item{variable}{One coefficient name present in both fits.}
  \item{truth}{Optional numeric reference vector, one value per location.
    Named vectors are aligned to fitted location identifiers.}
  \item{type}{Map type for plot(): "coefficients" or "context".}
  \item{ncol}{Number of columns in the coefficient facet layout.}
  \item{tile_width,tile_height}{Optional tile dimensions in coordinate units.
    If omitted, dimensions are inferred from the coordinate spacing.}
  \item{xlab,ylab}{Axis labels for the coordinate columns.}
  \item{...}{Additional arguments passed to the selected plotting function.}
}
\value{
A ggplot object.
}
\details{
The functions draw filled tiles centered on the coordinates stored in the
fitted object. Comparison plots require pooled and fixed-effects fits to
contain the same locations and coordinates. Context-effect maps are available
only for fixed-effects fits. Inferred tiles are suitable for regular grids;
irregular boundaries require supplied polygon geometry, which is not yet
supported.
}
\examples{
\dontrun{
plot_coefficients(fit, variables = c("x1", "x2"))
plot(fit, type = "context")
plot_mgwpr_comparison(fit_pmgwr, fit_mgwfer, variable = "x1")
}
}
\keyword{spatial}
\keyword{hplot}
