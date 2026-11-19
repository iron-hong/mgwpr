\name{mgwr}
\alias{mgwr}
\alias{coef.mgwr}
\alias{fitted.mgwr}
\alias{residuals.mgwr}
\alias{print.mgwr}
\alias{summary.mgwr}
\alias{print.summary.mgwr}
\alias{plot.mgwr}
\title{Fit Multiscale Geographically Weighted Regression}
\description{
Fits a cross-sectional MGWR model through \code{GWmodel::gwr.multiscale()}.
GWmodel selects an initial bandwidth for each model-matrix term and then
optimizes term-specific bandwidths during back-fitting. This selection can
take time.
}
\usage{
mgwr(formula, data, coords = NULL, spatial = NULL, id = NULL,
  kernel = "bisquare", adaptive = FALSE, approach = "AICc",
  criterion = "dCVR", max_iter = 2000L, tol = 1e-5, bws0 = NULL,
  bws_thresholds = NULL, bws_reopts = 5L, hatmatrix = TRUE,
  predictor_centered = NULL, nlower = 10L, longlat = NULL, crs = NA,
  parallel_method = FALSE, parallel_arg = NULL, verbose = FALSE, ...)

\method{coef}{mgwr}(object, ...)
\method{fitted}{mgwr}(object, ...)
\method{residuals}{mgwr}(object, ...)
\method{print}{mgwr}(x, ...)
\method{summary}{mgwr}(object, ...)
\method{print}{summary.mgwr}(x, ...)
\method{plot}{mgwr}(x, variables = NULL, ncol = 2,
  tile_width = NULL, tile_height = NULL, ...)
}
\arguments{
  \item{formula}{A regression formula with a numeric response.}
  \item{data}{A data frame or \code{sf} object with one row per spatial observation. Repeated panel rows are not supported.}
  \item{coords}{Names of two numeric coordinate columns in \code{data}. Supply either \code{coords} or \code{spatial}.}
  \item{spatial}{An \code{sf} object with one geometry per model-data row, in the same order.}
  \item{id}{Optional location-ID column name or vector. IDs must be unique.}
  \item{kernel}{One of \code{"bisquare"}, \code{"gaussian"}, \code{"exponential"}, \code{"tricube"}, or \code{"boxcar"}.}
  \item{adaptive}{If \code{TRUE}, bandwidths are nearest-neighbour counts; otherwise they are distances in coordinate units.}
  \item{approach}{Bandwidth-selection criterion: \code{"AICc"}, \code{"AIC"}, \code{"BIC"}, or \code{"CV"}.}
  \item{criterion}{Back-fitting convergence criterion: \code{"dCVR"} or \code{"CVR"}.}
  \item{max_iter}{Maximum number of back-fitting iterations.}
  \item{tol}{Convergence threshold passed to GWmodel.}
  \item{bws0}{Optional initial bandwidth, either a scalar or one value per model-matrix term.}
  \item{bws_thresholds}{Optional bandwidth-change thresholds, scalar or one per model-matrix term.}
  \item{bws_reopts}{Number of consecutive converged selections required before fixing a term bandwidth.}
  \item{hatmatrix}{If \code{TRUE}, calculate the hat matrix, standard errors, and t-values.}
  \item{predictor_centered}{Logical scalar or vector indicating which non-intercept predictors are centered internally.}
  \item{nlower}{Lower adaptive-neighbour bound for initial bandwidth selection.}
  \item{longlat}{Whether coordinates are longitude and latitude. If \code{NULL}, inferred from the \code{sf} CRS.}
  \item{crs}{CRS assigned when coordinate columns are converted to \code{sf} points.}
  \item{parallel_method}{Parallel approach passed to GWmodel.}
  \item{parallel_arg}{Optional parallel setting passed to GWmodel.}
  \item{verbose}{If \code{TRUE}, print bandwidth-selection progress.}
  \item{...}{Additional named arguments passed to \code{GWmodel::gwr.multiscale()}; for methods, arguments passed to the corresponding generic.}
  \item{object}{An MGWR fit.}
  \item{x}{An MGWR fit or summary.}
  \item{variables}{Coefficient names to plot. Defaults to all coefficient surfaces.}
  \item{ncol}{Number of facet columns.}
  \item{tile_width, tile_height}{Optional tile dimensions for coordinate-based maps.}
}
\details{
This is a cross-sectional MGWR interface. For panel data, select a single
period or aggregate to one observation per spatial unit before fitting. The
package's \code{mgwpr()} function provides panel estimators.
}
\value{
\code{mgwr()} returns an object of class \code{mgwr}. It contains location-
specific coefficients, bandwidths, fitted values, residuals, an \code{sf}
result layer, diagnostics, and the underlying GWmodel object in
\code{gwmodel_fit}. The S3 methods return or display these results;
\code{plot()} returns a coefficient-surface plot.
}
\references{
Gollini, I., Lu, B., Charlton, M., Brunsdon, C., and Harris, P. (2015).
GWmodel: An R Package for Exploring Spatial Heterogeneity Using Geographically
Weighted Models. \emph{Journal of Statistical Software}, 63(17), 1--50.
\url{https://doi.org/10.18637/jss.v063.i17}
}
\seealso{
\code{GWmodel::gwr.multiscale()}, \code{mgwpr()}, \code{plot_coefficients()}
}
\examples{
\dontrun{
# Requires the optional GWmodel package and one row per spatial location.
fit <- mgwr(y ~ x1 + x2, data = cross_section,
            coords = c("longitude", "latitude"), adaptive = TRUE)
summary(fit)
plot_coefficients(fit)
}
}
