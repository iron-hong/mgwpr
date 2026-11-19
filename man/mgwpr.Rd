\name{mgwpr}
\alias{mgwpr}
\alias{coef.mgwpr}
\alias{fitted.mgwpr}
\alias{residuals.mgwpr}
\alias{summary.mgwpr}
\alias{print.mgwpr}
\alias{print.summary.mgwpr}
\title{Multiscale Geographically Weighted Panel Regression}
\description{
Fits pooled or individual fixed-effects multiscale geographically weighted
panel regression by iterative backfitting. Gaussian and pooled binomial
responses are supported, with a separate spatial bandwidth for each model
matrix coefficient.
}
\usage{
mgwpr(formula, data, coords = NULL, spatial = NULL, id = NULL, time = NULL,
  method = c("fe", "pooled"), family = c("gaussian", "binomial"),
  bw = NULL, adaptive = TRUE,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"), selector = "cv",
  bw_min = NULL, bw_max = NULL, max_iter = 30L, tol = 1e-5,
  scale = FALSE, verbose = FALSE)
}
\arguments{
  \item{formula}{A model formula.}
  \item{data}{A data frame in long panel format.}
  \item{coords}{Names of the two numeric coordinate columns in data. Supply either coords or spatial.}
  \item{spatial}{Optional sf object containing one point or polygon geometry per panel unit and the ID column named by id. Polygon representative points are used for distance calculations; original geometries are retained for mapping.}
  \item{id}{Panel-unit identifier column name or vector. Required for method="fe".}
  \item{time}{Optional time column name or vector; used for validation only.}
  \item{method}{"fe" for individual fixed effects or "pooled".}
  \item{family}{"gaussian" for a continuous outcome or "binomial" for a binary outcome coded 0/1 or as a two-level factor. Binomial models currently require method="pooled".}
  \item{bw}{Optional positive scalar or vector of bandwidths, one per model-matrix coefficient.}
  \item{adaptive}{If TRUE, use neighbor-count bandwidths; otherwise coordinate distances.}
  \item{kernel}{Spatial kernel: "bisquare", "gaussian", "exponential", "tricube", or "boxcar".}
  \item{selector}{Bandwidth criterion; currently "cv" only.}
  \item{bw_min,bw_max}{Optional lower and upper bandwidth search bounds.}
  \item{max_iter}{Maximum number of backfitting iterations.}
  \item{tol}{Convergence tolerance on the largest coefficient change.}
  \item{scale}{Whether to center and scale the outcome and design matrix before fitting.}
  \item{verbose}{Whether to print iteration progress.}
  \item{object,x}{An object returned by mgwpr or summary.mgwpr.}
  \item{...}{Further arguments to methods (currently ignored).}
}
\value{
mgwpr returns an object of class "mgwpr" with location-level coefficients,
fixed-effects context effects when applicable, bandwidths, fitted values,
residuals, and convergence information.
coef, fitted, residuals, and summary extract or summarize these components.
}
\details{
Fixed effects are removed by individual within transformation. Binomial
responses are fitted using pooled geographically weighted logistic regression.
Bandwidths are selected by leave-one-location-out cross-validation. The sf
input is matched to panel units using id; locations must have one geometry
each. Moran's I, local F, Hausman, and Breusch-Pagan LM diagnostics are
available through the diagnose_* functions. Local F, Hausman, and LM
diagnostics are intended for Gaussian fits. The Hausman implementation uses a
local variance-components approximation, and diagnostics should be interpreted
with that limitation in mind.
}
\references{
Li, Z. and Fotheringham, A. S. (2026). Spatial Context as a Time-Invariant
Confounder: A Fixed-Effects Extension of MGWR. Annals of the American
Association of Geographers, 116(9), 2294-2318.
}
\seealso{summary.mgwpr, plot_coefficients, plot_context_effects,
plot_mgwpr_comparison, mgwpr_sf, diagnose_mgwpr}
\keyword{spatial}
\keyword{panel}
