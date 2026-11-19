\name{diagnose_moran}
\alias{diagnose_moran}
\alias{diagnose_local_f}
\alias{diagnose_hausman}
\alias{diagnose_lm}
\alias{diagnose_mgwpr}
\title{Spatial and panel diagnostics for MGWPR fits}
\usage{
diagnose_moran(object, spatial_weights = NULL, k = 4L,
  permutations = 999L, seed = NULL)

diagnose_local_f(object, bw = NULL)

diagnose_hausman(object, bw = NULL)

diagnose_lm(object, bw = NULL)

diagnose_mgwpr(object,
  tests = c("moran", "f_test", "hausman", "lm_test"),
  spatial_weights = NULL, k = 4L, permutations = 999L,
  seed = NULL, bw = NULL)
}
\arguments{
  \item{object}{An object returned by mgwpr().}
  \item{spatial_weights}{Optional non-negative location-by-location spatial weights matrix.}
  \item{k}{Number of nearest neighbours for default Moran weights.}
  \item{permutations}{Number of Moran permutation replicates; zero omits the permutation p-value.}
  \item{seed}{Optional seed for Moran permutations. The caller's RNG state is restored.}
  \item{bw}{Optional local diagnostic bandwidth; defaults to the median fitted bandwidth.}
  \item{tests}{Diagnostics to run: moran, f_test, hausman, and/or lm_test.}
}
\description{
Computes period-specific residual Moran's I with a permutation p-value;
location-specific weighted F tests comparing pooled and unit fixed effects;
local Hausman tests comparing within and random-effects slopes; and local
Breusch-Pagan LM tests for random intercepts.
}
\details{
Local F, Hausman, and LM tests require a Gaussian outcome. The Hausman
diagnostic uses a local Swamy-Arora-style variance-components approximation.
The LM formula assumes balanced local panels; unbalanced neighborhoods are
returned with a non-testable status. For binomial fits, Moran's I uses Pearson
residuals and is exploratory.
}
\value{
The individual diagnostics return data frames or a list (for Moran's I).
diagnose_mgwpr() returns a named list of requested diagnostic results.
}
\seealso{mgwpr, coefficient_inference, coefficient_significance}
\keyword{spatial}
\keyword{models}
