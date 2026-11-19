\name{compare_coefficient_surfaces}
\alias{compare_coefficient_surfaces}
\alias{coef_metrics}
\title{Compare Estimated and True Coefficient Surfaces}
\description{
Aligns pooled and fixed-effects coefficient estimates with known coefficient
surfaces by location, and calculates RMSE and correlation measures.
}
\usage{
compare_coefficient_surfaces(pooled, fixed, truth, id = "unit_id",
  variables = NULL, truth_vars = NULL)

coef_metrics(pooled, fixed, truth, id = "unit_id",
  variables = NULL, truth_vars = NULL)
}
\arguments{
  \item{pooled}{An MGWPR fit with method "pooled".}
  \item{fixed}{An MGWPR fit with method "fe".}
  \item{truth}{A data frame with a location ID and numeric true coefficient columns.}
  \item{id}{Name of the location ID column in truth.}
  \item{variables}{Coefficient names to compare. Defaults to common non-intercept terms.}
  \item{truth_vars}{Optional named character vector mapping model terms to truth columns,
    for example c(x1 = "beta1_true").}
}
\value{
compare_coefficient_surfaces returns a list of class "mgwpr_comparison" with
one row per location in coefficients and one row per coefficient in metrics.
coef_metrics returns only the metrics data frame.
}
\details{
The output metrics include pooled and fixed-effects RMSE and correlations.
Correlation is returned as NA for constant truth surfaces, such as a
coefficient that is zero at every location.
}
\examples{
\dontrun{
comparison <- compare_coefficient_surfaces(
  fit_pmgwr, fit_mgwfer, truth = units
)
comparison$coefficients
comparison$metrics
coef_metrics(fit_pmgwr, fit_mgwfer, truth = units)
}
}
\keyword{spatial}
\keyword{models}
