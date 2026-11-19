\name{mgwpr_sf}
\alias{mgwpr_sf}
\title{Join fitted MGWPR results to spatial geometries}
\usage{
mgwpr_sf(object, include_context = TRUE)
}
\arguments{
  \item{object}{An object returned by mgwpr().}
  \item{include_context}{Include recovered unit fixed effects when available.}
}
\description{
Returns an sf object with the fitted locations' original geometries, local
coefficient estimates, and optionally recovered context effects.
}
\value{An sf object keyed by location_id. Coefficient columns are prefixed
with coef_.}
\seealso{mgwpr, plot_coefficients, coefficient_significance}
\keyword{spatial}
