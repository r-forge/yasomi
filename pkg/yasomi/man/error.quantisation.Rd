\name{error.quantisation}
\alias{error.quantisation}
\title{Quantisation error for a Self-Organising Map}
\description{
  \code{error.quantisation} is a generic function for computing the
  quantisation error made by a fitted Self-Organising Map on a
  dataset. The exact definition depends on the data type, but should be
  compatible with the standard definition for vector data; in this case the quantisation error  is the mean of the Euclidean distance between each
observation and the prototype of its best matching unit.
}
\usage{
error.quantisation(som, newdata, \dots)
}
\arguments{
  \item{som}{an object of class derived from \code{"som"}}
  \item{newdata}{optional observations (see details)}
  \item{\dots}{additional arguments to be passed to methods}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the quantisation
  error made by the fitted som on the data used to fit it. 

  When \code{newdata} is specified, the function returns the
  quantisation error of the fitted som on the corresponding data. The
  representation of the new data depends on the actual SOM.  }
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus}}, \code{\link{som.tunecontrol}},
  \code{\link{error.quantisation.somnum}} for numerical SOM and
  \code{\link{error.quantisation.relationalsom}} for relational SOM}
%\examples{}
\keyword{cluster}
