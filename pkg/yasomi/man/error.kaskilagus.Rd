\encoding{UTF-8}
\name{error.kaskilagus}
\alias{error.kaskilagus}
\title{Kaski and Lagus' error measure for Self-Organising Maps}
\description{
Compute an error measure of a fitted Self-Organising Maps defined by
S. Kaski and K. Lagus.  This measure takes into account both the
quantisation error and a type of continuity (i.e., close observations
in the original space should be mapped to close units in the
grid). \code{error.quantisation} is the generic function for this
calculation. 
}
\usage{
error.kaskilagus(som, newdata, \dots)
}
\arguments{
  \item{som}{an object of class derived \code{"som"}}
  \item{newdata}{optional observations (see details)}
  \item{\dots}{additional arguments to be passed to methods}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the error made by
  the fitted som on the data used to fit it. Those data must have been
  saved in the som object (this is the default behavior of
  \code{\link{batchsom}}).

  When \code{newdata} is specified, the function returns the error made
  by the fitted som on the corresponding data. The representation of the
  new data depends on the actual SOM.
}
\references{
  Kaski, S. and Lagus, K. (1996) Comparing self-organizing maps, in:
  C. von der Malsburg, W. von Seelen, J. Vorbrüggen, B. Sendhoff (eds.),
  Proceedings of International Conference on Artificial Neural Networks
  (ICANN'96, Bochum, Germany), vol. 1112 of Lecture Notes in Computer
  Science, Springer, pp. 809--814.
}
\author{Fabrice Rossi}
\seealso{\code{\link{error.quantisation}},
  \code{\link{som.tunecontrol}},\code{\link{error.kaskilagus.somnum}}
  for numerical SOM and \code{\link{error.kaskilagus.relationalsom}} for
  relational SOM} 
%\examples{
%}
\keyword{cluster}
