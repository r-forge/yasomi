\encoding{UTF-8}
\name{error.kaskilagus}
\alias{error.kaskilagus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kaski and Lagus' error measure for Self-Organising Maps}
\description{
Compute an error measure of a fitted Self-Organising Maps defined by
S. Kaski and K. Lagus.  This measure takes into account both the
quantisation error and a type of continuity (i.e., close observations
in the original space should be mapped to close units in the grid).
}
\usage{
error.kaskilagus(som, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{data}{a matrix or a data frame of observations}
}
%\details{
%}
\value{
  The value of the error.
}
\references{S. Kaski, K. Lagus, Comparing self-organizing maps, in:
  C. von der Malsburg, W. von Seelen, J. Vorbrüggen, B. Sendhoff (eds.),
  Proceedings of International Conference on Artificial Neural Networks
  (ICANN'96), vol. 1112 of Lecture Notes in Computer Science, Springer,
  Berlin, Germany, 1996, pp. 809--814.}
\author{Fabrice Rossi}
%\note{   further notes   
%}
\seealso{\code{\link{error.quantisation}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
som <-
som.tune(sg,data,som.tunecontrol(sg,,criterion=error.kaskilagus))$best.som

# and display the prototypes
plot(som)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
