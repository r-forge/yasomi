\encoding{UTF-8}
\name{error.kaskilagus.somnum}
\alias{error.kaskilagus.somnum}
\title{Kaski and Lagus' error measure for Self-Organising Maps on vector
data}
\description{
Compute an error measure of a fitted Self-Organising Maps defined by
S. Kaski and K. Lagus. For vector data, the error combines for each
observation its quantisation error and a graph based distance in the
prior structure compute with the Euclidean distance between the
prototypes (in the data space).
}
\usage{
\method{error.kaskilagus}{somnum}(som, newdata,\dots)
}
\arguments{
  \item{som}{an object of class \code{"somnum"}}
  \item{newdata}{an optional matrix or a data frame of observations}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the error made by
  the fitted som on the data used to fit it. Those data must have been
  saved in the som object (this is the default behavior of
  \code{\link{batchsom}}). When \code{newdata} is specified, the
  function returns the error made by the fitted som on the corresponding
  data.
}
\references{
  Kaski, S. and Lagus, K. (1996) Comparing self-organizing maps, in:
  C. von der Malsburg, W. von Seelen, J. Vorbrüggen, B. Sendhoff (eds.),
  Proceedings of International Conference on Artificial Neural Networks
  (ICANN'96, Bochum, Germany), vol. 1112 of Lecture Notes in Computer
  Science, Springer, pp. 809--814.
}
\author{Fabrice Rossi}
\seealso{\code{\link{error.quantisation.somnum}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
st <- som.tune(data,sg,som.tunecontrol(sg,criterion=error.kaskilagus))
som <- st$best.som

# and display the prototypes
plot(som)
}
\keyword{cluster}
