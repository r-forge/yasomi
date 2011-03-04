\encoding{UTF-8}
\name{error.kaskilagus.relationalsom}
\alias{error.kaskilagus.relationalsom}
\title{Kaski and Lagus' error measure for Relational Self-Organising Maps}
\description{
Compute an error measure of a fitted Self-Organising Maps defined by
S. Kaski and K. Lagus. The relational formula is used both to compute
the quantisation error and dissimilarities between prototypes. 
}
\usage{
\method{error.kaskilagus}{relationalsom}(som, newdata, \dots)
}
\arguments{
  \item{som}{an object of class \code{"relationalsom"}}
  \item{newdata}{an optional object of class \code{"crossdist"}
    compatible with the fitted som (see details)}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the 
  error made by the fitted som on the data used to fit it. The
  dissimilarity between a data point and the prototype of its best
  matching unit is computed via the relationa formula. It is also used
  to compute dissimilarities between prototypes. Negative values
  that might occur in this formula are replaced by zero values (a
  warning is generated during this process). 

  When \code{newdata} is specified, the function returns the error of the
  fitted som on the corresponding data. The object must be of class
  \code{"crossdist"} as returned by \code{\link[proxy]{dist}} and must
  contain the dissimilarities between the original data (used to fit the
  SOM) and the new data (for which the error is to be
  computed). 
}
\references{
  Kaski, S. and Lagus, K. (1996) Comparing self-organizing maps, in:
  C. von der Malsburg, W. von Seelen, J. Vorbrüggen, B. Sendhoff (eds.),
  Proceedings of International Conference on Artificial Neural Networks
  (ICANN'96, Bochum, Germany), vol. 1112 of Lecture Notes in Computer
  Science, Springer, pp. 809--814.
}
\author{Fabrice Rossi}
\seealso{\code{\link{error.quantisation}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling and dissimilarity computation
data <- dist(scale(iris[1:4]))

# a medium hexagonal grid
sg <- somgrid(xdim=11,ydim=11,topo="hex")

# random initialisation (leads generally to medium quality results)
som <- batchsom(data,sg,init="random")
print(paste("Quantisation error:",error.quantisation(som)))
print(paste("Kaski and Lagus' error:",error.kaskilagus(som)))

# pca initialisation (leads generally to better results that random initialisation)
som <- batchsom(data,sg,init="pca")
print(paste("Quantisation error:",error.quantisation(som)))
print(paste("Kaski and Lagus' error:",error.kaskilagus(som)))

}
\keyword{cluster}
