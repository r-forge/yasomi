\name{error.quantisation.relationalsom}
\alias{error.quantisation.relationalsom}
\title{Quantisation error for a Relational Self-Organising Map}
\description{
Compute the quantisation error made by a fitted Self-Organising Map on
dissimilarity data: this is the mean of the dissimilarity between each
observation and the prototype of its best matching unit.
}
\usage{
\method{error.quantisation}{relationalsom}(som, newdata, \dots)
}
\arguments{
  \item{som}{an object of class \code{"somnum"}}
  \item{newdata}{an optional object of class \code{"crossdist"}
    compatible with the fitted som (see details)}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the quantisation
  error made by the fitted som on the data used to fit it. The
  dissimilarity between a data point and the prototype of its best
  matching unit is computed via the relationa formula. Negative values
  that might occur in this formula are replaced by zero values (a
  warning is generated during this process). 

  When \code{newdata} is specified, the function returns the
  quantisation error of the
  fitted som on the corresponding data. The object must be of class
  \code{"crossdist"} as returned by \code{\link[proxy]{dist}} and must
  contain the dissimilarities between the original data (used to fit the
  SOM) and the new data (for which the quantisation error is to be
  computed). 
}
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus.relationalsom}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling and dissimilarity computation
data <- dist(scale(iris[1:4]))

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

print(error.quantisation(som))

# a larger grid should in general give a lower quantisation error
sg <- somgrid(xdim=12,ydim=12,topo="hex")
som <- batchsom(data,sg)
print(error.quantisation(som))
}
\keyword{cluster}
