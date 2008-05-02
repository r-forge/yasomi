\name{print.som}
\alias{print.som}
\alias{print.somnum}
\alias{print.relationalsom}
\alias{summary.som}
\alias{summary.somnum}
\alias{summary.relationalsom}
\title{Describe a fitted Self-Organising Map}
\description{
  Print and summary methods for Self-Organising Maps. Gives some
  information on the prior structure of the SOM and the way it was
  fitted to a dataset.  
}
\usage{
\method{print}{som}(x, \dots)
\method{print}{somnum}(x, \dots)
\method{print}{relationalsom}(x, \dots)
\method{summary}{som}(object, \dots)
\method{summary}{somnum}(object, \dots)
\method{summary}{relationalsom}(object, \dots)
}
\arguments{
  \item{x,object}{an object of class \code{"som"}}
  \item{\dots}{not used}
}
\author{Fabrice Rossi}
\seealso{\code{\link{batchsom}} for fitting a SOM to some data,
  \code{\link{plot.som}} for a graphical representation.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=5,ydim=5,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# summary
print(som)
}
\keyword{cluster}

