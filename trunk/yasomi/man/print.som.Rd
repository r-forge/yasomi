\name{print.som}
\alias{print.som}
\alias{summary.som}
\title{Describe a fitted Self-Organising Map}
\description{
  Print and summary methods for Self-Organising Maps. Gives some
  information on the prior structure of the SOM and the way it was
  fitted to a dataset.  
}
\usage{
\method{print}{som}(x, ...)
\method{summary}{som}(object, ...)
}
\arguments{
  \item{x,object}{an object of class \code{"som"}}
  \item{\dots}{not used}
}
\author{Fabrice Rossi}
\seealso{\code{\link{batchsom}}}
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

