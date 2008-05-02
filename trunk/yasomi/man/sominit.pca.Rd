\name{sominit.pca}
\alias{sominit.pca}
\title{Initialise the prototypes of a SOM with PCA}
\description{
Initialise the prototypes of a Self-Organising Map with Principal Component
Analysis or a related method (e.g., classical multidimensional scaling
for a dissimilarity matrix). 
}
\usage{
sominit.pca(data, somgrid, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    e.g., a matrix or data frame of observations (which should be scaled), or a
    distance matrix}
  \item{somgrid}{a \code{somgrid} object}
  \item{\dots}{additional parameters to be passed to methods}
}
\value{
  A matrix containing appropriate initial prototypes
}
\author{Fabrice Rossi}
\seealso{\code{\link{sominit.pca.default}} and
  \code{\link{sominit.pca.dist}} for implementations,
and \code{\link{sominit.random}} for an alternative.}
\keyword{cluster}

