\name{sominit.pca.dist}
\alias{sominit.pca.dist}
\title{Initialise the prototypes of a dissimilarity SOM with classical
  MDS}
\description{
Initialise the prototypes of a dissimilarity Self-Organising Map with a
classical multidimensional scaling. The prototypes are regulary
positioned (according to the prior structure) in the two dimensional
space representation of the data according to a classical MDS. If the
dissimilarity matrix contains squared Euclidean distances between vector
data, then this method is equivalent to
\code{\link{sominit.pca.default}}.
}
\usage{
\method{sominit.pca}{dist}(data, somgrid, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted, i.e., a
    distance matrix}
  \item{somgrid}{a \code{somgrid} object}
  \item{\dots}{not used}
}
\details{
  The method first conducts an eigen analysis of the doubly centred
  squared dissimilarity matrix. It generates a regular grid in
  the two dimensional representation space spanned by the first two eigen
  vectors and returns the coordinates of the nodes of the grid as linear
  combinations of the original points. 
}
\value{
  A list with the following components
  \item{prototype}{a matrix containing appropriate initial prototypes}
  \item{D}{the matrix of squared dissimilarities}
  \item{D.c}{the doubly centered matrix obtained from D}
  \item{D.eigen}{the result of the eigen analysis of D.c}
}
\author{Fabrice Rossi}
\seealso{\code{\link{sominit.random}} and \code{\link{somgrid}}}
\examples{
data(iris)

X <- as.matrix(unique(iris[1:2]))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit.pca(dist(X),sg)$prototypes

proto.mapped <- proto\%*\%X

plot(X,pch="+",col="red",xlim=range(X[,1],proto.mapped[,1]),
     ylim=range(X[,2],proto.mapped[,2]))

points(proto.mapped,pch=20)
}
\keyword{cluster}

