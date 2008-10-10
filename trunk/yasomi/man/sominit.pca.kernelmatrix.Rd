\name{sominit.pca.kernelmatrix}
\alias{sominit.pca.kernelmatrix}
\title{Initialise the prototypes of a kernel SOM with kernel PCA}
\description{
Initialise the prototypes of a kernel Self-Organising Map with a
kernel Principal Component Analysis. The prototypes are regulary
positioned (according to the prior structure) in the hyperplane
spanned in the kernel induced feature space by the two first principal
components of the data mapped in this space.  If the
kernel matrix contains standard Euclidean inner products vector
data, then this method is equivalent to
\code{\link{sominit.pca.default}}.
}
\usage{
\method{sominit.pca}{kernelmatrix}(data, somgrid, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted, i.e., a object
    of class \code{"kernelmatrix"}}
  \item{somgrid}{a \code{somgrid} object}
  \item{\dots}{not used}
}
\details{
  The method first conducts an eigen analysis of the doubly centred
  kernel matrix. It generates a regular grid in
  the two dimensional representation space spanned by the first two eigen
  vectors and returns the coordinates of the nodes of the grid as linear
  combinations of the original points. 
}
\value{
  A list with the following components
  \item{prototype}{a matrix containing appropriate initial prototypes}
  \item{D.c}{the doubly centered kernel matrix}
  \item{D.eigen}{the result of the eigen analysis of D.c}
  \item{sdev}{standard deviations on each principal axis in kernel space}
}
%\references{Add reference to Kernel PCA}
\author{Fabrice Rossi}
\seealso{\code{\link{sominit.random}} and \code{\link{somgrid}}}
\examples{
data(iris)

X <- scale(as.matrix(unique(iris[1:2])))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit.pca(as.kernelmatrix(tcrossprod(X)),sg)$prototypes

proto.mapped <- proto\%*\%X

plot(X,pch="+",col="red",xlim=range(X[,1],proto.mapped[,1]),
     ylim=range(X[,2],proto.mapped[,2]))

points(proto.mapped,pch=20)
}
\keyword{cluster}

