\name{predict.kernelsom}
\alias{predict.kernelsom}
\title{Map new data to a fitted Kernel Self-Organising Map}
\description{
Map new data to a fitted Self-Organising Map, i.e., compute the best
matching unit for each given observation. 
}
\usage{
\method{predict}{kernelsom}(object, newdata, newdatanorms,
        with.secondwinner=FALSE,\dots)
}
\arguments{
  \item{object}{an object of class \code{"kernelsom"}}
  \item{newdata}{a matrix compatible with the fitted som (see details)}
  \item{newdatanorms}{a vector containing the squared norm in kernel space of
    the new data (see details)}
  \item{with.secondwinner}{switch to specify whether the second best
    matching unit should be computed and returned}
  \item{\dots}{not used}
}
\details{
  The \code{newdata} matrix must contain the values of the kernel used
  to fit the SOM evaluated between the original data and the new
  data. More precisely, \code{newdata[i,j]} contains the value of
  \eqn{K(x_i,nx_j)}, where \eqn{x_i} is the i-th original data point and
  \eqn{nx_j} is the j-th new data point.

  While this is not needed to compute the best matching unit, this
  function requires in addition the value of the squared norm in kernel
  space of the new data, that is a vector \code{newdatanorms} such that
  \code{newdatanorms[j]} contains \eqn{K(nx_j,nx_j)}. 
}
\value{
  A list with components
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{error}{the quantisation error of the observations by the
    prototypes of this Self-Organising Map (see
    \code{\link{error.quantisation.kernelsom}} for a definition of
    this value)}
  \item{rdist}{a matrix of squared dissimilarities between the new data
    and the prototypes (some values might be negative if the underlying
    kernel is not positive)}
  \item{winners}{if \code{with.secondwinner=TRUE} a two column matrix
    with a line per observation containing the index of the best
    matching unit and of the second best matching unit (in this order)}
}
\author{Fabrice Rossi}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# use only part of the data
train <- sample(1:nrow(data),100)
data.train <- data[train,]
K <- as.kernelmatrix(tcrossprod(data.train))

# a small hexagonal grid
sg <- somgrid(xdim=5,ydim=5,topo="hex")

# fit the SOM
som <- batchsom(K,sg)

# map remaining data
results <- predict(som,tcrossprod(data.train,data[-train,]),
                   diag(tcrossprod(data[-train,])))
print(paste("Learning quantisation error:",error.quantisation(som)))
print(paste("Test quantisation error:",results$error))
}
\keyword{cluster}
