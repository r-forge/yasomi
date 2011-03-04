\name{error.quantisation.kernelsom}
\alias{error.quantisation.kernelsom}
\title{Quantisation error for a Kernel Self-Organising Map}
\description{
Compute the quantisation error made by a Self-Organising Map fitted on
 data mapped in a kernel induced feature space: this is the mean of the
 distance between each observation and the prototype of its best
 matching unit, computed in the kernel space.
}
\usage{
\method{error.quantisation}{kernelsom}(som, newdata, newdatanorms, \dots)
}
\arguments{
  \item{som}{an object of class \code{"somnum"}}
  \item{newdata}{an optional matrix compatible with the fitted som (see
    details)}
  \item{newdatanorms}{an optional vector containing the squared norm in
    kernel space of the new data (see details)}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the quantisation
  error made by the fitted som on the data used to fit it. 

  When \code{newdata} is specified, the function returns the
  quantisation error of the fitted som on the corresponding data.
  The \code{newdata} matrix must contain the values of the kernel used
  to fit the SOM evaluated between the original data and the new
  data. More precisely, \code{newdata[i,j]} contains the value of
  \eqn{K(x_i,nx_j)}, where \eqn{x_i} is the i-th original data point and
  \eqn{nx_j} is the j-th new data point. In addition, when
  \code{newdata} is specified, then the function needs the value of the
  squared norm in kernel space of the new data, that is a vector
  \code{newdatanorms} such that \code{newdatanorms[j]} contains
  \eqn{K(nx_j,nx_j)}. 

  The distance between a data point and the prototype of its best
  matching unit is computed in the kernel space. If the kernel used to
  fit the SOM is not postive, negative values that might occur during
  the distance computation are replaced by zero values (a warning is
  generated during this process).
}
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus.kernelsom}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling and kernel computation
data <- as.kernelmatrix(tcrossprod(scale(iris[1:4])))

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
