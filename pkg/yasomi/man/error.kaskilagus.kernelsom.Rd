\encoding{UTF-8}
\name{error.kaskilagus.kernelsom}
\alias{error.kaskilagus.kernelsom}
\title{Kaski and Lagus' error measure for Kernel Self-Organising Maps}
\description{
Compute an error measure of a fitted Kernel Self-Organising Maps defined by
S. Kaski and K. Lagus. All the calculation are performed in the kernel
induced feature space.  
}
\usage{
\method{error.kaskilagus}{kernelsom}(som, newdata, newdatanorms, \dots)
}
\arguments{
  \item{som}{an object of class \code{"kernelsom"}}
  \item{newdata}{an optional matrix compatible with the fitted som (see
    details)}
  \item{newdatanorms}{an optional vector containing the squared norm in
    kernel space of the new data (see details)}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the 
  error made by the fitted som on the data used to fit it. 

  When \code{newdata} is specified, the function returns the
  error of the fitted som on the corresponding data.
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
# scaling and kernel computation
data <- as.kernelmatrix(tcrossprod(scale(iris[1:4])))

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
