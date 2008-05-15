\name{batchsom.dist}
\alias{batchsom.dist}
\title{Fit a Self-Organising Map to dissimilarity data}
\description{
\code{batchsom.dist} is used to fit a Self-Organising Map to
dissimilarity data.
}
\usage{
\method{batchsom}{dist}(data, somgrid, init=c("pca","random"), prototypes,
         assignment = c("single", "heskes"),radii=somradii(somgrid),
         maxiter = 75,kernel=c("gaussian","linear"),normalised,
         cut = 1e-07, verbose = FALSE, keepdata = TRUE,
         lowlevel=fastRelationalsom.lowlevel.R,\dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted represented by a
    dissimilarity matrix (an object of class \code{"dist"} from the
    \code{proxy} package) with all pairwise dissimilarities between the
    observations.}
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{prototypes}{a matrix of initial values for the
    prototypes. If missing chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{init}{the initialisation method (see details)}
  \item{assignment}{the assignment method used to compute the best
    matching unit of an observation during training}
  \item{radii}{the radii of the neighbourhoods used during training. See
  \code{\link{somradii}} for the default values}
  \item{maxiter}{the maximum number of of iteration performed for each
    radius (the minimum is one)}
  \item{kernel}{the kernel to use to compute neighbouring iterations}
  \item{normalised}{switch for normalising the neighbouring
    interactions. Has no influence with the \code{"single"} assignment method}
  \item{cut}{minimal value below wich neighbouring interactions are not
    take into account}
  \item{verbose}{switch for tracing the fitting process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{lowlevel}{the underlying low level implementation of the method
    (don't use this parameter, it's here only during the development process)}
  \item{\dots}{additional arguments to be passed to the initialisation method}
}
\details{
  This function implements the relational Self-Organising Map algorithm
  in which virtual linear combination of the original data are used to
  represent the prototypes. If the initial value of
\code{prototypes} is not provided, it is obtained by a call to
a function specified by the \code{init} parameter. If its value is
\code{"pca"} proprototypes are obtained by a call to 
\code{\link{sominit.pca.dist}} (this is also the case when
\code{init} is not specified), while
\code{\link{sominit.random.dist}} is called when \code{init} is
\code{"random"}. In both case, the additional parameters submitted to
the method are transmitted to the initialisation method. 
}
\value{
  An object of class \code{"som"} and of class \code{"relationalsom"}, a
  list with components including 
  \item{somgrid}{as in the arguments to \code{batchsom}}
  \item{prototypes}{a matrix containing the virtual coordinates of the
    prototypes: each row of the matrix sums to one and can be
    interpreted as the coefficients of a linear combination of the
    original observations.}
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{somradii}} for some reasonable default values for the
  radii, \code{\link{sominit.pca.dist}} and
  \code{\link{sominit.random.dist}} 
  for some control on the initial configuration of the prototypes, \code{\link{som.tune}} for the optimisation of some magic
  parameters (such as the radii), \code{\link{umatrix}} and
  \code{\link{distance.grid}} for visual analysis of the distances
  between the prototypes.}
\examples{
data(iris)
# scaling and dissimilarity computation
data <- dist(scale(iris[1:4]))

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM (random initialisation)
som <- batchsom(data,sg,init="random",method="cluster")

# and display the umatrix
umatrix(som)
}
\keyword{cluster}
