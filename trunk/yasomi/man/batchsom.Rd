\name{batchsom}
\alias{batchsom}
\title{Fit a Self-Organising Map to some data}
\description{
\code{batchsom} is used to fit a Self-Organising Map to some data, using a
batch algorithm.
}
\usage{
batchsom(somgrid, data, prototypes=somPCAInit(somgrid,data),
         assignment = c("single", "heskes"),radii,nbRadii=30,
         maxiter = 75,kernel=c("gaussian","linear"),normalised,
         cut = 1e-07, verbose = FALSE)
}
\arguments{
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{data}{a matrix or a data frame of observations}
  \item{prototypes}{a matrix of initial values for the
    prototypes. Defaults to the PCA based initialization provided by \code{\link{somPCAInit}}}
  \item{assignment}{the assignment method used to compute the best
    matching unit of an observation during training}
  \item{radii}{the radii of the neighbourhoods used during training. If
    missing, \code{nbRadii} radii are computed starting from a radius of 2/3
    of the diameter of the prior strucure and ending with a radius small
    enough to prevent any influence of a node on its neighbours}
  \item{nbRadii}{the number of radii to generate if 'radii' is missing}
  \item{maxiter}{the maximum number of of iteration performed for each
    radius (the minimum is one)}
  \item{kernel}{the kernel to use to compute neighbouring iterations}
  \item{normalised}{switch for normalising the neighbouring
    interactions. Has no influence with the \code{"single"} assignment method}
  \item{cut}{minimal value below wich neighbouring interactions are not
    take into account}
  \item{verbose}{switch for tracing the fitting process}
}
\details{
It is recommended to use the \code{"single"} assignment method which seems to
perform better in general. 
}
\value{
  An object of class \code{"som"}, a list with components
  \item{somgrid}{as in the arguments to \code{batchsom}}
  \item{prototypes}{a matrix containing the coordinates of the
    prototypes (also called the code book vectors)}
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{som.tune}} for the optimisation of some magic
  parameters (such as the radii). }
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(sg,data)

# and display the prototypes
plot(som)
}
\keyword{cluster}
