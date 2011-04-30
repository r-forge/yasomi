\name{batchsom.dist}
\alias{batchsom.dist}
\title{Fit a Self-Organising Map to dissimilarity data}
\description{
\code{batchsom.dist} is used to fit a Self-Organising Map to
dissimilarity data.
}
\usage{
\method{batchsom}{dist}(data, somgrid, init=c("pca","random"),
                prototypes,weights,
                mode = c("continuous","stepwise"), min.radius, max.radius,
                steps, decrease = c("power", "linear"), max.iter,
                kernel = c("gaussian", "linear", "zeroone"), normalised,
                assignment = c("single", "heskes"), cut = 1e-07,
                verbose = FALSE, keepdata = TRUE, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted represented by a
    dissimilarity matrix (an object of class \code{"dist"} from the
    \code{proxy} package) with all pairwise dissimilarities between the
    observations.}
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{prototypes}{a matrix of initial values for the
    prototypes. It contains linear coefficients that describe the
    prototypes as virtual linear combination of the initial data
    points. It has therefore one row for each prototype (as specified by
    the prior structure \code{somgrid}) and has many columns as data
    points. If missing chosen via 
    the method specified by the \code{init} parameter (see details)}
  \item{init}{the initialisation method (see details)}
  \item{weights}{optional weights for the data points}
  \item{mode}{annealing mode:
    \describe{
      \item{\code{"continuous"} (default)}{this is the standard
	annealing strategy for SOM: the influence of neighbours changes
	at each epoch of the algorithm, from \code{max.radius} to
	\code{min.radius} in exactly \code{step} steps.}
      \item{\code{"stepwise"}}{in this strategy, the algorithm performs
	several epochs (a maximum of \code{max.iter}) for each of the
	\code{step} radii (from \code{max.radius} to
	\code{min.radius}). The algorithm changes the neighbours
	influence only when the classification remains stable from one
	epoch to another. The \code{max.iter} parameter provides a
	safeguard against cycling behaviours.}
    }
  }
  \item{min.radius}{the minimum neighbourhood influence radius. If
    missing, the value might depend on the one of \code{kernel} but ensures
    in practice a local learning only (currently defaults to one)}
  \item{max.radius}{the maximal neighbourhood influence radius. If missing
    two third of the prior structure diameter plus one}
  \item{steps}{the number of radii to use during annealing}
  \item{decrease}{the radii generating formula (\code{"power"} or
    \code{"linear"}), i.e., the way the \code{steps} radii are generated
    from the extremal values given by \code{min.radius} and
    \code{max.radius}}
  \item{max.iter}{maximal number of epochs for one radius in the
    \code{"stepwise"} annealing mode (defaults to 75)}
  \item{kernel}{the kernel used to transform distances in the prior
    structure into influence coefficients}
  \item{normalised}{switch for normalising the neighbouring
    interactions. Has no influence with the \code{"single"} assignment method}
  \item{assignment}{the assignment method used to compute the best
    matching unit (BMU) of an observation during training:
    \describe{
      \item{\code{"single"} (default)}{this is the standard BMU
	calculation approach in which the best unit for an observation
	is the one of the closest prototype of this observation }
      \item{\code{"heskes"}}{Tom Heskes' variant for the BMU in which a
	weighted fit of all the prototypes to an observation is used to
	compute the best unit. The rationale is that the BMU's prototype
	and its neighbouring units' prototypes must be close to the
	observation.}
    }
  }
  \item{cut}{minimal value below wich neighbouring interactions are not
    take into account}
  \item{verbose}{switch for tracing the fitting process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
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
  \item{control}{a list containing all the parameters used to fit the SOM}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
  \item{weights}{the weights of the data points if the function is called with
    \code{keepdata = TRUE} and if the \code{weights} is given}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{sominit.pca.dist}} and
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
