\name{som.tunecontrol}
\alias{som.tunecontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Control parameters for the som.tune function}
\description{
Creates a list of parameters for the \code{\link{som.tune}} function.
}
\usage{
som.tunecontrol(somgrid, init = "pca", ninit = 1, assignment = "single",
                radii = c(2, 2/3 * somgrid$diam), nradii = 10,
                innernradii = 30, maxiter = 75, annealing = "power",
                kernel = "gaussian", criterion = error.quantisation)
}
\arguments{
  \item{somgrid}{an object of class \code{"somgrid"}}
  \item{init}{prototypes initialization method. Valid values are \code{"pca"}
    and \code{"random"}. The former 
    corresponds to principal component based initialization (see
    \code{\link{sominit.pca}}), while the latter uses randomly selected
  observations as initial values for the prototypes (see \code{\link{sominit.random}})}
  \item{ninit}{number of initial prototype values to test (only relevant
  for \code{init="random"})}
  \item{assignment}{assignment method with valid values \code{"single"} and
    \code{"heskes"} (see \code{\link{batchsom}})}
  \item{radii}{the range of radii to explore, i.e., a vector of length two
    containing a minimal and a maximal value of radii. The default
    minimum radius is 2 (almost purely local k-means like optimization)
    while the maximum is equal to two third of the diameter of the prior
    struture}
  \item{nradii}{number of radii to generate from the range specified in \code{radii}}
  \item{innernradii}{number of radii to use in the annealing scheme
    during the SOM fitting (see \code{\link{batchsom}})}
  \item{maxiter}{maximal number of iteration for each radius during
    fitting (see \code{\link{batchsom}})}
  \item{annealing}{annealing scheme with valid values \code{"power"}
    (exponential like annealing) and \code{"linear"} (linear scheme)}
  \item{kernel}{kernel chosen between \code{"gaussian"} and \code{"linear"}}
  \item{criterion}{an error criterion, i.e., a function that evaluate
    the quality of a fitted som on a dataset}
}
\details{
The parameters \code{init}, \code{assignment}, \code{annealing} and
\code{kernel} can contain a list of values rather than a single
value. In this case, all combinations are compared in
\code{\link{som.tune}}. For instance, if all parameters have default
values expect for \code{kernel=c("gaussian","linear")} then
\code{\link{som.tune}} will fit 10 SOM with a gaussian kernel as well as
10 SOM with a linear kernel and select the best one.

The error criterion is a function which takes a object returned by the
appropriate \code{\link{batchsom}} method (e.g. a \code{"somnum"} object
for standard vector data and a \code{"relationalsom"} object for
dissimilarity data) and returns a numerical value. \code{som.tune}
optimises for small values of this criterion.
}
\value{
  a list containing all the above parameters with one exception: 'radii'
  contains a vector of radii to test rather than the bounds specified in
  the function call. 
}
\author{Fabrice Rossi}
\note{this function is modelled after David Mayer's \link[e1071]{tune.control} function.
}
\seealso{\code{\link{som.tune}}}
\keyword{cluster}

