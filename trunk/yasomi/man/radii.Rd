\name{somradii}
\alias{somradii}
\title{Generate radii for SOM annealing}
\description{
This function generates a vector of decreasing radius values to be used
by \code{\link{batchsom}}. Loosely speaking, the radius is used by
\code{batchom} as a measure of the influence of a unit on another in the
SOM. More precisely, if the radius is k, only units at a distance
strictly lower than k have an influence on another unit. 
}
\usage{
somradii(somgrid, min, max, nb, annealing = c("power", "linear"), kernel = c("gaussian", "linear"))
}
\arguments{
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}}}
  \item{min}{the minimum radius of the values to be produced. If
    missing, the value depends on the one of \code{kernel} but ensures
    in practice a local learning only (see details)}
  \item{max}{the maximal radius of the values to be produced. If missing
  two third of the prior structure diameter plus one}
  \item{nb}{the number of radii to generate. If missing twice the value
    of max (with a minimum of two radii)}
  \item{annealing}{the annealing formula (\code{"power"} or
    \code{"linear"}), i.e., the way the \code{nb} radii are generated
    from the extremal values given by \code{min} and \code{max}}
  \item{kernel}{the kernel that will be used by the \code{batchsom}
    call. Used only when \code{min} is not specified.}
}
\details{
The function tries to use some reasonable default values for missing
parameters. In particular, it is generally considered that a SOM fitting
should end by a k-means like phase, that is by a final optimisation in
which the prior structure of the SOM is not used anymore. This is
obtained in \code{\link{batchsom}} by setting the radius to a value
small enough to discard the neighbour influence. This value is one for
the \code{"linear"} kernel and 0.5 (or lower) for the \code{"gaussian"}
kernel.

More generally, a radius of k ensures that a
unit will take into account only units at a distance strictly lower
than k from itself. The interpolation between a medium range influence
(here one plus two third of the diameter of the prior structure) and a
local influence only can be done in a linear way
(\code{annealing="linear"}), i.e with equally spaced radii or in an
exponential way (\code{annealing="power"}). In this latter case, the
difference between two consecutive radii decreases during the annealing
process to favor a longer "cooling" time in the final local phase of the
fitting process. 
}
\value{
  A vector of \code{nb} decreasing values for the radius parameter,
  starting from \code{max} and going down to \code{min}.
}
%\references{}
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%}
\seealso{\code{\link{som.tune}} for an optimisation of some of the
  parameters that define the radii list}
%\examples{
%}
\keyword{cluster}
