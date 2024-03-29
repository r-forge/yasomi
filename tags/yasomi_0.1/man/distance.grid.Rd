\name{distance.grid}
\alias{distance.grid}
\title{Build a surface that represents prototype distances in a fitted Self-Organising Map}
\description{
  Build a regularly sampled surface that represents the distances
  between prototypes of a fitted Self-Organising Map that are direct
  neighbours in the prior structure. The surface is generated via the
  U-matrix principle. 
}
\usage{
distance.grid(sompdist, mode = c("mean", "full"), nxo, nyo)
}

\arguments{
  \item{sompdist}{an object of class \code{"sompdist"} generated by \code{\link{prototype.distances}}}
  \item{mode}{specifies what distances are taken into account to
  generate the surface (see Details) }
  \item{nxo}{horizontal resolution of the generated surface}
  \item{nyo}{vertical resolution of the generated surface}
}
\details{
  There are two modes for generating the surface, \code{"mean"} and
  \code{"full"}. In the \code{"mean"} mode, the function first computes
  the mean distance between each prototype and all of its neighbours in
  the prior structure. Each value is associated to the position of the
  prototype in the two dimensional prior structure. In the
  \code{"full"}, individual distances between neighbouring prototypes
  are used. Each value is associated to the middle of the segment that
  links the positions of corresponding prototypes in the prior
  structure.

  Both prior topologies need some specific rules:
  \describe{
  \item{\code{"rectangular"}}{expect for those on the border, each prototype
  has eight neighbours. Segments linking diagonal neighbours of
  horizontal or vertical neighbours intersect in their middle. As a
  consequence, each corresponding intersection is used by two prototype
  pairs: the associated value is the mean of the distances between the
  involved prototypes. In other terms, only horizontal and vertical
  neighbours have exact distances mapped on the surface.}
\item{\code{"hexagonal"}}{hexagonal prior structures are not grids. The
  function must therefore interpolate the exact surface to evaluate it
  on a regular grid}
}

    When \code{nxo} and/or \code{nyo} is specified, the function first
  compute the original surface and then generate an interpolated surface
  that is evaluated on a regular grid with \code{nxo} vertical lines
  and \code{nyo} horizontal lines (when one of the value is missing, it is
  generated from the underlying \code{somgrid}. The interpolation method
  is based on Radial Basis Function and on a decomposition of the point
  set to avoid inverting large linear system. This a simplified version
  of the algorithm given in Pouderoux et al. article cited in the
  References section.

  The result of the function uses the conventions of standard surface
  plotting functions such as \code{\link{persp}}, \code{\link{contour}},
  \code{\link{filled.contour}}, etc.
}
\value{
  a list with components 
  \item{x}{locations of vertical lines of the grid}
  \item{y}{locations of horizontal lines of the grid}
  \item{z}{the value of the surface at the nodes of the grid}
}
\references{
  Ultsch, A. and Siemon, H. P. (1990) Kohonen's self organizing feature
  maps for exploratory data analysis, in: Proceedings of International
  Neural Network Conference (INNC'90).

 Pouderoux, J., Tobor, I., Gonzato, J.C. and Guitton, P. (2004) Adaptive
 hierarchical RBF interpolation for creating smooth digital elevation
 models, in: Proceedings of the Twelfth ACM International Symposium on
 Advances in Geographical Information System 2004, page 232--240.
}
\author{Fabrice Rossi}
%\note{
%}
\seealso{See \code{\link{prototype.distances}} for computing the
 distances and \code{\link{umatrix}} for a simpler visualisation of the
 same data.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a large hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# fit the SOM
som <- batchsom(sg,data)

# compute prototype distances
pdist <- prototype.distances(som)

# display them via interpolation
pdist.grid <- distance.grid(pdist,mode="mean",nxo=30,nyo=30)

# and contour lines
filled.contour(pdist.grid)
}
\keyword{cluster}

