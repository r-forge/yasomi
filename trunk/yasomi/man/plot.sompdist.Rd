\name{plot.sompdist}
\alias{plot.sompdist}
\title{Plot distances between prototypes of a fitted Self-Organising Map}
\description{
  Provide different type of visual representation of the distances
  between prototypes of neighbouring units of a fitted Self-Organising Map.
}
\usage{
\method{plot}{sompdist}(x, mode=c("mean","full"), ...)
}

\arguments{
  \item{x}{an object of class \code{"sompdist"} obtained via \code{\link{prototype.distances}}}
  \item{mode}{specifies how distances are aggregated for display (see
  below for details)}
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.somgrid}}}
}
\details{
This function provides a simple U-matrix like visualisation method for
  the distances between prototypes of direct neighbouring units of a
  fitted SOM. The main idea is to use colour coding of the cells of the
  prior structure to represent those distances. There are two modes for
  generating the picture:
  \describe{
    \item{\code{"mean"}}{in this mode, the visualisation grid is
  identical to the grid used as the prior structure. Each cell
  (rectangular or hexagonal, depending on the grid type) is filled with
  a colour chosen to represent the average distance between the
  corresponding prototype and the prototypes of its neighbour units in
  the prior structure.}
  \item{\code{"full"}}{in this mode, the function uses a grid
  approximately four time as large as the original one. The
  visualisation is based on the insertion of additional fake units
  between each unit and its direct neighbours in the original grid. As
  in the \code{"mean"} mode, each original unit displays via its colour
  the average of the distances between its prototype and the
  neighbouring one. In addition, fake units displays (again with a colour
  code) the actual distance between a prototype and its neighbour. In
  the case of an hexagonal grid, this can be done exactly. In the case
  of a rectangular grid, only horizontal and vertical neighbours can be
  represented exactly. Units added in diagonal are shared between two
  pair of prototypes and represent therefore the mean of the two
  corresponding distances (see also \code{\link{distance.grid}} for a
  similar solution to the same problem).}
  }
The additional parameters given to the function can be used to control
  the underlying \code{\link{plot.somgrid}} function, e.g. to change the
  default colour palette (\code{\link{heat.colors}}). 
}
\references{
  Ultsch, A. and Siemon, H. P. (1990) Kohonen's self organizing feature
  maps for exploratory data analysis, in: Proceedings of International
  Neural Network Conference (INNC'90).
}
\author{Fabrice Rossi}
\note{
}
\seealso{See \code{\link{prototype.distances}} to get the distance
  structure, \code{\link{umatrix}} for direct access to this type of
  display and \code{\link{distance.grid}} for possibly smoother plots.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
st <- som.tune(data,sg,som.tunecontrol(sg,,criterion=error.kaskilagus))
som <- st$best.som

# compute the distance
pdist <- prototype.distances(som)

# simple mean based umatrix
plot(pdist)

# more complete display
plot(pdist,mode="full")


}
\keyword{cluster}
\keyword{hplot}
