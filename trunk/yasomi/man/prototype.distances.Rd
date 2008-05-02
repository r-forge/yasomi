\encoding{UTF-8}
\name{prototype.distances}
\alias{prototype.distances}
\alias{as.matrix.sompdist}
\alias{as.dist.sompdist}
\title{Compute distances between neighbouring prototypes in a fitted Self-Organising Map}
\description{
Compute the distances between the prototypes of a fitted
Self-Organising Map (SOM) for neighbouring units in the prior structure.
}
\usage{
prototype.distances(som)
\method{as.matrix}{sompdist}(x, extended = TRUE, \dots)
\method{as.dist}{sompdist}(x, FUN = NULL)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{x}{an object of class \code{"sompdist"}}
  \item{extended}{logical. If \code{TRUE} the distance matrix or object
    build from neighbouring prototypes is extended to cover all pairs of
  prototypes via a graph based approach (see below for details). If
  \code{FALSE} only neighbouring prototypes are considered}
\item{\dots}{not used}
\item{FUN}{not used}
}
\details{
  The \code{as.matrix} method can be used to convert a distance matrix
  (an object of class\code{"sompdist"}) into a classical distance matrix
  as obtained via the \code{as.matrix} applied to objects of class
  \code{"dist"} (see \code{\link{dist}}). \strong{However}, the distance
  matrix obtained this way \strong{is not} the matrix of the Euclidean
  distances between the prototypes of the SOM: only
  the distances between direct neighbours in the grid are
  considered (see \code{\link{as.dist.somnum}} and
  \code{\link{as.dist.relationalsom}} for all pairwise
  distances between prototypes). There are used as edge lengths (or weights) for an
  undirected graph whose nodes are the units and which has edges only
  between direct neighbours. If the \code{as.matrix} method is called
  with \code{extended=TRUE} (the default value), then this graph is used
  to define a distance between any pair of nodes in the grid as the
  length of the shortest path between the considered node (this is the
  unit distance used in \code{\link{error.kaskilagus}}). If the
  method is used with \code{extended=FALSE}, then the missing distances
  (i.e. between non direct neighbours) are unknown (value \code{NA}).

  The \code{as.dist} method can be used to convert a distance matrix
  (an object of class \code{"sompdist"}) into an object of class
  \code{"dist"}. The extension procedure described above is always
  used. 

  Visualisation of the obtained distance matrix can be done directly
  with the \code{\link{plot.sompdist}} method or indirectly via the
  \code{\link{distance.grid}} function.
}
\value{
  \code{prototype.distances} returns an object of
  class\code{"sompdist"}, a list with components: 
  \item{pdist}{a matrix containing the distances between prototypes. It
  has one row for each prototype and 6 or 8 columns, 
  depending on the topology of the prior structure (6 for a hexagonal
  grid and 8 for a rectangular one). When a prototype has not the
  maximum number of neighbours because of its position in the grid,
  the corresponding entries are encoded with \code{NA}.
}
  \item{somgrid}{an object of class \code{"somgrid"}, the prior
    structure used for fitting the \code{som} object}
  \code{as.matrix} returns a square matrix and \code{as.dist} an object
  of class \code{"dist"}.
}
\references{
  Kaski, S. and Lagus, K. (1996) Comparing self-organizing maps, in:
  C. von der Malsburg, W. von Seelen, J. Vorbrüggen, B. Sendhoff (eds.),
  Proceedings of International Conference on Artificial Neural Networks
  (ICANN'96, Bochum, Germany), vol. 1112 of Lecture Notes in Computer
  Science, Springer, pp. 809--814.
}
\author{Fabrice Rossi}
\note{
  This function applies to standard numerical Self-Organising Map but
  also to relational SOM. In the first case, the distances between
  prototypes is the standard Euclidean distance. In the second case, the
  distance is obtained via the relational formula which can lead to
  negative values. When such values are obtained, they are replaced by
  zeros and a warning is generated. 
}
\seealso{See \code{\link{distance.grid}} for converting the prototypes
  distances into a grid of values that can be easily displayed,
  \code{\link{plot.sompdist}} and \code{\link{umatrix}} for simple
  visualisations of those distances and 
  \code{\link{error.kaskilagus}} for an error measure based on
  prototypes distances. See \code{\link{as.dist.somnum}} for all pairwise
  distances between prototypes.} 
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=8,ydim=8,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# compute the distance
pdist <- prototype.distances(som)

# crude visualisation
image(as.matrix(pdist))
}
\keyword{cluster}

