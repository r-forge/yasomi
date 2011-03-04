\name{hitMap}
\alias{hitMap}
\title{Clustering Hit Map}
\description{
Plot a hit map for the clustering performed by a Self-Organising Map.
The prior structure is represented by a regular grid of hexagons or
squares.  The surface of each element is proportional to the square root of the number of observations assigned to the corresponding unit in the SOM.
}
\usage{
hitMap(som,border=NA,with.grid=TRUE, \dots)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{border}{the color to draw the border.  The default, \code{NA},
    omits border. Use \code{border = NULL} to use \code{par("fg")}. See
    also \code{\link{polygon}} and \code{\link{rect}}}
  \item{with.grid}{logical. If \code{TRUE}, the hit map includes
    two cells for each unit: an enclosing cell of standard size (see
    \code{\link{plot.somgrid}}) and an internal cell whose surface
    size encodes the cluster size. If \code{FALSE}, the internal cell
    only is drawn} 
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.somgrid}}}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
\references{
  Vesanto, J. (2002) \emph{Data Exploration Process Based on the
    Self--Organizing Map}, PhD thesis, Helsinki University of Technology, Espoo (Finland), Acta Polytechnica Scandinavica, Mathematics and Computing Series No.
  115.
}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{plot.somgrid}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# display the hit map
hitMap(som)

# with another color and borders with default color
hitMap(som,col="blue",border=NULL)

}
\keyword{cluster}
\keyword{hplot}
