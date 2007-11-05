\name{hitMap}
\alias{hitMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering Hit Map}
\description{
Plot a hit map for the clustering performed by a Self-Organising Map.
The prior structure is represented by a regular grid of hexagons or
squares.  The surface of each element is proportional to the square root of the number of observations assigned to the corresponding unit in the SOM.
}
\usage{
hitMap(som,border=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{border}{the color to draw the border.  The default, \code{NA},
    omits border. Use \code{border = NULL} to use \code{par("fg")}. See
    also \code{\link{polygon}} and \code{\link{rect}}}
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.somgrid}}}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
\references{Vesanto's thesis}
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
som <- batchsom(sg,data)

# display the hit map
hitMap(som)

# with another color and borders with default color
hitMap(som,col="blue",border=NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{hplot}
