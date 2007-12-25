\name{plot.somnum}
\alias{plot.somnum}
\title{Plot a Self-Organising Map}
\description{
Plot a glyph based rendering of the prototypes of a Self-Organising Map.
Each glyph is positioned thanks to the prior structure. 
}
\usage{
\method{plot}{somnum}(x,y,...)
}
\arguments{
  \item{x}{an object of class \code{"somnum"}}
  \item{y}{a matrix or a data frame of data associated to the som's prototypes}
  \item{\dots}{additional parameters transmitted to the low level plot function}
}
%\details{
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
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# and display the prototypes
plot(som)

# display the class mapping
plot(som,iris[[5]])
}
\keyword{cluster}
\keyword{hplot}
