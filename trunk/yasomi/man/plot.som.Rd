\name{plot.som}
\alias{plot.som}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a Self-Organising Map}
\description{
Plot a glyph based rendering of the prototypes of a Self-Organising Map.
Each glyph is positioned thanks to the prior structure. 
}
\usage{
\method{plot}{som}(x,y,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"som"}}
  \item{y}{a matrix or a data frame of data associated to the som's prototypes}
  \item{\dots}{additional parameters transmitted to the low level plot function}
}
%\details{
%}
%\value{
%}
\references{Vesanto's thesis}
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
som <- batchsom(sg,data)

# and display the prototypes
plot(som)

# display the class mapping
plot(som,iris[[5]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{hplot}
