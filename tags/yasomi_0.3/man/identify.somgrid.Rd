\name{identify.somgrid}
\Rdversion{1.1}
\alias{identify.somgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify a grid cell in a grid display
}
\description{
This function reads the pointer position when the (first) mouse button
is pressed. The closest grid cell is identified and highlighted (if it
was not already selected). The indices of the identified cells are
returned as the value of the call.
}
\usage{
\method{identify}{somgrid}(x, n=x$size, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"somgrid"}}
  \item{n}{the maximum number of cells to be identified}
  \item{\dots}{not used}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An integer vector containing the indices of the identified cells.  
}
%\references{

%}
\author{Fabrice Rossi}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{identify}} for the generic function. 
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

