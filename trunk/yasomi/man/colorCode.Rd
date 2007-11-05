\name{colorCode}
\alias{colorCode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data coloring via a Self-Organising Map}
\description{
\code{colorcod}e maps observations to integers in an organised way:
close observations should be assigned close integers (and vice versa).
The coding is done by fitting a one dimensional Self-Organising Map to
the data.  The main use of this method is to provide color based linking
between views. 
}
\usage{
colorCode(data, nbcolor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix or a data frame of observations}
  \item{nbcolor}{the number of colors in the resulting mapping (i.e., the output will belong to [1,\code{nbcolor}]}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
\references{Vesanto's thesis}
\author{Fabrice Rossi}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# code calculation
code <- colorCode(data,50)

# display the result
pairs(data,bg=rainbow(50)[code],pch=21)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{dplot}
