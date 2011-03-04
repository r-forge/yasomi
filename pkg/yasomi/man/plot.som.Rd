\name{plot.som}
\alias{plot.som}
\title{Plot cell contents of a Self-Organising Map}
\description{
Plot a glyph based rendering of data associated to each cell of the
prior structure of a Self-Organising Map.
Each glyph is positioned thanks to the prior structure. 
}
\usage{
\method{plot}{som}(x, y, mode=c("prototype","data"),
type=c("parallel","stars","barplot"), 
with.grid=FALSE, barplot.align=c("bottom","center"), global.scale=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"som"}}
  \item{y}{a matrix or a data frame of data associated to the som's
    prototypes}
  \item{mode}{the mode of the display, either \code{"prototype"} for a
    single glyph per grid cell or \code{"data"} for multiple glyphs (see
    details)} 
  \item{type}{the type of the plot, chosen from \code{"parallel"},
    \code{"stars"} and \code{"barplot"} (see details)}
  \item{with.grid}{logical, with default to \code{FALSE}. If
    \code{TRUE}, the prior grid is drawn underneath the rest of the
    display}
  \item{barplot.align}{specifies the alignment used for the barplots
    when \code{type="barplot"}}
  \item{global.scale}{logical, with default to \code{FALSE}. If
    \code{TRUE}, the prototypes coordinates are scaled globally, while a
  per dimension scaling is used when \code{global.scale=TRUE}}
  \item{\dots}{additional parameters transmitted to the low level plot
    function (see details)}
}
\details{
  This function uses the prior structure of a fitted SOM to arrange some
  glyphs on a screen, either a single glyph per cell (when
  \code{mode="prototype"}) or an arbitrary number of glyphs per cell,
  including zero (when \code{mode="data"}). The primary usage of this
  function is to display either to the prototypes of a fitted SOM or to
  display the data points associated to each cluster.

  The glyph calculation depends on the \code{mode} parameter as follows:
  \describe{
    \item{\code{mode="prototype"}}{In this case, one glyph is computed
  per cell in the prior structure. When \code{y} is missing, the glyph
  corresponds to the prototype associated to the corresponding
  cluster. In this situation, the function is restricted to numerical
  SOMs (as obtained with \code{\link{batchsom.default}}). When \code{y}
  is specified, it must be a \code{matrix} or a \code{data.frame} with
  numerical content. Each row of \code{y} is used to compute a glyph.}
\item{\code{mode="data"}}{In this case, an arbitrary number of glyphs is
  associated to each cell. When \code{y} is missing, a glyph is computed
  for each original data point assigned to the corresponding
  cluster/cell. As in the previous case, the function is restricted to
  numerical SOMs. When \code{y} is specified, it must be a \code{list}
  with length equal to the number of cells in the prior structure (for
  instance, the result of \code{\link{mapToUnit}}). The
  list must be homogeneous and is restricted to two types of
  content. The first type corresponds to matrix/data.frame: all
  matrices/data.frames must have the same number of columns and each row
  defines a glyph. The second type corresponds to factors: each cell of
  \code{y} must be a factor defined on an unique and common set of levels.
    }}
  The actual glyphs depend on the \code{type} parameter:
  \describe{
    \item{\code{type="parallel"}}{this corresponds to a parallel coordinate
      like representation: each prototype is considered as a
      function mapping its coordinate index to the associated
  value. Functions are displayed as lines. 
  }
  \item{\code{type="stars"}}{in this case, the function uses
    \code{\link{stars}} to draw a star plot for each
    prototype, at the corresponding position in the prior structure}

  \item{\code{type="barplot"}}{this representation uses a barplot for each
  prototype. When \code{barplot.align="bottom"} (the default value),
  the prototypes coordinates are scaled to \eqn{[0,1]}{[0,1]} while they
  are scaled to \eqn{[-1,1]}{[-1,1]} when
  \code{barplot.align="center"}. In both cases, the baseline line of the
  barplots is at 0 which gives two different types of visualisation. In
  the case of \code{mode="data"}, bar plots are restricted to the factor
  valued case for \code{y}.
  }
}
Additional parameters can be transmitted to the underlying low level
  drawing functions. Parameters prefixed by \code{grid.} are used for
  the grid plotting when \code{with.grid=TRUE} (see
  \code{\link{plot.somgrid}}).  
}
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
\seealso{See \code{\link{mapToUnit}} for building \code{y} when \code{mode="data"}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# and display the prototypes as parallel coordinates
plot(som)

# display the prototypes as star plots, with a red grid
plot(som,type="stars",with.grid=TRUE,grid.border="red")

# display all the data as parallel coordinates
plot(som,mode="data",with.grid=TRUE)

# display an additional factor variable with a barplot
plot(som,mapToUnit(som,iris[[5]]),mode="data",type="barplot")

}
\keyword{cluster}
\keyword{hplot}
