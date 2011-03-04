testGridSizes <- function() {
    nb.grids <- 100
    xd <- ceiling(runif(min=0.1,max=50,nb.grids))
    yd <- ceiling(runif(min=0.1,max=50,nb.grids))
    for(i in seq_along(xd)) {
        sg <- somgrid(xd[i],yd[i],with.dist=FALSE)
        checkEquals(sg$xdim,xd[i])
        checkEquals(sg$ydim,yd[i])
        checkEquals(sg$size,xd[i]*yd[i])
        checkEquals(nrow(sg$pts),xd[i]*yd[i])
        checkEquals(ncol(sg$pts),2)
    }
}
