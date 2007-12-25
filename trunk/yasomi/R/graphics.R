grid2lines <- function(som,prototypes) {
    if(missing(prototypes)) {
        prototypes <- som$prototypes
    } else {
        if(nrow(prototypes)!=som$somgrid$size) {
            stop("the grid size is different from the number of prototypes")
        }
    }

    if(som$somgrid$topo=="rectangular") {
        result=matrix(NA,
        ncol=ncol(prototypes),
        nrow=(som$somgrid$xdim+1)*(som$somgrid$ydim)+(som$somgrid$xdim)*(som$somgrid$ydim+1))
    } else {
        result=matrix(NA,
        ncol=ncol(prototypes),
        nrow=(som$somgrid$xdim+1)*(som$somgrid$ydim)+(som$somgrid$xdim)*(som$somgrid$ydim+1)+(som$somgrid$xdim-1)*(som$somgrid$ydim+1))
    }
    ## first "horizontal" lines
    shift <- (som$somgrid$xdim+1)*som$somgrid$ydim
    targetpos <- 1:(shift-1)
    targetpos <- targetpos[-(seq(from=som$somgrid$xdim+1,to=length(targetpos),by=som$somgrid$xdim+1))]
    result[targetpos,] <- prototypes
    ## then "vertical" lines
    targetpos <- (shift+2):(shift+som$somgrid$xdim*(som$somgrid$ydim+1)+1)
    shift <- targetpos[length(targetpos)]
    targetpos <- targetpos[-(seq(from=som$somgrid$ydim+1,to=length(targetpos),by=som$somgrid$ydim+1))]
    result[targetpos,] <- prototypes[seq(from=1,by=som$somgrid$xdim,length.out=som$somgrid$ydim)+rep(0:(som$somgrid$xdim-1),each=som$somgrid$ydim),]
    ## other "vertical" lines in hexagonal grids
    if(som$somgrid$topo=="hexagonal") {
        targetpos <- (shift+1):(shift+(som$somgrid$xdim-1)*(som$somgrid$ydim+1))
        targetpos <- targetpos[-(seq(from=som$somgrid$ydim+1,to=length(targetpos),by=som$somgrid$ydim+1))]
        result[targetpos,] <- prototypes[seq(from=2,by=som$somgrid$xdim,length.out=som$somgrid$ydim)+(rep(c(0,-1),length.out=som$somgrid$ydim)+rep(0:(som$somgrid$xdim-2),each=som$somgrid$ydim)),]
    }
    result
}

componentPlane <- function(som,dim=1,...) {
    args <- list(...)
    if(is.null(args$main)) {
        main <- as.character(dim)
        if(!is.null(colnames(som$prototypes))[1]) {
            main <- colnames(som$prototypes)[dim]
        }
        args$main <- main
    }
    do.call("plot",c(list(x=som$somgrid,colorValues=som$prototypes[,dim],border=NA),args))
}

hitMap <- function(som,border=NA,with.cells=TRUE,...) {
    args <- list(...)
    if(is.null(args$col)) {
        args$col <- "red"
    }
    if(with.cells) {
        add <- !is.null(args$add) && args$add
        plot(som$somgrid,add=add)
        args$add <- TRUE
    }
    args$border <- border
    sizes <- table(factor(som$classif,levels=1:nrow(som$prototypes)))
    sizes <- sizes/max(sizes)
    do.call("plot",c(list(x=som$somgrid,size=sizes),args))
}

umatrix <- function(som,...) {
    args <- list(...)
    if(is.null(args$border)) {
        args$border <- NA
    }
    pdist <- prototype.distances(som)
    do.call("plot",c(list(x=pdist),args))
    invisible(pdist)
}

## very basic version
plot.somnum <- function(x,y,...) {
    if(missing(y)) {
        stars(x$prototypes,locations=x$somgrid$pts,len=0.5,radius=F,...)
    } else {
        if(is.factor(y)) {
            ymat <- mapFactorToUnit(x,y)
            stars(ymat,locations=x$somgrid$pts,len=0.5,draw.segments=T,col.segments=rainbow(ncol(ymat)),...)
        } else {
            stars(y,locations=x$somgrid$pts,len=0.5,...)
        }
    }
}
