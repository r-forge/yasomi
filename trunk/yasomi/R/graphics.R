grid2lines <- function(som,prototypes) {
    if(missing(prototypes)) {
        if(!inherits(som,"somnum")) {
            warning("'som' is not a standard numerical SOM, results might be meaningless")
        }
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
    if(!inherits(som,"somnum")) {
        stop("'som' is not a standard numerical SOM")
    }
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

hitMap <- function(som,border=NA,with.grid=TRUE,...) {
    args <- list(...)
    if(is.null(args$col)) {
        args$col <- "red"
    }
    if(with.grid) {
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

plot.som <- function(x,y,mode=c("prototype","data"),type=c("parallel","stars","barplot"),with.grid=FALSE,barplot.align=c("bottom","center"),global.scale=FALSE,...) {
    mode <- match.arg(mode)
    type <- match.arg(type)
    barplot.align <- match.arg(barplot.align)
    ## prepare optional arguments
    args <- list(...)
    grid.params <- split.on.prefix("grid",args)
    omar <- par()$mar
    on.exit(par(mar=omar))
    plane.mar <- c(0.5,0.5,4,0.5)
    if(is.null(args$main) | identical(args$main,"")) {
        plane.mar[3] <- 0.5
    }
    par(mar=plane.mar)
    if(x$somgrid$topo=="rectangular") {
        the.width <- 0.40
        the.height <- 0.45
    } else {
        the.width <- 0.4
        the.height <- 0.33
    }
    the.levels <- NULL
    if(missing(y)) {
        if(!inherits(x,"somnum")) {
            stop("'x' is not a standard numerical SOM")
        }
        if(mode=="prototype") {
            y <- x$prototypes
            classif <- 1:nrow(x$prototypes)
        } else {
            if(is.null(x$data)) {
                stop("cannot display data without saved data")
            }
            y <- x$data
            classif <- x$classif
        }
    } else {
        if(mode=="prototype") {
            if(nrow(y) != x$somgrid$size) {
                stop("'y' has not the same dimension has 'x$somgrid'")
            }
        } else {
            if(!is.list(y)) {
                stop("'y' must be a list when mode='data'")
            }
            if(length(y) !=  x$somgrid$size) {
                stop("'y' has not the same length has 'x$somgrid'")
            }
            data.e <- to.matrix(y)
            y <- data.e$m
            classif <- data.e$classif
            the.levels <- data.e$levels
        }
    }
    if(with.grid) {
        if(is.null(grid.params$with.prefix$border)) {
            grid.params$with.prefix$border="gray"
        }
        do.call("plot",c(list(x=x$somgrid),grid.params$with.prefix))
    } else {
        if(length(grid.params$with.prefix)>0) {
            warning("grid.* arguments are discarded")
        }
        plot(x$somgrid,type="n")
    }
    if(type=="stars") {
        if(is.null(grid.params$others$lwd)) {
            grid.params$others$lwd=1
        }
        y.scaled <- flex.scale(y,c(0,0.5),global.scale)
        do.call("stars",c(list(x=y.scaled,locations=x$somgrid$pts[classif,],radius=F,add=T,scale=F),grid.params$others))
    } else if(type=="parallel") { 
        y.scaled <- flex.scale(y,c(-the.height,the.height),global.scale)
        base.x <- seq(from=-the.width,to=the.width,length.out=ncol(y))
        x.dec <- rep(c(base.x,NA),times=nrow(y))
        y.dec <-
            as.vector(t(cbind(y.scaled,rep(NA,times=nrow(y)))))
        x.base <- rep(x$somgrid$pts[classif,1],each=ncol(y)+1)
        y.base <- rep(x$somgrid$pts[classif,2],each=ncol(y)+1)
        do.call("lines",c(list(x=x.dec+x.base,y=y.dec+y.base),grid.params$others))
    } else {
        if(mode=="data" && is.null(the.levels)) {
            stop("'barplot' is not supported in mode='data' for numerical data")
        }
        y.scaled <- flex.scale(y,c(-the.height,the.height),global.scale)
        bar.width <- 0.9*2*the.width/ncol(y)
        bar.space <- 0.1*2*the.width/(ncol(y)-1)
        base.x <-
            seq(from=-the.width,by=(bar.width+bar.space),length.out=ncol(y))
        x.left <- rep(base.x,times=nrow(y))+rep(x$somgrid$pts[,1],each=ncol(y))
        y.bottom <- rep(x$somgrid$pts[,2],each=ncol(y))
        if(barplot.align=="bottom") {
            y.bottom <- y.bottom - the.height
            y.scaled <- y.scaled + the.height
        }
        rect(x.left,y.bottom,x.left+bar.width,y.bottom+as.vector(t(y.scaled)),col=rainbow(ncol(y)))
    } 
}


