somgrid <- function(xdim,ydim,topo=c("rectangular", "hexagonal"),with.dist=TRUE) {
    topo <- match.arg(topo)
    if(xdim == 1) {
        tmp <- xdim
        xdim <- ydim
        ydim <- tmp
    }
    x <- seq(from=1,by=1,length.out=xdim)
    if(topo=="hexagonal" && ydim > 1) {
        y <- rev(seq(from=1,by=sqrt(3)/2,length.out=ydim))
    } else {
        y <- rev(seq(from=1,by=1,length.out=ydim))
    }
    pts <- as.matrix(expand.grid(x = x, y = y))
    if(topo=="hexagonal" && ydim > 1 ) {
        pts[,1] <- pts[,1]+rep(c(0,0.5),each=xdim,length.out=nrow(pts))
    }
    if(topo=="rectangular") {
        diam <- sqrt(xdim^2+ydim^2)
    } else {
        diam <- sqrt(sum((pts[1,]-pts[nrow(pts),])^2))
    }
    res <- list(pts = pts, xdim = xdim, ydim = ydim, topo = topo,
                size = xdim*ydim, diam = diam)
    if(with.dist) {
        res$dist <- as.matrix(dist(pts,method="Euclidean"),diag=0)
    }
    class(res) <- "somgrid"
    res
}

summary.somgrid <- function(object,...) {
    print(object)
}

print.somgrid <- function(x,...) {
    cat("\nPrior structure for a Self-Organising Map\n\n")
    cat("Parameters:\n")
    cat("    topology:",x$topo,"\n")
    cat("       units:",x$size,"\n")
    cat(" x dimension:",x$xdim,"\n")
    cat(" y dimension:",x$ydim,"\n")
    cat("    diameter:",x$diam,"\n")
    cat("\n")
    invisible(x)
}

plot.somgrid <- function(x,col=NA,size=NA,add=FALSE,border=NULL,
                         colorValues=NA,sizeValues=NA,colormap=50,
                         withkey=FALSE,color.palette=heat.colors,...) {
    ## some additional error handling is needed
    args <- list(...)
    if(!add) {
        omar <- par()$mar
        on.exit(par(mar=omar))
        plane.mar <- c(0.5,0.5,4,0.5)
        if(is.null(args$main) | identical(args$main,"")) {
            plane.mar[3] <- 0.5
        }
        par(mar=plane.mar)
    }
    if(withkey & add) {
        stop("'withkey=TRUE' is incompatible with 'add=TRUE'")
    }
    ## take care of border
    if(length(border)>0) {
        border <- rep(border,x$size)
    }
    filter <- rep(TRUE,x$size)
    if(!identical(colorValues,NA)) {
        ## remove cells with NA value
        filter <- !is.na(colorValues)
        ## let's compute values to display
        limits <- range(colorValues,na.rm=TRUE)
        if(is.numeric(colormap)) {
            ## minimum two colors
            nbcolors <- max(2,as.integer(colormap))
            colormap <- color.palette(nbcolors)
        } else {
            nbcolors <- length(colormap)
        }
        breaks <- seq(from=limits[1],to=limits[2],length.out=nbcolors+1)
        code <- cut(colorValues,breaks=breaks,labels=FALSE,include.lowest=TRUE)
        col <- colormap[code]
        if(withkey) {
            def.par <- par(no.readonly = TRUE)
            on.exit(par(def.par))
            layout(matrix(c(2,1),nrow=1),widths=c(0.85,0.15))
            key.mar <- c(1,0.25,max(plane.mar[3],1),3)
            par(mar=key.mar)
            plot.new()
            plot.window(xlim=c(0,1),ylim=limits,yaxs="r")
            rect(0,breaks[-length(breaks)],1,breaks[-1],col=colormap)
            axis(4)
            plane.mar <- plane.mar
            plane.mar[4] <- 0.25
            par(mar=plane.mar)
        }
    } else {
        if(withkey) {
            warning("no key to generate")
        }
        ## take care of colors
        if(!identical(col,NA)) {
            col <- rep(col,length.out=x$size)
        }
    }
    if(!identical(sizeValues,NA)) {
        ## let's compute values to display
        limits <- range(sizeValues,na.rm=TRUE)
        if(limits[1]!=limits[2]) {
            size <- sqrt((sizeValues-limits[1])/(limits[2]-limits[1]))
            ## remove cells with NA value
            filter <- filter & !is.na(sizeValues)
        } else {
            size=NA
        }
    }
    if(x$topo=="rectangular") {
        if(identical(size,NA)) {
            basesize <- 0.5
        } else {
            basesize <- rep(0,x$size)
            basesize[!is.na(size)] <- size[!is.na(size)]/2
        }
        xleft <- (x$pts[,1]-basesize)[filter]
        xright <- (x$pts[,1]+basesize)[filter]
        ybottom <- (x$pts[,2]-basesize)[filter]
        ytop <- (x$pts[,2]+basesize)[filter]
        if(!add) {
            plot(NA,type="n",xlim=range(xleft,xright),
                 ylim=range(ybottom,ytop),xlab="",ylab="",axes=FALSE,...)
        }
        rect(xleft,ybottom,xright,ytop,col=col[filter],border=border[filter])
    } else {
        if(identical(size,NA)) {
            edge <- 1/sqrt(3)
            dec <- 0.5
            hX <- rep(c(dec,   0,   -dec,  -dec,   0,    dec,NA),
                      times=x$size)
            hY <- rep(c(edge/2,edge,edge/2,-edge/2,-edge,-edge/2,NA),
                      times=x$size)
        } else {
            basesize <- rep(0,x$size)
            basesize[!is.na(size)] <- size[!is.na(size)]
            edge <- basesize/sqrt(3)
            dec <- basesize/2
            hX <- rep(dec,each=7)*rep(c(1,0,-1,-1,0,1,NA),times=x$size)
            hY <- rep(edge,each=7)*rep(c(1/2,1,1/2,-1/2,-1,-1/2,NA),times=x$size)
        }
        xpos <- (rep(x$pts[,1],each=7)+hX)[rep(filter,each=7)]
        ypos <- (rep(x$pts[,2],each=7)+hY)[rep(filter,each=7)]
        if(!add) {
            plot(NA,type="n",xlim=range(xpos,na.rm=TRUE),ylim=range(ypos,na.rm=TRUE),xlab="",ylab="",axes=FALSE,...)
        }
        polygon(xpos,ypos,col=col[filter],border=border[filter])
    }
}

radius.exp <- function(min,max,steps) {
    max*(min/max)^(seq(0,1,length.out=steps))
}

radius.lin <- function(min,max,steps) {
    seq(max,min,length.out=steps)
}

somradii <- function(somgrid,min,max,nb,annealing=c("power","linear"),
                     kernel=c("gaussian","linear")) {
    annealing <- match.arg(annealing,c("power","linear"))
    kernel <- match.arg(kernel,c("gaussian","linear"))
    if(missing(max)) {
        max <- 2/3*somgrid$diam+1
    }
    if(missing(min)) {
        min <- switch(kernel,"gaussian"=0.5,"linear"=1)
    }
    if(max<=min) {
        stop("max must be larger than min")
    }
    if(min<=0) {
        stop("min must be positive")
    }
    if(missing(nb)) {
        nb <- max(2,ceiling(2*max))
    }
    switch(annealing,
           "linear"=radius.lin(min,max,nb),
           "power"=radius.exp(min,max,nb))
}

