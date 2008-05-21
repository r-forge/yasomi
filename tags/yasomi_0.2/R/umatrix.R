protoDist.somnum <- function(som,i,j,k,l) {
    from <- som$prototypes[i+(j-1)*som$somgrid$xdim,]
    to <- som$prototypes[k+(l-1)*som$somgrid$xdim,]
    sqrt(sum((from-to)^2))
}

protoDist.relationalsom <- function(som,i,j,k,l) {
    from <- i+(j-1)*som$somgrid$xdim
    to <- k+(l-1)*som$somgrid$xdim
    relational.sqrt(c(som$prototypes[from,]%*%som$Dalpha[,to]-som$nf[from]-som$nf[to]))
}

prototype.distances <- function(som) {
    sg <- som$somgrid
    if(sg$topo=="rectangular") {
        ## column major
        distances <- matrix(NA,nrow=sg$size,ncol=8)
        pos <- 1
        for(j in 1:sg$ydim) {
            for(i in 1:sg$xdim) {
                ## right
                if(i<sg$xdim) {
                    distances[pos,1] <- protoDist(som,i,j,i+1,j)
                }
                ## lower right
                if(j<sg$ydim & i<sg$xdim) {
                    distances[pos,2] <- protoDist(som,i,j,i+1,j+1)
                }
                ## below neighbor
                if(j<sg$ydim) {
                    distances[pos,3] <- protoDist(som,i,j,i,j+1)
                }
                ## lower left
                if(j<sg$ydim & i>1) {
                    distances[pos,4] <- protoDist(som,i,j,i-1,j+1)
                }
                ## left
                if(i>1) {
                    distances[pos,5] <- distances[pos-1,1]
                }
                ## upper left
                if(i>1 & j>1) {
                    distances[pos,6] <- distances[pos-1-sg$xdim,2]
                }
                ## above
                if(j>1) {
                    distances[pos,7] <- distances[pos-sg$xdim,3]
                }
                ## upper right
                if(j>1 & i<sg$xdim) {
                    distances[pos,8] <- distances[pos+1-sg$xdim,4]
                }
                pos <- pos + 1
            }
        }
    } else {
        ## column major
        distances <- matrix(NA,nrow=sg$size,ncol=6)
        pos <- 1
        for(j in 1:sg$ydim) {
            for(i in 1:sg$xdim) {
                ## right
                if(i<sg$xdim) {
                    distances[pos,1] <- protoDist(som,i,j,i+1,j)
                }
                if(j%%2==1) {
                    if(j<sg$ydim) {
                        ## below right
                        distances[pos,2] <- protoDist(som,i,j,i,j+1)
                        ## below left
                        if(i>1) {
                            distances[pos,3] <- protoDist(som,i,j,i-1,j+1)
                        }
                    }
                    if(j>1) {
                        ## upper right
                        distances[pos,6] <- distances[pos-sg$xdim,3]
                        ## upper left
                        if(i>1) {
                            distances[pos,5] <- distances[pos-sg$xdim-1,2]
                        }
                    }
                } else {
                    if(j<sg$ydim) {
                        ## below left
                        distances[pos,3] <- protoDist(som,i,j,i,j+1)
                        ## below right
                        if(i<sg$xdim) {
                            distances[pos,2] <- protoDist(som,i,j,i+1,j+1)
                        }
                    }
                    if(j>1) {
                        ## upper left
                        distances[pos,5] <- distances[pos-sg$xdim,2]
                        ## upper right
                        if(i<sg$xdim) {
                            distances[pos,6] <- distances[pos-sg$xdim+1,3]
                        }
                    }
                }
                ## left
                if(i>1) {
                    distances[pos,4] <- distances[pos-1,1]
                }
                pos <- pos + 1
            }
        }
    }
    res <- list(pdist=distances,somgrid=sg)
    class(res) <- "sompdist"
    res
}

choose.depth <- function(size) {
    if(size>40) {
        ceiling(log(size/10,base=4))
    } else {
        1
    }
}

distance.grid <- function(sompdist,mode=c("mean","full"),nxo,nyo,outer=FALSE) {
    pdist <- sompdist$pdist
    sg <- sompdist$somgrid
    mode <- match.arg(mode)
    means <- rowMeans(pdist,na.rm=T)
    if(sg$topo=="rectangular") {
        if(mode=="mean") {
            reorder <-
                rep(1:sg$xdim,sg$ydim)+rep((sg$ydim-1):0,each=sg$xdim)*sg$xdim
            distances <- matrix(means[reorder],ncol=sg$ydim)
            x <- seq(from=1,by=1,length.out=sg$xdim)
            y <- seq(from=1,by=1,length.out=sg$ydim)
            radius <- 1.1
        } else {
            distances <- rep(NA,(2*sg$xdim-1)*(2*sg$ydim-1))
            pindex <- 1:sg$size
            ## first fill mean values
            vshift <- rep(seq(from=(2*sg$xdim-1)*(2*sg$ydim-2),
                              by=-2*(2*sg$xdim-1),length.out=sg$ydim),
                          each=sg$xdim)
            index.means <-
                rep(seq(1,by=2,length.out=sg$xdim),sg$ydim) + vshift
            distances[index.means] <- means
            ## then horizontal values
            vshift <- rep(seq(from=(2*sg$xdim-1)*(2*sg$ydim-2),
                              by=-2*(2*sg$xdim-1),length.out=sg$ydim),
                          each=sg$xdim-1)
            index.h <- rep(seq(2,by=2,length.out=sg$xdim-1),sg$ydim) + vshift
            distances[index.h] <-
                pdist[pindex[-seq(sg$xdim,by=sg$xdim,length.out=sg$ydim)],1]
            ## then vertical values
            vshift <- rep(seq(from=(2*sg$xdim-1)*(2*sg$ydim-3),
                              by=-2*(2*sg$xdim-1),length.out=sg$ydim-1),
                          each=sg$xdim)
            index.v <- rep(seq(1,by=2,length.out=sg$xdim),sg$ydim-1) + vshift
            distances[index.v] <- pdist[1:(sg$xdim*(sg$ydim-1)),3]
            ## then diagonal values
            vshift <- rep(seq(from=(2*sg$xdim-1)*(2*sg$ydim-3),
                              by=-2*(2*sg$xdim-1),length.out=sg$ydim-1),
                          each=sg$xdim-1)
            index.d <- rep(seq(2,by=2,length.out=sg$xdim-1),sg$ydim-1) + vshift
            pindex <- 1:(sg$xdim*(sg$ydim-1))
            distances[index.d] <- 0.5*(
                pdist[pindex[-seq(sg$xdim,by=sg$xdim,length.out=sg$ydim-1)],2]+
                    pdist[pindex[-seq(1,by=sg$xdim,length.out=sg$ydim-1)],4])
            distances <- matrix(distances,ncol=2*sg$ydim-1)
            x <- seq(from=1,to=sg$xdim,length.out=2*sg$xdim-1)
            y <- seq(from=1,to=sg$ydim,length.out=2*sg$ydim-1)
            radius <- 0.6
        }
        if(outer || !missing(nxo) || !missing(nyo)) {
            ## interpolation asked
            if(missing(nxo)) {
                xo <- x
            } else {
                xo <- seq(from=1,to=max(x),length.out=nxo)
            }
            if(missing(nyo)) {
                yo <- y
            } else {
                yo <- seq(from=1,to=max(y),length.out=nyo)
            }
            if(outer) {
                xo <- c(min(xo)-0.5,xo,max(xo)+0.5)
                yo <- c(min(yo)-0.5,yo,max(yo)+0.5)
            }
            depth <- choose.depth(nrow(distances))
            model <- terrainInterp(expand.grid(x,y),distances,depth,radius)
            distances <- matrix(predict(model,expand.grid(xo,yo)),
                                ncol=length(yo),nrow=length(xo))
            list(x=xo,y=yo,z=distances)
        } else {
            list(x=x,y=y,z=distances)
        }
    } else {
        ## interpolation is mandatory here
        if(mode=="mean") {
            depth <- choose.depth(sg$size)
            model <- terrainInterp(sg$pts,means,depth,1)
            if(missing(nxo)) {
                nxo <- 2*sg$xdim+1
            }
            if(missing(nyo)) {
                nyo <- 2*sg$ydim+1
            }
            distances <- means
            grid <- sg$pts
            radius <- 1.1
        } else {
            distances <- rep(NA,sg$size+(sg$xdim-1)*sg$ydim+
                             sg$xdim*(sg$ydim-1)+(sg$xdim-1)*(sg$ydim-1))
            grid <- matrix(NA,ncol=2,nrow=length(distances))
            ## first means
            grid[1:sg$size,] <- sg$pts
            distances[1:sg$size] <- means
            shift <- sg$size+1
            ## then right neighbors
            with.rn <- rep(1:(sg$xdim-1),sg$ydim)+
                rep((0:(sg$ydim-1))*sg$xdim,each=sg$xdim-1)
            indices <- shift:(shift+length(with.rn)-1)
            distances[indices] <- pdist[with.rn,1]
            grid[indices,] <- 0.5*(sg$pts[with.rn,]+sg$pts[with.rn+1,])
            shift <- max(indices)+1
            ## then vertical line (fixed x)
            thelines <- seq(from=1,by=sg$xdim,length.out=sg$ydim-1)+
                rep(0:(sg$xdim-1),each=sg$ydim-1)
            component <- rep(rep(2:3,length=sg$ydim-1),sg$xdim)
            indices <- shift:(shift+length(thelines)-1)
            distances[indices] <- pdist[cbind(thelines,component)]
            grid[indices,] <- 0.5*(sg$pts[thelines,]+sg$pts[thelines+sg$xdim,])
            shift <- max(indices)+1
            ## then remaining vertical lines (jagged x)
            otherlines <- (thelines+rep(1:0,length=length(thelines)))[1:(length(thelines)-sg$ydim+1)]
            component <- rep(rep(3:2,length=sg$ydim-1),sg$xdim-1)
            indices <- shift:(shift+length(otherlines)-1)
            distances[indices] <- pdist[cbind(otherlines,component)]
            grid[indices,] <- 0.5*(sg$pts[otherlines,]+
                                   sg$pts[otherlines+sg$xdim+
                                          rep(rep(c(-1,1),length.out=sg$ydim-1),times=sg$xdim-1),])
            if(missing(nxo)) {
                nxo <- 4*sg$xdim+1
            }
            if(missing(nyo)) {
                nyo <- 4*sg$ydim+1
            }
            radius <- 0.6
        }
        depth <- choose.depth(length(distances))
        model <- terrainInterp(grid,distances,depth,radius)
        xlim <- range(sg$pts[,1])
        ylim <- range(sg$pts[,2])
        x <- seq(xlim[1],xlim[2],length.out=nxo)
        y <- seq(ylim[1],ylim[2],length.out=nyo)
        if(outer) {
            x <- c(min(x)-0.5,x,max(x)+0.5)
            y <- c(min(y)-0.5,y,max(y)+0.5)
        }
        interpDist <- matrix(predict(model,expand.grid(x,y)),ncol=length(y))
        list(x=x,y=y,z=interpDist)
    }
}

plot.sompdist <- function(x,mode=c("mean","full"),...) {
    args <- list(...)
    args$mode <- match.arg(mode)
    if(is.null(args$border)) {
        args$border <- NA
    }
    pdist <- x$pdist
    sg <- x$somgrid
    if(sg$topo=="rectangular") {
        basevalues <- distance.grid(x,mode=args$mode)$z
        if(args$mode=="full") {
            plotsg <- somgrid(xdim=2*sg$xdim-1,ydim=2*sg$ydim-1,
                              topo="rectangular",with.dist=FALSE)
        } else {
            plotsg <- sg
        }
        values <- basevalues[cbind(rep(1:plotsg$xdim,plotsg$ydim),
                                   rep(plotsg$ydim:1,each=plotsg$xdim))]
    } else {
        if(args$mode=="mean") {
            values <- rowMeans(pdist,na.rm=T)
            plotsg <- sg
        } else {
            values <- rep(NA,2*sg$xdim*(2*sg$ydim-1))
            pindex <- 1:sg$size
            ## first fill mean values
            vshift <- seq(1,by=4*sg$xdim,length=sg$ydim)+
                rep(c(-1,0),length=sg$ydim)
            index.means <- rep(seq(1,by=2,length=sg$xdim),sg$ydim)+
                rep(vshift,each=sg$xdim)
            values[index.means] <- rowMeans(pdist,na.rm=T)
            ## then horizontal values
            index.h <- rep(seq(1,by=2,length=sg$xdim-1),sg$ydim)+
                rep(vshift,each=sg$xdim-1)+1
            values[index.h] <-
                pdist[pindex[-seq(sg$xdim,by=sg$xdim,length.out=sg$ydim)],1]
            if(sg$ydim>1) {
                ## then below right values
                vshift <- seq(2*sg$xdim,by=4*sg$xdim,length=sg$ydim-1) +
                    rep(c(0,1),length=sg$ydim-1)
                index.br <- rep(seq(1,by=2,length=sg$xdim),sg$ydim-1) +
                    rep(vshift,each=sg$xdim)
                index.br <- index.br[-seq(2*sg$xdim,by=2*sg$xdim,length=sg$ydim%/%2)]
                values[index.br] <- pdist[!is.na(pdist[,2]),2]
                ## then below left values
                vshift <- seq(2*sg$xdim,by=4*sg$xdim,length=sg$ydim-1) +
                    rep(c(1,0),length=sg$ydim-1)
                index.bl <- rep(seq(1,by=2,length=sg$xdim),sg$ydim-1) +
                    rep(vshift,each=sg$xdim)
                index.bl <- index.bl[-seq(sg$xdim,by=2*sg$xdim,length=sg$ydim%/%2)]
                values[index.bl] <- pdist[!is.na(pdist[,3]),3]

            }
            plotsg <- somgrid(2*sg$xdim,2*sg$ydim-1,topo="hexa",
                              with.dist=FALSE)
        }
    }
    args$mode <- NULL
    do.call("plot",c(list(x=plotsg,colorValues=values),args))
}

as.matrix.sompdist <- function(x, extended = TRUE,...) {
    sg <- x$somgrid
    pdist <- x$pdist
    result <- matrix(NA,ncol=sg$size,nrow=sg$size)
    diag(result) <- 0
    if(sg$topo=="rectangular") {
        delta <- c(1,sg$xdim+1,sg$xdim,sg$xdim-1)
        delta <- c(delta,-delta)
        baseindex <- 1:sg$size
        for(d in 1:8) {
            tofill <- !is.na(pdist[,d])
            index <- baseindex[tofill]
            result[cbind(index,index+delta[d])]=pdist[tofill,d]
        }
    } else {
        delta <- c(1,sg$xdim,sg$xdim-1,-1,-sg$xdim-1,-sg$xdim)
        baseindex <- 1:sg$size
        ## easy part
        for(d in c(1,4)) {
            tofill <- !is.na(pdist[,d])
            index <- baseindex[tofill]
            result[cbind(index,index+delta[d])]=pdist[tofill,d]
        }
        ## slightly more complex
        for(d in c(2,3,5,6)) {
            dec <- rep(delta[d],sg$size)+rep(rep(c(0,1),each=sg$xdim),length.out=sg$size)
            tofill <- !is.na(pdist[,d])
            index <- baseindex[tofill]
            result[cbind(index,index+dec[index])]=pdist[tofill,d]
        }
    }
    if(extended) {
        paths <- allShortestPaths(result)
        result <- paths$length
    }
    result
}

as.dist.sompdist <- function(x,FUN=NULL) {
    as.dist(as.matrix(x,extended=TRUE))
}
