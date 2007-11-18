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

distance.grid <- function(sompdist,mode=c("mean","full")) {
    pdist <- sompdist$pdist
    sg <- sompdist$somgrid
    mode <- match.arg(mode)
    if(sg$topo=="rectangular") {
        means <- rowMeans(pdist,na.rm=T)
        if(mode=="mean") {
            reorder <-
                rep(1:sg$xdim,sg$ydim)+rep((sg$ydim-1):0,each=sg$xdim)*sg$xdim
            distances <- matrix(means[reorder],ncol=sg$ydim)
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
        }
    }
    distances
}

