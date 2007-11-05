grid2lines <- function(som) {
  if(som$somgrid$topo=="rectangular") {
    result=matrix(NA,
      ncol=ncol(som$prototypes),
      nrow=(som$somgrid$xdim+1)*(som$somgrid$ydim)+(som$somgrid$xdim)*(som$somgrid$ydim+1))
  } else {
    result=matrix(NA,
      ncol=ncol(som$prototypes),
      nrow=(som$somgrid$xdim+1)*(som$somgrid$ydim)+(som$somgrid$xdim)*(som$somgrid$ydim+1)+(som$somgrid$xdim-1)*(som$somgrid$ydim+1))    
  }
  # first "horizontal" lines
  shift <- (som$somgrid$xdim+1)*som$somgrid$ydim
  targetpos <- 1:(shift-1)
  targetpos <- targetpos[-(seq(from=som$somgrid$xdim+1,to=length(targetpos),by=som$somgrid$xdim+1))]
  result[targetpos,] <- som$prototypes
  # then "vertical" lines
  targetpos <- (shift+2):(shift+som$somgrid$xdim*(som$somgrid$ydim+1)+1)
  shift <- targetpos[length(targetpos)]
  targetpos <- targetpos[-(seq(from=som$somgrid$ydim+1,to=length(targetpos),by=som$somgrid$ydim+1))]
  result[targetpos,] <- som$prototypes[seq(from=1,by=som$somgrid$xdim,length.out=som$somgrid$ydim)+rep(0:(som$somgrid$xdim-1),each=som$somgrid$ydim),]
  # other "vertical" lines in hexagonal grids
  if(som$somgrid$topo=="hexagonal") {
    targetpos <- (shift+1):(shift+(som$somgrid$xdim-1)*(som$somgrid$ydim+1))
    targetpos <- targetpos[-(seq(from=som$somgrid$ydim+1,to=length(targetpos),by=som$somgrid$ydim+1))]
    result[targetpos,] <- som$prototypes[seq(from=2,by=som$somgrid$xdim,length.out=som$somgrid$ydim)+(rep(c(0,-1),length.out=som$somgrid$ydim)+rep(0:(som$somgrid$xdim-2),each=som$somgrid$ydim)),]
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

hitMap <- function(som,border=NA,...) {
  args <- list(...)
  if(is.null(args$col)) {
    args$col <- "red"
  }
  args$border <- border
  sizes <- table(factor(som$classif,levels=1:nrow(som$prototypes)))
  sizes <- sizes/max(sizes)
  do.call("plot",c(list(x=som$somgrid,size=sizes),args))
}

protoDist <- function(som,i,j,k,l) {
  from <- som$prototypes[i+(j-1)*som$somgrid$xdim,]
  to <- som$prototypes[k+(l-1)*som$somgrid$xdim,]
  sqrt(sum((from-to)^2))
}

hexNeighbor <- function(somgrid,i,j) {
  decY <- c(0,1,1,0,-1,-1)
  if(j%%2==1) {
    decX <- c(1,0,-1,-1,-1,0)
  } else {
    decX <- c(1,1,0,-1,0,1)
  }
  icand <- i+decX
  icand[icand<=0 | icand>somgrid$xdim] <- NA
  jcand <- j+decY
  jcand[jcand<=0 | jcand>somgrid$ydim] <- NA
  result <- cbind(icand,jcand)
  result[!is.na(result[,1])&!is.na(result[,2]),]
}

umatrix <- function(som,...) {
  args <- list(...)
  if(is.null(args$border)) {
    args$border <- NA
  }
  sg <- som$somgrid
  if(sg$topo=="rectangular") {
    # "easy" case 
    distances <- matrix(NA,ncol=2*sg$xdim-1,nrow=2*sg$ydim-1)
    # first fill the matrix with prototype distances
    for(i in 1:sg$xdim) {
      for(j in 1:sg$ydim) {
        # below neighbor
	if(j<sg$ydim) {
          distances[2*i-1,2*j] <- protoDist(som,i,j,i,j+1)
	}
        # diagonal (lower right)
	if(j<sg$ydim & i<sg$xdim) {
	  distances[2*i,2*j] <- protoDist(som,i,j,i+1,j+1)
	}
        # diagonal (lower left)
	if(j<sg$ydim & i>1) {
	  distances[2*i-2,2*j] <- protoDist(som,i,j,i-1,j+1)
	}
	# right
	if(i<sg$xdim) {
	  distances[2*i,2*j-1] <- protoDist(som,i,j,i+1,j)
	}
      }
    }
    # then smooth the distances
    di <- c(-1,-1,-1,0,0,1,1,1)
    dj <- c(-1,0,1,-1,1,-1,0,1)
    for(i in 1:sg$xdim) {
      for(j in 1:sg$ydim) {
        ni <- (2*i-1)+di
     	nj <- (2*j-1)+dj
     	cond <- ni>0&ni<=2*sg$xdim-1&nj>0&nj<=2*sg$ydim-1
     	ni <- ni[cond]
     	nj <- nj[cond]
     	tmp <- 0
     	for(k in 1:length(ni)) {
          tmp <- tmp+distances[ni[k],nj[k]]
     	}
     	distances[2*i-1,2*j-1] <- tmp/length(ni)
      }
    } 
    ugrid <- somgrid(xdim=2*sg$xdim-1,ydim=2*sg$ydim-1,topo="rectangular")
    do.call("plot",c(list(x=ugrid,colorValues=as.vector(distances)),args))
    invisible(distances)
  } else {
    # hard case (hexagonal)
    distances <- matrix(NA,nrow=2*sg$xdim,ncol=2*sg$ydim-1)
    for(i in 1:sg$xdim) {
      for(j in 1:sg$ydim) {
	# right
	if(i<sg$xdim) {
	  distances[2*i+1-j%%2,2*j-1] <- protoDist(som,i,j,i+1,j)
	}
	if(j%%2==1) {
	  if(j<sg$ydim) {
            # below right neighbor 
            distances[2*i-1,2*j] <- protoDist(som,i,j,i,j+1)
	    # below left neighbor
	    if(i>1) {
              distances[2*i-2,2*j] <- protoDist(som,i,j,i-1,j+1)
	    }
	  }
	} else {
	  if(j<sg$ydim) {
	    # below left neighbor
              distances[2*i-1,2*j] <- protoDist(som,i,j,i,j+1)
	    if(i<sg$xdim) {
              # below right neighbor 
              distances[2*i,2*j] <- protoDist(som,i,j,i+1,j+1)
	    }
	  }
	}
      }
    }
    ugrid <- somgrid(2*sg$xdim,2*sg$ydim-1,topo="hexa")
    for(i in 1:sg$xdim) {
      for(j in 1:sg$ydim) {
        nn <- hexNeighbor(ugrid,2*i-j%%2,2*j-1)
	distances[2*i-j%%2,2*j-1] <- mean(distances[nn],na.rm=TRUE)
      }
    } 
    do.call("plot",c(list(x=ugrid,colorValues=as.vector(distances)),args))
    invisible(distances)
  }
}

# very basic version
plot.som <- function(x,y,...) {
  if(missing(y)) {
    stars(x$prototypes,locations=x$somgrid$pts,len=0.5*x$somgrid$basesize,radius=F,...)
  } else {
    if(is.factor(y)) {
      ymat <- mapFactorToUnit(x,y)
      stars(ymat,locations=x$somgrid$pts,len=0.5*x$somgrid$basesize,draw.segments=T,col.segments=rainbow(ncol(ymat)),...)
    } else {
      stars(y,locations=x$somgrid$pts,len=0.5*x$somgrid$basesize,...)
    }
  }
}
