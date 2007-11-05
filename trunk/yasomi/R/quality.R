som2graph <- function(som) {
  sg <- som$somgrid
  distances <- matrix(NA,ncol=sg$size,nrow=sg$size)
  if(sg$topo=="rectangular") {
    # first fill the matrix with prototype distances
    for(j in 1:sg$ydim) {
      base <- (j-1)*sg$xdim
      for(i in 1:sg$xdim) {
        current <- base + i
        # below neighbor
	if(j<sg$ydim) {
          distances[current,current+sg$xdim] <- protoDist(som,i,j,i,j+1)
	}
        # diagonal (lower right)
	if(j<sg$ydim & i<sg$xdim) {
	  distances[current,current+sg$xdim+1] <- protoDist(som,i,j,i+1,j+1)
	}
        # diagonal (lower left)
	if(j<sg$ydim & i>1) {
	  distances[current,current+sg$xdim-1] <- protoDist(som,i,j,i-1,j+1)
	}
	# right
	if(i<sg$xdim) {
	  distances[current,current+1] <- protoDist(som,i,j,i+1,j)
	}
      }
    }
    # then fill missing
    distances[lower.tri(distances)] <- t(distances)[lower.tri(distances)]
    matrix(distances,ncol=sg$size)
  } else {
    for(j in 1:sg$ydim) {
      base <- (j-1)*sg$xdim
      for(i in 1:sg$xdim) {
        current <- base + i
	# right
	if(i<sg$xdim) {
	  distances[current,current+1] <- protoDist(som,i,j,i+1,j)
	}
	if(j%%2==1) {
	  if(j<sg$ydim) {
            # below right neighbor 
            distances[current,current+sg$xdim] <- protoDist(som,i,j,i,j+1)
	    # below left neighbor
	    if(i>1) {
              distances[current,current+sg$xdim-1] <- protoDist(som,i,j,i-1,j+1)
	    }
	  }
	} else {
	  if(j<sg$ydim) {
	    # below left neighbor
              distances[current,current+sg$xdim] <- protoDist(som,i,j,i,j+1)
	    if(i<sg$xdim) {
              # below right neighbor 
              distances[current,current+sg$xdim+1] <- protoDist(som,i,j,i+1,j+1)
	    }
	  }
	}
      }
    }
    distances[lower.tri(distances)] <- t(distances)[lower.tri(distances)]
    matrix(distances,ncol=sg$size)    
  }
}

KaskiLagus <- function(som,data) {
  data <- as.matrix(data)
  winners <- secondBmu(som$prototypes,data)
  graph <- som2graph(som)
  paths <- allShortestPaths(graph)
  list(quant=sqrt(rowSums((data-som$prototypes[winners[,1],])^2)),path=paths$length[winners])
}

error.kaskilagus <- function(som,data) {
  pre <- KaskiLagus(som,data)
  mean(pre$quant+pre$path)
}

error.quantisation <- function(som,data) {
  bmu(som$prototypes,data)$error
}
