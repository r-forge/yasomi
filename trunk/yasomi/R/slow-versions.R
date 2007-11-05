neighborhood <- function(somgrid,T,kernel,normalised=TRUE) {
  raw <- kernel(somgrid$dist/T)
  if(normalised) {
    sweep(raw,1,rowSums(raw),"/")
  } else {
    raw
  }
}

kernel.gaussian <- function(x) {
  exp(-(3*x)^2/2)
}

kernel.linear <- function(x) {
  pre <- sapply(x,function(x) {
      if(x<1) {1-x} else {0}
    })
  if(is.matrix(x)) {
    matrix(pre,ncol=ncol(x),nrow=nrow(x))
  } else {
    pre
  }
}

bmu.R <- function(prototypes,data) {
  result=rep(0,nrow(data))
  for(i in 1:length(result)) {
    result[i]=which.min(rowSums(sweep(prototypes,2,data[i,],"-")^2))
  }
  result
}

bmu.heskes.R <- function(prototypes,data,nv) {
  result=rep(0,nrow(data))
  for(i in 1:length(result)) {
    result[i]=which.min(nv%*%rowSums(sweep(prototypes,2,data[i,],"-")^2))
  }
  result
}

batchsom.R <- function(somgrid,data,prototypes=somPCAInit(somgrid,data),
                       assignment=c("single","heskes"),radii,nbRadii,
                       maxiter=75,
                       kernel=c("gaussian","linear"),normalised,
                       cut=1e-7,verbose=FALSE) {
  # process parameters
  assignment <- match.arg(assignment)
  if(missing(normalised)) {
    normalised <- assignment=="heskes"
  }
  kernel <- match.arg(kernel)
  theKernel <- switch(kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear)
  # perform a few sanity checks 
  if(ncol(prototypes)!=ncol(data)) {
    stop("'prototypes' and 'data' have different dimensions")
  }
  if(class(somgrid)!="somgrid") {
    stop("'somgrid' is not of somgrid class")
  }
  # compute radii
  if(missing(radii)) {
    if(kernel=="gaussian") {
      minRadius <- 0.5
    } else {
      minRadius <- 1
    }
    radii <- radius.exp(minRadius,max(minRadius,somgrid$diam/3*2),nbRadii)
  }
  batchsom.lowlevel.R(somgrid,data,prototypes,assignment,radii,
                      maxiter,theKernel,normalised,cut,verbose)
}

batchsom.lowlevel.R <- function(somgrid,data,prototypes,
                                assignment,radii,maxiter,kernel,
                                normalised,cut,verbose) {
  data <- as.matrix(data)
  classif <- rep(NA,nrow(data))
  errors <- vector("list",length(radii))
  for(i in 1:length(radii)) {
    nv <- neighborhood(somgrid,radii[i],kernel,normalised=normalised)
    for(j in 1:maxiter) {
      if(assignment == "single") {
        bmus <- bmu(prototypes,data)
      } else {
        bmus <- bmu.heskes(prototypes,data,nv)
      }
      nclassif <- bmus$clusters
      noChange = identical(classif,nclassif)
      classif <- nclassif
      error <- bmus$error
      if(verbose) {
        print(paste(i,j,error))
      }
      errors[[i]] <- c(errors[[i]],error)
      weights <- nv[,classif]
      normed <- rowSums(weights)
      mask <- (1:length(normed))[(normed/length(normed))>cut]
      prototypes[mask,] <- sweep(weights%*%data,1,normed,"/")[mask,]
      if(noChange && j >=2) {
        if(verbose) {
          print(paste("radius:",radii[i],"iteration",j,"is stable, decreasing radius"))
	}	
        break;
      }      
    }
    if(!noChange && verbose) {
      print(paste("warning: can't reach a stable configuration with radius",i))
    }
  }
  res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,errors=unlist(errors))
  class(res) <- "som"
  res
  
}
