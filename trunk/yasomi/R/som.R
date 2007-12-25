sominit.default <- function(data,somgrid,...) {
    ## we don't scale the data
    data.pca <- prcomp(data)
    ## the more detailled axis is assigned to the first eigenvector
    if (somgrid$xdim>=somgrid$ydim) {
        x.ev <- 1
        y.ev <- 2
    } else {
        x.ev <- 2
        y.ev <- 1
    }
    if(somgrid$topo=="hexagonal") {
        xspan <- somgrid$xdim - 1
        if(somgrid$ydim>1) {
            xspan <- xspan+0.5
        }
        x <- seq(from=-2*data.pca$sdev[x.ev],by=4*data.pca$sdev[x.ev]/xspan,length.out=somgrid$xdim)
    } else {
        x <- seq(from=-2*data.pca$sdev[x.ev],to=2*data.pca$sdev[x.ev],length.out=somgrid$xdim)
    }
    y <- seq(from=2*data.pca$sdev[y.ev],to=-2*data.pca$sdev[y.ev],length.out=somgrid$ydim)
    base <- as.matrix(expand.grid(x = x, y = y))
    ## correction for hexagonal grids
    if(somgrid$topo=="hexagonal") {
        base[,1] <- base[,1]+rep(c(0,2*data.pca$sdev[x.ev]/xspan),each=somgrid$xdim,length.out=nrow(base))
    }
    ## map back the grid to the original space
    mapped <- tcrossprod(base,data.pca$rotation[,c(x.ev,y.ev)])
    ## decentering
    sweep(mapped,2,data.pca$center,"+")
}

somRandomInit <- function(data,somgrid) {
    data <- as.matrix(data)
    if(nrow(data)>=somgrid$size) {
        data[sample(1:nrow(data),size=somgrid$size),]
    } else {
        data[sample(1:nrow(data),size=somgrid$size,replace=TRUE),]
    }
}

bmu <- function(prototypes,data) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    result <- .C("bmu",
                 as.double(prototypes),
                 as.integer(nrow(prototypes)),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(prototypes)),
                 clusters=integer(nrow(data)),
                 error=as.double(0),
                 PACKAGE="yasomi")
    list(clusters=result$clusters+1,error=result$error)
}

secondBmu <- function(prototypes,data) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    matrix(.C("second_bmu",
              as.double(prototypes),
              as.integer(nrow(prototypes)),
              as.double(data),
              as.integer(nrow(data)),
              as.integer(ncol(prototypes)),
              clusters=integer(2*nrow(data)),
              PACKAGE="yasomi")$clusters,ncol=2)
}

bmu.heskes <- function(prototypes,data,nv) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    if(ncol(nv)!=nrow(nv)) {
        stop("'nv' is not a square matrix")
    }
    if(ncol(nv)!=nrow(prototypes)) {
        stop("'nv' and 'prototypes' have different dimensions")
    }
    ## nv must be in row major mode if normalised
    result <- .C("bmu_heskes",
                 as.double(prototypes),
                 as.double(t(nv)),
                 as.integer(nrow(prototypes)),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(prototypes)),
                 clusters=integer(nrow(data)),
                 error=as.double(0),
                 PACKAGE="yasomi")
    list(clusters=result$clusters+1,error=result$error)
}

radius.exp <- function(min,max,steps) {
    max*(min/max)^(seq(0,1,length.out=steps))
}

radius.lin <- function(min,max,steps) {
    seq(max,min,length.out=steps)
}

batchsom.default <- function(data,somgrid,prototypes=sominit(data,somgrid),
                     assignment=c("single","heskes"),radii,nbRadii=30,
                     maxiter=75,
                     kernel=c("gaussian","linear"),
                     normalised,
                     cut=1e-7,verbose=FALSE,...) {
    ## process parameters
    assignment <- match.arg(assignment)
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    theRule <- switch(assignment,"single"=0,"heskes"=1)
    kernel <- match.arg(kernel)
    kernelType <- switch(kernel,"gaussian"=0,"linear"=1)
    ## perform a few sanity checks
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    ## distances?
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    ## compute radii
    if(missing(radii)) {
        if(kernel=="gaussian") {
            minRadius <- 0.5
        } else {
            minRadius <- 1
        }
        radii <- radius.exp(minRadius,max(minRadius,somgrid$diam/3*2),nbRadii)
    }
    pre <- batchsom.lowlevel(somgrid,data,prototypes,theRule,radii,
                             maxiter,kernelType,normalised,cut,verbose)
    pre$assignment <- assignment
    pre$kernel <- kernel
    pre$normalised <- normalised
    pre$radii <- radii
    pre
}

batchsom.lowlevel <- function(somgrid,data,prototypes,
                              assignment,radii,maxiter,
                              kernelType,normalised,cut,verbose) {
    result <- .C("batch_som_optim",
                 proto=as.double(prototypes),
                 as.integer(somgrid$size),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(data)),
                 as.integer(assignment),
                 as.double(somgrid$dist),
                 as.double(radii),
                 as.integer(length(radii)),
                 as.integer(maxiter),
                 as.integer(kernelType)[1],
                 as.integer(normalised)[1],
                 as.double(cut)[1],
                 as.integer(verbose),
                 clusters=integer(nrow(data)),
                 errors=as.double(rep(-1,length(radii)*maxiter)),
                 PACKAGE="yasomi")
    prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                         dimnames=list(NULL,dimnames(data)[[2]]))
    res <- list(somgrid=somgrid,
                prototypes=prototypes,
                classif=result$cluster+1,
                errors=result$errors[result$errors>=0])
    class(res) <- c("som","somnum")
    res
}


colorCode <- function(data,nbcolor) {
    onedimgrid <- somgrid(xdim=nbcolor,ydim=1,topo="rectangular")
    colorsom <- batchsom(data,onedimgrid,nbRadii=20)
    colorsom$classif
}

mapToUnit <- function(som,values) {
    if(length(values)!=length(som$classif)) {
        stop("'values' is not of the same size as data use to fit the 'som'")
    }
    result <- vector("list",nrow(som$prototypes))
    if(is.null(dim(values))) {
        for(i in 1:length(result)) {
            result[[i]] <- values[som$classif==i]
        }
    } else {
        for(i in 1:length(result)) {
            result[[i]] <- values[som$classif==i,]
        }
    }
    result
}

mapFactorToUnit <- function(som,values) {
    if(length(values)!=length(som$classif)) {
        stop("'values' is not of the same size as data use to fit the 'som'")
    }
    if(!is.factor(values)) {
        stop("'values' is not a factor")
    }
    lv <- levels(values)
    result <- matrix(0,nrow=nrow(som$prototypes),ncol=length(lv),
                     dimnames=list(c(),lv))
    coded <- unclass(values)
    for(i in 1:length(som$classif)) {
        result[som$classif[i],coded[i]] <- result[som$classif[i],coded[i]] + 1
    }
    result
}


predict.somnum <- function(object,newdata,...) {
    som <- object
    newdata <- as.matrix(newdata)
    if(ncol(newdata)!=ncol(som$prototypes)) {
        stop("'newdata' and 'object$prototypes' have different dimensions")
    }
    pre <- bmu(som$prototypes,newdata)
    list(classif=pre$clusters,error=pre$error)
}

print.som <- function(x,...)
{
    cat("\nSelf-Organising Map\n\n")
    cat("Parameters:\n")
    cat("       grid: ",x$somgrid$topo," grid of size ",x$somgrid$xdim,"x",
        x$somgrid$ydim," with diameter ",x$somgrid$diam,"\n",sep="")
    cat("     kernel: ",x$kernel,
        if(!x$normalised) {" not normalised"} else {" normalised"},"\n",sep="")
    cat(" assignment:",x$assignment,"\n")
    if(length(x$radii)<6) {
        cat("      radii:",x$radii,"\n")
    } else {
        cat("      radii:",length(x$radii),"values from",max(x$radii),
            "down to",min(x$radii),"\n")
    }
    cat("\n")
    invisible(x)
}


summary.som <- function(object,...)
{
    cat("\nSelf-Organising Map\n\n")
    cat("Quantisation error: ",object$errors[length(object$errors)],"\n\n",
        sep="")
    cat("Parameters:\n")
    cat("       grid: ",object$somgrid$topo," grid of size ",
        object$somgrid$xdim,"x",object$somgrid$ydim," with diameter ",
        object$somgrid$diam,"\n",sep="")
    cat("     kernel: ",object$kernel,
        if(!object$normalised) {" not normalised"} else {" normalised"},
        "\n",sep="")
    cat(" assignment:",object$assignment,"\n")
    cat("      radii:",object$radii,"\n")
    invisible()
}
