partialSums.R <- function(cluster,nclust,diss) {
    fcluster <- factor(cluster,levels=1:nclust)
    pre <- apply(diss,2,tapply,fcluster,sum)
    pre[is.na(pre)] <- 0
    dimnames(pre) <- NULL
    bipre <- apply(pre,1,tapply,fcluster,sum)
    bipre[is.na(bipre)] <- 0
    dimnames(bipre) <- NULL
    list(ps=pre,bips=bipre)
}

weightedPartialSums.R <- function(cluster,nclust,diss,weights) {
    fcluster <- factor(cluster,levels=1:nclust)
    pre <- apply(diss,2,function(x) {tapply(x*weights,fcluster,sum)})
    pre[is.na(pre)] <- 0
    dimnames(pre) <- NULL
    bipre <- apply(pre,1,function(x) {tapply(x*weights,fcluster,sum)})
    bipre[is.na(bipre)] <- 0
    dimnames(bipre) <- NULL
    list(ps=pre,bips=bipre)
}


relationalsom.lowlevel.R <- function(somgrid,diss,prototypes,
                                     assignment,radii,weights,maxiter,kernel,
                                     normalised,cut,verbose,extended,
                                     data.norms) {
    classif <- rep(NA,nrow(diss))
    oldClassif <- rep(NA,nrow(diss))
    errors <- vector("list",length(radii))
    if(missing(weights) || is.null(weights)) {
        weights <- rep(1,nrow(diss))
    }
    nv <- neighborhood(somgrid,radii[1],kernel,normalised=normalised)
    for(i in 1:length(radii)) {
        for(j in 1:maxiter) {
            ## assignment
            if(extended) {
                if(assignment == "single") {
                    bmus <- extended.relationalbmu.R(prototypes,diss,data.norms,weights)
                } else {
                    stop(paste(assignment,"is not implemented for relational SOM"))
                }
                ## no more extended after the first iteration
                extended <- FALSE
            } else {
                if(assignment == "single") {
                    bmus <- relationalbmu.R(prototypes,diss,weights)
                } else {
                    stop(paste(assignment,"is not implemented for relational SOM"))
                }
            }
            nclassif <- bmus$clusters
            noChange <- identical(classif,nclassif)
            hasLoop <- identical(oldClassif,nclassif)
            oldClassif <- classif
            classif <- nclassif
            error <- bmus$error
            if(verbose) {
                print(paste(i,j,error))
            }
            errors[[i]] <- c(errors[[i]],error)
            ## prepare next iteration
            if(noChange || hasLoop) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i==length(radii)) {
                    ## the fitting is done
                    break;
                }
                ## preparing the loop with the next radius
                nv <- neighborhood(somgrid,radii[i+1],kernel,
                                   normalised=normalised)
            }
            ## always update the weights
            nvcl <-  sweep(nv[,classif],2,weights,"*")
            normed <- rowSums(nvcl)
            prototypes <- sweep(nvcl,1,normed,"/")
            ## break the loop if the partition is stable or when we have
            ## an oscillating behaviour
            if(noChange || hasLoop) {
                break;
            }
        }
        if(!noChange && verbose) {
            print(paste("warning: can't reach a stable configuration with radius",i))
        }
    }
    if(!noChange) {
        ## final assignment
        bmus <- relationalbmu.R(prototypes,diss,weights)
        classif <- bmus$clusters
        errors[[length(radii)]] <- c(errors[[length(radii)]],bmus$error)
    }
    Dalpha <- tcrossprod(diss,prototypes)
    nf <- 0.5*diag(prototypes%*%Dalpha)
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Dalpha=Dalpha,nf=nf)
    class(res) <- c("relationalsom","som")
    res
}

