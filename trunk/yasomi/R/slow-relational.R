partialSums.R <- function(cluster,nclust,diss) {
    fcluster <- factor(cluster,levels=1:nclust)
    pre <- apply(diss,2,tapply,fcluster,sum)
    pre[is.na(pre)] <- 0
    bipre <- apply(pre,1,tapply,fcluster,sum)
    bipre[is.na(bipre)] <- 0
    list(ps=pre,bips=bipre)
}

relationalsom.lowlevel.R <- function(somgrid,diss,prototypes,
                                     assignment,radii,maxiter,kernel,
                                     normalised,cut,verbose) {
    classif <- rep(NA,nrow(diss))
    oldClassif <- rep(NA,nrow(diss))
    errors <- vector("list",length(radii))
    nv <- neighborhood(somgrid,radii[1],kernel,normalised=normalised)
    for(i in 1:length(radii)) {
        for(j in 1:maxiter) {
            ## assignment
            if(assignment == "single") {
                bmus <- relationalbmu.R(prototypes,diss)
            } else {
                stop(paste(assignment,"is not implemented for relational SOM"))
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
            if(!noChange ||  j>=2 || j==maxiter || hasLoop) {
                ## representation is needed when:
                ##  - the partition has changed
                ##  - or with a new value of the radius if
                ##  -- the partition has not changed
                ##  -- a loop has been detected
                ##  -- the maximal number of iteration has been reached
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
                weights <- nv[,classif]
                normed <- rowSums(weights)
                prototypes <- sweep(weights,1,normed,"/")
            }
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
    Dalpha <- tcrossprod(diss,prototypes)
    nf <- 0.5*diag(prototypes%*%Dalpha)
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Dalpha=Dalpha,nf=nf)
    class(res) <- c("relationalsom","som")
    res
}

