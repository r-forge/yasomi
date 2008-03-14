som.tunecontrol <- function(somgrid,init="pca",ninit=1,
                            assignment="single",
                            radii=c(2,2/3*somgrid$diam),nradii=10,
                            innernradii=30,maxiter=75,
                            annealing="power",
                            kernel="gaussian",criterion=error.quantisation) {
    if(init=="pca" && ninit>1) {
        warning("PCA based initialization is deterministic")
    }
    list(init=match.arg(init,c("pca","random"),several.ok = TRUE),
         ninit=ninit,
         assignment=match.arg(assignment,c("single","heskes"),several.ok = TRUE),
         radii=seq(from=min(radii),to=max(radii),length.out=nradii),
         innernradii=innernradii,
         maxiter=maxiter,
         annealing=match.arg(annealing,c("power","linear"),several.ok = TRUE),
         kernel=match.arg(kernel,c("gaussian","linear"),several.ok = TRUE),
         criterion=criterion)
}


som.tune <- function(data,somgrid,control=som.tunecontrol(somgrid),
                     verbose=FALSE) {
    nbconf <- length(control$init)*control$ninit*length(control$assignment)*length(control$radii)*length(control$annealing)*length(control$kernel)
    dimensions <- character()
    if(length(control$init)>1) {
        dimensions <- c(dimensions,"Initialisation method")
    }
    if(control$ninit>1) {
        dimensions <- c(dimensions,"Initialisation")
    }
    if(length(control$assignment)>1) {
        dimensions <- c(dimensions,"Assignment method")
    }
    if(length(control$radii)>1) {
        dimensions <- c(dimensions,"Initial radius")
    }
    if(length(control$annealing)>1) {
        dimensions <- c(dimensions,"Annealing method")
    }
    if(length(control$kernel)>1) {
        dimensions <- c(dimensions,"Kernel")
    }
    bestSOM <- NULL
    performances <- rep(NA,nbconf)
    comp.init <- rep("",nbconf)
    comp.assign <- rep("",nbconf)
    comp.radii <- rep(NA,nbconf)
    comp.anneal <- rep("",nbconf)
    comp.kernel <- rep("",nbconf)
    if(!identical(control$criterion,error.quantisation)) {
        quantisation <- rep(NA,nbconf)
    } else {
        quantisation <- NULL
    }
    bestPerfSoFar <- Inf
    confIndex <- 1
    ## compute distances if they are missing
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    ## initialization method
    for(init in control$init) {
        ## assignment
        for(assignment in control$assignment) {
            ## kernel
            for(kernel in control$kernel) {
                if(kernel=="gaussian") {
                    minRadius <- 0.5
                } else {
                    minRadius <- 1
                }
                ## annealing
                for(annealing in control$annealing) {
                    ## initializations
                    for(i in 1:control$ninit) {
                        if(init=="pca") {
                            prototypes=sominit(data,somgrid)
                        } else {
                            prototypes=somRandomInit(data,somgrid)
                        }
                        ## radii
                        for(radius in control$radii) {
                            if(verbose) {
                                print(paste("Configuration ",confIndex,"/",nbconf,sep=""))
                            }
                            radii <- switch(annealing,
                                            "linear"=radius.lin(minRadius,radius,control$innernradii),
                                            "power"=radius.exp(minRadius,radius,control$innernradii))
                            som <- batchsom(data,somgrid,prototypes,assignment,radii,
                                            maxiter=control$maxiter,kernel=kernel)
                            performances[confIndex] <- control$criterion(som,data)
                            comp.init[confIndex] <- init
                            comp.assign[confIndex] <- assignment
                            comp.radii[confIndex] <- radius
                            comp.anneal[confIndex] <- annealing
                            comp.kernel[confIndex] <- kernel
                            if(performances[confIndex]<bestPerfSoFar) {
                                bestSOM <- som
                                bestPerfSoFar <- performances[confIndex]
                            }
                            if(!is.null(quantisation)) {
                                quantisation[confIndex] <- error.quantisation(som)
                            }
                            confIndex <- confIndex + 1
                        }
                    }
                }
            }
        }
    }
    if(is.null(quantisation)) {
        quantisation <- performances
        isquant <- TRUE
    } else {
        isquant <- FALSE
    }
    res <- list(best.som=bestSOM,quantisation=quantisation,errors=performances,
                isquant=isquant,control=control,dimensions=dimensions,
                init=comp.init,assignement=comp.assign,radii=comp.radii,
                annealing=comp.anneal,kernel=comp.kernel)
    class(res) <- "somtune"
    res
}

plot.somtune <- function(x,relative=TRUE,...) {
    args <- list(...)
    if(length(x$dimensions)==1) {
        if(is.null(args$xlab)) {
            args$xlab <- x$dimensions[1]
        }
        if(is.null(args$ylab)) {
            if(x$isquant) {
                args$ylab <- "Quantisation error"
            } else {
                args$ylab <- "Error measures"
            }
        }
        if(is.null(args$yaxs)) {
            args$yaxs <- "r"
        }
        if(x$isquant) {
            values <- x$quantisation
        } else {
            values <- rbind(x$quantisation,x$errors)
            args$beside <- TRUE
            if(is.null(args$legend.text)) {
                args$legend.text <- c("Quantisation","Error")
            }
        }
        if(relative) {
            baseline <- 0.9*min(values)
            values <- values-baseline
            args$offset <- baseline
        }
        config.names <- switch(x$dimensions[1],
                               "Initialisation method"=x$init,
                               "Initialisation"=1:length(x$quantisation),
                               "Assignment method"=x$assignement,
                               "Initial radius"=round(x$radii,digits=2),
                               "Annealing method"=x$annealing,
                               "Kernel"=x$kernel)
        do.call("barplot",
                c(list(height=values,names.arg=config.names,axis.lty=1),args))
    } else {
        stop("cannot plot an object of class 'somtune' with more than one dimension")
    }
}
