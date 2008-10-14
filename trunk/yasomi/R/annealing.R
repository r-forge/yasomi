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

batchsom.control.default <- function(data,somgrid,
                                     mode = c("continuous","stepwise"),
                                     min.radius, max.radius, steps,
                                     decrease = c("power", "linear"), max.iter,
                                     kernel = c("gaussian", "linear"),
                                     normalised,
                                     assignment = c("single", "heskes"),
                                     cut = 1e-07,...)
{
    mode <- match.arg(mode,c("stepwise","continuous"))
    decrease <- match.arg(decrease,c("power","linear"))
    kernel <- match.arg(kernel,c("gaussian","linear"))
    assignment <- match.arg(assignment,c("single", "heskes"))
    if(missing(max.radius)) {
        max.radius <- 2/3*somgrid$diam+1
    }
    if(missing(min.radius)) {
        min.radius <- switch(kernel,"gaussian"=0.5,"linear"=1)
    }
    if(max.radius<=min.radius) {
        stop("max.radius must be larger than min.radius")
    }
    if(min.radius<=0) {
        stop("min.radius must be positive")
    }
    if(missing(steps)) {
        steps <- switch(mode,
                        "stepwise"=max(2,ceiling(2*max.radius)),
                        "continuous"=max(20,ceiling(5*max.radius)))
    }
    if(missing(max.iter)) {
        max.iter <- switch(mode,
                           "stepwise"=75,
                           "continuous"=1)
    }
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    radii <- switch(decrease,
                    "linear"=radius.lin(min.radius,max.radius,steps),
                    "power"=radius.exp(min.radius,max.radius,steps))
    list(mode=mode,radii=radii,max.iter=max.iter,kernel=kernel,
         normalised=normalised,assignment=assignment,
         cut=cut)
}
