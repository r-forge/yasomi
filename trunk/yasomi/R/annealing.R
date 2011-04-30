radius.exp <- function(min,max,steps) {
    max*((min+1)/max)^(seq(0,1,length.out=steps))-1
}

radius.lin <- function(min,max,steps) {
    seq(max,min,length.out=steps)
}

batchsom.control.default <- function(data,somgrid,
                                     mode = c("continuous","stepwise"),
                                     min.radius, max.radius, steps,
                                     decrease = c("power", "linear"), max.iter,
                                     kernel = c("gaussian", "linear", "zeroone"),
                                     normalised,
                                     assignment = c("single", "heskes"),
                                     cut = 1e-07,...)
{
    mode <- match.arg(mode,c("continuous","stepwise"))
    decrease <- match.arg(decrease,c("power","linear"))
    kernel <- match.arg(kernel,c("gaussian","linear","zeroone"))
    assignment <- match.arg(assignment,c("single", "heskes"))
    if(missing(max.radius)) {
        max.radius <- 2/3*somgrid$diam+1
    }
    if(missing(min.radius)) {
        min.radius <- 0
    }
    if(max.radius<=min.radius) {
        stop("max.radius must be larger than min.radius")
    }
    if(min.radius<0) {
        stop("min.radius must be non negative")
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
