## in quality.R
error.quantisation <- function(som,newdata,...) {
    UseMethod("error.quantisation")
}

error.kaskilagus <- function(som,newdata,...) {
    UseMethod("error.kaskilagus")
}

## in umatrix.R
protoDist <- function(som,i,j,k,l) {
    UseMethod("protoDist")
}

## in som.R (slow-versions.R) and relational.R (slow-relational.R)

sominit.pca <- function(data,somgrid,...) {
    UseMethod("sominit.pca")
}

sominit.random <- function(data,somgrid,method=c("prototypes","random","cluster"),...) {
    UseMethod("sominit.random")
}

batchsom <- function(data,somgrid,init,prototypes,...) {
    UseMethod("batchsom")
}

## in kernel.R

as.kernelmatrix <- function(data,...) {
    UseMethod("as.kernelmatrix")
}
