## in quality.R
KaskiLagus <- function(som,data) {
    UseMethod("KaskiLagus")
}

## in umatrix.R
protoDist <- function(som,i,j,k,l) {
    UseMethod("protoDist")
}

## in som.R (slow-versions.R) and relational.R (slow-relational.R)

sominit <- function(data,somgrid,...) {
    UseMethod("sominit")
}

batchsom <- function(data,somgrid,...) {
    UseMethod("batchsom")
}
