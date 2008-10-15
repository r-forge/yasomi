library(yasomi)
## test data weights support
datadim <- ceiling(4+runif(1,min=1,max=16))
datasize <- datadim*100
X <- matrix(nrow=datasize,ncol=datadim,runif(datasize*datadim))

sg <- somgrid(5,6,topo="rectangular")

## random initial prototypes
pinit.choice <- sample(nrow(X),sg$size)
prototypes.init <- X[pinit.choice,]

radius.max <- runif(1,min=1/3,max=4/3)*sg$diam
nb.radii <- ceiling(runif(1,min=radius.max,max=3*radius.max))
mode <- sample(c("continuous","stepwise"),1)

## random integer weights
weights <- sample(5,size=nrow(X),replace=T,prob=exp(-seq(0,2,length=5)))

## with weight support
som <- batchsom(X,sg,prototypes=prototypes.init,weights=weights,verbose=TRUE,
                max.radius=radius.max,steps=nb.radii,mode=mode)

## weighting via replication
X.rep <- X[rep(1:nrow(X),times=weights),]
som.rep <- batchsom(X.rep,sg,prototypes=prototypes.init,verbose=TRUE,
                    max.radius=radius.max,steps=nb.radii,mode=mode)

stopifnot(all.equal(error.quantisation(som),error.quantisation(som.rep)))
stopifnot(all.equal(som.rep$errors,som$errors))
stopifnot(all.equal(error.kaskilagus(som),error.kaskilagus(som.rep)))
stopifnot(all.equal(som$classif[rep(1:nrow(X),times=weights)],som.rep$classif))
stopifnot(all.equal(som.rep$prototypes,som$prototypes))

## now relational SOM
dX <- dist(X)
rprototypes.init <- matrix(0,nrow=sg$size,ncol=nrow(X))
rprototypes.init[cbind(1:sg$size,pinit.choice)] <- 1

rsom <- batchsom(dX,sg,prototypes=rprototypes.init,weights=weights,verbose=TRUE,
                 max.radius=radius.max,steps=nb.radii,mode=mode)

rprototypes.final <- rsom$prototypes%*%X

stopifnot(all.equal(error.quantisation(som),error.quantisation(rsom)))
stopifnot(all.equal(som$errors,rsom$errors))
stopifnot(all.equal(error.kaskilagus(som),error.kaskilagus(rsom)))
stopifnot(all.equal(som$classif,rsom$classif))
stopifnot(all.equal(som$prototypes,rprototypes.final,
                    check.attributes = FALSE, check.names = FALSE))

## and finally kernel SOM
KX <- as.kernelmatrix(tcrossprod(X))

ksom <- batchsom(KX,sg,prototypes=rprototypes.init,weights=weights,verbose=TRUE,
                 max.radius=radius.max,steps=nb.radii,mode=mode)
stopifnot(all.equal(error.kaskilagus(som),error.kaskilagus(ksom)))
stopifnot(all.equal(error.quantisation(som),error.quantisation(ksom)))
stopifnot(all.equal(som$errors,ksom$errors))
stopifnot(all.equal(som$classif,ksom$classif))
stopifnot(all.equal(som$prototypes,rprototypes.final,
                    check.attributes = FALSE, check.names = FALSE))
