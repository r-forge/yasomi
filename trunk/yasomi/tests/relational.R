library(yasomi)
## test that the relational SOM gives identical results as the standard SOM
datadim <- ceiling(4+runif(1,min=1,max=16))
datasize <- datadim*100

X <- matrix(nrow=datasize,ncol=datadim,runif(datasize*datadim))
dX <- dist(X) # default to Euclidean distances

sg <- somgrid(datadim,ceiling(datadim+runif(1,min=-2,max=4)),topo="rectangular")

## random initialization by random prototypes
rprototypes.init <- sominit.random(dX,sg,method="prototypes")
## let's get back the chosen data points
init.points <- apply(rprototypes.init,1,which.max)
prototypes.init <- X[init.points,]

radius.max <- runif(1,min=1/3,max=4/3)*sg$diam
nb.radii <- ceiling(runif(1,min=radius.max,max=3*radius.max))
mode <- sample(c("continuous","stepwise"),1)

## let's fit the SOM
rsom <- batchsom(dX,sg,prototypes=rprototypes.init,verbose=TRUE,
                 max.radius=radius.max,steps=nb.radii,mode=mode)

rprototypes.final <- rsom$prototypes%*%X

## and now the same thing for the standard vector SOM
## cut=0 is needed as the dissimilarity SOM does not provide the cutting feature
som <- batchsom(X,sg,prototypes=prototypes.init,verbose=TRUE,
                max.radius=radius.max,steps=nb.radii,mode=mode,cut=0)

stopifnot(all.equal(error.kaskilagus(som),error.kaskilagus(rsom)))
stopifnot(all.equal(error.quantisation(som),error.quantisation(rsom)))
stopifnot(all.equal(som$classif,rsom$classif))
stopifnot(all.equal(som$prototypes,rprototypes.final,
                    check.attributes = FALSE, check.names = FALSE))
