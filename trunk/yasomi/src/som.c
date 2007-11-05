#include <R.h>
#include <Rmath.h>
#include "bmu.h"


/* Neighboorhood function calculation. The kernel functions are defined such
   that neighborhood(x,y) is approximatly zero when the distance between x and
   y is higher than the radius parameter.  
   When normalized, a row major approach is used to preserve memory locality. */

void neighborhood(double *distances,double *nv,int nbUnit,double radius,
		  int kernelType,int isNormalized) {
  int i,j,base;
  double tmp;
  int size = nbUnit*nbUnit;

  switch(kernelType) {
  case 0:
    /* Gaussian kernel */
    radius /= 3;
    for(i = 0 ; i < size; i++) {
      tmp = distances[i]/radius;
      nv[i] = exp (-0.5*tmp*tmp);
    }
    break;
  case 1:
    /* Linear kernel */
    for(i = 0 ; i < size; i++) {
      tmp = distances[i];
      nv[i] = tmp < radius ? (1-tmp/radius) : 0;
    }
    break;
  }
  if(isNormalized) {
    for(i = 0 ; i < nbUnit; i++) {
      tmp = 0;
      base = i*nbUnit;
      for(j = 0 ; j < nbUnit; j++) {
	tmp += nv[base+j];
      }
      for(j = 0 ; j < nbUnit; j++) {
	nv[base+j] /= tmp;
      }
    }
  }
}

/* proto (nproto,dim): initial values of the prototypes
   data (ndata,dim): dataset
   assign: assignment rule
   grid (nproto,nproto): distances between units (square matrix)
   radii (nradii): annealing scheme
   maxiter: maximal number of iterations per radius
   kernel: the kernel type
   normalized: normalization flag for the kernel
   cut: cuting value for updates
   verbose: verbose flag
   classif (ndata): clustering result
   errors (nradii,maxiter): evolution of the quantization error
 */

void batch_som(double *proto,int *nproto,double *data,int *ndata,
	       int *dim,int *assign,double *grid,double *radii,
	       int *nradii,int *maxiter,int *kernel,int *normalized,
	       double *cut,int *verbose,
	       int *classif,double *errors) {
  int dataSize=*ndata,protoSize=*nproto,dimension=*dim;
  int nbRadii=*nradii,maxIter=*maxiter;
  double cutValue=*cut;
  int kernelType=*kernel;
  int isNormalized=*normalized;
  int isVerbose=*verbose;

  int iRad,iter,indiv,j,k,l;
  int base;
  double *nv;
  double *preProto;
  double *tmpProto;

  int changed;
  double bestDist,dist,tmp,norm,neigh;
  int bestSoFar;

  double totalError;

  /* for heskes' rule */
  int assignmentRule=*assign;
  double *distances;

  if(assignmentRule == 1) {
    distances = (double *) R_alloc(protoSize, sizeof(double));
  }

  /* neighboorhood function */
  nv = (double *) R_alloc(protoSize*protoSize, sizeof(double));

  /* temporary prototype */
  tmpProto = (double *) R_alloc(dimension, sizeof(double));

  if(isVerbose) {
    Rprintf("Using");
    if(isNormalized) {
      Rprintf(" normalized");
    } else {
      Rprintf(" unnormalized");
    }
    switch(kernelType) {
    case 0:
      Rprintf(" Gaussian kernel");
      break;
    case 1:
      Rprintf(" linear kernel");
      break;
    }
    switch(assignmentRule) {
    case 0:
      Rprintf(" with simple BMU determination");
      break;
    case 1:
      Rprintf(" with Heskes' assignment rule");
      break;
    }
    Rprintf("\n");
  }

  /* annealing loop */
  for(iRad = 0; iRad < nbRadii; iRad++) {
    /* compute neighboorhood */
    neighborhood(grid,nv,protoSize,radii[iRad],kernelType,isNormalized);
    for(iter = 0; iter < maxIter; iter++) {
      /* assignement */
      if(assignmentRule == 1) {
	/* Heskes' rule */
	changed = bmu_heskes_ext_mem(proto,nv,nproto,data,ndata,dim,classif,
				     &totalError,distances);
      } else {
	changed = bmu(proto,nproto,data,ndata,dim,classif,&totalError);
      }
      errors[iter + iRad*maxIter] = totalError;
      if(isVerbose) {
	Rprintf("%i (%i) %g\n",iRad,iter,totalError);
      }
      /* representation */
      if(changed || iter == 0) {
	/* then compute prototypes */
	for(j = 0; j < protoSize; j++) {
	  base = j * protoSize;
	  norm = 0;
	  memset(tmpProto,0,dimension*sizeof(double));
	  for(indiv = 0 ; indiv < dataSize ; indiv++) {
	    l = classif[indiv];
	    neigh = nv[ base + l ];
	    norm += neigh;
	    /* loop on dimensions */
	    for(k = 0; k < dimension; k++) {
	      tmpProto[k] += neigh * data[indiv + k * dataSize];
	    }
	  }
	  /* normalization and storage */
	  if(norm>cutValue) {
	    for(k = 0; k < dimension; k++) {
	      proto[j + k * protoSize]=tmpProto[k]/norm;
	    }
	  } 
	}
      } else if(iter>=1) {
	if(isVerbose) {
	  Rprintf("Iteration %i for radius %g (%i) is stable, decreasing radius\n",iter,radii[iRad],iRad);
	}
	break;
      }
    }
  }
}

/* same calculation as batch_som but slightly optimized by first calculating 
   the center of mass of each cluster and then the updated prototypes */

void batch_som_optim(double *proto,int *nproto,double *data,int *ndata,
		     int *dim,int *assign,double *grid,double *radii,
		     int *nradii,int *maxiter,int *kernel,int *normalized,
		     double *cut,int *verbose,
		     int *classif,double *errors) {
  int dataSize=*ndata,protoSize=*nproto,dimension=*dim;
  int nbRadii=*nradii,maxIter=*maxiter;
  double cutValue=*cut;
  int kernelType=*kernel;
  int isNormalized=*normalized;
  int isVerbose=*verbose;

  int iRad,iter,indiv,j,k,l;
  int base,protoBase;
  double *nv;
  double *preProto;
  int *clusterSize;
  double *tmpProto;

  int changed;
  double bestDist,dist,tmp,norm,neigh;
  int bestSoFar;

  double totalError;

  /* for heskes' rule */
  int assignmentRule=*assign;
  double *distances;

  if(assignmentRule == 1) {
    distances = (double *) R_alloc(protoSize, sizeof(double));
  }

  /* neighboorhood function */
  nv = (double *) R_alloc(protoSize*protoSize, sizeof(double));

  /* class based prototypes */
  /* stored in row major mode */
  preProto = (double *) R_alloc(protoSize*dimension, sizeof(double));

  /* clusterSize */
  clusterSize = (int *) R_alloc(protoSize, sizeof(int));

  /* temporary prototype */
  tmpProto = (double *) R_alloc(dimension, sizeof(double));

  if(isVerbose) {
    Rprintf("Using");
    if(isNormalized) {
      Rprintf(" normalized");
    } else {
      Rprintf(" unnormalized");
    }
    switch(kernelType) {
    case 0:
      Rprintf(" Gaussian kernel");
      break;
    case 1:
      Rprintf(" linear kernel");
      break;
    }
    switch(assignmentRule) {
    case 0:
      Rprintf(" with simple BMU determination");
      break;
    case 1:
      Rprintf(" with Heskes' assignment rule");
      break;
    }
    Rprintf("\n");
  }

  /* annealing loop */
  for(iRad = 0; iRad < nbRadii; iRad++) {
    /* compute neighboorhood */
    neighborhood(grid,nv,protoSize,radii[iRad],kernelType,isNormalized);
    for(iter = 0; iter < maxIter; iter++) {
      /* assignement */
      if(assignmentRule == 1) {
	/* Heskes' rule */
	changed = bmu_heskes_ext_mem(proto,nv,nproto,data,ndata,dim,classif,
				     &totalError,distances);
      } else {
	changed = bmu(proto,nproto,data,ndata,dim,classif,&totalError);
      }
      errors[iter + iRad*maxIter] = totalError;
      if(isVerbose) {
	Rprintf("%i (%i) %g\n",iRad,iter,totalError);
      }
      /* representation */
      if(changed || iter == 0) {
	/* first compute class center of mass */
	/* zero */
	memset(preProto,0,protoSize*dimension*sizeof(double));
	memset(clusterSize,0,protoSize*sizeof(int));
	for(indiv = 0 ; indiv < dataSize ; indiv++) {
	  j = classif[indiv];
	  clusterSize[j]++;
	  /* loop on dimensions */
	  protoBase = j * dimension;
	  for(k = 0; k < dimension; k++) {
	    preProto[protoBase + k] += data[indiv + k * dataSize];
	  }
	}
	/* then compute prototypes */
	for(j = 0; j < protoSize; j++) {
	  norm = 0;
	  /* first compute normalization factor */
	  base = j * protoSize;
	  for(l = 0; l < protoSize; l++) {
	    if(clusterSize[l]>0) {
	      norm += nv[base + l] * clusterSize[l];
	    }
	  }
	  if(norm>cutValue) {
	    /* clear old value */
	    memset(tmpProto,0,dimension*sizeof(double));
	    /* compute new value */
	    for(l = 0; l < protoSize; l++) {
	      if(clusterSize[l]>0) {
		neigh = nv[base + l];
		protoBase = l * dimension;
		/* loop on dimensions */
		for(k = 0; k < dimension; k++) {
		  tmpProto[k] += neigh * preProto[protoBase + k];
		}
	      }
	    }
	    /* normalization and transfer */
	    for(k = 0; k < dimension; k++) {
	      proto[j + k * protoSize] = tmpProto[k]/norm;
	    }
	  }
	}
      } else if(iter>=1) {
	if(isVerbose) {
	  Rprintf("Iteration %i for radius %g (%i) is stable, decreasing radius\n",iter,radii[iRad],iRad);
	}
	break;
      }
    }
  }
}

