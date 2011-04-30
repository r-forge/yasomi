#include <R.h>
#include <Rmath.h>
#include "neighborhood.h"

/* Neighborhood function calculation. The kernel functions are defined such
   that neighborhood(x,y) is approximatly zero when the distance between x and
   y is higher than the radius parameter.  
   When normalized, a row major approach is used to preserve memory locality. */

void neighborhood(double *distances,double *nv,int nbUnit,double radius,
		  int kernelType,int isNormalized) 
{
    int i,j,base;
    double tmp;
    int size = nbUnit*nbUnit;
    kernel_type eKernelType = (kernel_type) kernelType;
    
    switch(eKernelType) {
    case gaussian:
	/* Gaussian kernel */
	radius /= 3;
	for(i = 0 ; i < size; i++) {
	    tmp = distances[i]/radius;
	    nv[i] = exp (-0.5*tmp*tmp);
	}
	break;
    case linear:
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

void neighborhood_single(double *distances,double *nv,int *nbUnit,
			 double *radius,int *kernelType,int *isNormalized) 
{
    int i;
    double tmp;
    kernel_type eKernelType = (kernel_type) *kernelType;
    int nb=*nbUnit;
    double theRadius=*radius;
    
    switch(eKernelType) {
    case gaussian:
	/* Gaussian kernel */
	theRadius /= 3;
	for(i = 0 ; i < nb; i++) {
	    tmp = distances[i]/theRadius;
	    nv[i] = exp (-0.5*tmp*tmp);
	}
	break;
    case linear:
	/* Linear kernel */
	for(i = 0 ; i < nb; i++) {
	    tmp = distances[i];
	    nv[i] = tmp < theRadius ? (1-tmp/theRadius) : 0;
	}
	break;
    }
    if(*isNormalized) {
	tmp = 0;
	for(i = 0 ; i < nb; i++) {
	    tmp += nv[i];
	}
	for(i = 0 ; i < nb; i++) {
	    nv[i] /= tmp;
	}
    }
}
