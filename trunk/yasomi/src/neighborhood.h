#ifndef YASOMI_NEIGHBORHOOD_H
#define YASOMI_NEIGHBORHOOD__H

typedef enum {
    gaussian = 0,
    linear = 1
} kernel_type;

void neighborhood(double *distances,double *nv,int nbUnit,double radius,
		  int kernelType,int isNormalized);

void neighborhood_single(double *distances,double *nv,int *nbUnit,
			 double *radius,int *kernelType,int *isNormalized);

#endif /* !YASOMI_NEIGHBORHOOD_H */
