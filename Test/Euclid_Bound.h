#ifndef EUCLID_BOUND_H
#define EUCLID_BOUND_H

#include "Library.h"

//No extern--> cannot compile
double ell_MBR(double*q,double**boundary,int dim);
double u_MBR(double*q,double**boundary,int dim);
double ell_MBR_sq(double*q, double**boundary, int dim);
double u_MBR_sq(double*q, double**boundary, int dim);
double u_tri(double*q,double*center,int dim,double radius,double& obt_dist);
double euclid_dist(double*q, double*p, int dim);
double euclid_dist_sq(double*q, double*p, int dim);

//Used in dual-group setting
double ell_MBR_MBR(double**Q_boundary, double**P_boundary, int dim);
double u_MBR_MBR(double**Q_boundary, double**P_boundary, int dim);

#endif