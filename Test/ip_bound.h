#ifndef IP_BOUND_H
#define IP_BOUND_H

double ell_ip_MBR(double*q,double**boundary,int dim);
double u_ip_MBR(double*q,double**boundary,int dim);
double ip_value(double*q,double*p,int dim);

double pow_term(double x,int degree);

#endif