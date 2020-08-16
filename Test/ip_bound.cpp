#include "ip_bound.h"

double ell_ip_MBR(double*q,double**boundary,int dim)
{
	double ell=0;
	for(int d=0;d<dim;d++)
	{
		if(q[d]>=0)
			ell+=q[d]*boundary[d][0];
		else //q[d]<0
			ell+=q[d]*boundary[d][1];
	}

	return ell;
}

double u_ip_MBR(double*q,double**boundary,int dim)
{
	double u=0;

	for(int d=0;d<dim;d++)
	{
		if(q[d]>=0)
			u+=q[d]*boundary[d][1];
		else //q[d]<0
			u+=q[d]*boundary[d][0];
	}

	return u;
}

double ip_value(double*q,double*p,int dim)
{
	double value=0;
	for(int d=0;d<dim;d++)
		value+=q[d]*p[d];

	return value;
}

double pow_term(double x,int degree)
{
	double value=1;

	for(int d=0;d<degree;d++)
		value=value*x;

	return value;
}