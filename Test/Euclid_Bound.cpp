#include "Euclid_Bound.h"

double ell_MBR(double*q,double**boundary,int dim)
{
	double ell2;
	double ell;

	ell2 = ell_MBR_sq(q, boundary, dim);
	ell = sqrt(ell2);

	return ell;
}

double u_MBR(double*q,double**boundary,int dim)
{
	double u2;
	double u;

	u2 = u_MBR_sq(q, boundary, dim);
	u = sqrt(u2);

	return u;
}

double ell_MBR_sq(double*q, double**boundary, int dim)
{
	double ell_sq = 0;
	double temp;
	for (int d = 0; d < dim; d++)
	{
		if (q[d] < boundary[d][0])
			temp = boundary[d][0] - q[d];
		else
		{
			if (q[d] > boundary[d][1])
				temp = q[d] - boundary[d][1];
			else
				temp = 0;
		}

		ell_sq += temp * temp;
	}

	return ell_sq;
}

double u_MBR_sq(double*q, double**boundary, int dim)
{
	double u_sq = 0;
	double temp;
	for (int d = 0; d < dim; d++)
	{
		temp = max(fabs(q[d] - boundary[d][0]), fabs(q[d] - boundary[d][1]));
		u_sq += temp * temp;
	}

	return u_sq;
}

double u_tri(double*q,double*center,int dim,double radius,double& obt_dist)
{
	double u;
	obt_dist=0;

	for(int d=0;d<dim;d++)
		obt_dist+=(q[d]-center[d])*(q[d]-center[d]);

	obt_dist=sqrt(obt_dist);

	u=obt_dist+radius;

	return u;

}

double euclid_dist(double*q, double*p, int dim)
{
	double dist;
	dist = euclid_dist_sq(q, p, dim);

	return sqrt(dist);
}

double euclid_dist_sq(double*q, double*p, int dim)
{
	double dist_sq = 0;
	for (int d = 0; d < dim; d++)
		dist_sq += (q[d] - p[d])*(q[d] - p[d]);

	return dist_sq;
}

double ell_MBR_MBR(double**Q_boundary, double**P_boundary, int dim)
{
	double MIN_Q_P;
	double ell;

	ell = 0;
	for (int d = 0; d < dim; d++)
	{
		if (Q_boundary[d][0] > P_boundary[d][1])
			MIN_Q_P = Q_boundary[d][0] - P_boundary[d][1];
		else
		{
			if (P_boundary[d][0] > Q_boundary[d][1])
				MIN_Q_P = P_boundary[d][0] - Q_boundary[d][1];
			else
				MIN_Q_P = 0;
		}
		ell += MIN_Q_P * MIN_Q_P;
	}

	return sqrt(ell);
}

double u_MBR_MBR(double**Q_boundary, double**P_boundary, int dim)
{
	double MAX_Q_P;
	double u;

	u = 0;

	for (int d = 0; d < dim; d++)
	{
		MAX_Q_P = max(fabs(Q_boundary[d][1] - P_boundary[d][0]), fabs(P_boundary[d][1] - Q_boundary[d][0]));
		u += MAX_Q_P * MAX_Q_P;
	}

	return sqrt(u);
}