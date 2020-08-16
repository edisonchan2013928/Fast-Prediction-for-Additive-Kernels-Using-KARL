#include "Validation.h"

bool validate_best(double L,double U,double rel_error,double& val_R)
{
	if (U - L <= rel_error * fabs(U + L))
	{
		val_R = 2 * L*U / (L + U);
		return true;
	}

	if (fabs(L) < eps && fabs(U) < eps)
	{
		val_R = 0;
		return true;
	}

	return false;
}