#include "slope_intercept.h"

inline double S_poly(double x_extreme, double x_extreme_pow_d, double t, int deg_INT)
{
	double S_value;
	double t_pow_d_minus_one;
	double t_pow_d;

	t_pow_d_minus_one = pow_term(t, deg_INT - 1);
	t_pow_d = t_pow_d_minus_one * t;

	S_value = (deg_INT - 1)*t_pow_d - deg_INT * t_pow_d_minus_one*x_extreme + x_extreme_pow_d;

	return S_value;
}

inline double S_sigmoid(double x_exterme, double tanh_x_extreme, double t)
{
	double S_value = 0;
	double tanh_t;
	double tanh_t_sq;

	tanh_t = tanh(t);
	tanh_t_sq = pow_term(tanh_t, 2);
	S_value = (x_exterme - t)*tanh_t_sq - tanh_t + t + (tanh_x_extreme - x_exterme);

	return S_value;
}

double interval_Bisection(double a, double b, double S_a, double S_b, double x_extreme, double y_extreme, int deg_INT, int k_type)
{
	double middle;
	double S_middle;
	if (S_b - S_a < eps)
		return (a + b) / 2.0;

	middle = (a + b) / 2.0;

	if (k_type == 1) //Polynomial kernel
		S_middle = S_poly(x_extreme, y_extreme, middle, deg_INT);
	if (k_type == 2) //Sigmoid kernel
		S_middle = S_sigmoid(x_extreme, y_extreme, middle);

	if (S_middle < 0)
		return interval_Bisection(middle, b, S_middle, S_b, x_extreme, y_extreme, deg_INT, k_type);
	else
		return interval_Bisection(a, middle, S_a, S_middle, x_extreme, y_extreme, deg_INT, k_type);
}

//Case x_min <= x <= x_max, otherwise: no need to use rotation
void poly_rotate_up(double x_min, double x_max, double deg, double& m, double& c)
{
	double y_min;
	double t;
	double t_pow_deg;
	int deg_INT=(int)deg;
	double S_x_max;
	double x_min_pow_d;

	y_min = pow_term(x_min, deg_INT);
	if (deg_INT == 3)
	{
		t = (-1 / 2.0)*x_min;
		t_pow_deg = -y_min / 8.0;
		//t_pow_deg = (-1 / 8.0)*pow_term(x_min, deg_INT);
	}
	else //other degrees, for example: 5,7,9...
	{
		x_min_pow_d = y_min;
		S_x_max = S_poly(x_min, x_min_pow_d, x_max, deg_INT);
		
		if (S_x_max <= 0) //S(x_max) <= 0 <==> x_max < t^*
		{
			t = x_max;
			t_pow_deg = pow_term(x_max, deg_INT);
		}
		else //S(x_max) > 0 <==> x_max > t^* (Numerical method)
		{
			t = interval_Bisection(0, x_max, x_min_pow_d, S_x_max, x_min, x_min_pow_d, deg_INT, 1);
			t_pow_deg = pow_term(t, deg_INT);
		}
	}

	m = (t_pow_deg - y_min) / (t - x_min);
	c = y_min - m * x_min;
}

//Case x_min <= x <= x_max, otherwise: no need to use rotation
void poly_rotate_down(double x_min, double x_max, double deg, double& m, double& c)
{
	double y_max;
	double t;
	double t_pow_deg;
	int deg_INT = (int)deg;
	double S_x_min;
	double x_max_pow_d;

	y_max = pow_term(x_max, deg_INT);
	if (deg_INT == 3)
	{
		t = -x_max / 2.0;
		t_pow_deg = -y_max / 8.0;
	}
	else //other degrees, for example: 5,7,9...
	{
		x_max_pow_d = y_max;
		S_x_min = S_poly(x_max, x_max_pow_d, x_min, deg_INT);
		if (S_x_min <= 0) //S(x_min) <= 0 <==> x_min < t^*
		{
			t = interval_Bisection(x_min, 0, S_x_min, x_max_pow_d, x_max, x_max_pow_d, deg_INT, 1);
			t_pow_deg = pow_term(t, deg_INT);
		}
		else //S(x_min) > 0 <==> x_min > t^* (Numerical method)
		{
			t = x_min;
			t_pow_deg = pow_term(t, deg_INT);
		}
	}

	m = (y_max - t_pow_deg) / (x_max - t);
	c = y_max - m * x_max;	
}

void poly_tangent(double*q, double*a_G, double b_G, model& our_model, double deg, double& m, double& c)
{
	double t;
	double t_pow_deg_minus_one;
	double t_pow_deg;
	double ip;
	int deg_INT = (int)deg;

	ip = ip_value(q, a_G, our_model.dim);
	t = (our_model.gamma*ip + our_model.beta*b_G) / b_G;
	t_pow_deg_minus_one = (deg - 1)* pow_term(t, deg_INT - 1);
	t_pow_deg = t_pow_deg_minus_one * t;
	m = deg * t_pow_deg_minus_one;
	c = -(deg - 1)*t_pow_deg;
}

void poly_chord(double x_min, double x_max, double deg, double& m, double& c)
{
	double y_min, y_max;
	//double t;
	//double t_pow_deg;
	int deg_INT = (int)deg;

	y_min = pow_term(x_min, deg_INT);
	y_max = pow_term(x_max, deg_INT);

	m = (y_max - y_min) / (x_max - x_min);
	c = y_max - m * x_max;
}

void sigmoid_chord(double x_min, double x_max, double& m, double& c)
{
	double y_min, y_max;
	
	y_min = tanh(x_min);
	y_max = tanh(x_max);

	m = (y_max - y_min) / (x_max - x_min);
	c = y_max - m * x_max;
}

void sigmoid_tangent(double*q, double*a_G, double b_G, model& our_model, double& m, double& c)
{
	double t;
	double tanh_t;
	double ip;

	ip = ip_value(q, a_G, our_model.dim);
	t = our_model.gamma*ip / b_G + our_model.beta;

	tanh_t = tanh(t);

	m = 1 - pow_term(tanh_t, 2);
	c = tanh_t - m * t;
}

void sigmoid_rotate_up(double x_min, double x_max, double& m, double& c)
{
	double y_max;
	double S_x_min;
	double t;
	double tanh_t;

	y_max = tanh(x_max);
	S_x_min = S_sigmoid(x_max, y_max, x_min);

	if (S_x_min <= 0) //x_min>=t^*
		sigmoid_chord(x_min, x_max, m, c);
	else //x_min<=t^*
	{
		t = interval_Bisection(0, x_min, y_max - x_max, S_x_min, x_max, y_max, 0, 2);
		tanh_t = tanh(t);
		m = 1 - pow_term(tanh_t, 2);
		c = y_max - m * x_max;
	}
}

void sigmoid_rotate_down(double x_min, double x_max, double& m, double& c)
{
	double y_min;
	double S_x_max;
	double t;
	double tanh_t;

	y_min = tanh(x_min);
	S_x_max = S_sigmoid(x_min, y_min, x_max);

	if (S_x_max >= 0)
		sigmoid_chord(x_min, x_max, m, c);
	else
	{
		t = interval_Bisection(x_max, 0, S_x_max, y_min - x_min, x_min, y_min, 0, 2);
		tanh_t = tanh(t);
		m = 1 - pow_term(tanh_t, 2);
		c = y_min - m * x_min;
	}
}