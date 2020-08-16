#include "kd_tree.h"

//kdNode
double kdNode::LB(int q_id, model& our_model)
{
	double lb_sq, ub_sq;
	double lb_ip, ub_ip;
	double lb, ub;
	double L = inf;
	double U = 0;
	double*q = our_model.queryMatrix[q_id];

	if (our_model.kernel_type == "rbf")
	{
		ub_sq = u_MBR_sq(q, this->boundary_twin[0], our_model.dim);
		L = sum_w_twin[0] * exp(-our_model.gamma*ub_sq);

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		//both positive and negative weights
		lb_sq = ell_MBR_sq(q, this->boundary_twin[1], our_model.dim);
		U = sum_w_twin[1] * exp(-our_model.gamma*lb_sq);
	}

	if (our_model.kernel_type == "poly")
	{
		lb_ip = ell_ip_MBR(q, this->boundary_twin[0], our_model.dim);
		L = sum_w_twin[0] * pow_term(our_model.gamma*lb_ip + our_model.beta, (int)our_model.deg);

		if (our_model.n_d == our_model.n_d_pos)
			return L;
		
		ub_ip = u_ip_MBR(q, this->boundary_twin[1], our_model.dim);
		U = sum_w_twin[1] * pow_term(our_model.gamma*ub_ip + our_model.beta, (int)our_model.deg);
	}

	if (our_model.kernel_type == "sig")
	{
		lb_ip = ell_ip_MBR(q, this->boundary_twin[0], our_model.dim);
		L = sum_w_twin[0] * tanh(our_model.gamma*lb_ip + our_model.beta);
		if (our_model.n_d == our_model.n_d_pos)
			return L;

		ub_ip = u_ip_MBR(q, this->boundary_twin[1], our_model.dim);
		U = sum_w_twin[1] * tanh(our_model.gamma*ub_ip + our_model.beta);
	}

	if (our_model.kernel_type == "Chi2")
	{
		lb = boundary_twin[0][0][0];
		if (lb < eps || q[0] < eps)
			L = 0;
		else
			L = sum_w_twin[0] * ((2 * q[0] * lb) / (q[0] + lb));

		if (our_model.n_d == our_model.n_d_pos)
			return L;

		ub = boundary_twin[1][0][1];
		if (ub < eps || q[0] < eps)
			U = 0;
		else
			U = sum_w_twin[1] * ((2 * q[0] * ub) / (q[0] + ub));
	}

	if (our_model.kernel_type == "JS")
	{
		//code here
	}

	if (our_model.kernel_type == "Hellinger")
	{
		lb = boundary_twin[0][0][0];
		if (lb < 0)
			L = 0;
		else
			L = sqrt(q[0] * lb)*sum_w_twin[0];

		if (our_model.n_d == our_model.n_d_pos)
			return L;

		ub = boundary_twin[1][0][1];
		if (ub < 0)
			U = 0;
		else
			U = sqrt(q[0] * ub)*sum_w_twin[1];
	}
	
	return (L - U);
}

double kdNode::UB(int q_id, model& our_model)
{
	double lb_sq, ub_sq;
	double lb_ip, ub_ip;
	double lb, ub;
	double L = inf;
	double U = 0;
	double*q = our_model.queryMatrix[q_id];

	if (our_model.kernel_type == "rbf")
	{
		lb_sq = ell_MBR_sq(q, this->boundary_twin[0], our_model.dim);
		U = sum_w_twin[0] * exp(-our_model.gamma*lb_sq);

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		//both positive and negative weights
		ub_sq = u_MBR_sq(q, this->boundary_twin[1], our_model.dim);
		L = sum_w_twin[1] * exp(-our_model.gamma*ub_sq);
	}

	if (our_model.kernel_type == "poly")
	{
		ub_ip = u_ip_MBR(q, this->boundary_twin[0], our_model.dim);
		U = sum_w_twin[0] * pow_term(our_model.gamma*ub_ip + our_model.beta, (int)our_model.deg);

		if (our_model.n_d == our_model.n_d_pos)
			return U;

		lb_ip = ell_ip_MBR(q, this->boundary_twin[1], our_model.dim);
		L = sum_w_twin[1] * pow_term(our_model.gamma*lb_ip + our_model.beta, (int)our_model.deg);
	}

	if (our_model.kernel_type == "sig")
	{
		ub_ip = u_ip_MBR(q, this->boundary_twin[0], our_model.dim);
		U = sum_w_twin[0] * tanh(our_model.gamma*ub_ip + our_model.beta);

		if (our_model.n_d == our_model.n_d_pos)
			return U;

		lb_ip = ell_ip_MBR(q, this->boundary_twin[1], our_model.dim);
		L = sum_w_twin[1] * tanh(our_model.gamma*lb_ip + our_model.beta);
	}

	if (our_model.kernel_type == "Chi2")
	{
		ub = boundary_twin[0][0][1];
		if (ub < eps || q[0] < eps)
			U = 0;
		else
			U = sum_w_twin[0] * ((2 * q[0] * ub) / (q[0] + ub));

		if (our_model.n_d == our_model.n_d_pos)
			return U;

		lb = boundary_twin[1][0][0];
		if (lb < eps || q[0] < eps)
			L = 0;
		else
			L = sum_w_twin[1] * ((2 * q[0] * lb) / (q[0] + lb));
	}

	if (our_model.kernel_type == "JS")
	{
		//code here
	}

	if (our_model.kernel_type == "Hellinger")
	{
		ub = boundary_twin[0][0][1];
		if (ub < 0)
			U = 0;
		else
			U = sqrt(q[0] * ub)*sum_w_twin[0];

		if (our_model.n_d == our_model.n_d_pos)
			return U;

		lb = boundary_twin[1][0][0];
		if (lb < 0)
			L = 0;
		else
			L = sqrt(q[0] * lb)*sum_w_twin[1];
	}

	return (U - L);
}

void kdNode::update_Aug(Node*node, Tree*t)
{
	update_Node(node, t);
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

//kdNode_KARL
double kdNode_KARL::LB(int q_id, model& our_model)
{
	double t_star;
	double ip;
	double exp_minus_t_star;
	double m, c;
	double L, U;
	double x_min, x_max;
	double exp_minus_x_min, exp_minus_x_max;
	double y_min, y_max;
	double q_ell;

	if (our_model.kernel_type == "rbf")
	{
		if (sum_w_twin[0] > eps) 
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);
			gamma_sum_twin[0] = our_model.gamma*(sum_w_twin[0] * our_model.qSquareNorm - 2 * ip + b_G_twin[0]);
			t_star = gamma_sum_twin[0] / sum_w_twin[0];

			exp_minus_t_star = exp(-t_star);
			m = -exp_minus_t_star;
			c = (1 + t_star)*exp_minus_t_star;

			L = m * gamma_sum_twin[0] + c * sum_w_twin[0];
		}
		else //sum of weights is 0
			L = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			x_min = our_model.gamma*ell_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim);
			x_max = our_model.gamma*u_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim);
			if (x_max - x_min < eps)
				return (L - sum_w_twin[1] * exp(-x_min));

			exp_minus_x_min = exp(-x_min);
			exp_minus_x_max = exp(-x_max);
			m = (exp_minus_x_max - exp_minus_x_min) / (x_max - x_min);
			c = (x_max*exp_minus_x_min - x_min * exp_minus_x_max) / (x_max - x_min);

			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);
			gamma_sum_twin[1] = our_model.gamma*(sum_w_twin[1] * our_model.qSquareNorm - 2 * ip + b_G_twin[1]);

			U = m * gamma_sum_twin[1] + c * sum_w_twin[1];
		}
		else
			U = 0;

		return (L - U);
	}

	if (our_model.kernel_type == "poly")
	{
		if (sum_w_twin[0] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);

			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				L = sum_w_twin[0] * pow_term(x_min, (int)our_model.deg);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0 || our_model.is_odd_degree == false) //Use tangent line
					poly_tangent(our_model.queryMatrix[q_id], a_G_twin[0], b_G_twin[0], our_model, our_model.deg, m, c);
				if (x_max <= 0 && our_model.is_odd_degree == true)//Use chord line
					poly_chord(x_min, x_max, our_model.deg, m, c);
				if (x_min < 0 && x_max > 0 && our_model.is_odd_degree == true) //Use rotation
					poly_rotate_up(x_min, x_max, our_model.deg, m, c);

				L = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[0];
			}
		}
		else
			L = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);
			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				U = sum_w_twin[1] * pow_term(x_min, (int)our_model.deg);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0 || our_model.is_odd_degree == false) //Use chord line
					poly_chord(x_min, x_max, our_model.deg, m, c);
				if (x_max <= 0 && our_model.is_odd_degree == true)//Use tangent line
					poly_tangent(our_model.queryMatrix[q_id], a_G_twin[1], b_G_twin[1], our_model, our_model.deg, m, c);
				if (x_min <= 0 && x_max >= 0 && our_model.is_odd_degree == true) //Use rotation
					poly_rotate_down(x_min, x_max, our_model.deg, m, c);

				U = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[1];
			}
		}
		else
			U = 0;

		return (L - U);
	}

	if (our_model.kernel_type == "sig")
	{
		if (sum_w_twin[0] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);

			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				L = sum_w_twin[0] * tanh(x_min);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0) //Use chord line
					sigmoid_chord(x_min, x_max, m, c);
				if (x_max <= 0) //Use tangent line
					sigmoid_tangent(our_model.queryMatrix[q_id], a_G_twin[0], b_G_twin[0], our_model, m, c);
				if (x_min <= 0 && x_max >= 0) //Use rotation
					sigmoid_rotate_up(x_min, x_max, m, c);

				L = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[0];
			}
		}
		else
			L = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);
			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				U = sum_w_twin[1] * tanh(x_min);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0) //Use tangent line
					sigmoid_tangent(our_model.queryMatrix[q_id], a_G_twin[1], b_G_twin[1], our_model, m, c);
				if (x_max <= 0)//Use chord line
					sigmoid_chord(x_min, x_max, m, c);
				if (x_min <= 0 && x_max >= 0) //Use rotation
					sigmoid_rotate_down(x_min, x_max, m, c);

				U = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[1];
			}
		}
		else
			U = 0;

		return (L - U);
	}

	if (our_model.kernel_type == "Chi2")
	{
		q_ell = our_model.queryMatrix[q_id][0];
		if (q_ell < eps)
			return 0;

		if (sum_w_twin[0] > eps)
		{
			x_min = boundary_twin[0][0][0];
			x_max = boundary_twin[0][0][1];
			//x_min = ell_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim);
			//x_max = u_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim);
			if (x_max - x_min < eps)
				L = sum_w_twin[0] * ((2 * q_ell * x_min) / (q_ell + x_min));
			else
			{
				y_min = (2 * q_ell * x_min) / (q_ell + x_min);
				y_max = (2 * q_ell * x_max) / (q_ell + x_max);

				m = (y_max - y_min) / (x_max - x_min);
				c = y_min - m * x_min;

				L = m * a_G_twin[0][0] + c * sum_w_twin[0];
			}
		}
		else
			L = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			t_star = a_G_twin[1][0] / sum_w_twin[1];
			m = (2 * q_ell * q_ell) / ((q_ell + t_star)*(q_ell + t_star));
			c = (2 * q_ell * t_star * t_star) / ((q_ell + t_star)*(q_ell + t_star));

			U = m * a_G_twin[1][0] + c * sum_w_twin[1];
		}
		else
			U = 0;

		return (L - U);
	}

	if (our_model.kernel_type == "JS")
	{
		//code here
	}

	if (our_model.kernel_type == "Hellinger")
	{
		q_ell = our_model.queryMatrix[q_id][0];
		if (q_ell < eps)
			return 0;

		if (sum_w_twin[0] > eps)
		{
			x_min = boundary_twin[0][0][0];
			x_max = boundary_twin[0][0][1];

			if (x_max - x_min < eps)
				L = sum_w_twin[0] * sqrt(q_ell * x_min);
			else
			{
				y_min = sqrt(q_ell * x_min);
				y_max = sqrt(q_ell * x_max);

				m = (y_max - y_min) / (x_max - x_min);
				c = y_min - m * x_min;

				L = m * a_G_twin[0][0] + c * sum_w_twin[0];
			}
		}
		else
			L = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return L;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			t_star = a_G_twin[1][0] / sum_w_twin[1];
			if (t_star < eps)
				U = 0;
			else
			{
				m = 0.5*sqrt(q_ell / t_star);
				c = 0.5*sqrt(q_ell*t_star);

				U = m * a_G_twin[1][0] + c * sum_w_twin[1];
			}
		}
		else
			U = 0;

		return (L - U);
	}
	
	return inf;
}

double kdNode_KARL::UB(int q_id, model& our_model)
{
	double t_star;
	double ip;
	double exp_minus_t_star;
	double m, c;
	double L, U;
	double x_min, x_max;
	double exp_minus_x_min, exp_minus_x_max;
	double y_min, y_max;
	double q_ell;

	if (our_model.kernel_type == "rbf")
	{
		if (sum_w_twin[0] > eps)
		{
			x_min = our_model.gamma*ell_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim);
			x_max = our_model.gamma*u_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim);
			if (x_max - x_min < eps)
				U = sum_w_twin[0] * exp(-x_min);
			else
			{
				exp_minus_x_min = exp(-x_min);
				exp_minus_x_max = exp(-x_max);
				m = (exp_minus_x_max - exp_minus_x_min) / (x_max - x_min);
				c = (x_max*exp_minus_x_min - x_min * exp_minus_x_max) / (x_max - x_min);

				ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);
				gamma_sum_twin[0] = our_model.gamma*(sum_w_twin[0] * our_model.qSquareNorm - 2 * ip + b_G_twin[0]);

				U = m * gamma_sum_twin[0] + c * sum_w_twin[0];
			}
		}
		else
			U = 0;
		
		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		if (sum_w_twin[1] > eps) //both positive and negative weights
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);
			gamma_sum_twin[1] = our_model.gamma*(sum_w_twin[1] * our_model.qSquareNorm - 2 * ip + b_G_twin[1]);
			t_star = gamma_sum_twin[1] / sum_w_twin[1];

			exp_minus_t_star = exp(-t_star);
			m = -exp_minus_t_star;
			c = (1 + t_star)*exp_minus_t_star;

			L = m * gamma_sum_twin[1] + c * sum_w_twin[1];
		}
		else
			L = 0;

		return (U - L);
	}

	if (our_model.kernel_type == "poly")
	{
		if (sum_w_twin[0] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);
			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				U = sum_w_twin[0] * pow_term(x_min, (int)our_model.deg);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0 || our_model.is_odd_degree == false) //Use chord line
					poly_chord(x_min, x_max, our_model.deg, m, c);
				if (x_max <= 0 && our_model.is_odd_degree == true)//Use tangent line
					poly_tangent(our_model.queryMatrix[q_id], a_G_twin[0], b_G_twin[0], our_model, our_model.deg, m, c);
				if (x_min <= 0 && x_max >= 0 && our_model.is_odd_degree == true) //Use rotation
					poly_rotate_down(x_min, x_max, our_model.deg, m, c);

				U = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[0];
			}
		}
		else
			U = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		if (sum_w_twin[1] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);

			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				L = sum_w_twin[1] * pow_term(x_min, (int)our_model.deg);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0 || our_model.is_odd_degree == false) //Use tangent line
					poly_tangent(our_model.queryMatrix[q_id], a_G_twin[1], b_G_twin[1], our_model, our_model.deg, m, c);
				if (x_max <= 0 && our_model.is_odd_degree == true)//Use chord line
					poly_chord(x_min, x_max, our_model.deg, m, c);
				if (x_min < 0 && x_max > 0 && our_model.is_odd_degree == true) //Use rotation
					poly_rotate_up(x_min, x_max, our_model.deg, m, c);

				L = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[1];
			}
		}
		else
			L = 0;

		return (U - L);
	}

	if (our_model.kernel_type == "sig")
	{
		if (sum_w_twin[0] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[0], our_model.dim);
			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[0], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				U = sum_w_twin[0] * tanh(x_min);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0) //Use tangent line
					sigmoid_tangent(our_model.queryMatrix[q_id], a_G_twin[0], b_G_twin[0], our_model, m, c);
				if (x_max <= 0)//Use chord line
					sigmoid_chord(x_min, x_max, m, c);
				if (x_min <= 0 && x_max >= 0) //Use rotation
					sigmoid_rotate_down(x_min, x_max, m, c);

				U = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[0];
			}
		}
		else
			U = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		if (sum_w_twin[1] > eps)
		{
			ip = ip_value(our_model.queryMatrix[q_id], a_G_twin[1], our_model.dim);

			x_min = our_model.gamma*ell_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;
			x_max = our_model.gamma*u_ip_MBR(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim) + our_model.beta;

			if (x_max - x_min < eps)
				L = sum_w_twin[1] * tanh(x_min);
			else
			{
				//obtain the slope m and the intercept c
				if (x_min >= 0) //Use tangent line
					sigmoid_chord(x_min, x_max, m, c);
				if (x_max <= 0)//Use chord line
					sigmoid_tangent(our_model.queryMatrix[q_id], a_G_twin[1], b_G_twin[1], our_model, m, c);
				if (x_min < 0 && x_max > 0) //Use rotation
					sigmoid_rotate_up(x_min, x_max, m, c);

				L = m * our_model.gamma*ip + (m*our_model.beta + c)*sum_w_twin[1];
			}
		}
		else
			L = 0;

		return (U - L);
	}

	if (our_model.kernel_type == "Chi2")
	{
		q_ell = our_model.queryMatrix[q_id][0];
		if (q_ell < eps)
			return 0;

		if (sum_w_twin[0] > eps)
		{
			t_star = a_G_twin[0][0] / sum_w_twin[0];
			m = (2 * q_ell * q_ell) / ((q_ell + t_star)*(q_ell + t_star));
			c = (2 * q_ell * t_star * t_star) / ((q_ell + t_star)*(q_ell + t_star));

			U = m * a_G_twin[0][0] + c * sum_w_twin[0];
		}
		else
			U = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		if (sum_w_twin[1] > eps)
		{
			x_min = boundary_twin[1][0][0];
			x_max = boundary_twin[1][0][1];
			//x_min = ell_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim);
			//x_max = u_MBR_sq(our_model.queryMatrix[q_id], boundary_twin[1], our_model.dim);
			if (x_max - x_min < eps)
				return (U - sum_w_twin[1] * ((2 * q_ell * x_min) / (q_ell + x_min)));

			y_min = (2 * q_ell * x_min) / (q_ell + x_min);
			y_max = (2 * q_ell * x_max) / (q_ell + x_max);

			m = (y_max - y_min) / (x_max - x_min);
			c = y_min - m * x_min;

			L = m * a_G_twin[1][0] + c * sum_w_twin[1];
		}
		else
			L = 0;

		return (U - L);
	}

	if (our_model.kernel_type == "JS")
	{
		//code here
	}

	if (our_model.kernel_type == "Hellinger")
	{
		q_ell = our_model.queryMatrix[q_id][0];
		if (q_ell < eps)
			return 0;

		if (sum_w_twin[0] > eps)
		{
			t_star = a_G_twin[0][0] / sum_w_twin[0];

			if (t_star < eps)
				U = 0;
			else
			{
				m = 0.5*sqrt(q_ell / t_star);
				c = 0.5*sqrt(q_ell*t_star);

				U = m * a_G_twin[0][0] + c * sum_w_twin[0];
			}
		}
		else
			U = 0;

		if (our_model.n_d == our_model.n_d_pos) //only positive weights
			return U;

		if (sum_w_twin[1] > eps)
		{
			x_min = boundary_twin[1][0][0];
			x_max = boundary_twin[1][0][1];
			
			if (x_max - x_min < eps)
				return (U - sum_w_twin[1] * sqrt(q_ell*x_min));

			y_min = sqrt(q_ell*x_min);
			y_max = sqrt(q_ell*x_max);

			m = (y_max - y_min) / (x_max - x_min);
			c = y_min - m * x_min;

			L = m * a_G_twin[1][0] + c * sum_w_twin[1];
		}
		else
			L = 0;

		return (U - L);
	}

	return -inf;
}

inline void update_node_twin_variables(int id, kdNode_KARL*node, model& our_model, int pn)
{
	double norm_square = 0;
	for (int d = 0; d < our_model.dim; d++)
	{
		node->a_G_twin[pn][d] += our_model.weightVector[id] * our_model.dataMatrix[id][d];
		norm_square += our_model.dataMatrix[id][d] * our_model.dataMatrix[id][d];
	}
	
	node->b_G_twin[pn] += our_model.weightVector[id] * norm_square;
}

void kdNode_KARL::update_KARL_Info(kdNode_KARL*node, Tree*t)
{
	int id;
	int pn;

	//obtain vec a_G and constant b_G
	node->a_G_twin[0] = new double[t->our_model.dim];
	node->a_G_twin[1] = new double[t->our_model.dim];
	for (int d = 0; d < t->our_model.dim; d++)
	{
		node->a_G_twin[0][d] = 0;
		node->a_G_twin[1][d] = 0;
	}

	node->b_G_twin[0] = 0;
	node->b_G_twin[1] = 0;

	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];

		if (id < t->our_model.n_d_pos)
			pn = 0;
		else
			pn = 1;
		//if (t->our_model.weightVector[id] > 0)
		//	pn = 0;
		//else
		//	pn = 1;

		update_node_twin_variables(id, node, t->our_model, pn);
	}

	if ((int)node->idList.size() <= t->our_model.leafCapacity) //this is the leaf node
		return;

	update_KARL_Info((kdNode_KARL*)node->childVector[0], t);
	update_KARL_Info((kdNode_KARL*)node->childVector[1], t);
}

void kdNode_KARL::update_Aug(Node*node, Tree*t)
{
	this->update_Node((kdNode_KARL*)node, t);
	this->update_KARL_Info((kdNode_KARL*)node, t);
}

kdNode_KARL*kdNode_KARL::createNode()
{
	return new kdNode_KARL();
}

//kdTree
kdTree::kdTree(model& our_model)
{
	this->our_model = our_model;
}

double kdTree::obtain_SplitValue(kdNode*node, int split_Dim)
{
	vector<double> tempVector;
	int id;
	int middle_left, middle_right, middle;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		tempVector.push_back(our_model.dataMatrix[id][split_Dim]);
	}

	sort(tempVector.begin(), tempVector.end());

	if ((int)tempVector.size() % 2 == 0)//even number
	{
		middle_right = (int)tempVector.size() / 2;
		middle_left = middle_right - 1;

		return ((tempVector[middle_left] + tempVector[middle_right]) / 2.0);
	}
	else
	{
		middle = ((int)tempVector.size() - 1) / 2;
		return tempVector[middle];
	}

	tempVector.clear();
}

void kdTree::KD_Tree_Recur(kdNode*node, int split_Dim)
{
	int id;
	int counter;
	//base case
	
	if ((int)node->idList.size() <= our_model.leafCapacity)
		return;

	//split_Dim=select_split_Dim(node);
	double splitValue = obtain_SplitValue(node, split_Dim); //code here

	//create two children
	kdNode*leftNode;
	kdNode*rightNode;

	leftNode = node->createNode();
	rightNode = node->createNode();

	counter = 0;
	int halfSize = ((int)node->idList.size()) / 2;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		if (our_model.dataMatrix[id][split_Dim] <= splitValue && counter <= halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	//getNode_Boundary(leftNode);
	//getNode_Boundary(rightNode);

	KD_Tree_Recur(leftNode, (split_Dim + 1) % our_model.dim);
	KD_Tree_Recur(rightNode, (split_Dim + 1) % our_model.dim);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree::build_kdTree()
{
	for (int i = 0; i < our_model.n_d; i++)
		rootNode->idList.push_back(i);

	KD_Tree_Recur((kdNode*)rootNode, 0);
	updateAugment((kdNode*)rootNode);
}

void kdTree::updateAugment(kdNode*node)
{
	node->update_Aug(node, this);
}