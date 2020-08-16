#include "SS.h"

bool compare_diff_SCAN(const diff_SCAN& lhs, const diff_SCAN& rhs) {
	return lhs.bound_diff > rhs.bound_diff;
}

double refine(double*q, Node*curNode, model& our_model)
{
	double incr_Value = 0;
	double sq_Value;
	double ip;
	int id;
	int node_size = (int)curNode->idList.size();

	for (int i = 0; i < node_size; i++)
	{
		id = curNode->idList[i];

		if (our_model.kernel_type == "Chi2" || our_model.kernel_type == "JS" || our_model.kernel_type == "Hellinger")
			incr_Value += our_model.weight_oriVector[id] * basic_function(q, our_model.dataMatrix[id], our_model);

		if (our_model.kernel_type == "rbf")
		{
			sq_Value = euclid_dist_sq(q, our_model.dataMatrix[id], our_model.dim);
			incr_Value += our_model.weight_oriVector[id] * exp(-our_model.gamma*sq_Value);
		}

		if (our_model.kernel_type == "poly" || our_model.kernel_type == "sig")
		{
			ip = ip_value(q, our_model.dataMatrix[id], our_model.dim);
			if (our_model.kernel_type == "poly")
				incr_Value += our_model.weight_oriVector[id] * pow_term(our_model.gamma*ip + our_model.beta, (int)our_model.deg);
			if (our_model.kernel_type == "sig")
				incr_Value += our_model.weight_oriVector[id] * tanh(our_model.gamma*ip + our_model.beta);
		}
	}

	return incr_Value;
}

void SS_iter(int q_id, model& our_model)
{
	double incr_Value = 0;
	double sq_Value;
	double ip;

	for (int i = 0; i < our_model.n_d; i++)
	{
		if (our_model.kernel_type == "Chi2" || our_model.kernel_type == "JS" || our_model.kernel_type == "Hellinger")
			incr_Value += our_model.weight_oriVector[i] * basic_function(our_model.queryMatrix[q_id], our_model.dataMatrix[i], our_model);

		if (our_model.kernel_type == "rbf")
		{
			sq_Value = euclid_dist_sq(our_model.queryMatrix[q_id], our_model.dataMatrix[i], our_model.dim);
			incr_Value += our_model.weight_oriVector[i] * exp(-our_model.gamma*sq_Value);
		}
		if (our_model.kernel_type == "poly" || our_model.kernel_type == "sig")
		{
			ip = ip_value(our_model.queryMatrix[q_id], our_model.dataMatrix[i], our_model.dim);
			if (our_model.kernel_type == "poly")
				incr_Value += our_model.weight_oriVector[i] * pow_term(our_model.gamma*ip + our_model.beta, (int)our_model.deg);
			if (our_model.kernel_type == "sig")
				incr_Value += our_model.weight_oriVector[i] * tanh(our_model.gamma*ip + our_model.beta);
		}
	}

	if (our_model.is_tau == 0)
		our_model.resultVector[q_id] = incr_Value;
	else
	{
		if (incr_Value > our_model.tau)
			our_model.resultVector[q_id] = 1;
		else
			our_model.resultVector[q_id] = -1;
	}
}

double basic_function(double*q, double*p, model& our_model)
{
	double value = 0;

	for (int d = 0; d < our_model.dim; d++)
	{
		if (our_model.kernel_type == "Chi2")
		{
			if (fabs(q[d]) > eps && fabs(p[d]) > eps)
				value += ((q[d] * p[d]) / (q[d] + p[d]));
		}
		if (our_model.kernel_type == "JS")
		{
			if (fabs(q[d]) > eps && fabs(p[d]) > eps)
				value += (q[d] * log2((q[d] + p[d]) / q[d]) + p[d] * log2((q[d] + p[d]) / p[d]));
		}
		if (our_model.kernel_type == "Hellinger")
			value += sqrt(q[d] * p[d]);
	}

	if (our_model.kernel_type == "Chi2")
		value *= 2;
	if (our_model.kernel_type == "JS")
		value *= 0.5;

	return value;
}

//Only for addictive kernel functions
void SS_one_D_additive(double q, int d, model& our_model, double& incr_Value_Pos, double& incr_Value_Neg)
{
	double**dataMatrix = our_model.dataMatrix;
	incr_Value_Pos = 0;
	incr_Value_Neg = 0;
	if (our_model.kernel_type == "Chi2")
	{
		for (int i = 0; i < our_model.n_d; i++)
		{
			if (fabs(q) > eps && fabs(dataMatrix[i][d]) > eps)
			{
				if (i < our_model.n_d_pos)
					incr_Value_Pos += our_model.weightVector[i] * (dataMatrix[i][d] / (q + dataMatrix[i][d]));
				else
					incr_Value_Neg += our_model.weightVector[i] * (dataMatrix[i][d] / (q + dataMatrix[i][d]));
			}
		}
		incr_Value_Pos *= (2 * q);
		incr_Value_Neg *= (2 * q);
	}

	if (our_model.kernel_type == "JS")
	{
		for (int i = 0; i < our_model.n_d; i++)
		{
			if (fabs(q) > eps && fabs(dataMatrix[i][d]) > eps)
			{
				if (i < our_model.n_d_pos)
					incr_Value_Pos += our_model.weightVector[i] * (q * log2((q + dataMatrix[i][d]) / q) + dataMatrix[i][d] * log2((q + dataMatrix[i][d]) / dataMatrix[i][d]));
				else
					incr_Value_Neg += our_model.weightVector[i] * (q * log2((q + dataMatrix[i][d]) / q) + dataMatrix[i][d] * log2((q + dataMatrix[i][d]) / dataMatrix[i][d]));
			}
		}
		incr_Value_Pos *= 0.5;
		incr_Value_Neg *= 0.5;
	}

	if (our_model.kernel_type == "Hellinger")
	{
		for (int i = 0; i < our_model.n_d; i++)
		{
			if (i < our_model.n_d_pos)
				incr_Value_Pos += our_model.weightVector[i] * sqrt(dataMatrix[i][d]);
			else
				incr_Value_Neg += our_model.weightVector[i] * sqrt(dataMatrix[i][d]);
		}
		incr_Value_Pos *= sqrt(q);
		incr_Value_Neg *= sqrt(q);
	}
}

bool check_termination(model& our_model, int q_id, double LB, double UB)
{
	if (our_model.is_tau == 0) //approximate-KAQ
	{
		double val_R;
		if (validate_best(LB, UB, our_model.epsilon, val_R) == true)
		{
			our_model.resultVector[q_id] = val_R;
			return true;
		}
	}
	else
	{
		if (LB > our_model.tau)
		{
			our_model.resultVector[q_id] = 1;
			return true;
		}
		if (UB < our_model.tau)
		{
			our_model.resultVector[q_id] = -1;
			return true;
		}
	}

	return false;
}

void sort_SS_iter(int q_id, double LB, double UB, vector<double>& LB_d, vector<double>& UB_d, model& our_model)
{
	vector<diff_SCAN> bound_diff_vec;
	diff_SCAN temp;
	double agg_pos, agg_neg;

	for (int d = 0; d < our_model.dim; d++)
	{
		temp.d = d;
		temp.bound_diff = UB_d[d] - LB_d[d];
		bound_diff_vec.push_back(temp);
	}

	sort(bound_diff_vec.begin(), bound_diff_vec.end(), compare_diff_SCAN);

	for (int d = 0; d < our_model.dim; d++)
	{
		LB -= LB_d[bound_diff_vec[d].d];
		UB -= UB_d[bound_diff_vec[d].d];

		SS_one_D_additive(our_model.queryMatrix[q_id][bound_diff_vec[d].d], bound_diff_vec[d].d, our_model, agg_pos, agg_neg);
		LB += (agg_pos - agg_neg);
		UB += (agg_pos - agg_neg);

		if (check_termination(our_model, q_id, LB, UB) == true)
			return;
	}
}