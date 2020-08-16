#pragma once
#ifndef SS_H
#define SS_H

#include "init_KARL.h"
#include "Tree.h"
#include "Euclid_Bound.h"
#include "Validation.h"

struct diff_SCAN
{
	int d;
	double bound_diff;

	//bool operator<(const diff_SCAN& rhs) const { return bound_diff < rhs.bound_diff; }
};

double refine(double*q, Node*curNode, model& our_model);
void SS_iter(int q_id, model& our_model);
double basic_function(double*q, double*p, model& our_model); //Used it for Chi2, JS and Hellinger kernel functions

//Only for addictive kernel functions
void SS_one_D_additive(double q, int d, model& our_model, double& incr_Value_Pos, double& incr_Value_Neg);
bool check_termination(model& our_model, int q_id, double LB, double UB);
void sort_SS_iter(int q_id, double LB, double UB, vector<double>& LB_d, vector<double>& UB_d, model& our_model);

#endif