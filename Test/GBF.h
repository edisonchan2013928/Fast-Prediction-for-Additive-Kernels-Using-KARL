#pragma once
#ifndef GBF_H
#define GBF_H

#include "init_KARL.h"
#include "kd_tree.h"
#include "ball_tree.h"
#include "Validation.h"
#include "SS.h"
#include "Cache.h"
#include "multiple_trees_GBF.h"

struct pqNode
{
	Node*node;
	double discrepancy;
	double node_L;
	double node_U;

	//int d; //used for additive kernels
};

//This is the maximum heap
struct comparePriority
{
	bool operator()(pqNode& p1, pqNode& p2)
	{
		return p1.discrepancy < p2.discrepancy;
	}
};

typedef priority_queue<pqNode, vector<pqNode>, comparePriority> PQ;

double computeSqNorm(double*q, int dim);
void GBF_iter(int q_id, Tree& t);
void KAQ_Algorithm(model& our_model);

#endif