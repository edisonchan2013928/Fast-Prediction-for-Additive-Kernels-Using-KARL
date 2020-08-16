#pragma once
#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"
#include "slope_intercept.h"

class kdNode : public Node
{
public:
	double LB(int q_id, model& our_model);
	double UB(int q_id, model& our_model);

	//void update_Aug(Node*node, Tree*t) {}
	void update_Aug(Node*node, Tree*t);

	kdNode*createNode();
};

class kdNode_KARL : public kdNode
{
public:
	double*a_G_twin[2];
	double b_G_twin[2];

	//facilitates online (sharing) computation
	double gamma_sum_twin[2];

	double LB(int q_id, model& our_model);
	double UB(int q_id, model& our_model);

	void update_KARL_Info(kdNode_KARL*node, Tree*t);

	void update_Aug(Node*node, Tree*t);
	kdNode_KARL*createNode();
};

class kdTree : public Tree
{
public:
	kdTree(model& our_model);
	double obtain_SplitValue(kdNode*node, int split_Dim);
	void KD_Tree_Recur(kdNode*node, int split_Dim);
	void build_kdTree();

	void updateAugment(kdNode*node);
};

#endif