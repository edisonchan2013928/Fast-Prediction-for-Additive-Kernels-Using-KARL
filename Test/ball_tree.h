#pragma once
#ifndef BALL_TREE_H
#define BALL_TREE_H

#include "kd_tree.h"

struct order_Entry
{
	int id;
	double dist;
	double norm;
	double ip;
};

class ballTree : public Tree
{
public:
	ballTree(model& our_model);
	void ball_Tree_Recur(kdNode*node, double*center);
	void build_ballTree();
	void updateAugment(kdNode*node);
};

#endif