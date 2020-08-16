#pragma once
#ifndef TREE_H
#define TREE_H

#include "init_KARL.h"
#include "Euclid_Bound.h"
#include "ip_bound.h"

class Node;
class Tree;

class Node
{
public:
	double sum_w_twin[2];
	vector<int> idList;
	double**boundary_twin[2];
	double**boundary;

	vector<Node*> childVector;

	virtual double LB(int q_id, model& our_model) = 0;
	virtual double UB(int q_id, model& our_model) = 0;
	//virtual double LB(double*q, int dim, model& our_model) = 0;
	//virtual double UB(double*q, int dim, model& our_model) = 0;

	virtual Node*createNode() = 0;
	virtual void update_Aug(Node*node, Tree*t) = 0;
	void update_Node(Node*node, Tree*t);
};

class Tree
{
public:
	model our_model;
	//int leafCapacity; //Set to be 20 in tKDC
	Node*rootNode;
};

#endif