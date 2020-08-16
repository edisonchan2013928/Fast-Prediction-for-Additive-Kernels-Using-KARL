#include "Tree.h"

void update_boundary(double**boundary, double*f_vector,int dim)
{
	for (int d = 0; d < dim; d++)
	{
		if (boundary[d][0] > f_vector[d])
			boundary[d][0] = f_vector[d];
		if (boundary[d][1] < f_vector[d])
			boundary[d][1] = f_vector[d];
	}
}

void Node::update_Node(Node*node, Tree*t)
{
	int n_d_pos = t->our_model.n_d_pos;
	int dim = t->our_model.dim;
	int id;
	int t_sel;

	//positive and negative boundaries initialization
	for (int pn = 0; pn <= 1; pn++) 
	{
		node->boundary_twin[pn] = new double*[dim];
		for (int d = 0; d < dim; d++)
			node->boundary_twin[pn][d] = new double[2];
	}

	for (int pn = 0; pn <= 1; pn++)
	{
		for (int d = 0; d < dim; d++)
		{
			node->boundary_twin[pn][d][0] = inf;
			node->boundary_twin[pn][d][1] = -inf;
		}
	}
	
	node->sum_w_twin[0] = 0;
	node->sum_w_twin[1] = 0;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		if (id < n_d_pos)
			t_sel = 0;
		else
			t_sel = 1;

		node->sum_w_twin[t_sel] += t->our_model.weightVector[id];
		update_boundary(node->boundary_twin[t_sel], t->our_model.dataMatrix[id], dim);
	}

	if ((int)node->idList.size() <= t->our_model.leafCapacity) //this is the leaf node
		return;

	for (int c = 0; c < (int)node->childVector.size(); c++)
		update_Node(node->childVector[c],t);
}