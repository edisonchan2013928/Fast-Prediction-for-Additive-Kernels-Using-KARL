#include "ball_tree.h"

bool compare_order_Entry(const order_Entry& o_Entry1, const order_Entry& o_Entry2)
{
	return o_Entry1.dist < o_Entry2.dist;
}

bool compare_order_Entry_decr(const order_Entry& o_Entry1, const order_Entry& o_Entry2)
{
	return o_Entry1.dist > o_Entry2.dist;
}

ballTree::ballTree(model& our_model)
{
	this->our_model = our_model;
}

void find_furthest(kdNode*node, double*center, model& our_model, int& id1, int& id2)
{
	double max_dist = -inf;
	double cur_dist;
	int cur_id;

	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		cur_id = node->idList[i];
		cur_dist = euclid_dist(center, our_model.dataMatrix[cur_id], our_model.dim);

		if (cur_dist > max_dist)
		{
			max_dist = cur_dist;
			id1 = cur_id;
		}
	}

	max_dist = -inf;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		cur_id = node->idList[i];
		if (cur_id == id1)
			continue;

		cur_dist = euclid_dist(our_model.dataMatrix[id1], our_model.dataMatrix[cur_id], our_model.dim);

		if (cur_dist > max_dist)
		{
			max_dist = cur_dist;
			id2 = cur_id;
		}
	}
}

void findCenter(kdNode*node, double*center, model& our_model)
{
	int id;

	for (int d = 0; d < our_model.dim; d++)
		center[d] = 0;

	//find positive and negative center
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		for (int d = 0; d < our_model.dim; d++)
			center[d] += our_model.dataMatrix[id][d];
	}

	for (int d = 0; d < our_model.dim; d++)
		center[d] /= (double)node->idList.size();
}

void balance_saparation_basic(int id1, int id2, vector<order_Entry>& o_EntryVector, kdNode*node, model& our_model)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector1;
	vector<order_Entry> o_EntryVector2;
	double euclid_1;
	double euclid_2;

	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		o_Entry.id = node->idList[i];

		if ((int)o_EntryVector1.size() > (int)(node->idList.size() / 2.0))
		{
			o_Entry.dist = euclid_dist(our_model.dataMatrix[id2], our_model.dataMatrix[node->idList[i]], our_model.dim);
			o_EntryVector2.push_back(o_Entry);
			continue;
		}

		if ((int)o_EntryVector2.size() > (int)(node->idList.size() / 2.0))
		{
			o_Entry.dist = euclid_dist(our_model.dataMatrix[id1], our_model.dataMatrix[node->idList[i]], our_model.dim);
			o_EntryVector1.push_back(o_Entry);
			continue;
		}

		euclid_1 = euclid_dist(our_model.dataMatrix[id1], our_model.dataMatrix[node->idList[i]], our_model.dim);
		euclid_2 = euclid_dist(our_model.dataMatrix[id2], our_model.dataMatrix[node->idList[i]], our_model.dim);
		if (euclid_1 < euclid_2)
		{
			o_Entry.dist = euclid_1;
			o_EntryVector1.push_back(o_Entry);
		}
		else
		{
			o_Entry.dist = euclid_2;
			o_EntryVector2.push_back(o_Entry);
		}
	}

	sort(o_EntryVector1.begin(), o_EntryVector1.end(), compare_order_Entry);
	sort(o_EntryVector2.begin(), o_EntryVector2.end(), compare_order_Entry_decr);

	o_EntryVector.insert(o_EntryVector.end(), o_EntryVector1.begin(), o_EntryVector1.end());
	o_EntryVector.insert(o_EntryVector.end(), o_EntryVector2.begin(), o_EntryVector2.end());
	o_EntryVector1.clear();
	o_EntryVector2.clear();
}

void half_division(vector<order_Entry>& o_EntryVector, kdNode*node1, kdNode*node2, model& our_model)
{
	int halfSize = (int)ceil((double)o_EntryVector.size() / 2.0);
	for (int i = 0; i < (int)o_EntryVector.size(); i++)
	{
		if (i < halfSize)
			node1->idList.push_back(o_EntryVector[i].id);
		else
			node2->idList.push_back(o_EntryVector[i].id);
	}

	o_EntryVector.clear();
}

void ballTree::ball_Tree_Recur(kdNode*node, double*center)
{
	int id1, id2;
	kdNode*leftNode;
	kdNode*rightNode;

	if ((int)node->idList.size() <= our_model.leafCapacity)
		return;

	leftNode = node->createNode();
	rightNode = node->createNode();

	vector<order_Entry> o_EntryVector;

	findCenter(node, center, our_model);
	find_furthest(node, center, our_model, id1, id2);
	balance_saparation_basic(id1, id2, o_EntryVector, node, our_model);
	half_division(o_EntryVector, leftNode, rightNode, our_model);

	ball_Tree_Recur(leftNode, center);
	ball_Tree_Recur(rightNode, center);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void ballTree::build_ballTree()
{
	double*center;
	center = new double[our_model.dim];

	for (int i = 0; i < our_model.n_d; i++)
		rootNode->idList.push_back(i);

	ball_Tree_Recur((kdNode*)rootNode, center);
	updateAugment((kdNode*)rootNode);
}

void ballTree::updateAugment(kdNode*node)
{
	node->update_Aug(node, this);
}