#include "multiple_trees_GBF.h"

inline void clearHeap_Multi_trees(vector<PQ>& pq_vec, int dim)
{
	int heapSize;

	for (int d = 0; d < dim; d++)
	{
		heapSize = (int)pq_vec[d].size();
		for (int h = 0; h < heapSize; h++)
			pq_vec[d].pop();
	}
}

void update_models(vector<model>& our_model_vec, model& our_model)
{
	model temp;
	temp.dim = 1;
	temp.epsilon = our_model.epsilon;
	//temp.queryMatrix = our_model.queryMatrix;
	if (our_model.method == 6)
		temp.is_tau = false;
	temp.kernel_type = our_model.kernel_type;
	temp.leafCapacity = our_model.leafCapacity;
	temp.method = our_model.method;
	temp.n_d = our_model.n_d;
	temp.n_d_pos = our_model.n_d_pos;
	temp.n_q = our_model.n_q;
	temp.tau = our_model.tau;
	temp.weightVector = our_model.weightVector;
	temp.weight_oriVector = our_model.weight_oriVector;

	for (int d = 0; d < our_model.dim; d++)
	{
		temp.queryMatrix = new double*[our_model.n_q];
		temp.dataMatrix = new double*[our_model.n_d];
		for (int q = 0; q < our_model.n_q; q++)
		{
			temp.queryMatrix[q] = new double[1];
			temp.queryMatrix[q][0] = our_model.queryMatrix[q][d];
		}
		for (int i = 0; i < our_model.n_d; i++)
		{
			temp.dataMatrix[i] = new double[1];
			temp.dataMatrix[i][0] = our_model.dataMatrix[i][d];
		}
		our_model_vec.push_back(temp);
	}
}

void build_multiple_trees(vector<Tree>& multi_trees, model& our_model)
{
	vector<model> our_model_vec;

	update_models(our_model_vec, our_model);
	for (int d = 0; d < our_model.dim; d++)
	{
		//method = 6: Cache_bound + linear bound for q_ell which is out-of-range
		//method = 7: linear bound
		//method = 8: SOTA bound (min/max bound)
		if (our_model.method == 6 || our_model.method == 7)
		{
			ballTree ball_Tree(our_model_vec[d]);
			multi_trees.push_back(ball_Tree);
			multi_trees[d].rootNode = new kdNode_KARL();
			//multi_trees[d].our_model = our_model_vec[d];
			((ballTree&)multi_trees[d]).build_ballTree();
		}
		if (our_model.method == 8) //SOTA
		{
			ballTree ball_Tree(our_model_vec[d]);
			multi_trees.push_back(ball_Tree);
			multi_trees[d].rootNode = new kdNode();
			((ballTree&)multi_trees[d]).build_ballTree();
		}
	}
}

int which_to_pop(vector<PQ>& pq_vec, int dim)
{
	double discrepancy;
	double max_discrepancy = -inf;
	int dimension = -100;
	for (int d = 0; d < dim; d++)
	{
		if ((int)pq_vec[d].size() == 0)
			continue;

		discrepancy = pq_vec[d].top().discrepancy;
		
		if (discrepancy > max_discrepancy)
		{
			max_discrepancy = discrepancy;
			dimension = d;
		}
	}

	return dimension;
}

void multiple_GBF_iter(int q_id, vector<Tree>& multi_trees, model& our_model)
{
	static vector<PQ> pq_vec;
	static int is_inside = 0;
	pqNode pq_entry;
	Node*curNode;
	double L, U;
	double temp_L, temp_U;
	double f_cur;
	double val_R;
	int dimension;
	int queue_vec_size = 0;
	L = 0; 
	U = 0;

	if (is_inside == 0)
	{
		PQ pq;
		for (int d = 0; d < our_model.dim; d++)
			pq_vec.push_back(pq);
		is_inside = 1;
	}

	for (int d = 0; d < our_model.dim; d++)
	{
		model& m = multi_trees[d].our_model;
		temp_L = multi_trees[d].rootNode->LB(q_id, m);
		temp_U = multi_trees[d].rootNode->UB(q_id, m);

		pq_entry.node = multi_trees[d].rootNode;
		pq_entry.node_L = temp_L;
		pq_entry.node_U = temp_U;
		pq_entry.discrepancy = temp_U - temp_L;

		L = L + temp_L;
		U = U + temp_U;

		pq_vec[d].push(pq_entry);
		queue_vec_size++;
	}

	while (queue_vec_size != 0)
	{
		if (our_model.is_tau == 0) //approximate-KAQ
		{
			//validation condition
			if (validate_best(L, U, our_model.epsilon, val_R) == true)
			{
				our_model.resultVector[q_id] = val_R;
				clearHeap_Multi_trees(pq_vec, our_model.dim);
				return;
			}
		}
		else //tau-KAQ
		{
			if (L > our_model.tau)
			{
				our_model.resultVector[q_id] = 1;
				clearHeap_Multi_trees(pq_vec, our_model.dim);
				return;
			}
			if (U < our_model.tau)
			{
				our_model.resultVector[q_id] = -1;
				clearHeap_Multi_trees(pq_vec, our_model.dim);
				return;
			}
		}

		dimension = which_to_pop(pq_vec, our_model.dim);
		pq_entry = pq_vec[dimension].top();
		pq_vec[dimension].pop();
		queue_vec_size--;

		L = L - pq_entry.node_L;
		U = U - pq_entry.node_U;

		curNode = pq_entry.node;
		model& m = multi_trees[dimension].our_model;

		//leaf Node
		if ((int)curNode->idList.size() <= our_model.leafCapacity)
		{
			f_cur = refine(m.queryMatrix[q_id], curNode, m);
			L = L + f_cur;
			U = U + f_cur;

			continue;
		}

		//Non-Leaf Node
		for (int c = 0; c < (int)curNode->childVector.size(); c++)
		{
			pq_entry.node_L = curNode->childVector[c]->LB(q_id, m);
			pq_entry.node_U = curNode->childVector[c]->UB(q_id, m);
			pq_entry.discrepancy = pq_entry.node_U - pq_entry.node_L;
			pq_entry.node = curNode->childVector[c];

			L = L + pq_entry.node_L;
			U = U + pq_entry.node_U;
			pq_vec[dimension].push(pq_entry);
			queue_vec_size++;
		}
	}

	if (our_model.is_tau == 0) //approximate-KAQ
		our_model.resultVector[q_id] = L;
	else //tau-KAQ
	{
		if (L > our_model.tau)
			our_model.resultVector[q_id] = 1;
		else
			our_model.resultVector[q_id] = -1;
	}
	clearHeap_Multi_trees(pq_vec, our_model.dim);
}