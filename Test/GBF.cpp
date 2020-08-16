#include "GBF.h"

inline void clearHeap(PQ& pq)
{
	int heapSize = (int)pq.size();
	for (int h = 0; h < heapSize; h++)
		pq.pop();
}

double computeSqNorm(double*q, int dim)
{
	double sqNorm = 0;
	for (int d = 0; d < dim; d++)
		sqNorm += q[d] * q[d];

	return sqNorm;
}

void GBF_iter(int q_id, Tree& t)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L, U;
	double f_cur;
	double val_R;
	//double sq_Euclid;

	model& m = t.our_model;

	Node*rootNode = t.rootNode;

	L = rootNode->LB(q_id, m);
	U = rootNode->UB(q_id, m);

	pq_entry.node = rootNode;
	pq_entry.node_L = L;
	pq_entry.node_U = U;
	pq_entry.discrepancy = U - L;

	pq.push(pq_entry);

	while (pq.size() != 0)
	{
		if (t.our_model.is_tau == 0) //approximate-KAQ
		{
			//validation condition
			if (validate_best(L, U, m.epsilon, val_R) == true)
			{
				if (m.method == 6) //Cache-based method
				{
					m.temp_LB = L;
					m.temp_UB = U;
				}
				else
					m.resultVector[q_id] = val_R;

				clearHeap(pq);
				return;
			}
		}
		else //tau-KAQ
		{
			if (L > m.tau)
			{
				m.resultVector[q_id] = 1;
				clearHeap(pq);
				return;
			}
			if (U < m.tau)
			{
				m.resultVector[q_id] = -1;
				clearHeap(pq);
				return;
			}
		}

		pq_entry = pq.top();
		pq.pop();

		L = L - pq_entry.node_L;
		U = U - pq_entry.node_U;

		curNode = pq_entry.node;

		//leaf Node
		if ((int)curNode->idList.size() <= m.leafCapacity)
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

			pq.push(pq_entry);
		}
	}

	//Case (L=exact=U):
	//m.resultVector.push_back(L);
	if (m.is_tau == 0) //approximate-KAQ
	{
		if (m.method == 6) //Cache-based method
		{
			m.temp_LB = L;
			m.temp_UB = U;
		}
		else
			m.resultVector[q_id] = L;
	}	
	else //tau-KAQ
	{
		if (L > m.tau)
			m.resultVector[q_id] = 1;
		else
			m.resultVector[q_id] = -1;
	}
	clearHeap(pq);
}

void updateResultVector(Tree& t, model& our_model)
{
	for (int r = 0; r < our_model.n_q; r++)
		our_model.resultVector[r] = t.our_model.resultVector[r];
}

void KAQ_Algorithm(model& our_model)
{
	double online_Time;
	Cache c;
	kdTree kd_Tree(our_model);
	ballTree ball_Tree(our_model);
	vector<Tree> multi_trees;

	//Used in cache bound methods
	//vector<Tree> tree_vec;

	#ifdef CACHE_STATISTICS
	c.one_D_refine_counter = 0;
	c.full_refine_counter = 0;
	#endif

	if (our_model.method == 1)//aKDE/ tKDC kd-tree + LB_MBR and UB_MBR
		kd_Tree.rootNode = new kdNode();
	if (our_model.method == 2)//kd-tree + KARL
		kd_Tree.rootNode = new kdNode_KARL();
	if (our_model.method == 3)
		ball_Tree.rootNode = new kdNode();
	if (our_model.method == 4)
		ball_Tree.rootNode = new kdNode_KARL();

	//Preprocessing: build cache/ augment tree(s)
	if (our_model.method == 1 || our_model.method == 2)
		kd_Tree.build_kdTree();
	if (our_model.method == 3 || our_model.method == 4)
		ball_Tree.build_ballTree();
	if (our_model.method == 5 || our_model.method == 6)
		build_uniform_cache(c, our_model);
	if (our_model.method == 6 || our_model.method == 7 || our_model.method == 8)
		build_multiple_trees(multi_trees, our_model);
		
	auto start_s = chrono::high_resolution_clock::now();

	for (int q_id = 0; q_id < our_model.n_q; q_id++)
	{
		//if (q_id == 94)
		//	cout << "HERE!" << endl;
		switch (our_model.method)
		{
			case 0:
				SS_iter(q_id, our_model);
				break;
			case 1:
				GBF_iter(q_id, kd_Tree);
				break;
			case 2:
				kd_Tree.our_model.qSquareNorm = computeSqNorm(our_model.queryMatrix[q_id], our_model.dim);
				GBF_iter(q_id, kd_Tree);
				break;
			case 3:
				GBF_iter(q_id, ball_Tree);
				break;
			case 4:
				ball_Tree.our_model.qSquareNorm = computeSqNorm(our_model.queryMatrix[q_id], our_model.dim);
				GBF_iter(q_id, ball_Tree);
				break;
			case 5:
			case 6:
				cache_bound(q_id, c, our_model, multi_trees);
				break;
			case 7:
			case 8:
				multiple_GBF_iter(q_id, multi_trees, our_model);
				break;
		}
	}

	auto end_s = chrono::high_resolution_clock::now();
	online_Time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;

	cout << "Method " << our_model.method << ": " << ((double)our_model.n_q / online_Time) << " Queries/sec" << endl;

	if (our_model.method == 1 || our_model.method == 2)
		updateResultVector(kd_Tree, our_model);
	if (our_model.method == 3 || our_model.method == 4)
		updateResultVector(ball_Tree, our_model);

	#ifdef CACHE_STATISTICS
	if (our_model.method == 5 || our_model.method == 6)
	{
		cout << "# of one_D_refine: " << c.one_D_refine_counter << endl;
		cout << "# of full_refine: " << c.full_refine_counter << endl;
	}
	#endif
}