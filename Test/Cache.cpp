#include "Cache.h"

void obtain_boundary(model& our_model, int d, double& x_min, double& x_max)
{
	x_max = -inf;
	x_min = inf;
	for (int i = 0; i < our_model.n_d; i++)
	{
		if (our_model.dataMatrix[i][d] > x_max)
			x_max = our_model.dataMatrix[i][d];
		if (our_model.dataMatrix[i][d] < x_min)
			x_min = our_model.dataMatrix[i][d];
	}

	if (x_min < eps)
		x_min = min_pos_value;
}

double estimate_uni_length(model& our_model, double p_ell_min, double p_ell_max)
{
	if (our_model.kernel_type == "Chi2")
		return our_model.epsilon*min_pos_value;
	if (our_model.kernel_type == "JS")
		return our_model.epsilon*min(min_pos_value, p_ell_min*log2(min_pos_value / p_ell_max + 1));
	if (our_model.kernel_type == "Hellinger")
		return our_model.epsilon*min(min_pos_value, p_ell_min);

	cout << "Error!" << endl;
	exit(0);
}

void uni_length_Cache(Cache& c, model& our_model)
{
	double x_min, x_max;
	int cache_size;
	double q;
	vector<double> d_vec;
	double incr_Value_Pos, incr_Value_Neg;
	double space_consumption = 0;

	if (load_Cache(our_model, c) == true)
		return;

	for (int d = 0; d < our_model.dim; d++)
	{
		c.uniform_cache_pos.push_back(d_vec);
		c.uniform_cache_neg.push_back(d_vec);
		//c.q_cache_vec.push_back(d_vec);

		obtain_boundary(our_model, d, x_min, x_max);
		c.uni_length_vec.push_back(estimate_uni_length(our_model, x_min, x_max));
		cache_size = (int)ceil((x_max - min_pos_value) / c.uni_length_vec[d]);

		space_consumption += cache_size * 8 * 2; //We have two vectors in the Cache
		if (space_consumption > 2.0 * 1024 * 1024 * 1024)
		{
			cout << "The space consumption is larger than 2GB!" << endl;
			exit(0);
		}

		c.q_min_cache.push_back(min_pos_value);
		for (int l = 0; l < cache_size; l++)
		{
			q = min_pos_value + l * c.uni_length_vec[d];
			SS_one_D_additive(q, d, our_model, incr_Value_Pos, incr_Value_Neg);

			//c.q_cache_vec[d].push_back(q);
			c.uniform_cache_pos[d].push_back(incr_Value_Pos);
			c.uniform_cache_neg[d].push_back(incr_Value_Neg);
		}
	}

	store_Cache(our_model, c);
}

/*void uni_size_Cache(Cache& c, model& our_model)
{
	//code here
}*/

void build_uniform_cache(Cache& c, model& our_model)
{
	if (our_model.method == 5 || our_model.method == 6) //Cache with uniform-length in each dimension
		uni_length_Cache(c, our_model);
	
	/*if (our_model.method == 6) //Cache with uniform-size in each dimension
	{
		c.uni_size = our_model.uni_size;
		uni_size_Cache(c, our_model);
	}*/
}

void cache_bound(int q_id, Cache& c, model& our_model, vector<Tree>& multi_trees)
{
	double*q = our_model.queryMatrix[q_id];
	double q_min;
	int q_index_low, q_index_up;
	double incr_Value_Pos, incr_Value_Neg;
	double temp;
	double LB = 0;
	double UB = 0;
	//Used in method 6
	vector<double> LB_d;
	vector<double> UB_d;

	for (int d = 0; d < our_model.dim; d++)
	{
		q_min = c.q_min_cache[d];
		//q_min = c.q_cache_vec[d][0];
		temp = (q[d] - q_min) / c.uni_length_vec[d];
		q_index_low = (int)floor(temp);
		q_index_up = (int)ceil(temp);

		if (q_index_low < 0 || q_index_up >= (int)c.uniform_cache_pos[d].size())
		{
			#ifdef CACHE_STATISTICS
			c.one_D_refine_counter++;
			#endif

			if (our_model.method == 5)
			{
				SS_one_D_additive(q[d], d, our_model, incr_Value_Pos, incr_Value_Neg);
				LB_d.push_back(incr_Value_Pos - incr_Value_Neg);
				UB_d.push_back(incr_Value_Pos - incr_Value_Neg);
			}
			if (our_model.method == 6)
			{
				GBF_iter(q_id, multi_trees[d]);
				LB_d.push_back(multi_trees[d].our_model.temp_LB);
				UB_d.push_back(multi_trees[d].our_model.temp_UB);
				//SS_one_D_additive(q[d], d, our_model, incr_Value_Pos, incr_Value_Neg);
			}

			LB += LB_d[d];
			UB += UB_d[d];
			continue;
		}

		LB_d.push_back(c.uniform_cache_pos[d][q_index_low] - c.uniform_cache_neg[d][q_index_up]);
		UB_d.push_back(c.uniform_cache_pos[d][q_index_up] - c.uniform_cache_neg[d][q_index_low]);
		LB += LB_d[d];
		UB += UB_d[d];
	}

	if (check_termination(our_model, q_id, LB, UB) == false)
	{
		#ifdef CACHE_STATISTICS
		c.full_refine_counter++;
		#endif

		if (our_model.method == 5 || our_model.method == 6)
			sort_SS_iter(q_id, LB, UB, LB_d, UB_d, our_model);
	}
}

void store_Cache(model& our_model, Cache& c)
{
	fstream cache_file;
	int size;

	cache_file.open(our_model.cache_fileName, ios::in | ios::out | ios::trunc);

	if (cache_file.is_open() == false)
	{
		cout << "Cannot open this cache file!" << endl;
		exit(0);
	}

	cache_file << our_model.dim << endl;
	for (int d = 0; d < our_model.dim; d++)
	{
		size = (int)c.uniform_cache_pos[d].size();
		cache_file << size << endl;
		if (size == 0)
			continue;
		cache_file << c.q_min_cache[d] << " " << c.uni_length_vec[d] << endl;
		for (int s = 0; s < size; s++)
			cache_file << c.uniform_cache_pos[d][s] << " " << c.uniform_cache_neg[d][s] << " ";
		cache_file << endl;
	}

	cache_file.close();
}

bool load_Cache(model& our_model, Cache& c)
{
	fstream cache_file;
	int size;
	int dim;
	double value;
	//double q_ell_min;
	double uni_length;
	vector<double> d_vec;

	cache_file.open(our_model.cache_fileName);

	if (cache_file.is_open() == false)
		return false;

	cache_file >> dim;
	if (dim != our_model.dim)
	{
		cout << "Incorrect dimension!" << endl;
		exit(0);
	}
	for (int d = 0; d < dim; d++)
	{
		//c.q_cache_vec.push_back(d_vec);
		c.uniform_cache_pos.push_back(d_vec);
		c.uniform_cache_neg.push_back(d_vec);

		cache_file >> size;
		if (size == 0)
			continue;

		cache_file >> value;
		c.q_min_cache.push_back(value);
		cache_file >> uni_length;
		c.uni_length_vec.push_back(uni_length);

		for (int i = 0; i < size; i++)
		{
			cache_file >> value;
			c.uniform_cache_pos[d].push_back(value);
			cache_file >> value;
			c.uniform_cache_neg[d].push_back(value);
			//c.q_cache_vec[d].push_back(q_ell_min + i * c.uni_length_vec[d]);
		}
	}

	return true;
}