#pragma once
#ifndef CACHE_H
#define CACHE_H

#include "init_KARL.h"
#include "SS.h"
#include "Validation.h"
#include "GBF.h"

#define CACHE_STATISTICS

const double min_pos_value = 0.01;

struct Cache
{
	//Cache with uniform condition
	//*********************************************//
	//vector< vector<double> > q_cache_vec;
	vector<double> q_min_cache;
	vector< vector<double> > uniform_cache_pos;
	vector< vector<double> > uniform_cache_neg;
	//double uni_length; //Cache with uniform-length in all dimension
	//int uni_size; //Cache with uniform-size in all dimension
	vector<double> uni_length_vec; //Cache with uniform-length in each dimension
	//*********************************************//

	#ifdef CACHE_STATISTICS
	double one_D_refine_counter;
	double full_refine_counter;
	#endif
};

void obtain_boundary(model& our_model, int d, double& x_min, double& x_max);
double estimate_uni_length(model& our_model, double p_ell_min, double p_ell_max);
void build_uniform_cache(Cache& c, model& our_model);
void uni_length_Cache(Cache& c, model& our_model);
//void uni_size_Cache(Cache& c, model& our_model);
void cache_bound(int q_id, Cache& c, model& our_model, vector<Tree>& multi_trees);
void store_Cache(model& our_model, Cache& c);
bool load_Cache(model& our_model, Cache& c);

#endif