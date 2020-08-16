#pragma once
#ifndef INIT_KARL_H
#define INIT_KARL_H

#include "Library.h"

const double inf = 99999999999;
const double eps = 0.000000001;

struct model
{
	int n_d;
	int n_d_pos; //The number of positive samples
	int n_q;
	int dim;
	int method;
	double**queryMatrix;
	double**dataMatrix;
	double*weightVector;
	double*weight_oriVector; //original weightvector

	int is_tau;

	//kernel functions
	//Gaussian kernel: exp(-gamma ||q-p||^2)
	//Polynomial kernel: (gamma q^T p + beta)^deg
	//Sigmoid kernel: tanh(gamma q^T p + beta)
	string kernel_type;
	double gamma;
	double beta;
	double deg;
	bool is_odd_degree;

	//regression/density estimation models
	double epsilon;

	//classification model
	double tau;

	//output parameters
	char*resultFileName;
	vector<double> resultVector;

	//Tree specification
	int leafCapacity; 
	int internalCapacity;
	//char*bulkLoad_TreeName; //Used for bulkLoad m-tree

	//faciliates online computation
	double qSquareNorm;

	//Information for Cache
	//double uni_length; //Cache with uniform-length in each dimension
	//int uni_size; //Cache with uniform-size in each dimension
	char*cache_fileName; //Store or load the cache file!
	double temp_LB;
	double temp_UB;
};

//void loadData(char*fileName, model& our_model, bool is_Query);
void init_model(int argc, char**argv, model& our_model);
void outputResultFile(model& our_model);

//debug (add)
void loadData(char*fileName, model& our_model, bool is_Query);

#endif