#pragma once
#ifndef FILE_CONVERTION_H
#define FILE_CONVERTION_H

#define _CRT_SECURE_NO_WARNINGS

#include "Library.h"

struct statistics
{
	int dim;
	double**queryMatrix;
	double**dataMatrix;
	double*weightArray;
	int n_d;
	int n_q;
	double threshold;

	double gamma; //Gaussian (rbf) kernel
};

struct KDE_stat : public statistics 
{
	double b; //bandwidth selection
	bool is_tau; //KDC or KDE
};

struct SVM_stat : public statistics
{
	string svm_type;
	string kernel_type;
	int nr_class;
	int*label;
	int*nr_sv;

	//Used in Polynomial or Sigmoid kernel
	double degree;
	double beta;
};

void LibSVM_to_KARL(char*LibSVM_fileName, char*KARL_fileName, SVM_stat& stat, bool isQuery);

//We currently only assume using rbf kernel function
void KDE_to_KARL(char*KDE_fileName, char*KARL_fileName, KDE_stat& stat, bool isQuery); 

#endif