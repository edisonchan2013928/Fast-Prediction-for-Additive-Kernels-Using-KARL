#include "GBF.h"

int main(int argc,char**argv)
{
	//string debug_string;
	model our_model;
	init_model(argc, argv, our_model);

	KAQ_Algorithm(our_model);
	outputResultFile(our_model);

	//cin >> debug_string;
}