#include "file_convertion.h"

int main(int argc, char**argv)
{
	SVM_stat stat;
	char*LibSVM_fileName = argv[1];
	char*KARL_fileName = argv[2];
	bool isQuery = (bool)atoi(argv[3]);
	
	stat.dim = atoi(argv[4]);
	if (isQuery == true)
		stat.n_q = atoi(argv[5]);

	//Example
	/*char*LibSVM_fileName = (char*)"../../../Datasets/ijcnn1/ijcnn1_test";
	char*KARL_fileName = (char*)"ijcnn1_our_test";
	bool isQuery = true;

	stat.dim = 22;
	if (isQuery == true)
		stat.n_q = 91701;*/

	LibSVM_to_KARL(LibSVM_fileName, KARL_fileName, stat, isQuery);
}