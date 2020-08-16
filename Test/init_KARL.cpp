#include "init_KARL.h"

void init_Matrix(fstream& file, double**& Matrix, model& our_model, int n, bool is_Query)
{
	int dim = our_model.dim;

	Matrix = new double*[n];
	for (int i = 0; i < n; i++)
		Matrix[i] = new double[dim];

	if (is_Query == false)
	{
		our_model.weightVector = new double[n];
		our_model.weight_oriVector = new double[n];
	}

	for (int i = 0; i < n; i++)
	{
		if (is_Query == false)
		{
			file >> our_model.weightVector[i];
			our_model.weight_oriVector[i] = our_model.weightVector[i];
		}

		for (int d = 0; d < dim; d++)
			file >> Matrix[i][d];
	}
}

void find_num_of_pos_samples(model& our_model)
{
	int n_d = our_model.n_d;
	our_model.n_d_pos = 0;

	for (int i = 0; i < n_d; i++)
	{
		if (our_model.weightVector[i] >= 0)
			our_model.n_d_pos++;
	}
}

void turn_weight_to_pos(model& our_model)
{
	int n_d = our_model.n_d;
	int n_d_pos = our_model.n_d_pos;

	//turn the weightVector to positive
	for (int i = n_d_pos; i < n_d; i++)
		our_model.weightVector[i] = fabs(our_model.weightVector[i]);
}

void loadData(char*fileName, model& our_model, bool is_Query)
{
	fstream file;
	file.open(fileName);
	//double**Matrix;

	if (file.is_open() == false)
	{
		if (is_Query == true)
			cout << "Cannot open query file!" << endl;
		if (is_Query == false)
			cout << "Cannot open model file!" << endl;
		exit(0);
	}

	if (is_Query == true)
	{
		file >> our_model.n_q;
		file >> our_model.dim;

		init_Matrix(file, our_model.queryMatrix, our_model, our_model.n_q, is_Query);
	}
	else
	{
		file >> our_model.n_d;
		file >> our_model.dim;
		file >> our_model.kernel_type;
		if (our_model.kernel_type == "rbf" || our_model.kernel_type == "poly" || our_model.kernel_type == "sig")
			file >> our_model.gamma;

		if (our_model.kernel_type == "poly")
		{
			file >> our_model.deg;
			file >> our_model.beta;

			if ((int)(our_model.deg) % 2 != 0)
				our_model.is_odd_degree = true;
			else
				our_model.is_odd_degree = false;
		}
		if (our_model.kernel_type == "sig")
			file >> our_model.beta;

		file >> our_model.is_tau;

		if (our_model.is_tau == 1)
			file >> our_model.tau;

		init_Matrix(file, our_model.dataMatrix, our_model, our_model.n_d, is_Query);
		find_num_of_pos_samples(our_model);
		turn_weight_to_pos(our_model);
	}

	file.close();
}

void init_model(int argc, char**argv, model& our_model)
{
	if (argc < 5)
	{
		cout << "The number of arguments is too few!" << endl;
		exit(0);
	}

	char*querysetFileName = argv[1];
	char*datasetFileName = argv[2];

	loadData(querysetFileName, our_model, true);
	loadData(datasetFileName, our_model, false);

	our_model.resultFileName = argv[3];
	our_model.method = atoi(argv[4]);
	our_model.leafCapacity = atoi(argv[5]);
	our_model.cache_fileName = argv[6];

	//if (our_model.is_tau == 0)
	our_model.epsilon = atof(argv[7]);

	for (int q = 0; q < our_model.n_q; q++)
		our_model.resultVector.push_back(-100);

	our_model.internalCapacity = 20; //only used for m-tree
}

void outputResultFile(model& our_model)
{
	fstream resultFile;
	resultFile.open(our_model.resultFileName, ios::in | ios::out | ios::trunc);

	if (resultFile.is_open() == false)
	{
		cout << "Cannot Open Result File!" << endl;
		return;
	}

	for (int r = 0; r < (int)our_model.resultVector.size(); r++)
		resultFile << our_model.resultVector[r] << endl;

	resultFile.close();
}
