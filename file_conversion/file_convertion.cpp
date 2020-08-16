#include "file_convertion.h"

void extract_Feature(fstream& file, SVM_stat&stat, bool is_Query)
{
	double**featureMatrix;
	int n;
	string lineString;
	char*token;
	int d;
	double d_Value;
	double temp_weight;

	if (is_Query == true)
	{
		stat.queryMatrix = new double*[stat.n_q];
		for (int i = 0; i < stat.n_q; i++)
			stat.queryMatrix[i] = new double[stat.dim];

		featureMatrix = stat.queryMatrix;
		n = stat.n_q;
	}
	else
	{
		stat.dataMatrix = new double*[stat.n_d];
		for (int i = 0; i < stat.n_d; i++)
			stat.dataMatrix[i] = new double[stat.dim];

		stat.weightArray = new double[stat.n_d];

		featureMatrix = stat.dataMatrix;
		n = stat.n_d;
	}

	//initialization of matrix to zero
	for (int i = 0; i < n; i++)
		for (int d = 0; d < stat.dim; d++)
			featureMatrix[i][d] = 0;

	for (int i = 0; i < n; i++)
	{
		getline(file, lineString);

		//get the first token from the lineString
		token = strtok((char*)lineString.c_str(), " :");

		if (is_Query == false) //training data
			stat.weightArray[i] = atof(token);
		else
			temp_weight = atof(token);

		//Use tokenization to obtain the feature Array
		while (true)
		{
			token = strtok(NULL, " :");

			//cout<<token<<endl;
			if (token != NULL)
			{
				d = atoi(token);
				if (d == 0)
					break;

				token = strtok(NULL, " :");
				if (token != NULL)
					d_Value = atof(token);
				else
				{
					cout << "Read File ERROR!" << endl;
					exit(1);
				}
				featureMatrix[i][d - 1] = d_Value;
			}
			else
				break;
		}
	}
}

void extract_Statistics(fstream& LibSVM_file, SVM_stat& stat)
{
	string temp;

	LibSVM_file >> temp;
	LibSVM_file >> stat.svm_type;
	LibSVM_file >> temp;
	LibSVM_file >> stat.kernel_type;

	if (stat.kernel_type == "polynomial")
	{
		LibSVM_file >> temp;
		LibSVM_file >> stat.degree;
	}

	if (stat.kernel_type == "rbf" || stat.kernel_type == "polynomial" || stat.kernel_type == "sigmoid")
	{
		LibSVM_file >> temp;
		LibSVM_file >> stat.gamma;
	}

	if (stat.kernel_type == "polynomial" || stat.kernel_type == "sigmoid")
	{
		LibSVM_file >> temp;
		LibSVM_file >> stat.beta;
	}

	LibSVM_file >> temp;
	LibSVM_file >> stat.nr_class;

	LibSVM_file >> temp;
	LibSVM_file >> stat.n_d;

	LibSVM_file >> temp;
	LibSVM_file >> stat.threshold;

	if (stat.svm_type == "c_svc")
	{
		LibSVM_file >> temp; LibSVM_file >> temp; LibSVM_file >> temp;
		LibSVM_file >> temp; LibSVM_file >> temp; LibSVM_file >> temp;
	}

	LibSVM_file >> temp;
	getline(LibSVM_file, temp);

	if (stat.kernel_type == "polynomial")
		stat.kernel_type = "poly";
	if (stat.kernel_type == "sigmoid")
		stat.kernel_type = "sig";
}

void output_to_KARL(fstream& KARL_file, double**featureMatrix, SVM_stat& stat, bool isQuery)
{
	int n;
	if (isQuery == true)
	{
		KARL_file << stat.n_q << " " << stat.dim << endl;
		n = stat.n_q;
	}
	else
	{
		KARL_file << stat.n_d << " " << stat.dim << " " << stat.kernel_type;

		if (stat.kernel_type == "rbf" || stat.kernel_type == "poly" || stat.kernel_type == "sig")
			KARL_file << " " << stat.gamma;
		if (stat.kernel_type == "poly")
			KARL_file << " " << stat.degree << " " << stat.beta;
		if (stat.kernel_type == "sig")
			KARL_file << " " << stat.beta;

		KARL_file << " " << 1 << " " << stat.threshold << endl;
		n = stat.n_d;
	}

	for (int i = 0; i < n; i++)
	{
		if (isQuery == false)
			KARL_file << stat.weightArray[i] << " ";
		for (int d = 0; d < stat.dim; d++)
			KARL_file << featureMatrix[i][d] << " ";
		KARL_file << endl;
	}
}

void LibSVM_to_KARL(char*LibSVM_fileName, char*KARL_fileName, SVM_stat& stat, bool isQuery)
{
	fstream LibSVM_file;
	fstream KARL_file;

	LibSVM_file.open(LibSVM_fileName);
	KARL_file.open(KARL_fileName, ios::in | ios::out | ios::trunc);

	if (LibSVM_file.is_open() == false || KARL_file.is_open() == false)
	{
		cout << "Cannot open file(s)!" << endl;
		exit(0);
	}

	if (isQuery == false)
	{
		extract_Statistics(LibSVM_file, stat);
		extract_Feature(LibSVM_file, stat, isQuery);
		output_to_KARL(KARL_file, stat.dataMatrix, stat, isQuery);
	}
	else
	{
		extract_Feature(LibSVM_file, stat, isQuery);
		output_to_KARL(KARL_file, stat.queryMatrix, stat, isQuery);
	}

	LibSVM_file.close();
	KARL_file.close();
}

double obtain_std(double**featureArray, int n, int d)
{
	double sq_sum = 0;
	double sum = 0;
	double mean;
	double std;

	for (int i = 0; i < n; i++)
	{
		sq_sum += (featureArray[i][d])*(featureArray[i][d]);
		sum += featureArray[i][d];
	}
	mean = sum / n;
	std = sqrt((sq_sum - n * mean*mean) / n);

	return std;
}

void preprocess_Data(double**featureArray, int n, int dim, bool isQuery, double b)
{
	double constant = pow((double)n, 1.0 / ((double)(dim + 4))) / (2.0*b);
	static double*coeff_vec;
	double std;

	if (isQuery == false)
	{
		coeff_vec = new double[dim];

		//obtain the coeff_vec
		for (int d = 0; d < dim; d++)
		{
			std = obtain_std(featureArray, n, d);
			coeff_vec[d] = sqrt(constant / std);
		}
	}

	for (int i = 0; i < n; i++)
		for (int d = 0; d < dim; d++)
			featureArray[i][d] *= coeff_vec[d];
}

void KDE_to_KARL(char*KDE_fileName, char*KARL_fileName, KDE_stat& stat, bool isQuery)
{
	fstream KDE_file;
	fstream KARL_file;
	double**featureMatrix;
	double n;
	int dim;

	KDE_file.open(KDE_fileName);
	KARL_file.open(KARL_fileName, ios::in | ios::out | ios::trunc);

	if (KDE_file.is_open() == false || KARL_file.is_open() == false)
	{
		cout << "Cannot open file(s)!" << endl;
		exit(0);
	}

	if (isQuery == true)
	{
		KDE_file >> stat.n_q;
		KDE_file >> stat.dim;
		stat.queryMatrix = new double*[stat.n_q];
		for (int i = 0; i < stat.n_q; i++)
			stat.queryMatrix[i] = new double[stat.dim];

		n = stat.n_q;
		dim = stat.dim;
		featureMatrix = stat.queryMatrix;
	}
	else
	{
		KDE_file >> stat.n_d;
		KDE_file >> stat.dim;

		stat.dataMatrix = new double*[stat.n_d];
		stat.weightArray = new double[stat.n_d];

		for (int i = 0; i < stat.n_d; i++)
		{
			stat.dataMatrix[i] = new double[stat.dim];
			stat.weightArray[i] = 1;
		}

		n = stat.n_d;
		dim = stat.dim;
		featureMatrix = stat.dataMatrix;
	}

	for (int i = 0; i < n; i++)
		for (int d = 0; d < dim; d++)
			KDE_file >> featureMatrix[i][d];

	preprocess_Data(featureMatrix, n, dim, isQuery, stat.b);

	if (isQuery == false)
	{
		stat.gamma = 1;
		KARL_file << stat.n_d << " " << stat.dim << " rbf " << stat.gamma;
		if (stat.is_tau == false)
			KARL_file << " " << 0 << endl;
		if (stat.is_tau == true)
			KARL_file << " " << 1 << " " << stat.threshold << endl;
	}
	else
		KARL_file << stat.n_q << " " << stat.dim << endl;

	for (int i = 0; i < n; i++)
	{
		if (isQuery == false)
			KARL_file << stat.weightArray[i] << " ";

		for (int d = 0; d < dim; d++)
			KARL_file << featureMatrix[i][d] << " ";
		KARL_file << endl;
	}
	
	KDE_file.close();
	KARL_file.close();
}