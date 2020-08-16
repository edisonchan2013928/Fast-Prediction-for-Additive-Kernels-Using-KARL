g++ -c file_convertion.cpp -w -o file_convertion.o -std=c++11
g++ LibSVM_binary_to_KARL.cpp -O3 -o LibSVM_binary_to_KARL file_convertion.o -std=c++11
#char*LibSVM_fileName = argv[1];
#char*KARL_fileName = argv[2];
#bool isQuery = (bool)atoi(argv[3]);
#
#stat.dim = atoi(argv[4]);
#if (isQuery == true)
#	stat.n_q = atoi(argv[5]);

#Once you have obtained the SVM model from the LIBSVM library, you need to use this code to convert both the testing dataset and the model file (in SVM format) to our format.
#The parameter stat.n_q denotes the number of queries (or testing data)
#Here, we show you the example for cod-rna testing file and model file. 

#convert testing file
./LibSVM_binary_to_KARL ../datasets/cod-rna_t_scale ../datasets/cod-rna_our_test 1 8 10000
#convert training file
./LibSVM_binary_to_KARL ../datasets/cod-rna_scale_Chi2_train ../datasets/cod-rna_Chi2_our_model 0 8