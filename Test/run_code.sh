#Before you use our code, please make sure you fully understand our paper.

#Compile our code
g++ -c Euclid_Bound.cpp -w -o Euclid_Bound.o -std=c++11
g++ -c GBF.cpp -w -o GBF.o -std=c++11
g++ -c init_KARL.cpp -w -o init_KARL.o -std=c++11
g++ -c ip_bound.cpp -w -o ip_bound.o -std=c++11
g++ -c kd_tree.cpp -w -o kd_tree.o -std=c++11
g++ -c ball_tree.cpp -w -o ball_tree.o -std=c++11
g++ -c slope_intercept.cpp -w -o slope_intercept.o -std=c++11
g++ -c SS.cpp -w -o SS.o -std=c++11
g++ -c Tree.cpp -w -o Tree.o -std=c++11
g++ -c Validation.cpp -w -o Validation.o -std=c++11
g++ -c Cache.cpp -w -o Cache.o -std=c++11
g++ -c multiple_trees_GBF.cpp -w -o multiple_trees_GBF.o -std=c++11

g++ KARL.cpp -O3 -o KARL Euclid_Bound.o GBF.o init_KARL.o ip_bound.o kd_tree.o ball_tree.o slope_intercept.o SS.o Tree.o Validation.o Cache.o multiple_trees_GBF.o -std=c++11

#These are the input parameters for our code:
#char*querysetFileName = argv[1];
#char*datasetFileName = argv[2];
#our_model.resultFileName = argv[3];
#our_model.method = atoi(argv[4]);
#our_model.leafCapacity = atoi(argv[5]);
#our_model.cache_fileName = argv[6];
#our_model.epsilon = atof(argv[7]);

#The parameters "querysetFileName" and "datasetFileName" denote the testing data and the SVM model respectively, which should follow our data format (see the folder "file_conversion").
#The parameter "our_model.resultFileName" denotes the result file.

#Our code supports different methods, which are specified by the parameter "our_model.method", where: 
#For Gaussian, polynomial and sigmoid kernels
#method = 1 (SOTA, Table 9)
#method = 2 (KARL, Table 9)

#For additive kernels
#method = 8 (SOTA, Table 10)
#method = 7 (KARL_{linear}, Table 10)
#method = 5 (KARL_{mono}, Table 10)
#method = 6 (KARL, Table 10)

#The parameter "our_model.leafCapacity" denotes the leaf capacity of the index structure, e.g., kd-tree.
#The parameter "our_model.cache_fileName" denotes the cache file, whcih stores the precomputed values, for our method KARL_{mono} (cf. Section 5.1).
#The parameter "our_model.epsilon" denotes the epsilon, i.e., tolerance, in Table 7.

#In this example, we test the support vector classification (SVC) in cod-rna dataset, using Chi2 kernel.
./KARL ../datasets/cod-rna_our_test ../datasets/cod-rna_Chi2_our_model ./Results/cod-rna_Chi2_M8 8 20 -1 -1
./KARL ../datasets/cod-rna_our_test ../datasets/cod-rna_Chi2_our_model ./Results/cod-rna_Chi2_M7 7 20 -1 -1
./KARL ../datasets/cod-rna_our_test ../datasets/cod-rna_Chi2_our_model ./Results/cod-rna_Chi2_M5 5 0 ./Cache/cod-rna_Chi2_Cache_e00125 0.0125
./KARL ../datasets/cod-rna_our_test ../datasets/cod-rna_Chi2_our_model ./Results/cod-rna_Chi2_M6 6 20 ./Cache/cod-rna_Chi2_Cache_e00125 0.0125