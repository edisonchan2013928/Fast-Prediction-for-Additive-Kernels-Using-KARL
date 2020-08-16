#Train different datasets with different kernel functions
#./svm-train -t 5 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_Chi2_train
#./svm-train -t 6 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_JS_train
#./svm-train -t 7 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_Hellinger_train
#./svm-train ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_train

./svm-train -t 5 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_Chi2_train
./svm-train -t 6 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_JS_train
./svm-train -t 7 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_Hellinger_train
./svm-train ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_train

./svm-train -t 5 ../../datasets/home/home_scale ../../datasets/home/home_scale_Chi2_train
./svm-train -t 6 ../../datasets/home/home_scale ../../datasets/home/home_scale_JS_train
./svm-train -t 7 ../../datasets/home/home_scale ../../datasets/home/home_scale_Hellinger_train
./svm-train ../../datasets/home/home_scale ../../datasets/home/home_scale_train
exit

#Train with Chi-squared kernel function
###################################################
#a9a
./svm-train -t 5 ../../datasets/a9a/a9a_scale ../../datasets/a9a/a9a_scale_Chi2_train

#ijcnn1
./svm-train -t 5 ../../datasets/ijcnn1/ijcnn1_scale ../../datasets/ijcnn1/ijcnn1_scale_Chi2_train

#Skin-non-skin
./svm-train -t 5 ../../datasets/skin_nonskin/skin_nonskin_scale ../../datasets/skin_nonskin/skin_nonskin_scale_Chi2_train

#covtype-b
./svm-train -t 5 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_Chi2_train

#cod-rna
./svm-train -t 5 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_Chi2_train

#home
./svm-train -t 5 ../../datasets/home/home_scale ../../datasets/home/home_scale_Chi2_train
###################################################

#Train with JS-squared kernel function
###################################################
#a9a
./svm-train -t 6 ../../datasets/a9a/a9a_scale ../../datasets/a9a/a9a_scale_JS_train

#ijcnn1
./svm-train -t 6 ../../datasets/ijcnn1/ijcnn1_scale  ../../datasets/ijcnn1/ijcnn1_scale_JS_train

#Skin-non-skin
./svm-train -t 6 ../../datasets/skin_nonskin/skin_nonskin_scale ../../datasets/skin_nonskin/skin_nonskin_scale_JS_train

#covtype-b
./svm-train -t 6 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_JS_train

#cod-rna
./svm-train -t 6 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_JS_train

#home
./svm-train -t 6 ../../datasets/home/home_scale ../../datasets/home/home_scale_JS_train
###################################################

#Train with Hellinger kernel function
###################################################
#a9a
./svm-train -t 7 ../../datasets/a9a/a9a_scale ../../datasets/a9a/a9a_scale_Hellinger_train

#ijcnn1
./svm-train -t 7 ../../datasets/ijcnn1/ijcnn1_scale  ../../datasets/ijcnn1/ijcnn1_scale_Hellinger_train

#Skin-non-skin
./svm-train -t 7 ../../datasets/skin_nonskin/skin_nonskin_scale ../../datasets/skin_nonskin/skin_nonskin_scale_Hellinger_train

#covtype-b
./svm-train -t 7 ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_Hellinger_train

#cod-rna
./svm-train -t 7 ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_Hellinger_train

#home
./svm-train -t 7 ../../datasets/home/home_scale ../../datasets/home/home_scale_Hellinger_train
###################################################

exit

#Train with Gaussian kernel function
###################################################
#a9a
./svm-train ../../datasets/a9a/a9a_scale ../../datasets/a9a/a9a_scale_train

#ijcnn1
./svm-train ../../datasets/ijcnn1/ijcnn1_scale ../../datasets/ijcnn1/ijcnn1_scale_train

#Skin-non-skin
./svm-train ../../datasets/skin_nonskin/skin_nonskin_scale ../../datasets/skin_nonskin/skin_nonskin_scale_train

#covtype-b
./svm-train ../../datasets/covtype/covtype.libsvm.binary.scale ../../datasets/covtype/covtype_scale_train

#cod-rna
./svm-train ../../datasets/cod-rna/cod-rna_scale ../../datasets/cod-rna/cod-rna_scale_train

#home
./svm-train ../../datasets/home/home_scale ../../datasets/home/home_scale_train
###################################################
