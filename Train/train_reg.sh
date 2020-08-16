#Train different datasets with different kernel functions

#Train with Chi-squared kernel function
###################################################
#3D_spatial_network
./svm-train -t 5 -s 3 -m 1024 ../../datasets_reg/3D_spatial_network/3D_spatial_network_scale ../../datasets_reg/3D_spatial_network/3D_spatial_network_Chi2_train

#cadata
./svm-train -t 5 -s 3 -m 1024 ../../datasets_reg/cadata/cadata_scale ../../datasets_reg/cadata/cadata_Chi2_train

#Wave
./svm-train -t 5 -s 3 -m 1024 ../../datasets_reg/Wave/Wave_scale ../../datasets_reg/Wave/Wave_Chi2_train
###################################################

#Train with Hellinger kernel function
###################################################
#3D_spatial_network
./svm-train -t 7 -s 3 -m 1024 ../../datasets_reg/3D_spatial_network/3D_spatial_network_scale ../../datasets_reg/3D_spatial_network/3D_spatial_network_Hellinger_train

#cadata
./svm-train -t 7 -s 3 -m 1024 ../../datasets_reg/cadata/cadata_scale ../../datasets_reg/cadata/cadata_Hellinger_train

#Wave
./svm-train -t 7 -s 3 -m 1024 ../../datasets_reg/Wave/Wave_scale ../../datasets_reg/Wave/Wave_Hellinger_train
###################################################

#Train with JS kernel function
###################################################
#3D_spatial_network
./svm-train -t 6 -s 3 -m 1024 ../../datasets_reg/3D_spatial_network/3D_spatial_network_scale ../../datasets_reg/3D_spatial_network/3D_spatial_network_JS_train

#cadata
./svm-train -t 6 -s 3 -m 1024 ../../datasets_reg/cadata/cadata_scale ../../datasets_reg/cadata/cadata_JS_train

#Wave
./svm-train -t 6 -s 3 -m 1024 ../../datasets_reg/Wave/Wave_scale ../../datasets_reg/Wave/Wave_JS_train
###################################################

