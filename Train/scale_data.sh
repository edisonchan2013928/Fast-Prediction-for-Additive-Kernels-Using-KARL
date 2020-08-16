#scale the datasets to [0,1]
#regression
#3D_spatial_network
./svm-scale -l 0 -u 1 -s ../../datasets_reg/3D_spatial_network/range ../../datasets_reg/3D_spatial_network/3D_spatial_network_train > ../../datasets_reg/3D_spatial_network/3D_spatial_network_scale
./svm-scale -r ../../datasets_reg/3D_spatial_network/range ../../datasets_reg/3D_spatial_network/3D_spatial_network_test > ../../datasets_reg/3D_spatial_network/3D_spatial_network_t_scale

#cadata
./svm-scale -l 0 -u 1 -s ../../datasets_reg/cadata/range ../../datasets_reg/cadata/cadata_train > ../../datasets_reg/cadata/cadata_scale
./svm-scale -r ../../datasets_reg/cadata/range ../../datasets_reg/cadata/cadata_test > ../../datasets_reg/cadata/cadata_t_scale

#Wave
./svm-scale -l 0 -u 1 -s ../../datasets_reg/Wave/range ../../datasets_reg/Wave/Wave_train > ../../datasets_reg/Wave/Wave_scale
./svm-scale -r ../../datasets_reg/Wave/range ../../datasets_reg/Wave/Wave_test > ../../datasets_reg/Wave/Wave_t_scale
exit

#classification 
#cod-rna
./svm-scale -l 0 -u 1 -s ../../datasets/cod-rna/range ../../datasets/cod-rna/cod-rna_train > ../../datasets/cod-rna/cod-rna_scale
./svm-scale -r ../../datasets/cod-rna/range ../../datasets/cod-rna/cod-rna_test > ../../datasets/cod-rna/cod-rna_t_scale

#home
./svm-scale -l 0 -u 1 -s ../../datasets/home/range ../../datasets/home/home_train > ../../datasets/home/home_scale
./svm-scale -r ../../datasets/home/range ../../datasets/home/home_test > ../../datasets/home/home_t_scale

#skin_non_skin
./svm-scale -l 0 -u 1 -s ../../datasets/skin_nonskin/range ../../datasets/skin_nonskin/skin_nonskin_train > ../../datasets/skin_nonskin/skin_nonskin_scale
./svm-scale -r ../../datasets/skin_nonskin/range ../../datasets/skin_nonskin/skin_nonskin_test > ../../datasets/skin_nonskin/skin_nonskin_t_scale

#a9a
./svm-scale -l 0 -u 1 -s ../../datasets/a9a/range ../../datasets/a9a/a9a > ../../datasets/a9a/a9a_scale
./svm-scale -r ../../datasets/a9a/range ../../datasets/a9a/a9a.t > ../../datasets/a9a/a9a_t_scale

#ijcnn1
./svm-scale -l 0 -u 1 -s ../../datasets/ijcnn1/range ../../datasets/ijcnn1/ijcnn1 > ../../datasets/ijcnn1/ijcnn1_scale
./svm-scale -r ../../datasets/ijcnn1/range ../../datasets/ijcnn1/ijcnn1.t > ../../datasets/ijcnn1/ijcnn1_t_scale
