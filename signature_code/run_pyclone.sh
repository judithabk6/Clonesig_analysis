#!/bin/bash

i=1
for PARAM in folder_path
# folder_path=20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done



source /data/users/jabecass/dl_tools_centos/miniconda3_activate


source activate my_python27_centos
purity=`cat $folder_path/purity.txt`
var_pyclone_1=`date +%s`
PyClone run_analysis_pipeline --in_files $folder_path/pyclone/input.tsv --working_dir $folder_path/pyclone --tumour_contents $purity 
var_pyclone_2=`date +%s`
source deactivate




source activate my_python36_centos
./signature_code/evaluate_pyclone.py $folder_path $var_pyclone_1 $var_pyclone_2


