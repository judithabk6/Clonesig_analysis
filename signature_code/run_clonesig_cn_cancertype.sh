#!/bin/bash

i=1
for PARAM in folder_path
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


./signature_code/run_clonesig_cn_cancertype.py $folder_path
