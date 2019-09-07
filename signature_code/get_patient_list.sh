#!/bin/bash

i=1
for PARAM in cancer_loc pythonenv
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos

./signature_code/get_patient_list.py $cancer_loc
