#!/bin/bash

i=1
for PARAM in patient_id cancer_loc cosmic_type
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

mkdir -p 20190704_TCGA/$patient_id

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


./signature_code/plot_tcga_restr.py $patient_id $cancer_loc $cosmic_type
