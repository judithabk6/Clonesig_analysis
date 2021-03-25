#!/bin/bash

i=1
for PARAM in folder_path
# samplename="Sim_500_19"
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate

source activate my_python36_centos
./signature_code/fix_evaluate.py $folder_path
source deactivate
