#!/bin/bash

i=1
for PARAM in samplename
# samplename="Sim_500_19"
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate


folder_path_cst=PhylogicNDT500/$samplename\_cst
folder_path_var=PhylogicNDT500/$samplename\_var


source activate my_python36_centos
./signature_code/evaluate_phylo500.py $folder_path_cst
./signature_code/evaluate_phylo500.py $folder_path_var
source deactivate


