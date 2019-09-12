#!/bin/bash

i=1
for PARAM in nb_clones nb_mut cancer_type perc_dip nb_seed
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


./signature_code/estimate_likelihood_test_statistics.py $nb_clones $nb_mut $cancer_type $perc_dip $nb_seed
