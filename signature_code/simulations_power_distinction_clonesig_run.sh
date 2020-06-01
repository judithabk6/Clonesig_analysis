#!/bin/bash

i=1
for PARAM in folder_path nb_mut
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


# run tracksig
export PATH=/bioinfo/local/build/R/R-3.3.2_centos/bin:$PATH
./signature_code/run_tracksig_power.R -f $folder_path

# run tracksigfreq
./signature_code/run_tracksigfreq_power.R -f $folder_path

# run palimpsest
./signature_code/run_palimpsest_power.R -f $folder_path

./signature_code/simulations_power_distinction_clonesig_run.py $folder_path
