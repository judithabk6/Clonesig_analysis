#!/bin/bash

i=1
for PARAM in patient_id cancer_loc cohort
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


grep "$patient_id" pcawg_download/snv/final_consensus_passonly.snv_mnv_indel.$cohort.public.maf > 20200510_pcawg/$patient_id.maf

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


./signature_code/run_clonesig_pcawg.py $patient_id $cancer_loc $cohort
bzip2 20200510_pcawg/*raw_results_restr_curated