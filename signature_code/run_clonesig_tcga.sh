#!/bin/bash

i=1
for PARAM in patient_id cancer_loc cosmic_type
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

mkdir -p 20190704_TCGA/$patient_id

grep "$patient_id" tmp/$cancer_loc\_useful_final_public_merge_cnv_purity.csv > 20190704_TCGA/$patient_id/public_snv.tsv
grep "$patient_id" tmp/$cancer_loc\_useful_final_qc_merge_cnv_purity.csv > 20190704_TCGA/$patient_id/protected_snv.tsv
grep "$patient_id" external_data/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv > 20190704_TCGA/$patient_id/sigprofiler_sig.tsv

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos


./signature_code/run_clonesig_tcga.py $patient_id $cancer_loc $cosmic_type
bzip2 20190704_TCGA/$patient_id/*raw_results_restr_curated