#!/bin/bash

i=1
for PARAM in folder_path
# folder_path=20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

nb_mut=`echo $folder_path | cut -d "t" -f 5`

source /data/users/jabecass/dl_tools_centos/miniconda3_activate

# run pyclone
if [ $nb_mut -lt 1001 ]; then 
source activate my_python27_centos
purity=`cat $folder_path/purity.txt`
var_pyclone_1=`date +%s`
PyClone run_analysis_pipeline --in_files $folder_path/pyclone/input.tsv --working_dir $folder_path/pyclone --tumour_contents $purity 
var_pyclone_2=`date +%s`
source deactivate
fi

# run sciclone
export PATH=/bioinfo/local/build/R/R-3.3.2_centos/bin:$PATH
var_sciclone_1=`date +%s`
./signature_code/run_sciclone.R -f $folder_path
var_sciclone_2=`date +%s`

# run ccube
var_ccube_1=`date +%s`
R --slave --no-restore --file=./signature_code/run_ccube.R --args -f  $folder_path
var_ccube_2=`date +%s`

# run deconstructsig
./signature_code/run_deconstructsigs.R -f $folder_path

# run tracksig
#if [ $nb_mut -gt 599 ]; then
./signature_code/run_tracksig.R -f $folder_path
#fi
# run palimpsest
./signature_code/run_palimpsest.R -f $folder_path


source activate my_python36_centos
./signature_code/evaluate_pyclone.py $folder_path $var_pyclone_1 $var_pyclone_2
./signature_code/evaluate_sciclone.py $folder_path $var_sciclone_1 $var_sciclone_2
./signature_code/evaluate_ccube.py $folder_path $var_ccube_1 $var_ccube_2
./signature_code/evaluate_deconstructsig.py $folder_path
./signature_code/evaluate_tracksig.py $folder_path
./signature_code/evaluate_palimpsest.py $folder_path
./signature_code/evaluate_clonesig_simu.py $folder_path
