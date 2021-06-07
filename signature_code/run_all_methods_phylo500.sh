#!/bin/bash

i=1
for PARAM in samplename
# samplename="Sim_500_19"
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate

source activate my_python36_centos
./signature_code/get_input_phylo500.py $samplename
folder_path_cst=PhylogicNDT500/$samplename\_cst
folder_path_var=PhylogicNDT500/$samplename\_var

./signature_code/run_clonesig_phylo500.py $folder_path_cst
./signature_code/run_clonesig_phylo500.py $folder_path_var
source deactivate

# run sciclone
export PATH=/bioinfo/local/build/R/R-3.3.2_centos/bin:$PATH
var_sciclone_1=`date +%s`
timeout 172800 ./signature_code/run_sciclone.R -f $folder_path_cst
var_sciclone_2=`date +%s`
echo $var_sciclone_1,$var_sciclone_2 > $folder_path_cst/sciclone_timing.txt

# run ccube
var_ccube_1=`date +%s`
timeout 172800 R --slave --no-restore --file=./signature_code/run_ccube.R --args -f  $folder_path_cst
var_ccube_2=`date +%s`
echo $var_ccube_1,$var_ccube_2 > $folder_path_cst/ccube_timing.txt

# run deconstructsig
timeout 172800 ./signature_code/run_deconstructsigs_phylo500.R -f $folder_path_cst
timeout 172800 ./signature_code/run_deconstructsigs_phylo500.R -f $folder_path_var

# run tracksig
#if [ $nb_mut -gt 599 ]; then
timeout 172800 ./signature_code/run_tracksig_phylo500.R -f $folder_path_cst
timeout 172800 ./signature_code/run_tracksig_phylo500.R -f $folder_path_var
#fi
# run tracksigfreq
timeout 172800 ./signature_code/run_tracksigfreq_phylo500.R -f $folder_path_cst
timeout 172800 ./signature_code/run_tracksigfreq_phylo500.R -f $folder_path_var

# run palimpsest
timeout 172800 ./signature_code/run_palimpsest_phylo500.R -f $folder_path_cst
timeout 172800 ./signature_code/run_palimpsest_phylo500.R -f $folder_path_var

# run dpclust
var_dpclust_1=`date +%s`
timeout 172800 R --vanilla --slave -q -f /bioinfo/users/jabecass/dl_tools_centos/dpclust/inst/example/dpclust_pipeline.R --args -r 1 -d $folder_path_cst/dpclust -o $folder_path_cst/dpclust -i $folder_path_cst/dpclust/info.tsv
var_dpclust_2=`date +%s`
echo $var_dpclust_1,$var_dpclust_2 > $folder_path_cst/dpclust_timing.txt

# run phylogicNDT
source activate my_python27_PhylogicNDT
purity=`cat $folder_path_cst/purity.txt`
var_phylogicndt_1=`date +%s`
cd $folder_path_cst/phylogicndt
timeout 172800 /data/users/jabecass/dl_tools_centos/PhylogicNDT/PhylogicNDT.py Cluster -i Test_Clust -s sample_01:input.maf::$purity:1 --maf_input_type calc_ccf
rm -r Test_Clust.phylogic_report.html Test_Clust_1d_mutation_plots Test_Clust_1d_cluster_plots
cd ../../..
var_phylogicndt_2=`date +%s`
source deactivate
echo $var_phylogicndt_1,$var_phylogicndt_2 > $folder_path_cst/phylogicndt_timing.txt

# run pyclone
source activate my_python27_centos
purity=`cat $folder_path_cst/purity.txt`
var_pyclone_1=`date +%s`
timeout 172800 PyClone run_analysis_pipeline --in_files $folder_path_cst/pyclone/input.tsv --working_dir $folder_path_cst/pyclone --tumour_contents $purity 
var_pyclone_2=`date +%s`
source deactivate
echo $var_pyclone_1,$var_pyclone_2 > $folder_path_cst/pyclone_timing.txt

source activate my_python36_centos
./signature_code/evaluate_phylo500.py $folder_path_cst
./signature_code/evaluate_phylo500.py $folder_path_var
source deactivate


