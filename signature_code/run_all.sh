#!/bin/bash

##############################################################
### step 1. calibrate features for clonesig on simulations ###
##############################################################
# establish a model selection rule
# get simulated samples
rm -f 20190528_simu_cn_cancertype.csv
j=1
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for cancer_type in {0..34}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20190528_simu_cn_cancertype.csv
j=$((j+1)); 
done
done
done
done
#5250
mkdir -p logs/20190528_simu_cn_cancertype
qsub -t 1-6300 -N 20190528_simu_cn_cancertype -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190528_simu_cn_cancertype -e logs/20190528_simu_cn_cancertype -v INPUT_FILE=20190528_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype.sh

# run clonesig up to 10 clones and get estimates for different model selection criteria
rm -f 20190528_simu_cn_cancertype_run.csv
j=1
for cancer_type in {0..34}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20190528_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}-true_c_mut100" >> 20190528_simu_cn_cancertype_run.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20190528_simu_cn_cancertype_run
qsub -t 1-6300 -N 20190528_simu_cn_cancertype_run -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190528_simu_cn_cancertype_run -e logs/20190528_simu_cn_cancertype_run -v INPUT_FILE=20190528_simu_cn_cancertype_run.csv signature_code/run_clonesig_cn_cancertype.sh


# calibration of the likelihood ratio test
rm -f 20190610_loglikelihood_test_ratio.csv
j=1
for ((nb_clones=2; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for cancer_type in {0..34}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20190610_loglikelihood_test_ratio.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20190610_loglikelihood_test_ratio
qsub -t 1-5250%59 -N 20190610_loglikelihood_test_ratio -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190610_loglikelihood_test_ratio -e logs/20190610_loglikelihood_test_ratio -v INPUT_FILE=20190610_loglikelihood_test_ratio.csv signature_code/estimate_likelihood_test_statistics.sh
# analysis of the results is in the corresponding jupyter notebook

################################################################
### step 2. Benchmark to evaluate clonesig and other methods ###
################################################################
# main benchmark with varying signatures in each clone
# use previously simulated samples (for the model selection criterion)
mkdir -p logs/20190718_simu_cn_cancertype
qsub -t 1-6300%58 -N 20190718_simu_cn_cancertype -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190718_simu_cn_cancertype -e logs/20190718_simu_cn_cancertype -v INPUT_FILE=20190528_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype.sh

# run all methods and evaluation
rm -f 20190718_simu_cn_cancertype_run.csv
j=1
for cancer_type in {0..34}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20190718_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20190718_simu_cn_cancertype_run.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20190718_simu_cn_cancertype_run
qsub -t 1-6300%52 -N 20190718_simu_cn_cancertype_run -q batch -d $PWD -l walltime=3:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20190718_simu_cn_cancertype_run -e logs/20190718_simu_cn_cancertype_run -v INPUT_FILE=20190718_simu_cn_cancertype_run.csv signature_code/evaluate_clonesig_simu.sh

# compress clonesig results (that are pickled)
qsub  -N 20190726_compression -q batch -d $PWD -l walltime=10:00:00,mem=10gb,nodes=1:ppn=1 signature_code/compress_clonesig_res_objects.sh

# adding run for pyclone and 5000 mutations for the main benchmark
rm -f 20190726_big_run_pyclone.csv
j=1
nb_mut=5000
for cancer_type in {0..34}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20190718_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20190726_big_run_pyclone.csv
j=$((j+1)); 
done
done
done

mkdir -p logs/20190726_big_run_pyclone
qsub -t 1-1260%100 -N 20190726_big_run_pyclone -q batch -d $PWD -l walltime=72:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20190726_big_run_pyclone -e logs/20190721_simu_cn_cancertype_run -v INPUT_FILE=20190726_big_run_pyclone.csv signature_code/run_pyclone.sh


# with constant signature activity across clones (without pyclone 5000 mutations, too long)
# get simulated samples
rm -f 20190729_simu_cn_cancertype.csv
j=1
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for cancer_type in {0..10}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20190729_simu_cn_cancertype.csv
j=$((j+1)); 
done
done
done
done

mkdir -p logs/20190729_simu_cn_cancertype_cst
qsub -t 1-1980%98 -N 20190729_simu_cn_cancertype_cst -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190729_simu_cn_cancertype_cst -e logs/20190729_simu_cn_cancertype_cst -v INPUT_FILE=20190729_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype_cst.sh

# run all methods and evaluation
rm -f 20190729_simu_cn_cancertype_run_cst.csv
j=1
for cancer_type in {0..10}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20190729_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20190729_simu_cn_cancertype_run_cst.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20190729_simu_cn_cancertype_run_cst
qsub -t 1-1980%98 -N 20190729_simu_cn_cancertype_run_cst -q batch -d $PWD -l walltime=72:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20190729_simu_cn_cancertype_run_cst -e logs/20190729_simu_cn_cancertype_run_cst -v INPUT_FILE=20190729_simu_cn_cancertype_run_cst.csv signature_code/evaluate_clonesig_simu.sh


############################################################
### step 3. Special evaluations of clonesig performances ###
############################################################
# simulations and run for assessing the separating power of clonesig
qsub  -N 20190803_simus -q batch -d $PWD -l walltime=10:00:00,mem=10gb,nodes=1:ppn=1 signature_code/simulations_power_distinction_clonesig.sh
rm -f 20190803_clonesig_power_eval.csv
j=1
for pi in {0..29}; do
for phi in {0..9}; do
for mut in {100,300,1000}; do
for perc_dip in 0.1 0.5 0.9; do
for depth in {100,500}; do
echo $j,"20190803_simulations_eval_clonesig_power/pi${pi}-phi${phi}-depth${depth}-percdip${perc_dip}-nb_mut${mut}" >> 20190803_clonesig_power_eval.csv
j=$((j+1)); 
done
done
done
done
done
mkdir -p logs/20190803_clonesig_power_eval
qsub -t 1-5400%98 -N 20190803_clonesig_power_eval -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20190803_clonesig_power_eval -e logs/20190803_clonesig_power_eval -v INPUT_FILE=20190803_clonesig_power_eval.csv signature_code/simulations_power_distinction_clonesig_run.sh


# simulations and run for evluating the sensitivity of the statistical test
rm -f 20190804_simus_sensitivity.csv
j=1
for nb_clones in {2..6}; do
echo $j,"${nb_clones}" >> 20190804_simus_sensitivity.csv
j=$((j+1)); 
done
mkdir -p logs/20190804_simus_sensitivity
qsub -t 1-5 -N 20190804_simus_sensitivity -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20190804_simus_sensitivity -e logs/20190804_simus_sensitivity -v INPUT_FILE=20190804_simus_sensitivity.csv signature_code/simulations_sensitivity_statistical_test.sh
rm -f 20190804_clonesig_power_eval.csv
j=1
for pi in {0..29}; do
for phi in {2..6}; do
for mut in {100,300,1000}; do
for perc_dip in 0.1 0.5 0.9; do
for depth in {100,500}; do
echo $j,"20190804_simulations_eval_clonesig_power/pi${pi}-phi${phi}-depth${depth}-percdip${perc_dip}-nb_mut${mut}" >> 20190804_clonesig_power_eval.csv
j=$((j+1)); 
done
done
done
done
done
mkdir -p logs/20190804_clonesig_power_eval
qsub -t 1-2700%118 -N 20190804_clonesig_power_eval -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20190804_clonesig_power_eval -e logs/20190804_clonesig_power_eval -v INPUT_FILE=20190804_clonesig_power_eval.csv signature_code/simulations_power_distinction_clonesig_run.sh


###############################################
### step 4. Download and organize TCGA data ###
###############################################
# download TCGA data (SNVs)
./signature_code/organise_dl_data_tcga.sh
# merge maf files, and organize
rm -f 20190629_get_snv_filters.csv
j=1
for cancer_loc in ACC BRCA CHOL DLBC GBM KICH KIRP LGG LUAD MESO PAAD PRAD SARC STAD THCA UCEC UVM BLCA CESC COAD ESCA HNSC KIRC LAML LIHC LUSC OV PCPG READ SKCM TGCT THYM UCS ; do
echo $j,$cancer_loc,$pythonenv >> 20190629_get_snv_filters.csv
j=$((j+1)); 
done
mkdir -p logs/20190629_tcga
qsub -t 1-33%1 -N 20190629_get_snv_filters -q batch -d $PWD -l walltime=10:00:00,mem=60gb,nodes=1:ppn=1 -o logs/20190629_tcga -e logs/20190629_tcga -v INPUT_FILE=20190629_get_snv_filters.csv ./signature_code/get_SNV_filters_and_data_exploration.sh
./merge_coad_read.py # because of clinical data
# build final patient list with complete info
qsub -t 1-33%1 -N get_patient_list -q batch -d $PWD -l walltime=1:00:00,mem=10gb,nodes=1:ppn=1 -o logs/20190629_tcga -e logs/20190629_tcga -v INPUT_FILE=20190629_get_snv_filters.csv ./signature_code/get_patient_list.sh

################################################
### step 5. run clonesig on the TCGA samples ###
################################################
# get curated cancer type-signature match table
./signature_code/curated_match_tcga_match_sig_list.py

# launch clonesig
mkdir -p logs/20190811_clonesig_TCGA
qsub -t 1-8958%48 -N 20190811_clonesig_TCGA -q batch -d $PWD -l walltime=3:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20190811_clonesig_TCGA -e logs/20190811_clonesig_TCGA -v INPUT_FILE=TCGA_pancancer_patient_list.csv signature_code/run_clonesig_tcga.sh

# plot results per sample in the cases with a significant change
mkdir -p logs/20190828_plot_TCGA
qsub -t 1-2836%58 -N 20190828_plot_TCGA -q batch -d $PWD -l walltime=3:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20190828_plot_TCGA -e logs/20190828_plot_TCGA -v INPUT_FILE=20190828_tcga_plots.csv signature_code/plot_tcga_restr.sh

#################################
### step 6. get result tables ###
#################################
# the following script was actually run by part along the steps of the analyses
# as some results are needed to adjust the parameters for the next step
# of analysis (model selection criterion and statisitcal test calibration)
./signature_code/get_result_tables.py

