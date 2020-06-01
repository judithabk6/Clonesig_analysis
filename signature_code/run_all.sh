#!/bin/bash

##############################################################
### step 1. calibrate features for clonesig on simulations ###
##############################################################
# establish a model selection rule
# get simulated samples
rm -f 20200213_get_simu_cn_cancertype.csv
j=1
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for cancer_type in {0..30}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20200213_get_simu_cn_cancertype.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200213_get_simu_cn_cancertype
qsub -t 1-7812 -N 20200213_get_simu_cn_cancertype -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20200213_get_simu_cn_cancertype -e logs/20200213_get_simu_cn_cancertype -v INPUT_FILE=20200213_get_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype.sh


# run clonesig up to 10 clones and get estimates for different model selection criteria
rm -f 20200213_simu_cn_cancertype_run.csv
j=1
for cancer_type in {0..30}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20200210_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20200213_simu_cn_cancertype_run.csv20200213_simu_cn_cancertype_run.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200213_simu_cn_cancertype_run_tmp
qsub -t 8-7812%57 -N 20200213_simu_cn_cancertype_run -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20200213_simu_cn_cancertype_run_tmp -e logs/20200213_simu_cn_cancertype_run_tmp -v INPUT_FILE=20200213_simu_cn_cancertype_run.csv signature_code/run_clonesig_cn_cancertype.sh

rm -f 20200213_simu_cn_cancertype_run_big.csv
j=1
for cancer_type in {0..30}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for perc_dip in {0,20,40,60,80,100} ; do
nb_mut=5000
echo $j,"20200210_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20200213_simu_cn_cancertype_run_big.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200213_simu_cn_cancertype_run_big
qsub -t 1-1116%57 -N 20200213_simu_cn_cancertype_run_big -q batch -d $PWD -l walltime=01:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20200213_simu_cn_cancertype_run_big -e logs/20200213_simu_cn_cancertype_run_big -v INPUT_FILE=20200213_simu_cn_cancertype_run_big.csv signature_code/run_clonesig_cn_cancertype.sh



# calibration of the likelihood ratio test
rm -f 20200217_loglikelihood_test_ratio.csv
j=1
for ((nb_clones=2; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for cancer_type in {0..30}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20200217_loglikelihood_test_ratio.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200217_loglikelihood_test_ratio
qsub -t 1-7812%57 -N 20200217_loglikelihood_test_ratio -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20200217_loglikelihood_test_ratio -e logs/20200217_loglikelihood_test_ratio -v INPUT_FILE=20200217_loglikelihood_test_ratio.csv signature_code/estimate_likelihood_test_statistics.sh
# analysis of the results is in the corresponding jupyter notebook

################################################################
### step 2. Benchmark to evaluate clonesig and other methods ###
################################################################
# main benchmark with varying signatures in each clone
# use previously simulated samples (for the model selection criterion)
mkdir -p logs/20190718_simu_cn_cancertype
qsub -t 1-6300%58 -N 20190718_simu_cn_cancertype -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20190718_simu_cn_cancertype -e logs/20190718_simu_cn_cancertype -v INPUT_FILE=20190528_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype.sh

# run all methods and evaluation
rm -f 20200218_simu_cn_cancertype_run.csv
j=1
for cancer_type in {0..30}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20200210_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20200218_simu_cn_cancertype_run.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200218_simu_cn_cancertype_run
qsub -t 1-7812%57 -N 20200218_simu_cn_cancertype_run -q batch -d $PWD -l walltime=3:00:00,mem=4gb,nodes=1:ppn=1 -o logs/20200218_simu_cn_cancertype_run -e logs/20200218_simu_cn_cancertype_run -v INPUT_FILE=20200218_simu_cn_cancertype_run.csv signature_code/evaluate_clonesig_simu.sh

# compress clonesig results (that are pickled)
qsub  -N 20190726_compression -q batch -d $PWD -l walltime=10:00:00,mem=10gb,nodes=1:ppn=1 signature_code/compress_clonesig_res_objects.sh


# with constant signature activity across clones (without pyclone 5000 mutations, too long)
# get simulated samples
rm -f 20200229_simu_cn_cancertype.csv
j=1
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for cancer_type in {0..10}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"${nb_clones},${nb_mut},${cancer_type},${perc_dip},${j}" >> 20200229_simu_cn_cancertype.csv
j=$((j+1)); 
done
done
done
done

mkdir -p logs/20200326_simu_cn_cancertype_cst
qsub -t 1-2772%18 -N 20200326_simu_cn_cancertype_cst -q batch -d $PWD -l walltime=01:00:00,mem=1gb,nodes=1:ppn=1 -o logs/20200326_simu_cn_cancertype_cst -e logs/20200326_simu_cn_cancertype_cst -v INPUT_FILE=20200229_simu_cn_cancertype.csv signature_code/simulations_clonesig_cn_cancertype_cst.sh

# adding run for pyclone and 5000 mutations for the main benchmark
rm -f 2020308_big_run_pyclone.csv
j=1
nb_mut=5000
for cancer_type in {0..30}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20200210_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 2020308_big_run_pyclone.csv
j=$((j+1)); 
done
done
done

mkdir -p logs/2020308_big_run_pyclone
qsub -t 1-1116%100 -N 2020308_big_run_pyclone -q batch -d $PWD -l walltime=96:00:00,mem=5gb,nodes=1:ppn=1 -o logs/2020308_big_run_pyclone -e logs/2020308_big_run_pyclone -v INPUT_FILE=2020308_big_run_pyclone.csv signature_code/run_pyclone.sh


# run all methods and evaluation
rm -f 20200229_simu_cn_cancertype_run_cst.csv
j=1
for cancer_type in {0..10}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for nb_mut in {20,50,100,300,600,1000,5000}; do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20200229_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 20200229_simu_cn_cancertype_run_cst.csv
j=$((j+1)); 
done
done
done
done
mkdir -p logs/20200229_simu_cn_cancertype_run_cst
qsub -t 1-2772%117 -N 20200229_simu_cn_cancertype_run_cst -q batch -d $PWD -l walltime=72:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20200229_simu_cn_cancertype_run_cst -e logs/20200229_simu_cn_cancertype_run_cst -v INPUT_FILE=20200229_simu_cn_cancertype_run_cst.csv signature_code/evaluate_clonesig_simu.sh

rm -f 2020420_big_run_pyclone_cst.csv
j=1
nb_mut=5000
for cancer_type in {0..10}; do
for ((nb_clones=1; nb_clones<=6; nb_clones++)); do
for perc_dip in {0,20,40,60,80,100} ; do
echo $j,"20200229_simulations_clonesig_cn_cancer_type/type${cancer_type}-perc_diploid${perc_dip}-nb_clones${nb_clones}-nb_mut${nb_mut}" >> 2020420_big_run_pyclone_cst.csv
j=$((j+1)); 
done
done
done

mkdir -p logs/2020420_big_run_pyclone_cst
qsub -t 1-396%90 -N 2020420_big_run_pyclone_cst -q batch -d $PWD -l walltime=96:00:00,mem=5gb,nodes=1:ppn=1 -o logs/2020420_big_run_pyclone_cst -e logs/2020308_big_run_pyclone -v INPUT_FILE=2020308_big_run_pyclone.csv signature_code/run_pyclone.sh

# evaluation on the dream dataset
rm -f 20200415_dream_run.csv
j=1
for i in {2..6}; do
for depth in {8,16,32,64,128}; do
echo $j,"T${i},${depth}X" >> 20200415_dream_run.csv
j=$((j+1)); 
done
done
mkdir -p logs/20200415_dream_run
qsub -t 1-25 -N 20200415_dream_run -q batch -d $PWD -l walltime=120:00:00,mem=10gb,nodes=1:ppn=1 -o logs/20200415_dream_run -e logs/20200415_dream_run -v INPUT_FILE=20200415_dream_run.csv signature_code/run_all_methods_dream.sh

############################################################
### step 3. Special evaluations of clonesig performances ###
############################################################
# simulations and run for assessing the separating power of clonesig
qsub  -N 20200305_simus -q batch -d $PWD -l walltime=10:00:00,mem=10gb,nodes=1:ppn=1 signature_code/simulations_power_distinction_clonesig.sh
rm -f 20200305_clonesig_power_eval.csv
j=1
for pi in {0..49}; do
for phi in {0..9}; do
for mut in {30,100,300,1000}; do
for perc_dip in 0.1 0.5 0.9; do
for depth in {100,500}; do
echo $j,"20200305_simulations_eval_clonesig_power/pi${pi}-phi${phi}-depth${depth}-percdip${perc_dip}-nb_mut${mut}" >> 20200305_clonesig_power_eval.csv
j=$((j+1)); 
done
done
done
done
done
mkdir -p logs/20200305_clonesig_power_eval
qsub -t 1-12000%117 -N 20200305_clonesig_power_eval -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20200305_clonesig_power_eval -e logs/20200305_clonesig_power_eval -v INPUT_FILE=20200305_clonesig_power_eval.csv signature_code/simulations_power_distinction_clonesig_run.sh
qsub -t 5001-12000%117 -N 20200305_clonesig_power_eval -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20200305_clonesig_power_eval -e logs/20200305_clonesig_power_eval -v INPUT_FILE=20200305_clonesig_power_eval.csv signature_code/simulations_power_distinction_clonesig_run.sh


# simulations and run for evluating the sensitivity of the statistical test
rm -f 20200430_simus_sensitivity.csv
j=1
for nb_clones in {2..6}; do
echo $j,"${nb_clones}" >> 20200430_simus_sensitivity.csv
j=$((j+1)); 
done
mkdir -p logs/20200430_simus_sensitivity
qsub -t 1-5 -N 20200430_simus_sensitivity -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20200430_simus_sensitivity -e logs/20200430_simus_sensitivity -v INPUT_FILE=20200430_simus_sensitivity.csv signature_code/simulations_sensitivity_statistical_test.sh
rm -f 20200501_clonesig_power_eval.csv
j=1
for pi in {0..29}; do
for phi in {2..6}; do
for mut in {30,100,300,1000}; do
for perc_dip in 0.1 0.5 0.9; do
for depth in {100,500}; do
echo $j,"20200430_simulations_eval_clonesig_power/pi${pi}-phi${phi}-depth${depth}-percdip${perc_dip}-nb_mut${mut}" >> 20200501_clonesig_power_eval.csv
j=$((j+1)); 
done
done
done
done
done
mkdir -p logs/20200501_clonesig_power_eval
qsub -t 1-3600%100 -N 20200501_clonesig_power_eval -q batch -d $PWD -l walltime=01:00:00,mem=3gb,nodes=1:ppn=1 -o logs/20200501_clonesig_power_eval -e logs/20200501_clonesig_power_eval -v INPUT_FILE=20200501_clonesig_power_eval.csv signature_code/simulations_power_distinction_clonesig_run.sh


#########################################################
### step 4. Download and organize TCGA and PCAWG data ###
#########################################################
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

##########################################################
### step 5. run clonesig on the TCGA and PCAWG samples ###
##########################################################
# get curated cancer type-signature match table
./signature_code/curated_match_tcga_match_sig_list.py

# launch clonesig
mkdir -p logs/20200503_clonesig_TCGA
qsub -t 1-8958%118 -N 20200503_clonesig_TCGA -q batch -d $PWD -l walltime=3:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20200503_clonesig_TCGA -e logs/20200503_clonesig_TCGA -v INPUT_FILE=TCGA_pancancer_patient_list.csv signature_code/run_clonesig_tcga.sh

# plot results per sample in the cases with a significant change
mkdir -p logs/20200519_plot_TCGA
qsub -t 1-8958%10 -N 20200519_plot_TCGA -q batch -d $PWD -l walltime=3:00:00,mem=5gb,nodes=1:ppn=1 -o logs/20200519_plot_TCGA -e logs/20200519_plot_TCGA -v INPUT_FILE=TCGA_pancancer_patient_list.csv signature_code/plot_tcga_restr.sh

mkdir -p logs/20200510_clonesig_pcawg
qsub -t 1-1950%117 -N 20200510_clonesig_pcawg -q batch -d $PWD -l walltime=6:00:00,mem=8gb,nodes=1:ppn=1 -o logs/20200510_clonesig_pcawg -e logs/20200510_clonesig_pcawg -v INPUT_FILE=20200510_pcawgs_final_patient_list.csv signature_code/run_clonesig_pcawg.sh

mkdir -p logs/20200510_clonesig_pcawg_tcga
qsub -t 5-828%85 -N 20200510_clonesig_pcawg_tcga -q batch -d $PWD -l walltime=6:00:00,mem=8gb,nodes=1:ppn=1 -o logs/20200510_clonesig_pcawg_tcga -e logs/20200510_clonesig_pcawg_tcga -v INPUT_FILE=20200510_pcawgs_tcga_final_patient_list.csv signature_code/run_clonesig_pcawg.sh


#################################
### step 6. get result tables ###
#################################
# the following script was actually run by part along the steps of the analyses
# as some results are needed to adjust the parameters for the next step
# of analysis (model selection criterion and statisitcal test calibration)
./signature_code/get_result_tables.py

