## Processed Results


#### 20200222_clonesig_model_selection.csv
The objective of this table is to calibrate a statistical test to assess the significance of a change in signatures between different clones. To that end, we compare the likelihood of a Clonesig model with $J$ clones as determined by the model selection criterion (loglikelihood column), and the likelihood of a model with the same clones but a single mixture of signatures common to all the clones (and found by fitting all observed mutations together) (loglikelihood_nopi column). The objective of the test is to determine whether the difference between the two likelihoods is significant. To that end, we have implemented a Likelihood-ratio test, and use this data to calibreate the parameter. MU is the full MU matrix with 49 signatures, and subMU only the subset corresponding to the given cancer type. (see the notebook [Clonesig_analysis/paper_notebooks/statistical_test_calibration.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/statistical_test_calibration.ipynb) for more details)

This table was used to generate Supplementary Figures 5, 6 and 7.


#### 20200222_clonesig_model_selection_criterion.csv
Table with results of CloneSig on simulated samples. Each row is a simulated sample (with varying characteristics cancer_type, perc_diploid, nb_clones, nb_mut, etc), and the result from a variety of criteria (bic, aic etc), and the likelihood of CloneSig model for a number of clones ranging from 1 to 10 (see the notebook [Clonesig_analysis/paper_notebooks/model_selection_criterion_learning.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/model_selection_criterion_learning.ipynb) for more details)

This table was used to generate Supplementary Figures 2, 3 and 4.



#### 20200512_eval_clonesig_power_statistical_test_sensitivity.csv
The objective of this table was to evaluate the sensitivity of the statistical test developed to assess the significance of signature change between subclones, and identify important variables influencing the sensitivity. Each row corresponds to a simulated sample, with some characteristics (number of clones, mutations, cosine distance between signature activities etc), and metrics associated with CloneSig's results (fitted number of clones, pvalue of the test etc) (see the notebook [Clonesig_analysis/paper_notebooks/sensitivity_eval_test.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/sensitivity_eval_test.ipynb) for more details)

This table was used to generate Supplementary Figures 8 and 9.



#### 20200513_eval_clonesig_power_2clones.csv
The objective of this table is to evaluate the ability of methods resolving both ITH structure and signature activity deconvolution (CloneSig, TrackSig, TrackSigFreq and Palimpsest) to distinguish clones with close frequencies depending on their difference in signature activity (see the notebook [Clonesig_analysis/paper_notebooks/clonesig_separation_power.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/clonesig_separation_power.ipynb) for more details). 

This table was used to generate Figure 4, and Supplementary Figures 110, 111 and 112.



#### 20200514_pcawg_results.csv
contains CloneSig results for each PCAWG sample (number of clones, statistic and p-value for the significance of signature activity change between subclones, some characteristics of clonal mutations, and the largest subclone in term of number of mutations, and signature activity for clonal mutations and the largest subclone)

This table was used to generate Supplementary Figures 69 to 102.



#### 20200514_tcga_results.csv
contains CloneSig results for each TCGA sample (number of clones, statistic and p-value for the significance of signature activity change between subclones, some characteristics of clonal mutations, and the largest subclone in term of number of mutations, and signature activity for clonal mutations and the largest subclone) (see the notebook [Clonesig_analysis/paper_notebooks/TCGA_results.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/TCGA_results.ipynb) for more details)

This table was used to generate Figure 6, and Supplementary Figures 36 to 68.


#### 20200214_tcga_results_survival_restr.csv
contains CloneSig results for each TCGA sample (number of clones, statistic and p-value for the significance of signature activity change between subclones, some characteristics of clonal mutations, and the largest subclone in term of number of mutations, and signature activity for clonal mutations and the largest subclone, also reported in table ```20200514_tcga_results.csv```), and some clinical variables (survival, age, sex, tumor stage) (see the notebook [Clonesig_analysis/paper_notebooks/TCGA_results.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/TCGA_results.ipynb) for more details)

This table was used to generate Figure 6, and Supplementary Figures 36 to 68.


#### 20200525_dream_results.csv
Benchmark results for the DREAM dataset (see the notebook [Clonesig_analysis/paper_notebooks/dream_results.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/dream_results.ipynb) for more details)

This table was used to generate Figure 2, and Supplementary Figures 28 to 33.

#### 20201011_phylo500_results.csv
Benchmark results for the PhylogicSim500 dataset (see the notebook [Clonesig_analysis/paper_notebooks/phylo500_results.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/phylo500_results.ipynb) for more details)

This table was used to generate Figure 2, 3 and 5.


#### 20201011_simclone1000_results.csv
Benchmark results for the PhylogicSim500 dataset (see the notebook [Clonesig_analysis/paper_notebooks/SimClone100_results.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/SimClone100_results.ipynb) for more details)

This table was used to generate Figure 2, and Supplementary Figures 34 and 35.

#### 20200520_eval_compare_simulations_new.csv and 20200520_eval_compare_simulations_cst.csv
Benchmark results for the CloneSigSim dataset for the varying and the constant settings (see the notebook [Clonesig_analysis/paper_notebooks/performance_analysis_simulations.ipynb](https://github.com/judithabk6/Clonesig_analysis/blob/master/paper_notebooks/performance_analysis_simulations.ipynb) for more details)

Those tables were used to generate Figure 2, and Supplementary Figures 10 to 27.