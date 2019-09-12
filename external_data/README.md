# External data
### TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv
downloaded from Synapse on July 8th (Alexandrov et al. 2018 supplementary material, syn11801497)

### sigProfiler_SBS_signatures_2018_03_28.csv and sigProfiler_exome_SBS_signatures.csv
downloaded from Synapse on May 2nd (Alexandrov et al. 2018 supplementary material, first file is not available anymore (syn11738319), and second can be found at syn11967913)

### signatures_probabilities.txt
downloaded from COSMIC and corresponds to the V2 of signatures.

### match_cancer_type_sig_v3.csv
obtained from the SigProfiler package using script ```Clonesig_analysis/signature_code/sigprofiler_alexandrov18_data.py```

### gdc_manifest.2019-05-19.txt
obtained from manually filled cart on the GDC portal, to download TCGA data

### curated_match_signature_cancertype_tcgawes_literature.csv
obtained using the file ```TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv``` and manual review of the literature, as detailed in the script ```Clonesig_analysis/signature_code/curated_match_tcga_match_sig_list.py```
