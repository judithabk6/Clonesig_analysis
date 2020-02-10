#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import os

try:
    rows, columns = os.popen('stty size', 'r').read().split()
    pd.set_option('display.width', int(columns))
    pd.set_option('display.max_columns', 600)
    pd.set_option('display.max_rows', 1000)
except:
    print("running on server, otherwise please investigate")

tcga = pd.read_csv('external_data/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv')
tcga = tcga.assign(patient_id=tcga['Sample Names'].str[:12])
patient_list = pd.read_csv('TCGA_pancancer_patient_list.csv', sep=',', index_col=0, header=None, names=['patient_id', 'tcga_name', 'pcawg_name'])
sig_cols = [c for c in tcga.columns if 'SBS' in c]
tcga_m = pd.merge(tcga, patient_list, on='patient_id', how='left')
tcga_m[sig_cols] = tcga_m[sig_cols].astype(bool)
final_table = (tcga_m.groupby('tcga_name')[sig_cols].sum() /tcga_m.groupby('tcga_name').patient_id.count()[:, np.newaxis] ).T
final_table_filter = (final_table>0.03).astype(int)

modif_dict = {'BRCA': ['SBS40', 'SBS18', 'SBS6', 'SBS20', 'SBS26', 'SBS17a', 'SBS17b', 'SBS8'],
              'CHOL': ['SBS6', 'SBS8', 'SBS16', 'SBS17a', 'SBS17b', 'SBS20', 'SBS22'],
              'COADREAD': ['SBS18', 'SBS6', 'SBS26'],
              'DLBC': ['SBS18', 'SBS9'],
              'GBM': ['SBS8'],
              'HNSC': ['SBS3', 'SBS18'],
              'KICH': ['SBS13'],
              'KIRC': ['SBS2', 'SBS3', 'SBS9', 'SBS13', 'SBS17a', 'SBS17b'],
              'LIHC': ['SBS3', 'SBS4', 'SBS6', 'SBS12', 'SBS17a', 'SBS17b', 'SBS18', 'SBS23', 'SBS35'],
              'LUSC': ['SBS18'],
              'MESO': ['SBS13'],
              'OV': ['SBS8', 'SBS16', 'SBS39'],
              'PAAD': ['SBS3', 'SBS17a', 'SBS17b', 'SBS18', 'SBS36'],
              'PCPG': ['SBS3'],
              'PRAD': ['SBS3', 'SBS8', 'SBS9', 'SBS16', 'SBS18', 'SBS34', 'SBS37'],
              'STAD': ['SBS3', 'SBS18'],
              'THCA': ['SBS18'],
              'UCEC': ['SBS3', 'SBS6', 'SBS18'],
              'UVM': ['SBS6', 'SBS12', 'SBS16']
              }

for cancer_loc, sig_list in modif_dict.items():
    for sig in sig_list:
        final_table_filter.loc[sig, cancer_loc] = 1

# exclude artefactual signatures SBS27 SBS43 SBS45 SBS46 SBS47 SBS48 SBS49 SBS50 SBS51 SBS52 SBS53 SBS54 SBS55 SBS56 SBS57 SBS58 SBS59 SBS60
artefactual_sigs = ['SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48',
                    'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54',
                    'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60']
for sig in artefactual_sigs:
    final_table_filter.loc[sig, :] = 0

final_table_filter.to_csv('external_data/curated_match_signature_cancertype_tcgawes_literature.csv', sep='\t')
'''
# adding some signatures if needed, with associated reference.
# ACC: no ref found
# BLCA: no missing signatures
# BRCA:
sig 40 and 18 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
sig 6, 20, 26, 17, 18, 8 Landscape of somatic mutations in 560 breast cancer whole-genome sequences
# CESC Cervical squamous cell carcinoma and endocervical adenocarcinoma: ok
# CHOL: Cholangiocarcinoma
sig 6, 8, 16, 17, 20, 22 Molecular genomic landscapes of hepatobiliary cancer
# COADREAD
sig 18 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
sig 6 26 Intra-tumour diversification in colorectal cancer at the single-cell level
# DLBC Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
sig 9, 18 Genetic landscape of hepatitis B virus-associated diffuse large B-cell lymphoma
# ESCA ok, évt 40 (COSMIC)
# GBM
sig 8 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
# HNSC
sig 18 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
sig 3
# KICH Kidney Chromophobe
S13 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
# KIRC 
sig 2,3,9,13,17a,17b, Genomic features of renal cell carcinoma with venous tumor thrombus
# KIRP ok
# LGG ok
# LIHC
sig 4, 6, 12, 17, 23, Mutational signatures reveal the dynamic interplay of risk factors and cellular processes during liver tumorigenesis
sig 35, 18 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
sig 3
# LUAD ok
# LUSC
sig 18 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
#MESO
S13 (there is S2....)
# OV
sig 8, 39 Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
sig 16 Copy number signatures and mutational processes in ovarian carcinoma
# PAAD
sig 3, 17, 18, 36 portraits et http://molecularcasestudies.cshlp.org/content/early/2019/03/01/mcs.a003681
https://www.cell.com/cancer-cell/pdf/S1535-6108(17)30299-4.pdf subclonal mutation of BRCA, so no sig 3 observed!!!
# PCPG
sig 3 (use of PARP inhibitors)
# PRAD
sig 3, 18, 37, 34 portraits
sig 16, 8, 9 The Evolutionary Landscape of Localized Prostate The Evolutionary Landscape of Localized Prostate Cancers Drives Clinical Aggression
Sequencing of prostate cancers identifies new cancer genes, routes of progression and drug targets
# SARC ok
# SKCM
# STAD
sig 18 portraits
sig 3 ArticleComparative Molecular Analysis of GastrointestinalAdenocarcinomas
# TCGT
# THCA
sig 18 portraits
# UCEC
sig 3, 18 portraits
sig 6 Analysis of mutational signatures in primary and metastatic endometrial cancer reveals distinct patterns of DNA repair defects and shifts during tumor progression
# UCS ok
# UVM
sig 12, 16, 6 Comprehensive Genetic Landscape of Uveal Melanoma by Whole-Genome Sequencing




'''