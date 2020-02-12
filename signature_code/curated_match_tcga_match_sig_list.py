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

modif_dict = {'BLCA': ['SBS8', 'SBS40', 'SBS29'],
              'BRCA': ['SBS40', 'SBS18', 'SBS6', 'SBS20', 'SBS26', 'SBS17a', 'SBS17b', 'SBS8', 'SBS9', 'SBS37', 'SBS41'],
              'CHOL': ['SBS6', 'SBS8', 'SBS16', 'SBS17a', 'SBS17b', 'SBS20', 'SBS22'],
              'COADREAD': ['SBS18', 'SBS6', 'SBS26', 'SBS10a', 'SBS10b', 'SBS28', 'SBS37', 'SBS44'],
              'DLBC': ['SBS18', 'SBS9', 'SBS3', 'SBS34', 'SBS36', 'SBS37'],
              'ESCA': ['SBS28', 'SBS40'],
              'GBM': ['SBS11'],
              'HNSC': ['SBS3', 'SBS18', 'SBS7a', 'SBS7b', 'SBS7d', 'SBS16', 'SBS17a', 'SBS17b', 'SBS33', 'SBS40'],
              'KICH': ['SBS13', 'SBS17b', 'SBS29'],
              'KIRC': ['SBS2', 'SBS3', 'SBS9', 'SBS13', 'SBS17a', 'SBS17b', 'SBS22', 'SBS29', 'SBS41'],
              'LIHC': ['SBS3', 'SBS4', 'SBS6', 'SBS12', 'SBS17a', 'SBS17b', 'SBS18', 'SBS23', 'SBS35', 'SBS14', 'SBS19', 'SBS26', 'SBS28', 'SBS29', 'SBS31', 'SBS35'],
              'LUAD': ['SBS3', 'SBS9', 'SBS17a', 'SBS17b', 'SBS18', 'SBS28']
              'LUSC': ['SBS8'],
              'MESO': ['SBS13'],
              'OV': ['SBS8', 'SBS16', 'SBS39', 'SBS18', 'SBS26', 'SBS35', 'SBS41'],
              'PAAD': ['SBS3', 'SBS6', 'SBS8', 'SBS17a', 'SBS17b', 'SBS18', 'SBS20', 'SBS26', 'SBS28', 'SBS30', 'SBS40'],
              'PCPG': ['SBS3'],
              'PRAD': ['SBS3', 'SBS8', 'SBS9', 'SBS16', 'SBS18', 'SBS2', 'SBS13', 'SBS33', 'SBS37', 'SBS41'],
              'SARC': ['SBS17a', 'SBS17b'],
              'SKCM': ['SBS2', 'SBS13', 'SBS38', 'SBS40', 'SBS17a', 'SBS17b'],
              'STAD': ['SBS3', 'SBS18', 'SBS21', 'SBS26', 'SBS28', 'SBS41', 'SBS44'],
              'UCEC': ['SBS6', 'SBS18'],
              'UCS': ['SBS3', 'SBS6', 'SBS10a', 'SBS10b', 'SBS14', 'SBS15', 'SBS26', 'SBS28', 'SBS40', 'SBS44'],
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
# BLCA: 
sig 8, 29, 40 The repertoire of mutational signatures in human cancer 2020 fig3
# BRCA:
sig 6, 20, 26, 17, 18, 8 Landscape of somatic mutations in 560 breast cancer whole-genome sequences
sig 9, 37, 41, 40, 18 The repertoire of mutational signatures in human cancer 2020 fig3
# CESC Cervical squamous cell carcinoma and endocervical adenocarcinoma: ok
# CHOL: Cholangiocarcinoma
sig 6, 8, 16, 17, 20, 22 Molecular genomic landscapes of hepatobiliary cancer
# COADREAD
sig 6 26 Intra-tumour diversification in colorectal cancer at the single-cell level
sig 10a, 10b, 28, 37, 44, 18 The repertoire of mutational signatures in human cancer 2020 fig3
# DLBC Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
sig 9, 18 Genetic landscape of hepatitis B virus-associated diffuse large B-cell lymphoma
sig 3, 34, 36, 37 The repertoire of mutational signatures in human cancer 2020 fig3
# ESCA ok, évt 40 (COSMIC)
sig 28, 40 The repertoire of mutational signatures in human cancer 2020 fig3
# GBM
sig 11 The repertoire of mutational signatures in human cancer 2020 fig3
# HNSC
sig 3
sig 7a, 7b, 7d, 16, 17a, 17b, 33, 40, 18 The repertoire of mutational signatures in human cancer 2020 fig3
# KICH Kidney Chromophobe
sig 17b, 29, 13  The repertoire of mutational signatures in human cancer 2020 fig3
# KIRC 
sig 2,3,9,13,17a,17b, Genomic features of renal cell carcinoma with venous tumor thrombus
sig 22, 29, 41  The repertoire of mutational signatures in human cancer 2020 fig3
# KIRP ok
# LGG ok
# LIHC
sig 4, 6, 12, 17, 23, Mutational signatures reveal the dynamic interplay of risk factors and cellular processes during liver tumorigenesis
sig 3 cosmic v3
sig 14, 19, 26, 28, 29, 31, 35, 18 The repertoire of mutational signatures in human cancer 2020 fig3
# LUAD
sig 3, 9, 17a, 17b, 18, 28 The repertoire of mutational signatures in human cancer 2020 fig3
# LUSC
sig 8 The repertoire of mutational signatures in human cancer 2020 fig3
#MESO
S13 (there is S2....)
# OV
sig 8, 39, 18, 26, 35, 41  The repertoire of mutational signatures in human cancer 2020 fig3
sig 16 Copy number signatures and mutational processes in ovarian carcinoma
# PAAD
https://www.cell.com/cancer-cell/pdf/S1535-6108(17)30299-4.pdf subclonal mutation of BRCA, so no sig 3 observed!!!
sig 3, 6, 8 17, 18, 20, 26, 28, 30, 40 The repertoire of mutational signatures in human cancer 2020 fig3
# PCPG
sig 3 (use of PARP inhibitors)
# PRAD
sig 16, 8, 9 The Evolutionary Landscape of Localized Prostate The Evolutionary Landscape of Localized Prostate Cancers Drives Clinical Aggression
Sequencing of prostate cancers identifies new cancer genes, routes of progression and drug targets
sig 2, 3, 13, 18, 33, 37, 41 The repertoire of mutational signatures in human cancer 2020 fig3
# SARC
sig 17a and b The repertoire of mutational signatures in human cancer 2020 fig3 (bone osteosarcoma and soft tissue leiomyosarcoma)
# SKCM
sig 2, 13, 17a, 17b, 38, 40 The repertoire of mutational signatures in human cancer 2020 fig3
# STAD
sig 18, 21, 26, 28, 41, 44  The repertoire of mutational signatures in human cancer 2020 fig3
sig 3 ArticleComparative Molecular Analysis of GastrointestinalAdenocarcinomas
# TCGT
# UCEC
sig 18  The repertoire of mutational signatures in human cancer 2020 fig3
sig 6 Analysis of mutational signatures in primary and metastatic endometrial cancer reveals distinct patterns of DNA repair defects and shifts during tumor progression
# UCS
sig 3, 6, 10a, 10b, 14, 15, 26, 28, 40, 44  The repertoire of mutational signatures in human cancer 2020 fig3
# UVM
sig 12, 16, 6 Comprehensive Genetic Landscape of Uveal Melanoma by Whole-Genome Sequencing




'''