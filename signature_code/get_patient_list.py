#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to conduct and exploratory analysis of the new tcga data.
"""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from util_functions import safe_mkdir
import matplotlib.pyplot as plt

import seaborn as sns
import requests


pd.set_option('display.max_columns', 200)

plt.style.use('seaborn-white')

cancer_loc = sys.argv[1]


# load ascat purity
all_ascat_purity = pd.read_csv('external_data/liftover_ASCAT_TCGA/filtered.combined.acf.ploidy.txt', sep='\t')
all_ascat_purity_relevant = all_ascat_purity[all_ascat_purity.cancer_type==cancer_loc]
all_ascat_purity_samples = [s[:16] for s in all_ascat_purity_relevant.barcodeTumour.unique()]
all_ascat_purity_patients = [s[:12] for s in all_ascat_purity_samples]


# load CNV data
all_cnv_hg38 = pd.read_csv('external_data/liftover_ASCAT_TCGA/cnasHg38.tsv',
                           sep='\t')
all_cnv_hg38 = all_cnv_hg38.assign(patient_id=all_cnv_hg38.participant.str[:12])
all_cnv_hg38 = all_cnv_hg38.assign(weight=all_cnv_hg38.stop-all_cnv_hg38.start)
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.type==cancer_loc]
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.participant.isin(all_ascat_purity['name'])]
all_cnv_hg38 = all_cnv_hg38.assign(key=all_cnv_hg38.patient_id+'_'+all_cnv_hg38.chromosome)

patient_list_list = list()
for file_type in ['qc', 'public']:
    # load snv data
    filename = 'tmp/useful_final_{}_{}.csv'.format(file_type, cancer_loc)
    useful_final = pd.read_csv(filename, sep='\t', index_col=0)
    useful_final = useful_final\
        .assign(patient_id=useful_final.sample_id.str[:12])
    useful_final = useful_final\
        .assign(cancer_loc=cancer_loc)
    try:
        annotations = pd.read_csv('data_tcga_wes/{}/annotations.txt'.format(cancer_loc),
                                  sep='\t')
        barcode_list = list()
        for ann_uuid in annotations['id'].tolist():
            r = requests.get('https://api.gdc.cancer.gov/annotations/{uuid}'.format(uuid=ann_uuid))
            barcode_list.append(r.json()['data']['entity_submitter_id'])
        annotations = annotations.assign(**{'Entity Barcode': barcode_list})
        exclude = [b[:12] for b in annotations['Entity Barcode'].unique().tolist()]
    except:
        exclude = []
    useful_final = useful_final[~useful_final.patient_id.isin(exclude)]

    useful_final = useful_final\
        .assign(sample_id_short=useful_final.sample_id.str[:-12])

    useful_final = useful_final.assign(key=useful_final.patient_id+'_'+useful_final['mutation_id.1'].str.split('_').str[1].str[3:])
    useful_final = useful_final.assign(position=useful_final['mutation_id.1'].str.split('_').str[2].astype(int))

    tmp_useful_final_cnv = pd.merge(useful_final, all_cnv_hg38, left_on='key', right_on='key', how='left')
    tmp_useful_final_cnv = tmp_useful_final_cnv[((tmp_useful_final_cnv.position >= tmp_useful_final_cnv.start) & (tmp_useful_final_cnv.position <= tmp_useful_final_cnv.stop))]
    useful_final_cnv = pd.merge(useful_final,
                                tmp_useful_final_cnv[['mutation_id.1', 'major', 'minor']],
                                left_on='mutation_id.1', right_on='mutation_id.1', how='left')
    useful_final_cnv = useful_final_cnv[(useful_final_cnv.sample_id_short.isin(all_ascat_purity_samples)) |
                                        (~useful_final_cnv.patient_id.isin(all_ascat_purity_patients))]
    useful_final_cnv_purity = pd.merge(useful_final_cnv,
                                       all_ascat_purity[['tissue', 'purity', 'ploidy']],
                                       left_on='patient_id', right_on='tissue', how='left')
    useful_final_cnv_purity = useful_final_cnv_purity.assign(total_cn=lambda x: x['minor'] + x['major'])
    useful_final_cnv_purity = useful_final_cnv_purity[useful_final_cnv_purity.total_cn > 0]

    patient_list_list.append(useful_final_cnv_purity.patient_id.unique().tolist())
    print(useful_final_cnv_purity.shape, useful_final_cnv_purity['mutation_id.1'].nunique())
    useful_final_cnv_purity.drop_duplicates(subset=['mutation_id.1'], inplace=True)
    useful_final_cnv_purity.to_csv('tmp/{}_useful_final_{}_merge_cnv_purity.csv'.format(cancer_loc, file_type), sep='\t', index=False)

intersection_patients = list(set(patient_list_list[0]).intersection(set(patient_list_list[1])))
udf = pd.DataFrame(intersection_patients, columns=['patient_id'], index=range(1, len(intersection_patients)+1))
udf = udf.assign(cancer_loc=cancer_loc)
udf.to_csv('{}_official_patient_list.csv'.format(cancer_loc), header=False)
