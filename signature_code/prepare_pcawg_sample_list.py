#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd

metadata = pd.read_csv('pcawg_download/WGS.metadata.tsv', sep='\t')
histo_data = pd.read_excel('pcawg_download/pcawg_specimen_histology_August2016_v9.xlsx')
match_cohort = pd.read_excel('pcawg_download/tumour_subtype_consolidation_map.tsv.xlsx', sheet_name=1)

match_cohort = match_cohort.dropna(subset=['Tier 2 (organ system)'])

match_cohort2 = match_cohort.join(
    match_cohort.pop('Contributing projects').str.split(',',expand=True))
cols = match_cohort.columns.to_list()
match_cohort_ext = pd.melt(match_cohort2, id_vars=cols, value_name='cohort').\
    dropna(subset=['cohort'])


TCGA_cohorts = match_cohort_ext[match_cohort_ext.cohort.str.contains('US')]
manual_additions = {'Abbreviation': ['Eso-AdenoCA', 'Panc-AdenoCA'],
                    'cohort': ['ESCA-US', 'PAAD-US']}
TCGA_cohorts_complete = pd.concat(
    [TCGA_cohorts,
     pd.DataFrame(manual_additions, columns=TCGA_cohorts.columns)])
missing_cohorts = set(match_cohort_ext.Abbreviation).\
    difference(TCGA_cohorts_complete.Abbreviation)
'''
{'Biliary-AdenoCA',
 'Bone-Benign',
 'Bone-Epith',
 'Bone-Osteosarc',
 'CNS-Medullo',
 'CNS-PiloAstro',
 'Lymph-CLL',
 'Myeloid-MPN',
 'Panc-Endocrine'}
'''

TCGA_signature_match = pd.read_csv('external_data/curated_match_signature_cancertype_tcgawes_literature.csv', sep='\t')
# From figure 3 of Alexandrov et al, Nature, 2020
modif_dict = {'Biliary-AdenoCA': ['SBS1', 'SBS2', 'SBS3', 'SBS5',
                                  'SBS9', 'SBS12', 'SBS13', 'SBS15', 'SBS17a',
                                  'SBS17b', 'SBS18', 'SBS21', 'SBS22', 'SBS24',
                                  'SBS32', 'SBS40', 'SBS44'],
              'Bone-Other': ['SBS1', 'SBS2', 'SBS5', 'SBS8', 'SBS13', 'SBS17a',
                             'SBS17b', 'SBS40'],
              'Bone-Osteosarc': ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8',
                                 'SBS13', 'SBS17a', 'SBS17b', 'SBS30', 'SBS40'],
              'CNS-Medullo': ['SBS1', 'SBS5', 'SBS8', 'SBS18', 'SBS39',
                              'SBS40'],
              'CNS-PiloAstro': ['SBS1', 'SBS5', 'SBS19', 'SBS23', 'SBS40'],
              'Lymph-CLL': ['SBS1', 'SBS5', 'SBS9', 'SBS40'],
              'Myeloid-AML': ['SBS1', 'SBS5', 'SBS18', 'SBS31'],
              'Myeloid-MDS/MPN': ['SBS1', 'SBS2', 'SBS5', 'SBS19', 'SBS23',
                              'SBS32'],
              'Panc-Endocrine': ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS6', 'SBS8',
                                 'SBS9', 'SBS11', 'SBS13', 'SBS26', 'SBS30',
                                 'SBS36', 'SBS39']}
TCGA_cohorts_complete = TCGA_cohorts_complete.assign(
    tcga_name=TCGA_cohorts_complete.cohort.str.split('-').str[0])
TCGA_cohorts_complete.loc[TCGA_cohorts_complete.tcga_name.isin(
    ['COAD', 'READ']), 'tcga_name'] = 'COADREAD'
TCGA_cohorts_complete = TCGA_cohorts_complete[TCGA_cohorts_complete.tcga_name!='LAML']
TCGA_cohorts_complete.loc[TCGA_cohorts_complete.Abbreviation.str.contains('Breast'), 'Abbreviation'] = 'Breast'
mini_match = TCGA_cohorts_complete.drop_duplicates(subset=['tcga_name'])
ok_names = mini_match.tcga_name.to_list()
pcawg_signature_match = TCGA_signature_match[['Unnamed: 0'] + ok_names]
pcawg_signature_match.columns = ['Unnamed: 0'] + mini_match.Abbreviation.to_list()
new_df = pd.DataFrame(index=TCGA_signature_match['Unnamed: 0'].values,
                      columns=modif_dict.keys())
for cancer_loc, sig_list in modif_dict.items():
    for sig in sig_list:
        new_df.loc[sig, cancer_loc] = 1
new_df = new_df.fillna(0).reset_index(drop=True)
final_pcawg_signature_match = pd.concat([pcawg_signature_match, new_df],
                                        axis=1)
final_pcawg_signature_match.to_csv('external_data/curated_match_signature_cancertype_pcawgs_literature.csv', sep='\t')

metadata_histo = pd.merge(metadata, histo_data[['icgc_sample_id', 'histology_abbreviation']], on='icgc_sample_id', how='left')
#!ls -lht pcawg_download/cnv/consensus.20170119.somatic.cna.icgc.public/ > icgc_cnv_file_list.tsv
icgc_plist = pd.read_csv('icgc_cnv_file_list.tsv', sep='\s+', skiprows=1,
                         header=None, names=['permissions', 'number', 'owner',
                                            'group', 'size', 'month', 'day',
                                            'year', 'filename'],
                         index_col=False)
file_list = icgc_plist.filename.to_list()
sample_list = [str(f).split('.')[0] for f in file_list]
metadata_histo_small = metadata_histo.drop(['read_group_info', 'read_length_r1', 'read_length_r2'], axis=1).drop_duplicates()
metadata_histo_small_relevant = metadata_histo_small[metadata_histo_small.aliquot_id.isin(sample_list)]
renaming_histo_dict = {'Bone-Benign': 'Bone-Other',
                       'Bone-Epith': 'Bone-Other',
                       'Breast-AdenoCA': 'Breast',
                       'Breast-DCIS': 'Breast',
                       'Breast-LobularCA': 'Breast',
                       'Myeloid-MDS': 'Myeloid-MDS/MPN',
                       'Myeloid-MPN': 'Myeloid-MDS/MPN'}
for k, v in renaming_histo_dict.items():
    metadata_histo_small_relevant.loc[metadata_histo_small_relevant.histology_abbreviation==k, 'histology_abbreviation'] = v
metadata_histo_small_relevant.index = range(1, len(metadata_histo_small_relevant)+1)
metadata_histo_small_relevant = metadata_histo_small_relevant.assign(cohort='icgc')
metadata_histo_small_relevant[['aliquot_id', 'histology_abbreviation', 'cohort']].to_csv('20200510_pcawgs_final_patient_list.csv', sep=',', header=False, index=True)

#!ls -lht pcawg_download/cnv/consensus.20170119.somatic.cna.tcga.public/ > tcga_pcawg_cnv_file_list.tsv
tcga_plist = pd.read_csv('tcga_pcawg_cnv_file_list.tsv', sep='\s+', skiprows=1,
                         header=None, names=['permissions', 'number', 'owner',
                                            'group', 'size', 'month', 'day',
                                            'year', 'filename'],
                         index_col=False)
file_list = tcga_plist.filename.to_list()
sample_list = [str(f).split('.')[0] for f in file_list]
metadata_histo_small = metadata_histo.drop(['read_group_info', 'read_length_r1', 'read_length_r2'], axis=1).drop_duplicates()
metadata_histo_small_relevant = metadata_histo_small[metadata_histo_small.aliquot_id.isin(sample_list)]
metadata_histo_small_relevant = metadata_histo_small_relevant.assign(cancer_loc=metadata_histo_small_relevant.dcc_project_code.str.split('-').str[0])
metadata_histo_small_relevant.loc[metadata_histo_small_relevant.cancer_loc.isin(['COAD', 'READ']), 'cancer_loc'] = 'COADREAD'
metadata_histo_small_relevant = metadata_histo_small_relevant.assign(cohort='tcga')
metadata_histo_small_relevant.index = range(1, len(metadata_histo_small_relevant)+1)
metadata_histo_small_relevant[['aliquot_id', 'cancer_loc', 'cohort']].to_csv('20200510_pcawgs_tcga_final_patient_list.csv', sep=',', header=False, index=True)

