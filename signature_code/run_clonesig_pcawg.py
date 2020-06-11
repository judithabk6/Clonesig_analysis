#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import sys
import pkg_resources
from clonesig import mixin_init_parameters
from clonesig.data_loader import PAT_LIST
import time
import signal
from clonesig.run_clonesig import run_clonesig
import scipy as sp
from clonesig.estimator import EV_DOF_THRESHOLD
from clonesig.evaluate import score_sig_1A_base

patient_id = sys.argv[1]
cancer_loc = sys.argv[2]
cohort = sys.argv[3]

"""
patient_id = '51800588-c622-11e3-bf01-24c6515278c0'
cancer_loc = 'Liver-HCC'
"""

def get_context(x):
    match_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if x['Tumor_Seq_Allele1'] in ('C', 'T'):
        context = x['ref_context'][9].upper() + '[' + x['Tumor_Seq_Allele1'] + '>' + x['Tumor_Seq_Allele2'] + ']' + \
            x['ref_context'][11].upper()
    else:
        context = match_dict[x['ref_context'][11].upper()] + '[' + match_dict[x['Tumor_Seq_Allele1']] +\
            '>' + match_dict[x['Tumor_Seq_Allele2']] + ']' + match_dict[x['ref_context'][9].upper()]
    return context

class TimeoutException(Exception):
    pass

def signal_handler(signum, frame):
    raise TimeoutException("Timed out!")

pcawg_folder = '20200510_pcawg'
cols = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
        'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele',
        'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
        'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode', 'Genome_Change', 'ref_context',
        'gc_content', 'i_1000genomes_AF', 'i_1000genomes_ID', 'i_Callers',
        'i_GERM1000G', 'i_GERMOVLP', 'i_LOWSUPPORT', 'i_NORMALPANEL',
        'i_NumCallers', 'i_OXOGFAIL', 'i_REMAPFAIL', 'i_SEXF', 'i_VAF',
        'i_bPcr', 'i_bSeq', 'i_qual', 'i_repeat_masker', 'i_signature_N3',
        'i_signature_R1', 'i_signature_R2', 'i_snv_near_indel', 't_alt_count',
        't_ref_count', 'i_model_score', 'i_n_vaf', 'Project_Code', 'Donor_ID']

filename = '{}/{}.maf'.format(pcawg_folder, patient_id)
maf_data = pd.read_csv(filename, sep='\t', header=None, names=cols)
cnv_folder = 'pcawg_download/cnv/consensus.20170119.somatic.cna.{}.public'.format(cohort)
## to modify and to be an argument to separate icgc and tcga
cnv_filename = '{}/{}.consensus.20170119.somatic.cna.txt'.format(cnv_folder, patient_id)
cnv_table = pd.read_csv(cnv_filename, sep='\t')
purity_filename = 'pcawg_download/cnv/consensus.20170217.purity.ploidy.txt'
purity_table = pd.read_csv(purity_filename, sep='\t')
purity = purity_table[purity_table.samplename==patient_id].purity.mean()
ploidy = purity_table[purity_table.samplename==patient_id].ploidy.mean()

signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
if cohort == 'icgc':
    match_cancertype = pd.read_csv(
      'external_data/curated_match_signature_cancertype_pcawgs_literature.csv',
      index_col=0, sep='\t')
    match_cancertype.index = match_cancertype['Unnamed: 0.1']
else:
    match_cancertype = pd.read_csv(
      'external_data/curated_match_signature_cancertype_tcgawes_literature.csv',
      index_col=0, sep='\t')


selected_sigs = match_cancertype[[cancer_loc]][match_cancertype[[cancer_loc]] > 0].dropna().index.tolist()
sub_sig = sig[[c for c in sig.columns if c in selected_sigs]]
prefit = False
sigprofiler = False

sig_matrix = sub_sig.values.astype(float).T
new_sig_matrix = sig_matrix + mixin_init_parameters.ZERO_PADDING \
    * (sig_matrix == 0)
MU = new_sig_matrix / new_sig_matrix.sum(axis=1)[:, np.newaxis]

final_cols_used = [c for c in sig.columns if c in selected_sigs]

model_selection_kws = {'factor':  0.065}

method = 'clonesig'
id_list = list()
metrics_list = list()
clonal_sigs = pd.DataFrame(index=[0], columns=sig.columns[2:])
subclonal_sigs = pd.DataFrame(index=[0], columns=sig.columns[2:])

maf_data_rel = maf_data[maf_data.Variant_Type=='SNP']
maf_data_rel = maf_data_rel.assign(
        pattern=pd.Categorical(
            maf_data_rel.apply(get_context, axis=1),
            categories=PAT_LIST, ordered=True))
maf_data_rel = maf_data_rel.assign(
        trinucleotide=maf_data_rel.apply(lambda x: PAT_LIST.index(x['pattern']),
                                         axis=1))
maf_data_rel = maf_data_rel.assign(
    t_depth=maf_data_rel.t_alt_count + maf_data_rel.t_ref_count)
maf_data_rel_cnv = pd.merge(maf_data_rel, cnv_table, left_on='Chromosome', right_on='chromosome')
maf_data_rel_cnv_filter = maf_data_rel_cnv[
    (maf_data_rel_cnv.start<=maf_data_rel_cnv.Start_position)&
    (maf_data_rel_cnv.end>=maf_data_rel_cnv.Start_position)&
    (maf_data_rel_cnv.major_cn>0.0)]

input_table = maf_data_rel_cnv_filter[['trinucleotide', 't_depth',
                                       't_alt_count', 'major_cn',
                                       'minor_cn', 'total_cn']]
input_table_nona = input_table.dropna(how='any', axis=0)


start = time.time()
signal.signal(signal.SIGALRM, signal_handler)
signal.alarm(48*3600)   # 48 hours
try:
    new_est, lr, pval, new_inputMU, cst_est, fitted_sigs = run_clonesig(
        input_table_nona.trinucleotide,
        input_table_nona.t_alt_count.values,
        input_table_nona.t_depth.values, np.array([2]*len(input_table_nona)),
        input_table_nona.total_cn.values,
        input_table_nona.minor_cn.values, purity, MU,
        return_sig_change_test=True, min_mut_clone=5, min_prop_sig=0.0,
        prefit_signatures=prefit, prefit_thresh=0.01,
        model_selection_kws=model_selection_kws, max_nb_clones=10)
except TimeoutException:
    print('too long')
end = time.time()
ev, _ = np.linalg.eig(
    1-sp.spatial.distance.squareform(
        sp.spatial.distance.pdist(new_inputMU, 'cosine')))
dof = sum(ev > EV_DOF_THRESHOLD)

prop_diploid = len(input_table_nona[(input_table_nona.minor_cn == 1) &
                                    (input_table_nona.major_cn == 1)]) / \
    len(input_table_nona)

id_list.append([patient_id, cancer_loc,
                prop_diploid, ploidy, purity, len(input_table_nona),
                MU.shape[0], new_inputMU.shape[0],
                input_table_nona.major_cn.mean(),
                input_table_nona.total_cn.mean(), method, str(prefit),
                str(sigprofiler), dof])
clonal_idx = np.argmax(new_est.phi)
clonal_phi = new_est.phi[clonal_idx]
clonal_xi = new_est.xi[clonal_idx]

if fitted_sigs is not None:
    clonal_sigs.loc[0, np.array(final_cols_used)[fitted_sigs]] = new_est.pi[clonal_idx, :]
    sq_dist = sp.spatial.distance.squareform(
        sp.spatial.distance.pdist(new_est.pi.dot(MU[fitted_sigs, :]),
                                  'cosine'))
else:
    clonal_sigs.iloc[0, :][final_cols_used] = new_est.pi[clonal_idx, :]
    sq_dist = sp.spatial.distance.squareform(
        sp.spatial.distance.pdist(new_est.pi.dot(MU), 'cosine'))
nb_mut_clonal = sum(new_est.qun.argmax(axis=1) == clonal_idx)
if new_est.J >= 2:
    sub_phi = new_est.phi.copy()
    sub_xi = new_est.xi.copy()
    sub_phi[clonal_idx] = sub_xi[clonal_idx] = 0
    largest_subclonal_idx = np.argmax(sub_xi)
    largest_subclonal_phi = new_est.phi[largest_subclonal_idx]
    largest_subclonal_xi = new_est.xi[largest_subclonal_idx]
    if fitted_sigs is not None:
        subclonal_sigs.loc[0, np.array(final_cols_used)[fitted_sigs]] = \
            new_est.pi[largest_subclonal_idx, :]
    else:
        subclonal_sigs.iloc[0, :][final_cols_used] = \
            new_est.pi[largest_subclonal_idx, :]
    nb_mut_largest_subclonal = \
        sum(new_est.qun.argmax(axis=1) == largest_subclonal_idx)
    clonal_largest_sub_pidist = sq_dist[clonal_idx, largest_subclonal_idx]
else:
    (largest_subclonal_phi, largest_subclonal_xi, nb_mut_largest_subclonal,
     clonal_largest_sub_pidist) = [np.nan] * 4
largest_pi_dist = np.max(sq_dist)

raw_data_distrib = np.zeros(96)
val, c = np.unique(input_table_nona.trinucleotide,
                   return_counts=True)
raw_data_distrib[val.astype(int)] = c
raw_data_distrib = raw_data_distrib / len(input_table_nona)
est_distrib = new_est.xi.dot(new_est.pi).dot(new_est.mu_matrix)
overall_profile_dist = score_sig_1A_base(raw_data_distrib, est_distrib)

metrics_list.append([new_est.J, lr, pval, clonal_phi, clonal_xi,
                     nb_mut_clonal, largest_subclonal_phi,
                     largest_subclonal_xi, nb_mut_largest_subclonal,
                     clonal_largest_sub_pidist, largest_pi_dist,
                     overall_profile_dist, end - start])


id_cols = ['patient_id', 'mutation_set', 
           'prop_diploid', 'ploidy', 'purity', 'nb_mut', 'nb_sigs',
           'nb_sigs_prefit', 'major_cn_mean', 'total_cn_mean', 'method',
           'prefit_bool', 'sigprofiler_bool', 'dof']
metrics_cols = ['nb_clones', 'lr', 'pval', 'clonal_phi', 'clonal_xi',
                'nb_mut_clonal', 'largest_subclonal_phi',
                'largest_subclonal_xi', 'nb_mut_largest_subclonal',
                'clonal_largest_sub_pidist', 'largest_pi_dist',
                'overall_profile_dist', 'runtime']

id_df = pd.DataFrame(id_list, columns=id_cols)
metrics_df = pd.DataFrame(metrics_list, columns=metrics_cols)

res_df = pd.concat([id_df, metrics_df], axis=1)
res_df.to_csv('{}/{}_clonesig_results_restr_curated.csv'
              .format(pcawg_folder, patient_id),
              index=False, sep='\t')
clonal_sigs.to_csv('{}/{}_clonesig_clonal_sigs_restr_curated.csv'
                   .format(pcawg_folder, patient_id),
                   index=False, sep='\t')
subclonal_sigs.to_csv('{}/{}_clonesig_subclonal_sigs_restr_curated.csv'
                      .format(pcawg_folder, patient_id),
                      index=False, sep='\t')


