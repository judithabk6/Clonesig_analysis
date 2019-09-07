#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import os
import sys
from clonesig.run_clonesig import run_clonesig
import pkg_resources
from clonesig import mixin_init_parameters
from clonesig.estimator import EV_DOF_THRESHOLD
import time
from clonesig.data_loader import get_context, PAT_LIST
import scipy as sp
import pickle
from clonesig.evaluate import score_sig_1A_base

try:
    rows, columns = os.popen('stty size', 'r').read().split()
    pd.set_option('display.width', int(columns))
    pd.set_option('display.max_columns', 600)
    pd.set_option('display.max_rows', 1000)
except:
    print("running on server, otherwise please investigate")

patient_id = sys.argv[1]
cancer_loc = sys.argv[2]
try:
    cosmic_type = sys.argv[3]
except IndexError:
    cosmic_type = ''

"""
patient_id = 'TCGA-OR-A5L4'
cancer_loc = 'ACC'
cosmic_type = ''

patient_id = 'TCGA-QG-A5YW'
cancer_loc = 'COADREAD'
cosmic_type = 'ColoRect-AdenoCA'
"""

tcga_folder = '20190704_TCGA'

cols = {'public': ['nb_caller', 'db_snp_status', 'exac_freq', '1000G_freq',
                   'COSMIC', 'PHENO', 'IMPACT', 'SIFT', 'PolyPhen',
                   'Hugo_Symbol', 'Gene', 'SYMBOL', 'SOMATIC', 'GDC_FILTER',
                   'CONTEXT', 't_depth', 't_ref_count', 't_alt_count',
                   'n_depth', 'n_ref_count', 'n_alt_count', 'mutation_id.1',
                   'sample_id', 'patient_id', 'cancer_loc', 'sample_id_short',
                   'key', 'position', 'major', 'minor', 'tissue', 'purity',
                   'ploidy', 'total_cn'],
        'protected': ['nb_caller', 'db_snp_status', 'exac_freq', '1000G_freq',
                      'COSMIC', 'PHENO', 'IMPACT', 'SIFT', 'PolyPhen',
                      'Hugo_Symbol', 'Gene', 'SYMBOL', 'SOMATIC', 'GDC_FILTER',
                      'CONTEXT', 't_depth', 't_ref_count', 't_alt_count',
                      'n_depth', 'n_ref_count', 'n_alt_count', 't_vaf',
                      'n_vaf', 'mutation_id.1', 'sample_id', 'patient_id',
                      'cancer_loc', 'sample_id_short', 'key', 'position',
                      'major', 'minor', 'tissue', 'purity', 'ploidy',
                      'total_cn']}
active_sigs_cols = ['Cancer Types', 'Sample Names', 'Accuracy', 'SBS1', 'SBS2',
                    'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c',
                    'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11',
                    'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16', 'SBS17a',
                    'SBS17b', 'SBS18', 'SBS19', 'SBS20', 'SBS21', 'SBS22',
                    'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28',
                    'SBS29', 'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34',
                    'SBS35', 'SBS36', 'SBS37', 'SBS38', 'SBS39', 'SBS40',
                    'SBS41', 'SBS42', 'SBS43', 'SBS44', 'SBS45', 'SBS46',
                    'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52',
                    'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58',
                    'SBS59', 'SBS60']

filenames = {'public': "{}/{}/public_snv.tsv".format(tcga_folder, patient_id),
             'protected': "{}/{}/protected_snv.tsv".format(tcga_folder,
                                                           patient_id)}
signature_filename = 'data/sigProfiler_exome_SBS_signatures.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
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

if prefit:
    model_selection_kws = {'factor': 0.048}
else:
    model_selection_kws = {'factor':  0.093}

method = 'clonesig'
id_list = list()
metrics_list = list()
clonal_sigs = pd.DataFrame(index=[0, 1], columns=sig.columns[2:])
subclonal_sigs = pd.DataFrame(index=[0, 1], columns=sig.columns[2:])
for i, folder_type in enumerate(['public', 'protected']):
    snv_table = pd.read_csv(filenames[folder_type], sep='\t', header=None,
                            names=cols[folder_type])

    snv_table = snv_table.assign(
            mutation_id=snv_table['mutation_id.1'].str[29:-2] +
            snv_table['mutation_id.1'].str[-1])
    snv_table = snv_table.assign(
            ref=snv_table.CONTEXT.str[5])
    snv_table = snv_table.assign(
            alt=snv_table.mutation_id.str[-1])
    snv_table = snv_table.assign(
            pattern=pd.Categorical(
                snv_table.apply(get_context, axis=1),
                categories=PAT_LIST, ordered=True))
    snv_table = snv_table.assign(
            trinucleotide=snv_table.apply(lambda x: PAT_LIST.index(x['pattern']), axis=1))
    snv_table_small = snv_table[['trinucleotide', 't_depth', 't_alt_count', 'major',
                                 'minor', 'total_cn']]
    snv_table_nona = snv_table_small.dropna(how='any', axis=0)

    purity = snv_table.purity.values[0]
    start = time.time()
    new_est, lr, pval, new_inputMU, cst_est, fitted_sigs = run_clonesig(
        snv_table_nona.trinucleotide,
        snv_table_nona.t_alt_count.values,
        snv_table_nona.t_depth.values, np.array([2]*len(snv_table_nona)),
        snv_table_nona.total_cn.values,
        snv_table_nona.minor.values, purity, MU,
        return_sig_change_test=True, min_mut_clone=5, min_prop_sig=0.0,
        prefit_signatures=prefit, prefit_thresh=0.01,
        model_selection_kws=model_selection_kws, max_nb_clones=10)
    end = time.time()
    ev, _ = np.linalg.eig(
        1-sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(new_inputMU, 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)

    prop_diploid = len(snv_table[(snv_table.minor == 1) &
                                 (snv_table.total_cn == 1)]) / len(snv_table)
    ploidy = snv_table.ploidy.values[0]

    id_list.append([patient_id, folder_type, cancer_loc, cosmic_type,
                    prop_diploid, ploidy, purity, len(snv_table_nona),
                    MU.shape[0], new_inputMU.shape[0],
                    snv_table_nona.major.mean(),
                    snv_table_nona.total_cn.mean(), method, str(prefit),
                    str(sigprofiler), dof])
    clonal_idx = np.argmax(new_est.phi)
    clonal_phi = new_est.phi[clonal_idx]
    clonal_xi = new_est.xi[clonal_idx]
    if fitted_sigs is not None:
        clonal_sigs.loc[i, np.array(final_cols_used)[fitted_sigs]] = new_est.pi[clonal_idx, :]
        sq_dist = sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(new_est.pi.dot(MU[fitted_sigs, :]),
                                      'cosine'))
    else:
        clonal_sigs.iloc[i, :][final_cols_used] = new_est.pi[clonal_idx, :]
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
            subclonal_sigs.loc[i, np.array(final_cols_used)[fitted_sigs]] = \
                new_est.pi[largest_subclonal_idx, :]
        else:
            subclonal_sigs.iloc[i, :][final_cols_used] = \
                new_est.pi[largest_subclonal_idx, :]
        nb_mut_largest_subclonal = \
            sum(new_est.qun.argmax(axis=1) == largest_subclonal_idx)
        clonal_largest_sub_pidist = sq_dist[clonal_idx, largest_subclonal_idx]
    else:
        (largest_subclonal_phi, largest_subclonal_xi, nb_mut_largest_subclonal,
         clonal_largest_sub_pidist) = [np.nan] * 4
    largest_pi_dist = np.max(sq_dist)

    raw_data_distrib = np.zeros(96)
    val, c = np.unique(snv_table_nona.trinucleotide,
                       return_counts=True)
    raw_data_distrib[val.astype(int)] = c
    raw_data_distrib = raw_data_distrib / len(snv_table_nona)
    est_distrib = new_est.xi.dot(new_est.pi).dot(new_est.mu_matrix)
    overall_profile_dist = score_sig_1A_base(raw_data_distrib, est_distrib)

    metrics_list.append([new_est.J, lr, pval, clonal_phi, clonal_xi,
                         nb_mut_clonal, largest_subclonal_phi,
                         largest_subclonal_xi, nb_mut_largest_subclonal,
                         clonal_largest_sub_pidist, largest_pi_dist,
                         overall_profile_dist, end - start])

    with open('{}/{}/{}_clonesig_raw_results_restr_curated'
              .format(tcga_folder, patient_id, folder_type), 'wb') as raw_res:
        my_pickler = pickle.Pickler(raw_res)
        my_pickler.dump([new_est, lr, pval, cst_est, fitted_sigs])

id_cols = ['patient_id', 'mutation_set', 'cancer_loc', 'cosmic_type',
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
res_df.to_csv('{}/{}/clonesig_results_restr_curated.csv'
              .format(tcga_folder, patient_id),
              index=False, sep='\t')
clonal_sigs.to_csv('{}/{}/clonesig_clonal_sigs_restr_curated.csv'
                   .format(tcga_folder, patient_id),
                   index=False, sep='\t')
subclonal_sigs.to_csv('{}/{}/clonesig_subclonal_sigs_restr_curated.csv'
                      .format(tcga_folder, patient_id),
                      index=False, sep='\t')








