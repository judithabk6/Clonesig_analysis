#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pkg_resources
import bz2
import pickle
import os
from clonesig.estimator import Estimator
from pandas.errors import EmptyDataError
from clonesig.evaluate import score_sig_1D_base



MIXTURE_THRESHOLD = 0.05


"""
folder_path = 'PhylogicNDT500/Sim_500_19_cst'
"""
folder_path = sys.argv[1]

signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
all_sigs = sig.columns[2:].to_list()

input_filename = '{}/input_t.tsv'.format(folder_path)
input_df = pd.read_csv(input_filename, sep='\t')
input_df = input_df.assign(signature=input_df.signature.astype(int))
true_signatures_1D_df = input_df[['mutation_id', 'signature']]

# tracksig
try:
    mixture_file = pd.read_csv('{}/tracksig/tracksig_mixtures_cancertype.csv'.
                               format(folder_path), sep=',')
except FileNotFoundError:
    print('non working')
try:
    changepoint_file = pd.read_csv(
        '{}/tracksig/tracksig_changepoints_cancertype.txt'.
        format(folder_path), header=None, sep=' ')
    changepoints_tracksig_list = changepoint_file.values[0]
except EmptyDataError:
    changepoints_tracksig_list = np.array(list())

input_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
with open('{}/purity.txt'.format(folder_path), 'r') as f:
    purity = float(f.read())
input_df = input_df.assign(mut_cn=1)
input_df = input_df.assign(vaf=input_df.var_counts /
                         (input_df.ref_counts + input_df.var_counts))
input_df = input_df.assign(
    total_cn=lambda x: x['minor_cn'] + x['major_cn'])
input_df = input_df.assign(
    vaf_cn=input_df.vaf * input_df['total_cn'] / input_df['mut_cn'])
input_df = input_df.assign(
    vaf_purity=input_df.apply(
        lambda x: x['vaf']/purity *
        ((1 - purity) * 2 + purity * x['total_cn']) /
        x['mut_cn'], axis=1))
input_df.sort_values(by='vaf_purity', inplace=True)
input_df.reset_index(inplace=True, drop=True)

input_df = input_df.assign(mutation_group=lambda x: x.index//100)
nbin = len(input_df)//100
input_df_filter = input_df[input_df.mutation_group <= nbin - 1]
cluster_id_list = np.zeros(input_df_filter.mutation_group.nunique())
i = 1
for chg_point in changepoints_tracksig_list:
    cluster_id_list[(chg_point - 1):] = i
    i += 1
input_df_filter = input_df_filter.assign(
    pred_cluster_id=input_df_filter.apply(
        lambda x: int(cluster_id_list[x['mutation_group']]), axis=1))

pred_signatures = np.zeros(len(all_sigs))
filename = '{}/subMU.csv'.format(folder_path)
sub_matrix = pd.read_csv(filename, sep='\t')
mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
sub_sigs = sub_matrix.columns[1:]
idx = [list(all_sigs).index(s) for s in sub_sigs]

est_sigs = mixture_file[mixture_file.columns[1:]].mean(axis=1).values
pred_signatures[idx] = est_sigs

pred_profile = est_sigs.dot(mu_mat_setting)

mut_sig = np.moveaxis(
    np.repeat([mixture_file.values[:, 1:].T[input_df_filter.mutation_group.astype(int)]],
              96, axis=0), [0, 1, 2], [2, 0, 1])
big_mu = np.repeat([mu_mat_setting], len(input_df_filter), axis=0)
big_everything = (mut_sig * big_mu / len(est_sigs))[np.arange(len(input_df_filter)), :, input_df_filter.trinucleotide]
signature_mut = np.argmax(big_everything, axis=1)
input_df_filter = input_df_filter.assign(signature=signature_mut)
tracksig_signatures_1D_df = input_df_filter[['mutation_id', 'signature']]

final_1D_df = pd.merge(true_signatures_1D_df, tracksig_signatures_1D_df, on='mutation_id', how='inner', suffixes=['_true', '_tracksig'])
score_1D_tracksig = score_sig_1D_base(final_1D_df.signature_true.values, final_1D_df.signature_tracksig.values.astype(int))

# tracksigfreq
try:
    mixture_file = pd.read_csv('{}/tracksigfreq/tracksigfreq_mixtures_cancertype.csv'.
                               format(folder_path), sep=',')
except FileNotFoundError:
    print('non working')
try:
    changepoint_file = pd.read_csv(
        '{}/tracksigfreq/tracksigfreq_changepoints_cancertype.txt'.
        format(folder_path), header=None, sep=' ')
    changepoints_tracksig_list = changepoint_file.values[0]
except EmptyDataError:
    changepoints_tracksig_list = np.array(list())

data_df = pd.read_csv('{}/tracksigfreq/vcaf.csv'.
    format(folder_path), sep='\t')

cluster_id_list = np.zeros(data_df.bin.nunique())
i = 1
for chg_point in changepoints_tracksig_list:
    cluster_id_list[(chg_point - 1):] = i
    i += 1
data_df = data_df.assign(
    pred_cluster_id=data_df.apply(lambda x: int(cluster_id_list[x['bin']-1]),
                             axis=1))
J_pred = len(changepoints_tracksig_list) + 1
weights_pred = data_df.groupby('pred_cluster_id').phi.count().values/len(data_df)
phi_pred_values = data_df.groupby('pred_cluster_id').phi.mean().values

est_clonal_idx = data_df.groupby('pred_cluster_id').phi.mean().idxmax()
data_df = data_df.assign(
    pred_subclonal=(data_df.pred_cluster_id != est_clonal_idx).astype(int))

pred_signatures = np.zeros(len(all_sigs))

filename = '{}/subMU.csv'.format(folder_path)
sub_matrix = pd.read_csv(filename, sep='\t')
mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
sub_sigs = sub_matrix.columns[1:]
idx = [list(all_sigs).index(s) for s in sub_sigs]
est_sigs = mixture_file[mixture_file.columns[1:]].mean(axis=1).values
pred_signatures[idx] = est_sigs

pred_profile = est_sigs.dot(mu_mat_setting)

mut_sig = np.moveaxis(
    np.repeat([mixture_file.values[:, 1:].T[data_df.bin.astype(int)-1]],
              96, axis=0), [0, 1, 2], [2, 0, 1])
big_mu = np.repeat([mu_mat_setting], len(data_df), axis=0)
big_everything = (mut_sig * big_mu / len(est_sigs))[np.arange(len(data_df)), :, data_df.trinucleotide]
signature_mut = np.argmax(big_everything, axis=1)

tracksigfreq_signatures_1D_df = input_df_filter[['mutation_id', 'signature']]

final_1D_df = pd.merge(true_signatures_1D_df, tracksigfreq_signatures_1D_df, on='mutation_id', how='inner', suffixes=['_true', '_tracksigfreq'])
score_1D_tracksigfreq = score_sig_1D_base(final_1D_df.signature_true.values, final_1D_df.signature_tracksigfreq.values.astype(int))

old_score_filename = [i for i in os.listdir(folder_path) if 'result_evaluation' in i][0]
old_score = pd.read_csv('{}/{}'.format(folder_path, old_score_filename), sep='\t')
old_score.loc[old_score.method=='tracksig', 'score_sig_1D'] = score_1D_tracksig
old_score.loc[old_score.method=='tracksigfreq', 'score_sig_1D'] = score_1D_tracksigfreq
new_name = old_score_filename.split('.')[0] + '_new.csv'
old_score.to_csv('{}/{}'.format(folder_path, new_name),
              sep='\t', index=False)



