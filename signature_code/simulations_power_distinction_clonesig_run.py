#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd
import os
import sys
from sklearn import preprocessing
from collections import Iterable
import numpy as np
from util_functions import safe_mkdir
import pickle
from clonesig import mixin_init_parameters
from clonesig.estimator import Estimator, log_binomial_coeff, _beta_binomial_logpmf, beta_binomial_pmf, EV_DOF_THRESHOLD
from clonesig.run_clonesig import run_clonesig, get_MU
from signature_code.evaluate_clonesig_simu import score2C as score2C_clonesig
from signature_code.evaluate_palimpsest import score2C as score2C_palimpsest
from signature_code.evaluate_tracksig import score2C as score2C_tracksig
from signature_code.evaluate_tracksigfreq import score2C as score2C_tracksigfreq
import time
import scipy as sp
from scipy.spatial.distance import cosine, pdist, euclidean
from scipy.stats import wasserstein_distance, pearsonr
import itertools
from sklearn import linear_model
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics.cluster import v_measure_score
import pkg_resources
from pandas.errors import EmptyDataError

folder_path = sys.argv[1]
MIXTURE_THRESHOLD = 0.05

'''
folder_path = '20190623_simulations_clonesig_cn_cancer_type/type2-perc_diploid100-nb_clones4-nb_mut100'
'''


nb_pi = int(folder_path.split('pi')[1].split('-')[0])
nb_phi = int(folder_path.split('phi')[1].split('-')[0])
depth = int(folder_path.split('depth')[1].split('-')[0])
nb_mut = int(folder_path.split('nb_mut')[1].split('-')[0])
perc_diploid = float(folder_path.split('percdip')[1].split('-')[0])
nb_clones = 2
cancer_type = None

MU = get_MU()

data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
with open('{}/purity.txt'.format(folder_path), 'r') as f:
    purity = float(f.read())

# get metrics from simulated data
with open('{}/sim_data'.format(folder_path), 'rb') as sim_pickle_file:
    sim_pickle = pickle.Unpickler(sim_pickle_file)
    sim_data_obj = sim_pickle.load()
dist_matrix = sp.spatial.distance.squareform(sp.spatial.distance.pdist(sim_data_obj.pi.dot(sim_data_obj.MU), 'cosine'))
if len(dist_matrix[dist_matrix > 0]):
    min_dist = np.min(dist_matrix[dist_matrix > 0])
    max_dist = np.max(dist_matrix[dist_matrix > 0])
    avg_dist = np.mean(dist_matrix[dist_matrix > 0])
else:
    min_dist, max_dist, avg_dist = np.nan, np.nan, np.nan
avg_major_cn = np.mean(sim_data_obj.C_tumor_tot - sim_data_obj.C_tumor_minor)
avg_tot_cn = np.mean(sim_data_obj.C_tumor_tot)
actual_perc_diploid = sum((sim_data_obj.C_tumor_tot == 2) & (sim_data_obj.C_tumor_minor == 1))/sim_data_obj.N
phi_dist = sim_data_obj.phi[0] - sim_data_obj.phi[1]
sim_data_obj.xi = np.array(sim_data_obj.xi)

# run and evaluate clonesig
method = 'clonesig'
id_list = list()
metrics_list = list()
for setting in ('prefit', 'all'):
    MU = get_MU()
    model_selection_kws = {'factor': 0.034}

    nuh = None
    pf = False
    if setting == 'prefit':
        model_selection_kws = {'factor': 0.065}
        pf = True
    start = time.time()
    new_est, lr, pval, new_inputMU, cst_est, fitted_sigs = run_clonesig(
        data_df.trinucleotide.values, data_df.var_counts.values,
        data_df.var_counts.values + data_df.ref_counts.values,
        data_df.normal_cn.values,
        data_df.minor_cn.values + data_df.major_cn.values,
        data_df.minor_cn.values, purity, MU, inputNu=None, nu_heuristics=nuh,
        return_sig_change_test=True, min_mut_clone=0, min_prop_sig=0.0,
        prefit_signatures=pf, prefit_thresh=0.01, model_selection_function=None,
        model_selection_kws=model_selection_kws, max_nb_clones=10)
    end = time.time()
    ev, _ = np.linalg.eig(1-sp.spatial.distance.squareform(sp.spatial.distance.pdist(new_inputMU, 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)

    id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                    MU.shape[0], new_inputMU.shape[0], min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, method, setting, dof,
                    phi_dist])
    sml = list()
    sml.append(new_est.J)
    sml.append(lr)
    sml.append(pval)
    auc, accuracy, sensitivity, specificity, precision = score2C_clonesig(sim_data_obj, new_est)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    metrics_list.append(sml)

# evaluate tracksig
mu_mat = {'all': MU}
for setting in ('all', 'prefit'):
    try:
        mixture_file = pd.read_csv('{}/tracksig/tracksig_mixtures_{}.csv'.
                                   format(folder_path, setting), sep=',')
    except FileNotFoundError:
        id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                        MU.shape[0], np.nan, min_dist, max_dist, avg_dist,
                        avg_major_cn, actual_perc_diploid, avg_tot_cn, "tracksig", setting, dof,
                        phi_dist])  
        sml = [np.nan] * 8
        metrics_list.append(sml)
        continue
    try:
        changepoint_file = pd.read_csv(
            '{}/tracksig/tracksig_changepoints_{}.txt'.
            format(folder_path, setting), header=None, sep=' ')
        changepoints_tracksig_list = changepoint_file.values[0]
    except EmptyDataError:
        changepoints_tracksig_list = np.array(list())
    if setting == 'prefit':
        premixture_file = pd.read_csv(
            '{}/tracksig/tracksig_premixtures_{}.csv'.
            format(folder_path, setting), sep=',')
        select_sig = premixture_file.x.values > MIXTURE_THRESHOLD
        mu_mat[setting] = MU[select_sig, :]
    ev, _ = np.linalg.eig(
        1-sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(mu_mat[setting], 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)
    data_df = sim_data_obj._get_data_df()
    data_df = data_df.assign(vaf=data_df.var_counts / (data_df.ref_counts + data_df.var_counts))
    data_df = data_df.assign(total_cn=lambda x: x['minor_cn'] + x['major_cn'])
    data_df = data_df.assign(vaf_cn=data_df.vaf * data_df['total_cn'] / data_df['mut_cn'])
    data_df = data_df.assign(vaf_purity=data_df.apply(lambda x: x['vaf']/sim_data_obj.purity * ((1 - sim_data_obj.purity) * 2 + sim_data_obj.purity * x['total_cn']) / x['mut_cn'], axis=1))
    data_df = data_df.assign(old_idx=data_df.index)
    data_df.sort_values(by='vaf_purity', inplace=True)
    data_df.reset_index(inplace=True, drop=True)

    data_df = data_df.assign(mutation_group=lambda x: x.index//100)
    data_df.sort_values(by='old_idx', inplace=True)
    cluster_id_list = np.zeros(data_df.mutation_group.nunique())
    i = 1
    for chg_point in changepoints_tracksig_list:
        cluster_id_list[(chg_point - 1):] = i
        i += 1
    data_df = data_df.assign(cluster_id=data_df.apply(lambda x: int(cluster_id_list[x['mutation_group']]), axis=1))

    if setting == 'prefit':
        nb_initial_sigs = 47
    else:
        nb_initial_sigs = mu_mat[setting].shape[0]
    id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                    MU.shape[0], mu_mat[setting].shape[0], min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, "tracksig", setting, dof,
                    phi_dist])

    sml = list()
    sml.append(len(changepoints_tracksig_list) + 1)
    sml.append(np.nan)
    sml.append(np.nan)
    auc, accuracy, sensitivity, specificity, precision = score2C_tracksig(sim_data_obj,
                                                                 data_df)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    metrics_list.append(sml)

# evaluate tracksigfreq
for setting in ('all', 'prefit'):
    try:
        mixture_file = pd.read_csv('{}/tracksigfreq/tracksigfreq_mixtures_{}.csv'.
                                   format(folder_path, setting), sep=',')
    except FileNotFoundError:
        id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                        MU.shape[0], np.nan, min_dist, max_dist, avg_dist,
                        avg_major_cn, actual_perc_diploid, avg_tot_cn, "tracksigfreq", setting, dof,
                        phi_dist])  
        sml = [np.nan] * 8
        metrics_list.append(sml)
        continue
    try:
        changepoint_file = pd.read_csv(
            '{}/tracksigfreq/tracksigfreq_changepoints_{}.txt'.
            format(folder_path, setting), header=None, sep=' ')
        changepoints_tracksig_list = changepoint_file.values[0]
    except EmptyDataError:
        changepoints_tracksig_list = np.array(list())
    if setting == 'prefit':
        premixture_file = pd.read_csv(
            '{}/tracksigfreq/tracksigfreq_premixtures_{}.csv'.
            format(folder_path, setting), sep=',')
        select_sig = premixture_file.x.values > MIXTURE_THRESHOLD
        mu_mat[setting] = MU[select_sig, :]
    ev, _ = np.linalg.eig(
        1-sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(mu_mat[setting], 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)

    data_df = pd.read_csv('{}/tracksigfreq/vcaf.csv'.
        format(folder_path), sep='\t')

    cluster_id_list = np.zeros(data_df.bin.nunique())
    i = 1
    for chg_point in changepoints_tracksig_list:
        cluster_id_list[(chg_point - 1):] = i
        i += 1
    data_df = data_df.assign(cluster_id=data_df.apply(lambda x: int(cluster_id_list[x['bin']-1]), axis=1))

    if setting == 'prefit':
        nb_initial_sigs = 47
    else:
        nb_initial_sigs = mu_mat[setting].shape[0]
    id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                    MU.shape[0], mu_mat[setting].shape[0], min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, "tracksigfreq", setting, dof,
                    phi_dist])

    sml = list()
    sml.append(len(changepoints_tracksig_list) + 1)
    sml.append(np.nan)
    sml.append(np.nan)
    auc, accuracy, sensitivity, specificity, precision = score2C_tracksigfreq(sim_data_obj,
                                                                 data_df)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    metrics_list.append(sml)

# evaluate palimpsest
sig_match_table = pd.read_csv(pkg_resources.resource_stream(
        'clonesig', 'data/curated_match_signature_cancertype_tcgawes_literature.csv'), sep='\t', index_col=0)
sig_names = sig_match_table[sig_match_table.sum(axis=1)>0].index.tolist()

for setting in ('all', 'prefit'):
    try:
        mixture_file = pd.read_csv('{}/palimpsest/palimpsest_mixtures_{}.csv'.
                                   format(folder_path, setting), sep='\t')
        ccf_file = pd.read_csv('{}/palimpsest/palimpsest_mut_data_{}.csv'
                               .format(folder_path, setting), sep='\t')
    except FileNotFoundError:
        id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                        MU.shape[0], np.nan, min_dist, max_dist, avg_dist,
                        avg_major_cn, actual_perc_diploid, avg_tot_cn, "palimpsest", setting, dof,
                        phi_dist])  
        sml = [np.nan] * 8
        metrics_list.append(sml)
        continue
    if setting == 'prefit':
        premixture_file = pd.read_csv(
            '{}/palimpsest/palimpsest_premixtures_{}.txt'.
            format(folder_path, setting), sep=' ')
        select_sig = premixture_file.columns.tolist()
        select_sig_idx = [sig_names.index(s) for s in select_sig]
        mu_mat[setting] = MU[select_sig_idx, :]
    ev, _ = np.linalg.eig(
        1-sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(mu_mat[setting], 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)

    if setting == 'prefit':
        nb_initial_sigs = 47
    else:
        nb_initial_sigs = mu_mat[setting].shape[0]
    id_list.append([nb_pi, nb_phi, depth, perc_diploid, nb_clones, nb_mut,
                    MU.shape[0], mu_mat[setting].shape[0], min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, "palimpsest", setting, dof,
                    phi_dist])

    sml = list()
    sml.append(2)
    sml.append(np.nan)
    sml.append(np.nan)
    ccf_file = ccf_file.assign(
        clonality_binary=ccf_file.apply(
            lambda x: 1 if x['Clonality'] == 'subclonal' else 0, axis=1))
    auc, accuracy, sensitivity, specificity, precision = score2C_palimpsest(sim_data_obj,
                                                                 ccf_file)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    metrics_list.append(sml)


sample_id_cols = ['nb_pi', 'nb_phi', 'depth', 'perc_diploid', 'nb_clones', 'nb_mut',
                  'nb_sig', 'nb_sig_fit', 'min_dist', 'max_dist', 'avg_dist',
                  'avg_major_cn', 'actual_perc_diploid', 'avg_tot_cn',
                  'method', 'setting', 'dof', 'phi_dist']
metrics_cols = ['fitted_nb_clones', 'll_ratio', 'pval', 'score2C_auc', 'score2C_accuracy',
                'score2C_sensitivity', 'score2C_specificity',
                'score2C_precision']
id_df = pd.DataFrame(id_list, columns=sample_id_cols)
metrics_df = pd.DataFrame(metrics_list, columns=metrics_cols)

res_df = pd.concat([id_df, metrics_df], axis=1)

res_df.to_csv('{}/eval_power.tsv'.format(folder_path), sep='\t',
              index=False)

