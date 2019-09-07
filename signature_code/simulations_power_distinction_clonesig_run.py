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
from clonesig.evaluate import score1B_base, score1C_base, score2A_base, score2C_base, score_sig_1A_base, score_sig_1B_base, score_sig_1C_base, score_sig_1D_base, score_sig_1E_base
import time
import scipy as sp
from scipy.spatial.distance import cosine, pdist, euclidean
from scipy.stats import wasserstein_distance, pearsonr
import itertools
from sklearn import linear_model
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics.cluster import v_measure_score
import pkg_resources

folder_path = sys.argv[1]


'''
folder_path = '20190623_simulations_clonesig_cn_cancer_type/type2-perc_diploid100-nb_clones4-nb_mut100'
'''

"""
1,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid0-nb_clones1-nb_mut100-true_c_mut100
2,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid20-nb_clones1-nb_mut100-true_c_mut100
3,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid40-nb_clones1-nb_mut100-true_c_mut100
4,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid60-nb_clones1-nb_mut100-true_c_mut100
5,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid80-nb_clones1-nb_mut100-true_c_mut100
6,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid100-nb_clones1-nb_mut100-true_c_mut100
7,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid0-nb_clones1-nb_mut300-true_c_mut100
8,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid20-nb_clones1-nb_mut300-true_c_mut100
9,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid40-nb_clones1-nb_mut300-true_c_mut100
10,20190528_simulations_clonesig_cn_cancer_type/type0-perc_diploid60-nb_clones1-nb_mut300-true_c_mut100
(my_python36_centos) -bash-4.2$ grep "^366," 20190618_eval_clonesig.csv
366,20190528_simulations_clonesig_cn_cancer_type/type2-perc_diploid100-nb_clones4-nb_mut100-true_c_mut100
20190528_simulations_clonesig_cn_cancer_type/type3-perc_diploid60-nb_clones1-nb_mut600-true_c_mut100
"""

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

def score1B(sim, est):
    J_true = len(sim.phi)
    J_pred = len(est.phi)
    return score1B_base(J_true, J_pred)


def score1C(sim, est):
    sim_w = np.unique(sim.U, return_counts=True)[1]/sim.N
    pre_est_w = np.zeros(est.J)
    pre_counts = np.unique(np.argmax(est.qun, axis=1),
                           return_counts=True)
    pre_est_w[pre_counts[0]] = pre_counts[1]
    est_w = pre_est_w/len(sim.U)
    return score1C_base(sim.phi, est.phi, sim_w, est_w)


def score2A(sim, est):
    return score2A_base(sim.U, np.argmax(est.qun, axis=1))


def score2C(sim, est):
    """
    this metric compares the sensitivity, specificity, precision of mutations
    classified as clonal vs subclonal. We take the convention that the clone
    with the largest phi is clonal and the rest is subclonal
    we take also the convention that subclonal are the positive examples, and
    clonal the negative ones.
    """
    sim_att = sim.U
    sim_clonal_idx = np.argmax(sim.phi)
    sim_binary = (sim_att != sim_clonal_idx).astype(int)

    est_att = np.argmax(est.qun, axis=1)
    est_clonal_idx = np.argmax(est.phi)
    est_binary = (est_att != est_clonal_idx).astype(int)
    return score2C_base(sim_binary, est_binary)


def score_sig_1A(sim, est):
    """
    euclidian norm between normalized trinucleotide context counts (empirical),
    and the reconstituted profile
    """
    raw_data_distrib = np.zeros(96)
    val, c = np.unique(sim.T, return_counts=True)
    raw_data_distrib[val.astype(int)] = c
    raw_data_distrib = raw_data_distrib / sim.N
    est_distrib = est.xi.dot(est.pi).dot(est.mu_matrix)
    return score_sig_1A_base(raw_data_distrib, est_distrib)


def score_sig_1B(sim, est):
    """
    euclidian norm between simulated and estimated signature profile (summed
    over all clones)
    """
    sim_distrib = sim.xi.dot(sim.pi).dot(sim.MU)
    est_distrib = est.xi.dot(est.pi).dot(est.mu_matrix)
    return score_sig_1B_base(sim_distrib, est_distrib)


def score_sig_1C(sim, est, cancer_type=None, fitted_sigs=None):
    """
    precision recall for detected signatures
    """
    sim_sig = sim.xi.dot(sim.pi)
    est_sig = est.xi.dot(est.pi)


    filter_filename = 'data/match_cancer_type_sig_v3.csv'
    cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
        'clonesig', filter_filename), index_col=0).values
    select = cancer_type_sig[cancer_type, :].astype(bool)
    if fitted_sigs is not None:
        selectb = fitted_sigs
    else:
        selectb = select
    if sim.MU.shape[0] == 65:
        big_pi_sim = sim_sig
    else:
        # ids of simulated signatures
        big_pi_sim = np.zeros(65)
        big_pi_sim[select] = sim_sig
    if est.mu_matrix.shape[0] == 65:
        big_pi_est = est_sig
    else:
        big_pi_est = np.zeros(65)
        big_pi_est[selectb] = est_sig
    sim_sig_idx = np.where(big_pi_sim > mixin_init_parameters.ZERO_PADDING)[0]

    # do the computation
    final_true_sig = np.zeros(65)
    final_true_sig[sim_sig_idx] = 1

    return score_sig_1C_base(final_true_sig, big_pi_est)


def score_sig_1D(sim, est, cancer_type=None, fitted_sigs=None):
    """
    percent of mutations with the right signature
    """
    est_sig = est.rnus[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1)
    sim_sig = sim.S.copy()
    if sim.MU.shape != est.mu_matrix.shape:
        filter_filename = 'data/match_cancer_type_sig_v3.csv'
        cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
            'clonesig', filter_filename), index_col=0).values
        select = cancer_type_sig[cancer_type, :].astype(bool)
        if fitted_sigs is not None:
            selectb = fitted_sigs
        else:
            selectb = select
        if sim.MU.shape[0] != 65:
            sim_sig = np.array([np.where(select)[0][int(i)] for i in sim_sig])
        if est.mu_matrix.shape[0] != 65:
            est_sig = np.array([np.where(selectb)[0][int(i)] for i in est_sig])
    return score_sig_1D_base(sim_sig, est_sig)


def score_sig_1E(sim, est):
    """
    cosine dist between the clone distrib that generated the mutation and the
    reconstituted one.
    """
    true_dist = sim.pi[sim.U, :].dot(sim.MU)
    est_dist = est.pi[est.qun.argmax(axis=1), :].dot(est.mu_matrix)
    return score_sig_1E_base(true_dist, est_dist)

method = 'clonesig'
id_list = list()
metrics_list = list()
for setting in ('prefit', 'all'):
    if setting == 'cancer_type':
        MU = get_MU(cancer_type=cancer_type)
        model_selection_kws = {'factor':  0.093}
    else:
        MU = get_MU()
        model_selection_kws = {'factor': 0.048}

    if setting == 'all_nuclonal':
        nuh = 'clonal'
    elif setting == 'all_minor':
        nuh = 'minor'
    else:
        nuh = None
    pf = False
    if setting == 'prefit':
        # model_selection_kws = {'factor': 0.022}
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
    sml.append(score1B(sim_data_obj, new_est))
    sml.append(score1C(sim_data_obj, new_est))
    sml.append(score2A(sim_data_obj, new_est))
    auc, accuracy, sensitivity, specificity, precision = score2C(sim_data_obj, new_est)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    sml.append(score_sig_1A(sim_data_obj, new_est))
    sml.append(score_sig_1B(sim_data_obj, new_est))
    auc, accuracy, sensitivity, specificity, precision = score_sig_1C(
        sim_data_obj, new_est, cancer_type=cancer_type,
        fitted_sigs=fitted_sigs)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    sml.append(score_sig_1D(sim_data_obj, new_est, cancer_type=cancer_type,
                            fitted_sigs=fitted_sigs))
    (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
        median_diff_distrib_mut, perc_dist_5, perc_dist_10) = score_sig_1E(
        sim_data_obj, new_est)
    for v in (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
              median_diff_distrib_mut, perc_dist_5, perc_dist_10):
        sml.append(v)
    sml.append(end-start)
    metrics_list.append(sml)
#    with open('{}/{}_clonesig_raw_results'
#              .format(folder_path, setting), 'wb') as raw_res:
#        if new_est.mu_matrix.shape[0] == 65:
#            new_est.mu_matrix = None
#            cst_est.mu_matrix = None
#        my_pickler = pickle.Pickler(raw_res)
#        my_pickler.dump([new_est, lr, pval, cst_est, fitted_sigs])


sample_id_cols = ['nb_pi', 'nb_phi', 'depth', 'perc_diploid', 'nb_clones', 'nb_mut',
                  'nb_sig', 'nb_sig_fit', 'min_dist', 'max_dist', 'avg_dist',
                  'avg_major_cn', 'actual_perc_diploid', 'avg_tot_cn',
                  'method', 'setting', 'dof', 'phi_dist']
metrics_cols = ['fitted_nb_clones', 'll_ratio', 'pval', 'score1B', 'score1C',
                'score2A', 'score2C_auc', 'score2C_accuracy',
                'score2C_sensitivity', 'score2C_specificity',
                'score2C_precision', 'score_sig_1A', 'score_sig_1B',
                'score_sig_1C_auc', 'score_sig_1C_accuracy',
                'score_sig_1C_sensitivity', 'score_sig_1C_specificity',
                'score_sig_1C_precision', 'score_sig_1D',
                'min_diff_distrib_mut', 'max_diff_distrib_mut',
                'std_diff_distrib_mut', 'median_diff_distrib_mut',
                'perc_dist_5', 'perc_dist_10', 'runtime']
id_df = pd.DataFrame(id_list, columns=sample_id_cols)
metrics_df = pd.DataFrame(metrics_list, columns=metrics_cols)

res_df = pd.concat([id_df, metrics_df], axis=1)

res_df.to_csv('{}/eval_clonesig.tsv'.format(folder_path), sep='\t',
              index=False)
print('+', '{}/eval_clonesig.tsv'.format(folder_path), '+', sep='')
