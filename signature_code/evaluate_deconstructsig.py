#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pickle
import scipy as sp
from clonesig.data_loader import SimLoader
from clonesig.evaluate import score_sig_1A_base, score_sig_1B_base, score_sig_1C_base, score_sig_1D_base, score_sig_1E_base
from clonesig.run_clonesig import get_MU
from clonesig.estimator import Estimator, EV_DOF_THRESHOLD
from clonesig import mixin_init_parameters
import pkg_resources

folder_path = sys.argv[1]


"""
folder_path = "20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300"
"""

cancer_type = int(folder_path.split('/')[1].split('type')[1].split('-')[0])
perc_diploid = int(folder_path.split('/')[1].split('perc_diploid')[1].split('-')[0])
nb_clones = int(folder_path.split('/')[1].split('nb_clones')[1].split('-')[0])
nb_mut = int(folder_path.split('/')[1].split('nb_mut')[1].split('-')[0])

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

MU = get_MU()
subMU = get_MU(cancer_type=cancer_type)


def score_sig_1A(sim, est_distrib):
    """
    euclidian norm between normalized trinucleotide context counts (empirical),
    and the reconstituted profile
    """
    raw_data_distrib = np.zeros(96)
    val, c = np.unique(sim.T, return_counts=True)
    raw_data_distrib[val.astype(int)] = c
    raw_data_distrib = raw_data_distrib / sim.N
    return score_sig_1A_base(raw_data_distrib, est_distrib)


def score_sig_1B(sim, est_distrib):
    """
    euclidian norm between simulated and estimated signature profile (summed
    over all clones)
    """
    sim_distrib = sim.xi.dot(sim.pi).dot(sim.MU)
    return score_sig_1B_base(sim_distrib, est_distrib)


def score_sig_1C(sim, est_sig, cancer_type=None):
    """
    precision recall for detected signatures
    """
    sim_sig = sim.xi.dot(sim.pi)

    filter_filename = 'data/match_cancer_type_sig_v3.csv'
    cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
        'clonesig', filter_filename), index_col=0).values
    select = cancer_type_sig[cancer_type, :].astype(bool)
    if sim.MU.shape[0] == 65:
        big_pi_sim = sim_sig
    else:
        # ids of simulated signatures
        big_pi_sim = np.zeros(65)
        big_pi_sim[select] = sim_sig
    if len(est_sig) == 65:
        big_pi_est = est_sig
    else:
        big_pi_est = np.zeros(65)
        big_pi_est[select] = est_sig
    sim_sig_idx = np.where(big_pi_sim > mixin_init_parameters.ZERO_PADDING)[0]

    # do the computation
    final_true_sig = np.zeros(65)
    final_true_sig[sim_sig_idx] = 1

    return score_sig_1C_base(final_true_sig, big_pi_est)


def score_sig_1D(sim, est_sig, inputMU, cancer_type=None):
    """
    percent of mutations with the right signature
    """
    data_df = sim._get_data_df()
    est = Estimator(data_df.trinucleotide.values, data_df.var_counts.values,
                    data_df.normal_cn.values,
                    data_df.minor_cn.values + data_df.major_cn.values,
                    data_df.minor_cn.values,
                    data_df.var_counts.values + data_df.ref_counts.values,
                    sim.purity, 1, inputMU=inputMU, pi=est_sig.reshape(1, -1))
    est_sig_att = est.rnus[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1)
    sim_sig = sim.S.copy()
    if sim.MU.shape != est.mu_matrix.shape:
        filter_filename = 'data/match_cancer_type_sig_v3.csv'
        cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
            'clonesig', filter_filename), index_col=0).values
        select = cancer_type_sig[cancer_type, :].astype(bool)
        if sim.MU.shape[0] != 65:
            sim_sig = np.array([np.where(select)[0][int(i)] for i in sim_sig])
        if est.mu_matrix.shape[0] != 65:
            est_sig_att = np.array([np.where(select)[0][int(i)] for i in est_sig_att])
    return score_sig_1D_base(sim_sig, est_sig_att)


def score_sig_1E(sim, est_sig, inputMU):
    """
    cosine dist between the clone distrib that generated the mutation and the
    reconstituted one.
    """
    true_dist = sim.pi[sim.U, :].dot(sim.MU)
    est_dist = np.repeat([est_sig], sim.N, axis=0).dot(inputMU)
    return score_sig_1E_base(true_dist, est_dist)


id_list = list()
metrics_list = list()
mu_mat = {'all': MU, 'cancertype': subMU}
for setting in ('all', 'cancertype'):
    result_file = pd.read_csv('{}/deconstructsigs/signatures_{}.csv'.format(folder_path, setting), sep=' ')
    ev, _ = np.linalg.eig(1-sp.spatial.distance.squareform(sp.spatial.distance.pdist(mu_mat[setting], 'cosine')))
    dof = sum(ev > EV_DOF_THRESHOLD)

    id_list.append([cancer_type, perc_diploid, nb_clones, nb_mut,
                    result_file.shape[1], result_file.shape[1], min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, 'deconstructsigs', setting, dof])

    sml = [np.nan] * 11
    sml.append(score_sig_1A(sim_data_obj, result_file.values.dot(mu_mat[setting])))
    sml.append(score_sig_1B(sim_data_obj, result_file.values.dot(mu_mat[setting])))
    auc, accuracy, sensitivity, specificity, precision = score_sig_1C(
        sim_data_obj, result_file.values[0], cancer_type=cancer_type)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    sml.append(score_sig_1D(sim_data_obj, result_file.values[0],
                            mu_mat[setting], cancer_type=cancer_type))
    (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
        median_diff_distrib_mut, perc_dist_5, perc_dist_10) = score_sig_1E(
        sim_data_obj, result_file.values[0], mu_mat[setting])
    for v in (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
              median_diff_distrib_mut, perc_dist_5, perc_dist_10):
        sml.append(v)
    fit_time = pd.read_csv('{}/deconstructsigs/deconstructsig_runtime_{}.csv'.format(folder_path, setting),
                           index_col=0).values[0][0]
    sml.append(fit_time)
    metrics_list.append(sml)

sample_id_cols = ['cancer_type', 'perc_diploid', 'nb_clones', 'nb_mut',
                  'nb_sig', 'nb_sig_fit', 'min_dist', 'max_dist', 'avg_dist',
                  'avg_major_cn', 'actual_perc_diploid', 'avg_tot_cn',
                  'method', 'setting', 'dof']
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

res_df.to_csv('{}/eval_deconstructsigs.tsv'.format(folder_path), sep='\t',
              index=False)

