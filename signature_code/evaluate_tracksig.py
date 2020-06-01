#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pickle
import scipy as sp
from clonesig.data_loader import SimLoader
from clonesig.evaluate import score1B_base, score1C_base, score2A_base, score2C_base, score_sig_1A_base, score_sig_1B_base, score_sig_1C_base, score_sig_1D_base, score_sig_1E_base
from clonesig.run_clonesig import get_MU
from clonesig.estimator import Estimator, EV_DOF_THRESHOLD
import pkg_resources
from pandas.errors import EmptyDataError
from clonesig import mixin_init_parameters


def score1B(sim, J_pred):
    J_true = len(sim.phi)
    return score1B_base(J_true, J_pred)


def score1C(sim, data_df):
    sim_w = np.unique(sim.U, return_counts=True)[1]/sim.N

    est_w = data_df.groupby('cluster_id').vaf_purity.count().values/len(data_df)
    est_phi = data_df.groupby('cluster_id').vaf_purity.mean().values
    return score1C_base(sim.phi, est_phi, sim_w, est_w)


def score2A(sim, data_df):
    true_df = sim._get_data_df()
    ordered_loci_table = pd.merge(true_df, data_df, on='mutation_id')
    return score2A_base(sim.U, ordered_loci_table.cluster_id)


def score2C(sim, data_df):
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

    true_df = sim._get_data_df()
    ordered_loci_table = pd.merge(true_df, data_df, on='mutation_id')
    est_att = ordered_loci_table.cluster_id
    est_clonal_idx = data_df.groupby('cluster_id').vaf_purity.mean().idxmax()
    est_binary = (est_att != est_clonal_idx).astype(int)
    return score2C_base(sim_binary, est_binary)


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


def score_sig_1C(sim, est_sig, cancer_type=None, fitted_sigs=None):
    """
    precision recall for detected signatures
    """
    sim_sig = sim.xi.dot(sim.pi)

    filter_filename = 'data/curated_match_signature_cancertype_tcgawes_literature.csv'
    cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
        'clonesig', filter_filename), sep='\t', index_col=0).values
    select = cancer_type_sig[:, cancer_type].astype(bool)
    big_select = cancer_type_sig.sum(axis=1).astype(bool)
    if fitted_sigs is not None:
        selectb = np.zeros(65)
        selectb[big_select] = fitted_sigs
        selectb = selectb.astype(bool)
    else:
        selectb = select
    if sim.MU.shape[0] == 65:
        big_pi_sim = sim_sig
    else:
        # ids of simulated signatures
        big_pi_sim = np.zeros(65)
        big_pi_sim[select] = sim_sig
    if len(est_sig) == 65:
        big_pi_est = est_sig
    elif len(est_sig) == 47:
        big_pi_est = np.zeros(65)
        big_pi_est[big_select] = est_sig
    else:
        big_pi_est = np.zeros(65)
        big_pi_est[selectb] = est_sig
    sim_sig_idx = np.where(big_pi_sim > mixin_init_parameters.ZERO_PADDING)[0]

    # do the computation
    final_true_sig = np.zeros(65)
    final_true_sig[sim_sig_idx] = 1

    return score_sig_1C_base(final_true_sig, big_pi_est)


def score_sig_1D(sim, est_sig_att, inputMU, cancer_type=None, fitted_sigs=None):
    """
    percent of mutations with the right signature
    """
    sim_sig = sim.S.copy()
    if sim.MU.shape != inputMU.shape:
        filter_filename = 'data/curated_match_signature_cancertype_tcgawes_literature.csv'
        cancer_type_sig = pd.read_csv(pkg_resources.resource_stream(
            'clonesig', filter_filename), sep='\t', index_col=0).values
        select = cancer_type_sig[:, cancer_type].astype(bool)
        big_select = cancer_type_sig.sum(axis=1).astype(bool)
        if fitted_sigs is not None:
            selectb = np.zeros(65)
            selectb[big_select] = fitted_sigs
            selectb = selectb.astype(bool)
        else:
            selectb = select
        if sim.MU.shape[0] != 65:
            sim_sig = np.array([np.where(select)[0][int(i)] for i in sim_sig])
        if inputMU.shape[0] == 47:
            est_sig_att = np.array([np.where(big_select)[0][int(i)] for i in est_sig_att])
        elif inputMU.shape[0] < 47:
            est_sig_att = np.array([np.where(selectb)[0][int(i)] for i in est_sig_att])
    return score_sig_1D_base(sim_sig, est_sig_att)


def score_sig_1E(sim, est_dist):
    """
    cosine dist between the clone distrib that generated the mutation and the
    reconstituted one.
    """
    true_dist = sim.pi[sim.U, :].dot(sim.MU)
    return score_sig_1E_base(true_dist, est_dist)


if __name__ == '__main__':
    folder_path = sys.argv[1]
    MIXTURE_THRESHOLD = 0.05

    """
    folder_path = "20190623_simulations_clonesig_cn_cancer_type/type0-perc_diploid100-nb_clones3-nb_mut1000"
    """

    cancer_type = int(folder_path.split('/')[1].split('type')[1].split('-')[0])
    perc_diploid = int(folder_path.split('/')[1].split('perc_diploid')[1].split('-')[0])
    nb_clones = int(folder_path.split('/')[1].split('nb_clones')[1].split('-')[0])
    nb_mut = int(folder_path.split('/')[1].split('nb_mut')[1].split('-')[0])

    # get metrics from simulated data
    with open('{}/sim_data'.format(folder_path), 'rb') as sim_pickle_file:
        sim_pickle = pickle.Unpickler(sim_pickle_file)
        sim_data_obj = sim_pickle.load()
    dist_matrix = sp.spatial.distance.squareform(
        sp.spatial.distance.pdist(sim_data_obj.pi.dot(sim_data_obj.MU), 'cosine'))
    if len(dist_matrix[dist_matrix > 0]):
        min_dist = np.min(dist_matrix[dist_matrix > 0])
        max_dist = np.max(dist_matrix[dist_matrix > 0])
        avg_dist = np.mean(dist_matrix[dist_matrix > 0])
    else:
        min_dist, max_dist, avg_dist = np.nan, np.nan, np.nan
    avg_major_cn = np.mean(sim_data_obj.C_tumor_tot - sim_data_obj.C_tumor_minor)
    avg_tot_cn = np.mean(sim_data_obj.C_tumor_tot)
    actual_perc_diploid = sum((sim_data_obj.C_tumor_tot == 2) &
                              (sim_data_obj.C_tumor_minor == 1)) / sim_data_obj.N

    MU = get_MU()
    subMU = get_MU(cancer_type=cancer_type)




    id_list = list()
    metrics_list = list()
    mu_mat = {'all': MU, 'cancertype': subMU}
    for setting in ('all', 'cancertype', 'prefit'):
        mixture_file = pd.read_csv('{}/tracksig/tracksig_mixtures_{}.csv'.
                                   format(folder_path, setting), sep=',')
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
            nb_initial_sigs = 65
        else:
            nb_initial_sigs = mu_mat[setting].shape[0]
        id_list.append([cancer_type, perc_diploid, nb_clones, nb_mut,
                        nb_initial_sigs, mu_mat[setting].shape[0], min_dist, max_dist,
                        avg_dist, avg_major_cn, actual_perc_diploid, avg_tot_cn,
                        'tracksig', setting, dof])

        sml = list()
        sml.append(len(changepoints_tracksig_list) + 1)
        sml.append(np.nan)
        sml.append(np.nan)
        sml.append(score1B(sim_data_obj, len(changepoints_tracksig_list) + 1))
        sml.append(score1C(sim_data_obj, data_df))
        sml.append(score2A(sim_data_obj, data_df))
        auc, accuracy, sensitivity, specificity, precision = score2C(sim_data_obj,
                                                                     data_df)
        for v in (auc, accuracy, sensitivity, specificity, precision):
            sml.append(v)

        est_distrib = mixture_file.values[:, 1:].T.dot(mu_mat[setting]).mean(axis=0)
        sml.append(score_sig_1A(sim_data_obj, est_distrib))
        sml.append(score_sig_1B(sim_data_obj, est_distrib))
        if setting == 'prefit':
            fitted_sigs = (premixture_file.x > 0.05).values
        else:
            fitted_sigs = None
        est_sigs = mixture_file.values[:, 1:].mean(axis=1)
        auc, accuracy, sensitivity, specificity, precision = score_sig_1C(
            sim_data_obj, est_sigs, cancer_type=cancer_type,
            fitted_sigs=fitted_sigs)
        for v in (auc, accuracy, sensitivity, specificity, precision):
            sml.append(v)
        mut_sig = np.moveaxis(
            np.repeat([mixture_file.values[:, 1:].T[data_df.mutation_group]],
                      96, axis=0), [0, 1, 2], [2, 0, 1])
        big_mu = np.repeat([mu_mat[setting]], sim_data_obj.N, axis=0)
        big_everything = (mut_sig * big_mu / len(est_sigs))[np.arange(sim_data_obj.N), :, data_df.trinucleotide]
        signature_mut = np.argmax(big_everything, axis=1)

        sml.append(score_sig_1D(sim_data_obj, signature_mut,  mu_mat[setting],
                                cancer_type, fitted_sigs))
        est_mut_profile = mixture_file.values[:, 1:].T[data_df.mutation_group].dot(mu_mat[setting])
        (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
            median_diff_distrib_mut, perc_dist_5, perc_dist_10) = score_sig_1E(
            sim_data_obj, est_mut_profile.astype(float))
        for v in (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
                  median_diff_distrib_mut, perc_dist_5, perc_dist_10):
            sml.append(v)
        fit_time = pd.read_csv('{}/tracksig/tracksig_runtime_{}.csv'.format(folder_path, setting),
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
    res_df.to_csv('{}/eval_tracksig.tsv'.format(folder_path), sep='\t',
                  index=False)
