#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pickle
import scipy as sp
from clonesig.data_loader import SimLoader
from clonesig.evaluate import score1B_base, score1C_base, score2A_base, score2C_base


def score1B(sim, cluster_table):
    J_true = len(sim.phi)
    J_pred = len(cluster_table)
    return score1B_base(J_true, J_pred)


def score1C(sim, cluster_table):
    sim_w = np.unique(sim.U, return_counts=True)[1]/sim.N

    est_w = cluster_table['no.of.mutations'] / cluster_table['no.of.mutations'].sum()
    est_phi = cluster_table['location']
    return score1C_base(sim.phi, est_phi, sim_w, est_w)


def score2A(sim, loci_table):
    true_df = sim._get_data_df()
    return score2A_base(sim.U, loci_table.cluster)


def score2C(sim, loci_table, cluster_table):
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
    est_att = loci_table.cluster
    est_clonal_idx = cluster_table.sort_values(by='location').iloc[-1]['cluster.no']
    est_binary = (est_att != est_clonal_idx).astype(int)
    return score2C_base(sim_binary, est_binary)


if __name__ == '__main__':
    folder_path = sys.argv[1]
    start_time = int(sys.argv[2])
    end_time = int(sys.argv[3])

    """
    folder_path = "20200210_simulations_clonesig_cn_cancer_type/type0-perc_diploid100-nb_clones1-nb_mut600"
    start_time = 100
    end_time = 234
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

    res_folder = '{}_DPoutput_2000iters_1000burnin_seed123'.format(folder_path.split('/')[1])
    loci_table = pd.read_csv('{}/dpclust/{}/{}_2000iters_1000burnin_bestConsensusAssignments.bed'.format(folder_path, res_folder, folder_path.split('/')[1]), sep='\t')
    cluster_table = pd.read_csv('{}/dpclust/{}/{}_2000iters_1000burnin_bestClusterInfo.txt'.format(folder_path, res_folder, folder_path.split('/')[1]), sep='\t')


    id_list = list()
    metrics_list = list()
    id_list.append([cancer_type, perc_diploid, nb_clones, nb_mut,
                    np.nan, np.nan, min_dist, max_dist, avg_dist,
                    avg_major_cn, actual_perc_diploid, avg_tot_cn, 'DPClust', "all", np.nan])
    sml = list()
    sml.append(len(cluster_table))
    sml.append(np.nan)
    sml.append(np.nan)
    sml.append(score1B(sim_data_obj, cluster_table[cluster_table['no.of.mutations']>1]))
    sml.append(score1C(sim_data_obj, cluster_table))
    sml.append(score2A(sim_data_obj, loci_table))
    auc, accuracy, sensitivity, specificity, precision = score2C(sim_data_obj, loci_table, cluster_table)
    for v in (auc, accuracy, sensitivity, specificity, precision):
        sml.append(v)
    # signature related metrics
    sml += [np.nan] * 14

    sml.append(end_time-start_time)
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
    res_df.to_csv('{}/eval_dpclust.tsv'.format(folder_path), sep='\t',
                  index=False)


