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


folder_path = sys.argv[1]
start_time = int(sys.argv[2])
end_time = int(sys.argv[3])

"""
folder_path = "20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300"
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

loci_table = pd.read_csv('{}/ccube/ssm_clusters.csv'.format(folder_path), sep='\t')
pre_cluster_table = loci_table.groupby('ccube_ccf_mean').rough_mult.count()
loci_table = loci_table.assign(cluster_id=loci_table.apply(lambda x: pre_cluster_table.index.tolist().index(x['ccube_ccf_mean']), axis=1))


def score1B(sim, loci_table):
    J_true = len(sim.phi)
    J_pred = loci_table.ccube_ccf_mean.nunique()
    return score1B_base(J_true, J_pred)


def score1C(sim, loci_table):
    sim_w = np.unique(sim.U, return_counts=True)[1]/sim.N
    pre_cluster_table = loci_table.groupby('ccube_ccf_mean').rough_mult.count()
    est_w = pre_cluster_table.values/len(loci_table)
    est_phi = pre_cluster_table.index.values
    return score1C_base(sim.phi, est_phi, sim_w, est_w)


def score2A(sim, loci_table):
    return score2A_base(sim.U, loci_table.cluster_id)


def score2C(sim, loci_table):
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

    est_att = loci_table.cluster_id
    est_clonal_idx = np.argmax(pre_cluster_table.index.tolist())
    est_binary = (est_att != est_clonal_idx).astype(int)
    return score2C_base(sim_binary, est_binary)

id_list = list()
metrics_list = list()
id_list.append([cancer_type, perc_diploid, nb_clones, nb_mut,
                np.nan, np.nan, min_dist, max_dist, avg_dist,
                avg_major_cn, actual_perc_diploid, avg_tot_cn, 'ccube', "all", np.nan])
sml = list()
sml.append(loci_table.ccube_ccf_mean.nunique())
sml.append(np.nan)
sml.append(np.nan)
sml.append(score1B(sim_data_obj, loci_table))
sml.append(score1C(sim_data_obj, loci_table))
sml.append(score2A(sim_data_obj, loci_table))
auc, accuracy, sensitivity, specificity, precision = score2C(sim_data_obj, loci_table)
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
res_df.to_csv('{}/eval_ccube.tsv'.format(folder_path), sep='\t',
              index=False)

