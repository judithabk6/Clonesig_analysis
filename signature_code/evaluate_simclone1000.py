#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pkg_resources
import bz2
import pickle
from signature_code.evaluate_phylo500 import format_clonesig, format_sciclone, format_pyclone, format_ccube, format_dpclust, format_deconstructsigs, format_palimpsest, format_tracksig, format_tracksigfreq
from clonesig.estimator import Estimator
from pandas.errors import EmptyDataError
from clonesig.evaluate import score1B_base, score1C_base, score2A_base, score2C_base, score_sig_1A_base, score_sig_1B_base, score_sig_1C_base, score_sig_1D_base, score_sig_1E_base



MIXTURE_THRESHOLD = 0.05


"""
folder_path = 'SimClone1000/sim07gxi7_cst'
"""

signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
all_sigs = sig.columns[2:].to_list()

def format_truth(folder_path):
    samplename = folder_path.split('/')[1].split('_')[0]
    truth_path = 'data_simu_pcawg/testing_SimClone1000_truth'
    cluster_filename = '{}/{}/simulated_0001/truth_tree/simulated_0001_subclonal_structure.txt'.format(truth_path, samplename)
    cluster_df = pd.read_csv(cluster_filename, sep='\t')
    J_true = cluster_df.shape[0]
    phi_true_values = cluster_df.simulated_0001_0001.values
    cluster_df = cluster_df.assign(weights=cluster_df.n_snvs/cluster_df.n_snvs.sum())
    weights_true = cluster_df.weights.values
    mut_assign_filename = '{}/{}/simulated_0001/truth_tree/simulated_0001_mutation_assignments.txt'.format(truth_path, samplename)
    final_df = pd.read_csv(mut_assign_filename, sep='\t')
    final_df = final_df.assign(mutation_id=final_df.chr.astype(str) + '_'
                               + final_df.pos.astype(str))
    final_df = final_df.assign(true_cluster_id=final_df.cluster)
    est_clonal_idx = cluster_df[cluster_df.simulated_0001_0001==1].cluster.unique()[0]
    final_df = final_df.assign(true_subclonal=
        (final_df.true_cluster_id != est_clonal_idx).astype(int))

    input_filename = '{}/input_t.tsv'.format(folder_path)
    input_df = pd.read_csv(input_filename, sep='\t')
    sig_profile_1A = np.zeros(96)
    val, c = np.unique(input_df.trinucleotide, return_counts=True)
    sig_profile_1A[val.astype(int)] = c
    sig_profile_1A = sig_profile_1A / len(input_df)

    pi_mat = np.loadtxt('{}/pi_matrix.csv'.format(folder_path), delimiter=',')
    mu_filename = '{}/subMU.csv'.format(folder_path)
    mu = pd.read_csv(mu_filename, sep='\t')
    mu_mat = mu.values[:, 1:].T
    if len(pi_mat.shape)>1:
        sig_profile_1B = weights_true.dot(pi_mat).dot(mu_mat)
    else:
        sig_profile_1B = weights_true.dot(pi_mat.reshape(1, -1)).dot(mu_mat)
    true_signatures_1C = np.zeros(len(all_sigs))
    relevant_sigs = mu.columns[1:]
    relevant_sigs_idx = [all_sigs.index(s) for s in relevant_sigs]
    true_signatures_1C[np.array(relevant_sigs_idx)] = pi_mat.sum(axis=0).astype(bool).astype(int)
    true_signatures_1D = input_df.signature.astype(int).values
    if len(pi_mat.shape)>1:
        true_profile_1E = pi_mat[input_df.clone.values-1, :].dot(mu_mat)
    else:
        true_profile_1E = pi_mat.reshape(1, -1)[input_df.clone.values-1, :].dot(mu_mat)

    return (J_true, phi_true_values, weights_true,
            final_df[['mutation_id', 'true_cluster_id']],
            final_df[['mutation_id', 'true_subclonal']], sig_profile_1A,
            sig_profile_1B, true_signatures_1C, true_signatures_1D,
            true_profile_1E)


def format_phylogicndt(folder_path):
    folder_path = folder_path.replace('var', 'cst')
    with open('{}/phylogicndt_timing.txt'.format(folder_path), 'r') as f:
        line = f.read()
        start, end = float(line.split(',')[0]), float(line.split(',')[1])
    runtime = end - start
    try:
        loci_table = pd.read_csv(
            '{}/phylogicndt/Test_Clust.mut_ccfs.txt'.format(folder_path),
            sep='\t')
        loci_table = loci_table.assign(chr_num=loci_table.Chromosome.str.replace('chr', ''))
        loci_table = loci_table.assign(mutation_id_short=loci_table.chr_num.astype(str) + '_' + loci_table.Start_position.astype(str))
        cluster_table = pd.read_csv(
            '{}/phylogicndt/Test_Clust.cluster_ccfs.txt'.format(folder_path),
            sep='\t')
        cluster_table = pd.merge(cluster_table,
                                 loci_table.Cluster_Assignment.value_counts().to_frame(),
                                 left_on='Cluster_ID', right_index=True)
    except FileNotFoundError:
        return [None] * 10 + [runtime]
    J_pred = len(cluster_table)
    weights_pred = cluster_table['Cluster_Assignment'] / cluster_table['Cluster_Assignment'].sum()
    phi_pred_values = cluster_table['postDP_ccf_mean']
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    data_df = data_df.assign(mutation_id_short=data_df.chromosome.astype(str) + '_' + data_df.position.astype(str))
    data_df_m = pd.merge(data_df, loci_table[['mutation_id_short', 'Cluster_Assignment']], on="mutation_id_short")
    data_df_m = data_df_m.assign(pred_cluster_id=data_df_m.Cluster_Assignment)
    est_clonal_idx = cluster_table.sort_values(by='postDP_ccf_mean').iloc[-1]['Cluster_ID']
    data_df_m = data_df_m.assign(
        pred_subclonal=(data_df_m.pred_cluster_id != est_clonal_idx).astype(int))
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df_m[['mutation_id', 'pred_cluster_id']],
            data_df_m[['mutation_id', 'pred_subclonal']],
            None, None, None, None, runtime)


if __name__ == '__main__':
    folder_path = sys.argv[1]
    method_list = ['pyclone', 'sciclone', 'ccube', 'dpclust', 'phylogicndt',
                   'clonesig', 'deconstructsigs', 'palimpsest', 'tracksig',
                   'tracksigfreq']
    method_function_dict = {'pyclone': format_pyclone, 'sciclone': format_sciclone,
                            'ccube': format_ccube, 'dpclust': format_dpclust,
                            'phylogicndt': format_phylogicndt,
                            'clonesig': format_clonesig,
                            'deconstructsigs': format_deconstructsigs,
                            'palimpsest': format_palimpsest,
                            'tracksig': format_tracksig,
                            'tracksigfreq': format_tracksigfreq}

    (J_true, phi_true_values, weights_true, true_cluster_assign, true_subclonal,
     sig_profile_1A, sig_profile_1B, true_signatures_1C, true_signatures_1D,
     true_profile_1E) = format_truth(folder_path)
    info_df = pd.read_csv('{}/info_df.csv'.format(folder_path), sep='\t')
    df_list = list()
    df_cols = ['sample', 'cancer_loc', 'nb_mut', 'true_nb_clones', 'true_purity',
               'perc_dip', 'var_or_cst', 'median_depth', 'fitted_nb_clones',
               'll_ratio', 'pval', 'score1B',
               'score1C', 'score2A', 'score2C_auc', 'score2C_accuracy',
               'score2C_sensitivity', 'score2C_specificity', 'score2C_precision',
               'score_sig_1A', 'score_sig_1B', 'score_sig_1C_auc',
               'score_sig_1C_accuracy', 'score_sig_1C_sensitivity',
               'score_sig_1C_specificity', 'score_sig_1C_precision',
               'score_sig_1D', 'min_diff_distrib_mut', 'max_diff_distrib_mut',
               'std_diff_distrib_mut', 'median_diff_distrib_mut', 'perc_dist_5',
               'perc_dist_10', 'runtime', 'method']
    for method in method_list:
        print(method)
        row_list = list()
        row_list.append(folder_path.split('/')[1].split('_')[0])
        row_list.append(info_df.cancer_loc.unique()[0])
        row_list.append(len(true_cluster_assign))
        row_list.append(J_true)
        input_filename = '{}/input_t.tsv'.format(folder_path)
        input_df = pd.read_csv(input_filename, sep='\t')
        row_list.append(input_df.purity.mean())
        row_list.append(info_df.perc_dip.mean())
        row_list.append(folder_path.split('_')[-1])
        row_list.append(info_df.median_depth.mean())
        (ll_ratio, pval, J_pred, phi_pred_values, weights_pred,
         pred_cluster_assign, pred_subclonal, pred_profile, pred_signatures,
         pred_signatures_mut, est_dist, runtime) = \
            method_function_dict[method](folder_path)
        row_list.append(J_pred)
        row_list.append(ll_ratio)
        row_list.append(pval)
        if J_pred is not None:
            row_list.append(score1B_base(J_true, J_pred))
        else:
            row_list.append(np.nan)
        if phi_pred_values is not None:
            row_list.append(score1C_base(phi_true_values, phi_pred_values,
                                         weights_true, weights_pred))
        else:
            row_list.append(np.nan)
        if pred_cluster_assign is not None:
            ordered_table = pd.merge(pred_cluster_assign, true_cluster_assign,
                                     on='mutation_id', how='inner')
            if len(true_cluster_assign)<20000:
                row_list.append(score2A_base(ordered_table.true_cluster_id,
                                             ordered_table.pred_cluster_id))
            else:
                row_list.append(np.nan)
            ordered_table = pd.merge(pred_subclonal, true_subclonal,
                                     on='mutation_id', how='inner')
            auc, accuracy, sensitivity, specificity, precision = \
                score2C_base(ordered_table.true_subclonal,
                             ordered_table.pred_subclonal)
            for v in (auc, accuracy, sensitivity, specificity, precision):
                row_list.append(v)
        else:
            for i in range(6):
                row_list.append(np.nan)
        if pred_profile is not None:
            row_list.append(score_sig_1A_base(sig_profile_1A, pred_profile))
            row_list.append(score_sig_1B_base(sig_profile_1B, pred_profile))
            auc, accuracy, sensitivity, specificity, precision = \
                score_sig_1C_base(true_signatures_1C, pred_signatures)
            for v in (auc, accuracy, sensitivity, specificity, precision):
                row_list.append(v)
            if method == 'deconstructsigs':
                nb_rows = min(est_dist.shape[0], true_profile_1E.shape[0])
                score_sig_1D = score_sig_1D_base(true_signatures_1D[0:nb_rows],
                                                 pred_signatures_mut[0:nb_rows])
                (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
                 median_diff_distrib_mut, perc_dist_5, perc_dist_10) = \
                    score_sig_1E_base(true_profile_1E[0:nb_rows, :].astype(float),
                                      est_dist[0:nb_rows, :].astype(float))
            else:
                ok_ids = ordered_table.mutation_id.values
                true_filter = (true_subclonal.mutation_id.isin(ok_ids)).values
                pred_filter = (pred_subclonal.mutation_id.isin(ok_ids)).values
                score_sig_1D = score_sig_1D_base(true_signatures_1D[true_filter],
                                                 pred_signatures_mut[pred_filter])
                (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
                 median_diff_distrib_mut, perc_dist_5, perc_dist_10) = \
                    score_sig_1E_base(true_profile_1E[true_filter, :].astype(float),
                                      est_dist[pred_filter, :].astype(float))
            for v in (score_sig_1D, min_diff_distrib_mut, max_diff_distrib_mut,
                      std_diff_distrib_mut, median_diff_distrib_mut,
                      perc_dist_5, perc_dist_10):
                row_list.append(v)
        else:
            for i in range(14):
                row_list.append(np.nan)
        row_list.append(runtime)
        row_list.append(method)
        df_list.append(row_list)
    res_df = pd.DataFrame(df_list, columns=df_cols)
    res_df.to_csv('{}/result_evaluation_simclone1000.csv'.format(folder_path),
                  sep='\t', index=False)
