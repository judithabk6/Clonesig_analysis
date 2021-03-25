#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import sys
from collections import Iterable
import numpy as np
import pkg_resources
import bz2
import pickle
from signature_code.evaluate_dream import format_sciclone as format_sciclone_old
from signature_code.evaluate_dream import format_ccube as format_ccube_old
from signature_code.evaluate_dream import format_pyclone as format_pyclone_old
from signature_code.evaluate_dream import format_dpclust as format_dpclust_old
from signature_code.evaluate_dream import format_phylogicndt as format_phylogicndt_old
from clonesig.estimator import Estimator
from pandas.errors import EmptyDataError
from clonesig.evaluate import score1B_base, score1C_base, score2A_base, score2C_base, score_sig_1A_base, score_sig_1B_base, score_sig_1C_base, score_sig_1D_base, score_sig_1E_base



MIXTURE_THRESHOLD = 0.05


"""
folder_path = 'PhylogicNDT500/Sim_500_19_var'
"""

signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
all_sigs = sig.columns[2:].to_list()

def format_truth(folder_path):
    num_tumor = folder_path.split('_')[2]
    truth_path = 'data_simu_pcawg/training_PhylogicNDT500_simulations_12_16'
    num_cluster_filename = '{}/Num_Clust/Sim_500_{}.number_of_clusters.txt'.\
        format(truth_path, num_tumor)
    num_cluster = pd.read_csv(num_cluster_filename, sep='\t')
    J_true = num_cluster.clusters.unique()[0]
    structure_filename = '{}/Subclonal_Structure/Sim_500_{}.subclonal_structure.txt'.\
        format(truth_path, num_tumor)
    structure_df = pd.read_csv(structure_filename, sep='\t')
    structure_df = structure_df.assign(weights=structure_df.n_ssms/structure_df.n_ssms.sum())
    phi_true_values = structure_df.ccf.values
    weights_true = structure_df.weights.values
    mut_assign_filename = '{}/Mut_Assign/Sim_500_{}.mutation_assignments.txt'.\
        format(truth_path, num_tumor)
    final_df = pd.read_csv(mut_assign_filename, sep='\t')
    final_df = final_df.assign(mutation_id=final_df.chr.astype(str) + '_'
                               + final_df.pos.astype(str))
    final_df = final_df.assign(true_cluster_id=final_df.cluster)
    est_clonal_idx = structure_df[structure_df.ccf==1].cluster.unique()[0]
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
        true_profile_1E = pi_mat[input_df.clone.values, :].dot(mu_mat)
    else:
        true_profile_1E = pi_mat.reshape(1, -1)[input_df.clone.values, :].dot(mu_mat)

    return (J_true, phi_true_values, weights_true,
            final_df[['mutation_id', 'true_cluster_id']],
            final_df[['mutation_id', 'true_subclonal']], sig_profile_1A,
            sig_profile_1B, true_signatures_1C, true_signatures_1D,
            true_profile_1E)


def format_clonesig(folder_path):
    try:
        raw_res = bz2.BZ2File('{}/cancertype_clonesig_raw_results.bz2'.format(folder_path), 'rb')
        new_est, ll_ratio, pval, cst_est, fitted_sigs, runtime = pickle.load(raw_res)
    except FileNotFoundError:
        return [None] * 12
    J_pred = new_est.J
    phi_pred_values = new_est.phi
    pre_est_w = np.zeros(new_est.J)
    pre_counts = np.unique(np.argmax(new_est.qun, axis=1),
                           return_counts=True)
    pre_est_w[pre_counts[0]] = pre_counts[1]
    weights_pred = pre_est_w/new_est.N
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    pred_cluster_assign = np.argmax(new_est.qun, axis=1)
    data_df = data_df.assign(pred_cluster_id=np.argmax(new_est.qun, axis=1))
    est_att = np.argmax(new_est.qun, axis=1)
    est_clonal_idx = np.argmax(new_est.phi)
    data_df = data_df.assign(
        pred_subclonal=(data_df.pred_cluster_id != est_clonal_idx).astype(int))

    pred_profile = new_est.xi.dot(new_est.pi).dot(new_est.mu_matrix)
    pred_signatures = np.zeros(len(all_sigs))
    est_sig = new_est.xi.dot(new_est.pi)
    mu_filename = '{}/subMU.csv'.format(folder_path)
    mu = pd.read_csv(mu_filename, sep='\t')
    mu_mat = mu.values[:, 1:].T
    sub_sigs = mu.columns[1:]
    idx = [list(all_sigs).index(s) for s in sub_sigs]
    pred_signatures[np.array(idx)] = est_sig
    est_sig = new_est.rnus[np.arange(new_est.N), new_est.qun.argmax(axis=1), :].argmax(axis=1)

    est_dist = new_est.pi[new_est.qun.argmax(axis=1), :].dot(new_est.mu_matrix)
    return (ll_ratio, pval, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']],
            pred_profile, pred_signatures, est_sig, est_dist, runtime)


def format_sciclone(folder_path):
    tmp_res = format_sciclone_old(folder_path.replace('var', 'cst'), None)
    new_res = list(tmp_res[:-1]) + [None, tmp_res[-1]]
    return new_res


def format_pyclone(folder_path):
    tmp_res = format_pyclone_old(folder_path.replace('var', 'cst'), None)
    new_res = list(tmp_res[:-1]) + [None, tmp_res[-1]]
    return new_res


def format_ccube(folder_path):
    tmp_res = format_ccube_old(folder_path.replace('var', 'cst'), None)
    new_res = list(tmp_res[:-1]) + [None, tmp_res[-1]]
    return new_res


def format_dpclust(folder_path):
    tmp_res = format_dpclust_old(folder_path.replace('var', 'cst'), None)
    new_res = list(tmp_res[:-1]) + [None, tmp_res[-1]]
    return new_res


def format_phylogicndt(folder_path):
    tmp_res = format_phylogicndt_old(folder_path.replace('var', 'cst'), None)
    new_res = list(tmp_res[:-1]) + [None, tmp_res[-1]]
    return new_res


def format_deconstructsigs(folder_path):
    res_filename = '{}/deconstructsigs/signatures_cancertype.csv'\
        .format(folder_path)
    result_file = pd.read_csv(res_filename, sep=' ')
    pred_signatures = np.zeros(len(all_sigs))
    filename = '{}/subMU.csv'.format(folder_path)
    sub_matrix = pd.read_csv(filename, sep='\t')
    mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
    sub_sigs = sub_matrix.columns[1:]
    idx = [list(all_sigs).index(s) for s in sub_sigs]
    pred_signatures[np.array(idx)] = result_file.iloc[0].values
    sig_profile = result_file.values.dot(mu_mat_setting)
    input_filename = '{}/deconstructsigs/pattern96.csv'.format(folder_path)
    pattern = pd.read_csv(input_filename, sep='\t')

    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    est = Estimator(data_df.trinucleotide.values, data_df.var_counts.values,
                    data_df.normal_cn.values,
                    data_df.minor_cn.values + data_df.major_cn.values,
                    data_df.minor_cn.values,
                    data_df.var_counts.values + data_df.ref_counts.values,
                    data_df.purity.mean(), 1,
                    inputMU=mu_mat_setting, pi=result_file.values.reshape(1, -1))
    est_sig_att = est.rnus[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1)


    nb_mut = pattern.sum().sum()
    pred_profile_1E = np.repeat([sig_profile], nb_mut, axis=0)
    runtime = pd.read_csv('{}/deconstructsigs/deconstructsig_runtime_cancertype.csv'
                          .format(folder_path),
                          index_col=0).values[0][0]
    return (None, None, None, None, None, None, None,
            sig_profile, pred_signatures, est_sig_att, pred_profile_1E, runtime)


def format_palimpsest(folder_path):
    try:
        mixture_file = pd.read_csv('{}/palimpsest/palimpsest_mixtures_cancertype.csv'.
                                   format(folder_path), sep='\t')
        ccf_file = pd.read_csv('{}/palimpsest/palimpsest_mut_data_cancertype.csv'
                               .format(folder_path), sep='\t')
    except FileNotFoundError:
        return [None] * 12
    J_pred = 2
    weights_pred = ccf_file.groupby('Clonality').CCF.count().values/len(ccf_file)
    phi_pred_values = ccf_file.groupby('Clonality').CCF.mean().values
    ccf_file = ccf_file.assign(
        clonality_binary=ccf_file.apply(
            lambda x: 1 if x['Clonality'] == 'subclonal' else 0, axis=1))
    ccf_file = ccf_file.reset_index()
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    data_df = data_df.assign(pred_cluster_id=ccf_file.clonality_binary)
    data_df = data_df.assign(pred_subclonal=ccf_file.clonality_binary)
    pred_signatures = np.zeros(len(all_sigs))
    filename = '{}/subMU.csv'.format(folder_path)
    sub_matrix = pd.read_csv(filename, sep='\t')
    mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
    sub_sigs = sub_matrix.columns[1:]
    idx = [list(all_sigs).index(s) for s in sub_sigs]
    pred_profile = (ccf_file.groupby('clonality_binary').CCF.count() /
                   len(ccf_file)).values.reshape(1, -1)\
        .dot(mixture_file).dot(mu_mat_setting)
    est_sigs = (ccf_file.groupby('clonality_binary').CCF.count() /
            len(ccf_file)).values.reshape(1, -1) \
        .dot(mixture_file)
    pred_signatures[np.array(idx)] = est_sigs[0]
    relevant_sigs = sub_matrix.columns[1:].to_list()
    signature_mut = np.array([relevant_sigs.index(s) for s in ccf_file.origin])
    est_dist = mixture_file.values[ccf_file.clonality_binary].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/palimpsest/palimpsest_runtime_cancertype.csv'.format(folder_path),
                           index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']], pred_profile,
            pred_signatures, signature_mut, est_dist, runtime)


def format_tracksig(folder_path):
    try:
        mixture_file = pd.read_csv('{}/tracksig/tracksig_mixtures_cancertype.csv'.
                                   format(folder_path), sep=',')
    except FileNotFoundError:
        return [None] * 12
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

    J_pred = len(changepoints_tracksig_list) + 1
    weights_pred = input_df_filter.groupby('pred_cluster_id').vaf_purity.count().values/len(input_df_filter)
    phi_pred_values = input_df_filter.groupby('pred_cluster_id').vaf_purity.mean().values
    est_clonal_idx = input_df_filter.groupby('pred_cluster_id').vaf_purity.mean().idxmax()
    input_df_filter = input_df_filter.assign(
        pred_subclonal=(input_df_filter.pred_cluster_id != est_clonal_idx).astype(int))


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

    est_dist = mixture_file.values[:, 1:].T[input_df_filter.mutation_group.astype(int)].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/tracksig/tracksig_runtime_cancertype.csv'.format(folder_path),
                          index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred, 
            input_df_filter[['mutation_id', 'pred_cluster_id']],
            input_df_filter[['mutation_id', 'pred_subclonal']], pred_profile,
            pred_signatures, signature_mut, est_dist, runtime)


def format_tracksigfreq(folder_path):
    try:
        mixture_file = pd.read_csv('{}/tracksigfreq/tracksigfreq_mixtures_cancertype.csv'.
                                   format(folder_path), sep=',')
    except FileNotFoundError:
        return [None] * 12
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

    est_dist = mixture_file.values[:, 1:].T \
        [data_df.bin.astype(int)-1].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/tracksigfreq/tracksigfreq_runtime_cancertype.csv'.format(folder_path),
                           index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred, 
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']], pred_profile, 
            pred_signatures, signature_mut, est_dist, runtime)


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
        row_list.append(folder_path.split('_')[2])
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
    res_df.to_csv('{}/result_evaluation_phylo500.csv'.format(folder_path),
                  sep='\t', index=False)


