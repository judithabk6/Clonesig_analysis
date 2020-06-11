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
import bz2


MIXTURE_THRESHOLD = 0.05
folder_path = sys.argv[1]

"""
folder_path = 'salcedo_dream_challenge/T2_8X'
"""


def format_truth(folder_path):
    truth_path = 'data/salcedo_dream_challenge/MuTect_truth'
    tumor = folder_path.split('/')[1].split('_')[0]
    depth = folder_path.split('/')[1].split('_')[1]
    filename_1A = '{}/MuTect_{}.T.{}.truth.1A.txt'\
        .format(truth_path, tumor, depth)
    with open(filename_1A, 'r') as f:
        true_purity = float(f.read())
    filename_1B = '{}/MuTect_{}.T.{}.truth.1B.txt'\
        .format(truth_path, tumor, depth)
    with open(filename_1B, 'r') as f:
        J_true = int(f.readline().strip())
    filename_1C = '{}/MuTect_{}.T.{}.truth.1C.txt'\
        .format(truth_path, tumor, depth)
    data_1C = pd.read_csv(filename_1C, sep='\t', header=None,
                          names=['cluster_id', 'nb_mut', 'ccf'])
    valid_1C = data_1C[data_1C.ccf > 0]
    phi_true_values = valid_1C.ccf.values / true_purity
    weights_true = valid_1C.nb_mut.values / sum(valid_1C.nb_mut)
    filename_2A = '{}/MuTect_{}.T.{}.truth.2A.txt'\
        .format(truth_path, tumor, depth)
    data_2A = pd.read_csv(filename_2A, header=None, names=['cluster_id'])
    unfiltered_df = pd.read_csv('{}/unrestricted_input_mut.csv'.format(folder_path), sep='\t')
    truth_vcf = pd.read_csv(
        '{}/MuTect_{}.T.{}.truth.scoring_vcf.vcf'.format(truth_path, tumor, depth),
        sep='\t', comment='#', index_col=False, header=None,
        names=['chromosome', 'position', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
               'INFO', 'FORMAT', 'normal', 'tumor', 'calling'])
    unfiltered_df = unfiltered_df.assign(calling=truth_vcf.calling)
    final_df = unfiltered_df[unfiltered_df.calling]
    final_df = final_df.reset_index()
    final_df = final_df.assign(true_cluster_id=data_2A.cluster_id.astype(int))
    est_clonal_idx = valid_1C.ccf.values.argmax() + 1
    final_df = final_df.assign(
        true_subclonal=(final_df.true_cluster_id != est_clonal_idx).astype(int))

    nb_mut = valid_1C.nb_mut.sum()
    # nb_mut à changer plus tard
    sig_prop_filename = 'data/salcedo_dream_challenge/input_trinucleotide_signatures.txt'
    sig_prop_data = pd.read_csv(sig_prop_filename, sep='\t')
    sig_prop_data = sig_prop_data[sig_prop_data.tumour==tumor]
    sig_matrix_filename = 'data/salcedo_dream_challenge/signatures.txt'
    sig_mat = pd.read_csv(sig_matrix_filename, sep='\t')
    relevant_sigs = sig_prop_data.signature.values
    sig_profile_1A = np.zeros(96)
    val, c = np.unique(unfiltered_df[unfiltered_df.calling].trinucleotide,
                       return_counts=True)
    sig_profile_1A[val.astype(int)] = c
    sig_profile_1A = sig_profile_1A / len(unfiltered_df[unfiltered_df.calling])
    sig_profile_1B = sig_prop_data.frequency.values.reshape(1,-1).dot(
        sig_mat[relevant_sigs].values.T)[0]
    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    true_signatures = np.isin(complete_mat.columns[1:], relevant_sigs).astype(int)
    true_profile_1E = np.repeat([sig_profile_1B], nb_mut, axis=0)
    return (J_true, phi_true_values, weights_true,
            final_df[['mutation_id', 'true_cluster_id']],
            final_df[['mutation_id', 'true_subclonal']], sig_profile_1A,
            sig_profile_1B, true_signatures, true_profile_1E)


def format_clonesig(folder_path, setting):
    try:
        raw_res = bz2.BZ2File('{}/{}_clonesig_raw_results.bz2'.format(folder_path, setting), 'rb')
        new_est, lr, pval, cst_est, fitted_sigs, runtime = pickle.load(raw_res)
    except FileNotFoundError:
        return [None] * 11
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
    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    complete_sigs = complete_mat.columns[1:]
    pred_signatures = np.zeros(len(complete_sigs))
    est_sig = new_est.xi.dot(new_est.pi)
    if setting in ('all', 'all_nuclonal'):
        pred_signatures = est_sig
    elif setting == 'cancertype':
        filename = '{}/sub_MU_matrix.csv'.format(folder_path)
        sub_matrix = pd.read_csv(filename, sep='\t')
        sub_sigs = sub_matrix.columns[1:]
        idx = [list(complete_sigs).index(s) for s in sub_sigs]
        pred_signatures[np.array(idx)] = est_sig
    elif setting == 'prefit':
        pred_signatures[fitted_sigs] = est_sig
    else:
        raise NameError('unknown setting for CloneSig')
    est_dist = new_est.pi[new_est.qun.argmax(axis=1), :].dot(new_est.mu_matrix)
    return (ll_ratio, pval, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']],
            pred_profile, pred_signatures, est_dist, runtime)


def format_pyclone(folder_path, setting):
    try:
        with open('{}/pyclone_timing.txt'.format(folder_path), 'r') as f:
            line = f.read()
            start, end = float(line.split(',')[0]), float(line.split(',')[1])
        runtime = end - start
    except:
        runtime = np.nan
    try:
        pyclone_res = '{}/pyclone/tables'.format(folder_path)
        cluster_table = pd.read_csv('{}/cluster.tsv'.format(pyclone_res),
                                    sep='\t')
        loci_table = pd.read_csv('{}/loci.tsv'.format(pyclone_res), sep='\t')
    except FileNotFoundError:
        return [None] * 10 + [runtime]
    J_pred = len(cluster_table[cluster_table['size'] > 1])
    weights_pred = cluster_table['size'] / cluster_table['size'].sum()
    phi_pred_values = cluster_table['mean']
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    ordered_loci_table = pd.merge(data_df, loci_table, on='mutation_id')
    ordered_loci_table = ordered_loci_table.assign(
        pred_cluster_id=ordered_loci_table.cluster_id)
    est_clonal_idx = cluster_table.sort_values(by='mean').iloc[-1].cluster_id
    ordered_loci_table = ordered_loci_table.assign(
        pred_subclonal=(ordered_loci_table.cluster_id != est_clonal_idx)
        .astype(int))
    runtime = end - start
    return (None, None, J_pred, phi_pred_values, weights_pred,
            ordered_loci_table[['mutation_id', 'pred_cluster_id']],
            ordered_loci_table[['mutation_id', 'pred_subclonal']],
            None, None, None, runtime)


def format_sciclone(folder_path, setting):
    with open('{}/sciclone_timing.txt'.format(folder_path), 'r') as f:
        line = f.read()
        start, end = float(line.split(',')[0]), float(line.split(',')[1])
    runtime = end - start
    try:
        loci_table = pd.read_csv('{}/sciclone/clusters1'.format(folder_path),
                                 sep='\t')
    except FileNotFoundError:
        return [None] * 10 + [runtime]
    J_pred = loci_table[loci_table.cluster > 0].cluster.nunique()
    weights_pred = loci_table[loci_table.cluster > 0].groupby('cluster')['tumor.vaf'].count()/len(loci_table)
    phi_pred_values = loci_table[loci_table.cluster > 0].groupby('cluster')['tumor.vaf'].mean()/100
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    data_df = data_df.assign(pred_cluster_id=loci_table.cluster)
    est_clonal_idx = (loci_table[loci_table.cluster > 0].groupby('cluster')['tumor.vaf'].mean()).idxmax()
    data_df = data_df.assign(
        pred_subclonal=(data_df.pred_cluster_id != est_clonal_idx).astype(int))
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']],
            None, None, None, runtime)


def format_dpclust(folder_path, setting):
    with open('{}/dpclust_timing.txt'.format(folder_path), 'r') as f:
        line = f.read()
        start, end = float(line.split(',')[0]), float(line.split(',')[1])
    runtime = end - start
    try:
        res_folder = '{}_DPoutput_2000iters_1000burnin_seed123'.format(folder_path.split('/')[1])
        loci_table = pd.read_csv('{}/dpclust/{}/{}_2000iters_1000burnin_bestConsensusAssignments.bed'.format(folder_path, res_folder, folder_path.split('/')[1]), sep='\t')
        cluster_table = pd.read_csv('{}/dpclust/{}/{}_2000iters_1000burnin_bestClusterInfo.txt'.format(folder_path, res_folder, folder_path.split('/')[1]), sep='\t')
    except FileNotFoundError:
        return [None] * 10 + [runtime]
    J_pred = len(cluster_table)
    weights_pred = cluster_table['no.of.mutations'] / cluster_table['no.of.mutations'].sum()
    phi_pred_values = cluster_table['location']
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    data_df = data_df.assign(pred_cluster_id=loci_table.cluster)
    est_clonal_idx = cluster_table.sort_values(by='location').iloc[-1]['cluster.no']
    data_df = data_df.assign(
        pred_subclonal=(data_df.pred_cluster_id != est_clonal_idx).astype(int))
    runtime = end - start
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']],
            None, None, None, runtime)


def format_phylogicndt(folder_path, setting):
    with open('{}/phylogicndt_timing.txt'.format(folder_path), 'r') as f:
        line = f.read()
        start, end = float(line.split(',')[0]), float(line.split(',')[1])
    runtime = end - start
    try:
        loci_table = pd.read_csv(
            '{}/phylogicndt/Test_Clust.mut_ccfs.txt'.format(folder_path),
            sep='\t')
        loci_table = loci_table.assign(chr_num=loci_table.Chromosome.str.replace('chr', '').astype(int))
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
            None, None, None, runtime)


def format_ccube(folder_path, setting):
    with open('{}/ccube_timing.txt'.format(folder_path), 'r') as f:
        line = f.read()
        start, end = float(line.split(',')[0]), float(line.split(',')[1])
    runtime = end - start
    try:
        loci_table = pd.read_csv(
            '{}/ccube/ssm_clusters.csv'.format(folder_path), sep='\t')
    except FileNotFoundError:
        return [None] * 10 + [runtime]
    pre_cluster_table = loci_table.groupby('ccube_ccf_mean').rough_mult.count()
    loci_table = loci_table.assign(
        cluster_id=loci_table.apply(
            lambda x:
            pre_cluster_table.index.tolist().index(x['ccube_ccf_mean']),
            axis=1))
    J_pred = loci_table.ccube_ccf_mean.nunique()
    weights_pred = pre_cluster_table.values/len(loci_table)
    phi_pred_values = pre_cluster_table.index.values
    data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
    data_df = data_df.assign(pred_cluster_id=loci_table.cluster_id)
    est_clonal_idx = np.argmax(pre_cluster_table.index.tolist())
    data_df = data_df.assign(
        pred_subclonal=(data_df.pred_cluster_id != est_clonal_idx).astype(int))
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']],
            None, None, None, runtime)


def format_deconstructsigs(folder_path, setting):
    res_filename = '{}/deconstructsigs/signatures_{}.csv'\
        .format(folder_path, setting)
    result_file = pd.read_csv(res_filename, sep=' ')
    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    complete_sigs = complete_mat.columns[1:]
    pred_signatures = np.zeros(len(complete_sigs))
    if setting == 'all':
        mu_mat_setting = complete_mat[complete_mat.columns[1:]].values.T
        pred_signatures = result_file.values[0]
    elif setting == 'cancertype':
        filename = '{}/sub_MU_matrix.csv'.format(folder_path)
        sub_matrix = pd.read_csv(filename, sep='\t')
        mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
        sub_sigs = sub_matrix.columns[1:]
        idx = [list(complete_sigs).index(s) for s in sub_sigs]
        pred_signatures[np.array(idx)] = result_file.iloc[0].values
    else:
        raise NameError('unknown setting for DeconstructSigs')
    sig_profile = result_file.values.dot(mu_mat_setting)
    input_filename = '{}/deconstructsigs/pattern96.csv'.format(folder_path)
    pattern = pd.read_csv(input_filename, sep='\t')
    nb_mut = pattern.sum().sum()
    pred_profile_1E = np.repeat([sig_profile], nb_mut, axis=0)
    runtime = pd.read_csv('{}/deconstructsigs/deconstructsig_runtime_{}.csv'
                          .format(folder_path, setting),
                          index_col=0).values[0][0]
    return (None, None, None, None, None, None, None,
            sig_profile, pred_signatures, pred_profile_1E, runtime)


def format_palimpsest(folder_path, setting):
    try:
        mixture_file = pd.read_csv('{}/palimpsest/palimpsest_mixtures_{}.csv'.
                                   format(folder_path, setting), sep='\t')
        ccf_file = pd.read_csv('{}/palimpsest/palimpsest_mut_data_{}.csv'
                               .format(folder_path, setting), sep='\t')
    except FileNotFoundError:
        return [None] * 11
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
    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    complete_sigs = complete_mat.columns[1:]
    pred_signatures = np.zeros(len(complete_sigs))
    if setting == 'all':
        mu_mat_setting = complete_mat[complete_mat.columns[1:]].values.T
        idx = np.where(pred_signatures==0)[0]
    elif setting == 'cancertype':
        filename = '{}/sub_MU_matrix.csv'.format(folder_path)
        sub_matrix = pd.read_csv(filename, sep='\t')
        mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
        sub_sigs = sub_matrix.columns[1:]
        idx = [list(complete_sigs).index(s) for s in sub_sigs]
    elif setting == 'prefit':
        premixture_file = pd.read_csv(
            '{}/palimpsest/palimpsest_premixtures_{}.txt'.
            format(folder_path, setting), sep=' ')
        sig_names = premixture_file.columns.to_list()
        sig_names = [s.replace('.', ' ') for s in sig_names]
        sig_names = [s.replace('PCAWG ', 'PCAWG-') for s in sig_names]
        idx = [list(complete_sigs).index(s) for s in sig_names]
        mu_mat_setting = complete_mat[sig_names].values.T
    else:
        raise NameError('unknown setting for Palimpsest')
    pred_profile = (ccf_file.groupby('clonality_binary').CCF.count() /
                   len(ccf_file)).values.reshape(1, -1)\
        .dot(mixture_file).dot(mu_mat_setting)
    est_sigs = (ccf_file.groupby('clonality_binary').CCF.count() /
            len(ccf_file)).values.reshape(1, -1) \
        .dot(mixture_file)
    pred_signatures[np.array(idx)] = est_sigs[0]
    est_dist = mixture_file.values[ccf_file.clonality_binary].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/palimpsest/palimpsest_runtime_{}.csv'.format(folder_path, setting),
                           index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred,
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']], pred_profile,
            pred_signatures, est_dist, runtime)


def format_tracksig(folder_path, setting):
    try:
        mixture_file = pd.read_csv('{}/tracksig/tracksig_mixtures_{}.csv'.
                                   format(folder_path, setting), sep=',')
    except FileNotFoundError:
        return [None] * 11
    try:
        changepoint_file = pd.read_csv(
            '{}/tracksig/tracksig_changepoints_{}.txt'.
            format(folder_path, setting), header=None, sep=' ')
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


    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    complete_sigs = complete_mat.columns[1:]
    pred_signatures = np.zeros(len(complete_sigs))

    if setting == 'all':
        mu_mat_setting = complete_mat[complete_mat.columns[1:]].values.T
        idx = np.where(pred_signatures==0)[0]
    elif setting == 'cancertype':
        filename = '{}/sub_MU_matrix.csv'.format(folder_path)
        sub_matrix = pd.read_csv(filename, sep='\t')
        mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
        sub_sigs = sub_matrix.columns[1:]
        idx = [list(complete_sigs).index(s) for s in sub_sigs]
    elif setting == 'prefit':
        sig_names = mixture_file[mixture_file.columns[0]].values
        sig_names = [s.replace('.', ' ').replace('PCAWG ', 'PCAWG-') for s in sig_names]
        idx = [list(complete_sigs).index(s) for s in sig_names]
        mu_mat_setting = complete_mat[sig_names].values.T
    else:
        raise NameError('unknown setting for TrackSig')

    est_sigs = mixture_file[mixture_file.columns[1:]].mean(axis=1).values
    pred_signatures[idx] = est_sigs

    pred_profile = est_sigs.dot(mu_mat_setting)

    est_dist = mixture_file.values[:, 1:].T[input_df_filter.mutation_group.astype(int)].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/tracksig/tracksig_runtime_{}.csv'.format(folder_path, setting),
                          index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred, 
            input_df_filter[['mutation_id', 'pred_cluster_id']],
            input_df_filter[['mutation_id', 'pred_subclonal']], pred_profile,
            pred_signatures, est_dist, runtime)

def format_tracksigfreq(folder_path, setting):
    try:
        mixture_file = pd.read_csv('{}/tracksigfreq/tracksigfreq_mixtures_{}.csv'.
                                   format(folder_path, setting), sep=',')
    except FileNotFoundError:
        return [None] * 11
    try:
        changepoint_file = pd.read_csv(
            '{}/tracksigfreq/tracksigfreq_changepoints_{}.txt'.
            format(folder_path, setting), header=None, sep=' ')
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

    complete_mat_filename = '{}/MU_matrix.csv'.format(folder_path)
    complete_mat = pd.read_csv(complete_mat_filename, sep='\t')
    complete_sigs = complete_mat.columns[1:]
    pred_signatures = np.zeros(len(complete_sigs))

    if setting == 'all':
        mu_mat_setting = complete_mat[complete_mat.columns[1:]].values.T
        idx = np.where(pred_signatures==0)[0]
    elif setting == 'cancertype':
        filename = '{}/sub_MU_matrix.csv'.format(folder_path)
        sub_matrix = pd.read_csv(filename, sep='\t')
        mu_mat_setting = sub_matrix[sub_matrix.columns[1:]].values.T
        sub_sigs = sub_matrix.columns[1:]
        idx = [list(complete_sigs).index(s) for s in sub_sigs]
    elif setting == 'prefit':
        sig_names = mixture_file[mixture_file.columns[0]].values
        sig_names = [s.replace('.', ' ').replace('PCAWG ', 'PCAWG-') for s in sig_names]
        idx = [list(complete_sigs).index(s) for s in sig_names]
        mu_mat_setting = complete_mat[sig_names].values.T
    else:
        raise NameError('unknown setting for TrackSigFreq')
    est_sigs = mixture_file[mixture_file.columns[1:]].mean(axis=1).values
    pred_signatures[idx] = est_sigs

    pred_profile = est_sigs.dot(mu_mat_setting)

    est_dist = mixture_file.values[:, 1:].T \
        [data_df.bin.astype(int)-1].dot(mu_mat_setting)
    runtime = pd.read_csv('{}/tracksigfreq/tracksigfreq_runtime_{}.csv'.format(folder_path, setting),
                           index_col=0).values[0][0]
    return (None, None, J_pred, phi_pred_values, weights_pred, 
            data_df[['mutation_id', 'pred_cluster_id']],
            data_df[['mutation_id', 'pred_subclonal']], pred_profile, 
            pred_signatures, est_dist, runtime)


method_setting_list = [('pyclone', None), ('sciclone', None), ('ccube', None),
                       ('dpclust', None), ('phylogicndt', None),
                       ('clonesig', 'all'), ('clonesig', 'cancertype'),
                       ('clonesig', 'prefit'), ('deconstructsigs', 'all'),
                       ('deconstructsigs', 'cancertype'),
                       ('palimpsest', 'all'), ('palimpsest', 'cancertype'),
                       ('palimpsest', 'prefit'), ('tracksig', 'all'),
                       ('tracksig', 'cancertype'), ('tracksig', 'prefit'),
                       ('tracksigfreq', 'all'), ('tracksigfreq', 'cancertype'),
                       ('tracksigfreq', 'prefit')]
method_function_dict = {'pyclone': format_pyclone, 'sciclone': format_sciclone,
                        'ccube': format_ccube, 'dpclust': format_dpclust,
                        'phylogicndt': format_phylogicndt,
                        'clonesig': format_clonesig,
                        'deconstructsigs': format_deconstructsigs,
                        'palimpsest': format_palimpsest,
                        'tracksig': format_tracksig,
                        'tracksigfreq': format_tracksigfreq}


(J_true, phi_true_values, weights_true, true_cluster_assign, true_subclonal,
 sig_profile_1A, sig_profile_1B, true_signatures, true_profile_1E) = \
    format_truth(folder_path)

tumor = folder_path.split('/')[1].split('_')[0]
depth = folder_path.split('/')[1].split('_')[1]

f = open('data/salcedo_dream_challenge/MuTect_truth/MuTect_{}.T.{}.truth.1A.txt'.format(tumor, depth), 'r')
true_purity = float(f.readline().strip())
perc_dip = 0 # jfkdjf

cnv_filename = 'data/salcedo_dream_challenge/MuTect_inputs/{}-{}_refit_subclones_noXY.txt'.format(tumor, depth)
cnv_table = pd.read_csv(cnv_filename, sep='\t')

def get_major_minor(x):
    if x.frac1_A==1:
        return pd.Series([x.nMaj1_A, x.nMin1_A])
    else:
        if x.frac1_A > x.frac2_A:
            return pd.Series([x.nMaj1_A, x.nMin1_A])
        else:
            return pd.Series([x.nMaj2_A, x.nMin2_A])

new_data = cnv_table.apply(get_major_minor, axis=1)
new_data.columns = ['major_cn', 'minor_cn']
cnv_final = pd.concat((cnv_table[['chr', 'startpos', 'endpos']],
                       new_data.astype(int)), axis=1)
cnv_final = cnv_final.assign(weight=cnv_final.endpos - cnv_final.startpos)
perc_dip = cnv_final[(cnv_final.major_cn==1)&(cnv_final.minor_cn==1)].weight.sum() / cnv_final.weight.sum()

df_list = list()
df_cols = ['tumor', 'depth', 'nb_mut', 'true_nb_clones', 'true_purity',
           'perc_dip', 'fitted_nb_clones', 'll_ratio', 'pval', 'score1B',
           'score1C', 'score2A', 'score2C_auc', 'score2C_accuracy',
           'score2C_sensitivity', 'score2C_specificity', 'score2C_precision',
           'score_sig_1A', 'score_sig_1B', 'score_sig_1C_auc',
           'score_sig_1C_accuracy', 'score_sig_1C_sensitivity',
           'score_sig_1C_specificity', 'score_sig_1C_precision',
           'min_diff_distrib_mut', 'max_diff_distrib_mut',
           'std_diff_distrib_mut', 'median_diff_distrib_mut', 'perc_dist_5',
           'perc_dist_10', 'runtime', 'method', 'setting']
for method, setting in method_setting_list:
    print(method, setting)
    row_list = list()
    row_list.append(tumor)
    row_list.append(depth)
    row_list.append(len(true_cluster_assign))
    row_list.append(J_true)
    row_list.append(true_purity)
    row_list.append(perc_dip)
    (ll_ratio, pval, J_pred, phi_pred_values, weights_pred,
     pred_cluster_assign, pred_subclonal, pred_profile, pred_signatures,
     est_dist, runtime) = \
        method_function_dict[method](folder_path, setting)
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
            score_sig_1C_base(true_signatures, pred_signatures)
        for v in (auc, accuracy, sensitivity, specificity, precision):
            row_list.append(v)
        if method == 'deconstructsigs':
            nb_rows = min(est_dist.shape[0], true_profile_1E.shape[0])
            (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
             median_diff_distrib_mut, perc_dist_5, perc_dist_10) = \
                score_sig_1E_base(true_profile_1E[0:nb_rows, :],
                                  est_dist[0:nb_rows, :])
        else:
            ok_ids = ordered_table.mutation_id.values
            true_filter = (true_subclonal.mutation_id.isin(ok_ids)).values
            pred_filter = (pred_subclonal.mutation_id.isin(ok_ids)).values
            (min_diff_distrib_mut, max_diff_distrib_mut, std_diff_distrib_mut,
             median_diff_distrib_mut, perc_dist_5, perc_dist_10) = \
                score_sig_1E_base(true_profile_1E[true_filter, :],
                                  est_dist[pred_filter, :].astype(float))
        for v in (min_diff_distrib_mut, max_diff_distrib_mut,
                  std_diff_distrib_mut, median_diff_distrib_mut,
                  perc_dist_5, perc_dist_10):
            row_list.append(v)
    else:
        for i in range(13):
            row_list.append(np.nan)
    row_list.append(runtime)
    row_list.append(method)
    row_list.append(setting)
    df_list.append(row_list)
res_df = pd.DataFrame(df_list, columns=df_cols)
res_df.to_csv('{}/result_evaluation_dream_new.csv'.format(folder_path),
              sep='\t', index=False)





