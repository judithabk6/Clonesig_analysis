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
from clonesig.estimator import Estimator, log_binomial_coeff, _beta_binomial_logpmf, beta_binomial_pmf
import time
import scipy as sp
from scipy.spatial.distance import cosine, pdist
import itertools
from sklearn import linear_model

folder_path = sys.argv[1]


'''
folder_path = '20190518_simulations_clonesig_cn_cancer_type/type0-perc_diploid20-nb_clones4-nb_mut100-true_c_mut100'
'''

cancer_type = int(folder_path.split('/')[1].split('type')[1].split('-')[0])
perc_diploid = int(folder_path.split('/')[1].split('perc_diploid')[1].split('-')[0])
nb_clones = int(folder_path.split('/')[1].split('nb_clones')[1].split('-')[0])
nb_mut = int(folder_path.split('/')[1].split('nb_mut')[1].split('-')[0])

sig_file_path = 'external_data/sigProfiler_SBS_signatures_2018_03_28.csv'
cancer_type_sig_filename = 'external_data/match_cancer_type_sig_v3.csv'

# open the matrix describing the signatures
SIG = pd.read_csv(sig_file_path)
SIG_MATRIX = SIG.values[:, 2:].astype(float).T
L, K = SIG_MATRIX.shape
NEW_SIG_MATRIX = SIG_MATRIX + 10**-20 * (SIG_MATRIX == 0)
MU = NEW_SIG_MATRIX / NEW_SIG_MATRIX.sum(axis=1)[:, np.newaxis]

# get the signatures specific of the cancer type
cancer_type_sig = pd.read_csv(cancer_type_sig_filename, index_col=0).values
select = cancer_type_sig[cancer_type, :]
subMU = MU[select.astype(bool), :]


def fit_model_special(T, B, C_normal, C_tumor_tot, C_tumor_minor, D, purity,
                      inputMU, nb_fits=1, seeds=None, max_nb_clones=6, extra=4):
    """
    possible metrics : F, loglikelihood, BIC, AIC, AICc, ICL_q, ICL_qn, SH.
    """
    L = inputMU.shape[0]
    if isinstance(seeds, Iterable):
        if len(seeds) != nb_fits:
            raise ValueError("Number of seeds is incompatible with number of required fits")
    if seeds is None:
        seeds = list(range(nb_fits))
    Fres = np.zeros((nb_fits, max_nb_clones+extra))
    bic = np.zeros((nb_fits, max_nb_clones+extra))
    loglikelihood = np.zeros((nb_fits, max_nb_clones+extra))
    aic = np.zeros((nb_fits, max_nb_clones+extra))
    aicc = np.zeros((nb_fits, max_nb_clones+extra))
    icl_q = np.zeros((nb_fits, max_nb_clones+extra))
    icl_qn = np.zeros((nb_fits, max_nb_clones+extra))
    bic_alt = np.zeros((nb_fits, max_nb_clones+extra))
    for i, s in enumerate(seeds):
        np.random.seed(s)
        for j, nb_clones in enumerate(range(1, max_nb_clones+1+extra)):
            print(j, i)
            if nb_clones >= 2:
                to_split = np.argmax(-(est.qun * np.log(est.qun)).sum(axis=0))
                mask = np.ones(nb_clones-1, dtype=bool)
                mask[to_split] = 0
                new_phi = np.zeros(nb_clones)
                new_phi[:nb_clones - 2] = est.phi[mask]
                new_phi[-2] = np.random.ranf() * 0.8 + 0.1
                new_phi[-1] = np.random.ranf() * 0.8 + 0.1
                new_xi = np.zeros(nb_clones)
                new_xi[:nb_clones - 2] = est.xi[mask]
                new_xi[-1], new_xi[-2] = [est.xi[to_split]] * 2
                new_pi = np.zeros((nb_clones, inputMU.shape[0]))
                new_pi[:nb_clones - 2, :] = est.pi[mask, :]
                new_pi[-1, :] = np.random.dirichlet(alpha=np.ones(inputMU.shape[0]))
                new_pi[-2, :] = np.random.dirichlet(alpha=np.ones(inputMU.shape[0]))
                est = Estimator(T, B, C_normal, C_tumor_tot,
                                C_tumor_minor, D, purity, nb_clones,
                                inputMU=inputMU, pi=new_pi, phi=new_phi, xi=new_xi)
            else:
                est = Estimator(T, B, C_normal, C_tumor_tot,
                                C_tumor_minor, D, purity, nb_clones,
                                inputMU=inputMU)
            est.fit()
            print(nb_clones, est.tau)
            Fres[i, j] = est.Fs[-1]
            bic[i, j] = est.get_bic()
            bic_alt[i, j] = np.nan
            loglikelihood[i, j] = est.get_loglikelihood
            aic[i, j] = est.get_aic()
            aicc[i, j] = est.get_aicc()
            icl_q[i, j] = est.get_icl()
            icl_qn[i, j] = est.get_icl(norm=True)
    dict_results = {'bic': np.argmax(bic.mean(axis=0)) + 1,
                    'aic': np.argmax(aic.mean(axis=0)) + 1,
                    'aicc': np.argmax(aicc.mean(axis=0)) + 1,
                    'icl_q': np.argmax(icl_q.mean(axis=0)) + 1,
                    'icl_qn': np.argmax(icl_qn.mean(axis=0)) + 1,
                    'bic_alt': np.argmax(bic_alt.mean(axis=0)) + 1}

    # compute SH estimate
    for mc in range(max_nb_clones-2, max_nb_clones + extra + 1):
        slopes = list()
        chpt = list()
        for end_p in range(0, mc-1):
            ransac = linear_model.LinearRegression()
            ransac.fit(((np.array(range(end_p+1, mc+1)))*(L+1)).reshape(-1, 1), loglikelihood.mean(axis=0)[end_p:mc])
            slopes.append(ransac.coef_)
            # print('pen', mc, end_p, loglikelihood[0][:mc] - np.arange(1, len(loglikelihood[0][:mc])+1) * (M+1) * 2 * ransac.estimator_.coef_)
            chpt.append(np.argmax(loglikelihood.mean(axis=0)[:mc] - np.arange(1, len(loglikelihood.mean(axis=0)[:mc])+1) * (L+1) * 2 * max(ransac.coef_, 0.0))+1)
        chpt = np.array(chpt)
        diff = chpt[1:] - chpt[0:-1]
        last_point = np.argmax(diff<0)
        if (last_point == 0) & (chpt[1] >= chpt[0]):
            last_point = mc
        counts = np.bincount(chpt[:last_point+1])
        # b = counts[::-1]
        # final_nb_clones = len(b) - np.argmax(b) - 1
        final_nb_clones = np.argmax(counts)
        dict_results['sh_{}'.format(mc)] = final_nb_clones

    ll = loglikelihood.mean(axis=0)
    dict_results['max_curvature'] = np.argmax(
        np.abs(ll[2:] + ll[0: -2] - 2 * ll[1: -1])) + 2
    return dict_results, ll


data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
with open('{}/purity.txt'.format(folder_path), 'r') as f:
    purity = float(f.read())

with open('{}/sim_data'.format(folder_path), 'rb') as sim_pickle_file:
    sim_pickle = pickle.Unpickler(sim_pickle_file)
    sim_data_obj = sim_pickle.load()
dist_matrix = sp.spatial.distance.squareform(sp.spatial.distance.pdist(sim_data_obj.true_pi.dot(sim_data_obj.MU), 'cosine'))
if len(dist_matrix[dist_matrix > 0]):
    min_dist = np.min(dist_matrix[dist_matrix > 0])
    max_dist = np.max(dist_matrix[dist_matrix > 0])
    avg_dist = np.mean(dist_matrix[dist_matrix > 0])
else:
    min_dist, max_dist, avg_dist = np.nan, np.nan, np.nan
avg_major_cn = np.mean(sim_data_obj.C_tumor_tot - sim_data_obj.C_tumor_minor)
avg_tot_cn = np.mean(sim_data_obj.C_tumor_tot)
actual_perc_diploid = sum((sim_data_obj.C_tumor_tot == 2) & (sim_data_obj.C_tumor_minor == 1))/sim_data_obj.N


# get clonesig results
start = time.time()
dict_results, ll = fit_model_special(data_df.trinucleotide.values, data_df.var_counts.values,
                                     data_df.normal_cn.values, data_df.minor_cn.values + data_df.major_cn.values,
                                     data_df.minor_cn.values,
                                     data_df.var_counts.values + data_df.ref_counts.values, purity,
                                     inputMU=MU)
sub_dict_results, sub_ll = fit_model_special(data_df.trinucleotide.values, data_df.var_counts.values,
                                             data_df.normal_cn.values, data_df.minor_cn.values + data_df.major_cn.values,
                                             data_df.minor_cn.values,
                                             data_df.var_counts.values + data_df.ref_counts.values, purity,
                                             inputMU=subMU)
total_time = time.time() - start

sample_id_cols = ['cancer_type', 'perc_diploid', 'nb_clones', 'nb_mut', 'nb_sig', 'min_dist', 'max_dist', 'avg_dist', 'avg_major_cn', 'actual_perc_diploid', 'avg_tot_cn']
id_df = pd.DataFrame([[cancer_type, perc_diploid, nb_clones, nb_mut,
                       MU.shape[0], min_dist, max_dist, avg_dist,
                       avg_major_cn, actual_perc_diploid, avg_tot_cn],
                      [cancer_type, perc_diploid, nb_clones, nb_mut,
                       subMU.shape[0], min_dist, max_dist, avg_dist,
                       avg_major_cn, actual_perc_diploid, avg_tot_cn]], columns=sample_id_cols)
ll_df = pd.DataFrame(np.vstack((ll, sub_ll)), columns=['ll_{}'.format(i+1) for i in range(len(ll))])
dict_df = pd.concat((pd.DataFrame.from_dict({k: list([v]) for k, v in dict_results.items()}), pd.DataFrame.from_dict({k:list([v]) for k, v in sub_dict_results.items()})))
dict_df = dict_df.reset_index()
res_df = pd.concat([id_df, dict_df, ll_df], axis=1)

res_df.to_csv('{}/clonesig_cn_cancer_type_features.tsv'.format(folder_path), sep='\t',
              index=False)
print('+', '{}/clonesig_cn_cancer_type_features.tsv'.format(folder_path), '+', sep='')
