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
import bz2
import signal

folder_path = sys.argv[1]
# folder_path="salcedo_dream_challenge/T2_16X"

tumor, depth = folder_path.split('/')[-1].split('_')

# add timeout.
class TimeoutException(Exception):
    pass

def signal_handler(signum, frame):
    raise TimeoutException("Timed out!")



signature_filename = 'data/salcedo_dream_challenge/signatures.txt'

signature_match = {'T2': ['Signature 1A', 'Signature 1B', 'Signature 4', # liver-HCC
                          'Signature 5', 'Signature 6', 'Signature 9',
                          'Signature 12', 'Signature 14', 'Signature 16',
                          'Signature 17', 'Signature 18', 'Signature 19'],
                   'T3': ['Signature 1A', 'Signature 1B', 'Signature 2', # lung-AdenoCA
                          'Signature 3', 'Signature 4', 'Signature 5',
                          'Signature 8', 'Signature 13', 'Signature 17',
                          'Signature 18'],
                   'T4': ['Signature 1A', 'Signature 1B', 'Signature 2', # Breast
                          'Signature 3', 'Signature 5', 'Signature 8',
                          'Signature 9', 'Signature 13', 'Signature 17',
                          'Signature 18'],
                   'T5': ['Signature 1A', 'Signature 1B', 'Signature 2', # Breast
                          'Signature 3', 'Signature 5', 'Signature 8',
                          'Signature 9', 'Signature 13', 'Signature 17',
                          'Signature 18'],
                   'T6': ['Signature PCAWG-1', 'Signature PCAWG-5A', # Colorect-AdenoCA
                          'Signature PCAWG-5B', 'Signature 10A',
                          'Signature 10B', 'Signature PCAWG-15',
                          'Signature PCAWG-17A', 'Signature PCAWG-17B',
                          'Signature PCAWG-18', 'Signature PCAWG-28',
                          'Signature PCAWG-N3']}

signature_df = pd.read_csv(signature_filename, sep='\t')
MU_cols = ['Signature 1A', 'Signature 1B', 'Signature 2', 'Signature 3',
           'Signature 4', 'Signature 5', 'Signature 6', 'Signature 7',
           'Signature 8', 'Signature 9', 'Signature 10', 'Signature 11',
           'Signature 12', 'Signature 13', 'Signature 14', 'Signature 15',
           'Signature 16', 'Signature 17', 'Signature 18', 'Signature 19',
           'Signature 20', 'Signature 21', 'Signature R1', 'Signature R2',
           'Signature R3', 'Signature U1', 'Signature U2']
if tumor == 'T6':
    MU_cols = ['Signature PCAWG-1', 'Signature PCAWG-2', 'Signature PCAWG-3',
               'Signature PCAWG-4', 'Signature PCAWG-5A', 'Signature PCAWG-5B',
               'Signature PCAWG-6', 'Signature PCAWG-7A', 'Signature PCAWG-7B',
               'Signature PCAWG-7C', 'Signature PCAWG-7D', 'Signature PCAWG-8',
               'Signature PCAWG-9', 'Signature 10A', 'Signature 10B',
               'Signature PCAWG-11', 'Signature PCAWG-12', 'Signature PCAWG-13',
               'Signature PCAWG-14', 'Signature PCAWG-15', 'Signature PCAWG-16',
               'Signature PCAWG-17A', 'Signature PCAWG-17B',
               'Signature PCAWG-18', 'Signature PCAWG-19', 'Signature PCAWG-20',
               'Signature PCAWG-21', 'Signature PCAWG-22', 'Signature PCAWG-24',
               'Signature PCAWG-26', 'Signature PCAWG-28', 'Signature PCAWG-30',
               'Signature PCAWG-R1', 'Signature PCAWG-R2', 'Signature SNPs',
               'Signature PCAWG-N1', 'Signature PCAWG-N2', 'Signature PCAWG-N3',
               'Signature PCAWG-N4', 'Signature PCAWG-N5', 'Signature PCAWG-N6',
               'Signature PCAWG-N7', 'Signature PCAWG-N8']

MU_df = signature_df[MU_cols]
MU_df.index = signature_df['Somatic Mutation Type']
SIG_MATRIX = MU_df.values
NEW_SIG_MATRIX = SIG_MATRIX + 10**-20 * (SIG_MATRIX == 0)
MU = NEW_SIG_MATRIX / NEW_SIG_MATRIX.sum(axis=0)[np.newaxis, :]
MU_df_norm = pd.DataFrame(MU, index=MU_df.index, columns=MU_df.columns)
sub_MU_df = MU_df_norm[signature_match[tumor]]
MU_df_norm.to_csv('{}/MU_matrix.csv'.format(folder_path), sep='\t', index=True)
sub_MU_df.to_csv('{}/sub_MU_matrix.csv'.format(folder_path), sep='\t', index=True)

MU_matrix = MU_df_norm.values
sub_MU_matrix = sub_MU_df.values

data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
with open('{}/purity.txt'.format(folder_path), 'r') as f:
    purity = float(f.read())

method = 'clonesig'
id_list = list()
metrics_list = list()
for setting in ('cancertype', 'prefit', 'all', 'all_nuclonal'):
    if setting == 'cancertype':
        MU = sub_MU_matrix.T
        model_selection_kws = {'factor':  0.093}
    else:
        MU = MU_matrix.T
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
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(48*3600)   # 48 hours
    try:
        new_est, lr, pval, new_inputMU, cst_est, fitted_sigs = run_clonesig(
            data_df.trinucleotide.values, data_df.var_counts.values,
            data_df.var_counts.values + data_df.ref_counts.values,
            data_df.normal_cn.values,
            data_df.minor_cn.values + data_df.major_cn.values,
            data_df.minor_cn.values, purity, MU, inputNu=None, nu_heuristics=nuh,
            return_sig_change_test=True, min_mut_clone=0, min_prop_sig=0.0,
            prefit_signatures=pf, prefit_thresh=0.01, model_selection_function=None,
            model_selection_kws=model_selection_kws, max_nb_clones=10)
    except TimeoutException:
        print('too long')
        continue
    end = time.time()
    raw_res = bz2.BZ2File('{}/{}_clonesig_raw_results.bz2'.format(folder_path, setting), 'wb')
    my_pickler = pickle.Pickler(raw_res)
    my_pickler.dump([new_est, lr, pval, cst_est, fitted_sigs, end-start])


