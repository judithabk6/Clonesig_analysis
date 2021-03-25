#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import sys
import numpy as np
import pickle
from clonesig.run_clonesig import run_clonesig
import bz2
import signal
import time

# add timeout.
class TimeoutException(Exception):
    pass

def signal_handler(signum, frame):
    raise TimeoutException("Timed out!")


folder_path = sys.argv[1]
# folder_path = 'PhylogicNDT500/Sim_500_19_var'


MU_df = pd.read_csv('{}/subMU.csv'.format(folder_path), sep='\t', index_col=0)
MU_matrix = MU_df.values

MU = MU_matrix.T
model_selection_kws = {'factor':  0.093}

data_df = pd.read_csv('{}/input_t.tsv'.format(folder_path), sep='\t')
with open('{}/purity.txt'.format(folder_path), 'r') as f:
    purity = float(f.read())
    
start = time.time()
signal.signal(signal.SIGALRM, signal_handler)
signal.alarm(48*3600)   # 48 hours
try:
    new_est, lr, pval, new_inputMU, cst_est, fitted_sigs = run_clonesig(
        data_df.trinucleotide.values, data_df.var_counts.values,
        data_df.var_counts.values + data_df.ref_counts.values,
        data_df.normal_cn.values,
        data_df.minor_cn.values + data_df.major_cn.values,
        data_df.minor_cn.values, purity, MU, inputNu=None, nu_heuristics=None,
        return_sig_change_test=True, min_mut_clone=0, min_prop_sig=0.0,
        prefit_signatures=False, prefit_thresh=0.01, model_selection_function=None,
        model_selection_kws=model_selection_kws, max_nb_clones=10)
except TimeoutException:
    print('too long')
end = time.time()

raw_res = bz2.BZ2File('{}/cancertype_clonesig_raw_results.bz2'.format(folder_path), 'wb')
my_pickler = pickle.Pickler(raw_res)
my_pickler.dump([new_est, lr, pval, cst_est, fitted_sigs, end-start])
raw_res.close()