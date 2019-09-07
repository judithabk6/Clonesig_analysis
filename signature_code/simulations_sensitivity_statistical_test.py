#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import uniform
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from clonesig.data_loader import SimLoader
from util_functions import safe_mkdir
from scipy.stats import beta
import sys
import scipy as sp


sig_file_path = 'external_data/sigProfiler_SBS_signatures_2018_03_28.csv'
cancer_type_sig_filename = 'external_data/match_cancer_type_sig_v3.csv'

# open the matrix describing the signatures
SIG = pd.read_csv(sig_file_path)
SIG_MATRIX = SIG.values[:, 2:].astype(float).T
L, K = SIG_MATRIX.shape
NEW_SIG_MATRIX = SIG_MATRIX + 10**-20 * (SIG_MATRIX == 0)
MU = NEW_SIG_MATRIX / NEW_SIG_MATRIX.sum(axis=1)[:, np.newaxis]

nb_clones = int(sys.argv[1])
# get pi values (haha) that cover the spectrum
np.random.seed(7)

min_dist_list = list()
max_dist_list = list()
avg_dist_list = list()
for i in range(10000):
    nb_active_sig = np.random.poisson(5) + 2
    active_signatures = np.random.choice(L, nb_active_sig, replace=False)
    pi = np.zeros((nb_clones, L))
    for i in range(nb_clones):
        pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    dist_matrix = sp.spatial.distance.squareform(sp.spatial.distance.pdist(pi.dot(MU), 'cosine'))
    min_dist = np.min(dist_matrix[dist_matrix > 0])
    max_dist = np.max(dist_matrix[dist_matrix > 0])
    avg_dist = np.mean(dist_matrix[dist_matrix > 0])
    min_dist_list.append(min_dist)
    max_dist_list.append(max_dist)
    avg_dist_list.append(avg_dist)
d = sm.nonparametric.KDEUnivariate(max_dist_list)
d.fit()

M_r = 100
valid_pi = list()
pi_dist_kept = list()
while len(valid_pi) < 30:
    u = np.random.random()
    nb_active_sig = np.random.poisson(5) + 2
    active_signatures = np.random.choice(L, nb_active_sig, replace=False)
    pi = np.zeros((nb_clones, L))
    for i in range(nb_clones):
        pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    dist_matrix = sp.spatial.distance.squareform(sp.spatial.distance.pdist(pi.dot(MU), 'cosine'))

    max_dist = np.max(dist_matrix[dist_matrix > 0])

    if u < beta.pdf(max_dist, 1.5, 8)/(M_r * d.evaluate(max_dist)):
        valid_pi.append(pi)
        pi_dist_kept.append(max_dist)


nb_phi = nb_clones
# simulate data
expname = '20190804_simulations_eval_clonesig_power'
safe_mkdir(expname)
for nb_pi, pi in enumerate(valid_pi):
    for nb_mut in (100, 300, 1000):
        for perc_dip in (0.1, 0.5, 0.9):
            for depth in (100, 500):
                foldername = ('{}/pi{}-phi{}-depth{}-percdip{}-nb_mut{}'.
                              format(expname, nb_pi, nb_phi, depth,
                                     perc_dip, nb_mut))
                uu = SimLoader(nb_mut, nb_clones, inputMU=MU, xi_param=np.array([1.0/nb_clones] * nb_clones),
                               pi_param=pi, phi_param=np.linspace(0.1, 1, nb_clones),
                               rho_param=100, cn=True, D_param=depth,
                               purity_param=0.8, dip_prop=perc_dip)
                uu._get_unobserved_nodes()
                uu._get_observed_nodes()
                uu.write_object(foldername)
                uu.write_clonesig(foldername)
                uu.write_pyclone_sciclone_ccube(foldername)
                uu.write_deconstructsig(foldername)
                uu.write_tracksig(foldername)
                uu.write_palimpsest(foldername)



