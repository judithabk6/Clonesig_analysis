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
from clonesig.run_clonesig import get_MU



MU = get_MU()
L, K = MU.shape

# get pi values (haha) that cover the spectrum
np.random.seed(7)

dist_list = list()
for i in range(10000):
    nb_active_sig = np.random.poisson(7) + 1
    active_signatures = np.random.choice(L, nb_active_sig, replace=False)
    pi = np.zeros((2, L))
    pi[0, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    pi[1, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    dist_list.append(cosine(pi[0, :].dot(MU), pi[1, :].dot(MU)))


d = sm.nonparametric.KDEUnivariate(dist_list)
d.fit()

M_r = 20
valid_pi = list()
pi_dist_kept = list()
while len(valid_pi) < 50:
    u = np.random.random()
    nb_active_sig = np.random.poisson(7) + 1
    active_signatures = np.random.choice(L, nb_active_sig, replace=False)
    pi = np.zeros((2, L))
    pi[0, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    pi[1, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
    y = cosine(pi[0, :].dot(MU), pi[1, :].dot(MU))
    if u < 1.0/(M_r * d.evaluate(y)):
        valid_pi.append(pi)
        pi_dist_kept.append(y)

# sns.distplot(dist_list)
# sns.distplot(pi_dist_kept)
# plt.show()


# simulate data
expname = '20200305_simulations_eval_clonesig_power'
safe_mkdir(expname)
for nb_pi, pi in enumerate(valid_pi):
    for nb_phi, subclone_phi in enumerate((1-np.logspace(-0.5, 0.95, 10)/10)):
        for nb_mut in (30, 100, 300, 1000):
            for perc_dip in (0.1, 0.5, 0.9):
                for depth in (100, 500):
                    foldername = ('{}/pi{}-phi{}-depth{}-percdip{}-nb_mut{}'.
                                  format(expname, nb_pi, nb_phi, depth,
                                         perc_dip, nb_mut))
                    uu = SimLoader(nb_mut, 2, inputMU=MU, xi_param=[0.5, 0.5],
                                   pi_param=pi, phi_param=[1, subclone_phi],
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
                    uu.write_dpclust(foldername)
                    uu.write_phylogicndt(foldername)
                    uu.write_tracksigfreq(foldername)


