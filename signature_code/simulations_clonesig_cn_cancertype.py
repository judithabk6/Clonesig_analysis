#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import sys
import itertools
import pandas as pd
from scipy.spatial.distance import cosine


from util_functions import safe_mkdir
from clonesig.data_loader import SimLoader

nb_clones = int(sys.argv[1])
nb_mut = int(sys.argv[2])
cancer_type = int(sys.argv[3])
perc_dip = float(sys.argv[4]) / 100.0
seed_job = int(sys.argv[5])

'''
nb_clones = 2
nb_mut = 1000
cancer_type = 5
perc_dip = 20 / 100.0
seed_job = 4
'''

sig_file_path = 'external_data/sigProfiler_SBS_signatures_2018_03_28.csv'
cancer_type_sig_filename = 'external_data/curated_match_signature_cancertype_tcgawes_literature.csv'

# open the matrix describing the signatures
SIG = pd.read_csv(sig_file_path)
SIG_MATRIX = SIG.values[:, 2:].astype(float).T
L, K = SIG_MATRIX.shape
NEW_SIG_MATRIX = SIG_MATRIX + 10**-20 * (SIG_MATRIX == 0)
MU = NEW_SIG_MATRIX / NEW_SIG_MATRIX.sum(axis=1)[:, np.newaxis]

# get the signatures specific of the cancer type
cancer_type_sig = pd.read_csv(cancer_type_sig_filename, sep='\t', index_col=0).values
select = cancer_type_sig[:, cancer_type]
subMU = MU[select.astype(bool), :]

THRESHOLD = 0.05


expname = '20200210_simulations_clonesig_cn_cancer_type'
safe_mkdir(expname)
np.random.seed(seed_job)

# get pi !! careful not the tracksig setting anymore
steady_pi = np.zeros((nb_clones, subMU.shape[0]))
nb_active_sig = np.min((np.random.poisson(7) + 1, subMU.shape[0]))
active_signatures = np.random.choice(subMU.shape[0], nb_active_sig, replace=False)
steady_pi[0, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))
for i in range(1, nb_clones):
    steady_pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))
    while cosine(steady_pi[i, :].dot(subMU), steady_pi[i-1, :].dot(subMU)) < THRESHOLD:
        steady_pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))
# get phi
steady_phi = np.zeros(nb_clones)
steady_phi[0] = 1.0
# get steady_phi in decreasing order
for i in range(1, nb_clones):
    steady_phi[i] = np.random.uniform(
        low=0.1 + 0.1 * (nb_clones - i - 1),
        high=steady_phi[i-1] - 0.1)

# get xi
steady_xi = np.random.dirichlet(alpha=np.ones(nb_clones))
while min(steady_xi) < 0.1:
    steady_xi = np.random.dirichlet(alpha=np.ones(nb_clones))


foldername = ('{}/type{}-perc_diploid{}-nb_clones{}-nb_mut{}'.
              format(expname, cancer_type, int(perc_dip * 100), nb_clones, nb_mut))
uu = SimLoader(nb_mut, nb_clones, inputMU=subMU, xi_param=steady_xi,
               pi_param=steady_pi, phi_param=steady_phi,
               rho_param=100, cn=True, dip_prop=perc_dip)
uu._get_unobserved_nodes()
uu._get_observed_nodes()
uu.write_clonesig(foldername)
uu.write_pyclone_sciclone_ccube(foldername)
uu.write_deconstructsig(foldername)
uu.write_tracksig(foldername)
uu.write_tracksigfreq(foldername)
uu.write_palimpsest(foldername)
uu.write_dpclust(foldername)
uu.write_phylogicndt(foldername)
uu.write_object(foldername)






