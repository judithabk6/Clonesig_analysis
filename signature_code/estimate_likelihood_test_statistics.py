#!/usr/bin/env python
# -*- coding:utf-8 -*-
import pandas as pd
import sys
from collections import Iterable
import numpy as np
from util_functions import safe_mkdir
from clonesig.estimator import Estimator
from clonesig.data_loader import SimLoader
from sklearn import linear_model

nb_clones = int(sys.argv[1])
nb_mut = int(sys.argv[2])
cancer_type = int(sys.argv[3])
perc_dip = float(sys.argv[4]) / 100
nb_seed = int(sys.argv[5])

output_dir = '20200216_estimate_likelihood_test_statistics'
safe_mkdir(output_dir)

'''
nb_mut = 100
nb_clones = 4
nb_seed = 5
cancer_type=2
perc_dip=0.2
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
MU = MU[cancer_type_sig.sum(axis=1).astype(bool), :]

# simulate pi
nb_active_sig = np.min((np.random.poisson(7) + 1, subMU.shape[0]))
active_signatures = np.random.choice(subMU.shape[0], nb_active_sig, replace=False)
pi = np.zeros((nb_clones, subMU.shape[0]))
pi[0, active_signatures] = np.random.dirichlet(alpha=np.ones(nb_active_sig))
for i in range(1, nb_clones):
    pi[i, :] = pi[0, :]

# simulate phi
steady_phi = np.zeros(nb_clones)
steady_phi[0] = 1.0
# get steady_phi in decreasing order
for i in range(1, nb_clones):
    steady_phi[i] = np.random.uniform(
        low=0.1 + 0.1*(nb_clones - i - 1),
        high=steady_phi[i-1] - 0.1)

# get xi
steady_xi = np.random.dirichlet(alpha=np.ones(nb_clones))
while min(steady_xi) < 0.1:
    steady_xi = np.random.dirichlet(alpha=np.ones(nb_clones))

np.random.seed(20190610 + nb_seed)
uu = SimLoader(nb_mut, nb_clones, inputMU=subMU,
               pi_param=pi, phi_param=steady_phi, xi_param=steady_xi,
               rho_param=100, cn=True, dip_prop=perc_dip)
uu._get_unobserved_nodes()
uu._get_observed_nodes()


sc_est_subMU = Estimator(uu.T, uu.B, uu.C_normal, uu.C_tumor_tot,
                         uu.C_tumor_minor, uu.D, uu.purity, 1,
                         inputMU=subMU)
sc_est_subMU.fit()
est_subMU = Estimator(uu.T, uu.B, uu.C_normal, uu.C_tumor_tot,
                      uu.C_tumor_minor, uu.D, uu.purity, nb_clones,
                      inputMU=subMU)
est_subMU.fit()

sc_est_MU = Estimator(uu.T, uu.B, uu.C_normal, uu.C_tumor_tot,
                      uu.C_tumor_minor, uu.D, uu.purity, 1,
                      inputMU=MU)
sc_est_MU.fit()
est_MU = Estimator(uu.T, uu.B, uu.C_normal, uu.C_tumor_tot,
                   uu.C_tumor_minor, uu.D, uu.purity, nb_clones,
                   inputMU=MU)
est_MU.fit()

loglikelihood_subMU = est_subMU.get_loglikelihood
alt_pi_subMU = np.repeat(sc_est_subMU.pi, nb_clones, axis=0)
est_subMU.pi = alt_pi_subMU
loglikelihood_nopi_subMU = est_subMU.get_loglikelihood

loglikelihood_MU = est_MU.get_loglikelihood
alt_pi_MU = np.repeat(sc_est_MU.pi, nb_clones, axis=0)
est_MU.pi = alt_pi_MU
loglikelihood_nopi_MU = est_MU.get_loglikelihood


cols = ['nb_mut', 'nb_clones', 'cancer_type', 'nb_sig_subMU', 'seed', 'perc_dip',
        'loglikelihood_MU', 'loglikelihood_nopi_MU', 'loglikelihood_subMU',
        'loglikelihood_nopi_subMU']

data = [nb_mut, nb_clones, cancer_type, subMU.shape[0], nb_seed, perc_dip,
        loglikelihood_MU, loglikelihood_nopi_MU, loglikelihood_subMU,
        loglikelihood_nopi_subMU]
res_df = pd.DataFrame([data], columns=cols)

res_df.to_csv('{}/cancertype{}_nbmut{}_nbclones{}_nbseed{}_percdip{}.csv'
              .format(output_dir, cancer_type, nb_mut, nb_clones, nb_seed, int(perc_dip*100)),
              sep='\t', index=False)

