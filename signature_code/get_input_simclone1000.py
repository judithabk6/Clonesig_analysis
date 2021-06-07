#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import pkg_resources
from clonesig import mixin_init_parameters
from clonesig.data_loader import DataWriter
from scipy.spatial.distance import cosine
import pathlib
import sys

THRESHOLD = 0.05

samplename = sys.argv[1]
samplenumber = int(sys.argv[2])

"""
samplename = "sim0h786c"
samplenumber = 8
sed 's/^\#\#/\&/g' data_simu_pcawg/testing_SimClone1000_truth/sim0h786c/simulated_0001/simulated_0001_0001/snv/icgc_consensus/simulated_0001_0001.vcf >  data_simu_pcawg/testing_SimClone1000_truth/sim0h786c/simulated_0001/simulated_0001_0001/snv/icgc_consensus/simulated_0001_0001_clean.vcf
"""

np.random.seed(samplenumber * 12)
data_folder_path = 'data_simu_pcawg/testing_SimClone1000_input'
truth_folder_path = 'data_simu_pcawg/testing_SimClone1000_truth'
vcf_filename = '{}/{}/simulated_0001/simulated_0001_0001/snv/icgc_consensus/simulated_0001_0001_clean.vcf'.format(truth_folder_path, samplename, samplename)
cnv_filename = '{}/{}/{}_segments.txt'.format(data_folder_path, samplename, samplename)
mut_assign_filename = '{}/{}/simulated_0001/truth_tree/simulated_0001_mutation_assignments.txt'.format(truth_folder_path, samplename)
purity_filename = '{}/{}/{}_purity_ploidy.txt'.format(truth_folder_path, samplename, samplename)


vcf_data = pd.read_csv(vcf_filename, sep='\t', comment='&')
vcf_data = vcf_data.assign(chromosome=vcf_data['#CHROM'])
vcf_data = vcf_data.assign(position=vcf_data['POS'].astype(int))
vcf_data = vcf_data.assign(mutation_id=vcf_data.chromosome.astype(str) + '_' +
                           vcf_data.position.astype(str))
vcf_data = vcf_data.assign(ref_counts=vcf_data.INFO.str.split('t_ref_count=').str[1].str.split(';').str[0].astype(int))
vcf_data = vcf_data.assign(var_counts=vcf_data.INFO.str.split('t_alt_count=').str[1].str.split(';').str[0].astype(int))

signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')

match_cancertype = pd.read_csv(
  'external_data/curated_match_signature_cancertype_pcawgs_literature.csv',
  index_col=0, sep='\t')
match_cancertype.index = match_cancertype['Unnamed: 0.1']

cancer_loc = np.random.choice(match_cancertype.columns.to_list()[1:])

selected_sigs = match_cancertype[[cancer_loc]][match_cancertype[[cancer_loc]] > 0].dropna().index.tolist()
sub_sig = sig[[c for c in sig.columns if c in selected_sigs]]

sig_matrix = sub_sig.values.astype(float).T
new_sig_matrix = sig_matrix + mixin_init_parameters.ZERO_PADDING \
    * (sig_matrix == 0)
MU = new_sig_matrix / new_sig_matrix.sum(axis=1)[:, np.newaxis]

final_cols_used = [c for c in sig.columns if c in selected_sigs]
MU_df_norm = pd.DataFrame(MU.T, index=sig.Type + '_' + sig.SubType,
                          columns=final_cols_used)

# get true PI matrix
mut_assign = pd.read_csv(mut_assign_filename, sep='\t')
mut_assign = mut_assign.assign(mutation_id=mut_assign.chr.astype(str) + '_' +
                           mut_assign.pos.astype(str))
nb_clones = mut_assign.cluster.nunique()

steady_pi = np.zeros((nb_clones, MU.shape[0]))
nb_active_sig = np.min((np.random.poisson(7) + 1, MU.shape[0]))
active_signatures = np.random.choice(MU.shape[0], nb_active_sig, replace=False)
steady_pi[0, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))
for i in range(1, nb_clones):
    steady_pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))
    while cosine(steady_pi[i, :].dot(MU), steady_pi[i-1, :].dot(MU)) < THRESHOLD:
        steady_pi[i, active_signatures] = np.random.dirichlet(alpha=np.ones(len(active_signatures)))

vcf_data_cluster = pd.merge(vcf_data, mut_assign, on='mutation_id')

# get signature and trinucleotide
N = vcf_data_cluster.shape[0]
varS = np.zeros(N)
cstS = np.zeros(N)
J = vcf_data_cluster.cluster.nunique()
U = vcf_data_cluster.cluster.values
L, K = MU.shape
for i in range(J):
    varS[U == i] = np.random.choice(L, sum(U == i),
                                    replace=True,
                                    p=steady_pi[i, :])
    cstS[U == i] = np.random.choice(L, sum(U == i),
                                    replace=True,
                                    p=steady_pi[0, :])

varT = np.zeros(N)
cstT = np.zeros(N)
for i in range(L):
    varT[varS == i] = np.random.choice(K, sum(varS == i),
                                       replace=True,
                                       p=MU[i, :])
    cstT[cstS == i] = np.random.choice(K, sum(cstS == i),
                                       replace=True,
                                       p=MU[i, :])

vcf_data_cluster = vcf_data_cluster.assign(varSig=varS,
                                           varTrinucleotide=pd.Categorical(varT, categories=list(range(K)), ordered=True),
                                           cstSig=cstS,
                                           cstTrinucleotide=pd.Categorical(cstT, categories=list(range(K)), ordered=True))

# CNV data
cnv_table = pd.read_csv(cnv_filename, sep='\t')
vcf_data_cnv = pd.merge(vcf_data_cluster, cnv_table, left_on='chromosome',
                        right_on='chromosome', how='left')
vcf_data_cnv_filter = vcf_data_cnv[(vcf_data_cnv.position>=vcf_data_cnv.start)&(vcf_data_cnv.position<=vcf_data_cnv.end)]
vcf_data_cnv_filter = vcf_data_cnv_filter.assign(normal_cn=2)

purity_table = pd.read_csv(purity_filename, sep='\t')
purity = purity_table['rho.purityOut'].mean()

vcf_data_cnv_filter = vcf_data_cnv_filter.assign(purity=purity)

var_df = vcf_data_cnv_filter[['chromosome', 'position', 'ref_counts',
                                'var_counts', 'cluster', 'varSig',
                                'varTrinucleotide', 'major_cn', 'minor_cn',
                                'normal_cn', 'purity', 'mutation_id']]
var_df.columns = ['chromosome', 'position', 'ref_counts', 'var_counts',
                    'clone', 'signature', 'trinucleotide', 'major_cn',
                    'minor_cn', 'normal_cn', 'purity', 'mutation_id']

cst_df = vcf_data_cnv_filter[['chromosome', 'position', 'ref_counts',
                                'var_counts', 'cluster', 'cstSig',
                                'cstTrinucleotide', 'major_cn', 'minor_cn',
                                'normal_cn', 'purity', 'mutation_id']]
cst_df.columns = ['chromosome', 'position', 'ref_counts', 'var_counts',
                    'clone', 'signature', 'trinucleotide', 'major_cn',
                    'minor_cn', 'normal_cn', 'purity', 'mutation_id']

suffixes = ['var', 'cst']
for ii, df in enumerate([var_df, cst_df]):
    data = DataWriter()
    data.data = df.copy()
    data.cn_profile = cnv_table.copy()
    data.purity = purity

    foldername = 'SimClone1000/{}_{}'.format(samplename, suffixes[ii])
    pathlib.Path(foldername).mkdir(parents=True, exist_ok=True)
    data.write_clonesig(foldername)
    data.data = data.data.assign(total_cn=data.data.major_cn + data.data.minor_cn)
    data.write_pyclone_sciclone_ccube(foldername)
    data.write_deconstructsig(foldername)
    data.data = data.data.assign(mut_cn=1)
    data.write_tracksig(foldername)
    data.write_tracksigfreq(foldername)
    data.cn_profile = data.cn_profile.assign(chromosome=data.cn_profile.chromosome,
                                             start=data.cn_profile.start,
                                             end=data.cn_profile.end,
                                             minor=data.cn_profile.minor_cn,
                                             major=data.cn_profile.major_cn)
    data.write_palimpsest(foldername) # need column chromosome in cn data
    data.write_dpclust(foldername)
    data.write_phylogicndt(foldername) # mut_cn
    # write pi matrix
    if ii==0:
        np.savetxt('{}/pi_matrix.csv'.format(foldername), steady_pi,
                   delimiter=',')
    else:
        np.savetxt('{}/pi_matrix.csv'.format(foldername),
                   np.repeat([steady_pi[0, :]], J, axis=0), delimiter=',')

    # write relevant info (nb mut, perc_diploid, median depth, tumor type, cst, var)
    cnv_table = cnv_table.assign(weight=cnv_table.end - cnv_table.start)
    perc_dip = cnv_table[(cnv_table.major_cn==1)&(cnv_table.minor_cn==1)].weight.sum() / cnv_table.weight.sum()
    info_df = pd.DataFrame([[len(df), perc_dip, np.median(df.var_counts + df.ref_counts), cancer_loc, suffixes[ii]]],
                            columns=['nb_mut', 'perc_dip', 'median_depth', 'cancer_loc', 'var_cst'])
    info_df.to_csv('{}/info_df.csv'.format(foldername), index=False, sep='\t')
    MU_df_norm.to_csv('{}/subMU.csv'.format(foldername), index=True, sep='\t')



