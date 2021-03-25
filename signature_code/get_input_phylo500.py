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

"""
samplename = "Sim_500_19"
"""

np.random.seed(int(samplename.split('_')[2]))
data_folder_path = 'data_simu_pcawg/training_PhylogicNDT500_simulations_12_16'
maf_filename = '{}/Mafs/{}.no_real_info.maf.tsv'.format(data_folder_path, samplename)
cnv_filename = '{}/BB/{}.BB_file.txt'.format(data_folder_path, samplename)
mut_assign_filename = '{}/Mut_Assign/{}.mutation_assignments.txt'.format(data_folder_path, samplename)
purity_filename = '{}/pp_table.txt'.format(data_folder_path)
cancer_type_filename = '{}/cancer_type.txt'.format(data_folder_path)

ctype = pd.read_csv(cancer_type_filename, sep='\t')
pre_loc = ctype[ctype['sample']==samplename].cancer_type.unique()[0]

loc_match = {'Bone-Cart': 'Bone-Other',
             'Bone-Epith': 'Bone-Other',
             'Bone-Leiomyo': 'Bone-Other',
             'Breast-AdenoCa': 'Breast',
             'Breast-LobularCa': 'Breast',
             'Cervix-AdenoCA': 'Cervix-SCC',
             'Myeloid-MDS': 'Myeloid-MDS/MPN',
             'Myeloid-MPN': 'Myeloid-MDS/MPN'}
if pre_loc in loc_match:
    cancer_loc = loc_match[pre_loc]
else:
    cancer_loc = pre_loc


maf_df = pd.read_csv(maf_filename, sep='\t')

maf_df = maf_df.assign(chromosome=maf_df['Chromosome'].astype(int))
maf_df = maf_df.assign(position=maf_df['Start_position'].astype(int))
maf_df = maf_df.assign(mutation_id=maf_df.chromosome.astype(str) + '_' +
                           maf_df.position.astype(str))
maf_df = maf_df.assign(ref_counts=maf_df.t_ref_count)
maf_df = maf_df.assign(var_counts=maf_df.t_alt_count)

# get cancer type MU matrix
signature_filename = 'data/sigProfiler_SBS_signatures_2018_03_28.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')

match_cancertype = pd.read_csv(
  'external_data/curated_match_signature_cancertype_pcawgs_literature.csv',
  index_col=0, sep='\t')
match_cancertype.index = match_cancertype['Unnamed: 0.1']


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

maf_df_cluster = pd.merge(maf_df, mut_assign, on='mutation_id')

# get signature and trinucleotide
N = maf_df_cluster.shape[0]
varS = np.zeros(N)
cstS = np.zeros(N)
J = maf_df_cluster.cluster.nunique()
U = maf_df_cluster.cluster.values
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

maf_df_cluster = maf_df_cluster.assign(varSig=varS,
                                       varTrinucleotide=pd.Categorical(varT, categories=list(range(K)), ordered=True),
                                       cstSig=cstS,
                                       cstTrinucleotide=pd.Categorical(cstT, categories=list(range(K)), ordered=True))

# CNV data
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

vcf_data_cnv = pd.merge(maf_df_cluster, cnv_final, left_on='chromosome',
                        right_on='chr', how='left')
vcf_data_cnv_filter = vcf_data_cnv[(vcf_data_cnv.position>=vcf_data_cnv.startpos)&(vcf_data_cnv.position<=vcf_data_cnv.endpos)]
vcf_data_cnv_filter = vcf_data_cnv_filter.assign(normal_cn=2)

purity_table = pd.read_csv(purity_filename, sep='\t')
purity = purity_table[purity_table['sample']==samplename].purity.mean()
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
    data.cn_profile = cnv_final.copy()
    data.purity = purity

    foldername = 'PhylogicNDT500/{}_{}'.format(samplename, suffixes[ii])
    pathlib.Path(foldername).mkdir(parents=True, exist_ok=True)
    data.write_clonesig(foldername)
    data.data = data.data.assign(total_cn=data.data.major_cn + data.data.minor_cn)
    data.write_pyclone_sciclone_ccube(foldername)
    data.write_deconstructsig(foldername)
    data.data = data.data.assign(mut_cn=1)
    data.write_tracksig(foldername)
    data.write_tracksigfreq(foldername)
    data.cn_profile = data.cn_profile.assign(chromosome=data.cn_profile.chr,
                                             start=data.cn_profile.startpos,
                                             end=data.cn_profile.endpos,
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
    cnv_final = cnv_final.assign(weight=cnv_final.endpos - cnv_final.startpos)
    perc_dip = cnv_final[(cnv_final.major_cn==1)&(cnv_final.minor_cn==1)].weight.sum() / cnv_final.weight.sum()
    info_df = pd.DataFrame([[len(df), perc_dip, np.median(df.var_counts + df.ref_counts), cancer_loc, suffixes[ii]]],
                            columns=['nb_mut', 'perc_dip', 'median_depth', 'cancer_loc', 'var_cst'])
    info_df.to_csv('{}/info_df.csv'.format(foldername), index=False, sep='\t')
    MU_df_norm.to_csv('{}/subMU.csv'.format(foldername), index=True, sep='\t')
