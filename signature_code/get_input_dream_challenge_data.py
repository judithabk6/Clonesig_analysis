#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import sys
from pyfaidx import Faidx
from clonesig.data_loader import PAT_LIST, DataWriter
import pathlib
from scipy.stats import beta
import numpy as np


tumor = sys.argv[1]
seq_depth = sys.argv[2]

"""
tumor = "T2"
seq_depth = "8X"
#sed 's/^\#\#/\&/g' data/salcedo_dream_challenge/MuTect_inputs/mutect_filtered_T2.T.16XnoXY.vcf >  data/salcedo_dream_challenge/MuTect_inputs/mutect_filtered_T2.T.16XnoXY_clean.vcf
"""

np.random.seed(int(seq_depth[:-1])*int(tumor[1]))
def get_context(x):
    match_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if x['REF'] in ('C', 'T'):
        context = x['triplet'][0] + '[' + x['REF'] + '>' + x['ALT'] + ']' + \
            x['triplet'][2]
    else:
        context = match_dict[x['triplet'][2]] + '[' + match_dict[x['REF']] +\
            '>' + match_dict[x['ALT']] + ']' + match_dict[x['triplet'][0]]
    return context

vcf_filename = 'data/salcedo_dream_challenge/MuTect_inputs/mutect_filtered_{}.T.{}noXY_clean.vcf'.format(tumor, seq_depth)
cnv_filename = 'data/salcedo_dream_challenge/MuTect_inputs/{}-{}_refit_subclones_noXY.txt'.format(tumor, seq_depth)
purity_filename = 'data/salcedo_dream_challenge/MuTect_inputs/{}-{}_refit_cellularity_ploidy.txt'.format(tumor, seq_depth)
ref_fasta = 'external_data/hs37d5.fa'
fa = Faidx(ref_fasta)
# downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence

vcf_data = pd.read_csv(vcf_filename, sep='\t', comment='&')

vcf_data = vcf_data.assign(chromosome=vcf_data['#CHROM'].astype(int))
vcf_data = vcf_data.assign(position=vcf_data['POS'].astype(int))
vcf_data = vcf_data.assign(mutation_id=vcf_data.chromosome.astype(str) + '_' +
                           vcf_data.position.astype(str) + '_' +
                           vcf_data.REF + vcf_data.ALT)
vcf_data = vcf_data.assign(ref_counts=vcf_data.tumor.str.split(':').str[1].
                           str.split(',').str[0].astype(int))
vcf_data = vcf_data.assign(var_counts=vcf_data.tumor.str.split(':').str[1].
                           str.split(',').str[1].astype(int))
vcf_data = vcf_data.assign(triplet=vcf_data.apply(
    lambda x: fa.fetch(str(x.chromosome), x.position-1, x.position+1).seq,
    axis=1))
vcf_data = vcf_data.assign(pattern=pd.Categorical(
            vcf_data.apply(get_context, axis=1),
            categories=PAT_LIST, ordered=True))
vcf_data = vcf_data.assign(
        trinucleotide=pd.Categorical(
            vcf_data.apply(lambda x: PAT_LIST.index(x['pattern']),
                                     axis=1),
            categories=list(range(len(PAT_LIST))), ordered=True))
foldername = 'salcedo_dream_challenge/{}_{}'.format(tumor, seq_depth)
pathlib.Path(foldername).mkdir(parents=True, exist_ok=True)
vcf_data.to_csv(
    '{}/unrestricted_input_mut.csv'.format(foldername), sep='\t')
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

vcf_data_cnv = pd.merge(vcf_data, cnv_final, left_on='chromosome',
                        right_on='chr', how='left')
vcf_data_cnv_filter = vcf_data_cnv[(vcf_data_cnv.position>=vcf_data_cnv.startpos)&(vcf_data_cnv.position<=vcf_data_cnv.endpos)]
vcf_data_cnv_filter = vcf_data_cnv_filter.assign(normal_cn=2)



vcf_data_cnv_filter = vcf_data_cnv_filter[vcf_data_cnv_filter.major_cn>0]
purity_data = pd.read_csv(purity_filename, sep='\t')
purity = purity_data.cellularity.values[0]

data = DataWriter()
data.data = vcf_data_cnv_filter.copy()
data.cn_profile = cnv_final.copy()
data.purity = purity

foldername = 'salcedo_dream_challenge/{}_{}'.format(tumor, seq_depth)
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




