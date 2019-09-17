#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import os
import sys
from clonesig.run_clonesig import run_clonesig
import pkg_resources
from clonesig import mixin_init_parameters
from clonesig.estimator import EV_DOF_THRESHOLD
import time
from clonesig.data_loader import LE, get_context, PAT_LIST
import scipy as sp
import pickle
from clonesig.evaluate import score_sig_1A_base
import seaborn as sns
import matplotlib.pyplot as plt
from seaborn.distributions import _freedman_diaconis_bins
import matplotlib
import bz2
import collections

try:
    rows, columns = os.popen('stty size', 'r').read().split()
    pd.set_option('display.width', int(columns))
    pd.set_option('display.max_columns', 600)
    pd.set_option('display.max_rows', 1000)
except:
    print("running on server, otherwise please investigate")

patient_id = sys.argv[1]
cancer_loc = sys.argv[2]
try:
    cosmic_type = sys.argv[3]
except IndexError:
    cosmic_type = ''

"""
patient_id = 'TCGA-OR-A5JO'
cancer_loc = 'ACC'
cosmic_type = ''

patient_id = 'TCGA-A8-A06X'
cancer_loc = 'BRCA'
cosmic_type = 'Breast'


patient_id = 'TCGA-4H-AAAK'
cancer_loc = 'BRCA'
cosmic_type = 'Breast'

"""

color_dict = collections.OrderedDict()
color_dict['SBS1']= '#F2786D'
color_dict['SBS2']= '#D86009'
color_dict['SBS13']= '#D86009'
color_dict['SBS3']= '#DA8C05'
color_dict['SBS4']= '#A9A307'
color_dict['SBS5']= '#36B608'
color_dict['SBS6']= '#07BC7C'
color_dict['SBS7a']= '#f2b67c'
color_dict['SBS7b']= '#f2b67c'
color_dict['SBS7c']= '#f2b67c'
color_dict['SBS7d']= '#f2b67c'
color_dict['SBS8']= '#02BFC0'
color_dict['SBS9']= '#00AEF7'
color_dict['SBS10a']= '#66A428'
color_dict['SBS10b']= '#66A428'
color_dict['SBS11']= '#938CFF'
color_dict['SBS12']= '#E071EC'
color_dict['SBS14']= '#F566BE'
color_dict['SBS15']= '#7CCA7C'
color_dict['SBS16']= '#356AB2'
color_dict['SBS17a']= '#C15B00'
color_dict['SBS17b']= '#C15B00'
color_dict['SBS18']= '#666666'
color_dict['SBS19']= '#756EB5'
color_dict['SBS20']= '#E7AC00'
color_dict['SBS21']= '#A5CEE4'
color_dict['SBS22']= '#2CA121'
color_dict['SBS23']= '#335a9c'
color_dict['SBS24']= '#FEC068'
color_dict['SBS25']= '#a5f0b1'
color_dict['SBS26']= '#6A399C'
color_dict['SBS28']= '#E71408'
color_dict['SBS29']= '#4e88a6'
color_dict['SBS30']= '#994BA5'
color_dict['SBS31']= '#0cf531'
color_dict['SBS33']= '#A75620'
color_dict['SBS34']= '#62C3A4'
color_dict['SBS35']= '#E988C4'
color_dict['SBS36']= '#E6C591'
color_dict['SBS37']= '#BEADD5'
color_dict['SBS38']= '#F3007E'
color_dict['SBS39']= '#089F76'
color_dict['SBS40']= '#A77709'
color_dict['SBS41']= '#a80526'
color_dict['SBS42']= 'black'
color_dict['SBS43']= '#6b088a'
color_dict['SBS44']= '#fce428'
color_dict['SBS45']= '#f0d0f2'
color_dict['SBS46']= '#78aaff'
color_dict['SBS47']= '#ffb108'
color_dict['SBS49']= '#c108ff'
color_dict['SBS52']= '#c27289'
color_dict['SBS54']= '#1320d1'
color_dict['SBS56']= '#e8dba2'
color_dict['SBS58']= '#b0eb31'
color_dict['SBS59']= '#b8b8b8'

tcga_folder = '20190704_TCGA'


signature_filename = 'data/sigProfiler_exome_SBS_signatures.csv'
sig = pd.read_csv(
    pkg_resources.resource_stream(
        'clonesig', signature_filename),
    sep=',')
match_cancertype = pd.read_csv(
  'external_data/curated_match_signature_cancertype_tcgawes_literature.csv',
  index_col=0, sep='\t')

selected_sigs = match_cancertype[[cancer_loc]][match_cancertype[[cancer_loc]] > 0].dropna().index.tolist()
sub_sig = sig[[c for c in sig.columns if c in selected_sigs]]
prefit = False
sigprofiler = False

sig_matrix = sub_sig.values.astype(float).T
new_sig_matrix = sig_matrix + mixin_init_parameters.ZERO_PADDING \
    * (sig_matrix == 0)
MU = new_sig_matrix / new_sig_matrix.sum(axis=1)[:, np.newaxis]

final_cols_used = [c for c in sig.columns if c in selected_sigs]


clonal_sigs = pd.DataFrame(index=[0, 1], columns=sig.columns[2:])
subclonal_sigs = pd.DataFrame(index=[0, 1], columns=sig.columns[2:])
for i, folder_type in enumerate(['public', 'protected']):
    with bz2.BZ2File('{}/{}/{}_clonesig_raw_results_restr_curated.bz2'.format(tcga_folder, patient_id, folder_type), 'rb') as est_pickle_file:
        est_pickle = pickle.Unpickler(est_pickle_file)
        new_est, lr, pval, cst_est, fitted_sigs = est_pickle.load()

    clonal_idx = np.argmax(new_est.phi)
    clonal_phi = new_est.phi[clonal_idx]
    clonal_xi = new_est.xi[clonal_idx]
    if fitted_sigs is not None:
        actual_fitted = np.array(final_cols_used)[fitted_sigs]
        clonal_sigs.iloc[i, :][actual_fitted] = new_est.pi[clonal_idx, :]
        sq_dist = sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(new_est.pi.dot(MU[fitted_sigs, :]),
                                      'cosine'))
    else:
        clonal_sigs.iloc[i, :][final_cols_used] = new_est.pi[clonal_idx, :]
        sq_dist = sp.spatial.distance.squareform(
            sp.spatial.distance.pdist(new_est.pi.dot(MU), 'cosine'))
    nb_mut_clonal = sum(new_est.qun.argmax(axis=1) == clonal_idx)
    if new_est.J >= 2:
        sub_phi = new_est.phi.copy()
        sub_xi = new_est.xi.copy()
        sub_phi[clonal_idx] = sub_xi[clonal_idx] = 0
        largest_subclonal_idx = np.argmax(sub_xi)
        largest_subclonal_phi = new_est.phi[largest_subclonal_idx]
        largest_subclonal_xi = new_est.xi[largest_subclonal_idx]
        if fitted_sigs is not None:
            actual_fitted = np.array(final_cols_used)[fitted_sigs]
            subclonal_sigs.iloc[i, :][actual_fitted] = \
                new_est.pi[largest_subclonal_idx, :]
        else:
            subclonal_sigs.iloc[i, :][final_cols_used] = \
                new_est.pi[largest_subclonal_idx, :]
        nb_mut_largest_subclonal = \
            sum(new_est.qun.argmax(axis=1) == largest_subclonal_idx)
        clonal_largest_sub_pidist = sq_dist[clonal_idx, largest_subclonal_idx]
    else:
        (largest_subclonal_phi, largest_subclonal_xi, nb_mut_largest_subclonal,
         clonal_largest_sub_pidist) = [np.nan] * 4
    largest_pi_dist = np.max(sq_dist)

    if pval < 0.05:
        with bz2.BZ2File('{}/{}/{}_clonesig_raw_results_restr_curated.bz2'.format(tcga_folder, patient_id, folder_type), 'rb') as est_pickle_file:
            est_pickle = pickle.Unpickler(est_pickle_file)
            est_list = est_pickle.load()
            est = est_list[0]
        clonesig_res_patient = pd.read_csv('{}/{}/clonesig_results_restr_curated.csv'.format(tcga_folder, patient_id), sep='\t')

        cancer_loc = clonesig_res_patient.cancer_loc.unique()[0]
        cosmic_type = clonesig_res_patient.cosmic_type.unique()[0]

        mut_table = pd.DataFrame({'trinucleotide': est.T, 'var_counts': est.B,
                                  'minor_cn': est.C_tumor_minor,
                                  'major_cn': est.C_tumor_major,
                                  'total_cn': est.C_tumor_tot, 'depth': est.D,
                                  'clone': est.qun.argmax(axis=1),
                                  'signature': np.array(selected_sigs)[est.rnus[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1)],
                                  'mult': est.vmnu[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1) +1})
        mut_table = mut_table.assign(vaf=mut_table.var_counts / mut_table.depth)
        mut_table = mut_table.assign(vaf_cn=mut_table.vaf * mut_table['total_cn'] / mut_table['mult'])
        mut_table = mut_table.assign(vaf_purity=mut_table.apply(lambda x: x['vaf']/est.p * ((1 - est.p) * 2 + est.p * x['total_cn']) / x['mult'], axis=1))
        mut_table = mut_table.assign(trinucleotide=pd.Categorical(mut_table.trinucleotide, ordered=True, categories=range(96)))

        nb_bins = min(_freedman_diaconis_bins(mut_table.vaf_purity) * 2, 50)
        final_bins = np.linspace(min(mut_table.vaf_purity), max(mut_table.vaf_purity), nb_bins)
        # fig, ax = plt.subplots(1, figsize=(8, 28), sharex=False, gridspec_kw={'hspace': 0.08, 'wspace': 0, 'height_ratios': [1, 6, 1]})
        fig, ax = plt.subplots(1, figsize=(8, 5))
        clone_cols = sns.husl_palette(mut_table.clone.nunique(), l=0.8, s=.7)
        est_sigs = [s for s in selected_sigs if s in mut_table.signature.unique()]
        mylist = [color_dict[s] for s in est_sigs]
        my_palette = sns.color_palette(mylist)
        #cols = sns.color_palette("Set2", len(est_sigs))
        cols = sns.color_palette(my_palette, len(est_sigs))
        for i in range(est.J):
            sns.distplot(mut_table[mut_table.clone==i].vaf_purity, kde=False, label='clone {}'.format(i), rug=True, ax=ax, color=clone_cols[i], bins=final_bins)
        ax.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0., fontsize=14)
        plt.ylabel('count')
        ax.invert_xaxis()
        ax.set_xlim([1.5, 0])
        ax.set_xlabel('Cancer Cell Fraction', fontsize=17)
        ax.set_ylabel('Number of Mutations', fontsize=17)
        ax.set_title(patient_id, fontsize=20)
        plt.savefig('20190828_tcga_plots/{}_{}_clonesig_panel1.pdf'.format(patient_id, folder_type), bbox_inches='tight')
        fig, ax = plt.subplots(1, figsize=(8, 18))
        sns.scatterplot(x='vaf_purity', y='trinucleotide', data=mut_table,
                             hue='signature',
                             hue_order=[s for s in selected_sigs if s in mut_table.signature.unique()],
                             palette=cols, s=110, marker='d', ax=ax)
                             #hue='clone', palette=clone_cols, s=70, marker='d', ax=ax)
        ax.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0., fontsize=14)
        for i in range(96):
            ax.axhline(i, ls=':', lw=0.5, color='gray')
        for i in np.arange(-0.5, 96.5, 16):
            ax.axhline(i, color='gray', lw=1.2)
        ax.invert_xaxis()
        ax.set_xlim([1.5, 0])
        ax.set_ylim([-1, 96])

        ax.set_yticks(range(96))
        ax.get_yaxis().set_tick_params(direction='out', pad=18)
        ax.set_yticklabels(['{}{}{}'.format(c[0], c[2], c[-1]) for c in PAT_LIST], fontsize=8, color='gray', ha='left')

        axes2 = ax.twinx()   # mirror them
        axes2.set_ylim([-1, 96])
        axes2.set_yticks(np.arange(7.5, 96, 16))
        axes2.get_yaxis().set_tick_params(length=0, pad=6)
        _ = axes2.set_yticklabels(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'], fontsize=15, color='gray', ha='left', rotation=90, va='center')
        ax.set_xlabel('Cancer Cell Fraction', fontsize=17)
        ax.set_ylabel('Trinucleotide context', fontsize=17)
        ax.set_xlim([1.5, 0])
        ax.set_title(patient_id, fontsize=20)
        plt.savefig('20190828_tcga_plots/{}_{}_clonesig_panel2.pdf'.format(patient_id, folder_type), bbox_inches='tight')
        
        fig, ax = plt.subplots(1, figsize=(8, 5))
        for i, present_sig in enumerate(est_sigs):
            g=sns.kdeplot(mut_table[mut_table.signature==present_sig].vaf_purity, label=present_sig, ax=ax, color=cols[i])
        for i, present_sig in enumerate(est_sigs):
            dat = g.lines[i].get_xydata()
            dat[:,1] = dat[:,1]*len(mut_table[mut_table.signature==present_sig])/len(mut_table)
            g.lines[i].set_data(dat.T)
        ax.set_ylim([0, 2])

        for i in range(est.J):
            g=sns.distplot(mut_table[mut_table.clone==i].vaf_purity, kde=False, label='clone {}'.format(i), ax=ax, color=clone_cols[i], bins=final_bins)
        for jj in range(len(g.get_children())):
            if isinstance(g.get_children()[jj], matplotlib.patches.Rectangle):
                g.get_children()[jj].set_height(g.get_children()[jj].get_height()/(len(mut_table))*8)
        ax.set_xlabel('Cancer Cell Fraction', fontsize=17)
        ax.set_ylabel('Density', fontsize=17)
        ax.legend(bbox_to_anchor=(1.05, 1.05), loc=2, borderaxespad=0., fontsize=14)
        ax.invert_xaxis()
        ax.set_xlim([1.5, 0])
        ax.set_ylim([0, 2])
        ax.set_title(patient_id, fontsize=20)
        plt.savefig('20190828_tcga_plots/{}_{}_clonesig_panel3.pdf'.format(patient_id, folder_type), bbox_inches='tight')


        clone_cols = sns.husl_palette(mut_table.clone.nunique(), l=0.8, s=0.7)
        fig = plt.figure(figsize=(23, 10), dpi= 80)
        grid = plt.GridSpec(4, 8, hspace=0.1, wspace=0.1)

        # Define the axes
        ax_main = fig.add_subplot(grid[:-1, :-1])
        ax_right = fig.add_subplot(grid[:-1, -1], xticklabels=[], yticklabels=[])
        ax_bottom = fig.add_subplot(grid[-1, 0:-1], xticklabels=[], yticklabels=[])
        # Scatterplot on main ax
        sns.scatterplot(y='vaf_purity', x='trinucleotide', data=mut_table[mut_table.vaf_purity<1.5],
                        hue='signature', palette=cols,
                        hue_order=[s for s in selected_sigs if s in mut_table.signature.unique()],
                        s=250, marker='d', ax=ax_main,
                        legend=False)
        #ax_main.invert_yaxis()
        for i in range(96):
            ax_main.axvline(i, ls=':', lw=0.5, color='gray')
        for i in np.arange(-0.5, 96.5, 16):
            ax_main.axvline(i, color='gray', lw=1.2)
        ax_main.set_xlim([-1, 96])
        ax_main.set_ylim([1.5, 0])

        # histogram on the right
        #sns.distplot(mut_table.vaf_purity, kde=False, rug=True, ax=ax_right, color='gray', vertical=True, bins=final_bins)
        for i, present_sig in enumerate(est_sigs):
            
            g=sns.kdeplot(mut_table[mut_table.signature==present_sig].vaf_purity,
                          label=present_sig, ax=ax_right, color=cols[i], vertical=True,
                         lw=4)
        for i, present_sig in enumerate(est_sigs):
            dat = g.lines[i].get_xydata()
            dat[:,0] = dat[:,0]*len(mut_table[mut_table.signature==present_sig])/len(mut_table)
            g.lines[i].set_data(dat.T)
            if len(mut_table[mut_table.signature==present_sig]) < 4:
                dat = g.lines[i].get_xydata()
                dat[:,0] = dat[:,0] / 10
                g.lines[i].set_data(dat.T)

        ax_right.set_xlim([0, 1])

        for i in range(est.J):
            g=sns.distplot(mut_table[mut_table.clone==i].vaf_purity, kde=False,
                           label='clone {}'.format(i), ax=ax_right,
                           color=clone_cols[i], bins=final_bins, vertical=True)
        for jj in range(len(g.get_children())):
            if isinstance(g.get_children()[jj], matplotlib.patches.Rectangle):
                g.get_children()[jj].set_width(g.get_children()[jj].get_width()/(len(mut_table))*6)
        ax_right.set_xlabel('', fontsize=17)
        ax_right.legend(bbox_to_anchor=(1.05, 0.8), loc=2, borderaxespad=0., fontsize=14)
        #ax_right.set_ylim([1.5, 0])


        #ax_right.invert_yaxis()
        ax_right.set_ylim([1.5, 0])
        sns.countplot(x="trinucleotide", data=mut_table, ax=ax_bottom, color='gray')
        ax_bottom.invert_yaxis()
        ax_bottom.set_xlim([-1, 96])
        ax_bottom.set_xticks(range(96))
        ax_bottom.get_xaxis().set_tick_params(direction='out', pad=25, rotation=90)
        _=ax_bottom.set_xticklabels(['{}{}{}'.format(c[0], c[2], c[-1]) for c in PAT_LIST], fontsize=10, color='gray', ha='center', va='bottom')
        _=ax_main.set_xticklabels([])
        _=ax_main.set_xticks([])
        ax_main.set_ylabel('Cancer Cell Fraction', fontsize=23)
        ax_bottom.set_xlabel('Trinucleotide context', fontsize=23)


        ax_main.set_xlim([-1, 96])
        ax_main.set_xticks(range(96))
        ax_main.get_xaxis().set_tick_params(direction='out', pad=25, rotation=90)
        _=ax_main.set_xticklabels(['' for c in PAT_LIST], fontsize=10, color='gray', ha='center', va='bottom')
        ax_main.set_ylabel('Cancer Cell Fraction', fontsize=23)
        ax_main.set_xlabel('', fontsize=23)
        ax_main.tick_params(axis="y", labelsize=15)
        axes2 = ax_main.twiny()   # mirror them
        axes2.set_xlim([-1, 96])
        axes2.set_xticks(np.arange(7.5, 96, 16))
        axes2.get_xaxis().set_tick_params(length=0, pad=10)
        _ = axes2.set_xticklabels(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'], fontsize=18, color='gray', ha='left', rotation=0, va='center')

        _=ax_main.set_xticks([])
        _=ax_right.set_xticks([])
        _=ax_right.set_yticks([])
        _=ax_bottom.set_yticks([])

        axes3 = ax_right.twiny()   # mirror them
        axes3.set_xlabel('Density', fontsize=23)
        ax_right.set_ylabel('', fontsize=23)
        ax_bottom.set_ylabel('Mutation\nCount', fontsize=23)
        _=axes3.set_xticks([])

        plt.savefig('20190828_tcga_plots/{}_{}_clonesig_global.pdf'.format(patient_id, folder_type), bbox_inches='tight')



# # black and white figure

# fig = plt.figure(figsize=(28, 10), dpi= 80)
# grid = plt.GridSpec(5, 5, hspace=0.5, wspace=0.2)

# # Define the axes
# ax_main = fig.add_subplot(grid[:-1, :-1])
# ax_right = fig.add_subplot(grid[:-1, -1], xticklabels=[], yticklabels=[])
# ax_bottom = fig.add_subplot(grid[-1, 0:-1], xticklabels=[], yticklabels=[])

# # Scatterplot on main ax
# sns.scatterplot(y='vaf_purity', x='trinucleotide', data=mut_table,
#                      hue='signature',
#                      hue_order=[s for s in selected_sigs if s in mut_table.signature.unique()],
#                      palette=['gray'] * mut_table.signature.nunique(), s=110, marker='d', ax=ax_main,
#                      legend=False)
# ax_main.invert_yaxis()
# ax_main.set_ylim([1.5, 0])
# # histogram on the right
# sns.distplot(mut_table.vaf_purity, kde=False, rug=True, ax=ax_right, color='gray', vertical=True, bins=40)
# ax_right.invert_yaxis()
# ax_right.set_ylim([1.5, 0])

# #ax_bottom.hist(df.displ, 40, histtype='stepfilled', orientation='vertical', color='deeppink')
# sns.countplot(x="trinucleotide", data=mut_table, ax=ax_bottom)
# ax_bottom.invert_yaxis()
# plt.show(block=False)


# Decorations
# ax_main.set(title=patient_id, xlabel='displ', ylabel='hwy')
# ax_main.title.set_fontsize(20)
# for item in ([ax_main.xaxis.label, ax_main.yaxis.label] + ax_main.get_xticklabels() + ax_main.get_yticklabels()):
#     item.set_fontsize(14)

# xlabels = ax_main.get_xticks().tolist()
# ax_main.set_xticklabels(xlabels)






