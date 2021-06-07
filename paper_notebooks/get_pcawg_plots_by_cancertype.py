import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import collections
from scipy.stats import chi2_contingency
from pathlib import Path
import warnings

pd.options.display.max_columns = 200
output_path = '20200515_paper_figures/PCAWG_results/heatmaps'
Path(output_path).mkdir(parents=True, exist_ok=True)

warnings.filterwarnings("ignore")

clonesig_res = pd.read_csv('20200514_pcawg_results.csv', sep='\t')

def plot_function_PCAWG(loc, input_df):

    letters = list('abcde')
    letter_idx = 0
    sbs_names = [c.split('_')[1] for c in input_df if 'subclonal_SBS' in c]
    bla_protected = pd.DataFrame(input_df[[c for c in input_df if 'subclonal_SBS' in c]].fillna(0).values - input_df[[c for c in input_df if ('clonal_SBS' in c) and ("sub" not in c)]].fillna(0).values, columns=sbs_names)
    relevant_idx = input_df[input_df.pval<=0.05].index
    bla_protected[~bla_protected.index.isin(relevant_idx)] = 0
    relevant_sigs = [c for i, c in enumerate(bla_protected.columns) if np.abs(bla_protected).sum(axis=0)[c]>0]
    relevant_sigs_clonal = [c for i, c in enumerate(bla_protected.columns) if input_df['clonal_{}'.format(c)].sum()>0]

    sns.set_context('poster', font_scale=0.15*len(relevant_sigs)+0.3)

    bla_protected = bla_protected.assign(patient_id = input_df.patient_id)
    #row_colors.append(bla_protected_s.clust.map(kmeans_colordict).fillna('#ababab'))

    g=sns.clustermap(bla_protected[relevant_sigs], cmap="RdBu_r", vmin=-1, vmax=1,
                     yticklabels=False,
                     row_cluster=True, col_cluster=False, figsize=(2 * len(relevant_sigs), 2 * len(relevant_sigs)))
    #g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([0.95, .2, .03, .45])
    g.cax.set_ylabel('signature change intensity', position=(0, 0.5), x=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=70)
    g.ax_heatmap.set_xlabel('signatures')
    g.ax_heatmap.set_ylabel('patients', labelpad=25*len(relevant_sigs))
    g.ax_heatmap.yaxis.set_label_position("left")
    g.ax_heatmap.text(-1.5, -1, letters[letter_idx], weight='bold')
    letter_idx = letter_idx + 1

    plt.savefig('{}/{}_heatmap.pdf'.format(output_path, loc.replace('/', '_')), bbox_inches='tight')

    final_figure_list = list()
    pval_list = list()

    clonal_sigs = ['clonal_{}'.format(c) for c in relevant_sigs_clonal]
    subclonal_sigs = ['subclonal_{}'.format(c) for c in relevant_sigs_clonal]

    sub_sigs = input_df[clonal_sigs + subclonal_sigs]
    sub_sigs_novar = sub_sigs[(input_df.pval>0.05)&(input_df.nb_clones>=2)]
    if len(sub_sigs_novar)==1:
        sub_sigs[(input_df.pval>0.05)&(input_df.nb_clones>=2)] = np.tile((sub_sigs_novar[clonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].clonal_xi.values.reshape(-1,1)], len(clonal_sigs), axis=0).T + \
            sub_sigs_novar[subclonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].largest_subclonal_xi.values.reshape(-1,1)], len(subclonal_sigs), axis=0).T), 2)
    else:
        sub_sigs[(input_df.pval>0.05)&(input_df.nb_clones>=2)] = np.tile((sub_sigs_novar[clonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].clonal_xi], len(clonal_sigs), axis=0).T + \
            sub_sigs_novar[subclonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].largest_subclonal_xi], len(subclonal_sigs), axis=0).T), 2)

    sub_sigs = sub_sigs.assign(patient_id = input_df.patient_id)

    g = sns.clustermap(sub_sigs[clonal_sigs + subclonal_sigs].fillna(0), cmap="Greens", vmin=0, vmax=1,
                       yticklabels=False,
                       row_cluster=True, col_cluster=False, figsize=(4 * len(relevant_sigs_clonal), 2 * len(relevant_sigs_clonal)))

    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([0.95, .2, .03, .45])
    g.cax.set_ylabel('signature intensity', position=(0, 0.5), x=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=70, ha='right', rotation_mode="anchor")
    g.ax_heatmap.set_xticklabels([item.get_text().split('_')[1] for item in g.ax_heatmap.get_xticklabels()])
    g.ax_heatmap.set_xlabel('signatures')
    g.ax_heatmap.set_ylabel('patients', labelpad=35*len(relevant_sigs_clonal))
    g.ax_heatmap.yaxis.set_label_position("left")
    g.ax_heatmap.axvline(x=len(clonal_sigs), color='k')
    g.ax_heatmap.text(len(clonal_sigs)/3, -2, 'clonal', fontsize=6*len(relevant_sigs_clonal))
    g.ax_heatmap.text(len(clonal_sigs)*(1+1/3), -2, 'subclonal', fontsize=6*len(relevant_sigs_clonal))
    g.ax_heatmap.text(-1.5, -1, letters[letter_idx], fontsize=6*len(relevant_sigs_clonal), weight="bold")
    letter_idx = letter_idx + 1
    plt.savefig('{}/{}_heatmap_absolute_values.pdf'.format(output_path, loc.replace('/', '_')), bbox_inches='tight')

    nb_patients = len(input_df)
    nb_patients_sig = len(input_df[input_df.pval<=0.05])

    print(r'\begin{figure}')
    print(r'\centering')
    print(r'\includegraphics[height=0.4\textheight]{{{}/{}_heatmap.pdf}}'.format(output_path, loc.replace('/', '_')))
    cap_text = '\\textbf{{Overview of CloneSig results for the PCAWG sub-cohort {} (n={}).}} Panel a: Stratification of patients depending on their pattern of signature change for {} patients ({} patients, including {} with a significant signature change). The heatmap represents the difference between the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF).'.format(loc, nb_patients, loc, nb_patients, nb_patients_sig)
    # for i, fig_name in enumerate(final_figure_list):
    #     print(r'\includegraphics[height=0.15\textheight]{{figures/{}_{}.pdf}}'.format(loc, fig_name))
    #     cap_text += text_dict[fig_name].format(letters[i+1], np.round(pval_list[i], 4))
    # print('\n')
    print(r'\includegraphics[height=0.35\textheight]{{{}/{}_heatmap_absolute_values.pdf}}'.format(output_path, loc.replace('/', '_')))   
    cap_text += 'Panel {}: Stratification of patients depending on their complete pattern of signature exposure. The heatmap represents the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF).'.format(letters[len(final_figure_list)+1], loc, nb_patients, nb_patients_sig)
    print('\caption{{{}}}'.format(cap_text))
    print('\label{{{}}}'.format(loc.lower()))
    print('\end{figure}\n')


#\caption{Stratification of patients depending on their pattern of signature change for UVM patients (80 patients, including 3 with a significant signature change). The heatmap represents the difference between the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF). The two lower panels represent the repartition of tumor size classes (left, Chi-square test of independence p=0.55), and clinical tumor stage (right, Chi-square test of independence p=0.58) in the different clusters.}

match_cohort = {'BRCA': 'Breast', 'COADREAD': 'ColoRect-AdenoCA',
                'LIHC': 'Liver-HCC', 'THCA': 'Thy-AdenoCA',
                'LUSC': 'Lung-SCC', 'KICH': 'Kidney-ChRCC',
                'UCEC': 'Uterus-AdenoCA', 'HNSC': 'Head-SCC',
                'OV': 'Ovary-AdenoCA', 'KIRC': 'Kidney-RCC',
                'SKCM': 'Skin-Melanoma', 'STAD': 'Stomach-AdenoCA',
                'KIRP': 'Kidney-PRCC', 'LUAD': 'Lung-AdenoCA',
                'GBM': 'CNS-GBM', 'SARC': 'Bone-Leiomyo',
                'BLCA': 'Bladder-TCC', 'PRAD': 'Prost-AdenoCA',
                'CESC': 'Cervix-SCC', 'LGG': 'CNS-Oligo',
                'DLBC': 'Lymph-BNHL'}

clonesig_res = clonesig_res.assign(new_mutation_set=clonesig_res.apply(lambda x: match_cohort[x['mutation_set']] if x['mutation_set'] in match_cohort else x['mutation_set'], axis=1))


sns.set_style('white')
sns.set_context('poster')
for ms in sorted(clonesig_res.new_mutation_set.unique()):
    sub_res = clonesig_res[clonesig_res.new_mutation_set==ms].reset_index()
    plot_function_PCAWG(ms, sub_res)


