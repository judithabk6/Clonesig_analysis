#!/usr/bin/env python
# coding: utf-8

# # Analysis of Clonesig results on the TCGA
# ## Pancancer analysis


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

pd.options.display.max_columns = 200
phd_folder_path = '/Users/JudithAbecassis/Documents/PhD'
os.chdir('{}/TCGA_signatures'.format(phd_folder_path))
import warnings
warnings.filterwarnings("ignore")



clonesig_res = pd.read_csv('20190824_tcga_results_survival_restr.csv', sep='\t')




sex_colordict = {}
sex_colordict['Female'] = "#e65050"
sex_colordict['Male'] = "#508ae6"

stage_pt_colordict = {}
stage_pt_colordict['T0'] = '#b5f5c6'
stage_pt_colordict['T1'] = "#81f7a0"
stage_pt_colordict['T2'] = "#69c982"
stage_pt_colordict['T3'] = "#498a5a"
stage_pt_colordict['T4'] = "#22402a"

stage_global_colordict = {}
stage_global_colordict['I'] = "#ff52b4"
stage_global_colordict['II'] = "#c23a87"
stage_global_colordict['III'] = "#87295e"
stage_global_colordict['IV'] = "#541b3b"
'#ababab'
'#ffffff'

age_colordict = {}
age_colordict['<40'] = '#cfcfcf'
age_colordict['40-50'] = "#a1a1a1"
age_colordict['50-60'] = "#757575"
age_colordict['60-70'] = "#4f4f4f"
age_colordict['>70'] = "#000000"

kmeans_colordict = {}
colorblind=["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC",
            "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"]
for i in range(len(colorblind)):
    kmeans_colordict[i] = colorblind[i]
# In[31]:


def plot_function_TCGA(loc, input_df, n_clusters):

    letters = list('abcde')
    letter_idx = 0
    sbs_names = [c.split('_')[1] for c in input_df if 'subclonal_SBS' in c]
    bla_protected = pd.DataFrame(input_df[[c for c in input_df if 'subclonal_SBS' in c]].fillna(0).values - input_df[[c for c in input_df if ('clonal_SBS' in c) and ("sub" not in c)]].fillna(0).values, columns=sbs_names)
    relevant_idx = input_df[input_df.pval<=0.05].index
    bla_protected[~bla_protected.index.isin(relevant_idx)] = 0

    relevant_sigs = [c for i, c in enumerate(bla_protected.columns) if np.abs(bla_protected).sum(axis=0)[c]>0.2]

    relevant_sigs_clonal = [c for i, c in enumerate(bla_protected.columns) if input_df['clonal_{}'.format(c)].sum()>0]
    print(len(relevant_sigs), len(relevant_sigs_clonal))

    sns.set_context('poster', font_scale=0.15*len(relevant_sigs)+0.3)

    from sklearn.cluster import KMeans
    clust = KMeans(n_clusters=n_clusters, random_state=0).fit_predict(bla_protected.loc[relevant_idx, relevant_sigs]) + 2

    bla_protected = bla_protected.assign(clust = 0)
    bla_protected.loc[input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].index, 'clust'] = 1
    bla_protected.loc[relevant_idx, 'clust'] = clust
    bla_protected.clust.value_counts()
    bla_protected = bla_protected.assign(SEX = input_df.SEX)
    bla_protected = bla_protected.assign(patient_id = input_df.patient_id)
    if 'age_group' in input_df.columns:
        bla_protected = bla_protected.assign(age_group = input_df.age_group)
    if 'staging_global' in input_df.columns:
        bla_protected = bla_protected.assign(staging_global = input_df.staging_global)
    if 'staging_pt' in input_df.columns:
        bla_protected = bla_protected.assign(staging_pt = input_df.staging_pt)
    bla_protected_s = bla_protected.sort_values('clust')
    
    row_colors = list()
    leg_list = list()
    row_colors.append(bla_protected_s.SEX.map(sex_colordict).fillna('#ababab'))
    if 'age_group' in input_df.columns:
        row_colors.append(bla_protected_s.age_group.map(age_colordict).fillna('#ffffff'))
    if 'staging_global' in input_df.columns:
        row_colors.append(bla_protected_s.staging_global.map(stage_global_colordict).fillna('#ababab'))

    if 'staging_pt' in input_df.columns:
        row_colors.append(bla_protected_s.staging_pt.map(stage_pt_colordict).fillna('#ababab'))

    #row_colors.append(bla_protected_s.clust.map(kmeans_colordict).fillna('#ababab'))


    g=sns.clustermap(bla_protected_s[relevant_sigs], cmap="RdBu_r", vmin=-1, vmax=1,
                     yticklabels=False, row_colors=row_colors,
                    row_cluster=True, col_cluster=False, figsize=(2 * len(relevant_sigs), 2 * len(relevant_sigs)))
    #g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([0.95, .2, .03, .45])
    g.cax.set_ylabel('signature change intensity', position=(0, 0.5), x=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=70)
    g.ax_heatmap.set_xlabel('signatures')
    g.ax_heatmap.set_ylabel('patients', labelpad=40*len(relevant_sigs))
    g.ax_heatmap.yaxis.set_label_position("left")
    g.ax_heatmap.text(-1.5, -1, letters[letter_idx], weight='bold')
    letter_idx = letter_idx + 1

    patchList = []
    patchList.append(Patch(color='white', label='Gender:'))
    for key in sorted(list(set(bla_protected_s.SEX.unique()).intersection(set(sex_colordict.keys())))):
        data_key = Patch(color=sex_colordict[key], label=key)
        patchList.append(data_key)
    patchList.append(Patch(color='#ababab', label='missing'))
    
    my_own_order = ['<40', '40-50', '50-60', '60-70', '>70', 'nan', np.nan]
    age_order = {key: i for i, key in enumerate(my_own_order)}
    if 'age_group' in input_df.columns:
        patchList.append(Patch(color='white', label='Age:'))
        for key in sorted(list(set(bla_protected_s.age_group.unique()).intersection(set(age_colordict.keys()))), key=lambda d: age_order[d]):
            data_key = Patch(color=age_colordict[key], label=key)
            patchList.append(data_key)
        patchList.append(Patch(color='#ffffff', label='missing', edgecolor='#ababab'))

    if 'staging_global' in input_df.columns:
        patchList.append(Patch(color='white', label='Tumor staging:'))
        for key in sorted(list(set(bla_protected_s.staging_global.unique()).intersection(set(stage_global_colordict.keys())))):
            data_key = Patch(color=stage_global_colordict[key], label=key)
            patchList.append(data_key)
        patchList.append(Patch(color='#ababab', label='missing'))
        
    if 'staging_pt' in input_df.columns:
        patchList.append(Patch(color='white', label='Tumor size:'))
        for key in sorted(list(set(bla_protected_s.staging_pt.unique()).intersection(set(stage_pt_colordict.keys())))):
            data_key = Patch(color=stage_pt_colordict[key], label=key)
            patchList.append(data_key)
        patchList.append(Patch(color='#ababab', label='missing'))
        
    # patchList.append(Patch(color='white', label='Cluster:'))
    # for key in sorted(list(set(bla_protected_s.clust.unique()).intersection(set(kmeans_colordict.keys())))):
    #     data_key = Patch(color=kmeans_colordict[key], label=key)
    #     patchList.append(data_key)
    leg_list.append(g.ax_heatmap.legend(handles=patchList, loc=2, bbox_to_anchor=(1.4, 1), fontsize=2.7*len(relevant_sigs)))
    plt.savefig('20190801_paper_figures/{}_heatmap.pdf'.format(loc), bbox_inches='tight')
    
    final_figure_list = list()
    pval_list = list()
    
    # if 'staging_pt' in input_df.columns:
    #     pivot_df = bla_protected_s.pivot_table(index='clust', columns='staging_pt', values='SBS1', aggfunc='count').fillna(0)
    #     chi2, p, dof, ex = chi2_contingency(pivot_df)
    #     pval_list.append(p)
    #     pivot_df.div(pivot_df.sum(axis=1), axis=0).plot.bar(stacked=True, colors=[stage_pt_colordict.get(k, '#ababab') for k in sorted(bla_protected_s.staging_pt.unique())], legend=False)
    #     plt.gca().text(-0.5, 1.1, letters[letter_idx], weight="bold")
    #     letter_idx = letter_idx + 1
    #     plt.savefig('20190801_paper_figures/{}_staging_pt.pdf'.format(loc), bbox_inches='tight')
    #     final_figure_list.append('staging_pt')


    # if 'staging_global' in input_df.columns:
    #     pivot_df = bla_protected_s.pivot_table(index='clust', columns='staging_global', values='SBS1', aggfunc='count').fillna(0)
    #     chi2, p, dof, ex = chi2_contingency(pivot_df)
    #     pval_list.append(p)
    #     pivot_df.div(pivot_df.sum(axis=1), axis=0).plot.bar(stacked=True, colors=[stage_global_colordict.get(k, '#ababab') for k in sorted(bla_protected_s.staging_global.unique())], legend=False)
    #     plt.gca().text(-0.5, 1.1, letters[letter_idx], weight="bold")
    #     letter_idx = letter_idx + 1
    #     plt.savefig('20190801_paper_figures/{}_staging_global.pdf'.format(loc), bbox_inches='tight')
    #     final_figure_list.append('staging_global')


    # if 'age_group' in input_df.columns:
    #     pivot_df = bla_protected_s.pivot_table(index='clust', columns='age_group', values='SBS1', aggfunc='count').fillna(0)
    #     chi2, p, dof, ex = chi2_contingency(pivot_df)
    #     pval_list.append(p)
    #     pivot_df.div(pivot_df.sum(axis=1), axis=0).plot.bar(stacked=True, colors=[age_colordict.get(k, '#ffffff') for k in sorted(bla_protected_s.age_group.unique(), key=lambda d: age_order[d])], legend=False)
    #     plt.gca().text(-0.5, 1.1, letters[letter_idx],weight="bold")
    #     letter_idx = letter_idx + 1
    #     plt.savefig('20190801_paper_figures/{}_age.pdf'.format(loc), bbox_inches='tight')
    #     final_figure_list.append('age')

        
    clonal_sigs = ['clonal_{}'.format(c) for c in relevant_sigs_clonal]
    subclonal_sigs = ['subclonal_{}'.format(c) for c in relevant_sigs_clonal]

    sub_sigs = input_df[clonal_sigs + subclonal_sigs]
    sub_sigs_novar = sub_sigs[(input_df.pval>0.05)&(input_df.nb_clones>=2)]
    sub_sigs[(input_df.pval>0.05)&(input_df.nb_clones>=2)] = np.tile((sub_sigs_novar[clonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].clonal_xi], len(clonal_sigs), axis=0).T + \
        sub_sigs_novar[subclonal_sigs].values * np.repeat([input_df[(input_df.pval>0.05)&(input_df.nb_clones>=2)].largest_subclonal_xi], len(subclonal_sigs), axis=0).T), 2)
    sub_sigs = sub_sigs.assign(SEX = input_df.SEX)
    sub_sigs = sub_sigs.assign(patient_id = input_df.patient_id)
    if 'age_group' in input_df.columns:
        sub_sigs = sub_sigs.assign(age_group = input_df.age_group)
    if 'staging_global' in input_df.columns:
        sub_sigs = sub_sigs.assign(staging_global = input_df.staging_global)
    if 'staging_pt' in input_df.columns:
        sub_sigs = sub_sigs.assign(staging_pt = input_df.staging_pt)
    sub_sigs = pd.merge(sub_sigs, bla_protected[['patient_id', 'clust']])
    sub_sigs.sort_values('clust', inplace=True)
    
    row_colors = list()
    leg_list = list()
    row_colors.append(sub_sigs.SEX.map(sex_colordict).fillna('#ababab'))
    if 'age_group' in input_df.columns:
        row_colors.append(sub_sigs.age_group.map(age_colordict).fillna('#ffffff'))
    if 'staging_global' in input_df.columns:
        row_colors.append(sub_sigs.staging_global.map(stage_global_colordict).fillna('#ababab'))
    if 'staging_pt' in input_df.columns:
        row_colors.append(sub_sigs.staging_pt.map(stage_pt_colordict).fillna('#ababab'))
    #row_colors.append(bla_protected_s.clust.map(kmeans_colordict).fillna('#ababab'))
        
    g=sns.clustermap(sub_sigs[clonal_sigs + subclonal_sigs].fillna(0), cmap="Greens", vmin=0, vmax=1,
                     yticklabels=False, row_colors=row_colors,
                    row_cluster=True, col_cluster=False, figsize=(4 * len(relevant_sigs_clonal), 2 * len(relevant_sigs_clonal)))
    #g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([0.95, .2, .03, .45])
    g.cax.set_ylabel('signature intensity', position=(0, 0.5), x=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=70, ha='right', rotation_mode="anchor")
    g.ax_heatmap.set_xticklabels([item.get_text().split('_')[1] for item in g.ax_heatmap.get_xticklabels()])
    g.ax_heatmap.set_xlabel('signatures')
    g.ax_heatmap.set_ylabel('patients', labelpad=70*len(relevant_sigs_clonal))
    g.ax_heatmap.yaxis.set_label_position("left")
    g.ax_heatmap.axvline(x=len(clonal_sigs), color='k')
    g.ax_heatmap.text(len(clonal_sigs)/3, -2, 'clonal', fontsize=6*len(relevant_sigs_clonal))
    g.ax_heatmap.text(len(clonal_sigs)*(1+1/3), -2, 'subclonal', fontsize=6*len(relevant_sigs_clonal))
    g.ax_heatmap.text(-1.5, -1, letters[letter_idx], fontsize=6*len(relevant_sigs_clonal), weight="bold")
    letter_idx = letter_idx + 1
    plt.savefig('20190801_paper_figures/{}_heatmap_absolute_values.pdf'.format(loc), bbox_inches='tight')
    
    nb_patients = len(input_df)
    nb_patients_sig = len(input_df[input_df.pval<=0.05])
    text_dict = {'age': 'Panel {} represents the repartition of age groups in the different clusters (Chi-square test of independence p={}).',
         'staging_global': 'Panel {} represents the repartition of the clinical tumor stage in the different clusters (Chi-square test of independence p={}).',
         'staging_pt': 'Panel {} represents the repartition of tumor size classes in the different clusters (Chi-square test of independence p={}).'}

    print(r'\begin{figure}')
    print(r'\centering')
    print(r'\includegraphics[height=0.4\textheight]{{figures/{}_heatmap.pdf}}'.format(loc))
    cap_text = 'Panel a: Stratification of patients depending on their pattern of signature change for {} patients ({} patients, including {} with a significant signature change). The heatmap represents the difference between the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF).'.format(loc, nb_patients, nb_patients_sig)
    # for i, fig_name in enumerate(final_figure_list):
    #     print(r'\includegraphics[height=0.15\textheight]{{figures/{}_{}.pdf}}'.format(loc, fig_name))
    #     cap_text += text_dict[fig_name].format(letters[i+1], np.round(pval_list[i], 4))
    # print('\n')
    print(r'\includegraphics[height=0.35\textheight]{{figures/{}_heatmap_absolute_values.pdf}}'.format(loc))   
    cap_text += 'Panel {}: Stratification of patients depending on their complete pattern of signature exposure. The heatmap represents the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF).'.format(letters[len(final_figure_list)+1], loc, nb_patients, nb_patients_sig)
    print('\caption{{{}}}'.format(cap_text))
    print('\label{{{}}}'.format(loc.lower()))
    print(print('\end{figure}\n'))
    
    

#\caption{Stratification of patients depending on their pattern of signature change for UVM patients (80 patients, including 3 with a significant signature change). The heatmap represents the difference between the signature activity in the largest subclone (in terms of number of mutations) and the clonal mutations (defined as belonging to the clone of highest CCF). The two lower panels represent the repartition of tumor size classes (left, Chi-square test of independence p=0.55), and clinical tumor stage (right, Chi-square test of independence p=0.58) in the different clusters.}



# In[32]:


clonesig_res_acc = clonesig_res[clonesig_res.cancer_loc_x=='ACC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/ACC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.PATH_T_STAGE.str[:2])
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
sub_protected_chg = clonesig_res_acc[(clonesig_res_acc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[33]:


sns.set_style('white')
sns.set_context('poster')
plot_function_TCGA('ACC', sub_protected_chg_m, 2)



# ### BLCA

# In[209]:


clonesig_res_blca = clonesig_res[clonesig_res.cancer_loc_x=='BLCA'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/BLCA/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
#clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])

clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
sub_protected_chg = clonesig_res_blca[(clonesig_res_blca.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[210]:


plot_function_TCGA('BLCA', sub_protected_chg_m, 3)




# ### BRCA

# In[211]:


clonesig_res_brca = clonesig_res[clonesig_res.cancer_loc_x=='BRCA'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/BRCA/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_brca[(clonesig_res_brca.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[212]:


plot_function_TCGA('BRCA', sub_protected_chg_m, 3)



# ### CESC

# In[221]:


clonesig_res_cesc = clonesig_res[clonesig_res.cancer_loc_x=='CESC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/CESC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.CLINICAL_STAGE.str.replace('B', '').str.replace('Stage ', '').str.replace('A', '').str.replace('1', '').str.replace('2', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_cesc[(clonesig_res_cesc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')



plot_function_TCGA('CESC', sub_protected_chg_m, 3)



# ### CHOL

# In[219]:


clonesig_res_chol = clonesig_res[clonesig_res.cancer_loc_x=='CHOL'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/CHOL/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_chol[(clonesig_res_chol.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[220]:


plot_function_TCGA('CHOL', sub_protected_chg_m, 1)






# ### COADREAD

# In[223]:


clonesig_res_coadread = clonesig_res[clonesig_res.cancer_loc_x=='COADREAD'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/COADREAD/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_coadread[(clonesig_res_coadread.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[224]:


plot_function_TCGA('COADREAD', sub_protected_chg_m, 4)






# ### DLBC

# In[225]:


clonesig_res_dlbc = clonesig_res[clonesig_res.cancer_loc_x=='DLBC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/DLBC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
#clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.CLINICAL_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_dlbc[(clonesig_res_dlbc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[226]:


plot_function_TCGA('DLBC', sub_protected_chg_m, 2)




# ### ESCA

# In[227]:


clonesig_res_esca = clonesig_res[clonesig_res.cancer_loc_x=='ESCA'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/ESCA/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_esca[(clonesig_res_esca.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[228]:


plot_function_TCGA('ESCA', sub_protected_chg_m, 4)






# ### GBM

# In[229]:


clonesig_res_gbm = clonesig_res[clonesig_res.cancer_loc_x=='GBM'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/GBM/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_gbm[(clonesig_res_gbm.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')
sub_protected_chg = clonesig_res_gbm[(clonesig_res_gbm.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[230]:


plot_function_TCGA('GBM', sub_protected_chg_m, 3)




# ### HNSC

# In[232]:

clonesig_res_hnsc = clonesig_res[clonesig_res.cancer_loc_x=='HNSC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/HNSC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_hnsc[(clonesig_res_hnsc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[233]:


plot_function_TCGA('HNSC', sub_protected_chg_m, 4)







# ### KICH

# In[59]:


clonesig_res_kich = clonesig_res[clonesig_res.cancer_loc_x=='KICH'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/KICH/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_kich[(clonesig_res_kich.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[60]:


plot_function_TCGA('KICH', sub_protected_chg_m, 3)







# ### KIRC

# In[235]:


clonesig_res_kirc = clonesig_res[clonesig_res.cancer_loc_x=='KIRC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/KIRC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE, bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_kirc[(clonesig_res_kirc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[236]:


plot_function_TCGA('KIRC', sub_protected_chg_m, 3)







# ### KIRP

# In[61]:


clonesig_res_kirp = clonesig_res[clonesig_res.cancer_loc_x=='KIRP'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/KIRP/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_kirp[(clonesig_res_kirp.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[62]:


plot_function_TCGA('KIRP', sub_protected_chg_m, 3)




# ### LGG

# In[240]:


clonesig_res_lgg = clonesig_res[clonesig_res.cancer_loc_x=='LGG'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/LGG/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
#clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
#clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_lgg[(clonesig_res_lgg.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[241]:


plot_function_TCGA('LGG', sub_protected_chg_m, 1)


# In[ ]:





# ### LIHC

# In[242]:


clonesig_res_lihc = clonesig_res[clonesig_res.cancer_loc_x=='LIHC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/LIHC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_lihc[(clonesig_res_lihc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[243]:


plot_function_TCGA('LIHC', sub_protected_chg_m, 4)






# ### LUAD

# In[244]:


clonesig_res_luad = clonesig_res[clonesig_res.cancer_loc_x=='LUAD'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/LUAD/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_luad[(clonesig_res_luad.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[245]:


plot_function_TCGA('LUAD', sub_protected_chg_m, 4)







# ### LUSC

# In[63]:


clonesig_res_lusc = clonesig_res[clonesig_res.cancer_loc_x=='LUSC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/LUSC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_lusc[(clonesig_res_lusc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[64]:


plot_function_TCGA('LUSC', sub_protected_chg_m, 4)







# ### MESO

# In[251]:


clonesig_res_meso = clonesig_res[clonesig_res.cancer_loc_x=='MESO'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/MESO/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_meso[(clonesig_res_meso.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[252]:


plot_function_TCGA('MESO', sub_protected_chg_m, 1)






# ### OV

# In[257]:


clonesig_res_ov = clonesig_res[clonesig_res.cancer_loc_x=='OV'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/OV/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
#clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.CLINICAL_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_ov[(clonesig_res_ov.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[258]:


plot_function_TCGA('OV', sub_protected_chg_m, 4)






# ### PAAD

# In[259]:


clonesig_res_paad = clonesig_res[clonesig_res.cancer_loc_x=='PAAD'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/PAAD/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_paad[(clonesig_res_paad.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[260]:


plot_function_TCGA('PAAD', sub_protected_chg_m, 2)






# ### PCPG

# In[65]:


clonesig_res_pcpg = clonesig_res[clonesig_res.cancer_loc_x=='PCPG'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/PCPG/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
#clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
#clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_pcpg[(clonesig_res_pcpg.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[66]:


plot_function_TCGA('PCPG', sub_protected_chg_m, 1)


# ### PRAD

# In[278]:


clonesig_res_prad = clonesig_res[clonesig_res.cancer_loc_x=='PRAD'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/PRAD/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.CLIN_T_STAGE.str[:2])
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_prad[(clonesig_res_prad.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[279]:


plot_function_TCGA('PRAD', sub_protected_chg_m, 2)


# In[ ]:





# ### SARC

# In[6]:


clonesig_res_sarc = clonesig_res[clonesig_res.cancer_loc_x=='SARC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/SARC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_sarc[(clonesig_res_sarc.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[10]:


plot_function_TCGA('SARC', sub_protected_chg_m, 3)


# In[ ]:





# ### SKCM

# In[11]:


clonesig_res_skcm = clonesig_res[clonesig_res.cancer_loc_x=='SKCM'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/SKCM/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_skcm[(clonesig_res_skcm.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[13]:


plot_function_TCGA('SKCM', sub_protected_chg_m, 4)






# ### STAD

# In[14]:


clonesig_res_stad = clonesig_res[clonesig_res.cancer_loc_x=='STAD'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/STAD/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_stad[(clonesig_res_stad.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[15]:


plot_function_TCGA('STAD', sub_protected_chg_m, 3)



# ### TGCT
# 

# In[19]:


clonesig_res_tgct = clonesig_res[clonesig_res.cancer_loc_x=='TGCT'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/TGCT/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_tgct[(clonesig_res_tgct.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[21]:


plot_function_TCGA('TGCT', sub_protected_chg_m, 1)


# In[ ]:





# ### THCA

# In[25]:


clonesig_res_thca = clonesig_res[clonesig_res.cancer_loc_x=='THCA'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/THCA/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_thca[(clonesig_res_thca.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[26]:


plot_function_TCGA('THCA', sub_protected_chg_m, 1)



# ### THYM

# In[29]:


clonesig_res_thym = clonesig_res[clonesig_res.cancer_loc_x=='THYM'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/THYM/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_thym[(clonesig_res_thym.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[31]:


plot_function_TCGA('THYM', sub_protected_chg_m, 2)



# ### UCEC

# In[34]:


clonesig_res_ucec = clonesig_res[clonesig_res.cancer_loc_x=='UCEC'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/UCEC/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_global=clinical_data.CLINICAL_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_ucec[(clonesig_res_ucec.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[35]:


plot_function_TCGA('UCEC', sub_protected_chg_m, 3)



# ### UCS

# In[43]:


clonesig_res_ucs = clonesig_res[clonesig_res.cancer_loc_x=='UCS'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/UCS/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_global=clinical_data.CLINICAL_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', '').str.replace('2', '').str.replace('1', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_ucs[(clonesig_res_ucs.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[44]:


plot_function_TCGA('UCS', sub_protected_chg_m, 1)



# ### UVM

# In[45]:


clonesig_res_uvm = clonesig_res[clonesig_res.cancer_loc_x=='UVM'].reset_index()
clinical_data = pd.read_csv('data_tcga_wes/UVM/clinical/data_bcr_clinical_data_patient.txt', sep='\t', skiprows=4)
clinical_data = clinical_data.assign(staging_pt=clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
clinical_data = clinical_data.assign(staging_global=clinical_data.AJCC_PATHOLOGIC_TUMOR_STAGE.str.replace('A', '').str.replace('Stage ', '').str.replace('B', '').str.replace('C', ''))
clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
sub_protected_chg = clonesig_res_uvm[(clonesig_res_uvm.mutation_set=='protected')]#&(clonesig_res.nb_mut>200)&(clonesig_res.purity>0.4)]
sub_protected_chg_m = pd.merge(sub_protected_chg, clinical_data, how='left', left_on='patient_id', right_on='PATIENT_ID')


# In[48]:


plot_function_TCGA('UVM', sub_protected_chg_m, 1)

clonesig_res.pivot_table(index='cancer_loc_x', columns='mutation_set', values='nb_mut', agg_func='mean')
clonesig_res.groupby('cancer_loc_x').agg({'a':['sum', 'max'], 
                         'b':'mean', 
                         'c':'sum', 
                         'd': lambda x: x.max() - x.min()})
print(clonesig_res.pivot_table(index='cancer_loc_x', columns='mutation_set', values='nb_mut',aggfunc=['mean', 'std', 'count']).to_latex(float_format="{:0.2f}".format)) 

clonesig_res = clonesig_res.assign(binary_vital_status_15y=clonesig_res.binary_vital_status)
clonesig_res = clonesig_res.assign(survival_months_15y=(clonesig_res.survival_days/30.5))

clonesig_res.loc[(clonesig_res.survival_days>=15*365)&(clonesig_res.binary_vital_status==1), 'binary_vital_status_15y'] = 0
clonesig_res.loc[clonesig_res.survival_days>=15*365, 'survival_months_15y'] = 15*365/30.5

print(clonesig_res.pivot_table(index='cancer_loc_x', columns='mutation_set', values='survival_months_15y',aggfunc=['median']).to_latex(float_format="{:0.2f}".format)) 
print(clonesig_res.pivot_table(index='cancer_loc_x', columns='mutation_set', values='binary_vital_status_15y',aggfunc=['sum', 'count']).to_latex(float_format="{:0.2f}".format)) 


clonesig_res.groupby('cancer_loc_x').survival_days.count()

fig, ax = plt.subplots(figsize=(13, 8))
sns.boxplot(x='cancer_loc_x', y='survival_months_15y', data=clonesig_res[(clonesig_res.nb_mut>0)&(clonesig_res.mutation_set=='protected')], ax=ax, palette='gist_ncar') 
ax.xaxis.set_tick_params(rotation=70)
ax.set_xlabel('cancer type')
ax.set_ylabel('followup (in months)')
plt.savefig('20190801_paper_figures/followup_time_distribution.pdf', bbox_inches='tight')
