import os
import pandas as pd
import time
import numpy as np
from datetime import datetime 


# get table from the model selection criteria analysis
inp = pd.read_csv('20200213_simu_cn_cancertype_run.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    filename = '{}/clonesig_cn_cancer_type_features.tsv'.format(row.folderpath)
    try:
        timestamp = os.path.getmtime(filename)
        df = pd.read_csv(filename, sep='\t')
        df = df.assign(timestamp=timestamp)
        tables.append(df)
        #print(filename, 'ok')
    except FileNotFoundError:
        print(index, filename)
        non_working.append(index + 1)

big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200222_clonesig_model_selection_criterion.csv', sep='\t', index=False)


# get the table from the calibration of the statistical test
inp = pd.read_csv('20200217_loglikelihood_test_ratio.csv', header=None,
                  names=['idx', 'nb_clones', 'nb_mut', 'cancer_type', 'perc_dip', 'nb_seed'])
tables = list()
output_dir = "20200216_estimate_likelihood_test_statistics"
non_working = list()
for index, row in inp.iterrows():
    pass
    filename = '{}/cancertype{}_nbmut{}_nbclones{}_nbseed{}_percdip{}.csv'.\
        format(output_dir, row['cancer_type'],  row['nb_mut'],
               row['nb_clones'], row['nb_seed'], row['perc_dip'])
    #print(f, filename)
    try:
        timestamp = os.path.getmtime(filename)
        df = pd.read_csv(filename, sep='\t')
        df = df.assign(timestamp=timestamp)
        tables.append(df)
    except FileNotFoundError:
        print(filename, ' is missing')
        non_working.append(index+1)

big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200222_clonesig_model_selection.csv', sep='\t', index=False)


# get table for the main benchmark
inp = pd.read_csv('20200218_simu_cn_cancertype_run.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for method in ('clonesig', 'pyclone', 'sciclone', 'ccube', 'deconstructsigs', 'tracksig', 'palimpsest', 'dpclust', 'phylogicndt', 'tracksigfreq'):
        filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
        try:
            timestamp = os.path.getmtime(filename)
            df = pd.read_csv(filename, sep='\t')
            df = df.assign(timestamp=timestamp)
            tables.append(df)
            #print(filename, 'ok')
        except FileNotFoundError:
            if 'eval_pyclone' in filename:
                print(index, filename)
                non_working.append(row.folderpath)
non_working_df = pd.DataFrame(non_working, columns=['folderpath'])
non_working_df.index = list(range(1, len(non_working) + 1))
non_working_df.to_csv('20200504_to_rerun_evalbug.csv', index=True, header=False)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200520_eval_compare_simulations_new.csv', sep='\t', index=False)

# get table from the main benchmark (constant signature activity)
inp = pd.read_csv('20200229_simu_cn_cancertype_run_cst.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for method in ('clonesig', 'pyclone', 'sciclone', 'ccube', 'deconstructsigs', 'tracksig', 'palimpsest', 'dpclust', 'phylogicndt', 'tracksigfreq'):
        filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
        try:
            timestamp = os.path.getmtime(filename)
            df = pd.read_csv(filename, sep='\t')
            df = df.assign(timestamp=timestamp)
            tables.append(df)
            #print(filename, 'ok')
        except FileNotFoundError:
            if 'eval_clonesig' in filename:
                print(index, filename)
                non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200520_eval_compare_simulations_cst.csv', sep='\t', index=False)

# get table for the analysis of the statistical test sensitivity
inp = pd.read_csv('20200501_clonesig_power_eval.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
inp2 = pd.read_csv('20200505_clonesig_power_eval.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
inp3 = pd.concat([inp, inp2]).reset_index()
tables = list()
non_working = list()
for index, row in inp3.iterrows():
    filename = '{}/eval_power.tsv'.format(row.folderpath)
    try:
        timestamp = os.path.getmtime(filename)
        df = pd.read_csv(filename, sep='\t')
        df = df.assign(timestamp=timestamp)
        tables.append(df)
        #print(filename, 'ok')
    except FileNotFoundError:
        print(index, filename)
        non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200512_eval_clonesig_power_statistical_test_sensitivity.csv', sep='\t', index=False)


# get table for the power of separation of clonesig
inp = pd.read_csv('20200305_clonesig_power_eval.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    method = 'clonesig'
    filename = '{}/eval_power.tsv'.format(row.folderpath)
    try:
        timestamp = os.path.getmtime(filename)
        df = pd.read_csv(filename, sep='\t')
        df = df.assign(timestamp=timestamp)
        tables.append(df)
        #print(filename, 'ok')
    except FileNotFoundError:
        print(index, filename)
        non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200513_eval_clonesig_power_2clones.csv', sep='\t', index=False)


# get table of results on the TCGA
tcga_folder = '20200503_TCGA'
inp = pd.read_csv('TCGA_pancancer_patient_list.csv', header=None, names=['idx', 'patient_id', 'cancer_loc', 'cosmic_type'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    try:
        res_df = pd.read_csv('{}/{}/clonesig_results_restr_curated.csv'.format(tcga_folder, row['patient_id']), sep='\t')
        clonal_sigs = pd.read_csv('{}/{}/clonesig_clonal_sigs_restr_curated.csv'.format(tcga_folder, row['patient_id']), sep='\t')
        subclonal_sigs = pd.read_csv('{}/{}/clonesig_subclonal_sigs_restr_curated.csv'.format(tcga_folder, row['patient_id']), sep='\t')
        clonal_sigs.columns = ['clonal_{}'.format(c) for c in clonal_sigs.columns]
        subclonal_sigs.columns = ['subclonal_{}'.format(c) for c in subclonal_sigs.columns]
        df = pd.concat([res_df, clonal_sigs, subclonal_sigs], axis=1)
        tables.append(df)
        #print(filename, 'ok')
    except FileNotFoundError:
        print(index, row['patient_id'])
        non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)



CANCER_LOCS = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COADREAD', 'DLBC',
               'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG',
               'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD',
               'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS',
               'UVM']
LOC_STAGE_COL = {loc: 'AJCC_PATHOLOGIC_TUMOR_STAGE' for loc in CANCER_LOCS}
LOC_STAGE_COL['CESC'] = 'CLINICAL_STAGE'
LOC_STAGE_COL['DLBC'] = 'CLINICAL_STAGE'
LOC_STAGE_COL['GBM'] = None
LOC_STAGE_COL['LAML'] = None
LOC_STAGE_COL['LGG'] = None
LOC_STAGE_COL['OV'] = 'CLINICAL_STAGE'
LOC_STAGE_COL['PCPG'] = None
LOC_STAGE_COL['PRAD'] = None
LOC_STAGE_COL['SARC'] = None
LOC_STAGE_COL['THYM'] = None
LOC_STAGE_COL['UCEC'] = 'CLINICAL_STAGE'
LOC_STAGE_COL['UCS'] = 'CLINICAL_STAGE'

other_cols_merge = {'age_at_initial_pathologic_diagnosis', 'gender', 'race',
                    'ajcc_pathologic_tumor_stage', 'clinical_stage',
                    'histological_type', 'histological_grade',
                    'initial_pathologic_dx_year', 'menopause_status',
                    'birth_days_to', 'vital_status', 'tumor_status',
                    'last_contact_days_to', 'death_days_to', 'cause_of_death',
                    'new_tumor_event_type', 'new_tumor_event_site',
                    'new_tumor_event_site_other', 'new_tumor_event_dx_days_to',
                    'treatment_outcome_first_course', 'margin_status',
                    'residual_tumor'}


clinical_data_dict = dict()
clinical_cols_dict = dict()
for cancer_loc in CANCER_LOCS:
    # load clinical data
    clinical_data = pd.read_csv(
        'data_tcga_wes/{}/clinical/data_bcr_clinical_data_patient.txt'.format(cancer_loc),
        sep='\t', skiprows=4)
    clinical_data = clinical_data.assign(
        binary_vital_status=clinical_data.apply(
            lambda x: 1 if x.OS_STATUS == 'DECEASED' else
            (0 if x.OS_STATUS == 'LIVING' else np.nan), axis=1))
    try:
        clinical_data = clinical_data[clinical_data.OS_MONTHS != '[Not Available]']
        #clinical_data = clinical_data[clinical_data.DAYS_TO_BIRTH != '[Not Available]']
        clinical_data = clinical_data.assign(survival_days=clinical_data.apply(lambda x: 30.5*float(x['OS_MONTHS']), axis=1))
        #clinical_data = clinical_data.assign(age_at_diagnosis=clinical_data.apply(lambda x: int(np.round(-float(x.DAYS_TO_BIRTH)/365)), axis=1))
        clinical_data = clinical_data.assign(cancer_loc=cancer_loc)
    except ValueError:
        print(cancer_loc, clinical_data.shape)
        clinical_data = clinical_data.assign(cancer_loc=cancer_loc)
        clinical_data = clinical_data.assign(survival_days=np.nan)
    if LOC_STAGE_COL[cancer_loc]:
        clinical_data = clinical_data.assign(
            stage=pd.Categorical(
                clinical_data[LOC_STAGE_COL[cancer_loc]].str.replace('A', '')
                .str.replace('Stage ', '').str.replace('B', '')
                .str.replace('C', ''), categories=['I', 'II', 'III', 'IV'],
                ordered=True))
    else:
        clinical_data = clinical_data.assign(stage=pd.Categorical(
            [np.nan]*len(clinical_data), categories=['I', 'II', 'III', 'IV'],
                ordered=True))
    clinical_data.loc[clinical_data.AGE=='[Not Available]', 'AGE'] = np.nan
    clinical_data = clinical_data.assign(age_group=pd.cut(clinical_data.AGE.astype(float), bins=[0, 39, 49, 59, 69, 150], labels=['<40', '40-50', '50-60', '60-70', '>70']))
    clinical_data_dict[cancer_loc] = clinical_data
    clinical_filling = clinical_data.replace('\[Not .*\]', np.nan, regex=True)
    filled_enough_columns = clinical_filling\
        .columns[clinical_filling.count() /
                 clinical_filling.shape[0] > 0.8]
    clinical_cols_dict[cancer_loc] = set(filled_enough_columns.to_list())
    # pour essayer de merger avec les variants germline, fichier 1-s2.0-S0092867418303635-mmc2.xlsx
    # mais en fait je vais essayer de récupérer des données directement avec
    # l'espoir d'avoir le TCGA barcode pour merger facilement
    # print(cancer_loc, {c.lower() for c in clinical_data.columns}.intersection(other_cols_merge))


common_cols = ['PATIENT_ID', 'binary_vital_status', 'survival_days', 'AGE', 'cancer_loc', 'age_group', 'SEX', 'stage']
clinical_data_merge = pd.concat([clinical_data_dict[loc][common_cols] for loc in CANCER_LOCS], axis=0)
big_merge = pd.merge(big_table2, clinical_data_merge, left_on='patient_id', right_on='PATIENT_ID', how='left')
big_merge.to_csv('20200214_tcga_results_survival_restr.csv', sep='\t', index=False)

big_merge.to_csv('20200205_tcga_results_survival_restr.csv', sep='\t', index=False)


pcawg_folder = '20200510_pcawg'
inp = pd.read_csv('20200510_pcawgs_final_patient_list.csv', header=None, names=['idx', 'patient_id', 'cancer_loc', 'cohort'])
inp2 = pd.read_csv('20200510_pcawgs_tcga_final_patient_list.csv', header=None, names=['idx', 'patient_id', 'cancer_loc', 'cohort'])
inp3 = pd.concat([inp, inp2]).reset_index()
tables = list()
non_working = list()
for index, row in inp3.iterrows():
    try:
        res_df = pd.read_csv('{}/{}_clonesig_results_restr_curated.csv'.format(pcawg_folder, row['patient_id']), sep='\t')
        clonal_sigs = pd.read_csv('{}/{}_clonesig_clonal_sigs_restr_curated.csv'.format(pcawg_folder, row['patient_id']), sep='\t')
        subclonal_sigs = pd.read_csv('{}/{}_clonesig_subclonal_sigs_restr_curated.csv'.format(pcawg_folder, row['patient_id']), sep='\t')
        clonal_sigs.columns = ['clonal_{}'.format(c) for c in clonal_sigs.columns]
        subclonal_sigs.columns = ['subclonal_{}'.format(c) for c in subclonal_sigs.columns]
        df = pd.concat([res_df, clonal_sigs, subclonal_sigs], axis=1)
        tables.append(df)
        #print(filename, 'ok')
    except FileNotFoundError:
        print(index, row['patient_id'])
        non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200514_pcawg_results.csv', sep='\t', index=False)


salcedo_path = 'salcedo_dream_challenge'
inp = pd.read_csv('20200415_dream_run.csv', header=None, names=['idx', 'tumor', 'depth'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    try:
        res_df = pd.read_csv('{}/{}_{}/result_evaluation_dream_new.csv'.format(salcedo_path, row.tumor, row.depth), sep='\t')
        tables.append(res_df)
    except FileNotFoundError:
        print(index+1, row.tumor, row.depth)
        non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20200525_dream_results.csv', sep='\t', index=False)      


phylo500_path = 'PhylogicNDT500'
inp = pd.read_csv('20200830_phylo500_run.csv', header=None, names=['idx', 'samplename'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for suffix in ('cst', 'var'):
        try:
            res_df = pd.read_csv('{}/{}_{}/result_evaluation_phylo500_new.csv'.format(phylo500_path, row.samplename, suffix), sep='\t')
            tables.append(res_df)
        except FileNotFoundError:
            print(index+1, row.samplename, suffix)
            non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20201011_phylo500_results.csv', sep='\t', index=False)      


simclone1000_path = 'SimClone1000'
inp = pd.read_csv('20200916_simclone1000_run.csv', header=None, names=['idx', 'samplename'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for suffix in ('cst', 'var'):
        #input_table = pd.read_csv('{}/{}_{}/input_t.tsv'.format(simclone1000_path, row.idx, suffix), sep='\t')
        try:
            res_df = pd.read_csv('{}/{}_{}/result_evaluation_simclone1000_new.csv'.format(simclone1000_path, row.idx, suffix), sep='\t')
            tables.append(res_df)
        except FileNotFoundError:
            print(index+1, row.idx, suffix)
            non_working.append(index + 1)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20201011_simclone1000_results.csv', sep='\t', index=False)      




