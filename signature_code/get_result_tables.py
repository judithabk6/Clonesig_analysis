import os
import pandas as pd
import time

# get table from the model selection criteria analysis
inp = pd.read_csv('20190528_simu_cn_cancertype_run.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
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
# fix some little error in code
cancer_type_sig_filename = 'external_data/match_cancer_type_sig_v3.csv'
cancer_type_sig = pd.read_csv(cancer_type_sig_filename)
nb_sig_df = cancer_type_sig.drop(cancer_type_sig.columns[0], 1)\
                                .sum(axis=1).to_frame(name='corrected_nb_sig')
big_table2_corr = pd.merge(big_table2,nb_sig_df,
                           left_on='cancer_type', right_index=True)
big_table2_corr.loc[big_table2_corr.index==1, 'nb_sig'] = big_table2_corr.loc[big_table2_corr.index==1, 'corrected_nb_sig'] 
big_table2_corr_tosave = big_table2_corr.drop(big_table2_corr.columns[-1], 1)

big_table2_corr_tosave.to_csv('20190611_clonesig_model_selection_criterion.csv', sep='\t', index=False)


# get the table from the calibration of the statistical test
inp = pd.read_csv('20190610_loglikelihood_test_ratio.csv', header=None,
                  names=['idx', 'nb_clones', 'nb_mut', 'cancer_type', 'perc_dip', 'nb_seed'])
tables = list()
output_dir = "20190610_estimate_likelihood_test_statistics"
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
big_table2.to_csv('20190613_clonesig_model_selection.csv', sep='\t', index=False)


# get table for the main benchmark
inp = pd.read_csv('20190718_simu_cn_cancertype_run.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for method in ('clonesig', 'pyclone', 'sciclone', 'ccube', 'deconstructsigs', 'tracksig', 'palimpsest'):
        filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
        try:
            timestamp = os.path.getmtime(filename)
            df = pd.read_csv(filename, sep='\t')
            df = df.assign(timestamp=timestamp)
            tables.append(df)
            #print(filename, 'ok')
        except FileNotFoundError:
            if 'pyclone' in filename:
                print(index, filename)
                non_working.append(row.folderpath)
non_working_df = pd.DataFrame(non_working, columns=['folderpath'])
non_working_df.index = list(range(1, len(non_working) + 1))
non_working_df.to_csv('20190824_to_rerun_pyclone.csv', index=True, header=False)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('20190824_eval_compare_simulations_new.csv', sep='\t', index=False)

# get table from the main benchmark (constant signature activity)
inp = pd.read_csv('20190729_simu_cn_cancertype_run_cst.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    for method in ('clonesig', 'pyclone', 'sciclone', 'ccube', 'deconstructsigs', 'tracksig', 'palimpsest'):
        filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
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
big_table2.to_csv('201907824_eval_compare_simulations_cst.csv', sep='\t', index=False)

# get table for the analysis of the statistical test sensitivity
inp = pd.read_csv('20190804_clonesig_power_eval.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    method = 'clonesig'
    filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
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
big_table2.to_csv('201907824_eval_clonesig_power_statistical_test_sensitivity.csv', sep='\t', index=False)


# get table for the power of separation of clonesig
inp = pd.read_csv('20190803_clonesig_power_eval.csv', header=None, names=['idx', 'folderpath', 'nb_mut', 'fito'])
tables = list()
non_working = list()
for index, row in inp.iterrows():
    method = 'clonesig'
    filename = '{}/eval_{}.tsv'.format(row.folderpath, method)
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
big_table2.to_csv('201907805_eval_clonesig_power_2clones.csv', sep='\t', index=False)


# get table of results on the TCGA
tcga_folder = '20190704_TCGA'
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
#big_table2.to_csv('20190707_tcga_results_sigprofiler.csv', sep='\t', index=False)


CANCER_LOCS = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COADREAD', 'DLBC',
               'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG',
               'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD',
               'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS',
               'UVM']
clinical_data_dict = dict()
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

    clinical_data_dict[cancer_loc] = clinical_data

common_cols = ['PATIENT_ID', 'binary_vital_status', 'survival_days', 'AGE', 'cancer_loc']
clinical_data_merge = pd.concat([clinical_data_dict[loc][common_cols] for loc in CANCER_LOCS], axis=0)
big_merge = pd.merge(big_table2, clinical_data_merge, left_on='patient_id', right_on='PATIENT_ID', how='left')
big_merge.to_csv('20190824_tcga_results_survival_restr.csv', sep='\t', index=False)












