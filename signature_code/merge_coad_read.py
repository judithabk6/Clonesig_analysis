import pandas as pd
import numpy as np
import os

try:
    rows, columns = os.popen('stty size', 'r').read().split()
    pd.set_option('display.width', int(columns))
    pd.set_option('display.max_columns', 600)
    pd.set_option('display.max_rows', 1000)
except:
    print("running on server, otherwise please investigate")

read_public = pd.read_csv('tmp/READ_useful_final_public_merge_cnv_purity.csv', sep='\t')
read_protected = pd.read_csv('tmp/READ_useful_final_qc_merge_cnv_purity.csv', sep='\t')
coad_public = pd.read_csv('tmp/COAD_useful_final_public_merge_cnv_purity.csv', sep='\t')
coad_protected = pd.read_csv('tmp/COAD_useful_final_qc_merge_cnv_purity.csv', sep='\t')

coadread_public = pd.concat([read_public, coad_public], axis=0)
coadread_public.to_csv('tmp/COADREAD_useful_final_public_merge_cnv_purity.csv', sep='\t', index=False)
coadread_protected = pd.concat([read_protected, coad_protected], axis=0)
coadread_protected.to_csv('tmp/COADREAD_useful_final_qc_merge_cnv_purity.csv', sep='\t', index=False)


patient_list_public = coadread_public.patient_id.unique().tolist()
patient_list_protected = coadread_protected.patient_id.unique().tolist()
intersection_patients = list(set(patient_list_public).intersection(set(patient_list_protected)))
udf = pd.DataFrame(intersection_patients, columns=['patient_id'], index=range(1, len(intersection_patients)+1))
udf = udf.assign(cancer_loc='COADREAD')
udf.to_csv('COADREAD_official_patient_list.csv', header=False)

# in bash console
# cd data_tcga_wes/
# mkdir COADREAD
# mv ../COAD/clinical/* clinical/
