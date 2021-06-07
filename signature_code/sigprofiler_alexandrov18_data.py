#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np

'''
download code (2019/05/18)
synapse get -r syn11804040
'''

input_path = 'external_data/sig_profiler_results'

pcawg_inp_filename = 'PCAWG_sigProfiler_SBS_signatures_in_samples.csv'
tcga_inp_filename = 'TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv'
other_wes_inp_filename = 'nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv'
other_wgs_inp_filename = 'nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv'

pcawg_inp = pd.read_csv(input_path + '/' + pcawg_inp_filename)
tcga_inp = pd.read_csv(input_path + '/' + tcga_inp_filename)
other_wes_inp = pd.read_csv(input_path + '/' + other_wes_inp_filename)
other_wgs_inp = pd.read_csv(input_path + '/' + other_wgs_inp_filename)

"""
it is annoying to make the mapping of all the cancer types... I have written
to cosmic to know if they could provide the data in a better format... We will see
For simulations I can use the current table, so let's move forward with that.
"""


"""
other part: get temporary data from the MAthworks Sigprofiler package
downloaded on May 16th 2019
"""
sig_profiler_path = '/Users/JudithAbecassis/Documents/PhD/dl_tools/'
cancertype_sig = sp.io.loadmat(sig_profiler_path + 'SigProfiler_2_5_1_7/SigProfiler_2_5_1_7/tools/SigProfilerSingleSample/input/common/signatures_in_samples_and_cancer_types.mat')
sig_cancer_type = pd.DataFrame(data=cancertype_sig['signaturesInCancerTypes'],
                               index=[i[0] for i in cancertype_sig['cancerTypes'][:,0]],
                               columns=[i[0] for i in cancertype_sig['sigNames'][:,0]])

sig_cancer_type.to_csv('/Users/JudithAbecassis/Documents/PhD/TCGA_signatures/external_data/match_cancer_type_sig_v3.csv')

"""
download of signatures
WGS  signatures syn11738319
WES signatures syn11967914
"""