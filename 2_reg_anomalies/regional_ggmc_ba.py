"""
regional_ggmc.py

Author: M. Zemp
Date: 16 May 2020
Last changes: 19 May 2020

Scripted for Python 3.7

Description:
This script reads glacier-wide mass-balance data from the WGMS FoG database
and provides functions for related analysis and plots.

Input: ggmc_bw-bs-ba.csv (UTF-8 encoding)

Return: tbd.svg
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions_ggmc import *
import scipy  # usage: scipy.stats.funct()

"""main code"""
# define input

path = os.path.dirname(os.path.abspath(__file__))

in_data_file = path + '\\in_data\\' + 'ggmc_bw-bs-ba.csv'
# in_data_file_short = path + '\\in_data\\' + 'vei-on-ggmc_bw-bs-ba_acn.csv'

# read glacier data from csv files into dataframe
input_df = pd.read_csv(in_data_file, delimiter=',', header=0)
print('Input file ({}), data fields: \n {}'.format(in_data_file, list(input_df.columns)))

# number crunching: create unique list of glacier ids and years with data
wgms_id_lst = input_df['WGMS_ID'].unique().tolist()
yr_lst = list(range(min(input_df['YEAR']), max(input_df['YEAR']), 1))
reg_lst = input_df['GLACIER_REGION_CODE'].unique().tolist()

# number crunching: create & export data frame with mass-balance data from input file
ba_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '.csv'
ba_unc_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '_unc.csv'

# create mass-balance data csv if it has not been produced before
ba_df = create_mb_dataframe(input_df, wgms_id_lst, yr_lst, 'ANNUAL_BALANCE')
ba_df.to_csv(ba_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')

ba_unc_df = create_mb_dataframe(input_df, wgms_id_lst, yr_lst, 'ANNUAL_BALANCE_UNC')
ba_unc_df.to_csv(ba_unc_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')
exit()

# read mass-balance data from csv into dataframe (if produced)
ba_df = pd.read_csv(ba_file, delimiter=',', header=0, index_col=0)
ba_df.columns = ba_df.columns.map(int)  # make columns names great again


# number crunching: create list of glacier ids for glacier region categories
# add neighbouring glaciers for regions not represented in the present period.

glac_ALA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ALA'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_WNA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'WNA'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ACN = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ACN'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ACS = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ACS') | (input_df['GLACIER_REGION_CODE'] == 'ACN')), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_GRL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'GRL'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ISL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ISL'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SJM = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'SJM'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SCA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'SCA'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_RUA = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'RUA') | (input_df['GLACIER_REGION_CODE'] == 'SJM')), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASN = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ASN') | (input_df['WGMS_ID'] == 853)), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_CEU = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'CEU'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_CAU = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'CAU'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASC = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ASC'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASE = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ASE'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASW = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ASW'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_TRP = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'TRP'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SAN = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'SAN'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_NZL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'NZL'), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ANT = input_df.loc[(input_df['GLACIER_REGION_CODE'] == 'ANT'), ['GLACIER_REGION_CODE', 'WGMS_ID']] # incl. Ant. Peninsula
# print(glac_ANT['WGMS_ID'].unique().tolist())

# number crunching: select mass-balance data for glacier region groups
# annual mass-balance data
ba_ALA_df = ba_df.loc[:, list(glac_ALA['WGMS_ID'].unique().tolist())]
ba_WNA_df = ba_df.loc[:, list(glac_WNA['WGMS_ID'].unique().tolist())]
ba_ACN_df = ba_df.loc[:, list(glac_ACN['WGMS_ID'].unique().tolist())]
ba_ACS_df = ba_df.loc[:, list(glac_ACS['WGMS_ID'].unique().tolist())]
ba_GRL_df = ba_df.loc[:, list(glac_GRL['WGMS_ID'].unique().tolist())]
ba_ISL_df = ba_df.loc[:, list(glac_ISL['WGMS_ID'].unique().tolist())]
ba_SJM_df = ba_df.loc[:, list(glac_SJM['WGMS_ID'].unique().tolist())]
ba_SCA_df = ba_df.loc[:, list(glac_SCA['WGMS_ID'].unique().tolist())]
ba_RUA_df = ba_df.loc[:, list(glac_RUA['WGMS_ID'].unique().tolist())]
ba_ASN_df = ba_df.loc[:, list(glac_ASN['WGMS_ID'].unique().tolist())]
ba_CEU_df = ba_df.loc[:, list(glac_CEU['WGMS_ID'].unique().tolist())]
ba_CAU_df = ba_df.loc[:, list(glac_CAU['WGMS_ID'].unique().tolist())]
ba_ASC_df = ba_df.loc[:, list(glac_ASC['WGMS_ID'].unique().tolist())]
ba_ASW_df = ba_df.loc[:, list(glac_ASW['WGMS_ID'].unique().tolist())]
ba_ASE_df = ba_df.loc[:, list(glac_ASE['WGMS_ID'].unique().tolist())]
ba_TRP_df = ba_df.loc[:, list(glac_TRP['WGMS_ID'].unique().tolist())]
ba_SAN_df = ba_df.loc[:, list(glac_SAN['WGMS_ID'].unique().tolist())]
ba_NZL_df = ba_df.loc[:, list(glac_NZL['WGMS_ID'].unique().tolist())]
ba_ANT_df = ba_df.loc[:, list(glac_ANT['WGMS_ID'].unique().tolist())]


# number crunching: calculate global mean series of regional averages
ba_reg_means_concat_df = pd.concat((ba_ALA_df.mean(axis=1), ba_WNA_df.mean(axis=1), ba_ACN_df.mean(axis=1), ba_ACS_df.mean(axis=1), ba_GRL_df.mean(axis=1),
                                    ba_ISL_df.mean(axis=1), ba_SJM_df.mean(axis=1), ba_SCA_df.mean(axis=1), ba_RUA_df.mean(axis=1), ba_ASN_df.mean(axis=1),
                                    ba_CEU_df.mean(axis=1), ba_CAU_df.mean(axis=1), ba_ASC_df.mean(axis=1), ba_ASW_df.mean(axis=1), ba_ASE_df.mean(axis=1),
                                    ba_TRP_df.mean(axis=1), ba_SAN_df.mean(axis=1), ba_NZL_df.mean(axis=1), ba_ANT_df.mean(axis=1)))
ba_reg_means_by_yr_index_df = ba_reg_means_concat_df.groupby(ba_reg_means_concat_df.index)
ba_global_moravgs_df = round(ba_reg_means_by_yr_index_df.mean(), 0)

# number crunching: calculate anomalies over reference period
year_ini = 2009
year_fin = 2018

#read years of interest
# dict = [{'year': year_ini, 'vei': 5}, {'year': year_fin, 'vei': 5}]
# yr_df = pd.DataFrame(yr_dict)
# yr_df.set_index('year', inplace=True)

reference_period = range(year_ini, year_fin + 1)
glac_anom = calc_anomalies(ba_SAN_df, reference_period)
# glac_anom['mean'] = glac_anom.mean(axis=1)
reg_anom = glac_anom.mean(axis=1)

print(ba_SAN_df)
print(glac_anom)
# print(reg_anom)
glac_anom.plot()
plt.show()


yr_lst_anom=list(np.arange(1960,2020))
anom_df = pd.DataFrame(index=yr_lst_anom)
anom_df.index.name='year'
anom_df = anom_df.reset_index()

pd.concat([anom_df,new_df],axis=1)
print(anom_df)
exit()




# figures: plot mass-balance time series & volcanic eruptions, no data: [9999]
ba_df_list = [ba_ALA_df.mean(axis=1), ba_WNA_df.mean(axis=1), ba_ACN_df.mean(axis=1), ba_ACS_df.mean(axis=1), ba_GRL_df.mean(axis=1),
                    ba_ISL_df.mean(axis=1), ba_SJM_df.mean(axis=1), ba_SCA_df.mean(axis=1), ba_RUA_df.mean(axis=1), ba_ASN_df.mean(axis=1),
                    ba_CEU_df.mean(axis=1), ba_CAU_df.mean(axis=1), ba_ASC_df.mean(axis=1), ba_ASW_df.mean(axis=1), ba_ASE_df.mean(axis=1),
                    ba_TRP_df.mean(axis=1), ba_SAN_df.mean(axis=1), ba_NZL_df.mean(axis=1), ba_ANT_df.mean(axis=1), ba_global_moravgs_df]


plot_time_series(ba_df_list, bw_df_list, bs_df_list, yr_df)

# figures: plot epoch analysis
plot_epoch_analysis(glac_anom, reg_anom, year_ini, 'Pinatubo & Cerro Hudson')

# tables: export regional mean balances to csv
ba_reg_file = path + '\\out_data\\' + 'ba_regions' + '.csv'


ba_out = pd.concat([ba_ALA_df.mean(axis=1), ba_WNA_df.mean(axis=1), ba_ACN_df.mean(axis=1), ba_ACS_df.mean(axis=1), ba_GRL_df.mean(axis=1),
                    ba_ISL_df.mean(axis=1), ba_SJM_df.mean(axis=1), ba_SCA_df.mean(axis=1), ba_RUA_df.mean(axis=1), ba_ASN_df.mean(axis=1),
                    ba_CEU_df.mean(axis=1), ba_CAU_df.mean(axis=1), ba_ASC_df.mean(axis=1), ba_ASW_df.mean(axis=1), ba_ASE_df.mean(axis=1),
                    ba_TRP_df.mean(axis=1), ba_SAN_df.mean(axis=1), ba_NZL_df.mean(axis=1), ba_ANT_df.mean(axis=1), ba_global_moravgs_df], axis=1)
ba_out = ba_out.round(decimals=0)
ba_out.columns = reg_lst
ba_out.to_csv(ba_reg_file, sep=',', encoding='utf-8', index=True, index_label='Year')



print('.........................................................................................')
print('"If they are well cleaned out, volcanoes burn slowly and steadily, without any eruptions."')
print('Antoine de Saint-Exupéry')
print('.........................................................................................')
