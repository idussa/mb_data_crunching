"""
vei_on_ggmc.py

Author: M. Zemp
Date: 16 May 2020
Last changes: 19 May 2020

Scripted for Python 3.7

Description:
This script reads glacier-wide mass-balance data from the WGMS FoG database
and provides functions for related analysis and plots.

Input: vei-on-ggmc_bw-bs-ba.csv (UTF-8 encoding)

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

# number crunching: create & export data frame with mass-balance data from input file
ba_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '.csv'
bw_file = path + '\\in_data\\' + 'ggmc_' + 'bw' + '.csv'
bs_file = path + '\\in_data\\' + 'ggmc_' + 'bs' + '.csv'

# create mass-balance data csv if it has not been produced before
ba_df = create_mb_dataframe(input_df, wgms_id_lst, 'ANNUAL_BALANCE')
ba_df.to_csv(ba_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')
# bw_df = create_mb_dataframe(input_df, wgms_id_lst, 'WINTER_BALANCE')
# bw_df.to_csv(bw_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')
# bs_df = create_mb_dataframe(input_df, wgms_id_lst, 'SUMMER_BALANCE')
# bs_df.to_csv(bs_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')

# read mass-balance data from csv into dataframe (if produced)
ba_df = pd.read_csv(ba_file, delimiter=',', header=0, index_col=0)
ba_df.columns = ba_df.columns.map(int)  # make columns names great again
bw_df = pd.read_csv(bw_file, delimiter=',', header=0, index_col=0)
bw_df.columns = bw_df.columns.map(int)
bs_df = pd.read_csv(bs_file, delimiter=',', header=0, index_col=0)
bs_df.columns = bs_df.columns.map(int)

# number crunching: create list of glacier ids for glacier region categories

glac_ALA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ALA), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_WNA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == WNA), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ACN = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ACN), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ACS = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ACS), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_GRL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == GRL), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ISL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ISL), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SJM = input_df.loc[(input_df['GLACIER_REGION_CODE'] == SJM), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SCA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == SCA), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_RUA = input_df.loc[(input_df['GLACIER_REGION_CODE'] == RUA), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASN = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ASN), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_CEU = input_df.loc[(input_df['GLACIER_REGION_CODE'] == CEU), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_CAU = input_df.loc[(input_df['GLACIER_REGION_CODE'] == CAU), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASC = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ASC), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASE = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ASE), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ASW = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ASW), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_TRP = input_df.loc[(input_df['GLACIER_REGION_CODE'] == TRP), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_SAN = input_df.loc[(input_df['GLACIER_REGION_CODE'] == SAN), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_NZL = input_df.loc[(input_df['GLACIER_REGION_CODE'] == NZL), ['GLACIER_REGION_CODE', 'WGMS_ID']]
glac_ANT = input_df.loc[(input_df['GLACIER_REGION_CODE'] == ANT), ['GLACIER_REGION_CODE', 'WGMS_ID']]  # incl. Ant. Peninsula
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
# winter mass-balance data
bw_ALA_df = bw_df.loc[:, list(glac_ALA['WGMS_ID'].unique().tolist())]
bw_WNA_df = bw_df.loc[:, list(glac_WNA['WGMS_ID'].unique().tolist())]
bw_ACN_df = bw_df.loc[:, list(glac_ACN['WGMS_ID'].unique().tolist())]
bw_ACS_df = bw_df.loc[:, list(glac_ACS['WGMS_ID'].unique().tolist())]
bw_GRL_df = bw_df.loc[:, list(glac_GRL['WGMS_ID'].unique().tolist())]
bw_ISL_df = bw_df.loc[:, list(glac_ISL['WGMS_ID'].unique().tolist())]
bw_SJM_df = bw_df.loc[:, list(glac_SJM['WGMS_ID'].unique().tolist())]
bw_SCA_df = bw_df.loc[:, list(glac_SCA['WGMS_ID'].unique().tolist())]
bw_RUA_df = bw_df.loc[:, list(glac_RUA['WGMS_ID'].unique().tolist())]
bw_ASN_df = bw_df.loc[:, list(glac_ASN['WGMS_ID'].unique().tolist())]
bw_CEU_df = bw_df.loc[:, list(glac_CEU['WGMS_ID'].unique().tolist())]
bw_CAU_df = bw_df.loc[:, list(glac_CAU['WGMS_ID'].unique().tolist())]
bw_ASC_df = bw_df.loc[:, list(glac_ASC['WGMS_ID'].unique().tolist())]
bw_ASW_df = bw_df.loc[:, list(glac_ASW['WGMS_ID'].unique().tolist())]
bw_ASE_df = bw_df.loc[:, list(glac_ASE['WGMS_ID'].unique().tolist())]
bw_TRP_df = bw_df.loc[:, list(glac_TRP['WGMS_ID'].unique().tolist())]
bw_SAN_df = bw_df.loc[:, list(glac_SAN['WGMS_ID'].unique().tolist())]
bw_NZL_df = bw_df.loc[:, list(glac_NZL['WGMS_ID'].unique().tolist())]
bw_ANT_df = bw_df.loc[:, list(glac_ANT['WGMS_ID'].unique().tolist())]
# summer mass-balance data
bs_ALA_df = bs_df.loc[:, list(glac_ALA['WGMS_ID'].unique().tolist())]
bs_WNA_df = bs_df.loc[:, list(glac_WNA['WGMS_ID'].unique().tolist())]
bs_ACN_df = bs_df.loc[:, list(glac_ACN['WGMS_ID'].unique().tolist())]
bs_ACS_df = bs_df.loc[:, list(glac_ACS['WGMS_ID'].unique().tolist())]
bs_GRL_df = bs_df.loc[:, list(glac_GRL['WGMS_ID'].unique().tolist())]
bs_ISL_df = bs_df.loc[:, list(glac_ISL['WGMS_ID'].unique().tolist())]
bs_SJM_df = bs_df.loc[:, list(glac_SJM['WGMS_ID'].unique().tolist())]
bs_SCA_df = bs_df.loc[:, list(glac_SCA['WGMS_ID'].unique().tolist())]
bs_RUA_df = bs_df.loc[:, list(glac_RUA['WGMS_ID'].unique().tolist())]
bs_ASN_df = bs_df.loc[:, list(glac_ASN['WGMS_ID'].unique().tolist())]
bs_CEU_df = bs_df.loc[:, list(glac_CEU['WGMS_ID'].unique().tolist())]
bs_CAU_df = bs_df.loc[:, list(glac_CAU['WGMS_ID'].unique().tolist())]
bs_ASC_df = bs_df.loc[:, list(glac_ASC['WGMS_ID'].unique().tolist())]
bs_ASW_df = bs_df.loc[:, list(glac_ASW['WGMS_ID'].unique().tolist())]
bs_ASE_df = bs_df.loc[:, list(glac_ASE['WGMS_ID'].unique().tolist())]
bs_TRP_df = bs_df.loc[:, list(glac_TRP['WGMS_ID'].unique().tolist())]
bs_SAN_df = bs_df.loc[:, list(glac_SAN['WGMS_ID'].unique().tolist())]
bs_NZL_df = bs_df.loc[:, list(glac_NZL['WGMS_ID'].unique().tolist())]
bs_ANT_df = bs_df.loc[:, list(glac_ANT['WGMS_ID'].unique().tolist())]

# number crunching: calculate global mean of regional averages
ba_reg_means_concat_df = pd.concat((ba_ALA_df.mean(axis=1), ba_WNA_df.mean(axis=1), ba_ACN_df.mean(axis=1), ba_ACS_df.mean(axis=1), ba_GRL_df.mean(axis=1),
                                    ba_ISL_df.mean(axis=1), ba_SJM_df.mean(axis=1), ba_SCA_df.mean(axis=1), ba_RUA_df.mean(axis=1), ba_ASN_df.mean(axis=1),
                                    ba_CEU_df.mean(axis=1), ba_CAU_df.mean(axis=1), ba_ASC_df.mean(axis=1), ba_ASW_df.mean(axis=1), ba_ASE_df.mean(axis=1),
                                    ba_TRP_df.mean(axis=1), ba_SAN_df.mean(axis=1), ba_NZL_df.mean(axis=1), ba_ANT_df.mean(axis=1)))
ba_reg_means_by_yr_index_df = ba_reg_means_concat_df.groupby(ba_reg_means_concat_df.index)
ba_global_moravgs_df = round(ba_reg_means_by_yr_index_df.mean(), 0)

bw_reg_means_concat_df = pd.concat((bw_ALA_df.mean(axis=1), bw_WNA_df.mean(axis=1), bw_ACN_df.mean(axis=1), bw_ACS_df.mean(axis=1), bw_GRL_df.mean(axis=1),
                                    bw_ISL_df.mean(axis=1), bw_SJM_df.mean(axis=1), bw_SCA_df.mean(axis=1), bw_RUA_df.mean(axis=1), bw_ASN_df.mean(axis=1),
                                    bw_CEU_df.mean(axis=1), bw_CAU_df.mean(axis=1), bw_ASC_df.mean(axis=1), bw_ASW_df.mean(axis=1), bw_ASE_df.mean(axis=1),
                                    bw_TRP_df.mean(axis=1), bw_SAN_df.mean(axis=1), bw_NZL_df.mean(axis=1), bw_ANT_df.mean(axis=1)))
bw_reg_means_by_yr_index_df = bw_reg_means_concat_df.groupby(bw_reg_means_concat_df.index)
bw_global_moravgs_df = round(bw_reg_means_by_yr_index_df.mean(), 0)

bs_reg_means_concat_df = pd.concat((bs_ALA_df.mean(axis=1), bs_WNA_df.mean(axis=1), bs_ACN_df.mean(axis=1), bs_ACS_df.mean(axis=1), bs_GRL_df.mean(axis=1),
                                    bs_ISL_df.mean(axis=1), bs_SJM_df.mean(axis=1), bs_SCA_df.mean(axis=1), bs_RUA_df.mean(axis=1), bs_ASN_df.mean(axis=1),
                                    bs_CEU_df.mean(axis=1), bs_CAU_df.mean(axis=1), bs_ASC_df.mean(axis=1), bs_ASW_df.mean(axis=1), bs_ASE_df.mean(axis=1),
                                    bs_TRP_df.mean(axis=1), bs_SAN_df.mean(axis=1), bs_NZL_df.mean(axis=1), bs_ANT_df.mean(axis=1))) # without tropics
bs_reg_means_by_yr_index_df = bs_reg_means_concat_df.groupby(bs_reg_means_concat_df.index)
bs_global_moravgs_df = round(bs_reg_means_by_yr_index_df.mean(), 0)

# number crunching: calculate anomalies over reference period
year_ini = 2010
year_fin = 2018

#read years of interest
dict = [{'year': year_ini, 'vei': 5}, {'year': year_fin, 'vei': 5}]
yr_df = pd.DataFrame(yr_dict)
yr_df.set_index('year', inplace=True)

reference_period = range(year_ini, year_fin + 1)
glac_anom = calc_anomalies(ba_ALA_df, reference_period)
reg_anom = glac_anom.mean(axis=1)

# figures: plot mass-balance time series & volcanic eruptions, no data: [9999]
ba_df_list = [ba_ALA_df.mean(axis=1), ba_WNA_df.mean(axis=1), ba_ACN_df.mean(axis=1), ba_ACS_df.mean(axis=1), ba_GRL_df.mean(axis=1),
                    ba_ISL_df.mean(axis=1), ba_SJM_df.mean(axis=1), ba_SCA_df.mean(axis=1), ba_RUA_df.mean(axis=1), ba_ASN_df.mean(axis=1),
                    ba_CEU_df.mean(axis=1), ba_CAU_df.mean(axis=1), ba_ASC_df.mean(axis=1), ba_ASW_df.mean(axis=1), ba_ASE_df.mean(axis=1),
                    ba_TRP_df.mean(axis=1), ba_SAN_df.mean(axis=1), ba_NZL_df.mean(axis=1), ba_ANT_df.mean(axis=1), ba_global_moravgs_df]
bw_df_list = [bw_ALA_df.mean(axis=1), bw_WNA_df.mean(axis=1), bw_ACN_df.mean(axis=1), bw_ACS_df.mean(axis=1), bw_GRL_df.mean(axis=1),
                    bw_ISL_df.mean(axis=1), bw_SJM_df.mean(axis=1), bw_SCA_df.mean(axis=1), bw_RUA_df.mean(axis=1), bw_ASN_df.mean(axis=1),
                    bw_CEU_df.mean(axis=1), bw_CAU_df.mean(axis=1), bw_ASC_df.mean(axis=1), bw_ASW_df.mean(axis=1), bw_ASE_df.mean(axis=1),
                    bw_TRP_df.mean(axis=1), bw_SAN_df.mean(axis=1), bw_NZL_df.mean(axis=1), bw_ANT_df.mean(axis=1), bw_global_moravgs_df]
bs_df_list = [bs_ALA_df.mean(axis=1), bs_WNA_df.mean(axis=1), bs_ACN_df.mean(axis=1), bs_ACS_df.mean(axis=1), bs_GRL_df.mean(axis=1),
                    bs_ISL_df.mean(axis=1), bs_SJM_df.mean(axis=1), bs_SCA_df.mean(axis=1), bs_RUA_df.mean(axis=1), bs_ASN_df.mean(axis=1),
                    bs_CEU_df.mean(axis=1), bs_CAU_df.mean(axis=1), bs_ASC_df.mean(axis=1), bs_ASW_df.mean(axis=1), bs_ASE_df.mean(axis=1),
                    bs_TRP_df.mean(axis=1), bs_SAN_df.mean(axis=1), bs_NZL_df.mean(axis=1), bs_ANT_df.mean(axis=1), bs_global_moravgs_df]

plot_time_series(ba_df_list, bw_df_list, bs_df_list, yr_df)

# figures: plot epoch analysis
plot_epoch_analysis(glac_anom, reg_anom, year_ini, 'Pinatubo & Cerro Hudson')

# tables: export regional mean balances to csv
reg_short_lst = ['Arctic', 'NorthMidLat', 'LowLat', 'SouthMidLat', 'Antarctic', 'Global']
reg_short_excltrp_lst = ['Arctic', 'NorthMidLat', 'SouthMidLat', 'Antarctic', 'Global']
ba_reg_file = path + '\\out_data\\' + 'ba_regions' + '.csv'
bw_reg_file = path + '\\out_data\\' + 'bw_regions' + '.csv'
bs_reg_file = path + '\\out_data\\' + 'bs_regions' + '.csv'

ba_out = pd.concat([ba_arctic_df.mean(axis=1), ba_nmidlats_df.mean(axis=1), ba_lowlats_df.mean(axis=1),
                  ba_smidlats_df.mean(axis=1), ba_antarctic_df.mean(axis=1), ba_global_moravgs_df], axis=1)
ba_out = ba_out.round(decimals=0)
ba_out.columns = reg_short_lst
ba_out.to_csv(ba_reg_file, sep=',', encoding='utf-8', index=True, index_label='Year')

bw_out = pd.concat([bw_arctic_df.mean(axis=1), bw_nmidlats_df.mean(axis=1),
                  bw_smidlats_df.mean(axis=1), bw_antarctic_df.mean(axis=1), bw_global_moravgs_df], axis=1)
bw_out = bw_out.round(decimals=0)
bw_out.columns = reg_short_excltrp_lst
bw_out.to_csv(bw_reg_file, sep=',', encoding='utf-8', index=True, index_label='Year')

bs_out = pd.concat([bs_arctic_df.mean(axis=1), bs_nmidlats_df.mean(axis=1),
                  bs_smidlats_df.mean(axis=1), bs_antarctic_df.mean(axis=1), bs_global_moravgs_df], axis=1)
bs_out = bs_out.round(decimals=0)
bs_out.columns = reg_short_excltrp_lst
bs_out.to_csv(bs_reg_file, sep=',', encoding='utf-8', index=True, index_label='Year')

print('.........................................................................................')
print('"If they are well cleaned out, volcanoes burn slowly and steadily, without any eruptions."')
print('Antoine de Saint-Exupéry')
print('.........................................................................................')
