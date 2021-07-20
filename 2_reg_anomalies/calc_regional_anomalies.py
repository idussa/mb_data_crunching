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

## Separate Andes in two regions:
reg_lst.remove('SAN')
reg_lst= reg_lst + ['SA1','SA2']


# number crunching: create & export data frame with mass-balance data from input file
ba_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '.csv'

# read mass-balance data from csv into dataframe (if produced)
ba_df = pd.read_csv(ba_file, delimiter=',', header=0, index_col=0)
ba_df.columns = ba_df.columns.map(int)  # make columns names great again

# # number crunching: create list of glacier ids for glacier region categories
reg_anom_df = pd.DataFrame(index=yr_lst)
reg_anom_df.index.name = 'YEAR'
reg_anom_lst=[]

# Define reference period
year_ini = 2009
year_fin = 2018

reference_period = range(year_ini, year_fin + 1)
print(list(reference_period))


for region in reg_lst:
    # region='CEU'
    print('working on region... ', str(region))

    ## number crunching:  create list of glacier ids for glacier region categories
    # add neighbouring glaciers for regions not represented in the present period
    if region == 'RUA':
        add_reg_lst = ['RUA', 'SJM']
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'].isin(add_reg_lst))), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    elif region == 'ASN':
        add_id_lst = [853] # Urumqi (ASC)
        rem_id = 897  # Hamagury yuki
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ASN') | (input_df['WGMS_ID'].isin(add_id_lst))), ['GLACIER_REGION_CODE', 'WGMS_ID']]
        glac = glac.drop(glac[glac['WGMS_ID'] == rem_id].index)

    elif region == 'ASC':
        rem_id_lst = [1511, 1512]  # Urumqi East and west branches
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ASC')), ['GLACIER_REGION_CODE', 'WGMS_ID']]
        glac = glac.drop(glac[glac['WGMS_ID'].isin(rem_id_lst)].index)

    elif region == 'SA1':
        rem_id_lst = [3902,3903,3904,3905,1344] # keep Martial Este only
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'SAN')), ['GLACIER_REGION_CODE', 'WGMS_ID']]
        glac = glac.drop(glac[glac['WGMS_ID'].isin(rem_id_lst)].index)

    elif region == 'SA2': # keep Echaurren Norte only
        rem_id_lst = [3902,3903,3904,3905,2000]
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'SAN')), ['GLACIER_REGION_CODE', 'WGMS_ID']]
        glac = glac.drop(glac[glac['WGMS_ID'].isin(rem_id_lst)].index)

    elif region == 'ASW':
        add_id_lst = [817, 853] # Ts. Tuyuksuyskiy (ASC), Urumqi (ASC)
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ASW') | (input_df['WGMS_ID'].isin(add_id_lst))), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    elif region == 'ASE':
        add_id_lst = [817, 853] # Ts. Tuyuksuyskiy (ASC), Urumqi (ASC)
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ASE') | (input_df['WGMS_ID'].isin(add_id_lst))), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    elif region == 'ACS':
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'ACS') | (input_df['GLACIER_REGION_CODE'] == 'ACN')), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    else:
        glac = input_df.loc[(input_df['GLACIER_REGION_CODE'] == str(region)), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    ## number crunching:   select mass-balance data for glacier region groups
    ba_reg_df = ba_df.loc[:, list(glac['WGMS_ID'].unique().tolist())]
    print(ba_reg_df)

    ## plot regional time series
    # ax= ba_reg_df.plot(title='Temporal coverage Glaciological obs. region: ' + str(region) + '\n Nb glaciers measured: ' + str(len(ba_reg_df.columns)), fontsize='medium')
    # ax.set_xlabel('Year', fontsize='medium')
    # ax.set_xlim([1900, 2020])
    # ax.set_ylabel('MB (mm w.e. yr-1)', fontsize='medium')
    # ax.legend(ncol=3, title= 'RGIId', loc ='best', fontsize=6)
    # out_fig= path + '\\out_data\\Reg_glac_series_gla_add\\'+str(region)+'_glaciologial_obs.png'
    # plt.savefig(out_fig)
    # print('Plot saved as {}.'.format(out_fig))
    # # plt.show()

    ## plot regional anomalies and anomaly mean
    glac_anom = calc_anomalies(ba_reg_df, reference_period, region)
    print(glac_anom)
    glac_anom.to_csv(path + '\\out_data\\Reg_anomaly_ref_period_' + str(year_ini) + '-' + str(year_fin)+'\\'+region+'_anomalies.csv')

    ## observe correlation matrix between glacier anomalies
    # corr_df=glac_anom.corr(method='pearson')
    # print(corr_df.round(decimals=2))

    ax = glac_anom.plot(colormap= 'Paired', legend=False)
    plt.title(region, fontsize=20)
    ax.set_xlim([1995, 2020])
    ax.set_ylim([-3500, 4000])
    ax.axhline(0, color='Grey', linewidth=1)
    ax.set_ylabel('Annomaly (mm w.e.)', fontsize='medium')
    ax.text(2014, 3500, 'N = '+str(len(glac_anom.columns))+' glaciers')

    reg_anom = glac_anom.mean(axis=1)
    reg_anom.plot(ax=ax, color='black', alpha=0.7, linewidth= 3)
    out_fig= path + '\\out_data\\plot_reg_anomaly_ref_period_'+str(year_ini)+'-'+str(year_fin)+'_final\\Anomaly_reg_'+region+'_ref_period_'+str(year_ini)+'-'+str(year_fin)+'.png'
    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))
    # reg_anom.plot()
    # plt.show()

    anom_df=pd.DataFrame(reg_anom, columns= [str(region)])
    reg_anom_lst.append(anom_df)

reg_anom_df=pd.concat(reg_anom_lst, axis='columns')

print(reg_anom_df)
# reg_anom_df.plot()
# plt.show()

reg_anom_df.to_csv(path + '\\out_data\\Mean_regional_anomalies_ref_period_'+str(year_ini)+'-'+str(year_fin)+'.csv')

exit()









print('.........................................................................................')
print('"If they are well cleaned out, volcanoes burn slowly and steadily, without any eruptions."')
print('Antoine de Saint-Exupéry')
print('.........................................................................................')
