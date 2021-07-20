"""
regional_ggmc.py

Author: M. Zemp
Date: 16 May 2020
Last changes: 19 May 2020

Scripted for Python 3.7

Description:
This script reads glacier-wide mass-balance data from the WGMS FoG database
and provides functions for related analysis and plots.

Input: C3S Glacier mass balance series
    ggmc_bw-bs-ba.csv (UTF-8 encoding)

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
reg_lst.remove('SAN')
reg_lst= reg_lst + ['SA1','SA2'] # Separate Andes in two regions:

# number crunching: create & export data frame with mass-balance data from input file
ba_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '.csv'
ba_unc_file = path + '\\in_data\\' + 'ggmc_' + 'ba' + '_unc.csv'

# create mass-balance data csv if it has not been produced before
# ba_df = create_mb_dataframe(input_df, wgms_id_lst, yr_lst, 'ANNUAL_BALANCE')
# ba_df.to_csv(ba_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')
#
# ba_unc_df = create_mb_dataframe(input_df, wgms_id_lst, yr_lst, 'ANNUAL_BALANCE_UNC')
# ba_unc_df.to_csv(ba_unc_file, sep=',', encoding='utf-8', index=True, index_label='YEAR')

# read mass-balance data and uncertainty from csv into dataframe
ba_df = pd.read_csv(ba_file, delimiter=',', header=0, index_col=0)
ba_df.columns = ba_df.columns.map(int)  # make columns names great again

ba_unc_df = pd.read_csv(ba_unc_file, delimiter=',', header=0, index_col=0)
ba_unc_df.columns = ba_unc_df.columns.map(int)  # make columns names great again

# # create empty dataframes for regional anomalies and uncertainties
reg_anom_df = pd.DataFrame(index=yr_lst)
reg_anom_df.index.name = 'YEAR'
reg_anom_lst=[]

reg_anom_err_df = pd.DataFrame(index=yr_lst)
reg_anom_err_df.index.name = 'YEAR'
sig_reg_anom_lst=[]

# Define reference period
year_ini = 2009
year_fin = 2018

reference_period = range(year_ini, year_fin + 1)
# print(list(reference_period))

plt_year_min = 1950 # starting year for Regional series plot
plt_year_max = 2020 # end year for Regional series plot

axes = 'eq_axes' # all region with same Y axes, visualizes best the contribution between regions
# axes = 'tight' # fits Y axes to each region, visualizes best the temporal variability of each region


for region in reg_lst:
    region='CEU'
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

    elif region == 'NZL':
        add_id_lst = [2000]  # Martial Este (SAN-01)
        glac = input_df.loc[((input_df['GLACIER_REGION_CODE'] == 'NZL') | (input_df['WGMS_ID'].isin(add_id_lst))), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    else:
        glac = input_df.loc[(input_df['GLACIER_REGION_CODE'] == str(region)), ['GLACIER_REGION_CODE', 'WGMS_ID']]

    ## number crunching:   select mass-balance data for glacier region groups
    ba_reg_df = ba_df.loc[:, list(glac['WGMS_ID'].unique().tolist())]
    ba_unc_reg_df = ba_unc_df.loc[:, list(glac['WGMS_ID'].unique().tolist())]
    # print(ba_unc_reg_df)

    gla_dir = path + '\\out_data\\'+str(region)+'_glaciological_sample_figures\\'
    if not os.path.exists(gla_dir):
        os.mkdir(gla_dir)

    # # PLOT:  All selected glacier mass balance series by region
    # ba_reg_df = ba_reg_df/ 1000
    # wgms_id = 491
    #
    #
    # ax= ba_reg_df.plot(color='Silver', alpha=0.7, linewidth= 0.8,fontsize='large', legend=False)
    # ax.axhline(0, color='Silver', linewidth=1)
    # ax = ba_reg_df[wgms_id].plot(color='rebeccapurple', alpha=0.3, linewidth= 1.5,fontsize='large', legend=False)
    #
    # ax.set_xlabel('Year', fontsize='large')
    # ax.set_xlim([1950, 2020])
    # ax.set_ylim([-4.5, 3.0])
    # ax.set_ylabel(r'B$_{glac}$ (m w.e. yr$^{-1}$)', fontsize='large')
    # # ax.legend(ncol=3, title= 'RGIId', loc ='best', fontsize=6)
    #
    # ax.text(2003, 2.5, 'N = ' + str(len(ba_reg_df.columns)) + ' glaciers')
    #
    # out_fig= gla_dir +str(region)+'_all_glaciers_Ba_series.png'
    # plt.savefig(out_fig)
    # plt.close()
    #
    # print('Plot saved as {}.'.format(out_fig))
    # # plt.show()
    # exit()
    #
    # # PLOT:  ONE glacier mass balance series
    # wgms_id = 491
    #
    # ax = ba_reg_df[wgms_id].plot(color='rebeccapurple', alpha=0.3, linewidth= 1.5,fontsize='large', legend=False)
    # ax.axhline(0, color='Silver', linewidth=1)
    # ax.set_xlabel('Year', fontsize='large')
    #
    # ax.set_xlim([plt_year_min, plt_year_max])
    # ax.set_ylim([-3.5, 2.0])
    # ax.set_ylabel(r'B$_{glac}$ (m w.e. yr$^{-1}$)', fontsize='large')
    #
    # out_fig = gla_dir + str(region) + '_glacier_'+str(wgms_id)+'_Ba_series'+str(plt_year_min)+'-'+str(plt_year_min)+'.png'
    # plt.savefig(out_fig)
    #
    # # plt.close()
    # print('Plot saved as {}.'.format(out_fig))
    # # plt.show()
    # exit()


    ## CALCULATE:  individual glacier anomalies by region
    glac_anom = calc_anomalies(ba_reg_df, reference_period, region)
    glac_mba_unc = calc_anomalies_unc(ba_reg_df,ba_unc_reg_df, reference_period, region)

    # if no uncertainty measurement use the regional annual mean uncertainty of the glaciological sample (SA1 and SA2)
    if glac_mba_unc.isnull().sum().sum():
        for id in glac_mba_unc.columns.tolist():
            year_min = glac_anom[id].first_valid_index()
            yrs = list(range(1885, year_min))
            glac_mba_unc[id].fillna(np.nanmean(ba_unc_reg_df), inplace=True)
            glac_mba_unc[id].mask(glac_mba_unc.index.isin(yrs), np.nan, inplace=True)
    else:
        continue

    # glac_anom.to_csv(path + '\\out_data\\Reg_anomaly_ref_period_' + str(year_ini) + '-' + str(year_fin)+'\\'+region+'_gla_anomalies.csv')
    # glac_mba_unc.to_csv(path + '\\out_data\\Reg_valid_gla_mba_uncertainties\\' + region + '_gla_mba_uncertainties.csv')

    ## CALCULATE:  regional anomaly mean and uncertainty
    ## Regional mean
    reg_anom = glac_anom.mean(axis=1)
    anom_df=pd.DataFrame(reg_anom, columns= [str(region)])
    reg_anom_lst.append(anom_df)

    ## Uncertainty
    mean_gla_mba_err = glac_mba_unc.mean(axis=1)
    mean_gla_mba_err_df=pd.DataFrame(mean_gla_mba_err, columns= [str(region)])
    # print(mean_gla_mba_err_df)

    nb_gla=(len(glac_anom.columns))
    if nb_gla >1:
        gla_anom_std = glac_anom.std(axis=1)
        var_gla_anom_err = 1.96 * gla_anom_std

    elif nb_gla == 1:
        var_gla_anom_err = mean_gla_mba_err_df*0

    var_gla_anom_err_df = pd.DataFrame(var_gla_anom_err, columns=[str(region)])
    # print(var_gla_anom_err_df)

    #fill nans with mean error
    min_year = anom_df.first_valid_index()
    yrs = list(range(1885, min_year))
    var_gla_anom_err_df.fillna(var_gla_anom_err_df.mean(), inplace=True)
    var_gla_anom_err_df.loc[var_gla_anom_err_df.index.isin(yrs), [str(region)]] = np.nan

    sig_anom_df = np.sqrt(var_gla_anom_err_df**2 + mean_gla_mba_err_df**2)
    sig_reg_anom_lst.append(sig_anom_df)

    ## observe correlation matrix between glacier anomalies
    # corr_df=glac_anom.corr(method='pearson')
    # print(corr_df.round(decimals=2))

    # ## PLOT glacier anomalies and Regional mean

    if len(glac_anom.columns) > 1:
        plot_anom_df = anom_df.mul(1)/1000
        plot_anom_df['UNC'] = var_gla_anom_err_df /1000
        glac_anom = glac_anom/1000

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax = glac_anom.plot(color= 'Silver',linewidth= 0.5, alpha=0.8, legend=False)

        # plt.title(region, fontsize=20)
        ax.set_xlim([plt_year_min, plt_year_max])
        ax.set_ylim([-3.5, 4.0])
        ax.axhline(0, color='Grey', linewidth=1)
        ax.set_ylabel('\u03B2 (m w.e.)', fontsize=11)
        ax.tick_params(labelsize=11)
        # ax.text(1984, 3.5, 'N = '+str(len(glac_anom.columns))+' glaciological observations')
        ax.text(1952, -3.1, 'b - \u03B2$_R$ (B$_{gla}$)', fontsize=22)
        # ax.text(1952, -3.2, 'reference period = 2009-2018')

        wgms_id=491
        glac_anom[wgms_id].plot(ax=ax, color='rebeccapurple', alpha=0.3, linewidth= 1.5)
        plot_anom_df[region].plot(ax=ax, color='black', alpha=0.7, linewidth=1.5)
        plt.fill_between(plot_anom_df.index, plot_anom_df[region] + plot_anom_df['UNC'], plot_anom_df[region] - plot_anom_df['UNC'], color='grey', alpha=0.3, linewidth= 0)
        out_fig= gla_dir+ 'Anomaly_reg_'+region+'_ref_period_'+str(year_ini)+'-'+str(year_fin)+'_'+str(plt_year_min)+'-'+str(plt_year_max)+'.pdf'
        plt.savefig(out_fig)
        print('Plot saved as {}.'.format(out_fig))
        # reg_anom.plot()
        plt.close()
    else:
        print('................... Region with only one glacier anomaly:', region, '............................')
    exit()



print('.........................................................................................')
print('"The End"')
print('.........................................................................................')
