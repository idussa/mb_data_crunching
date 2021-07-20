"""
Calculate the regional mass loss

calc_reg_mass_loss.py

Author: idussa
Date: Feb 2021
Last changes: Feb 2021

Scripted for Python 3.7

Description:
This script reads glacier-wide mass balance data edited from WGMS FoG database
and regional glacier anomalies produced by calc_regional_anomalies_and_error.py
and provides the observational consensus estimate for every individual glacier
with available geodetic observations WGMS Id

Input:  C3S_GLACIER_DATA_20200824.csv
        OCE_files_by_region\\
        (UTF-8 encoding)

Return: tbd.svg
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions_ggmc import *
import scipy  # usage: scipy.stats.funct()

rgi_region= {'ACN' : 'ArcticCanadaNorth',
             'WNA' : 'WesternCanadaUS',
             'ALA' : 'Alaska',
             'ACS' : 'ArcticCanadaSouth',
             'TRP' : 'LowLatitudes',
             'SCA' : 'Scandinavia',
             'SJM' : 'Svalbard',
             'CEU' : 'CentralEurope',
             'CAU' : 'CaucasusMiddleEast',
             'ASC' : 'CentralAsia',
             'ASN' : 'NorthAsia',
             'ASE' : 'SouthAsiaEast',
             'NZL' : 'NewZealand',
             'ASW' : 'SouthAsiaWest',
             'GRL' : 'GreenlandPeriphery',
             'ANT' : 'AntarcticSubantarctic',
             'ISL' : 'Iceland',
             'RUA' : 'RussianArctic',
             'SA1' : 'SouthernAndes',
             'SA2' : 'SouthernAndes'}

rgi_code= {
           'ALA' : '01',
           'WNA' : '02',
           'ACN' : '03',
           'ACS' : '04',
           'GRL' : '05',
           'ISL' : '06',
           'SJM' : '07',
           'SCA' : '08',
           'RUA' : '09',
           'ASN' : '10',
           'CEU' : '11',
           'CAU' : '12',
           'ASC' : '13',
           'ASW' : '14',
           'ASE' : '15',
           'TRP' : '16',
           'SA1' : '17',
           'SA2' : '17',
           'NZL' : '18',
           'ANT' : '19'}

reg_lst= ['ALA', 'WNA', 'ACN', 'ACS', 'GRL', 'ISL', 'SJM', 'SCA', 'RUA', 'ASN', 'CEU', 'CAU', 'ASC', 'ASW', 'ASE', 'TRP', 'SA1', 'SA2', 'NZL', 'ANT']


##########################################
##########################################
"""main code"""
##########################################
##########################################

## Define input paths
path = os.path.dirname(os.path.abspath(__file__))


in_ba_file = path + '\\in_data\\Regional_B_series_AreaWeighted.csv'
in_sig_ba_file = path + '\\in_data\\Regional_B_series_uncertainty.csv'
in_area_zemp_file = path+ '\\in_data\\Regional_area_change_Zemp.csv'
in_rgi = 'C:\\Users\\idussail2\\Documents\\PROJECTS\\G3P_project\\data\\00_rgi60\\00_rgi60_attribs\\' # path to rgi attribute file

ba_df = pd.read_csv(in_ba_file, encoding='latin1', delimiter=',', header=0, index_col='YEAR')
sig_ba_df = pd.read_csv(in_sig_ba_file, encoding='latin1', delimiter=',', header=0, index_col='YEAR')
area_zemp_df = pd.read_csv(in_area_zemp_file, encoding='latin1', delimiter=',', header=0, index_col='YEAR')

# Define Variables
PoR = list(range(2006, 2017)) # period to calculate the cumulative mass loss

S_ocean = 362.5 * 10**6 # Cogley et al 2012
sig_S_ocean = 362.5 * 10**6 # Cogley et al 2012

sig_area = 0.05 # Paul et al 2015

# Define Output parameters
out_dir = path + '\\out_data\\'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

plt_year_min = 1950 # starting year for Regional series plot
plt_year_max = 2020 # end year for Regional series plot

axes = 'eq_axes' # all region with same Y axes, visualizes best the contribution between regions
# axes = 'tight' # fits Y axes to each region, visualizes best the temporal variability of each region

############################################################################################################################

cols = ['region', 'area_mean_km2', 'nb_gla_rgi' ,'DM_Gt_yr', 'sigma_DM_Gt_yr', 'SLE_mm_yr', 'sigma_SLE_mm_yr']
glob_cum_df = pd.DataFrame(index=reg_lst, columns=cols)

for region in reg_lst:
    # region='CEU'
    out_DM_series = out_dir +'regional_mass_loss_series\\'
    if not os.path.exists(out_DM_series):
        os.mkdir(out_DM_series)

    in_rgi_area = in_rgi + rgi_code[region] +'_rgi60_'+ rgi_region[region] +'.csv'
    rgi_area_df = pd.read_csv(in_rgi_area, encoding='latin1', delimiter=',', header=0, usecols= ['RGIId', 'O1Region', 'O2Region', 'Area'], index_col='RGIId').sort_index()

    if region == 'SA1':
        rgi_area_df= rgi_area_df.loc[rgi_area_df['O2Region']== 1]

    if region == 'SA2':
        rgi_area_df= rgi_area_df.loc[rgi_area_df['O2Region']== 2]

    ba_mmwe = ba_df[region]
    ba_kmwe = ba_mmwe/10**6

    sig_ba_mmwe = sig_ba_df[region]
    sig_ba_kmwe = sig_ba_mmwe/10**6

    area = area_zemp_df[region]

    dm_Gt = ba_kmwe * area

    sig_dm_sum = np.sqrt((sig_ba_kmwe/ba_kmwe)**2 + sig_area**2)
    sig_dm = np.abs(dm_Gt) * sig_dm_sum

    reg_file =pd.DataFrame()
    reg_file['Aw_mwe'] = ba_df[region]
    reg_file['sig_tot_mwe'] = sig_ba_df[region]
    reg_file['area_tot_km2'] = area_zemp_df[region]
    reg_file['DM_Gt'] = dm_Gt
    reg_file['sig_tot_DM'] = sig_dm

    reg_file.to_csv(out_DM_series + 'results_region_' + rgi_code[region] + '_' + region + '_' + rgi_region[region] + '.csv')

    # # Calculate cumulative mass loss for PoR

    mean_area = reg_file['area_tot_km2'].loc[(reg_file.index.isin(PoR))].mean()
    DM_Gt_yr = reg_file['DM_Gt'].loc[(reg_file.index.isin(PoR))].sum() / len(PoR)
    sigma_DM_Gt_yr = np.sqrt(reg_file['sig_tot_DM'].loc[(reg_file.index.isin(PoR))].pow(2).sum()) / len(PoR)
    SLE = (-DM_Gt_yr / S_ocean) * 10**6
    Sigma_SLE = (sigma_DM_Gt_yr / S_ocean) * 10**6

    for index, row in glob_cum_df.iterrows():
        if index == region:
            row['region']= rgi_region[region]
            row['area_mean_km2'] = "{:.0f}".format(mean_area)
            row['nb_gla_rgi'] = len(rgi_area_df.index)
            row['DM_Gt_yr'] = "{:.2f}".format(DM_Gt_yr)
            row['sigma_DM_Gt_yr'] = "{:.2f}".format(sigma_DM_Gt_yr)
            row['SLE_mm_yr'] = "{:.3f}".format(SLE)
            row['sigma_SLE_mm_yr'] = "{:.3f}".format(Sigma_SLE)

    # ## PLOT regional mass loss
    reg_file['DM_Gt'].plot(color='black', alpha=0.7, linewidth=1)
    plt.title(region +'glacier mass loss (Gt)', fontsize=18)
    plt.xlim([plt_year_min, plt_year_max])
    if axes == 'eq_axes':
        plt.ylim([-200, 100])
    plt.axhline(0, color='Grey', linewidth=1)
    plt.ylabel('Mass loss (Gt)', fontsize='medium')
    plt.xlabel('Year', fontsize='medium')

    plt.fill_between(reg_file['DM_Gt'].index, reg_file['DM_Gt'] + reg_file['sig_tot_DM'],
                     reg_file['DM_Gt'] - reg_file['sig_tot_DM'], color='purple', alpha=0.2, linewidth=0)

    out_plot_dir = out_dir + 'plot_reg_mass_loss_Gt_'+str(plt_year_min)+'-'+str(plt_year_max)+'_'+axes+'\\'
    if not os.path.exists(out_plot_dir):
        os.mkdir(out_plot_dir)

    out_fig = out_plot_dir + region + '_dM_series_Gt.pdf'
    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))
    plt.close()

glob_cum_df.to_csv(out_dir + 'Cum_DM_PoR_'+str(min(PoR))+'_'+str(max(PoR))+'.csv')

print('.........................................................................................')
print('"The End"')
print('.........................................................................................')
exit()


