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

##    Define input paths
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

# axes = 'eq_axes' # all region with same Y axes, visualizes best the contribution between regions
axes = 'tight' # fits Y axes to each region, visualizes best the temporal variability of each region

############################################################################################################################
reg_dm_file = pd.DataFrame()
sig_reg_dm_file = pd.DataFrame()

for region in reg_lst:
    # region='CEU'

    ba_mmwe = ba_df[region]
    ba_kmwe = ba_mmwe/10**6

    sig_ba_mmwe = sig_ba_df[region]
    sig_ba_kmwe = sig_ba_mmwe/10**6

    area = area_zemp_df[region]

    dm_Gt = ba_kmwe * area

    sig_dm_sum = np.sqrt((sig_ba_kmwe/ba_kmwe)**2 + sig_area**2)
    sig_dm = np.abs(dm_Gt) * sig_dm_sum

    reg_dm_file[region] = dm_Gt
    sig_reg_dm_file[region] = sig_dm




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



print('.........................................................................................')
print('"The End"')
print('.........................................................................................')
exit()


