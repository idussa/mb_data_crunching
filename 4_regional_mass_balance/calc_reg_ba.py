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
path_rgi = 'C:\\Users\\idussail2\\Documents\\PROJECTS\\G3P_project\\data\\00_rgi60\\00_rgi60_attribs\\' # path to rgi attribute file
path_oce = path + '\\in_data\\OCE_files_by_region\\'# path to regional OCE files

##    Define output parameters

out_dir = path + '\\out_data\\'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

DM_series_min_yr = 1950 # starting year of DM series saved in the .csv files

# end_plot= '_OCE_and_B_reg_' # plots Regional Mass balance and the individual region OCE
end_plot = '_B_reg_and_err_' # plots Regional Mass balance and the relative mass balance uncertainty

# gla_area = 'area_CS3' # Uses CS3 survey year individual glacier areas
gla_area = 'area_FoG-Latest' # Uses Latest FoG individual glacier areas

plt_year_min = 1950 # starting year for Regional series plot
plt_year_max = 2020 # end year for Regional series plot

axes = 'eq_axes' # all region with same Y axes, visualizes best the contribution between regions
# axes = 'tight' # fits Y axes to each region, visualizes best the temporal variability of each region

##    Define glacier areas for area weighting
if gla_area == 'area_FoG-Latest':
    in_fog_area = path + '\\in_data\\FoG_20200824_latest-S.csv'
    fog_area_df= pd.read_csv(in_fog_area, encoding='latin1', delimiter=',', header=0, usecols=['WGMS_ID','S'], index_col='WGMS_ID').sort_index()
    area_gla_df = fog_area_df.dropna()

elif gla_area == 'area_CS3':
    in_fog_area = path + '\\in_data\\GEO_MASS_BALANCE_DATA_20200824.csv'
    fog_area_df = pd.read_csv(in_fog_area, encoding='latin1', delimiter=',', header=0, usecols=['WGMS_ID','AREA_SURVEY_YEAR'],index_col='WGMS_ID').sort_index()
    fog_area_df = fog_area_df.rename(columns = {'AREA_SURVEY_YEAR': 'S'})
    area_gla_df = fog_area_df.groupby(fog_area_df.index).mean()
    area_gla_df= area_gla_df.dropna()

print(len(area_gla_df))
exit()
# Define correlation factor for error propagation
corr_factor = 50 # from Zemp et al. 2019

############################################################################################################################

###### Calculate Total glacier mass loss by region ######

Reg_mb_df = pd.DataFrame()
Reg_sig_mb_df = pd.DataFrame()

for region in reg_lst:
    # region='TRP'
    print('working on region, ', region)

    ## Define and read input:   regional OCE series and Uncertainty
    in_oce_file = path_oce + region + '_regional_OCEs.csv'
    in_sig_oce_file = path_oce + region + '_regional_sigma_OCEs.csv'

    oce_df = pd.read_csv(in_oce_file, encoding='latin1', delimiter=',', header=0, index_col='YEAR')
    sig_oce_df = pd.read_csv(in_sig_oce_file, encoding='latin1', delimiter=',', header=0, index_col='YEAR')

    in_rgi_area = path_rgi + rgi_code[region] +'_rgi60_'+ rgi_region[region] +'.csv'
    rgi_area_df = pd.read_csv(in_rgi_area, encoding='latin1', delimiter=',', header=0, usecols= ['RGIId', 'O1Region', 'O2Region', 'Area'], index_col='RGIId').sort_index()

    if region == 'SA1':
        rgi_area_df= rgi_area_df.loc[rgi_area_df['O2Region']== 1]

    if region == 'SA2':
        rgi_area_df= rgi_area_df.loc[rgi_area_df['O2Region']== 2]

    nb_gla_reg = len(rgi_area_df.index)
    tot_area_rgi_reg = rgi_area_df['Area'].sum()

    ## select areas of wgms_ids belonging to the region group
    wgms_id_lst = oce_df.columns.to_list()
    wgms_id_lst = [int(i) for i in wgms_id_lst]


    ## Calculate total area of observed glaciers presenting an area value in FoG

    id_lst=[]
    for id in wgms_id_lst:
        if id in area_gla_df.index:
            id_lst.append(id)
        else:
            pass

    area_gla_obs= area_gla_df.loc[id_lst]
    nb_gla_obs= len(area_gla_obs.index)
    tot_area_obs = area_gla_obs['S'].sum()

    print('area glaciers observed', area_gla_obs)
    print('number glaciers with OCE and area', nb_gla_obs)
    print('total area observed', tot_area_obs)

    ####### Calculate unobserved glaciers time series and uncertainties ##########

    ## 1. Calculate OCE series for unobserved glaciers as the Weigthed mean from the regional glacier sample with observations
    rel_mb_df = pd.DataFrame()
    rel_sig_mb_df = pd.DataFrame()

    for id in id_lst:
        print('working on glacier, ', id)

        area= area_gla_obs.loc[id, 'S']
        mb_oce= oce_df[str(id)]
        mb_oce_rel = (mb_oce * area) / tot_area_obs
        rel_mb_df[id] = mb_oce_rel

        sig_oce = sig_oce_df[str(id)]
        sig_rel = (sig_oce * area)/tot_area_obs
        rel_sig_mb_df[id] = sig_rel

    # Result: OCE for unobserved glaciers
    Aw_oce_obs_df = rel_mb_df.sum(axis=1, min_count=1)

    ## 2. Calculate OCE uncertainties for unobserved glaciers

    # Weighted mean Sigma OCE
    Aw_sig_obs_df = rel_sig_mb_df.sum(axis=1, min_count=1)
    # print(Aw_sig_obs_df)

    # Regional OCE Arithmetic mean
    Arith_oce = oce_df.mean(axis=1)

    # Error Variability of intepolation method (Area Weighted mean oce’s vs Arithmetic mean oce’s)
    mb_int_df = pd.concat([Aw_oce_obs_df, Arith_oce], axis=1)
    mb_int_std = mb_int_df.std(axis=1)
    mb_int_var_err = 1.96 * mb_int_std
    # print(mb_int_var_err)

    # Result: sigma OCE for unobserved glaciers
    Sig_oce_unobs_gla = np.sqrt(mb_int_var_err ** 2 + Aw_sig_obs_df ** 2)
    # print('-------------------------sig unobs: ',Sig_oce_unobs_gla)

    ## 3. Add OCE series and uncertainties for unobserved glaciers
    # Id -9999 for unobserved glaciers, OCE is the area weighthed average of the regional observed series
    oce_df[-9999]= Aw_oce_obs_df
    out_oce = out_dir + 'OCE_all_gla_by_region\\'
    if not os.path.exists(out_oce):
        os.mkdir(out_oce)
    oce_df.to_csv(out_oce + region +'_all_gla_OCEs.csv')

    sig_oce_df[-9999] = Sig_oce_unobs_gla
    sig_oce_df.to_csv(out_dir + 'OCE_all_gla_by_region\\' + region + '_all_gla_sigma_OCEs.csv')

    ####### Calculate Regional specific mass balance time series ##########

    Reg_mb_df[region] = Aw_oce_obs_df
    nb_unobs_gla = nb_gla_reg - nb_gla_obs

    rel_Aw_reg_sig_df = pd.DataFrame()
    new_id_lst= id_lst+[-9999]

    for id in new_id_lst:
        # id=0
        print('working on glacier, ', id)
        if id == -9999:
            area_unobs = tot_area_rgi_reg - tot_area_obs
            sig_rel_unobs = (Sig_oce_unobs_gla * area_unobs) / tot_area_rgi_reg
            rel_Aw_reg_sig_df[-9999]= sig_rel_unobs
        else:
            area = area_gla_obs.loc[id, 'S']
            sig_reg = sig_oce_df[str(id)]
            sig_rel = (sig_reg * area) / tot_area_rgi_reg
            rel_Aw_reg_sig_df[id] = sig_rel

    # rel_Aw_reg_sig_df.to_csv(out_dir + region +'_sig_test.csv')

    if nb_gla_obs < corr_factor:
        reg_sig = np.sqrt(rel_Aw_reg_sig_df.pow(2).sum(axis=1, min_count= 1))

    elif nb_gla_obs >= corr_factor:
        reg_sig_sqsum = rel_Aw_reg_sig_df.pow(2).sum(axis=1, min_count= 1)
        reg_sig_corr = (1 / (nb_gla_obs / corr_factor)) * reg_sig_sqsum
        reg_sig = np.sqrt(reg_sig_corr)

    Reg_sig_mb_df[region] = reg_sig

    # ## PLOT regional mass balance
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot zero balance line
    ax.axhline(0, color='Grey', linewidth=1)

    if end_plot == '_OCE_and_B_reg_':
        # # plot regional individual OCEs
        for item in oce_df.columns:
            ax.plot(oce_df.index, oce_df[item], color='Silver', alpha=0.5, linewidth=0.5)

    elif end_plot == '_B_reg_and_err_':
        # plot regional B uncertainty
        plt.fill_between(Reg_mb_df.index, Reg_mb_df[region] + Reg_sig_mb_df[region],
                         Reg_mb_df[region] - Reg_sig_mb_df[region], color='purple', alpha=0.5, linewidth=0)

    # plot regional B
    ax.plot(Reg_mb_df[region], color='black', label='area weighted mean', alpha=0.7, linewidth= 1)


    per_gla_obs= "{:.1f}".format(nb_gla_obs/nb_gla_reg * 100)
    ax.set_title('region '+ region + ' \nN OCE='+str(nb_gla_obs), fontsize=12)
    ax.axhline(0, color='Grey', linewidth=1)
    ax.set_xlabel('Year')
    ax.set_ylabel(r'B$_{R}$ (m w.e. yr$^{-1}$) (mm w.e.)', fontsize='medium')

    plt.xlim([plt_year_min, plt_year_max])

    if axes == 'eq_axes':
        plt.ylim((-4500, 4500))
    plt.legend(loc=3)
    # plt.tight_layout()

    out_plot_dir = out_dir + 'plot_'+end_plot+'_'+gla_area+'_'+str(plt_year_min)+'-'+str(plt_year_max)+'_'+axes+'\\'
    if not os.path.exists(out_plot_dir):
        os.mkdir(out_plot_dir)

    out_fig= out_plot_dir + region +'_Ba_series_Aw.png'

    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))
    # plt.show()
    plt.close()

Reg_mb_df = Reg_mb_df.loc[(Reg_mb_df.index >= 1950)]
Reg_sig_mb_df = Reg_sig_mb_df.loc[(Reg_sig_mb_df.index >= 1950)]

### Save regional Mass balance series
Reg_mb_df.to_csv(out_dir + 'Regional_B_series_AreaWeighted.csv')
Reg_sig_mb_df.to_csv(out_dir + 'Regional_B_series_uncertainty.csv')


print('.........................................................................................')
print('"The End"')
print('.........................................................................................')
exit()
