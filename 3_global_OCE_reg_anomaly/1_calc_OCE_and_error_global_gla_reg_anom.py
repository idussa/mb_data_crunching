"""
Calculate the observational consensus estimate for every individual glacier

calc_OCE_and_error_global_gla_reg_anom.py

Author: idussa
Date: Feb 2021
Last changes: Feb 2021

Scripted for Python 3.7

Description:
This script reads glacier-wide mass balance data edited from WGMS FoG database
and regional glacier anomalies produced by calc_regional_anomalies_and_error.py
and provides the observational consensus estimate for every individual glacier
with available geodetic observations WGMS Id

Input:  GEO_MASS_BALANCE_DATA_20200824.csv
        Regional_anomalies_ref_period_2009-2018.csv
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

##########################################
##########################################
"""main code"""
##########################################
##########################################

# Define input
##### 1. PATH TO FILES ######

path = os.path.dirname(os.path.abspath(__file__))
# path = os.path.dirname(os.path.abspath(__file__))+'\\codes\\mb_data_crunching'

# read file with global geodetic mass-balance data and regional anomalies
in_data_geo = path + '\\in_data\\GEO_MASS_BALANCE_DATA_20200824.csv'
in_data_anom = path + '\\in_data\\Regional_anomalies_ref_period_2009-2018.csv'
in_data_anom_err = path + '\\in_data\\Regional_anomalies_ERROR_ref_period_2009-2018.csv'

out_dir = path + '\\out_data\\OCE_files_by_region\\'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

##### 2. DEFINE GEODETIC DATA ######

# # read global geodetic mass-balance data from csv into dataframe
geo_df= pd.read_csv(in_data_geo, encoding='latin1', delimiter=',', header=0, index_col='WGMS_ID').sort_index()
geo_df.reset_index(inplace=True)

##### 3. DEFINE VARIABLES ######

#  Make a list of all RGI first order regions
reg_lst = geo_df['GLACIER_REGION_CODE'].unique().tolist()
reg_lst.remove('SAN')
reg_lst= reg_lst + ['SA1','SA2']
print(reg_lst)


# all_wgms_id_lst = geo_df['WGMS_ID'].unique().tolist()
# print('Nb glaciers with geodetic obs C3S 2020: '+str(len(all_wgms_id_lst)))

#  Make a list of the full period of interest 1885 to present
yr_df = pd.read_csv(in_data_anom, sep=',', na_values='NaN', index_col=['YEAR'])
yr_lst = yr_df.index.to_list()
print(yr_lst)


min_year_geo_obs = 0 # start year of the period of interest, only geodetic from this year on will be considered, 0 if to use the start of the anomaly period
min_lenght_geo= 5 # minimum lenght in years of the geodetic observations accounted for anomaly calibration
run = 'reg_anomaly'  # define run using the regional anomaly

############################################################################################################################

###### DATA CRUNCHING: Slice glaciological and geodetic dataframe for a region of interest ######

for region in reg_lst:
    region='ANT'
    print('working on region, ', region)

    ## create regional directory for regional glaciers OCE
    out_reg_dir= path + '\\out_data\\' + str(region) + '_oce_by_gla\\'
    if not os.path.exists(out_reg_dir):
        os.mkdir(out_reg_dir)

    ## create regional OCE and sigma OCE empty dataframes
    reg_oce_df = pd.DataFrame(index=yr_lst)
    reg_oce_df.index.name = 'YEAR'
    reg_sig_oce_df = pd.DataFrame(index=yr_lst)
    reg_sig_oce_df.index.name = 'YEAR'

    ## number crunching: select geodetic data for glacier region group
    if region == 'SA1':
        geo_reg_df = geo_df.loc[(geo_df['GLAC_REG2'] == 'SAN-01')]
        # print(geo_reg_df)
        # exit()

    elif region == 'SA2':
        geo_reg_df = geo_df.loc[(geo_df['GLAC_REG2'] == 'SAN-02')]
        # print(geo_reg_df)
        # exit()

    else:
        geo_reg_df = geo_df.loc[(geo_df['GLACIER_REGION_CODE'] == str(region))]
    # print(geo_reg_df)

    ## create a list of wgms_ids belonging to the region group
    reg_wgms_id_lst= geo_reg_df['WGMS_ID'].unique().tolist()
    # print(reg_wgms_id_lst)
    # print('Nb glaciers in region ' + str(region) + ' with geodetic obs C3S 2020: ' + str(len(reg_wgms_id_lst)))

    ############################################################################################################################

    ###### CALCULATING OCE: Loop through all glaciers in the region with available geodetic estimates ######

    for wgms_id in reg_wgms_id_lst:
        # wgms_id=491
        print('working on glacier Id, ', wgms_id)

        # # create individual glacier directory
        out_gla_dir = out_reg_dir + 'wgmsId_' + str(wgms_id) + '_oce\\'
        if not os.path.exists(out_gla_dir):
            os.mkdir(out_gla_dir)

        # # read regional anomaly data and uncertainties from csv files into dataframe
        reg_anom_df = pd.read_csv(in_data_anom, sep=',', usecols=['YEAR', region], na_values='NaN')
        reg_anom_df = reg_anom_df.rename(columns={region: wgms_id}).set_index('YEAR')

        reg_anom_err_df = pd.read_csv(in_data_anom_err, sep=',', usecols=['YEAR', region], na_values='NaN')
        reg_anom_err_df = reg_anom_err_df.rename(columns={region: wgms_id}).set_index('YEAR')
        # print(reg_anom_df)
        # print(reg_anom_err_df)

        # # Define period of the complete anomaly series
        val = reg_anom_df[wgms_id].loc[reg_anom_df.first_valid_index()]

        if min_year_geo_obs == 0:
            if reg_anom_df.loc[reg_anom_df[wgms_id]==val].index[0] > 2000: ## For anomalies startig after 2000, use 2000 in order to use the geodetic estimates that start 2000
                min_year = 2000
            else:
                min_year = reg_anom_df.loc[reg_anom_df[wgms_id] == val].index[0]
        else:
            if reg_anom_df.loc[reg_anom_df[wgms_id]==val].index[0] > 2000:
                min_year = 2000
            else:
                min_year = min_year_geo_obs

        max_year = reg_anom_df.index.max()

        # # create geodetic mass balance series and geodetic Dataframe for selected glacier
        geo_mb_gla_df = geo_reg_df.loc[(geo_reg_df['WGMS_ID'] == wgms_id)]
        # print(geo_mb_gla_df)

        # # Select geodetic estimates inside the period of interest and longer than min_lenght_geo

        geo_mb_gla_df = geo_mb_gla_df[['WGMS_ID', 'ini_date', 'fin_date', 'mb_chg_rate', 'sigma_tot_mb_chg']]
        geo_ind_gla_sel_df = geo_mb_gla_df.loc[(geo_mb_gla_df['ini_date'] >= min_year - 2) & (geo_mb_gla_df['fin_date'] <= max_year + 1)]
        # print(geo_ind_gla_sel_df)

        #create empty dataframes for calibrated series, calibrated series uncertainty, sigma geodetic uncertainty and distance to observation period
        cal_series_df = pd.DataFrame(index=yr_lst)
        cal_series_df.index.name='YEAR'

        sig_cal_series_df = pd.DataFrame(index=yr_lst)
        sig_cal_series_df.index.name='YEAR'

        sigma_geo_df = pd.DataFrame(index=yr_lst)
        sigma_geo_df.index.name= 'YEAR'

        dist_geo_df = pd.DataFrame(index=yr_lst)
        dist_geo_df.index.name= 'YEAR'
        dist_geo_df = dist_geo_df.reset_index()

        ###### Calculate the calibrated series #####
        for index, row in geo_ind_gla_sel_df.iterrows():
            if (int(row['fin_date'])-int(row['ini_date'])) >= min_lenght_geo:
                ref_period_geo_obs = range(int(row['ini_date']), int(row['fin_date']))
                ref_anom_period_df = reg_anom_df.loc[reg_anom_df.index.isin(ref_period_geo_obs)]
                avg_ref_anom = ref_anom_period_df.mean()
                ref_anom= reg_anom_df - avg_ref_anom
                # print(ref_anom)

                cal_val = row['mb_chg_rate'] + ref_anom
                cal_series_df['serie_'+str(index)] = cal_val[wgms_id]

                sig_cal_series_df['sig_serie_'+str(index)] = np.sqrt(row['sigma_tot_mb_chg']**2 + reg_anom_err_df**2)
                # print(cal_series_df)
                # print(sig_cal_series_df)

                # # Create geodetic estimate uncertainty dataframe
                sigma_geo_df['serie_' + str(index)]=row['sigma_tot_mb_chg']
                i_date=int(row['ini_date'])
                f_date=int(row['fin_date'])

                # # Create Distance to geodetic observation period dataframe
                dist_geo_df['serie_' + str(index)]= dist_geo_df['YEAR'].apply(lambda row: dis_fil(row, i_date, f_date))
            else:
                pass

        # if cal_series_df.empty == True:
        #     print('No calibrated series for glacier ' + str(wgms_id))

        # print(cal_series_df)

        if cal_series_df.empty == True:
            os.rmdir(out_gla_dir)

        ###### Apply weights to calculate the mean calibrated series #####

        ## Calculate weight related to geodetic estimate uncertainty
        if cal_series_df.empty == False:

            if math.isnan(geo_ind_gla_sel_df['sigma_tot_mb_chg'].max()) == True:
                fill_sigma = 1.0
            else:
                fill_sigma = geo_ind_gla_sel_df['sigma_tot_mb_chg'].max()

            weight_dir = out_gla_dir + 'wgmsId_' + str(wgms_id) + '_weights\\'
            if not os.path.exists(weight_dir):
                os.mkdir(weight_dir)

            sigma_geo_df.fillna(fill_sigma, inplace=True)
            # sigma_geo_df.to_csv(weight_dir + 'Sigma_series_WGMS_ID_' + str(wgms_id) + '.csv', index=False)
            sigma_ratio_df=sigma_geo_df.apply(lambda x: (1 / x))
            wgt1_sigma_df=sigma_ratio_df.div(sigma_ratio_df.sum(axis=1), axis=0) # pass to percentage
            # wgt1_sigma_df.to_csv(weight_dir + 'Weight_percent_Sigma_series_WGMS_ID_' +str(wgms_id)+ '.csv', index=False)

            ## Calculate weight related to distance to the geodetic estimate survey period
            p=2 ## P value for inverse distance weighting
            dist_geo_df=dist_geo_df.set_index('YEAR')
            # dist_geo_df.to_csv(weight_dir + 'Distance_series_WGMS_ID_' + str(wgms_id) + '.csv', index=False)
            inv_dist_df=dist_geo_df.apply(lambda x: (1 / x) ** p)
            wgt2_dist_df=inv_dist_df.div(inv_dist_df.sum(axis=1), axis=0) # pass to percentage
            # wgt2_dist_df.to_csv(weight_dir + 'Weight_percent_Distance_series_WGMS_ID_' +str(wgms_id)+ '.csv', index=False)

            ## Calculate weight related to uncertaity and distance combined
            W1_W2_comb_df=wgt1_sigma_df.add(wgt2_dist_df)

            ##### Calculate MEANS of calibrated series: Artihmetic, Weight_1, Weight_2, Weight_combined #####

            cal_mean_df = pd.DataFrame(index=yr_lst)
            cal_mean_df.index.name='YEAR'
            # print(cal_series_df)

            ### Apply the weights to the calibrated series
            cal_series_W1_df = cal_series_df.mul(wgt1_sigma_df)
            cal_series_W2_df = cal_series_df.mul(wgt2_dist_df)
            cal_series_W1_W2_df=(cal_series_df.mul(W1_W2_comb_df))/2

            ## calibrated series means
            cal_mean_df['MEAN']=cal_series_df.mean(axis=1)
            cal_mean_df['MEAN_sigma_W']=cal_series_W1_df.sum(axis=1, min_count=1)
            cal_mean_df['MEAN_dist_W']=cal_series_W2_df.sum(axis=1, min_count=1)
            cal_mean_df['MEAN_combined_W']=cal_series_W1_W2_df.sum(axis=1, min_count=1)

            ## Plot the different means of the calibrated series
            # # print(cal_mean_df)
            # fig=cal_mean_df.plot()
            # plt.savefig(out_gla_dir + 'wgmsId_' + str(wgms_id) + '_weights\\Fig_diff_mean_cal_series_glacier_id_' + str(wgms_id) + '_' + region + '_Geo_+' + str(min_lenght_geo) + 'years.png')
            # plt.close()

            ## cumulative series of the different means
            cal_mean_cum_df = pd.DataFrame(index=yr_lst)
            cal_mean_cum_df.index.name='YEAR'

            cal_mean_cum_df['Cum_MEAN']=cal_mean_df['MEAN'].cumsum()
            cal_mean_cum_df['Cum_MEAN_sigma_W']=cal_mean_df['MEAN_sigma_W'].cumsum(skipna=True)
            cal_mean_cum_df['Cum_MEAN_dist_W']=cal_mean_df['MEAN_dist_W'].cumsum(skipna=True)
            cal_mean_cum_df['Cum_MEAN_combined_W']=cal_mean_df['MEAN_combined_W'].cumsum(skipna=True)

            ## Plot the cumulative values from the different means of the calibrated series
            # # print(cal_mean_cum_df)
            # fig=cal_mean_cum_df.plot()
            # plt.savefig(out_gla_dir + 'wgmsId_' + str(wgms_id) + '_weights\\Fig_cum_series_diff_mean_glacier_id_' + str(wgms_id) + '_Geo_+' + str(min_lenght_geo) + 'years.png')
            # plt.close()

            ############################################################################################################################

            ### Calculate an save individual glacier OCE uncertainty

            nb_cal= (len(cal_series_df.columns))
            # print(nb_cal)
            # print(sig_cal_series_df)
            # sig_cal_series_df.to_csv(out_gla_dir + 'Sigma_calibrated_series' + str(wgms_id) + '.csv',index=False)

            if nb_cal > 1:
                cal_series_std = cal_series_df.std(axis=1)
                cal_series_var_err = 1.96 * cal_series_std
                # print(cal_series_var_err)

                sig_cal_mean = sig_cal_series_df.mean(axis=1)
                # sig_cal_sqsum=sig_cal_series_df.pow(2).sum(axis=1)
                # sig_cal_series_sum = np.sqrt(sig_cal_sqsum)/nb_cal
                # print(sig_cal_mean)

                # print(sig_cal_series_sum)
                sig_oce = np.sqrt(cal_series_var_err**2 + sig_cal_mean**2)
                # print(sig_oce)
                # exit()

            elif nb_cal == 1:
                sig_oce = sig_cal_series_df

            # if cal_series_df.empty == False:
            reg_sig_oce_df[wgms_id] = sig_oce
            sig_oce.to_csv(out_gla_dir + 'sigma_OCE_wgmsId_' + str(wgms_id) + '_' + region + '.csv')

            ### Plot and save individual glacier OCE

            if nb_cal >= 1:
                cal_series_df['normal_MEAN']=cal_series_df.mean(axis=1)
                cal_series_df['weighted_MEAN']=cal_series_W1_W2_df.sum(axis=1, min_count=1)

                plot_gla_oce(cal_series_df, geo_ind_gla_sel_df, wgms_id, min_year, max_year, min_lenght_geo, region, run, out_gla_dir)
                cal_series_df = cal_series_df.rename(columns={'weighted_MEAN':wgms_id})
                oce_df= cal_series_df[wgms_id]

                oce_df.to_csv(out_gla_dir + 'OCE_wgmsId_' + str(wgms_id) + '_' + region + '.csv')
                # print(oce_df)
                plot_gla_oce_and_unc(cal_series_df, oce_df, geo_ind_gla_sel_df, reg_sig_oce_df, wgms_id, min_year, max_year, min_lenght_geo, region, run,
                             out_gla_dir)

                reg_oce_df[wgms_id] = oce_df


    # print(reg_oce_df)
    # print(reg_sig_oce_df)


    ### Save regional OCEs
    reg_oce_df.to_csv(out_dir + region + '_regional_OCEs.csv')
    reg_sig_oce_df.to_csv(out_dir + region + '_regional_sigma_OCEs.csv')
    exit()
print('.........................................................................................')
print('"The End"')
print('.........................................................................................')
exit()
