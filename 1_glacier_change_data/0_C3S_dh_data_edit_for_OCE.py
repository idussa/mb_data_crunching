"""
C3S_dh_data_edit_for_OCE.py

Author: idussail2
Date: 2/22/2021
Last changes: Enter new date here

Scripted for Python 3.7

Description:
This script reads glacier-wide elevation change data from the WGMS FoG database
and transform it into specific glacier mass balance change
Mass balance uncertainties include dh/dt uncertainty + density conversion uncertainty

Input: C3S Glacier elevation change series
    _C3S_ELEVATION_CHANGE_DATA_G3P_20200824.csv

Return: GEO_MASS_BALANCE_DATA_20200824.csv
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
path = os.path.dirname(os.path.abspath(__file__))

# read file  global geodetic elevation change data (FoG format)

c3s_geo_data = path + '\\in_data\\_C3S_ELEVATION_CHANGE_DATA_20200824.csv'
geo_data_df = pd.read_csv(c3s_geo_data, delimiter=',', header=0)

c3s_geo_series = path + '\\in_data\\_C3S_ELEVATION_CHANGE_SERIES_20200824.csv'
geo_series_df = pd.read_csv(c3s_geo_series, delimiter=',', header=0, usecols=[2,6])

in_data_geo= geo_data_df.merge(geo_series_df.drop_duplicates(subset=['WGMS_ID']), how='left')

#################################
#################################
#### 1. EDIT GEODETIC DATA
#################################
#################################

# # read global geodetic mass-balance data from csv into dataframe
input_geo= pd.read_csv(in_data_geo, encoding='latin1', delimiter=',', header=0, index_col='WGMS_ID').sort_index()
# print(input_geo.columns)

## Get reference date (ini)
input_geo['str_ref_date']= input_geo['REFERENCE_DATE'].astype(str)
input_geo['ref_year']=input_geo['str_ref_date'].str.slice(0, 4)
input_geo['ref_year']=pd.to_numeric(input_geo['ref_year'], downcast="float")
input_geo['ref_month']=input_geo['str_ref_date'].str.slice(4, 6)
input_geo['ref_month']=pd.to_numeric(input_geo['ref_month'], downcast="float")

## Get survey date (fin)
input_geo['str_sur_date']= input_geo['SURVEY_DATE'].astype(str)
input_geo['sur_year']=input_geo['str_sur_date'].str.slice(0, 4)
input_geo['sur_year']=pd.to_numeric(input_geo['sur_year'], downcast="float")
input_geo['sur_month']=input_geo['str_sur_date'].str.slice(4, 6)
input_geo['sur_month']=pd.to_numeric(input_geo['sur_month'], downcast="float")

## Change date into decimal format
input_geo['ini_date']=input_geo.apply(lambda x: date_format(x['ref_month'], x['ref_year']), axis=1)
input_geo['fin_date']=input_geo.apply(lambda x: date_format(x['sur_month'], x['sur_year']), axis=1)

## Transform cumulative elevation changes to rate
input_geo['elevation_chg_rate']=input_geo.apply(lambda x: cum_to_rate(x['ELEVATION_CHANGE'], x['fin_date'], x['ini_date']), axis=1)
input_geo['sigma_elevation_chg']=input_geo.apply(lambda x: cum_to_rate(x['ELEVATION_CHANGE_UNC'], x['fin_date'], x['ini_date']), axis=1)

## Transform elevation change to specific mass balance: Apply density conversion factor
f_dens = 0.85
sig_dens = 0.06

input_geo['mb_chg_rate']= input_geo['elevation_chg_rate'] * f_dens
input_geo['sigma_obs_mb_chg']= input_geo['sigma_elevation_chg'] * f_dens

## Calculate combined mass balance uncertainty: sigma dh + density transformation

mb_rate = input_geo['mb_chg_rate']
sigma_mb = input_geo['sigma_obs_mb_chg']
sigma_mean = input_geo['sigma_obs_mb_chg'].mean()
sigma_mb = sigma_mb.fillna(sigma_mean)

input_geo['sigma_tot_mb_chg'] = abs(mb_rate) * np.sqrt((sigma_mb/mb_rate)**2 +(sig_dens/f_dens)**2)

geo_df =input_geo.drop(['str_ref_date', 'ref_year', 'ref_month', 'str_sur_date', 'sur_year', 'sur_month',
                        'ELEVATION_CHANGE', 'ELEVATION_CHANGE_UNC', 'sigma_obs_mb_chg'], axis=1)

geo_df.reset_index(inplace=True)

geo_df.to_csv(path+'\\in_data\\GEO_MASS_BALANCE_DATA_20200824.csv', index=False)
