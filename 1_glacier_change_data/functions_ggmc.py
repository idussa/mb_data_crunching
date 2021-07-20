"""functions"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def create_mb_dataframe(in_df, id_lst, yr_lst, mb_field):
    """This function selects mass-balance data from the input file into a dataframe."""
    print('Creating dataframe of {}...'.format(mb_field))

    mb_dict = {}

    for id in id_lst:
        mb_lst = []
        for yr in yr_lst:
            if yr in in_df.loc[in_df['WGMS_ID'] == id]['YEAR'].values:
                mb_lst.append(
                    in_df[(in_df['WGMS_ID'] == id) & (in_df['YEAR'] == yr)].iloc[0][mb_field])
            else:
                mb_lst.append(np.nan)
        mb_dict.update({id: mb_lst})

    mb_df = pd.DataFrame(mb_dict, index=yr_lst)

    print('..done.')
    return mb_df

def calc_anomalies(in_df, ref_period, region):
    """This function calculates the anomalies of glacier mass balances over a defined reference period."""
    print('Calculating anomalies for reference period from {} to {}...'.format(min(ref_period), max(ref_period)))

    # create subset over reference period
    in_ref_df = in_df.loc[in_df.index.isin(ref_period)]

    # create subset with minimal data
    ref_period_df = in_df.loc[in_df.index.isin(ref_period)]
    if region == 'ASN':
        ref_period_ids = ref_period_df.count() > 3
    else:
        ref_period_ids = ref_period_df.count() > 7
    ref_period_id_lst = list(ref_period_ids[ref_period_ids].index)

    # create subset of glacier ids with good data over reference period
    good_ids_in_ref_df = in_ref_df[ref_period_id_lst]
    good_ids_in_df=in_df[ref_period_id_lst]

    # calculate anomaly (x_i-x_avg) for data over reference period
    avg_ref_df = good_ids_in_ref_df.mean()
    anomaly_ref_df = round(good_ids_in_df - avg_ref_df, 0)
    print(anomaly_ref_df)
    print('done.')
    return anomaly_ref_df

def calc_anomalies_unc(in_df, in_unc_df, ref_period, region):
    """This function calculates the uncertainties of glacier mass balances series"""

    # create subset with minimal data
    ref_period_df = in_df.loc[in_df.index.isin(ref_period)]
    if region == 'ASN':
        ref_period_ids = ref_period_df.count() > 3
    else:
        ref_period_ids = ref_period_df.count() > 7
    ref_period_id_lst = list(ref_period_ids[ref_period_ids].index)

    # calculate mb uncertainty of glacier ids with good data over reference period
    unc_ref_df = in_unc_df[ref_period_id_lst]
    reg_unc_mean= np.nanmean(unc_ref_df)

    for id in ref_period_id_lst:
        year_min = in_df[id].first_valid_index()
        yrs= list(range(1885,year_min))
        if unc_ref_df[id].isnull().all():
            unc_ref_df[id].fillna(reg_unc_mean, inplace=True)
        else:
            unc_ref_df[id].fillna(unc_ref_df[id].mean(), inplace=True)
        # unc_ref_df.loc[unc_ref_df.index.isin(yrs), [id]] = np.nan
        unc_ref_df[id].mask(unc_ref_df.index.isin(yrs), np.nan, inplace=True)

    print('done.')
    return unc_ref_df


def plot_reg_anomalies(in_ba_df_lst, path):
    """This function plots the regional mass-balance series and volcanic eruptions over time."""
    print('Plotting mass-balance time series...')

    # set output format
    file_type = '.svg'

    # define name of figure
    out_fig = path + '\\out_data\\Fig_regional_mass-balance_timeseries_' + file_type

    # define glacier region labels
    reg_lst = ['(a) ALA', '(b) WNA', '(c) ACN', '(d) ACS', '(e) GRL',
               '(f) ISL', '(g) SJM', '(h) SCA', '(i) RUA', '(j) ASN',
               '(k) CEU', '(l) CAU', '(m) ASC', '(n) ASW', '(o) ASE',
               '(p) TRP', '(q) SAN', '(r) NZL', '(s) ANT', '(t) Global mean of regional averages']

    # initialize plot
    # fig, ax = plt.subplots(20, 3, sharex='col', figsize=(8, 24))
    fig, ax = plt.subplots(20, 1, sharex='col', figsize=(24, 24))

    # plot regional and global mass-balance series
    for i in range(0, len(in_ba_df_lst)):
        ax[i].axhline(0, color='black', linewidth=0.5, zorder=1)

        # plot annual values
        ax[i].bar(in_ba_df_lst[i].index, in_ba_df_lst[i] / 1000, color='black', zorder=4)

        # set labels
        ax[i].set_ylim(-5, 4.5)
        ax[i].set_xlim(1960, 2020)
        ax[i].set_yticks(range(-3, 4, 3))
        ax[i].text(0.005, 0.85, str('{}'.format(reg_lst[i])), fontsize='small', ha='left', va='center',
                   transform=ax[i].transAxes)
        ax[i].set_ylabel('m w.e.', fontsize='small')


    # save plot
    # plt.tight_layout()
    plt.savefig(out_fig)

    print('Plot saved as {}.'.format(out_fig))

    # clear plot
    plt.clf()

    return

def plot_gla_geo(in_cal_series_df, in_geo_df, glacier_id, min_year, max_year, min_lenght_geo, region, run, gla_dir):
    #add in_geo_df to arguments
    """This function plots the anomalies of glacier mass balances over a defined reference period."""
    print('Plotting calibrated series for glacier = {}'.format(glacier_id))

    # in_geo_df = in_geo_df/1000
    # in_cal_series_df = in_cal_series_df/1000

    # set output format
    file_type = '.pdf'

    # define name of figure
    out_fig = gla_dir + 'Fig1_'+str(region) +'_glacier_'+str(glacier_id)+'_geodetic_sample'+ file_type

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot zero balance line
    ax.axhline(0, color='Grey', linewidth=1)

    # plot geodetic mass change trends
    for index, row in in_geo_df.iterrows():
        color = 'b' if row['mb_chg_rate'] > 0 else 'r'
        x1 = row['ini_date']
        x2 = row['fin_date']
        y1 = row['mb_chg_rate']
        y2 = row['sigma_tot_mb_chg']
        ax.fill([x1, x1, x2, x2], [0, y1, y1, 0], color, alpha=0.1)
        ax.plot([x1, x2], [y1, y1], color, linewidth=0.75, alpha=0.8, solid_capstyle='butt')

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('B$_{geo}$ (m w.e.)', fontsize='large')
    ax.tick_params(labelsize='large')
    ax.text(1952, -3.2, 'c - B$_{geo}$', fontsize=24)

    # save plot
    plt.xlim(1950, 2020)
    plt.ylim((-3.5, 2.0))
    # plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))
    # plt.show()

    # clear plot
    plt.close()

    return


def plot_gla_oce_and_unc(in_cal_series_df, in_oce_df, in_geo_df, in_oce_unc_df,  glacier_id, min_year, max_year, min_lenght_geo, region, run, gla_dir):
    #add in_geo_df to arguments
    """This function plots the anomalies of glacier mass balances over a defined reference period."""
    print('Plotting calibrated series for glacier = {}'.format(glacier_id))

    # set output format
    file_type = '.png'

    # define name of figure
    out_fig = gla_dir + 'Fig2_'+str(region)+'_glacier_'+str(glacier_id)+'_OCE_v2'+ file_type

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot zero balance line
    ax.axhline(0, color='Grey', linewidth=1)

    # plot geodetic mass change trends
    for index, row in in_geo_df.iterrows():
        color = 'b' if row['mb_chg_rate'] > 0 else 'r'
        x1 = row['ini_date']
        x2 = row['fin_date']
        y1 = row['mb_chg_rate']
        y2 = row['sigma_tot_mb_chg']
        # ax.fill([x1, x1, x2, x2], [y1+y2, y1-y2, y1-y2, y1+y2], color='grey', alpha=0.1)
        # ax.plot([x1, x2], [y1, y1], color, linewidth=0.75, alpha=0.8, solid_capstyle='butt')

        # ax.fill([x1, x1, x2, x2], [0, y1, y1, 0], color, alpha=0.1)
        ax.plot([x1, x2], [y1, y1], color, linewidth=0.75, alpha=0.8, solid_capstyle='butt')

    # plot calibrated glaciological series
    for item in in_cal_series_df.columns:
        ax.plot(in_cal_series_df.index, in_cal_series_df[item], color='Silver', alpha=0.7, linewidth=0.7)

    plt.fill_between(in_oce_df.index, in_oce_df + in_oce_unc_df[glacier_id],
                     in_oce_df - in_oce_unc_df[glacier_id], color='grey', alpha=0.3, linewidth=0)

    # plot average calibrated series
    ax.plot(in_oce_df, color='black', alpha=0.7, linewidth=1.5)

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('B$_{OCE}$ (m w.e.)', fontsize='large')
    # ax.set_ylabel('B$_{calibrated}$ (m w.e.)', fontsize='large')
    ax.tick_params(labelsize='large')
    ax.text(1952, -3.2, 'd - B$_{OCE}$ ', fontsize=24)

    # save plot
    plt.xlim(1950, 2020)
    plt.ylim((-3.5, 2.0))
    # plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))
    # plt.show()

    # clear plot
    plt.close()

    return

def date_format(month,year):
    default_month = 6.0
    month=default_month if month > 12.0 else month
    offset = year + (1.0 / 12 * (month - 1))
    return offset

def cum_to_rate(elev_chg_cum,sur_date,ref_date):
    if (sur_date-ref_date) != 0:
        elev_chg_rate = elev_chg_cum/(sur_date-ref_date) if elev_chg_cum != 0 else elev_chg_cum
        return elev_chg_rate
    else:
        return elev_chg_cum

def dis_fil(row, ini_date, fin_date):
    if row > fin_date:
        return row - fin_date + 1
    if row <= ini_date:
        return ini_date + 2 - row
    else:
        return 1.0


