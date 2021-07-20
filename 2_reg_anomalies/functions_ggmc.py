"""functions"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

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
    # print(anomaly_ref_df)
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
        # unc_ref_df.loc[unc_ref_df.index.isin(yrs), id] = np.nan
        # unc_ref_df.loc[unc_ref_df.index.isin(yrs), str(id)] = np.nan
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


def plot_epoch_analysis(in_glac_anomaly, in_reg_anomaly, event_year, event_name):
    """This function plots the anomalies of glacier mass balances over a defined reference period."""
    print('Plotting compositing epoch analysis for the eruption of {} ({})'.format(event_name, event_year))

    # set output format
    file_type = '.svg'

    # define name of figure
    out_fig = path + '\\out_data\\' + 'Fig_CEA_' + event_name + '_' + str(event_year) + file_type

    # initialize plot
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot zero balance line, and event line
    ax.axhline(0, color='Grey', linewidth=1)
    ax.axvline(event_year, color='LightGrey', linewidth=1)

    # plot individual glacier anomalies
    for item in glac_anom.columns:
        ax.plot(glac_anom.index, glac_anom[item], color='Grey', linewidth=0.75)

    # plot average regional anomay
    ax.plot(reg_anom.index, reg_anom.values, color='Black', linewidth=4)

    # set labels
    ax.set_title('Glacier mass-change anomaly around the eruption of {} ({})'.format(event_name, event_year),
                 fontsize='medium')
    ax.set_xlabel('Year')
    ax.set_ylabel('Specific mass-change anomaly (m w.e.)')
    ax.set_xticks(np.arange(min(reg_anom.index), max(reg_anom.index) + 1, step=1))

    # save plot
    plt.tight_layout()
    plt.savefig(out_fig)
    print('Plot saved as {}.'.format(out_fig))

    # clear plot
    plt.clf()

    return

