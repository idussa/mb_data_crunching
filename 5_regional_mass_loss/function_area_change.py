# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# glacier_area.py
#
# Author: M. Zemp
# Date: 07 Feb 2019
# Last changes: 23 July 2019
#
# Scripted for Python 2.7
#
# Description:
# This script contains a functions to calculate the glacier area
# for all first-order regions and year based on RGI 6.0 and area
# change rates compiled and used by Zemp et al. (2019)
#
# Input: data lists
#
# Return: data dictionary
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# import required packages


def read_change_rate_data():
    """Read in regional glacier area and current change rates"""

    # for each region, define area and reference year (from RGI 6.0) and change rates (%) from Zemp et al. 2019
    area_data = [
        {'REGION': 'ACN', 'AREA_km2': 105111, 'REF_YEAR': 2000, 'CHANGE_RATE_%': -0.07},
        {'REGION': 'ACS', 'AREA_km2': 40888, 'REF_YEAR': 2000, 'CHANGE_RATE_%': -0.08},
        {'REGION': 'ALA', 'AREA_km2': 86725, 'REF_YEAR': 2009, 'CHANGE_RATE_%': -0.48},
        {'REGION': 'ANT', 'AREA_km2': 132867, 'REF_YEAR': 1989, 'CHANGE_RATE_%': -0.27},
        {'REGION': 'ASC', 'AREA_km2': 49303, 'REF_YEAR': 2003, 'CHANGE_RATE_%': -0.18},
        {'REGION': 'ASE', 'AREA_km2': 14734, 'REF_YEAR': 2003, 'CHANGE_RATE_%': -0.47},
        {'REGION': 'ASN', 'AREA_km2': 2410, 'REF_YEAR': 2011, 'CHANGE_RATE_%': -0.43},
        {'REGION': 'ASW', 'AREA_km2': 33568, 'REF_YEAR': 2003, 'CHANGE_RATE_%': -0.36},
        {'REGION': 'CAU', 'AREA_km2': 1307, 'REF_YEAR': 2001, 'CHANGE_RATE_%': -0.53},
        {'REGION': 'CEU', 'AREA_km2': 2092, 'REF_YEAR': 2003, 'CHANGE_RATE_%': -0.93},
        {'REGION': 'GRL', 'AREA_km2': 89717, 'REF_YEAR': 2001, 'CHANGE_RATE_%': -0.82},
        {'REGION': 'ISL', 'AREA_km2': 11060, 'REF_YEAR': 2000, 'CHANGE_RATE_%': -0.36},
        {'REGION': 'NZL', 'AREA_km2': 1162, 'REF_YEAR': 1978, 'CHANGE_RATE_%': -0.69},
        {'REGION': 'RUA', 'AREA_km2': 51592, 'REF_YEAR': 2006, 'CHANGE_RATE_%': -0.08},
        {'REGION': 'SAN', 'AREA_km2': 29429, 'REF_YEAR': 2000, 'CHANGE_RATE_%': -0.18},
        {'REGION': 'SCA', 'AREA_km2': 2949, 'REF_YEAR': 2002, 'CHANGE_RATE_%': -0.27},
        {'REGION': 'SJM', 'AREA_km2': 33959, 'REF_YEAR': 2001, 'CHANGE_RATE_%': -0.26},
        {'REGION': 'TRP', 'AREA_km2': 2341, 'REF_YEAR': 2000, 'CHANGE_RATE_%': -1.19},
        {'REGION': 'WNA', 'AREA_km2': 14524, 'REF_YEAR': 2006, 'CHANGE_RATE_%': -0.54},
    ]

    return area_data


def calc_current_area_li(sample_year):
    """Calculate regional glacier area for a given year based on linear equation"""
    print '..............................'
    print 'Calculating regional glacier area for {} based on linear equation.'.format(sample_year)

    # read in area data
    regional_area = read_change_rate_data()

    print regional_area[0]
    print sample_year
    # create list to hold results
    regional_area_sample_year = []

    # loop through regions
    for reg in regional_area:
        rgi_area = int(reg['AREA_km2'])
        change_rate_km2 = reg['CHANGE_RATE_%'] / 100 * rgi_area

        # calculate glacier area for sample year
        area_sample_year = int(round(rgi_area + (sample_year - reg['REF_YEAR']) * change_rate_km2, 0))

        # add results as dict to list
        entry = {
            'REGION': reg['REGION'],
            'SAMPLE_YEAR': sample_year,
            'SAMPLE_YEAR_AREA_km2': area_sample_year
        }
        regional_area_sample_year.append(entry)

    return regional_area_sample_year


def calc_current_area_ci(sample_year):
    """Calculate regional glacier area for a given year based on compound interest equation"""
    print '..............................'
    print 'Calculating regional glacier area for {} based on compound interest equation.'.format(sample_year)

    # read in area data
    regional_area = read_change_rate_data()

    print regional_area[0]
    print sample_year
    # create list to hold results
    regional_area_sample_year = []

    # loop through regions
    for reg in regional_area:
        rgi_area = int(reg['AREA_km2'])
        # test if sample year equals RGI reference year
        if sample_year == reg['REF_YEAR']:
            # print '{} =='.format(reg['REF_YEAR'])
            area_sample_year = rgi_area
        # test if sample year is after RGI reference year
        elif sample_year > reg['REF_YEAR']:
            # print '{} >'.format(reg['REF_YEAR'])
            n_yrs = sample_year - reg['REF_YEAR']
            change_rate = reg['CHANGE_RATE_%']
            area_sample_year = calc_compound_interest(rgi_area, change_rate, n_yrs)
        # then sample year must be before RGI reference year
        else:
            # print '{} <'.format(reg['REF_YEAR'])
            n_yrs = reg['REF_YEAR'] - sample_year
            change_rate = -1 * reg['CHANGE_RATE_%']
            area_sample_year = calc_compound_interest(rgi_area, change_rate, n_yrs)
        # add results as dict to list
        entry = {
            'REGION': reg['REGION'],
            'SAMPLE_YEAR': sample_year,
            'SAMPLE_YEAR_AREA_km2': area_sample_year
        }
        regional_area_sample_year.append(entry)

    return regional_area_sample_year


def calc_compound_interest(ko, p, n):
    """Calculate compound interest for n years and interest p"""
    kz = int(round((ko * (1 + float(p) / 100)**n), 0))
    return kz
