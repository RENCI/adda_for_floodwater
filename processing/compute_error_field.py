#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

# Class to manage the construction of the ADCIRC-OBSERVATONS error fields
# Check the station df_final error measures. If z-score > (user specified) 3 remove them 

import os, sys
import numpy as np
import pandas as pd
import datetime as dt

from argparse import ArgumentParser
from utilities.utilities import utilities
from scipy import stats

import json

def interpolate_and_sample( diurnal_range, df_in )-> pd.DataFrame:
    """
    This function take an input DataFrame (usually sampled on the hour) and transforms it to a semidiurnal 
    time series based on the desired indexes assembled in durnal_range. The new diurnal_range entries are
    combined with the given DataFrame entries and assigned a NaN value.  Data are sorted and duplicates removed.
    The data are then interpolated using a default model of Linear. The desired diurnal_range of values are
    then returned to the caller.

    Parameters:
        df_in: A DataFrame times (datetime64) x stations
        diurnal_range:  list of datetime64
       
    Returns:
        df_out: DataFrame. (datetime64 x station) on a semidiurnal time sequence

    """
    df_x = pd.DataFrame(diurnal_range, columns=['TIME'])
    df_x.set_index('TIME',inplace=True)
    df_out = df_in.append(df_x) # This merges the two indexes
    df_out = df_out.loc[~df_out.index.duplicated(keep='first')] #Keep first because real data is first on the list
    df_out.sort_index(inplace=True) # this is sorted with intervening nans that need to be imputted
    df_out_int = df_out.interpolate(method='linear')
    df_out_int = df_out_int.loc[diurnal_range]
    df_out_int.index.name='TIME' 
    return df_out_int

def combine_data_to_dict(in_adc,in_obs,in_err, product='WL')->dict:
    """
    Combines the three data frames into a single multiindex object
    This is an optional procedure that maintains compatability with presently
    executing data pipelines

    Parameters:
        in_adc,in_obs,in_diff:  DataFrames (datetime64 x station)
        product: (str) A user selected product header name to add to final dict object

    Returns:
        df_merged_dict: A Dict of the merged 3 data source DataFrame suitable for JSON storage  
    """
    adc,obs,err = in_adc.copy(), in_obs.copy(), in_err.copy()
    adc['SRC']='ADC'
    obs['SRC']='OBS'
    err['SRC']='ERR'
    # Specify current index name
    adc.index.name='TIME'
    obs.index.name='TIME'
    err.index.name='TIME'
    # Do the work
    adc.reset_index(inplace=True)
    obs.reset_index(inplace=True)
    err.reset_index(inplace=True)
    adc.set_index(['TIME','SRC'], inplace=True)
    obs.set_index(['TIME','SRC'], inplace=True)
    err.set_index(['TIME','SRC'], inplace=True)
    df_merged = pd.concat([adc,obs,err])
    utilities.log.debug('Combined data sets into multi-index form')
    data_dict = generate_merged_dict_data(df_merged,product='WL')
    return data_dict

def generate_merged_dict_data(df_merged, product='WL')->dict:
    """
    Reformat the df_merged data into a dict with Stations as the main key
    Create the dict: self.df_merged_dict
    For this class we can expect ADS/OBS/ERR data to all be available.
    Must convert timestamp index to Strings YYYYMMDD HH:MM:SS

    Parameters:
        df_merged: DataFrame. A multi-index,combined 3-data source (adc,obs,diff) object 
        product: (str) A user selected product name to add to final dict object

    Returns:
        df_merged_dict: A Dict of the merged 3 data source DataFrame suitable for JSON storage  
    """
    utilities.log.debug('Begin processing DICT data format')
    variables = ['ADC','OBS','ERR']
    df = df_merged
    df.reset_index(inplace=True) # Remove SRC from the multiindex
    df.set_index(['TIME'], inplace=True)
    df.index = df.index.strftime('%Y-%m-%d %H:%M:%S')
    dictdata = {}
    for variable in variables: 
        df_all = df[df['SRC']==variable]
        dataall = df_all.drop('SRC',axis=1).T
        stations = dataall.index
        cols = dataall.columns.to_list()
        for station in stations:
            val = dataall.loc[station].to_list()
            if station in dictdata.keys():
                dictdata[station].update({variable: {'TIME': cols, product:val}})
            else:
                dictdata[station]={variable: {'TIME': cols, product:val}}
    utilities.log.debug('Constructed DICT time series data')
    return dictdata 

class compute_error_field(object):
    """ 
    This class supports computing station errors between observations and adcirc results.
    Some facility for alignment of times, station exclusins, transformation to a semidiurnal time are supported.
    The final results (errors) are packaged up for passage to an interpolator to build adcirc offset scaler felds.

    Additional time-range seldction may be performed with bound_lo,bound_hi values. If no bounds are supplied, then the intersection of values between adcirc
    and observations will be selected

    The input data structures are generally expected to arise from using the AST/Harvester codes
    
    Parameters:
        obs: A dataframe of index (datetime x stations) from observations (get_obs_stations)
       meta: A dataframe of meta data for the stations generally from calling the harvesterr code (get_obs_stations)
        adc: A dataframe of index (datetime x stations) from, ASGS/ADCIRC (get_adc_stations)

        All time series data are expected to come in as hourly sampling, though other sampling is supported

        Case uses: Generally for ADDA or APSVIZ2, researchers want to assemble data into 12-hour periods and
            examine the average error over each period. In this case one would set parameters as n_hours_per_periods=12, 
            n_averging_periods=4. 
        For example, if you want to take a years worth of hourly data and return the daily averages for all data
        one would: set n_hours_per_periods = 24, and n_averging_periods = n_yearly_length(hours) // 24 
    """

    def __init__(self, obs, adc, meta, n_hours_per_period=12, n_hours_per_tide=12.42, n_pad=1, zthresh=3):
        """
        The set of input parameters are not neccesarily all used all the time. But all are specified
        here as an example

        Parameters:
            obs: DataFrame time series of OBS station.
            adc: DataFrame time series of ADC station.
            meta: DataFrame meta data of OBS station. (names, lat/lon/etc)
            n_hours_per_period: (int). The number of hours (entries) defining a period
            n_hours_per_tide: (int) The semidiurnal period (if needed)
            n_pad: (int) If requestring diurnal tidals, then this adds n_pad of flanking entries to the time series
            z_thresh: (int) If exclusion stations based on Z-scoring, this is the threshold
            exclude_outlier: (bool). Exclude (or not) stations based on Z-score
        """

        self.obs=obs
        self.adc=adc
        self.meta=meta

        self.n_hours_per_period=n_hours_per_period
        self.n_hours_per_tide=n_hours_per_tide # This never needs to change
        self.n_pad=n_pad

        self.z_thresh=zthresh
        utilities.log.debug('Averaging: Hours per period {}'.format(n_hours_per_period))
        utilities.log.debug('Tidal: Hours per tidal period {}, num of padding hours {}'.format(self.n_hours_per_tide, self.n_pad))
        utilities.log.debug('Exclusions: If chosen use Z-thresh of {}'.format(self.z_thresh))

    def _intersection_stations(self):
        """ 
        Here we remove stations not shared. Update inplace
        """
        common_stations = self.obs.columns & self.adc.columns
        self.obs = self.obs[common_stations]
        self.adc = self.adc[common_stations]
        utilities.log.debug('OBS and ADC DataFrames reduced (inplace) to common stations of len {}'.format(len(common_stations)))

    def _intersection_times(self):
        """
        Align/intersect the times between OBS and ADC. This will update OBS,ADC data inplace.
        Index times are expected to be datetime64 (or equiv) entries
        """
        common_times = self.obs.index & self.adc.index
        obs = self.obs.loc[common_times]
        adc = self.adc.loc[common_times]
        self.obs=obs
        self.adc=adc
        utilities.log.debug('OBS and ADC DataFrames reduced (inplace) to common times of len {}'.format(len(common_times)))

#
# The tidal transform step used here is a bit tricky. We begin with the input hourly data (starttime,endtime) from the adc/obs sources.
# Then, beginning at starttime, we build an estimated tidal time series my incrementing by 12/12.42 to the endtime (plus any hourly n_pad if set)
# Then the new times are interpolated and the hourly indexing is removed (except the initial time). 
# Then we seek (in reverse) the number of time steps that comprise full tidal periods (12.42). We take the int( total time steps/12.42 )
# to compute the number of full periods (num_periods). Then selecting rows backwards from the final time we step back num_periods * 12.
# (inclusive). We can step back 12 rows for each tidal period since the new row indexing is implied to be on a 12.42 hour scale
# Then we keep only the full tidal period data.
# You would expect the number of tidal periods to be <  input range/12 unless user specified n_pad > 0

    def _tidal_transform_data(self):
        """
        For ADC and OBS. transform data (inplace). Interpolate on a diurnal time of period n_hours_per_tide = 12.42
        hours (3726 secs). Once done, remove the input hourly data keeping only the newly interpolated 
        results
        Expects intersection of times and bounds to have been applied. No checks are performed
        NOTE timein, timeout are interpolation ranges and not actual data fetches. So extending timeout
        beyond the ADCIRC 'now' state is okay.

        The tidal transform is skipped if not enough data exist to complete a full 12.42 hour period
        """
        n_range = self.adc.index.tolist()
        n_range.sort()
        timein, timeout = n_range[0], n_range[-1]
        normalRange = pd.date_range(str(timein), str(timeout), freq='3600S') 
        n_hours_per_period = self.n_hours_per_period
        n_hours_per_tide = self.n_hours_per_tide
        n_pad = self.n_pad # This is used to push inteprlation end-nans to outside the time bounds

        time_step =  int(3600*n_hours_per_tide/n_hours_per_period) # Always scale to an hour (3600s)
        diurnal_range = pd.date_range(timein, timeout+np.timedelta64(n_pad,'h'), freq=str(time_step)+'S').to_list()
    
        self.adc = interpolate_and_sample( diurnal_range, self.adc )
        self.obs = interpolate_and_sample( diurnal_range, self.obs )

        # Truncate the time ranges to only retain full tidal periods. Update data inplace. 
        num_periods = int(len(self.adc)/12.42)
        if num_periods >= 1:
            num_rows_keep = num_periods * 12
            self.adc = self.adc.tail(num_rows_keep)
            self.obs = self.obs.tail(num_rows_keep)
            utilities.log.debug('Tidal transform: Num retained periods {}, num rows kept {}'.format(num_periods, num_rows_keep))
            utilities.log.debug('Tidal transform: adc data index {}'.format(self.adc.index))
            # Now only keep a complete set of 12.42 ur full-tidal periods whos last time <=timeout
            ttimestart, ttimeend = diurnal_range[0],diurnal_range[-1]
            utilities.log.debug('inplace tidal_transform_data: n_pad {}, n_hours_per_period {}, n_hours_per_tide {}'.format(n_pad, n_hours_per_period, n_hours_per_tide))
            utilities.log.debug('inplace tidal_transform_data: timein {}, timeout {}, trans_time_start, trans_time_end {}'.format(timein, timeout, ttimestart, ttimeend))
        else:
            utilities.log.debug('inplace tidal_transform_data: Not enough data to perform a tidal transform. Skip')

# MIGHT need something special for Hurricanes ?
    def _apply_time_bounds(self, time_range):
        """
        Take the input time Tuple (time_start,time_end) and apply to the OBS,ADC data sets
        Retains only the time ranges that are bounded by the Tuple range (inclusive)
        Update data sets inplace

        Parameters:
            time_range: Tuple of (str,str) ('%Y-%m-%d %H:%M:%S'): DataFrame time series of OBS station.

        Returns:
             Update OBS,ADC (inplace)
        """
        dformat='%Y-%m-%d %H:%M:%S'
        try:
            bound_lo = dt.datetime.strptime(time_range[0], dformat)
            bound_hi = dt.datetime.strptime(time_range[1], dformat)
        except (ValueError,TypeError):
            try:
                outid = int(test_val)
                bound_lo = time_range[0]
                bound_hi = time_range[1]
            except ValueError:
                utilities.log.error('_apply_time_bounds Input time range is wrong. Must be %Y-%m-%d %H:%M:%S or a hurricane advisory number.  Got {}: Abort'.format(time_range))
                sys.exit(1)
        self.adc=self.adc.loc[ (self.adc.index >= bound_lo) & (self.adc.index <= bound_hi) ]
        self.obs=self.obs.loc[ (self.obs.index >= bound_lo) & (self.obs.index <= bound_hi) ]
        self._intersection_times() # Should not be needed here but just in case
        utilities.log.debug('New adc time lo {}, New adc time hi {}'.format( min(self.adc.index).strftime(dformat), max(self.adc.index.strftime(dformat))))

# Statistical stuff

    def _remove_station_outliers(self):
        """ 
        Process the Summary error data and look for stations with 
        outlier status. This may imply an outside forcing function 
        and warrents removal of the station. 
        Use self.z_thresh to remove stations based on Z-score.

        Removal of station from the summary data are performed inplace

        Return
            droplist: (list). List of stations to remove from OBS/ADC
        """
        z = np.abs(stats.zscore(self.df_final['mean'].dropna(axis=0)))> self.z_thresh
        droplist = self.df_final.dropna(axis=0)[z].index.tolist()
        self.df_final.drop(droplist,axis=0,inplace=True)
        utilities.log.debug('Requested check for station error outlier status. {} stations removed using a zthresh of {}.'.format(len(droplist),self.z_thresh))

    def _compute_and_average_errors(self):
        """
        Average over periods in multiples of n_hours_per_period steps
        This method does not require a semidiurnal treansformation to have been perfrormed

        For nowcast-type data Period 0 is the closest period to the the input timemark. For forecast data
        the time bounds are flipped and so the FIRST Period (most negative) is the closest to the time timemark

        We expect that time bounds have been set
        """
        # Get the raw diffs, NaNs and all. Some applications want this

        self.diff = self.adc - self.obs

        ## Prep diff for some averaging
        self.diff.reset_index(inplace=True)
        # Groupup perdios and get group averages
        self.df_cycle_aves = self.diff.groupby(self.diff.index // self.n_hours_per_period).mean().T
        # Invert columns labels. Period 0 is the newest, then -1, -2, etc.
        num_cols = len(self.df_cycle_aves.columns)-1
        self.df_cycle_aves.columns = ['Period_'+str(x) for x in self.df_cycle_aves.columns-num_cols]
        num_periods = self.df_cycle_aves.shape[1]
        self.df_cycle_aves = pd.concat([self.meta[['LON','LAT','NAME']], self.df_cycle_aves],axis=1)
        self.df_cycle_aves.index.name='STATION'
        # Get overall mean and average as well but not for the TIME entry
        self.df_final = pd.DataFrame([self.diff.drop('TIME',axis=1).mean(), self.diff.drop('TIME',axis=1).std()]).T
        self.df_final.columns=['mean','std']
        self.df_final = pd.concat([self.meta[['LON','LAT','NAME']], self.df_final],axis=1)
        self.df_final.index.name='STATION'
        # Reset diff format back
        self.diff.set_index('TIME',inplace=True)
        utilities.log.debug('Completed Cycle averaging. Num of periods {}'.format(num_periods))

##
## Simple test based on data previously stored into files
##

def main(args):

    config = utilities.init_logging(subdir=None, config_file='../config/main.yml')

    meta = args.obsmeta
    obsf = args.obsdata
    adcf = args.adcdata

    adc = pd.read_pickle(args.adcdata)
    obs = pd.read_pickle(args.obsdata) 
    meta = pd.read_pickle(args.obsmeta)

    comp_err = compute_error_field(obs, adc, meta) # All default params
    comp_err._intersection_stations()
    comp_err._intersection_times()
    comp_err._tidal_transform_data()
    # Demonstrate limiting the time_range
    dformat='%Y-%m-%d %H:%M:%S'
    time_start='2022-02-17 00:00:00'
    time_stop='2022-02-19 00:00:00'
    time_range = (time_start, time_stop)
    comp_err._apply_time_bounds(time_range)
    comp_err._compute_and_average_errors()

##
## Dump pkl files to local disks
##
    comp_err.diff.to_pickle('tideTimeErrors.pkl')
    comp_err.df_final.to_csv('stationSummaryAves.csv')
    comp_err.df_cycle_aves.to_csv('stationPeriodAves.csv')
    data_dict = combine_data_to_dict(comp_err.adc,comp_err.obs,comp_err.diff, product='WL')
    with open('adc_obs_error_merged.json', 'w') as fp:
        json.dump(data_dict, fp)
    print('Finished')
 
if __name__ == '__main__':
    parser = ArgumentParser(description=main.__doc__)
    parser.add_argument('--obsmeta', action='store', dest='obsmeta',default=None, help='FQFN to obs_wl_metadata.pkl', type=str)
    parser.add_argument('--obsdata', action='store', dest='obsdata',default=None, help='FQFN to obs_wl_smoothed.pkl', type=str)
    parser.add_argument('--adcdata', action='store', dest='adcdata',default=None, help='FQFN to adc_wl.pkl', type=str)
    parser.add_argument('--extraExpDir', action='store', dest='extraExpDir', default=None, help='Subdir to store files', type=str)
    args = parser.parse_args()
    sys.exit(main(args))

