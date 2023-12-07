#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

#
# The get_stations methods for NOAA,NDBC and CONTRAILS are basically all the possible stations
# that we might want to get. These are slowly varying lists. Station ids that do not or no longer
# exist are quietly ignored.
#
import os,sys
import numpy as np
import pandas as pd
import datetime as dt

import harvester.fetch_data as fetch_data

from utilities.utilities import utilities as utilities

class get_obs_stations(object):
    """ 
    Class to establish connection to the lower level code to access noaa, ncbc, or contrails servers and acquire a range of product levels for the set
    of stations input stations

    Input station IDs must be stationids and they are treated as string values.

    This is the Harvester layer that would nornmally be called by ADDA,APSVIZ,Reanalysis work
    """

    # Currently supported sources and products

    SOURCES = ['NOAAWEB', 'NOAA','CONTRAILS','NDBC','NDBC_HISTORIC']
    NOAA_PRODUCTS = ['water_level','hourly_height','predictions','air_pressure','wind_speed']
    CONTRAILS_PRODUCTS = ['river_stream_elevation','river_water_level','river_flow_volume','coastal_water_level','air_pressure']
    NDBC_PRODUCTS=['wave_height','air_pressure','wind_speed']

    # Default to NOAA/NOS WL run

    def __init__(self, source='NOAA',product='water_level', 
                contrails_yamlname=None,knockout_dict=None, station_list_file=None):
        """
        get_obs_stations constructor

        Parameters: 
            contrails_yamlname: <str> Location of the login info for Contrails (optional) Not used if source=NOAA)
            source: <str> Named source. For now either NOAA or CONTRAILS
            product: <str>, product type desired. Values are source-specific
            knockout: <dict> used to remove ranges of time(s) for a given station
            station_list_file: <str> file containing station id data 

         The input OBS.YML file is used only as a placeholder for future considerations
   
        """
        self.source = source.upper()

        if self.source=='NOAA':
            selected_products = self.NOAA_PRODUCTS

        elif self.source=='NOAAWEB':
            selected_products = self.NOAA_PRODUCTS

        elif self.source=='CONTRAILS':
            if contrails_yamlname is None:
                utilities.log.error('For Contrails work an authentication yaml is required. It was missing: Abort')
                sys.exit(1)
            selected_products = self.CONTRAILS_PRODUCTS
        elif self.source=='NDBC' or self.source=='NDBC_HISTORIC':
            selected_products = self.NDBC_PRODUCTS
            utilities.log.info(f'NDBC request source is {self.source}')
        else:
            utilities.log.error(f'No valid source specified {self.source}')
            sys.exit(1)

        # Setup master list of stations. Can override with a list of stationIDs later if you must
        if self.source=='NOAA':
            self.station_list = fetch_data.get_noaa_stations(station_list_file)

        if self.source=='NOAAWEB':
            self.station_list = fetch_data.get_noaa_stations(station_list_file)

        # It is up to the caller to differentiate river from coastal stations 
        if self.source=='CONTRAILS':
            self.station_list=fetch_data.get_contrails_stations(station_list_file)
            print(f'get stations {self.station_list}')

        if self.source=='NDBC':
            self.station_list = fetch_data.get_ndbc_buoys(station_list_file)

        if self.source=='NDBC_HISTORIC':
            self.station_list = fetch_data.get_ndbc_buoys(station_list_file)
       
        utilities.log.debug(f'Station list is {self.station_list}')

        # Specify the desired products
        self.product = product.lower()
        if self.product not in selected_products:
            utilities.log.error(f'Requested product not available {self.product}, Possible choices {selected_products}') 
            sys.exit(1)

        # May Need to get Contrails secrets 
        if self.source=='CONTRAILS':
            self.contrails_yamlname=contrails_yamlname  
            print(f'Contrail secrets {self.contrails_yamlname}')
            utilities.log.info('Got Contrails access information')

        # Specify up any knockout filename based exclusions as a Dict()
        self.knockout_dict = knockout_dict
        if self.knockout_dict is not None:
            utilities.log.info('A knockout dict has been specified')

        utilities.log.info(f'SOURCE to fetch from {self.source}')
        utilities.log.info(f'PRODUCT to fetch is {self.product}')

    def remove_knockout_stations(self, df_station) -> pd.DataFrame:
        """
        Occasionally, a station can report poor results for large ranges of times. This
        method allows the caller to nan those results
        Input should be a dataframe of time indexing and stations as columns.
        The time ranges specified in the args.knockout will be set to Nans inclusively.
        This method can be useful when the indicated statin has historically shown 
        poor unexplained performance over a large window of time

        Uses the class variables:
        self.knockout

        Parameters:
           dataframe: time x stations.
           knockout_dict: <dict> Class variable instantiated at class invocation
        Returns:
           dataframe (time x stations) with self.knockout stations removed
        """
        stations=list(self.knockout_dict.keys()) # How many. We anticipate only one but multiple stations can be handled
        cols=list(map(str,df_station.columns.to_list()))
        # 1 Do any stations exist? If not quietly leave
        if not bool(set(stations).intersection(cols)):
            return df_station
        # 2 Okay for each station loop over time range(s) and NaN away
        utilities.log.debug(f'Stations is {stations}')
        utilities.log.debug(f'dict {self.knockout_dict}')
        print(df_station)
        for station in stations:
            for key, value in self.knockout_dict[station].items():
                df_station.loc[value[0]:value[1]][station]=np.nan
        utilities.log.info('Knockout dict has been applied')
        return df_station

    def override_station_IDs(self, station_list):
        """
        This method allows the user to override the current list of stations to process
        They will be subject to the usual validity testing. We overwrite any station data fetched in the class
    
        Parameters:
            stationlist: list(str) of stationIDs. Overrides any existing list.
        Returns:
            self.station_list: list(str) of stations in the class variable
        """
        if isinstance(station_list, list):
            self.station_list=station_list
            utilities.log.info(f'Manually resetting value of station_list {self.station_list}')
        else:
            utilities.log.error(f'Manual station list can only be a list of ids {self.station_list}')
            sys.exit(1)

    def remove_missingness_stations(self, df_in, max_nan_percentage_cutoff=100)-> pd.DataFrame:
        """
        maxmum percentage of allowable nans in any station. 
        max_nan_percentage_cutoff indicates the MINIMUM percent valid data per station. Corresponds to
        (100-max_nan_percentage_cutoff) as the rate of nans per station
        Note: This should generally only apply to high resolution data (eg 6 min NOAA)

        Parameters:
            df_in: Input dataframe (time x stations).
            max_nan_max_nan_percentage_cutoff: <float>. Percent of allowable nans. (dafault=100)
        Returns:
            df_out: dataframe (time x stationID) with stations kept
        """
        num_times=df_in.shape[0]
        df_counts = pd.DataFrame()
        df_counts['VALUES'] = df_in.count()
        df_counts['PERCENTAGE_VALUES']=100*df_counts/num_times
        df_counts['KEEP_STATION']=df_counts['PERCENTAGE_VALUES']>=100-max_nan_percentage_cutoff
        utilities.log.debug('Station thresholding status {}'.format(df_counts))
        # Now filter away
        df_out = df_in[ df_counts[ df_counts['KEEP_STATION'] == True].index.tolist()]
        exclude_station_list = df_in[ df_counts[ df_counts['KEEP_STATION'] != True].index.tolist()].columns.tolist()
        utilities.log.info(f"Remaining {df_out.shape[1]} stations based on a percent cutoff of {max_nan_percentage_cutoff}")
        utilities.log.debug(f'Removing the following stations because of Max Nan %:{exclude_station_list}')
        return df_out

# Choosing a sampling min is non-trivial and depends on the data product selected. Underestimating is better than overestimating
# Since we will do a rolling averager later followed by a final resampling at 1 hour freq.

    def fetch_station_product(self, time_range, return_sample_min=0, interval=None):
        """
        Fetch the data. The main information is part of the class (sources, products, etc.). However, one must still specify the return_sample_minutes
        to sample the data. This harvesting code will read the raw data for the selected product. Perform an interpolation (it doesn't pad nans), and then
        resample the data at the desired freq (in minutes)

        Parameters:
            time_range: <Tuple> (str,str). Starttime,endtime,inclusive. format='%Y-%m-%d %H:%M:%S' 
            return_sample_min: <int> sampling frequency of the returned, interpolated, data set
            interval: <str>: Values of None or 'h'. Only applied to NOAA. None gives full avail freq
        Returns:
            data: <dataframe> Sampled data of dims (time x stations)
            meta: <dataframe> associated metadata
        """
        starttime = time_range[0]
        endtime=time_range[1]
        utilities.log.info(f'Retrieving data for the time range {starttime} to {endtime}')
        template = "An exception of type {0} occurred."

        interval=interval

        data=np.nan
        meta=np.nan

        time_range=(starttime,endtime)
        if self.source.upper()=='NOAA':
            excludedStations=list()
            # Use default station list
            noaa_stations=self.station_list
            noaa_metadata=f"_{endtime.replace(' ','T')}" 
            try:
                data, meta = fetch_data.process_noaa_stations(time_range, noaa_stations, data_product=self.product, interval=interval, resample_mins=return_sample_min)
            except Exception as ex:
                utilities.log.error(f'NOAA error {template.format(type(ex).__name__, ex.args)}')
                #sys.exit(1)

        if self.source.upper()=='NOAAWEB':
            excludedStations=list()
            # Use default station list
            noaa_stations=self.station_list
            noaa_metadata=f"_{endtime.replace(' ','T')}"
            try:
                data, meta = fetch_data.process_noaaweb_stations(time_range, noaa_stations, data_product=self.product, interval=interval, resample_mins=return_sample_min)
            except Exception as ex:
                utilities.log.error(f'NOAA error {template.format(type(ex).__name__, ex.args)}')
                #sys.exit(1)


        if self.source.upper()=='CONTRAILS':
            contrails_config = utilities.load_config(self.contrails_yamlname)['DEFAULT']
            excludedStations=list()
            meta='_RIVERS' if self.product=='river_water_level' or self.product=='river_stream_elevation' else '_COASTAL'
            contrails_stations=self.station_list
            try:
                # Get default station list
                contrails_metadata=f"{meta}_{endtime.replace(' ','T')}" 
                data, meta = fetch_data.process_contrails_stations(time_range, contrails_stations, contrails_config, data_product=self.product, resample_mins=return_sample_min)
            except Exception as ex:
                utilities.log.error(f'CONTRAILS error {template.format(type(ex).__name__, ex.args)}')
                #sys.exit(1)

        if self.source.upper()=='NDBC':
            excludedStations=list()
            # Use default station list
            ndbc_stations=self.station_list
            ndbc_metadata=f"_{endtime.replace(' ','T')}"
            try:
                data, meta = fetch_data.process_ndbc_buoys(time_range, ndbc_stations, data_product=self.product, resample_mins=return_sample_min)
                utilities.log.info(f'NDBC data {data}')
            except Exception as ex:
                utilities.log.error(f'NDBC process error {template.format(type(ex).__name__, ex.args)}')
                #sys.exit(1)

        if self.source.upper()=='NDBC_HISTORIC':
            excludedStations=list()
            # Use default station list
            ndbc_stations=self.station_list
            ndbc_metadata=f"_{endtime.replace(' ','T')}"
            try:
                data, meta = fetch_data.process_ndbc_historic_buoys(time_range, ndbc_stations, data_product=self.product, resample_mins=return_sample_min)
                utilities.log.info(f'NDBC_HISTORIC data {data}')
            except Exception as ex:
                utilities.log.error(f'NDBC_HISTORIC process error {template.format(type(ex).__name__, ex.args)}')
                #sys.exit(1)

        utilities.log.info(f'Finished with data source {self.source}')

        if self.knockout_dict is not None:
            data = self.remove_knockout_stations(data)
            utilities.log.info('Removing knockouts from acquired observational data')

        #utilities.log.info('Data file {}, meta file {}'.format(data,meta))
        utilities.log.info('Finished')
        return data, meta

        if self.source.upper()=='NDBC_HISTORIC':
            excludedStations=list()
            # Use default station list
            ndbc_stations=self.station_list
            ndbc_metadata='_'+endtime.replace(' ','T')
            try:
                data, meta = fetch_data.process_ndbc_historic_buoys(time_range, ndbc_stations, data_product=self.product, resample_mins=return_sample_min)
                utilities.log.info(f'NDBC_HISTORIC data {data}')
            except Exception as ex:
                utilities.log.error(f'NDBC_HISTORIC process error {template.format(type(ex).__name__, ex.args)}')
                sys.exit(1)
        utilities.log.info(f'Finished with data source {self.source}')

        if self.knockout_dict is not None:
            data = self.remove_knockout_stations(data)
            utilities.log.info('Removing knockouts from acquired observational data')

        utilities.log.info(f'Data file {data}, meta file {meta}')
        utilities.log.info('Finished')
        return data, meta

    def fetch_smoothed_station_product(self, df_in, return_sample_min=60, window=11) -> pd.DataFrame:
        """
        Takes the provided input df, smooths using the indicated window and resamples on the
        input return_sample_min (usually set to an hourly). Lastly, interpolated using a linear model 

        CENTERED window rolling average.

        Parameters:
            df_in: <dataframe> of times x stations.
            window: <int> (def=11) width of window.
            return_sample_min: <int> (def=60) Number of mins to sample on the output product
                upsampling will pad with Nans. Values <=0 indicates no sampling returning
                raw averaged data
        Returns:
            df_smoothed: dataframe (time x stations) smoothed and possibly resampled.
        """
        utilities.log.info(f'Smoothing requested. Window of {window}')
        df_smooth = df_in.rolling(window=window, center=True).mean()
        # Double check if completely empty stations persist
        indlist = df_smooth.loc[df_smooth.isnull().all(1)].index # Only if ALL columns are nan
        df_smooth.loc[indlist] = df_in.loc[indlist] 
        # Optional Resample
        if return_sample_min > 0:
            timesampling = f'{return_sample_min}min'
            df_smooth.interpolate(method='polynomial', order=1, limit=1, inplace=True) 
            df_smooth=df_smooth.resample(timesampling).asfreq()
            df_smooth = self.remove_columns_with_onevalue(df_smooth)
            utilities.log.debug(f'Averaged data has been resampled and then interpolated to {return_sample_min}mins')
        return df_smooth

    def remove_columns_with_onevalue(self, df):
        """
        Depending on the starttime,stoptime and the selection of the missinglness thresholds (esp==100%)
        It is possible a station may only retain a single value. When this happens, if that station is passed to an interpolator
        it will cause a failure. Instead of try/except trapping the interpolation failure, we do a check here to remove
        the offending station as a station with a single value is of little benefit
        Alternatively, the caller could explore changing the time range OR tightening up the threshold
        """
        delete_stations = list()
        for stationid in df.columns.to_list():
            if df[stationid].count() <=1:
                delete_stations.append(stationid)
        if len(delete_stations) > 0:
            utilities.log.warning(f'remove_columns_with_onevalue: Some stations have too little data to interpolate: Removing them {delete_stations}')
        return df.drop(delete_stations,axis=1)

def main(args):
    """
    A simple main method to demonstrate the use of this class
    It assumes the existance of a proper main.yml to get IO information
    It assumes the existance of a proper obs.yml 
    It assumes the existance of a proper contrails.yml (if needed) for accessing contrails
    """

    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'../config','main.yml'))

    # Set up IO env
    utilities.log.info(f"Product Level Working in {os.getcwd()}")

    interval = args.interval

    # Set up time ranges
    if args.starttime is None and args.stoptime is None:
        starttime='2022-02-12 00:00:00' # '%Y-%m-%d %H:%M:%S'
        endtime='2022-02-14 00:00:00'
    else:
        starttime=args.starttime
        endtime=args.stoptime

    # Build informative metadata for filenames
    timein=dt.datetime.strptime(starttime, "%Y-%m-%d %H:%M:%S")
    timeout=dt.datetime.strptime(endtime, "%Y-%m-%d %H:%M:%S")
    iometadata = '_'+timein.strftime('%Y%m%d%H%M')+'_'+timeout.strftime('%Y%m%d%H%M')
    #iometadata=''

    station_list = args.station_list

    # Invoke the class
    # No longer use the obs.yml file inside the class. keep for future considerations
    if args.data_source.upper() == 'NOAA':
        if station_list is None:
            noaa_stations=os.path.join(os.path.dirname(__file__), '../supporting_data', 'CERA_NOAA_HSOFS_stations_V3.1.csv')
        else:
            noaa_stations=station_list
        rpl = get_obs_stations(source=args.data_source, product=args.data_product,
                    contrails_yamlname=None,
                    knockout_dict=None, station_list_file=noaa_stations)
    elif args.data_source.upper() == 'NOAAWEB':
        if station_list is None:
            noaa_stations=os.path.join(os.path.dirname(__file__), '../supporting_data', 'CERA_NOAA_HSOFS_stations_V3.1.csv')
        else:
            noaa_stations=station_list
        rpl = get_obs_stations(source=args.data_source, product=args.data_product,
                    contrails_yamlname=None,
                    knockout_dict=None, station_list_file=noaa_stations)
    elif args.data_source.upper() == 'CONTRAILS':
        if station_list is None: 
            contrails_stations=os.path.join(os.path.dirname(__file__), '../supporting_data','contrails_stations.csv')
        else:
            contrails_stations=station_list
        contrails_yamlname=os.path.join(os.path.dirname(__file__),'../secrets','contrails.yml')
        rpl = get_obs_stations(source=args.data_source, product=args.data_product,
                    contrails_yamlname=contrails_yamlname,
                    knockout_dict=None, station_list_file=contrails_stations)
    elif args.data_source.upper() == 'NDBC' or args.data_source.upper() == 'NDBC_HISTORIC':
        if station_list is None:
            ndbc_stations=os.path.join(os.path.dirname(__file__), '../supporting_data', 'ndbc_buoys.csv')
        else:
            ndbc_stations=station_list
        rpl = get_obs_stations(source=args.data_source, product=args.data_product,
                    contrails_yamlname=None,
                    knockout_dict=None, station_list_file=ndbc_stations)
    else:
        print('No source specified')
        sys.exit(1)

    # Fetch best resolution and no resampling

    data,meta=rpl.fetch_station_product((starttime,endtime), return_sample_min=0, interval=interval )

    # Revert Harvester filling of nans to -99999 back to nans
    data.replace('-99999',np.nan,inplace=True)
    meta.replace('-99999',np.nan,inplace=True)

    # Remove stations with too many nans ( Note Harvester would have previously removed stations that are ALL NANS)
    data_thresholded = rpl.remove_missingness_stations(data, max_nan_percentage_cutoff=100) # (maxmum percentage of allowable nans)

    # Because of the large number of ways for stationsa to get removed/have missingness, etc, AND weith pandas deprecations,
    # This following intersection is required
    meta_list = set(data_thresholded.columns.tolist()).intersection(meta.index.to_list())
    meta_thresholded = meta.loc[meta_list]

    # Apply a moving average (smooth) the data performed the required resampling to the desire rate followed by interpolating
    data_smoothed = rpl.fetch_smoothed_station_product(data_thresholded, return_sample_min=60, window=11)

    # Write the data to disk in a way that mimics ADDA
    # Write selected in Pickle data 
    metapkl = f'./obs_wl_metadata%s.pkl'%iometadata
    detailedpkl = f'./obs_wl_detailed%s.pkl'%iometadata
    smoothpkl = f'./obs_wl_smoothed%s.pkl'%iometadata
    meta_thresholded.to_pickle(metapkl)
    data_thresholded.to_pickle(detailedpkl)
    data_smoothed.to_pickle(smoothpkl)

    # Write selected in JSON format

    # Convert and write selected JSON data
    #metajson = utilities.writePickle(meta_thresholded.index.strftime('%Y-%m-%d %H:%M:%S'),rootdir=rpl.rootdir,subdir=rpl.iosubdir,fileroot='obs_wl_metadata',iometadata=rpl.iometadata)
    data_thresholded.index = data_thresholded.index.strftime('%Y-%m-%d %H:%M:%S')
    data_smoothed.index = data_smoothed.index.strftime('%Y-%m-%d %H:%M:%S')

    metajson = f'./obs_wl_metadata%s.json'%iometadata
    detailedjson = f'./obs_wl_detailed%s.json'%iometadata
    smoothjson = f'./obs_wl_smoothed%s.json'%iometadata
    meta_thresholded.to_json(metajson)
    data_thresholded.to_json(detailedjson)
    data_smoothed.to_json(smoothjson)
    print('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--starttime', action='store', dest='starttime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='NOAA', type=str,
                        help='choose supported data source (case independant) eg NOAA or CONTRAILS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product eg water_level')
    parser.add_argument('--interval', action='store', dest='interval', default=None, type=str,
                        help='Interval request to the fetcher (h): Only used by NOAA/NOAAWEB')
    parser.add_argument('--station_list', action='store', dest='station_list', default=None, type=str,
                        help='Choose a non-default location/filename for a stationlist')
    args = parser.parse_args()
    sys.exit(main(args))
