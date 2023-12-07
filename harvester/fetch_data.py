#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

import os,sys
import pandas as pd
import datetime as dt
import math

#from harvester.fetch_station_data import noaa_web_fetch_data, noaanos_fetch_data, contrails_fetch_data, ndbc_fetch_data, ndbc_fetch_historic_data
from harvester.fetch_station_data import noaa_web_fetch_data, noaanos_fetch_data, contrails_fetch_data
from utilities.utilities import utilities as utilities

##
## Some basic functions that wrap the underlying calls
##

# Currently supported sources
SOURCES = ['NOAA','CONTRAILS','NDBC','NDBC_HISTORIC','NOAAWEB']


def intersect_stations(df_data,df_meta):
    """
    Occasionally, (esp for NDBC data) the returned list of data is not the same as the
    returned list of metadata. This will cause weird plotting scenarios for apsviz
    So here we check end align as needed
    """
    data_stations=df_data.columns.tolist()
    meta_stations=df_meta.index.tolist()
    stations=list(set(data_stations) & set(meta_stations))
    utilities.log.debug(f' Total number of OBS intersected stations {len(stations)}')
    return df_data[stations], df_meta.loc[stations]


#def get_noaa_stations(fname='./config/noaa_stations.txt'):
def get_noaa_stations(fname=None):
    """
    Simply read a CSV file containing stations under the header of stationid
    Expected format is
        serial_nr, stationid

    Parameters:
        fname: <str> full path to a valid stationid file

    Returns:
        noaa_stations: list(str). List of valid noaa station ids 
    """
    if fname is None:
        utilities.log.error('No NOAA station file assigned: Abort')
        sys.exit(1)
    df = pd.read_csv(fname, index_col=0, header=0, skiprows=[1],dtype=str)
    noaa_stations_list = df["stationid"].to_list() 
    noaa_stations=[word.rstrip() for word in noaa_stations_list] 
    return noaa_stations

def get_contrails_stations(fname=None):
    """
    Simply read a CSV file containing stations under the header of stationid
    A convenience method to fetch river guage lists. 
    Contrails data

    Expected format is
        serial_nr, stationid
    Parameters:
        fname: <str> full path to a valid stationid file

    Returns:
        contrails_stations: list(str). List of valid contrails station ids 
    """
    if fname is None:
        utilities.log.error('No Contrails station file assigned: Abort')
        sys.exit(1)
    df = pd.read_csv(fname, index_col=0, header=0, skiprows=[1],dtype=str)
    contrails_stations_list = df["stationid"].to_list()
    contrails_stations=[word.rstrip() for word in contrails_stations_list] 
    return contrails_stations

def get_ndbc_buoys(fname=None):
    """
    Read a list of buoy data stations. These data are more complicated because
    the NDBC reader doesnt easily provide the location and state information. Thus
    we expect this input file to carry that infomration.
    
    Expected format is
        serial_nr, stationid, location, state
    
    Parameters:
        fname: <str> full path to a valid stationid file

    Returns:
     a list of tuples:
        station_tuples: [(station,location,state), (station, location, state), etc]

    """
    if fname is None:
        utilities.log.error('No NDBC station file assigned: Abort')
        sys.exit(1)
    df_buoys = pd.read_csv(fname,index_col=0, header=0, skiprows=[1],dtype=str)
    # Many buoys have NO STATE affiliation to replace nans with NONE
    df_buoys['state']=df_buoys['state'].fillna('NONE')
    list_stations = df_buoys["stationid"].tolist()
    list_locations = df_buoys["location"].tolist()
    list_states = df_buoys["state"].tolist()
    station_tuples = tuple(zip(list_stations, list_locations, list_states))
    return station_tuples

def choose_common_header_name(product):
    """
    For harvesting, we only want select common data names in the final data time series
    This is complicated by the fact that different sources use different product names. So here
    we manually construct a dictionary of current harvester supported data products

    Input:
        product: <str> input product name

    Returns:
        name: <str >selected common name
    """
    product_name_maps={
        'water_level': 'water_level',
        'wave_height': 'wave_height',
        'predictions': 'water_level', 
        'hourly_height': 'water_level',
        'river_water_level': 'water_level',
        'river_stream_elevation': 'water_level',
        'river_flow_volume': 'flow_volume',
        'coastal_water_level': 'water_level',
        'air_pressure':'air_pressure',
        'pressure':'pressure',
        'wind_speed':'wind_speed'
        }

    if product in product_name_maps.keys():
        name = product_name_maps[product]
        return name.upper()
    else:
        utilities.log.error(f'choose_common_header_name. No such product name {product}')
        sys.exit(1)

def format_data_frames(df, product) -> pd.DataFrame:
    """
    A Common formatting used by all sources. 
    """
    df.index = df.index.strftime('%Y-%m-%dT%H:%M:%S')
    df.reset_index(inplace=True)
    df_out=pd.melt(df, id_vars=['TIME'])
    df_out.columns=('TIME','STATION',choose_common_header_name(product))
    df_out.set_index('TIME',inplace=True)
    return df_out

##
## End functions
##

##
## Globals
##

dformat='%Y-%m-%d %H:%M:%S'
GLOBAL_TIMEZONE='gmt' # Every source is set or presumed to return times in the zone

#TODO abstract this
#PRODUCT='water_level' # Used by all sources regardless of specifici data product selected
# wave_height
# pressure"
# wind speed"

##
## Run stations
##

##
## Added a data aligner. I have seens some cases (esp in the NDBC), where we can find data but no metadata
## Here we now keep data common to both,lest we get plots with no metadata
##

def process_noaa_stations(time_range, noaa_stations, interval=None, data_product='water_level', resample_mins=15 ):
    """
    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set

    Parameters:
        time_range: <tuple> (<str>,<str>). Input time range ('%Y-%m-%dT%H:%M:%S)
        noaa_stations: list(str). List of desired NOAA stations
        interval: <str> (def=None) A NOAA specific interval setting
        data_product: <str >(def water_level). A generic AST named data product ( Not the True NOAA data product name) 
        resample_mins: <int> Returned time series with a sampling of resample_mins

    Returns:
        df_noaa_data: DataFrame (time x station)
        df_noaa_meta: DataFrame (station x metadata)
    """
    # Fetch the data
    noaa_products=['water_level', 'predictions', 'hourly_height', 'air_pressure', 'wind_speed']
    try:
        if not data_product in noaa_products:
            utilities.log.error(f'NOAA: data product can only be {noaa_products}')
            #sys.exit(1)
        noaanos = noaanos_fetch_data(noaa_stations, time_range, product=data_product, interval=interval, resample_mins=resample_mins)
        df_noaa_data = noaanos.aggregate_station_data()
        df_noaa_meta = noaanos.aggregate_station_metadata()
        df_noaa_data,df_noaa_meta = intersect_stations(df_noaa_data.copy(),df_noaa_meta.copy())
        df_noaa_meta.index.name='STATION'
    except Exception as e:
        utilities.log.error(f'Error: NOAA: {e}')
    return df_noaa_data, df_noaa_meta

def process_noaaweb_stations(time_range, noaa_stations, interval=None, data_product='water_level', resample_mins=15 ):
    """
    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set

    Parameters:
        time_range: <tuple> (<str>,<str>). Input time range ('%Y-%m-%dT%H:%M:%S)
        noaa_stations: list(str). List of desired NOAA stations
        interval: <str> (def=None) A NOAA specific interval setting
        data_product: <str >(def water_level). A generic AST named data product ( Not the True NOAA data product name) 
        resample_mins: <int> Returned time series with a sampling of resample_mins

    Returns:
        df_noaa_data: DataFrame (time x station)
        df_noaa_meta: DataFrame (station x metadata)
    """
    # Fetch the data
    noaa_products=['water_level', 'predictions', 'hourly_height', 'air_pressure', 'wind_speed']
    try:
        if not data_product in noaa_products:
            utilities.log.error(f'NOAA WEB: data product can only be {noaa_products}')
            #sys.exit(1)
        noaanos = noaa_web_fetch_data(noaa_stations, time_range, product=data_product, interval=interval, resample_mins=resample_mins)
        df_noaa_data = noaanos.aggregate_station_data()
        df_noaa_meta = noaanos.aggregate_station_metadata()
        df_noaa_data,df_noaa_meta = intersect_stations(df_noaa_data.copy(),df_noaa_meta.copy())
        df_noaa_meta.index.name='STATION'
    except Exception as e:
        utilities.log.error(f'Error: NOAA WEB: {e}')
    return df_noaa_data, df_noaa_meta

def process_contrails_stations(time_range, contrails_stations, authentication_config, data_product='river_water_level', resample_mins=15 ):
    """
    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set

    Parameters:
        time_range: <tuple> (<str>,<str>). Input time range ('%Y-%m-%dT%H:%M:%S)
        contrails_stations: list(str). List of desired Contrails stations
        authentication_config: <dict>. A Contrails specific authorization dict
        data_product: <str> (def river_water_level). A generic AST named data product ( Not the True Contrails data product name) 
        resample_mins: <int> Returned time series with a sampling of resample_mins

    Returns:
        df_contrails_data: DataFrame (time x station)
        df_contrails_meta: DataFrame (station x metadata)
    """
    # Fetch the data
    contrails_product=['river_flow_volume','river_water_level','coastal_water_level', 'air_pressure', 'river_stream_elevation']
    try:
        if data_product not in contrails_product:
            utilities.log.error(f'Contrails data product can only be: {contrails_product} was {data_product}')
            #sys.exit(1)
        contrails = contrails_fetch_data(contrails_stations, time_range, authentication_config, product=data_product, owner='NCEM', resample_mins=resample_mins)
        df_contrails_data = contrails.aggregate_station_data()
        df_contrails_meta = contrails.aggregate_station_metadata()
        df_contrails_data,df_contrails_meta = intersect_stations(df_contrails_data.copy(),df_contrails_meta.copy())

        df_contrails_meta.index.name='STATION'
    except Exception as e:
        utilities.log.error(f'Error: CONTRAILS: {e}')
    return df_contrails_data, df_contrails_meta

#def process_ndbc_buoys(time_range, ndbc_buoys, data_product='wave_height', resample_mins=15 ):
#    """
#    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set
#    
#    Parameters:
#        time_range: <tuple> (<str>,<str>). Input time range ('%Y-%m-%dT%H:%M:%S)
#        ndbc_buoys: list(str). List of desired NDBC buoys
#        data_product: <str> (def wave_height). An AST named data product ( Not the True NDBC data product name) 
#        resample_mins: (int) Returned time series with a sampling of resample_mins
#
#    Returns:
#        df_ndbc_data: DataFrame (time x station)
#        df_ndbc_meta: DataFrame (station x metadata)
#    """
#    # Fetch the data
#    ndbc_products=['wave_height', 'air_pressure', 'wind_speed']
#    try:
#        if not data_product in ndbc_products:
#            utilities.log.error(f'NDBC: data product can only be {ndbc_products}')
#            #sys.exit(1)
#        ndbc = ndbc_fetch_data(ndbc_buoys, time_range, product=data_product, resample_mins=resample_mins)
#        df_ndbc_data = ndbc.aggregate_station_data()
#        df_ndbc_meta = ndbc.aggregate_station_metadata()
#        df_ndbc_data,df_ndbc_meta = intersect_stations(df_ndbc_data.copy(),df_ndbc_meta.copy())
#        df_ndbc_meta.index.name='STATION'
#    except Exception as e:
#        utilities.log.error(f'Error: NEW NDBC: {e}')
#    return df_ndbc_data, df_ndbc_meta
#
#def process_ndbc_historic_buoys(time_range, ndbc_buoys, data_product='wave_height', resample_mins=15 ):
#    """
#    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set
#    
#    Parameters:
#        time_range: <tuple> (<str>,<str>). Input time range ('%Y-%m-%dT%H:%M:%S)
#        ndbc_buoys: list(str). List of desired NDBC buoys
#        data_product: <str> (def wave_height). An AST named data product ( Not the True NDBC data product name) 
#        resample_mins: <int> Returned time series with a sampling of resample_mins
#
#    Returns:
#        df_ndbc_data: DataFrame (time x station)
#        df_ndbc_meta: DataFrame (station x metadata)
#    """
#    # Fetch the data
#    ndbc_products=['wave_height', 'air_pressure', 'wind_speed']
#    try:
#        if not data_product in ndbc_products:
#            utilities.log.error(f'NDBC: data product can only be {ndbc_products}')
#            #sys.exit(1)
#        ndbc = ndbc_fetch_historic_data(ndbc_buoys, time_range, product=data_product, resample_mins=resample_mins)
#        df_ndbc_data = ndbc.aggregate_station_data()
#        df_ndbc_meta = ndbc.aggregate_station_metadata()
#        df_ndbc_data,df_ndbc_meta = intersect_stations(df_ndbc_data.copy(),df_ndbc_meta.copy())
#        df_ndbc_meta.index.name='STATION'
#    except Exception as e:
#        utilities.log.error(f'Error: NEW NDBC HISTORIC: {e}')
#    return df_ndbc_data, df_ndbc_meta

def main(args):
    """
    Generally we anticipate inputting a STOPTIME
    Then the STARTTIME is ndays on the past
    """

    #main_config = utilities.init_logging(subdir=None, config_file='../config/main.yml')
    #main_name = args.main_config if args.main_config is not None else os.path.join(os.path.dirname(__file__),'../config','main.yml')
    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'../config','main.yml'))

    if args.sources:
         print('Return list of sources')
         return SOURCES
         sys.exit(0)
    data_source = args.data_source
    data_product = args.data_product

    if data_source.upper() in SOURCES:
        utilities.log.info('Found selected data source {}'.format(data_source))
    else:
        utilities.log.error('Invalid data source {}'.format(data_source))
        sys.exit(1)

    utilities.log.info(f'Input product {data_product}')

    # Setup times and ranges
    if args.stoptime is not None:
        time_stop=dt.datetime.strptime(args.stoptime,dformat)
    else:
        time_stop=dt.datetime.now()

    # Do we want a hard starttime?
    if args.starttime is None:
        time_start=time_stop+dt.timedelta(days=args.ndays) # How many days BACK
    else:
        time_start=dt.datetime.strptime(args.starttime,dformat)

    starttime=dt.datetime.strftime(time_start, dformat)
    endtime=dt.datetime.strftime(time_stop, dformat)
    #starttime='2021-12-08 12:00:00'
    #endtime='2021-12-10 00:00:00'

    utilities.log.info(f'Selected time range is {starttime} to {endtime}, ndays is {args.ndays}')

    # metadata are used to augment filename
    #NOAA/NOS
    if data_source.upper()=='NOAA':
        excludedStations=list()
        time_range=(starttime,endtime) # Can be directly used by NOAA 
        # Use default station list
        noaa_stations=get_noaa_stations(args.station_list) if args.station_list is not None else get_noaa_stations(fname=os.path.join(os.path.dirname(__file__),'../supporting_data','noaa_stations.csv'))
        noaa_metadata=f"_{data_product}_{endtime.replace(' ','T')}"  # +'_'+starttime.replace(' ','T')
        data, meta = process_noaa_stations(time_range, noaa_stations, data_product = data_product)
        df_noaa_data = format_data_frames(data, data_product) # Melt the data :s Harvester default format
        # Output
        # If choosing non-default locations BOTH variables must be specified
        try:
            if args.ofile is not None:
                dataf=f'%s/noaa_stationdata%s.csv'% (args.ofile,noaa_metadata)
                metaf=f'%s/noaa_stationdata_meta%s.csv'% (args.ometafile,noaa_metadata)
            else:
                dataf=f'./noaa_stationdata%s.csv'%noaa_metadata
                metaf=f'./noaa_stationdata_meta%s.csv'%noaa_metadata
            df_noaa_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'NOAA data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: NOAA: Failed Write {e}')
            sys.exit(1)

    # metadata are used to augment filename
    #NOAA/NOS WEB API
    if data_source.upper()=='NOAAWEB':
        excludedStations=list()
        time_range=(starttime,endtime) # Can be directly used by NOAA 
        # Use default station list
        noaa_stations=get_noaa_stations(args.station_list) if args.station_list is not None else get_noaa_stations(fname=os.path.join(os.path.dirname(__file__),'../supporting_data','noaa_stations.csv'))
        noaa_metadata=f"_{data_product}_{endtime.replace(' ','T')}"  # +'_'+starttime.replace(' ','T')
        data, meta = process_noaaweb_stations(time_range, noaa_stations, data_product = data_product)
        df_noaa_data = format_data_frames(data, data_product) # Melt the data :s Harvester default format
        # Output
        # If choosing non-default locations BOTH variables must be specified
        try:
            if args.ofile is not None:
                dataf=f'%s/noaa_stationdata%s.csv'% (args.ofile,noaa_metadata)
                metaf=f'%s/noaa_stationdata_meta%s.csv'% (args.ometafile,noaa_metadata)
            else:
                dataf=f'./noaa_stationdata%s.csv'%noaa_metadata
                metaf=f'./noaa_stationdata_meta%s.csv'%noaa_metadata
            df_noaa_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'NOAA WEB data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: NOAA WEB: Failed Write {e}')
            sys.exit(1)

    #Contrails
    if data_source.upper()=='CONTRAILS':
        # Load contrails secrets
        conf_name = args.config_name if args.config_name is not None else os.path.join(os.path.dirname(__file__),'../secrets','contrails.yml')
        contrails_config = utilities.load_config(conf_name)['DEFAULT']
        utilities.log.info('Got Contrails access information')
        template = "An exception of type {0} occurred."
        excludedStations=list()
        if data_product=='river_water_level' or data_product=='river_flow_volume' or data_product=='river_stream_elevation':
            fname=os.path.join(os.path.dirname(__file__),'../supporting_data','contrails_stations_rivers.csv')
            meta='RIVERS'
        else:
            fname=os.path.join(os.path.dirname(__file__),'../supporting_data','contrails_stations_coastal.csv')
            meta='COASTAL'
        try:
            # Build ranges for contrails ( and noaa/nos if you like)
            time_range=(starttime,endtime) 
            # Get default station list
            contrails_stations=get_contrails_stations(args.station_list) if args.station_list is not None else get_contrails_stations(fname)
            contrails_metadata=f"_{data_product}_{meta}_{endtime.replace(' ','T')}" # +'_'+starttime.replace(' ','T')
            data, meta = process_contrails_stations(time_range, contrails_stations, contrails_config, data_product = data_product )
            df_contrails_data = format_data_frames(data, data_product) # Melt: Harvester default format
        except Exception as ex:
            utilities.log.error(f'CONTRAILS error {type(ex).__name__}, {ex.args}')
            sys.exit(1)
        # If choosing non-default locations BOTH variables must be specified
        try:
            if args.ofile is not None:
                dataf=f'%s/contrails_stationdata%s.csv'% (args.ofile,contrails_metadata)
                metaf=f'%s/contrails_stationdata_meta%s.csv'% (args.ometafile,contrails_metadata)
            else:
                dataf=f'./contrails_stationdata%s.csv'%contrails_metadata
                metaf=f'./contrails_stationdata_meta%s.csv'%contrails_metadata
            df_contrails_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'CONTRAILS data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: CONTRAILS: Failed Write {e}')
            sys.exit(1)

    #NDBC
    if data_source.upper()=='NDBC':
        time_range=(starttime,endtime) # Can be directly used by NDBC
        # Use default station list
        ndbc_stations=get_ndbc_buoys(args.station_list) if args.station_list is not None else get_ndbc_buoys(fname=os.path.join(os.path.dirname(__file__),'../supporting_data','ndbc_buoys.csv'))
        ndbc_metadata=f"_{data_product}_{endtime.replace(' ','T')}" # +'_'+starttime.replace(' ','T')
        data, meta  = process_ndbc_buoys(time_range, ndbc_stations, data_product = data_product)
        df_ndbc_data = format_data_frames(data, data_product) # Melt the data :s Harvester default format
        # Output
        # If choosing non-default locations BOTH variables must be specified
        try:
            if args.ofile is not None:
                dataf=f'%s/ndbc_stationdata%s.csv'% (args.ofile,ndbc_metadata)
                metaf=f'%s/ndbc_stationdata_meta%s.csv'% (args.ometafile,ndbc_metadata)
            else:
                dataf=f'./ndbc_stationdata%s.csv'%ndbc_metadata
                metaf=f'./ndbc_stationdata_meta%s.csv'%ndbc_metadata
            df_ndbc_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'NDBC data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: NDBC: Failed Write {e}')
            sys.exit(1)

    if data_source.upper()=='NDBC_HISTORIC':
        time_range=(starttime,endtime) # Can be directly used by NDBC
        # Use default station list
        ndbc_stations=get_ndbc_buoys(args.station_list) if args.station_list is not None else get_ndbc_buoys(fname=os.path.join(os.path.dirname(__file__),'../supporting_data','ndbc_buoys.csv'))
        ndbc_metadata=f"_{data_product}_{endtime.replace(' ','T')}" # +'_'+starttime.replace(' ','T')
        data, meta  = process_ndbc_historic_buoys(time_range, ndbc_stations, data_product = data_product)
        df_ndbc_data = format_data_frames(data, data_product) # Melt the data :s Harvester default format
        # Output
        # If choosing non-default locations BOTH variables must be specified
        try:
            if args.ofile is not None:
                dataf=f'%s/ndbc_stationdata%s.csv'% (args.ofile,ndbc_metadata)
                metaf=f'%s/ndbc_stationdata_meta%s.csv'% (args.ometafile,ndbc_metadata)
            else:
                dataf=f'./ndbc_stationdata%s.csv'%ndbc_metadata
                metaf=f'./ndbc_stationdata_meta%s.csv'%ndbc_metadata
            df_ndbc_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'NDBC data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: NDBC: Failed Write {e}')
            sys.exit(1)

    utilities.log.info(f'Finished with data source {data_source}')
    utilities.log.info('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--ndays', action='store', dest='ndays', default=-2, type=int,
                        help='Number of look-back days from stoptime (or now): default -2')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS. Default=now')
    parser.add_argument('--starttime', action='store', dest='starttime', default=None, type=str,
                        help='Desired starttime YYYY-mm-dd HH:MM:SS. Default=None')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default=None, type=str,
                        help='choose supported data source (case independant) eg NOAA or CONTRAILS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product eg river_water_level: Only required for Contrails')
    parser.add_argument('--station_list', action='store', dest='station_list', default=None, type=str,
                        help='Choose a non-default location/filename for a stationlist')
    parser.add_argument('--config_name', action='store', dest='config_name', default=None, type=str,
                        help='Choose a non-default contrails auth config_name')
    parser.add_argument('--ofile', action='store', dest='ofile', default=None, type=str,
                        help='Choose a non-default data product output directory')
    parser.add_argument('--ometafile', action='store', dest='ometafile', default=None, type=str,
                        help='Choose a non-default metadata output directory')
    args = parser.parse_args()
    sys.exit(main(args))
