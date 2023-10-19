#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

#
# A file suitable for use by ADDA,APSVIZ2,Reanalysis to fetch ADCIRC water levels from the TDS
# The TDS input to this class is a list of URLs. If you require running this by specifying TIMES, 
# then you must preprocess the data into a list of URLs.
#
# TODO Check into the case where ADCIRC returns completely empty stations. This filtering may have been 
# turned off n the Harvester codes.

import os,sys
import pandas as pd
import datetime as dt
import math

from harvester.fetch_station_data import adcirc_fetch_data
from utilities.utilities import utilities as utilities

import time as tm

# Currently supported sources
SOURCES = ['TDS','FILE']

def get_adcirc_stations_fort63_style(fname=None)->pd.DataFrame:
    """
    This style of data input is required when using the fort63_style methods. Simply read a list of stations from a csv file.
    This gets read into a DataFrame. File MUST contain at least Node and stationid columns

    Parameters:
        fname: <str> full path to a valid stationid file
    Returns:
        DataFrame: [stationid, Node]

    """
    if fname is None:
        utilities.log.error(f'No Fort63_style ADCIRC station file assigned {fname}: Abort with status 0')
        sys.exit(0)
    df = pd.read_csv(fname, index_col=0, header=0, skiprows=[1], sep=',',dtype=str)
    try:
        df["Node"]=df["Node"].astype(int)
    except IndexError as e:
        utilities.log.error(f'Unsuccessful fort63_station read. Perhaps no Node column ? {e}')
    return df

def get_adcirc_stations_fort61_style(fname=None):
    """
    Simply read a list of stations from a csv file. File MUST contain a stationid column
    Generally, we simply combine NOAA and Contrails into a single list. It is okay to include stations not likely to exist since
    the processing stage will simply remove them

    Parameters:
        fname: <str> full path to a valid stationid file
    Returns:
        DataFrame: stationid

    """
    if fname is None:
        utilities.log.error(f'No Fort61_style ADCIRC station file assigned {fname}: Abort with status 0')
        sys.exit(0)
    df = pd.read_csv(fname, index_col=0, header=0, skiprows=[1],dtype=str)
    adcirc_stations_list = df["stationid"].to_list()
    adcirc_stations=[word.rstrip() for word in adcirc_stations_list]
    return adcirc_stations

def format_data_frames(df):
    """
    A Common formatting used by all sources
    """
    df.index = df.index.strftime('%Y-%m-%dT%H:%M:%S')
    df.reset_index(inplace=True)
    df_out=pd.melt(df, id_vars=['TIME'])
    df_out.columns=('TIME','STATION',PRODUCT.upper())
    df_out.set_index('TIME',inplace=True)
    return df_out

# The hurricane methods are for the future
def check_advisory(value, dformat='%Y%m%d%H'):
    """
    Try to determine if an advisory number was passed instead of a time value

    Parameters:
        value: <str> The time/adv word extracted form a url
        dformat: <str> format of the input str time 
    Returns:
        state_hurricane: <bool> True if hurricane
    """
    state_hurricane=False
    utilities.log.debug(f'Check advisory {value}')
    try:
        test=dt.datetime.strptime(value,dformat) # '%Y%m%d%H')
        utilities.log.info(f'A timestamp data was found: Not a Hurricane advisory ? {test}')
    except ValueError:
        try:
            outid = int(value)
            state_hurricane=True
        except ValueError:
            utilities.log.error(f'Expected an Advisory value but could not convert to int {value}')
            sys.exit(1)
    utilities.log.info(f'URL state_hurricane is {state_hurricane}')
    return state_hurricane

def check_if_hurricane(urls):
    """
    Very simple procedure but requires having the TDS nomenclature. This will not work 
    for generic locally stored files
    Only need to check one valid url from the list. This presumnes they all have the same
    grid, instance, class, etc

    Parameters:
        urls: list(str) A list of valid urls
    Returns:
        state_hurricane: <bool>  True if hurricane
    """
    if not isinstance(urls, list):
        utilities.log.error('time: URLs must be in list form')
        sys.exit(1)
    for url in urls:
        try:
            words=url.split('/')
            state_hurricane = check_advisory(words[-6])
            break
        except IndexError as e:
            utilities.log.error(f'check_if_hurricane Uexpected failure try next:{e}')
    return state_hurricane

# This is retained for backward compatability
# and is only used by this method's main
def convert_input_url_to_nowcast(urls):
    """
    Though one could call this method using a nowcast url, occasionally we want to be able to
    only pass a forecast type url and, from that, figure out what the corresponding nowcast url might be.
    This assume a proper TDS formatted url and makes no attempts to validate the 
    the constructed url. Either it exists or this methiod exits(1)

    To use this feature:
    We mandate that the url is used to access TDS data. The "ensemble" information must be in position .split('/')[-2]

    Parameters:
        urls: list(str) A list of valid urls
    Returns:
        urls: list(str) with all ensmble values set to "nowcast"

    """
    if not isinstance(urls, list):
        utilities.log.error('nowcast: URLs must be in list form')
        sys.exit(1)
    newurls=list()
    for url in urls:
        urlwords=url.split('/')
        urlwords[-2]='nowcast'
        newurls.append('/'.join(urlwords))
    utilities.log.info('Modified input URL to be a nowcast type')
    return newurls

# This is also retained for backward compatability
# and is only used by this method's main
def combine_metadata_with_station_tuples(df_meta, station_nodes, fort63_style):
    """
    Grabs the list of station tuples [(stationid,adcirc node)] and adds a 
    column to the df_meta object

    Parameters:
        df_meta: Input metadata dataframe stations x metadata
        station_nodes: list(str) list of stationids
        fort63_style: <bool> Value modifies the associated station Node value for prepending 63_ or 61_
        
    Returns:
        df_meta: An improved ddata ataframe of entries stations x metadata
    """
    df_stations=pd.DataFrame(station_nodes, columns=['STATION','Node'])
    df_stations.set_index('STATION',inplace=True)
    df_stations['Node']=df_stations['Node'].astype('str')
    typ = '63_' if fort63_style else '61_'
    df_stations['Node']=typ+df_stations['Node'] # Create something like 61-12345
    df_meta = pd.concat([df_meta,df_stations],axis=1)
    return df_meta

##
## Globals
##
dformat='%Y-%m-%d %H:%M:%S'
GLOBAL_TIMEZONE='gmt' # Every source is set or presumed to return times in the zone
PRODUCT='water_level'

##
## Run stations. These functions generaly get called by outside callers
##

def process_adcirc_stations(urls, adcirc_stations, gridname, ensemble, sitename, data_product='water_level', resample_mins=0, fort63_style=False, variable_name='zeta', keep_earliest_url=True):
    """
    Helper function to take an input list of times, stations, and product and return a data set and associated metadata set

    Parameters:
        urls: list<str>. Previously generated list of urls that span the desired time range 
        stations: list(str). List of desired stations
        gridname: <str> gridname of the urls (opt)
        ensemble: <str> ensemble (opt)
        sitename <str>: sitename (opt)
        data_product: <str> (def data_product). An AST named data product (Not the True data product name) 
        resample_mins: <int> Returned time series with a sampling of resample_mins
        fort63_style: <bool> Request fetching water levels from fort.63.nc/swan_HS.63.nc. Requires compatible station dataframe 
        variable_nme: <str>  If wanting swan results,must pass in the correct NC4 data varable name of 'swan_HS'
    Returns:
        df_adcirc_data: DataFrame (time x station)
        df_adcirc_meta: DataFrame (station x metadata)
    """

    # Fetch the data
    try:
        if data_product != 'water_level':
            utilities.log.error('ADCIRC data product can only be: water_level')
            sys.exit(1)
        adcirc = adcirc_fetch_data(adcirc_stations, urls, data_product, sitename=sitename, gridname=gridname, castType=ensemble.rstrip(), resample_mins=resample_mins, fort63_style=fort63_style, variable_name=variable_name, keep_earliest_url=keep_earliest_url)
        df_adcirc_data = adcirc.aggregate_station_data()
        df_adcirc_meta = adcirc.aggregate_station_metadata()
        df_adcirc_meta.index.name='STATION'
        # Add insitu found nodes for metadata
        station_nodes = adcirc.available_stations_tuple
        df_adcirc_meta=combine_metadata_with_station_tuples(df_adcirc_meta, station_nodes, fort63_style)
    except Exception as e:
        utilities.log.error(f'Failed fetching ADCIRC: {e}')
        raise
    return df_adcirc_data, df_adcirc_meta 

def first_true(iterable, default=False, pred=None):
    """
    itertools recipe found in the Python 3 docs
    Returns the first true value in the iterable.
    If no true value is found, returns *default*
    If *pred* is not None, returns the first item
    for which pred(item) is true.

    first_true([a,b,c], x) --> a or b or c or x
    first_true([a,b], x, f) --> a if f(a) else b if f(b) else x

    """
    return next(filter(pred, iterable), default)

def strip_time_from_url(urls)->str:
    """
    We mandate that the URLs input to this fetcher are those used to access the TDS data. The "time" information will be in position .split('/')[-6]
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
         time: <str>: in either TDS formatted '%Y%m%d%H' or possibly as a hurricane advisory string (to be checked later)
    """
    url = grab_first_url_from_urllist(urls)
    try:
        words = url.split('/')
        ttime=words[-6] # Always count from the back. NOTE if a hurrican this could be an advisory number.
    except IndexError as e:
        utilities.log.error(f'strip_time_from_url Uexpected failure try next:{e}')
    return ttime

def strip_ensemble_from_url(urls)->str:
    """
    We mandate that the URLs input to this fetcher are those used to access the TDS data. The "ensemble" information will be in position .split('/')[-2]
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
        Ensemble: <str>
    """
    url = grab_first_url_from_urllist(urls)
    try:
        words = url.split('/')
        ensemble=words[-2] # Usually nowcast,forecast, etc 
    except IndexError as e:
        utilities.log.error(f'strip_ensemble_from_url Uexpected failure try next:{e}')
    return ensemble

def strip_instance_from_url(urls)->str:
    """
    We mandate that the URLs input to this fetcher are those used to access the TDS data. The "instance" information will be in position .split('/')[-3]
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
        Instance: <str>
    """
    url = grab_first_url_from_urllist(urls)

    try:
        words = url.split('/')
        instance=words[-3] 
    except IndexError as e:
        utilities.log.error(f'strip_instance_from_url Uexpected failure try next:{e}')
    return instance 

def strip_storm_number_from_url(urls, fill='NONE')->str:
    """
    We mandate that the Hurricane URLs input to this fetcher are those used to access the TDS/ECFlow data. The "storm number" (eg al09) 
    information will be in position .split('/')[-7]
    eg. 'http://tds.renci.org/thredds/dodsC/2021/al09/11/hsofs/hatteras.renci.org/hsofs-al09-bob/nhcOfcl/fort.61.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
        Storm name: <str>
    """
    url = grab_first_url_from_urllist(urls)
    try:
        words = url.split('/')
        stormnumber=words[-7]
    except IndexError as e:
        utilities.log.error(f'strip_storm_number_from_url Uexpected failure try next:{e}')
        stormnumber=fill
    return stormnumber

def strip_sitename_from_url(urls, fill='NoSite')->str:
    """
    Here we attempt to find which site the url was computed. We read the
    machine name from the url and lookup in the local dict for the predefined canonical site
    name. If no such name use the fill value
    The "machine name" information will be in position .split('/')[-4]. It may consist of multiple words.
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
        fill: (str) manually specify the value for the sitename metadata
    Returns:
        canonical site name: (str) eg RENCI,PSC
    """
    known_sites= {'hatteras.renci.org':'RENCI', 
                  'ht-ncfs.renci.org': 'RENCI',
                  'sapelo2': 'UGA',
                  'bridges2.psc.edu': 'PSC'}

    url = grab_first_url_from_urllist(urls)
    try:
        words = url.split('/')
        machine=words[-4] 
    except IndexError as e:
        utilities.log.error(f'strip_sitename_from_url Uexpected failure try next:{e}')
    site = known_sites[machine] if machine in known_sites.keys() else fill
    return site

def grab_gridname_from_url(urls)->str:
    """
    We mandate that the URLs input to this fetcher are those used to access the TDS data. The "grid" information will be in position .split('/')[-2]
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
        grid.upper(): <str>
    """
    url = grab_first_url_from_urllist(urls)
    try:
        words = url.split('/')
        grid=words[-5] # Usually nowcast,forecast, etc 
    except IndexError as e:
        utilities.log.error(f'strip_gridname_from_url Uexpected failure try next:{e}')
    return grid.upper()

def grab_first_url_from_urllist(urls)->str:
    """
    eg. 'http://tds.renci.org/thredds/dodsC/2021/nam/2021052318/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.63.nc'
    
    Parameters:
        urls: list(str). list of valid urls
    Returns:
        url: <str> . Fetch first available, valid url in the list
    """
    if not isinstance(urls, list):
        utilities.log.error('first url: URLs must be in list form')
        sys.exit(1)
    url = first_true(urls)
    return url

def main(args):
    """
    We require the provided URL is using the typical TDS nomenclature and that the timestamp is in ('/') position -6
    Moreover, This time stamp behaves a little differently if fetching a nowcast versus a forecast (pre vs post). For now, we will
    annotate final .csv files with _TIME_ corresponding to the actual url data starttime.
    """

    #main_config = utilities.init_logging(subdir=None, config_file='../config/main.yml')
    main_name = args.main_config if args.main_config is not None else os.path.join(os.path.dirname(__file__),'../config','main.yml')
    main_config = utilities.init_logging(subdir=None,config_file=main_name)

    if args.sources:
         print('Return list of sources')
         return SOURCES
         sys.exit(0)

    data_source = args.data_source

    if data_source.upper() in SOURCES:
        utilities.log.info(f'Found selected data source {data_source}')
    else:
        utilities.log.error(f'Invalid data source {data_source}')
        sys.exit(1)

    urls = args.urls
    if urls==None:
        utilities.log.error('No URL was specified: Abort')
        sys.exit(1)

    if not isinstance(urls, list):
        utilities.log.error('urls: URLs must be in list form: Converting')
        urls = [urls]

    if args.convertToNowcast:
        utilities.log.info('Requested conversion to Nowcast')
        urls = convert_input_url_to_nowcast(urls)

    data_product = args.data_product
    if data_product != 'water_level':
        utilities.log.error('ADCIRC: The Only available data product is water_level: {data_product}')
        sys.exit(1)
    else:
        utilities.log.info(f'Chosen data source {data_source}')

    variable_name = args.variable_name
    utilities.log.info(f'Selected variable name is {variable_name}')

    # If reading in a LOCAL url then we must skip the following steps
    if not args.raw_local_url:
        # Check if this is a Hurricane
        if not check_if_hurricane(urls):
            utilities.log.info('URL is not a Hurricane advisory')
            #sys.exit(1)
            urltimeStr = strip_time_from_url(urls)
            urltime = dt.datetime.strptime(urltimeStr,'%Y%m%d%H')
            runtime=dt.datetime.strftime(urltime, dformat)
        else:
            utilities.log.info('URL is a Hurricane')
            urladvisory = strip_time_from_url(urls)
            runtime=urladvisory
        ensemble = strip_ensemble_from_url(urls)  # Only need to check on of them
        gridname = grab_gridname_from_url(urls)   # ditto
        sitename = strip_sitename_from_url(urls)
        #starttime='2021-12-08 12:00:00'
        utilities.log.info(f'Selected run time/Advisory range is {runtime}')
    else:
        print('Raw url data structure specified')
        utilities.log.info('Running a RAW url type')
        ensemble = 'NONE'
        gridname = 'NONE'
        sitename = 'NONE'

    if args.fort63_style:
        utilities.log.info('Fort_63 style station inputs specified')

    ##
    ## Start the processing
    ##

    t0 = tm.time()
    if data_source.upper()=='TDS':
        excludedStations=list()
        # Use default station list
        fname_stations = args.station_list if args.station_list is not None else '../supporting_data/CERA_NOAA_HSOFS_stations_V3.1.csv'
        try:
            if args.fort63_style:
                adcirc_stations=get_adcirc_stations_fort63_style(fname=fname_stations)
            else:
                adcirc_stations=get_adcirc_stations_fort61_style(fname=fname_stations)
            # Need to specify something for the metadata 
            if not args.raw_local_url:
               adcirc_metadata='_'+ensemble+'_'+gridname.upper()+'_'+runtime.replace(' ','T')
               adcirc_metadata=f"_{ensemble}_{gridname.upper()}_{runtime.replace(' ','T')}"
            else:
                adcirc_metadata='Raw_data'

            data, meta = process_adcirc_stations(urls, adcirc_stations, gridname, ensemble, sitename, data_product, resample_mins=0, fort63_style=args.fort63_style, variable_name=variable_name)
        except Exception as e:
            utilities.log.info('Failed to process the ADCIRC data: {e}')
            raise            
        ## df_adcirc_data = format_data_frames(data)

        ## Skip the melt: Note this is only for the written file
        df_adcirc_data = data
        # Output 
        try:
            if args.ofile is not None:
                dataf=f'%s/adcirc_stationdata%s.csv'% (args.ofile,adcirc_metadata)
                metaf=f'%s/adcirc_stationdata_meta%s.csv'% (args.ometafile,adcirc_metadata)
            else:
                dataf=f'./adcirc_stationdata%s.csv'%adcirc_metadata
                metaf=f'./adcirc_stationdata_meta%s.csv'%adcirc_metadata
            df_adcirc_data.to_csv(dataf)
            meta.to_csv(metaf)
            utilities.log.info(f'ADCIRC data has been stored {dataf},{metaf}')
        except Exception as e:
            utilities.log.error(f'Error: ADCIRC: Failed Write {e}')
            raise
            #sys.exit(1)
    total_time = tm.time() - t0
    utilities.log.info(f'Finished with data source {data_source}. Time was {total_time}')
    utilities.log.info('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='TDS', type=str,
                        help='choose supported data source: default = TDS')
    parser.add_argument('--urls', nargs='+', action='store', dest='urls', default=None, type=str,
                        help='TDS url to fetcb ADCIRC data')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product: default is water_level')
    parser.add_argument('--convertToNowcast', action='store_true',
                        help='Attempts to force input URL into a nowcast url assuming normal TDS conventions')
    parser.add_argument('--fort63_style', action='store_true', 
                        help='Boolean: Will inform Harvester to use fort.63.methods to get station nodesids')
    parser.add_argument('--station_list', action='store', dest='station_list', default=None, type=str,
                        help='Choose a non-default location/filename for a stationlist')
    parser.add_argument('--ofile', action='store', dest='ofile', default=None, type=str,
                        help='Choose a non-default data product output directory')
    parser.add_argument('--ometafile', action='store', dest='ometafile', default=None, type=str,
                        help='Choose a non-default metadata output directory')
    parser.add_argument('--raw_local_url', action='store_true',
                        help='Specify input url is locally stored')
    parser.add_argument('--variable_name', action='store', dest='variable_name', default='zeta', type=str,
                        help='Choose a non-default netCDF4 variable name')
    parser.add_argument('--main_config', action='store', dest='main_config', default=None, type=str,
                        help='Choose a non-default main.yml')
    args = parser.parse_args()
    sys.exit(main(args))
