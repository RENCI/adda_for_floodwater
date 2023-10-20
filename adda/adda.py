#!/usr/bin/env python
#
import os,sys
import glob
import numpy as np
import pandas as pd
import datetime as dt
print('\n'.join(sys.path))

import harvester.fetch_adcirc_data as fetch_adcirc_data
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.get_observations_stations as get_obs_stations
import processing.compute_error_field as compute_error_field
import processing.interpolate_scaled_offset_field as interpolate_scaled_offset_field
import gridmap.grid_to_station_maps as grid_to_station_maps
import adda_visualization_plots as adda_visualization_plots

from utilities.utilities import utilities as utilities
import utilities.io_utilities as io_utilities
from argparse import ArgumentParser

import joblib

def main(args):
    """
    """

    config_file=args.da_config_file
    config = utilities.init_logging(subdir=None, config_file=config_file)
    print(f'config={config}')
    map_file=config['mapfile']
    Ndays=int(-config['max_lookback_days'])
    rootdir=config['rundir'] 
    dwlc_filename=config['dwlc_filename'] 

    # Basic checks
    if args.gridname is None:
        utilities.log.error('Input gridname cannot be None')
        sys.exit(1)
    if args.fw_arch_dir is None:
        utilities.log.error('Floodwater Arch Dir cannot be None, empty, or "./".')
        sys.exit(1)

    fort63_style=args.fort63_style
    fw_arch_dir=args.fw_arch_dir

    # Set up IO env
    #utilities.log.info("Product Level Working in {}.".format(os.getcwd()))

    # Set up the times information No need to worry about values for hh:mm:ssZ Subsequent resampling cleans that up

#    if args.timeout is None: # Set to a default of now()
#        tnow = dt.datetime.now()
#        stoptime = tnow.strftime('%Y-%m-%d %H:%M:%S')
#    else:
#        stoptime=args.timeout
    tnow = dt.datetime.now()
    stoptime = tnow.strftime('%Y-%m-%d %H:%M:%S')

    #dt_starttime = dt.datetime.strptime(stoptime,'%Y-%m-%d %H:%M:%S')+dt.timedelta(days=args.ndays)
    total_hours = 0 if args.ndays==0 else abs(args.ndays*24) - 1 # Ensures the final list length INCLUDES the now time as a member and is a multiple of 24 hours

    if args.input_url is None:
        # generate list of urls
        utilities.log.info('Generating list of urls based on fw_arch_dir')
        utilities.log.info(f'fw_arch_dir={fw_arch_dir}')
        temp=fw_arch_dir+"/archive/"
        n=2*4+1
        urls=glob.glob(temp + "/*/*/adcirc/analysis/fort.63.nc", recursive = True)
        urls=urls[-(n+1):-1]
        #dt_starttime = dt.datetime.strptime(stoptime,'%Y-%m-%d %H:%M:%S')+np.sign(args.ndays)*dt.timedelta(hours=total_hours) # Should support lookback AND look forward
        #starttime = dt_starttime.strftime('%Y-%m-%d %H:%M:%S')
        #print('Total_hours, Starttime, Stoptime and ndays {}, {}. {}'.format(total_hours, starttime, stoptime,args.ndays))
    else:
        # assume "input_url" is a file of urls.
        if not os.path.exists(args.input_url):
            utilities.log.error(f'Input_url file {args.input_url} not found.') 
            sys.exit(1)
        with open(args.input_url) as f:
            urls = f.read().splitlines()

    utilities.log.info(f'URL list is {urls}') 
##
## get the adcirc station data
##

    #Note specifying the map_file REQUIRES the listed file to have a fullpathname
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)   # ,mapfile=map_file)
    file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)
    #print(f'station_file={station_file}')

    rpl = get_adcirc_stations.get_adcirc_stations(source='TDS', product=args.data_product,
                station_list_file=station_file, 
                knockout_dict=None, fort63_style=fort63_style )

    # Convert URLs to desired fort type
    #if fort63_style:
    #    urls=get_adcirc_stations.convert_urls_to_63style(urls)
    #else:
    #    urls=get_adcirc_stations.convert_urls_to_61style(urls)
    #utilities.log.info(f'Set of URLs to process {urls}')

    #print(f'fort63_style={fort63_style}')

    # Fetch best resolution and no resamplingA

    data_adc,meta_adc=rpl.fetch_station_product(urls, return_sample_min=args.return_sample_min, fort63_style=fort63_style  )

    # Revert Harvester filling of nans to -99999 back to nans
    data_adc.replace('-99999',np.nan,inplace=True)
    meta_adc.replace('-99999',np.nan,inplace=True)

    # Grab the grid coordinates for the url 
    #urls_63 = get_adcirc_stations.convert_urls_to_63style(urls)
    #print(urls_63)
    adc_coords = get_adcirc_stations.extract_adcirc_grid_coords( urls )

    #lons = adc_coords['LON']
    #lats = adc_coords['LAT']
    #adc_coords = {'lon':lons, 'lat':lats} # This created for backward compatibility

    #print(f'Grid name {rpl.gridname}')
    #print(f'Instance name {rpl.instance}')

    # Get a last piece of metadata for first url in iterable grabs either the time (%Y%m%s%H) or hurricane advisory (int)
    # Grab the adcirc time ranges for calling the observations code
    
    # Since ADCIRC start 6 hours early, using the ADCRC starttime can throw off the error averages at the end
    obs_starttime = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
    obs_endtime = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
    #obs_starttime = dt.datetime.strftime(starttime, '%Y-%m-%d %H:%M:%S')
    #obs_stoptime = dt.datetime.strftime(stoptime, '%Y-%m-%d %H:%M:%S')

    # Currently being used: Derived from ndays and timeout specification
    # If a special case input_url, need to grab times form the input url incase it was a hurricane
#    if args.input_url is None:
#        obs_starttime=starttime
#        obs_endtime=stoptime

    #

    # If equal fetch the min/max values of the adcirc data.
    if obs_starttime == obs_endtime:
            obs_starttime = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
            obs_endtime = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')

    # Now double check that the time orders properly
    (obs_starttime,obs_endtime)= (obs_endtime,obs_starttime) if obs_starttime > obs_endtime else (obs_starttime,obs_endtime)

    print(f'Time ranges {obs_starttime} and {obs_endtime}')

    # Cnstruct iometadata to update all filename - want diff format
    #io_start = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y%m%d%H%M')
    #io_end = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y%m%d%H%M')
    #iometadata = io_start+'_'+io_end 

#    if args.use_iometadata:
#        rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra='ADDA'+'_'+iometadata)
#    else:
#        rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra='')
#        iometadata='' 
    iometadata='' 

    # Write the data to disk
    iosubdir='adcpkl'
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)

    # Write selected in Pickle data 
    metapkl = io_utilities.write_pickle(meta_adc,rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_metadata')
    detailedpkl = io_utilities.write_pickle(data_adc, rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_detailed')

    # Grab the adcirc time ranges for calling the observations code
    ##obs_starttime = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
    ##obs_endtime = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')

    # Optionally Write selected in JSON format Convert and write selected JSON data
    #data_adc.index = data_adc.index.strftime('%Y-%m-%d %H:%M:%S')
    data_adc_4json = data_adc.copy()
    data_adc_4json.index = data_adc_4json.index.strftime('%Y-%m-%d %H:%M:%S')
    metajson = io_utilities.write_json(meta_adc,rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_metadata')
    detailedjson = io_utilities.write_json(data_adc_4json, rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_detailed')
   
    #print('Write out coordinates')
    #ADCfilecoords=utilities.write_csv(df_adc_coords, rootdir=rootdir,subdir=iosubdir,fileroot='adc_coord',iometadata=iometadata)
    #utilities.log.info('Wrote grid coords to {}'.format(ADCfilecoords))
    #print('Write out coordinates dict in JSON format')
    ADCfilecoordsJson = io_utilities.write_dict_to_json(adc_coords, rootdir=rootdir,subdir=iosubdir,fileroot='adc_coord')
    #ADCfilecoordsJson = rootdir+'/adc_coord'+iometadata+'.json'
    io_utilities.write_json_file(adc_coords, ADCfilecoordsJson)
    utilities.log.info('Wrote grid coords to {}'.format(ADCfilecoordsJson))

##
## set up and fetch the NOAA observations - using a time range specified by the available adcirc data set
##

    obs = get_obs_stations.get_obs_stations(source='NOAAWEB', product='water_level',
                contrails_yamlname='None',
                knockout_dict=None, station_list_file=station_file)
    # Get data at highest resolution
    data_obs,meta_obs=obs.fetch_station_product((obs_starttime,obs_endtime), return_sample_min=0)
    data_obs.replace('-99999',np.nan,inplace=True)
    meta_obs.replace('-99999',np.nan,inplace=True)
    # Remove stations with too many nans ( Note Harvester would have previously removed stations that are ALL NANS)
    data_thresholded = obs.remove_missingness_stations(data_obs, max_nan_percentage_cutoff=10)  # (Maximum allowable nans %)
    # meta_thresholded = meta_obs.loc[data_thresholded.columns.tolist()]
    meta = set(data_thresholded.columns.tolist()).intersection(meta_obs.index.to_list())
    meta_thresholded = meta_obs.loc[meta]
    # Apply a moving average (smooth) the data performed the required resampling to the desired rate followed by interpolating
    data_obs_smoothed = obs.fetch_smoothed_station_product(data_thresholded, return_sample_min=60, window=11)

    # Write the data to disk in a way that mimics ADDA
    iosubdir='obspkl'
    # Write selected in Pickle data 
    metapkl = io_utilities.write_pickle(meta_thresholded,rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_metadata')
    detailedpkl = io_utilities.write_pickle(data_thresholded, rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_detailed')
    smoothpkl = io_utilities.write_pickle(data_obs_smoothed, rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_smoothed')
    print('Finished OBS')

##
## compute errors
##
    comp_err = compute_error_field.compute_error_field(data_obs_smoothed, data_adc, meta_obs, n_pad=1) # All default params
    comp_err._intersection_stations()
    comp_err._intersection_times()
    comp_err._tidal_transform_data()
    #print(f'COMPERR input times {obs_starttime} and {obs_endtime}')
    comp_err._apply_time_bounds((obs_starttime,obs_endtime)) # redundant but here for illustration
    comp_err._compute_and_average_errors()

    # Set up IO env
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))

    iosubdir='errorfield'

    # Write Pickle files 
    tideTimeErrors = io_utilities.write_pickle(comp_err.diff,rootdir=rootdir,subdir=iosubdir,fileroot='tideTimeErrors')
    utilities.log.info('Wrote out Pickle {}'.format(tideTimeErrors))

    # Write selected in JSON format. Basically the Merged_dict data multi-index set
    data_dict = compute_error_field.combine_data_to_dict(comp_err.adc,comp_err.obs,comp_err.diff, product='WL')
    dict_json = io_utilities.write_dict_to_json(data_dict, rootdir=rootdir,subdir=iosubdir,fileroot='adc_obs_error_merged')
    utilities.log.info('Wrote out JSON {}'.format(dict_json))

    # Write the CSV files which carry averages
    station_summary_aves=io_utilities.write_csv(comp_err.df_final, rootdir=rootdir,subdir=iosubdir,fileroot='stationSummaryAves')
    station_period_aves=io_utilities.write_csv(comp_err.df_cycle_aves, rootdir=rootdir,subdir=iosubdir,fileroot='stationPeriodAves')
    #utilities.log.info('Wrote out CSV cyclic averages {}'.format(station_period_aves))
    utilities.log.info('Wrote out CSV summary averages {}'.format(station_summary_aves))
##
## Perform the interpolation. Must use the LinearRBF approach
##
    # map_file=os.path.join(os.path.dirname(__file__), './supporting_data', 'grid_to_stationfile_maps.yml')
    file_stations = station_summary_aves # Better to reread using "get_station_values"
    #file_land_controls = grid_to_station_maps.read_specified_filename_from_map(gridname=args.gridname, datatype='LANDCONTROL')   #  mapfile=map_file
    #file_water_controls = grid_to_station_maps.read_specified_filename_from_map(gridname=args.gridname, datatype='WATERCONTROL') #  mapfile=map_file
    #
    df_stations = interpolate_scaled_offset_field.get_station_values(fname=file_stations)
    df_land_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_land_controls)
    df_water_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_water_controls)

    set_of_dfs=list()
    set_of_clamp_dfs=list()

    set_of_dfs.append(df_stations) # Always need this
    if df_water_controls is not None:
        set_of_dfs.append(df_water_controls)
        set_of_clamp_dfs.append(df_water_controls)
    if df_land_controls is not None:
        set_of_dfs.append(df_land_controls)
        set_of_clamp_dfs.append(df_land_controls)
        utilities.log.info(f'Final Land control results {df_land_controls}')

    df_ClampControl = pd.concat(set_of_clamp_dfs,axis=0)
    # This can only be used with LinearRBF. More testing is required
    pathpoly = interpolate_scaled_offset_field.build_polygon_path(df_ClampControl)

    utilities.log.info('construct_interpolation_model: Number of dfs to combine for interpolation is {}'.format(len(set_of_dfs)))

    df_combined=interpolate_scaled_offset_field.combine_datasets_for_interpolation(set_of_dfs)
    print('df_combined')
    print(df_combined)
    model = interpolate_scaled_offset_field.interpolation_model_fit(df_combined, interpolation_type='LinearRBF')

    # Build new grid
    df_extrapolated_ADCIRC_GRID = interpolate_scaled_offset_field.interpolation_model_transform(adc_coords, model=model, input_grid_type='points',pathpoly=pathpoly)

    # do A TEST FIT
#    if args.cv_testing:
#        print('TEST FIT')
#        full_scores, best_scores = interpolate_scaled_offset_field.test_interpolation_fit(df_source=df_stations, df_land_controls=df_land_controls.copy(), df_water_controls=df_water_controls, cv_splits=5, nearest_neighbors=3)
#        print('ADDA')
#        print(best_scores)
#        print(full_scores)

##
## Write out datafiles 
##
    gridfile = io_utilities.write_ADCIRC_formatted_gridfield_to_Disk(df_extrapolated_ADCIRC_GRID, value_name='VAL', rootdir='../',subdir='',filename=dwlc_filename)
    utilities.log.info('Wrote ADCIRC offset field to {}'.format(gridfile))
    adcirc_extrapolated_pkl = io_utilities.write_pickle(df_extrapolated_ADCIRC_GRID, rootdir=rootdir,subdir=iosubdir,fileroot='interpolated_wl')
    utilities.log.info('Wrote ADCIRC offset field PKL {}'.format(adcirc_extrapolated_pkl))
    #print('Finished')

    # Test using the generic grid and plot to see the generated offset surface. Use the same model as previously generated
    adc_plot_grid = interpolate_scaled_offset_field.generic_grid()
    df_plot_transformed = interpolate_scaled_offset_field.interpolation_model_transform(adc_plot_grid, model=model, input_grid_type='grid',pathpoly=pathpoly) 

    # Write out the model for posterity
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(rootdir, iosubdir, 'interpolate_linear_model.h5')
    try:
        joblib.dump(model, newfilename)
        status = True
        utilities.log.info('Saved model file '+str(newfilename))
    except:
        utilities.log.error('Could not dump model file to disk '+ newfilename)

##
## Apply the model to a 400x500 grid and plot the surface, stations, clamps
##
    iosubdir=''
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(rootdir, iosubdir, 'surface.png')
    adda_visualization_plots.save_plot_model( adc_plot_grid=adc_plot_grid, df_surface=df_plot_transformed, df_stations=df_stations, df_land_control=df_land_controls, df_water_control=df_water_controls, filename=newfilename, plot_now=False)
    utilities.log.info('Saved IMAGE file to {}'.format(newfilename))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product eg water_level')
    parser.add_argument('--return_sample_min', action='store', dest='return_sample_min', default=60, type=int,
                        help='return_sample_min is the time stepping in the final data objects. (mins)')
    parser.add_argument('--ndays', default=-2, action='store', dest='ndays',help='Day lag (usually < 0)', type=int)
    #parser.add_argument('--timeout', default=None, action='store', dest='timeout', help='YYYY-mm-dd HH:MM:SS. Latest day of analysis def to now()', type=str)
    parser.add_argument('--da_config_file', action='store', dest='da_config_file', default=None,
                        help='String: yml config for DA info')
    parser.add_argument('--fort63_style', action='store_true',default=True, 
                        help='Boolean: Will inform Harvester to use fort.63.methods to get station nodesids')
    parser.add_argument('--gridname', action='store',dest='gridname', default=None,
                        help='str: Select appropriate gridname Default is hsofs')
    parser.add_argument('--map_file', action='store',dest='map_file', default=None,
                        help='str: Select appropriate map_file yml; for grid lookup')
    parser.add_argument('--cv_testing', action='store_true', dest='cv_testing', default=False,
                        help='Boolean: Invoke a CV procedure for post model CV testing')
    #parser.add_argument('--use_iometadata', action='store_true', dest='use_iometadata', default=False,
                        #help='Boolean: Include the iometadata time range to all output files and dirs')
    parser.add_argument('--input_url', action='store', dest='input_url', default=None, type=str,
                        help='TDS url to fetch ADCIRC data')
    parser.add_argument('--fw_arch_dir', action='store', dest='fw_arch_dir', default=None, type=str,
                        help='Location of FloodWater archive dir for the suite.')
    args = parser.parse_args()
    sys.exit(main(args))
