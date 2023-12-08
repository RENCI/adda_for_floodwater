#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

import sys, os
import pandas as pd
import numpy as np
from utilities.utilities import utilities
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import mean_squared_error
from scipy import interpolate as sci
import matplotlib.path as mpltPath

############################################################################
# Helper functions

def get_lon_lat_control_values(fname=None, header_data='val')->pd.DataFrame:
    """
    Simply read the CSV file 
    Ensure the columns are all have the names: LON,LAT,VALUE

    land/water control data are always lon,lat,val. 

    Parameters:
        fname: (str) Filename of the desired land/water control points
        header_data: (str) Alternative name for the clamp/control node data header if not the default
            Will be changed to VAL on exit

    Returns:
        df: DataFrame with the headers: LON,LAT,VAL 

    """
    if fname is None:
        utilities.log.error('No lon_lat_val data file assigned: Abort interpolation')
        sys.exit(1)
    df = pd.read_csv(fname, header=0)
    df=df.rename(columns = {header_data:'val'})
    df.columns = df.columns.str.replace(' ', '') 
    df.columns = [x.upper() for x in df.columns]
    # Check columns: MUST include LON,LAT,VAL
    col_headers=list()
    for item in df.columns.to_list(): # Many times weve seen leading spaces in thes 
        if item=='LON' or item=='LAT' or item=='VAL': 
            col_headers.append(item)
    if len(col_headers)!=3:
        utilities.log.info('get_lon_lat_control_values: control data file {} doesnt have all the needed headers. {}: Abort'.format(fname,df.columns)) 
        sys.exit(1)
    return df[['LON','LAT','VAL']]

def get_station_values(fname=None, header_data='mean')->pd.DataFrame:
    """
    Simply read a proper CSV file 

    station (averaging) files can be more complicated than control point files
    but we assume stationids are on the index and that each has lon,lat,mean (or header_data)

    The INPUT headers are very specific to AST usage
    examples of typical header_data could be 'mean','std'

    Parameters:
        fname: (str) Filename of the desired station (averaging) file
        header_data: (str) Alternative name for the input desired values. Eg if not mean it might be "Period_-0"
            which by default is the mean of the most last period
            Will be changed to VAL on exit

    Returns:
        df: In the format STATIONID(index),LON,LAT,VAL 
    """
    if fname is None:
        utilities.log.error('No station (averaging) data file assigned: Abort interpolation')
        sys.exit(1)
    df = pd.read_csv(fname, header=0)
    df=df.rename(columns = {header_data:'val'})
    df['STATION'] = df['STATION'].astype(str)   # Need to be sure
    df.columns = df.columns.str.replace(' ', '')
    df.columns = [x.upper() for x in df.columns]
    df=df.rename(columns = {'STATION':'STATIONID'}) # Just a catch-all. No harm, no foul
    col_headers=list()
    for item in df.columns.to_list(): # Many times weve seen leading spaces in thes 
        if item=='STATIONID' or item=='LON' or item=='LAT' or item=='VAL':
            col_headers.append(item)
    if len(col_headers)!=4:
        utilities.log.info('get_station_values: station data file {} doesnt have all the needed headers. {}: Abort'.format(fname,df.columns))
        sys.exit(1)
    df=df[['STATIONID','LON','LAT','VAL']]
    return df.set_index('STATIONID')

def randomly_select_dataframe_rows(df, frac=1)->pd.DataFrame:
    """
    Intended for, but not restricted to, station (averaging) data
    This is primarily used if the user wants to perform statistical testing on the
    resulting interpolated fields. Setting frac==1 (default) returns the entire data set 
    with rows shuffled (potentially removing some order biases)
    Setting to < 1, returns a subset of the data randomly selected.

    This method is a placeholder for more sophisticated splitting/subsampling
    
    Parameters:
        df: DataFrame to reshuffle/resample
        frac: (int) fraction of randomly selected rows to retain

    Returns:
        dataFrame
    """
    return df.sample(frac=frac).reset_index(drop=True)

def knn_fit_control_points(df_source, df_controls, nearest_neighbors=3)->pd.DataFrame:
    """
    Fit the values of input control_points based on the KNN(k) of the input df_source DataFrame
    indata errors. Then apply that model to predict the control point values

    The fitting of control points using KNN has been generally suspended but still usable if needed.
  . Now land_controls are generally expected to simply be clamps(=0)

    This is setup to use a kd_tree method with uniform weights. The default metric is minkowski
         with a default of using manhattan_distance (l1),

    The input df_source data are double checked for missingness which must be removed

    Parameters:
        df_controls: DataFrame with the headers ('LON','LAT','VAL')
            Usually read by get_lon_lat_control_values()

        df_source: DataFrame with (at least) the headers ('LON','LAT','VAL')
            Usually read by get_station_values()
        nearest_neighbors: (int)  Number of neighbors to perform regression fit

    Returns:
        df_fit_controls: DatraFrame with the columns ('LON','LAT','VAL')

    """
    # Must remove entries with nans as VAL
    before_missing = df_source.shape[0]
    df_drop = df_source.dropna()
    after_missing = df_drop.shape[0]
    utilities.log.info('knn_fit_control_points: Input source rows {}, after missingness check rows {}'.format(before_missing, after_missing))

    weight_value = 'uniform' # '=distance' #='uniform'
    # Grab the Lon,Lat from the source (station) data
    source_coords = df_drop[['LON','LAT']].to_numpy()
    source_values =  df_drop[['VAL']].to_numpy()
    control_coords = df_controls[['LON','LAT']].to_numpy()
    # Build the KNN model
    neigh = KNeighborsRegressor(n_neighbors=nearest_neighbors, algorithm='kd_tree', weights=weight_value)
    neigh.fit(source_coords, source_values)
    # predict on the control points
    control_values = neigh.predict(control_coords)
    df_controls['VAL']=control_values
    #
    utilities.log.info('KNN estimated control points using {} and a nearest_neighbors of {} Total source points were {}'.format(weight_value, nearest_neighbors, len(source_values)))
    utilities.log.info('Mean value of predicted control point values is {}'.format(np.mean(control_values)))
    utilities.log.info('Standard Deviation of predicted control point values is {}'.format(np.std(control_values)))
    return df_controls

def combine_datasets_for_interpolation(list_dfs) ->pd.DataFrame:
    """
    This method simply combines the provided list of DataFrames into a single DataFrame for passage to
    the interpolation method. We use a list here, since, the caller may vary what data actually gets included
    for interpolation. The list is unbound. Each dataframe MUST contain the headers:
    ('LON','LAT','VAL')
    
    Parameters:
        list_dfs: list(dfs)

    Returns:
        df_combined: DataFrame containing only ('LON','LAT','VAL')
    """ 
    utilities.log.debug('combine_datasets_for_interpolation: Combining {} dataframes'.format(len(list_dfs)))
    filtered_dataframes = list()
    template = "An exception of type {0} occurred."
    try:
        for df in list_dfs:
            filtered_dataframes.append(df[['LON','LAT','VAL']])
    except Exception as ex:
        message = template.format(type(ex).__name__, ex.args)
        utilities.log.warn('Failed df concat. Perhaps non-conforming header names?: msg {}'.format(message))
        sys.exit(1)
    df_combined=pd.concat(filtered_dataframes, ignore_index=True)
    return df_combined

def build_polygon_path(df_clampcontrol):
    """
    For the input dataset construct a polygon object and generate a path. This is used as a boundary
    contour for the final interpolaton field. This EXPECTS the clamps positions to be ORDERED such that
    the polygon generation behaves well. This constraint is required to be performed by the calling program
    If the inteprolation procedure is not RBF, then the polygon method can be skipped

    Parameters:
        df_clampcontrol: A data frame the contains (at least) 'LON','LAT' columns. These should
            only include rows of clamps and/or control points. Not station data. Ideally, these points
            have been properly ORDERED by the user to generate a clean,singly clustered  polygon 
    """
    path = mpltPath.Path(df_clampcontrol[['LON','LAT']].to_numpy())
    return path

def interpolation_model_fit(df_combined, fill_value=0.0, interpolation_type='LinearRBF'):
    """
    Method to manage building interpolations of gauge (ADCIRC-Observational) errors (residuals)
    The model is constructed by interpolation of the input error points
    usually constrained with both zero-nodes positioned offcoast and land control points to better control
    fluctuations at the coast. 
    The types LINEATNDINTERPOLATOR, NEARESTNDINTERPOLATOR,  and CLOUGHTOCHER2DINTERPOLATOR are both based on a QHull generation.

    The types of interpolation models that are supported are:
        LinearNDInterpolator (*)
        CloughTocher2DInterpolator (*)
        NearestNDInterpolator (*)
        LinearRBF (default)
        (*) deprecated

    Parameters:
        df_combined (DatraFrame) is the DataFrame with columns: LON,LAT,VAL that contains all 
        source and clamp points. (DataFrame) It is expected the caller, if req, has already KNN processed the 
        land_control points and concatenated the water_control points

    Returns:
        model: The actual model use for subsequent extrapolations
    """
    model_choices=['LINEARNDINTERPOLATOR','CLOUGHTOCHER2DINTERPOLATOR','NEARESTNDINTERPOLATOR','LINEARRBF'] # version 2, March 2022 
    if interpolation_type.upper() in model_choices:
        interpolation_type = interpolation_type.upper()
    else:
        utilities.log.error('No valid interpolation_type specified. Got {}, Avail are {}'.format(interpolation_type,model_choices))
        sys.exit(1)

    # Must remove entries with nans as VAL
    before_missing = df_combined.shape[0]
    df_drop = df_combined.dropna() # Will drop a row for ANY value. TO ensure inut df_sources doesn;t have more than 'LON','LAT','VAL' columns

    after_missing = df_drop.shape[0]
    utilities.log.info('interpolation_model_fit: Input rows {}, after missingness check rows {}'.format(before_missing, after_missing))
   
    X = df_drop['LON'].to_list()
    Y = df_drop['LAT'].to_list()
    V = df_drop['VAL'].to_list()
    template = "An exception of type {0} occurred."
    try:
        if interpolation_type=='LINEARRBF':
            model = sci.Rbf(X,Y,V,function='linear', smooth=0)
        elif interpolation_type=='LINEARNDINTERPOLATOR':
            model = sci.LinearNDInterpolator((X,Y), V, fill_value=fill_value)
        elif interpolation_type=='NEARESTNDINTERPOLATOR':
            model = sci.NearestNDInterpolator((X,Y), V)
        elif interpolation_type=='CLOUGHTOCHER2DINTERPOLATOR':
            model = sci.CloughTocher2DInterpolator( (X,Y),V,fill_value=fill_value)
        else:
            utilities.log.info(f'Incorrect interpolator specified {interpolation_type}')
            sys.exit(1)
    except Exception as ex:
        message = template.format(type(ex).__name__, ex.args)
        utilities.log.warn('Interpolation failed: msg {}'.format(message))
        sys.exit(1)
    utilities.log.info('Finished building interpolation model')
    return model
    
# TODO Verify behavior if a point has a VAL=nan
def interpolation_model_transform(adc_combined, model=None, input_grid_type='points', pathpoly=None)-> pd.DataFrame: 
    """
    Apply the precomputed model to the LON/LAT input points. These inputs usually
    reflect the final grid such as an ADCIRC grid. The inputs can be of types points or grid.
    If the LON/LAT lists are meant to be actual points to extrapolate on, then points. Else, if
    that are INDEXES in the LON/LAT directions, then choose GRID.

    Parameters:
        adc_combined: dict with keys 'LON':list,'LAT':list
        input_grid_type: (str) either points or grid.
        model: Resuting from a prior 'fit'
        pathpoly: An mplpath path object that, if provided, is used to set the surface to zero outside the polygon 
            If none is provided, then this extra boundary step is skipped

    Returns:
        DataFrame: inteprolated grid
    """
    if model is None:
        utilities.log.error('interpolation_model_transform: No model was supplied: Abort')
        sys.exit(1)
    gridx = adc_combined['LON']
    gridy = adc_combined['LAT']
    d = []
    xl = len(gridx)
    yl = len(gridy)
    if input_grid_type=='grid':
        if pathpoly is not None:
            (x, y) = np.meshgrid(gridx, gridy)
            xxyy=np.array([x.flatten(), y.flatten()]).T
            notinpolyg = ~ pathpoly.contains_points(xxyy).reshape(x.shape)
        for y in range(0,yl):
            gy = gridy[y]
            for x in range(0,xl):
                gx=gridx[x]
                zval=model(gx,gy).item(0) # Assume only a single item
                d.append([gx, gy, zval]) 
        df = pd.DataFrame(d,columns=['LON','LAT','VAL'])
        if pathpoly is not None: df.loc[notinpolyg.flatten(),'VAL']=0
    else:
        if pathpoly is not None:
            xxyy_adc=np.array([gridx, gridy]).T
            notinpolyg = ~ pathpoly.contains_points(xxyy_adc).reshape(len(gridx)) # .shape)
        for y in range(0,yl): # Performing like this ensures proper Fortran style ordering
            gy = gridy[y]
            gx = gridx[y]
            zval=model(gx,gy).item(0) # ONLY ONE VALUE
            d.append([gx, gy, zval])
        df = pd.DataFrame(d,columns=['LON','LAT','VAL'])
        if pathpoly is not None: df.loc[notinpolyg,'VAL']=0
    return df
##
## Potential use methods for quick testing of new interpolation models
##

def generic_grid():
    """
    Build 2D grid for testing the interpolation_model_transform.
    this is a style=grid type data set

    Returns:
        df: DataFrame with columns ['LON','LAT']
    """
    utilities.log.debug('Generating points for surface figure.')
    lowerleft_x = -100
    lowerleft_y =5 
    res = .1  # resolution in deg
    nx = 500
    ny = 500
    x = np.arange(lowerleft_x, lowerleft_x + nx * res, res)
    y = np.arange(lowerleft_y, lowerleft_y + ny * res, res)
    coords = {'LON':x.tolist(), 'LAT':y.tolist()}
    return coords 

def test_interpolation_fit(df_source, df_land_controls=None, df_water_controls=None, cv_splits=5, nearest_neighbors=3, stratified_header_name=None):
    """
    df_source is the set of stationids and their values (errors)

    kNN fitting land control points has been deprecated.

    If df_water_controls are provided, they will be combined with the training set+(opt) land_controls prior to 
    the interpolation fit

    The MEAN Squared Error computation is only performed using the test set of stations. If we included the water_controls
    that would down-bias the statistics (they are always zero and always correct). The land_controls are also leftout but 
    tested separately

    NOTE For test _interpolation we do not need the polygon path because we are only comparing station values to predictions

    Parameters:
        df_source (DataFrame) stationds x errors
        df_land_controls: (DataFrame) stationsids x values(==0). Will be updated using KNN
        df_water_controls: (DataFrame) stationsids x values(==0) 
        cv_splits: (int) Type of Cross-validation splitting
        nearest_neighbors: (int) How many gauges should be applied to KNN compute land_control nodes 
        stratified_header_name: (str) Name iof an existing header form which to stratify the KFolds. No checking is made
            to enforce a minimum number of members to train and test

    Returns:
        kf_dict: (dict) Statistical values for the full data set. Returns aveMSE, for each CV fold and overall best_cnt
        kf_cntr_dict: (dict) Statistical values for the only the land control nodes. Returns aveCnrlMSE, for each CV fold and overall best_cnt
    """
    utilities.log.info('Initiating the CV testing: Using station drop outs to check overfitting. cv_splits {}, knn {}'.format(cv_splits, nearest_neighbors))
    df_source_drop = df_source.dropna()
    utilities.log.info('test_interpolation_fit: Nan Removal. Init {}. after {}'.format(len(df_source), len(df_source_drop)))

    if stratified_header_name is None:
        kf = KFold(n_splits=cv_splits)
        kfs = kf.split(df_source_drop)
        folds = kf.get_n_splits(df_source_drop)
    else:
        try:
            kf = StratifiedKFold(n_splits=cv_splits)
            folds=kf.get_n_splits()
            kfs = kf.split(df_source_drop,  df_source_drop[stratified_header_name])
        except Exception as ex:
            print('Stratified header_id error {}'.format(ex.args))
            sys.exit(1)

    kf_dict = dict([("fold_%s" % i,[]) for i in range(1, folds+1)])
    kfcntl_dict = dict([("Cnrl_fold_%s" % i,[]) for i in range(1, folds+1)])
    fold = 0
    for train_index, test_index in kfs: # kf.split(df_source_drop):
        print(f'TRAIN {train_index}, TEST {test_index}')
        train_list = list()
        fold += 1
        print('Fold: {}'.format(fold))
        df_source_train, df_source_test = df_source_drop.iloc[train_index], df_source_drop.iloc[test_index]
        train_list.append(df_source_train)
        if df_water_controls is not None:
            train_list.append(df_water_controls)
        if df_land_controls is not None:
        #    df_land_controls_train = knn_fit_control_points(df_source_train, df_land_controls, nearest_neighbors=nearest_neighbors)
        #    ktest = len(df_source_test) if len(df_source_test) < nearest_neighbors else nearest_neighbors
        #    df_land_controls_test = knn_fit_control_points(df_source_test, df_land_controls, nearest_neighbors=ktest)
            utilities.log.info('No KNN performed on the land control points')
            train_list.append(df_land_controls) 
        utilities.log.info('test_interpolation_fit: Num of DFs to combined for fit {}'.format(len(train_list)))
        df_Train=pd.concat(train_list, ignore_index=True)
        df_Test=df_source_test # pd.concat(test_list, ignore_index=True)

        # Start the computations
        model = interpolation_model_fit(df_Train, interpolation_type='LinearRBF')
        test_coords = {'LON': df_Test['LON'].to_list(), 'LAT': df_Test['LAT'].to_list()}
        df_fit_test_coords = interpolation_model_transform(test_coords, model=model, input_grid_type='points') # ['LON','LAT','VAL'])
        test_mse = mean_squared_error(df_Test['VAL'], df_fit_test_coords['VAL'])
        kf_dict["fold_%s" % fold].append(test_mse)

        # Now check the fitting of just the land_controls using the main model: Exclude the stations and water clamps
        if df_land_controls is not None:
            cntl_x = df_land_controls['LON']
            cntl_y = df_land_controls['LAT']
            land_coords = {'LON':cntl_x[:].tolist(),'LAT':cntl_y[:].tolist()}
            df_fit_land_coords = interpolation_model_transform(land_coords, model=model, input_grid_type='points') 
            testcntl_mse = mean_squared_error(df_land_controls['VAL'],df_fit_land_coords['VAL'])
        else:
             testcntl_mse =-99999
        kfcntl_dict["Cnrl_fold_%s" % fold].append(testcntl_mse)

    best_score = min(kf_dict.values())
    bestcntl_score = min(kfcntl_dict.values())
    #
    kf_dict["best"]=best_score
    kf_dict["avgMSE"] = 0.0
    for i in range(1, folds+1):
        kf_dict["avgMSE"] += kf_dict["fold_%s" % i][0]
    kf_dict["avgMSE"] /= float(folds)
    kfcntl_dict["best_cntl"]=bestcntl_score
    kfcntl_dict["avgCnrlMSE"] = 0.0
    for i in range(1, folds+1):
        kfcntl_dict["avgCnrlMSE"] += kfcntl_dict["Cnrl_fold_%s" % i][0]
    kfcntl_dict["avgCnrlMSE"] /= float(folds)
    return {**kf_dict, **kfcntl_dict}, best_score

def main(args):
    """
    This reads pre-generated files containing station errors, land controls and water clamps. The land controls
    have not been processed usng KNN.
    """

    config = utilities.init_logging(subdir=None, config_file='../config/main.yml')

    try:
        df_source = pd.read_pickle(args.source_file)
        df_land_controls = pd.read_csv(args.land_control_file,header=0)
        df_water_controls = pd.read_csv(args.water_control_file,header=0)
        df_land_controls.rename(columns={'lat':'LAT','lon':'LON','val':'VAL'}, inplace=True)
        df_water_controls.rename(columns={'lat':'LAT','lon':'LON','val':'VAL'}, inplace=True)
    except Exception as ex:
        print('LOAD error {}'.format(ex.args))
        sys.exit(1)

    full_scores, best_scores = test_interpolation_fit(df_source=df_source, df_land_controls=df_land_controls, df_water_controls=df_water_controls, cv_splits=5, nearest_neighbors=3)
    
    # Example of using generic grid to get an interpolation surface
    df_combined=combine_datasets_for_interpolation([df_source, df_water_controls, df_land_controls])
    model = interpolation_model_fit(df_combined, fill_value=0.0, interpolation_type='LinearRBF')
    adc_plot_grid = generic_grid()
    df_plot_transformed = interpolation_model_transform(adc_plot_grid, model=model, input_grid_type='grid')

    # Save the full scores for future reference
    utilities.log.debug('Save scores')
    df_full_scores = pd.DataFrame.from_dict(full_scores)
    df_full_scores.to_pickle('full_scores.pkl')
    utilities.log.debug('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--source_file', action='store', dest='source_file', default='../interpolate_testing/station_errors_test_interpolation.pkl',
                        type=str, help='PKL file with station errors')
    parser.add_argument('--land_control_file', action='store', dest='land_control_file', default='../interpolate_testing/land_controls_interpolation.pkl',
                        type=str, help='PKL file with land_control points - no KNN applied')
    parser.add_argument('--water_control_file', action='store', dest='water_control_file', default='../interpolate_testing/water_controls_interpolation.pkl',
                        type=str, help='PKL file with water_control points')
    parser.add_argument('--stratified_header_name', action='store', dest='stratified_header_name', default=None,
                        type=str, help='Header for performing stratified sampling')
    args = parser.parse_args()

if __name__ == '__main__':
    sys.exit(main(args))

