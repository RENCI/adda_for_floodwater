#!/usr/bin/env python

import sys, os
import numpy as np
import pandas as pd
from utilities.utilities import utilities
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

vmax=0.5
vmin=-vmax
xlim=[-100, -50]
ylim=[5, 55]
N=20
#base_cmap='tab20c' # Or simply use None tab20c is also okay
#cmap= plt.cm.get_cmap('jet', N)

def save_plot_model(plot_grid=None, df_surface=None, df_stations=None, df_land_control=None, df_water_control=None, 
                    filename=None, plot_now=False, title=None):
    """
    Wrapper to take generated plt file and save to disk
    """
    if filename is None:
        utilities.log.error('save_plot_model: No filename provided to save the plot')
        return

    plt_data = plot_model(plot_grid=plot_grid, 
                          df_surface=df_surface, df_stations=df_stations, 
                          df_land_control=df_land_control, df_water_control=df_water_control, 
                          plot_now=plot_now, title=title)
    plt_data.savefig(filename, bbox_inches='tight')

def plot_model(plot_grid=None, df_surface=None, df_stations=None, df_land_control=None, df_water_control=None, plot_now=True, title=None):
    """
    Basic plotter to display the error field. 

    adc_plot_grid carries the (x,y) values spanning the data: df_surface
    So if, len(LAT) = 300, len(LON)=400, len(dim(df_surface)) is (400*300)

    Parameters:
        plot_grid: dict with Keys of 'LON','LAT'
        df_surface: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        df_stations: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        df_land_control: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        df_water_control: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        plot_now: bool. True display the plot. Else not
    Results:
        A plot 
    """
    coastline=np.loadtxt(os.path.join(os.path.dirname(__file__), "misc", "coarse_us_coast.dat"))

    #
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,10), dpi=144) #, sharex=True)
    cmap=discrete_cmap(N,'jet')

    # draw surface
    if df_surface is not None:
        x = plot_grid['LON']
        y = plot_grid['LAT']
        v = df_surface['VAL'].values
        v = v.reshape(-1,len(x)) # Reshapes to LAT column-major format
        mesh = ax.pcolormesh(x, y, v, shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)

    # add in points
    if df_stations is not None:
        stations_X=df_stations['LON'].values
        stations_Y=df_stations['LAT'].values
        stations_V=df_stations['VAL'].values
        ax.scatter(stations_X, stations_Y, s=50, marker='o',
                   c=stations_V, cmap=cmap, edgecolor='black',
                   vmin=vmin, vmax=vmax)
    
    if df_land_control is not None:
        land_X=df_land_control['LON'].values
        land_Y=df_land_control['LAT'].values
        land_V=df_land_control['VAL'].values
        ax.scatter(land_X, land_Y, s=30, marker='o', edgecolor='black')
                   #c=land_V, edgecolor='black',
                   #vmin=vmin, vmax=vmax)
    
    if df_water_control is not None:
        water_X=df_water_control['LON'].values
        water_Y=df_water_control['LAT'].values
        water_V=df_water_control['VAL'].values
        ax.scatter(water_X, water_Y, s=30, marker='x',c='black')
                   
    ax.plot(coastline[:,0],coastline[:,1],color='black',linewidth=1.5)

    ax.axis('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylabel('Latitude', fontsize=16)
    ax.set_xlabel('Longitude', fontsize=16)
    plt.grid(True)

    if title is not None:
        ax.set_title(title, fontsize=20)

    v1 = np.linspace(vmin, vmax, 11, endpoint=True)
    cbar=plt.colorbar(mesh,ticks=v1,orientation='vertical')
    cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in v1]) # add the labels

    if (plot_now):
        plt.show()
    return plt

def discrete_cmap(N, base_cmap='jet_r'):
    """
    Create an N-bin discrete colormap from the specified input map
    """

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def main(args):
    # Grab the test data
    f_land_control = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_land_controls.csv'
    f_water_control = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_water_controls.csv'
    f_surface = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_surface.csv'
    f_grid = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/adc_plot_grid.json'
    print("HERE HERE HERE HERE") 

    df_land_control=pd.read_csv(f_land_control,header=0,index_col=0)
    df_water_control=pd.read_csv(f_water_control,header=0, index_col=0)
    df_surface=pd.read_csv(f_surface,header=0, index_col=0)
    adc_coords = utilities.read_json_file(f_grid)
    plot_model(adc_plot_grid=adc_coords, df_surface=df_surface ,df_land_control=None,df_water_control=None)
    #plot_model(adc_plot_grid=adc_coords, df_surface=df_surface ,df_land_control=df_land_control,df_water_control=df_water_control)
    plt.show()

