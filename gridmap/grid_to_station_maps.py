#!/usr/bin/env python
#

import os,sys
import pandas as pd
from utilities.utilities import utilities as utilities
from argparse import ArgumentParser

##
## Support methods to find an appropriate dataframe formatted station list  
##

# SOURCES=['NOAA_STATIONS','CONTRAILS_RIVERS','CONTRAILS_COASTAL','NDBC_BUOYS']
# The use may define anyt new key word to the grid yaml and access that data through here

def read_specified_filename_from_map(gridname=None, mapfile=None, datatype=None):
    """
    Parameters:
        datatype: (str) is one of SOURCES" Eg 'NOAA_STATIONS','CONTRAILS_RIVERS','CONTRAILS_COASTAL','NDBC_BUOYS'
        gridname: (str) a gridname such as hsofs
        a caller speciific gridmap file. If None choose the default one
            If default, then must build a filename path as well
    Return: fullpath name of desitred file or None
    """
    if mapfile is None:
        map_yamlname=os.path.join(os.path.dirname(__file__), '../supporting_data', 'grid_to_stationfile_maps.yml')
    else:
        map_yamlname=mapfile
    print(map_yamlname)

    # Does mapfile exists. if not abort
    try:
        map_config = utilities.load_config(map_yamlname)
    except OSError as e:
        utilities.log.error('grid_to_station_maps: No such map file {}'.format(map_yamlname))
        list_mapfile_contents(map_yamlname)
        sys.exit(1)
    utilities.log.info('Found a grid map file {}'.format(map_yamlname))

    # Does Gridname exist. If not abort
    try:
        test_grid_name = map_config['GRIDMAP'][gridname.upper()]
    except Exception as e:
        utilities.log.error('grid_to_station_maps: Error finding grid. {}'.format(e,gridname))
        print('Failed: choices for selected gridnames are')
        list_mapfile_contents(map_yamlname)
        sys.exit(1)

    # Soft Fetch the requested value. return None as an option (usually land/water control files)
    dataFile = map_config.get('GRIDMAP',{}).get(gridname.upper(),{}).get(datatype.upper())
    if mapfile is None and dataFile is not None:
        dataFile=os.path.join(os.path.dirname(__file__), '../supporting_data', dataFile)
    #print('data file is {}'.format(dataFile))
    return dataFile

def find_station_list_from_map( gridname=None, datatype='NOAA_STATIONS', mapfile=None ):
    """
    This method abstracts finding appropriate station lists for use for EDS applications.
    Each dataframe MUST contains a column with header "stationid"
    Each dataframe MAY also contain a column with the header "Node"

    If the calling routine is scheduled to fetch ADCIRC water levels using a fort63_style 
    approach, then the Node column MUST be populated. It is not required if selecting a fort61_style
    approach
    
    The default mapfile, points to known map files in the ../supporting_data directory
    User specified mapfiles, should probably include full pathnames in the yml file
    """
    stationFile = read_specified_filename_from_map(gridname, mapfile, datatype=datatype)
    print('station file is {}'.format(stationFile))
    fort63_compliant=False
    if stationFile is not None:
        fort63_compliant = check_for_fort63_compliance(stationFile)
    return stationFile, fort63_compliant

def find_land_control_points_from_map( gridname=None, mapfile=None ):
    """
    This method abstracts finding appropriate land_control lists for use for EDS applications.
    Each dataframe MUST contains the columns: lon,lat,val, where val==0. These get KNN specified
    as part ofa subsequenty interpolation scheme
    The default mapfile, points to known map files in the ../supporting_data directory
    User specified mapfiles, should probably include full pathnames in the yml file
    """
    landcontrolFile = read_specified_filename_from_map(gridname, mapfile, datatype='LANDCONTROL')
    print('land control file is {}'.format(landcontrolFile))
    return landcontrolFile

def find_water_control_points_from_map( gridname=None, mapfile=None ):
    """
    This method abstracts finding appropriate water_control lists for use for EDS applications.
    Each dataframe MUST contains the columns: lon,lat,val, where val==0. These are always zero clamps
    The default mapfile, points to known map files in the ../supporting_data directory
    User specified mapfiles, should probably include full pathnames in the yml file
    """
    watercontrolFile = read_specified_filename_from_map(gridname, mapfile, datatype='WATERCONTROL')
    print('water control file is {}'.format(watercontrolFile))
    return watercontrolFile

def find_secondary_water_control_points_from_map( gridname=None, mapfile=None ):
    """
    This method abstracts finding appropriate secondary water_control lists for use for EDS applications.
    Each dataframe MUST contains the columns: lon,lat,val, where val==0. These are always zero clamps
    The default mapfile, points to known map files in the ../supporting_data directory
    User specified mapfiles, should probably include full pathnames in the yml file
    """
    watercontrolFile = read_specified_filename_from_map(gridname, mapfile, datatype='SECONDARY_WATERCONTROL')
    print('secondary water control file is {}'.format(watercontrolFile))
    return watercontrolFile

def check_for_fort63_compliance(stationfile):
    """
    Open the station file and check for the existance of the "Node" column.
    If not, return False"
    """
    try:
        df = pd.read_csv(stationfile, index_col=0, header=0, skiprows=[1])
    except Exception as e:
        utilities.log.error('grid_to_station_maps: compliance: {}'.format(e))
    avail_columns = df.columns.to_list()
    compliant = True if 'Node' in avail_columns else False
    return compliant

def list_mapfile_contents(mapfile):
    """
    Read and output all entries in the provided mapfile
    """
    if mapfile is None:
        map_yamlname=os.path.join(os.path.dirname(__file__), '../supporting_data', 'grid_to_stationfile_maps.yml')
    else:
        map_yamlname=mapfile
    map_config = utilities.load_config(map_yamlname)
    print('Grid choices')
    for item in map_config['GRIDMAP'].items():
        print(item)

def main(args):
    """
    This method abstracts finding appropriate station lists for use for EDS applications.
    Each dataframe MUST contains a column with header "stationid"
    Each dataframe MAY also contain a column with the header "Node"

    If the calling routine is scheduled to fetch ADCIRC water levels using a fort63_style 
    approach, then the Node column MUST be populated. It is not required if selecting a fort61_style
    approach
    
    The default mapfile, points to known map files in the ../supporting_data directory
    User specified mapfiles, should probably include full pathnames in the yml file
    """
    # Basic checks
    print('Search for station list: Selecting the grid {}'.format(args.gridname))

    stationfile,fort_63_compliant = find_station_list_from_map( args.gridname, args.mapfile)
    #fort_63_compliant=check_for_fort63_compliance(stationfile)
    list_mapfile_contents(args.mapfile)

    print('Obtained Stationfile {} for gridname {}. fort63 compliant {}'.format(stationfile,args.gridname,fort_63_compliant))
    print('Success. found station map file')

    # Repeat station list
    stationFile = find_station_list_from_map( args.gridname, args.mapfile)
    print('stationFile is {}'.format(stationFile))

    # TEST finding a landcontrol file
    landcontrolFile = find_land_control_points_from_map( args.gridname, args.mapfile)
    print('main: Landcontrol is {}'.format(landcontrolFile))

    watercontrolFile = find_water_control_points_from_map( args.gridname, args.mapfile)
    print('main: Watercontrol is {}'.format(watercontrolFile))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--gridname', action='store',dest='gridname', default='hsofs',
                        help='str: Select appropriate gridname Default is hsofs')
    parser.add_argument('--mapfile', action='store',dest='mapfile', default=None,
                        help='str: FQFN to find an alternative mapfilename.')
    args = parser.parse_args()
    sys.exit(main(args))
