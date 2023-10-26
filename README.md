## The ADCIRC Data Assimilator (ADDA) in the TWI Floodwater/ECFLOW environment

25 Oct 2023

Blanton/UNC-RENCI <br>
Bunya/UNC-CRC

ADDA generates smooth surfaces of an error field, computed between NOAA / NOS observations and corresponding ADCIRC output.  It does this by computing the average errors from time series of errors at each specified NOAA gague location and building a nearest-neighbor-based model to calculate the surface values at all ADCIRC grid nodes.  The surface is constrained by including offshore and land "control points".  NOAA / NOS gauge observations are retrieved using the noaa-coops python package (git@github.com:GClunies/noaa_coops.git).  

TODO: Describe surface generation

Currently, ADDA requires a specific python environment, but will eventually be embedded in the "thirdparty" section of the Floodwater package as a submodule.

### Installation / Virtual Env:

- Clone this repo.  This locations is referred to as **PATHTO** below.
- Make a python virtual environment (venv) with the requirements.txt file, called **adda**.  The actual name of the venv does not matter, as long as the name is speficied correctly in the **data_assimilation.yaml** file (see below).
  - conda create --name adda --file requirements.txt
- Several required packages are not available through conda.  These need to be pip-installed:
  - conda activate adda
  - pip install -r pip.reqs.txt
- Add ADDA paths to the virtual environment.  Either:
  - Make a conda.pth in "envs/adda/lib/python3.8/site-packages/" that contains:
    - **PATHTO**/adda_for_floodwater
    - **PATHTO**/adda_for_floodwater/adda
    - Note that the path to the venv's site-packages location may different than that above.  Modify as needed. 
- or (probably better):
  - Add the PYTHONPATH variable to the conda env. Activate the **adda** environment, add path:
<pre>conda activate adda
conda env config vars set PYTHONPATH=PATHTO/adda_for_floodwater:PATHTO/adda_for_floodwater/adda</pre>
    

### Connecting to Floodwater

Floodwater has been instrumented with hooks to include a dynamic water level correction (DWLC) surface based on error analysis over previous analysis cycles.  To turn this capability on, include the following in the suite config yaml file: 

<pre>
...
...
models:

    #...ADCIRC module configuration
    adcirc:
    ...
    ...
      #...Set the data assimilation options
      data_assimilation:
        enabled: true
        configuration_file: PATHTO/ecflow_configs/da/data_assimilation.yaml
        run_non_da_forecast: true
...
...
</pre>

### The **data_assimilation.yaml** File

The configuration_file **data_assimilation.yaml** needs to contain the following, suitably adjusted to the local environment (path_to):
<pre>
LOGGING: true
LOGLEVEL: DEBUG
rundir: "./adda"
max_lookback_days: 2
venv: adda
dwlc_filename: "da_error_surface.dat.1"
addahome: "PATHTO/adda_for_floodwater/"
mapfile: "PATHTO/adda_for_floodwater/gridmap/grid_to_stationfile_maps.yml"
</pre>

### The grid_to_stationfile_maps.yml

This yaml file contains pointers to files ADDA needs in order to retrieve specific NOAA/NOS station observations and corresponding ADCIRC grid nodes at which to compute the errors, and specify "control" points in open water and on land to constraint the resulting surface evaluation.  For each ADCIRC grid, it specifies locations of 3 files:

- **NOAA_STATIONS** - csv file containing NOAA / NOS ids and ADCIRC grid nodes.  The file can have lots of info in it, but must have stationid and Node as columns in the first line, and "units" second line.  Other columns are ignored. E.g., 
<pre>
serial_nr,stationid,stationname,state,vertical_datum,NOAA lon,NOAA lat,navd_to_msl [ft],navd_to_msl [m],Datum Source,lon,lat,bathy,Element,Node
-,-,-,-,-,deg,deg,ft,m,-,deg,deg,m MSL,-,-
1,8410140,Eastport Passamaquoddy Bay,ME,MSL,-66.982315,44.903300,0.2300,0.0701,NOAA gage,-66.982315,44.903300,36.7010,2034577,1034613
2,8413320,Bar Harbor,ME,MSL,-68.205000,44.391700,0.2760,0.0841,NOAA gage,-68.203000,44.393000,4.5782,2287004,1162475 
...
</pre>
- **LANDCONTROL** - csv file containing "land" lon/lat locations and surface values to force the surface toward:
<pre>lon,lat,val
-60,50,0.
-80,55,0.
-100,50,0.
-110,30,0.
</pre>
- **WATERCONTROL** - csv file containing open-water lon/lat locations and surface values, same format as the LANDCONTROL file

The ADDA repo and **grid_to_stationfile_maps.yml** contain files for the HSOFS and ec95d grids.  Users can add grids as needed, following the above formatting and information. 
<pre>
GRIDMAP: &gridmap
 HSOFS:
    NOAA_STATIONS: "PATHTO/adda_for_floodwater/gridmap/HSOFS_stations_V2.csv"
    LANDCONTROL: "PATHTO/adda_for_floodwater/gridmap/hsofs_land_control_list.dat"
    WATERCONTROL: "PATHTO/adda_for_floodwater/gridmap/hsofs_water_control_list.dat"
 EC95D:
    NOAA_STATIONS: "PATHTO/adda_for_floodwater/gridmap/ec95d_stations.V1.csv"
    LANDCONTROL: "PATHTO/adda_for_floodwater/gridmap/ec95d_land_control_list.dat"
    WATERCONTROL: "PATHTO/adda_for_floodwater/gridmap/ec95d_water_control_list.dat"
</pre>
Gridnames in the grid_to_stationfile_maps.yml file must be in **ALL_CAPS**, and must match (modulo the case) the grid name as defined in the grid's config json file, as in (e.g., PATHTO/ecflow_configs/ec95d/config.json) 
<pre>
{
   "name": "ec95d",
   "dt": 10,
   ...
   ...
</pre>
  
#### TODO:
- mv gridmap files to main config dir
- send logging data to the main floodwater logging stream
- put ADDA into the "thirdparty" structure in Floodwater

