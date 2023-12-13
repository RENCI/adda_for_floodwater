
#################################### 
**grid_to_stationfile_maps.yml**
#################################### 

This yaml file contains pointers to files ADDA needs in order to retrieve specific NOAA/NOS station observations and corresponding ADCIRC grid nodes at which to compute the errors, and specify "control" points in open water and on land to constraint the resulting surface evaluation.  For each ADCIRC grid, it specifies locations of 3 files:

* **NOAA_STATIONS** - csv file containing NOAA / NOS ids and ADCIRC grid nodes.  The file can have lots of info in it, but must have stationid and Node as columns in the first line, and "units" second line.  Other columns are ignored. E.g., 

.. code-block:: text

  serial_nr,stationid,stationname,state,Element,Node
  -,-,-,-,-,-
  1,8410140,Eastport Passamaquoddy Bay,ME,2034577,1034613
  2,8413320,Bar Harbor,ME,2287004,1162475
  3,8418150,Portland,ME,2314208,1176249
  4,8419317,Wells,ME,2327314,1182888
  ...

* **LANDCONTROL** - csv file containing "land" lon/lat locations and surface values to force the surface toward:

.. code-block:: text

  lon,lat,val
  -60,50,0.
  -80,55,0.
  -100,50,0.
  -110,30,0.

* **WATERCONTROL** - csv file containing open-water lon/lat locations and surface values, same format as the LANDCONTROL file

The ADDA repo and **grid_to_stationfile_maps.yml** contain files for the HSOFS and ec95d grids.  Users can add grids as needed, following the above formatting and information.   This is the "default" file content; **PATHTO** will need to be set by the user.

.. code-block:: yaml

  GRIDMAP: &gridmap
    HSOFS:
      NOAA_STATIONS: "PATHTO/adda_for_floodwater/gridmap/HSOFS_stations_V2.csv"
      LANDCONTROL: "PATHTO/adda_for_floodwater/gridmap/hsofs_land_control_list.dat"
      WATERCONTROL: "PATHTO/adda_for_floodwater/gridmap/hsofs_water_control_list.dat"
   EC95D:
      NOAA_STATIONS: "PATHTO/adda_for_floodwater/gridmap/ec95d_stations.V1.csv"
      LANDCONTROL: "PATHTO/adda_for_floodwater/gridmap/ec95d_land_control_list.dat"
      WATERCONTROL: "PATHTO/adda_for_floodwater/gridmap/ec95d_water_control_list.dat"


Gridnames in the grid_to_stationfile_maps.yml file must be in **ALL_CAPS**, and must match (modulo the case) the grid name as defined in the grid's config json file, as in (e.g., PATHTO/ecflow_configs/hsofs/config.json) 

.. code-block:: yaml

   name: hsofs
   dt: 2
  
