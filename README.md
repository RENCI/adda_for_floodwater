### Installing / running the ADCIRC Data Assimilator (ADDA)

19 Oct 2023
Blanton/RENCI

#### Installation:

- Clone this repo
- Make a python virtual environment with the requirements.txt file, called adcirc_DA

#### Connecting to Floodwater

Floodwater has been instrumented with hooks to include a dynamic water level surface based on error analysis over previous analysis cycles.  To turn this capability on, include the following in the suite config yaml file: 

<pre>
models:

    #...ADCIRC module configuration
    adcirc:
    ...
    ...
      #...Set the data assimilation options
      data_assimilation:
        enabled: true
        configuration_file: /home/bblanton/ecflow_configs/da/data_assimilation.yaml
        run_non_da_forecast: true

</pre>

The configuration_file *data_assimilation.yaml* needds to contain the following, suitably adjusted to the local environment:
<pre>
    LOGGING: true
    LOGLEVEL: DEBUG
    rundir: "./adda"
    max_lookback_days: 2
    venv: adcirc_DA
    addahome: "<path to>/adda_for_floodwater/"
    mapfile: "<path to>/adda_for_floodwater/gridmap/grid_to_stationfile_maps.yml"
</pre>

#### The grid_to_stationfile_maps.yml
This file contains file locations that ADDA needs in order to 
- retrieve specific NOAA/NOS station observations
- know the correspondence between the NOAA stations and the ADCIRC grid node at which to compute the errors
- specify "control" points in open water and on land
  
#### TODO:
- mv gridmap files to main config dir

