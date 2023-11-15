========================
Connecting to Floodwater
========================

`Floodwater <https://waterinstitute.github.io/floodwater/index.html>`_ has been instrumented with hooks to include a dynamic water level correction (DWLC) surface based on error analysis over previous analysis cycles.  To turn this capability on, include the following in a `Floodwater suite config yaml file <https://waterinstitute.github.io/floodwater/configuration_files.html#suite-configuration-file>`_: 

.. code-block:: yaml

    models:
    #...ADCIRC module configuration
    adcirc:
      #...Set the data assimilation options
      data_assimilation:
        enabled: true
        configuration_file: PATHTO/ecflow_configs/da/data_assimilation.yaml
        run_non_da_forecast: true

To turn off the DA, set both enabled and run_non_da_forecast to false.

===================================
**data_assimilation.yaml** File
===================================

The configuration_file **data_assimilation.yaml** needs to contain the following, adjusted to the local environment (**PATHTO**):

.. code-block:: yaml

   LOGGING: true
   LOGLEVEL: DEBUG
   rundir: "./adda"
   max_lookback_cycles: 8
   min_lookback_cycles: 2
   venv: adda
   dwlc_filename: "da_error_surface.dat"
   addahome: "PATHTO/adda_for_floodwater/"
   mapfile: "PATHTO/adda_for_floodwater/gridmap/grid_to_stationfile_maps.yml"

Configuration details: 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* LOGGING
   turn logging on (recommended) or off (not recommended)0
* LOGLEVEL
    log output detail.  INFO or DEBUG
* rundir
    output directory, relative to FLOODWATER archive directory, for the current cycle's ADDA results
* max_lookback_cycles
    maximum number of analysys (nowcast) cycles to use for error calculations
* min_lookback_cycles
    minimum number of analysys (nowcast) cycles to use for error calculations
* venv
    virtual python environment name built for ADDA
* dwlc_filename
    name of the file to contain the error surface.  Currently hardwired in Floodwater to "da_error_surface.dat".  I.e., do not change.
* addahome
    path to the installation of adda_for_floodwater
* mapfile
    path to the grid-to-station mapping file

