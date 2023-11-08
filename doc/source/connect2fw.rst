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
The **data_assimilation.yaml** File
===================================

The configuration_file **data_assimilation.yaml** needs to contain the following, suitably adjusted to the local environment (path_to):

.. code-block:: yaml

   LOGGING: true
   LOGLEVEL: DEBUG
   rundir: "./adda"
   max_lookback_days: 2
   venv: adda
   dwlc_filename: "da_error_surface.dat.1"
   addahome: "PATHTO/adda_for_floodwater/"
   mapfile: "PATHTO/adda_for_floodwater/gridmap/grid_to_stationfile_maps.yml"

