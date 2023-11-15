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

Do not change the value of dwlc_filename.  This is currently hardwired in Floodwater.

LOGLEVEL can be set to "INFO" to decrease the verbosity of the ADDA logging output.

min_lookback_cycles and max_lookback_days control how many previous nowcast cycles are used in the error analysis.  

It is probably best to leave others alone, except for setting **PATHTO** appropriately.
