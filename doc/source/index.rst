
.. image:: ../img/logo.png
|
|
===================================================================================================================================
The ADCIRC Data Assimilator (ADDA) in the ECFLOW/`Floodwater <https://waterinstitute.github.io/floodwater/index.html>`_ environment
===================================================================================================================================

ADDA generates smooth surfaces of an error field, computed between NOAA / NOS observations and corresponding ADCIRC output.  It does this by computing the average errors from time series of errors at each specified NOAA gauge location and building a nearest-neighbor-based model to calculate the surface values at all ADCIRC grid nodes.  The surface is constrained by including offshore and land "control points".  NOAA / NOS gauge observations are retrieved using the `noaa-coops python package <git@github.com:GClunies/noaa_coops.git>`_ .  

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    installation.rst
    connect2fw.rst
    gridmap.rst
    acknowledgements.rst
