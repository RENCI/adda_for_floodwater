## The ADCIRC Data Assimilator (ADDA) in the ECFLOW/Floodwater environment

Version 0.2

Blanton/UNC-RENCI <br>
Bunya/UNC-CRC

ADDA generates smooth surfaces of an error field, computed between NOAA / NOS observations and corresponding ADCIRC output.  It does this by computing the average errors from time series of errors at each specified NOAA gague location and building a nearest-neighbor-based model to calculate the surface values at all ADCIRC grid nodes.  The surface is constrained by including offshore and land "control points".  NOAA / NOS gauge observations are retrieved using the noaa-coops python package (git@github.com:GClunies/noaa_coops.git).  

Currently, ADDA requires a specific python environment, but will eventually be embedded in the "thirdparty" section of the Floodwater package as a submodule.

ADDA documented is here: https://renci.github.io/adda_for_floodwater

#### TODO:
- Mv gridmap files to main config dir
- Send logging data to the main floodwater logging stream
- Put ADDA into the "thirdparty" structure in Floodwater
- Describe surface generation

