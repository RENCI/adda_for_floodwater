##
## Set of the codes
##

# 1) Fetch the AST base codes. 
git clone https://github.com/RENCI/AST.git

# 2) Fetch the EDS aps which include adda and reanalysis
git clone https://github.com/RENCI/EDSASTAPPS.git

##
## Set up you working env
## RUNTIMEDIR is created under your working dir and is the top level dir where all results get stored
##

export PYTHONPATH="$PWD/AST:$PWD/EDSASTAPPS"
export RUNTIMEDIR=./MYRUNS
SRC="$PWD/EDSASTAPPS/adda"

##
## Setup a calculation (on ht1) To perform the full offset file computation for 1 or more years
##

mkdir adda_testing
cd adda_testing

##
## Fetch the station/clamps/grid map data
## I only clone into this directory to obviate the need to manually change the contents of the
## grid_to_stationfile_maps.yml internal pointers. 
##

git clone https://github.com/RENCI/AST_gridstations.git

##
## Execute the pipeline
##

python $SRC/adda.py  --instance_name 'hsofs-nam-bob-2021' --ndays -2  --timeout '2022-02-20 00:00:00' --return_sample_min 60  --instance_name 'hsofs-nam-bob-2021' --gridname 'hsofs' --ensemble='nowcast' --fort63_style --config_name '../EDSASTAPPS/adda/secrets/url_framework.yml' --map_file './AST_gridstations/interpolation_stationlist/grid_to_stationfile_maps.yml' --main_yamlname '../EDSASTAPPS/adda/config/main.yml'



