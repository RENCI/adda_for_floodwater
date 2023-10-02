##
## Setup the codes. This is a procedure one might use to build a functional ADDA/AST system without downloading
## all of the EDSASTAPPS codes. Many steps can be modified such as where to clone AST_gridstations but this should
## be a decent start
##

# 1) A rather newer version of git is required. 1.8 doesn't work, git version 2.39.2 does

# 2) I needed to update my Centos7. YMMV
sudo yum -y install https://packages.endpointdev.com/rhel/7/os/x86_64/endpoint-repo.x86_64.rpm
yum remove git
sudo yum install git

# 3) Decide where to build the packages
mkdir /projects/sequence_analysis/vol1/prediction_work/PARTIALADDA
cd /projects/sequence_analysis/vol1/prediction_work/PARTIALADDA

#4) Will need AST proper
git clone https://github.com/RENCI/AST.git

#5) Will also need AST_gridstations
git clone https://github.com/RENCI/AST_gridstations.git


# 6) Fetch only the ADDA bits and pieces required for adda
#    into the target directory of choice

git clone --filter=blob:none --no-checkout --depth 1 --sparse https://github.com/RENCI/EDSASTAPPS.git
cd EDSASTAPPS
git sparse-checkout add adda io_utilities gridmap
git checkout

# 7) From within PARTIALADDA, do a test of adda

# Set up the env to find AST
export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/PARTIALADDA/AST:/projects/sequence_analysis/vol1/prediction_work/PARTIALADDA/EDSASTAPPS
export RUNTIMEDIR=.

# Run a job

SRC='/projects/sequence_analysis/vol1/prediction_work/PARTIALADDA/EDSASTAPPS/adda'

# Example: Pass in a valid URL
python $SRC/adda.py --input_url "http://tds.renci.org/thredds/dodsC/2022/al09/12/hsofs/hatteras.renci.org/hsofs-al09-bob/nowcast.nodwlc/fort.63.nc" --ndays -4 --gridname 'hsofs' --map_file ./AST_gridstations/interpolation_stationlist/grid_to_stationfile_maps.yml --ensemble 'nowcast.nodwlc' --fort63_style


# Another example: Using YAML to construct proper URLS
# NOTE: This assumes you have a properly formatted url_framework.yml file
python $SRC/adda.py  --instance_name 'hsofs-nam-bob-2021' --ndays -2  --timeout '2022-02-20 00:00:00' --return_sample_min 60  --instance_name 'hsofs-nam-bob-2021' --gridname 'hsofs' --ensemble='nowcast' --fort63_style --config_name './EDSASTAPPS/adda/secrets/url_framework.yml' --map_file './AST_gridstations/interpolation_stationlist/grid_to_stationfile_maps.yml' --main_yamlname './EDSASTAPPS/adda/config/main.yml'
