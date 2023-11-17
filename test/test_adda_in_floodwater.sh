#!/bin/bash

fw_dir="/scratch/bblanton/ecflow_output/dev/ec95d_dev_da_idalia/"
#fw_dir="/home/bblanton/storm_surge/ec95d_dev_da_idalia.test1"

ADCIRC_DA_CONFIG_FILE="/home/bblanton/ecfc/da/data_assimilation.yaml"

#...Retrieve config variables from yaml
venv=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value venv`
addahome=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value addahome`
dwlcfile=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value dwlc_filename`
path2adda="$addahome/adda/"

echo "venv=$venv"
echo "adda=$path2adda/adda.py"

RUN="conda run --no-capture-output"

current_time="--current_time advisory_017"

grid="--gridname ec95d"
meteo="--met nhc"
com="$RUN -n $venv $path2adda/adda.py  $current_time $meteo --da_config_file $ADCIRC_DA_CONFIG_FILE $grid --fw_arch_dir $fw_dir"
echo $com
$com

