#!/bin/bash


#fw_dir="/scratch/bblanton/ecflow_output/dev/ec95d_dev_da_idalia/"
#fw_dir="/home/bblanton/storm_surge/ec95d_dev_da_idalia.test1"
fw_dir="/scratch/sbunya/ecflow_output/ncsc123_gfs_da/"

ADCIRC_DA_CONFIG_FILE="/home/bblanton/ecfc/da/data_assimilation.yaml"
#ADCIRC_DA_CONFIG_FILE="/home/sbunya/ecflow_configs/da/data_assimilation.yaml"

#...Retrieve config variables from yaml
venv=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value venv`
addahome=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value addahome`
dwlcfile=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value dwlc_filename`
path2adda="$addahome/adda/"

echo "venv=$venv"
echo "adda=$path2adda/adda.py"

RUN="conda run --no-capture-output"

#current_time="--current_time advisory_017"
current_time="--current_time 20231118/hour_12"

#grid="--gridname ec95d"
grid="--gridname NCSC_SAB_V1.23"
#meteo="--met nhc"
meteo="--met gfs"
com="$RUN -n $venv $path2adda/adda.py  $current_time $meteo --da_config_file $ADCIRC_DA_CONFIG_FILE $grid --fw_arch_dir $fw_dir"
echo $com
$com

