#!/bin/bash


#fw_dir="/scratch/bblanton/ecflow_output/dev/ec95d_dev_da_idalia/"
#fw_dir="/scratch/bblanton/ecflow_output/dev/ec95d_dev_da/"
#fw_dir="/home/bblanton/storm_surge/ec95d_dev_da_idalia.test1"
fw_dir="/scratch/sbunya/ecflow_output/dev/ec95d_dev_da/"

#ADCIRC_DA_CONFIG_FILE="/home/bblanton/ecfc/da/data_assimilation_test.yaml"
ADCIRC_DA_CONFIG_FILE="/home/bblanton/ecfc/da/data_assimilation.yaml"
#ADCIRC_DA_CONFIG_FILE="/home/sbunya/ecflow_configs/da/data_assimilation.yaml"

#...Retrieve config variables from yaml
venv=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value venv`
venv="addaTest"
addahome=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value addahome`
dwlcfile=`cat $ADCIRC_DA_CONFIG_FILE | shyaml get-value dwlc_filename`
path2adda="$addahome/adda/"

echo "venv=$venv"
echo "adda=$path2adda/adda.py"

RUN="conda run --no-capture-output"

#meteo="--met nhc"
#current_time="--current_time advisory_017"

meteo="--met gfs"
current_time="--current_time 20231207/hour_00"

grid="--gridname ec95d"
#grid="--gridname NCSC_SAB_V1.23"
da_config_file="--da_config_file $ADCIRC_DA_CONFIG_FILE"
fw_arch_dir="--fw_arch_dir $fw_dir"
com="$RUN -n $venv $path2adda/adda.py  $current_time $meteo $da_config_file  $grid $fw_arch_dir"
echo $com
$com

