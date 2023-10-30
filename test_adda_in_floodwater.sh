
export ASTHOME="$HOME/GitHub/RENCI/adda_for_floodwater"
export ADDAHOME="$HOME/GitHub/RENCI/adda_for_floodwater/adda"
export PYTHONPATH="$ADDAHOME:$ASTHOME"

fw_dir="/scratch/bblanton/ecflow_output/dev/ec95d_dev_da/"

da_config_file="--da_config_file /home/bblanton/ecfc/da/data_assimilation.yaml"
ndays="--ndays -2" 
venv="adda" 
inputFile="--input_url adda.filelist"
adda="/home/bblanton/GitHub/RENCI/adda_for_floodwater/adda/adda.py"
RUN="conda run --no-capture-output"

gridname="--gridname ec95d"

#$RUN -n $venv $adda $da_config_file $inputFile $ndays   $gridname --fw_arch_dir $fw_dir
$RUN -n $venv $adda $da_config_file   $gridname --fw_arch_dir $fw_dir


