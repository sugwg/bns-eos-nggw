#! /bin/bash

set -x

export PYCBC_DETECTOR_CONFIG=/home/abandopa/projects/github-bns-eos-population/eos_population/network_asharp_only/custom_detector_network/detectors.ini

pycbc_config_file=${1}
pycbc_inj_file=${2}
pycbc_output_file=${3}
nprocs=62

echo "Using ${pycbc_config_file} as configuration file"
echo "Using ${pycbc_inj_file} as injection file"
echo "Writing output to ${pycbc_output_file}"

# expand EOS data
tar -xzvf eos_nsat.tar.gz
pwd
hostname
ls -al

pycbc_seed=11185

echo "Using ${nprocs} processors"
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed ${pycbc_seed} \
    --config-file ${pycbc_config_file} \
    --output-file ${pycbc_output_file} \
    --processing-scheme cpu \
    --nprocesses ${nprocs}
