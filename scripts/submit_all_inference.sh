#!/bin/bash

CONFIGDIR=${1}
CONFIGDIR=${CONFIGDIR%/}
TAG=${CONFIGDIR#*configs_}

echo "Tag is ${TAG}"

for INIFILE in ${CONFIGDIR}/config_run*[0-9].ini; do
    CONFIGFILE=${INIFILE#${CONFIGDIR}/}
    INJID=${INIFILE#*config_run}
    INJID=${INJID%.ini}
    OUTFILE=out${INJID}.hdf
    INJFILE=injection_${INJID}.hdf
    echo "Submitting ${CONFIGFILE}"
    echo "Transferring ${INJFILE}"
    condor_submit submitfile.sub config_file=${CONFIGFILE} injfile=${INJFILE} outfile=${OUTFILE} 
done