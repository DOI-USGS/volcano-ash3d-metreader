#!/bin/bash

# This bash script is run for the target 'make check'

rc=0
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi

# Check to see if the NCEP data for 1980 is present
ls -1r ${WINDROOT}/NCEP/1980/air.1980.nc
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Warning: Could not find NCEP data for 1980"
  echo "Only the Sonde test cases will be run."
  echo "To download the NCEP data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 1980"
  RUNNCEP=F
  nwin=1
else
  RUNNCEP=T
  nwin=2
fi

pushd tests
sh clean.sh
sh run_tests_Sonde.sh
sh clean.sh
if [[ "$RUNNCEP" -eq T ]] ; then
  sh run_tests_NCEP.sh
  sh clean.sh
fi
popd

echo "Note: Any test failures should be inspected."
echo "      This script only uses the cmp command for comparison."
echo "      Values can be off by 0.001% using different compilation flags."


