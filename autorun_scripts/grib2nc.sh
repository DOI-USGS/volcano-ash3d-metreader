#!/bin/bash

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

# Shell script that converts a grib2 file to netcdf using netcdf-java
# This script takes one command-line argument, the name of the grib file
rc=0
# Search for required packages
echo "Looking for latest netcdfAll in ~/ncj/"
ls -1r ~/ncj/netcdfAll*.jar
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find netcdfAll in ~/ncj/ rc=$rc"
  echo "Please make sure to put netcdfAll-[].jar in ~/ncj or specify the path."
  echo "The latest version can be downloaded from"
  echo "  https://www.unidata.ucar.edu/downloads/netcdf-java/"
  exit 1
fi
NCJv=`ls -1r ~/ncj/netcdfAll*.jar | head -n 1`
#NCJv="${HOME}/ncj/netcdfAll-4.5.jar"
echo "Found $NCJv"

echo "Looking for java"
which java
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find java in your path rc=$rc"
  exit 1
fi
JAVA=`which java`
#JAVA="/usr/bin/java"
echo "Found ${JAVA}"

echo "------------------------------------------------------------"
echo "running grib2nc.sh $1"
echo `date`
echo "------------------------------------------------------------"

FILEROOT=$1
GribFile="${FILEROOT}"
NCFile="${FILEROOT}.nc"

${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset \
     -in ${GribFile} -out ${NCFile} -IsLargeFile

rm ${FILEROOT}.ncx2 ${FILEROOT}.gbx9

