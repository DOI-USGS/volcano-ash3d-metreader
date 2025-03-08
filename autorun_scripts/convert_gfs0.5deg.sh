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

# Script that converts the gfs grib files to netcdf using netcdf-java.
# This script is called from autorun_gfs0.5deg.sh and takes two command-line arguments
#   get_gfs0.5deg.sh YYYYMMDD HR

# Check environment variables WINDROOT and USGSROOT
#  WINDROOT = location where the downloaded windfiles will be placed.
#  USGSROOT = location where the MetReader tools and scripts were placed.
# Please edit these to suit your system or ensure WINDROOT/USGSROOT are set as environment
# variables in ${HOME}/.bash_profile or ${HOME}/.bashrc
if [ -z ${WINDROOT} ];then
 # default location
 WINDROOT="/data/WindFiles"
fi
if [ -z ${USGSROOT} ];then
 # default location
 USGSROOT="/opt/USGS"
fi

# Select the netcdf version to write
#NCv=3
NCv=4
rc=0
yearmonthday=$1
FChour=$2
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


GFS="0p50"
HourMax=198
HourStep=3

echo "------------------------------------------------------------"
echo "running convert_gfs0.5deg.sh ${yearmonthday} ${FChour}"
echo `date`
echo "------------------------------------------------------------"

rc=0
GFSDATAHOME="${WINDROOT}/gfs"
if [[ -d ${GFSDATAHOME} ]] ; then
   echo "Error:  Download directory ${GFSDATAHOME} does not exist"
   rc=$((rc + 1))
   exit $rc
fi

#name of directory containing current files
FC_day=gfs.${yearmonthday}${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
echo "going to ${GFSDATAHOME}/${FC_day}"
cd ${GFSDATAHOME}/${FC_day}

#Convert to NetCDF
t=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  gfsfile="gfs.t${FChour}z.pgrb2.${GFS}.f${hour}"
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  if test -r ${gfsfile}
  then
     echo "Converting ${gfsfile} to ${netcdffile}"
     echo "java -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${gfsfile} -out ${netcdffile} -IsLargeFile"
     ${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${gfsfile} -out ${netcdffile} -IsLargeFile
     if [[ $? -ne 0 ]]; then
          exit 1
     fi
     t=$((t+${HourStep}))
   else
     echo "error: ${gfsfile} does not exist."
     exit 1
   fi
done

#Make sure the netcdf files all exist
echo "making sure all netcdf files exist"
t=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  if test -r ${netcdffile}
  then
     echo "${netcdffile} exists"
     t=$((t+${HourStep}))
   else
     echo "error: ${netcdffile} does not exist."
     exit 1
   fi
done

echo "all done with windfiles"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished convert_gfs0.5deg.sh ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
