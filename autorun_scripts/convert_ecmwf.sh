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

# Script that converts the ecmwf grib files to netcdf using netcdf-java.
# This script is called from autorun_ecmwf.sh and takes three command-line arguments
#   convert_ecmwf.sh RES YYYYMMDD HR

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
if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: convert_ecmwf.sh Resolution YYYYMMDD FCpackage"
  echo "       where Resolution = 1p00, 0p50, or 0p25"
  echo "             YYYYMMDD   = year, month, day"
  echo "             FCpackage  = 00, 06, 12, or 18"
  exit
fi

RES=$1
yearmonthday=$2
FChour=$3

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

# Checking input parameters
case ${RES} in
 0p25)
  echo "RES 0.25 degree"
  ;;
 *)
  echo "RES product not recognized"
  echo "Valid values: 0p25 "
  exit
esac

echo "------------------------------------------------------------"
echo "running convert_ecmwf.sh ${RES} ${yearmonthday} ${FChour}"
echo `date`
echo "------------------------------------------------------------"

case ${RES} in
 0p25)
  # ECMWF 0.25 degree
  HourMax=99
  HourStep=3
  #        20250211000000-0h-oper-fc.grib2
  FilePre="${yearmonthday}${FChour}0000-"
  FilePost="h-oper-fc.grib2"
  iwf=22
  ;;
 *)
  echo "ECMWF product not recognized"
  echo "Valid values: 0p25"
  exit
esac

rc=0
validlist="valid_files.txt"
ECMWFDATAHOME="${WINDROOT}/ecmwf"
#if [[ -d ${ECMWFDATAHOME} ]] ; then
#   echo "Error:  Download directory ${ECMWFDATAHOME} does not exist"
#   rc=$((rc + 1))
#   exit $rc
#fi

#name of directory containing current files
FC_day=ecmwf.${yearmonthday}${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
echo "going to ${ECMWFDATAHOME}/${FC_day}"
cd ${ECMWFDATAHOME}/${FC_day}

#Convert to NetCDF
t=0
rm -f ${ECMWFDATAHOME}/${FC_day}/${validlist}
touch ${ECMWFDATAHOME}/${FC_day}/${validlist}
vcount=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  ecmwffile=${FilePre}${hour}
  ncmlfile="ecmwf.t${FChour}z.f${hour}.ncml"
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  netcdf4file="${yearmonthday}${FChour}.f${hour}.nc4"
  if test -r ${ecmwffile}
  then
     #echo "making ${ncmlfile}"
     #${USGSROOT}/bin/makegfsncml ${ecmwffile} ${ncmlfile}
     #if [[ $? -ne 0 ]]; then
     #     exit 1
     #fi
     echo "Converting ${ecmwffile} to ${netcdffile}"
     # Note: nc4 is significantly smaller, but the direct conversion to nc4 from netcdf-java via the flag -netcdf4
     #       results in incompatible files.  Here, we use netcdf-java to product a temporary file, then convert
     #       to nc4 with nccopy
     echo "java -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${ecmwffile} -out ${netcdffile} -IsLargeFile"
     if [ $NCv -eq 4 ]
     then
       ${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${ecmwffile} -out tmp.nc -IsLargeFile
       nccopy -k 4 -d 5 tmp.nc ${netcdffile}
       rm tmp.nc
     else
       ${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${ecmwffile} -out ${netcdffile} -IsLargeFile
     fi
     if [[ $? -ne 0 ]]; then
          exit 1
     fi
     # Check converted file for valid values
     echo "checking ${netcdffile} for corrupt values"
     ${USGSROOT}/bin/MetCheck ${iwf} 2 ${netcdffile}
     if [[ $? -eq 0 ]]; then
       cat MetCheck_log.txt >> ${ECMWFDATAHOME}/${FC_day}/${validlist}
       vcount=$((vcount+1))
     fi
   else
     echo "Warning: ${ecmwffile} does not exist. Skipping."
   fi
   t=$((t+${HourStep}))
done

#Make sure the netcdf files all exist
echo "making sure all netcdf files exist"
t=0        # time index
numfiles=0 # file index
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
     numfiles=$((numfiles+1))
   else
     echo "error: ${netcdffile} does not exist."
     exit 1
   fi
done

#set soft links in "latest" directory
echo "creating soft links in latest directory"
mkdir -p ${ECMWFDATAHOME}/latest
echo "rm ${ECMWFDATAHOME}/latest/*"
rm ${ECMWFDATAHOME}/latest/*
echo "rm ${ECMWFDATAHOME}/ecmwflist.txt"
rm ${ECMWFDATAHOME}/ecmwflist.txt
echo "4" > ${ECMWFDATAHOME}/ecmwflist.txt
echo "$numfiles" >> ${ECMWFDATAHOME}/ecmwflist.txt
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
  linkfile="latest.f${hour}.nc"
  echo "creating soft link for ${netcdffile}"
  ln -s ${ECMWFDATAHOME}/${FC_day}/${netcdffile} ${ECMWFDATAHOME}/latest/${linkfile}
  echo "latest/${linkfile}" >> ${ECMWFDATAHOME}/ecmwflist.txt
  t=$((t+3))
done

echo "removing *.ncml, *.ncx2, and *.gbx9 files"
echo "rm ${ECMWFDATAHOME}/${FC_day}/*.ncml ${ECMWFDATAHOME}/${FC_day}/*.gbx9"
rm -f ${ECMWFDATAHOME}/${FC_day}/*.ncml ${ECMWFDATAHOME}/${FC_day}/*.gbx9 ${ECMWFDATAHOME}/${FC_day}/*.ncx2

echo "writing last_downloaded.txt"
echo ${yearmonthday}${FChour} > ${ECMWFDATAHOME}/last_downloaded.txt

echo "all done with windfiles"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished convert_ecmwf.sh ${RES} ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
