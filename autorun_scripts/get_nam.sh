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

# Shell script that downloads nam data files (091, 181, 196) for the date supplied
# on the command line.
# This script is called from autorun_nam.sh and takes three command-line arguments
#   get_ecmwf.sh NAM YYYYMMDD HR

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

NAM=$1
yearmonthday=$2
FChour=$3
#SERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod"
SERVER="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod"

echo "------------------------------------------------------------"
echo "running get_nam.sh script for ${NAM} $yearmonthday ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

case ${NAM} in
 181)
  # Caribbean 0.108 degrees
  HourMax=36
  HourStep=3
  #        nam.t00z.hawaiinest.hiresf00.tm00.grib2
  FilePre="nam.t${FChour}z.afwaca"
  FilePost=".tm00.grib2"
  ;;
 196)
  # HI 2.5 km
  HourMax=36
  HourStep=1
  #        nam.t00z.hawaiinest.hiresf00.tm00.grib2
  FilePre="nam.t${FChour}z.hawaiinest.hiresf"
  FilePost=".tm00.grib2"
  ;;
 091)
  # AK 2.95 km
  HourMax=36
  HourStep=1
  #        nam.t06z.alaskanest.hiresf00.tm00.grib2
  FilePre="nam.t${FChour}z.alaskanest.hiresf"
  FilePost=".tm00.grib2"
  ;;
 *)
  echo "NAM product not recognized"
  echo "Valid values: 091, 196"
  exit
esac

NAMDATAHOME="${WINDROOT}/nam/${NAM}"
install -d ${NAMDATAHOME}
if [[ $? -ne 0 ]] ; then
   echo "Error:  Download directory ${NAMDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi

#name of directory containing current files
FC_day=${yearmonthday}_${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
cd $NAMDATAHOME
mkdir -p $FC_day
cd $FC_day

t=0
while [ "$t" -le ${HourMax} ]; do
  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  INFILE=${FilePre}${hour}${FilePost}
  fileURL=${SERVER}/nam.${yearmonthday}/$INFILE
  time wget ${fileURL}
  ${USGSROOT}/bin/gen_GRIB_index $INFILE
  ${USGSROOT}/bin/autorun_scripts/grib2nc.sh $INFILE
  t=$(($t+${HourStep}))
done

mkdir -p $NAMDATAHOME/latest
cd $NAMDATAHOME/latest
rm nam.*
ln -s ../$FC_day/* .

t1=`date`
echo "download start: $t0"
echo "download   end: $t1"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished get_nam.sh ${NAM} ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
