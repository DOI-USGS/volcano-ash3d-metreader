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

#  Script that get NCEP wind files for the year specified.  
#  This script is run at the beginning of the new year, to complete the NCEP record of
#  the previous year's wind files.  But it can also be run to get the wind files for
#  any previous year.  This script is invoked by typing:
#    ./get_NCEP_50YearReanalysis.sh $year
#  where $year is the year you want to get.

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

echo "------------------------------------------------------------"
echo "running get_NCEP_50YearReanalysis.sh for year $1"
echo `date`
echo "------------------------------------------------------------"

starttime=`date`                  #record when we're starting the download
echo "starting get_NCEP_getyear.sh at $starttime"

rc=0
NCEPDATAHOME="${WINDROOT}/NCEP"
#if [[ -d ${NCEPDATAHOME} ]] ; then
#   echo "Error:  Download directory ${NCEPDATAHOME} does not exist"
#   rc=$((rc + 1))
#   exit $rc
#fi
echo "NCEPDATAHOME=$NCEPDATAHOME"

#get year
if [[ $1 -eq '' ]]; then
   echo "Error: you need to provide your desired year as an argument"
   exit 1
else
   y=$1
   echo "year=$y"
fi

#Make sure the destination directory exists
#echo "making sure the directory for year ${y} exists"
#if [ ! -r "${NCEPDATAHOME}/${y}" ]         
#then
#   echo "Error: Directory ${NCEPDATAHOME}/${y} does not exist"
#   exit 1
#else
#   echo "Good.  It does."
#   echo "Copying contents of ${NCEPDATAHOME}/${y} to ${NCEPDATAHOME}/backup"
#   cp -v ${NCEPDATAHOME}/${y}/* ${NCEPDATAHOME}/backup
#   echo "all done copying files"
#fi

#if test -r ${NCEPDATAHOME}/dbuffer
#then
#    echo "moving to ${NCEPDATAHOME}/dbuffer"
    cd ${NCEPDATAHOME}/dbuffer
#  else
#    echo "error: ${NCEPDATAHOME}/dbuffer does not exist."
#    exit 1
#fi

var=( air hgt omega shum uwnd vwnd )
for (( i=0;i<=5;i++))
do
 echo "wget http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep/${var[i]}.${y}.nc"
 wget http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis/pressure/${var[i]}.${y}.nc
 mv ${NCEPDATAHOME}/dbuffer/${var[i]}.${y}.nc ${NCEPDATAHOME}/${y}/${var[i]}.${y}.nc       #overwrite older file
done

echo "Finished downloading."
endtime=`date`

echo "download started at $starttime"
echo "download ended at $endtime"
echo "Checking files for valid values and time steps"
${USGSROOT}/bin/MetCheck 25 2 ${NCEPDATAHOME} ${y}
if [[ $? -eq 0 ]]; then
   cat MetCheck_log.txt >> ${NCEPDATAHOME}/${y}/${validlist}
fi
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished get_NCEP_50YearReanalysis.sh"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
