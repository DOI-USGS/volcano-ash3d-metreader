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

# Shell script that manages the probing of gfs data at volcano locations in volc_NOVAC.txt.
#
# Requires that the following be installed:
#  ncks
# 
# This script expects a command line argument indicating which forecast package to probe.
# The date of the forecast package is determined at run-time by the UTC time, but can be
# overridden by setting the variable 'yearmonthday'
#   probe_volc.sh gfs0p50 0   for the 0.5 degree 00 forecast package

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

VOLCFILE="${USGSROOT}/share/volc_NOVAC.txt"

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: probe_volc.sh product FCpackage"
  echo "       where Resolution = gfs1p00, gfs0p50, gfs0p25, or nam091"
  echo "             FCpackage  = 0, 6, 12, or 18"
  exit
fi

PROD=$1
FC=$2

case ${PROD} in
 gfs0p25)
  echo "GFS 0.25 degree"
  ;;
 gfs0p50)
  echo "GFS 0.50 degree"
  HourMax=198
  HourStep=3
  ;;
 gfs1p00)
  echo "GFS 1.00 degree"
  ;;
 nam091)
  echo "NAM 2.95 km"
  HourMax=34
  #HourMax=1
  HourStep=1
  ;;
 *)
  echo "Product not recognized"
  echo "Valid values: gfs0p25, gfs0p50, gfs1p00, nam092"
  exit
esac

case ${FC} in
 0)
  FChour="00"
  FChourR="0.0"
  ;;
 6)
  FChour="06"
  FChourR="6.0"
  ;;
 12)
  FChour="12"
  FChourR="12.0"
  ;;
 18)
  FChour="18"
  FChourR="18.0"
  ;;
 *)
  echo "Forecast package not recognized"
  echo "Valid values: 0, 6, 12, 18"
  exit
esac

yearmonthday=`date -u +%Y%m%d`
# Here you can over-ride the date if need be
#yearmonthday="20210128"

echo "------------------------------------------------------------"
echo "running probe_volc ${PROD} ${yearmonthday} ${FChour} script"
echo "------------------------------------------------------------"

SCRIPTDIR="${USGSROOT}/bin/autorun_scripts"

#####  This is where the output files will be written
SONDEDIR="${WINDROOT}/sonde"
SONDEDIR="/data/www/vsc-ash.wr.usgs.gov/sonde"
SONDEDIR="/var/www/html/sonde"
mkdir -p ${SONDEDIR}
# First, we copy the master list from /opt/USGS/share/ to the sonde directory
vfile=${SONDEDIR}/volcs.dat
# all further parsing of the volcano list will apply to this new file.
cp ${VOLCFILE} ${vfile}
nvolc=`wc -l ${vfile} | cut -d' ' -f1`
echo $nvolc

for (( iv=1;iv<=${nvolc};iv++ ))
do
  volc=`head -n ${iv} ${vfile} | tail -1 | cut -f1 -d':'`
  lon=`head -n ${iv} ${vfile} | tail -1 | cut -f2 -d':'`
  lat=`head -n ${iv} ${vfile} | tail -1 | cut -f3 -d':'`
  echo "$volc $lon $lat"

  mkdir -p ${SONDEDIR}/${volc}/${yearmonthday}
  cd ${SONDEDIR}/${volc}/${yearmonthday}
  #mkdir -p ${SONDEDIR}/${volc}/temp
  #cd ${SONDEDIR}/${volc}/temp
  
  rm -f probe.log
  touch probe.log
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
  
    case ${PROD} in
     gfs0p50)
      ARGS="1 0 $lon $lat T 4 1 2 3 5 4 20 2"
      WINDPATH=${WINDROOT}/gfs/gfs.${yearmonthday}${FChour}
      WINDFILE=${yearmonthday}${FChour}.f${hour}.nc
      CLOUD="ncks -C -d lon,${lon} -d lat,${lat} -H -Q -s '%f ' -v Total_cloud_cover_isobaric ${WINDFILE} > Cloud.dat"
      ;;
     nam091)
      hour=`echo "${hour}" | cut -c2-3`
      ARGS="1 1 $lon $lat F 4 1 2 3 5 4 13 3"
      WINDPATH=${WINDROOT}/nam/091/${yearmonthday}_${FChour}
      WINDFILE=nam.t00z.alaskanest.hiresf${hour}.tm00.loc.grib2
      CLOUD="touch Cloud.dat"
      ;;
     *)
      echo "Product not recognized"
      echo "Valid values: gfs0p25, gfs0p50, gfs1p00, nam092"
      exit
    esac
  
  
    ln -s ${WINDPATH}/${WINDFILE} .
    echo "${USGSROOT}/bin/probe_Met ${WINDFILE} ${ARGS}"
    ${USGSROOT}/bin/probe_Met ${WINDFILE} ${ARGS}
    HourOffset=`echo "${FChour} + ${t}"  | bc`
    NewYYYYMMDD=`date -d"${yearmonthday} +${HourOffset} hour" -u +%Y%m%d`
    Newhour=`date -d"${yearmonthday} +${HourOffset} hour" -u +%H`
    mkdir -p ${SONDEDIR}/${volc}/${NewYYYYMMDD}
    mv NWP_prof.dat ${SONDEDIR}/${volc}/${NewYYYYMMDD}/${volc}_${PROD}_phuvt_${NewYYYYMMDD}_${Newhour}.dat
    ncks -C -d lon,${lon} -d lat,${lat} -H -Q -s '%f ' -v Total_cloud_cover_isobaric ${WINDFILE} > Cloud.dat
    mv Cloud.dat ${SONDEDIR}/${volc}/${NewYYYYMMDD}/${volc}_${PROD}_cloud_${NewYYYYMMDD}_${Newhour}.dat
    rm ${WINDFILE}
    #echo "$t : Mapping ${WINDFILE} to ${volc}_gfs_phuvt_${NewYYYYMMDD}_${Newhour}.dat" >> probe.log
    t=$((t+${HourStep}))
  done
  
  cd ${SONDEDIR}/
done

ls -l */* > tree.txt

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished probe_volc script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

