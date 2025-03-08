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

yearmonthday=$1
FChour=$2
SERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod"
#SERVER="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod"

echo "------------------------------------------------------------"
echo "running get_nam091.sh script for $yearmonthday ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

HourMax=36
HourStep=1

NAMDATAHOME="${WINDROOT}/nam/091"

#name of directory containing current files
FC_day=${yearmonthday}_${FChour}

#******************************************************************************
#START EXECUTING

# Note: Since grid 91 files are so big (~800 Mb/hourly-file), we will follow the
#       instructions on
#       http://nomads.ncep.noaa.gov/txt_descriptions/fast_downloading_grib.shtml
#       and first download the idx file, process it, and grab only the grib
#       layers needed.
#       The following extra utilities are required:
# ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_inv.pl
# ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_grib.pl
#       These files will be downloaded if they are not on the path

#go to correct directory
mkdir -p ${NAMDATAHOME}/${FC_day}
cd ${NAMDATAHOME}/${FC_day}

# Make sure we have the needed perl scripts for processing the index files
which get_inv.pl > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
        echo "Can not fine file get_inv.pl"
        echo "Downloading a fresh copy"
        wget ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_inv.pl
        chmod 775 get_inv.pl
        my_get_inv=${NAMDATAHOME}/${FC_day}/get_inv.pl
else
        my_get_inv=get_inv.pl
fi
which get_grib.pl > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
        echo "Can not fine file get_grib.pl"
        echo "Downloading a fresh copy"
        wget ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_grib.pl
        chmod 775 get_grib.pl
        my_get_grib=${NAMDATAHOME}/${FC_day}/get_grib.pl
else
        my_get_grib=get_grib.pl
fi




####################################################################
#  Set up parameters describing the variables we intend to download.
####################################################################

# Variables as a function of isobaric2
#  Geopotential_height_isobaric
#  u-component_of_wind_isobaric
#  v-component_of_wind_isobaric
#  Vertical_velocity_pressure_isobaric
#  Temperature_isobaric
#  Relative_humidity_isobaric
#  Specific_humidity_isobaric
niso2=42
iso2=("10 mb" "20 mb" "30 mb" "50 mb" "75 mb" "100 mb" "125 mb" "150 mb" "175 mb" "200 mb" "225 mb" "250 mb" "275 mb" "300 mb" "325 mb" "350 mb" "375 mb" "400 mb" "425 mb" "450 mb" "475 mb" "500 mb" "525 mb" "550 mb" "575 mb" "600 mb" "625 mb" "650 mb" "675 mb" "700 mb" "725 mb" "750 mb" "775 mb" "800 mb" "825 mb" "850 mb" "875 mb" "900 mb" "925 mb" "950 mb" "975 mb" "1000 mb")
nvar_iso2=7
var_iso2=("HGT" "UGRD" "VGRD" "VVEL" "TMP" "RH" "SPFH")
# Variables as a function of isobaric3
#  Cloud_mixing_ratio_isobaric
#  Snow_mixing_ratio_isobaric
niso3=40
iso3=("30 mb" "50 mb" "75 mb" "100 mb" "125 mb" "150 mb" "175 mb" "200 mb" "225 mb" "250 mb" "275 mb" "300 mb" "325 mb" "350 mb" "375 mb" "400 mb" "425 mb" "450 mb" "475 mb" "500 mb" "525 mb" "550 mb" "575 mb" "600 mb" "625 mb" "650 mb" "675 mb" "700 mb" "725 mb" "750 mb" "775 mb" "800 mb" "825 mb" "850 mb" "875 mb" "900 mb" "925 mb" "950 mb" "975 mb" "1000 mb")
nvar_iso3=2
var_iso3=("CLWMR" "SNMR")
# Variables as a function of height_above_ground7
#  u-component_of_wind_height_above_ground
#  v-component_of_wind_height_above_ground
nhag7=2
hag7=("10 m above ground" "80 m above ground")
nvar_hag7=2
var_hag7=("UGRD" "VGRD")
# Variables as a function of depth_below_surface_layer
#  Volumetric_Soil_Moisture_Content_depth_below_surface_layer
ndbsl=4
dbsl=("0-0.1 m below ground" "0.1-0.4 m below ground" "0.4-1 m below ground" "1-2 m below ground")
nvar_dbsl=1
var_dbsl=("SOILW")
# Variables as a function of surface
#  Planetary_Boundary_Layer_Height_surface
#  Frictional_Velocity_surface
#  Snow_depth_surface
#  Surface_roughness_surface
#  Wind_speed_gust_surface
#  Categorical_Rain_surface
#  Categorical_Snow_surface
#  Categorical_Freezing_Rain_surface
#  Categorical_Ice_Pellets_surface
#  Precipitation_rate_surface
nsurf=1
surf=("surface")
nvar_surf=10
var_surf=("HPBL" "FRICV" "SNOD" "SFCR" "GUST" "CRAIN" "CSNOW" "CFRZR" "CICEP" "PRATE")
# Variables as a function of some special 2d variable
#  Pressure_cloud_base
#  Pressure_cloud_tops
#  Total_cloud_cover_entire_atmosphere
nmisc2d=3
misc2d=("cloud base" "cloud top" "entire atmosphere (considered as a single layer)")
nvar_misc2d=("PRES" "PRES" "TCDC")

####################################################################
#  Now start the loop over all time steps
####################################################################
t=0
while [ "$t" -le ${HourMax} ]; do
  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  # Set up file names and get the index file listing locations of
  # the grib records.
  INFILE=nam.t${FChour}z.alaskanest.hiresf${hour}.tm00.grib2
  INFILEx=${INFILE}.idx
  OUTFILE=nam.t${FChour}z.alaskanest.hiresf${hour}.tm00.loc.grib2
  URL=${SERVER}/nam.${yearmonthday}/${INFILE}
  echo "$t of ${HourMax} $URL"
  # This script needs to be current (12/2018 or later) to handle https
  ${my_get_inv} $URL.idx > my_inv

  # Get all variables that are a function of isobaric2
  rm -f iso2.grib2
  touch iso2.grib2
  for (( iv=0;iv<$nvar_iso2;iv++))
  do
    for (( id=0;id<$niso2;id++))
    do
      echo "${t}/${HourMax} isobaric2 ${iv}/${nvar_iso2} ${var_iso2[iv]} ${iso2[id]}"
      grep ":${var_iso2[iv]}:" my_inv | grep ":${iso2[id]}:" > rec.tmp
      cat rec.tmp | ${my_get_grib} $URL tmp.grib2 > /dev/null 2>&1
      cat iso2.grib2 tmp.grib2 >> iso2.grib2 2>/dev/null
    done
  done
  # Get all variables that are a function of isobaric3
  rm -f iso3.grib2
  touch iso3.grib2
  for (( iv=0;iv<$nvar_iso3;iv++))
  do
    for (( id=0;id<$niso3;id++))
    do
      echo "${t}/${HourMax} isobaric3 ${iv}/${nvar_iso3} ${var_iso3[iv]} ${iso3[id]}"
      grep ":${var_iso3[iv]}:" my_inv | grep ":${iso3[id]}:" > rec.tmp
      cat rec.tmp | ${my_get_grib} $URL tmp.grib2  > /dev/null 2>&1
      cat iso3.grib2 tmp.grib2 >> iso3.grib2 2>/dev/null
    done
  done
  # Get all variables that are function of height_above_ground7
  rm -f hag7.grib2
  touch hag7.grib2
  for (( iv=0;iv<$nvar_hag7;iv++))
  do
    for (( id=0;id<$nhag7;id++))
    do
      echo "${t}/${HourMax} ght_above_ground7 ${iv}/${nvar_hag7} ${var_hag7[iv]} ${hag7[id]}"
      grep ":${var_hag7[iv]}:" my_inv | grep ":${hag7[id]}:" > rec.tmp
      cat rec.tmp | ${my_get_grib} $URL tmp.grib2  > /dev/null 2>&1
      cat hag7.grib2 tmp.grib2 >> hag7.grib2 2>/dev/null
    done
  done
  # Get all variables that are function of depth_below_surface_layer
  rm -f dbsl.grib2
  touch dbsl.grib2
  for (( iv=0;iv<$nvar_dbsl;iv++))
  do
    for (( id=0;id<$ndbsl;id++))
    do
      echo "${t}/${HourMax} depth_below_surface_layer ${iv}/${nvar_dbsl} ${var_dbsl[iv]} ${dbsl[id]}"
      grep ":${var_dbsl[iv]}:" my_inv | grep ":${dbsl[id]}:" > rec.tmp
      cat rec.tmp | ${my_get_grib} $URL tmp.grib2  > /dev/null 2>&1
      cat dbsl.grib2 tmp.grib2 >> dbsl.grib2 2>/dev/null
    done
  done
  # Get all variables that are function of surface
  rm -f surf.grib2
  touch surf.grib2
  for (( iv=0;iv<$nvar_surf;iv++))
  do
    for (( id=0;id<$nsurf;id++))
    do
      echo "${t}/${HourMax} surface ${iv}/${nvar_surf} ${var_surf[iv]} ${surf[id]}"
      grep ":${var_surf[iv]}:" my_inv | grep ":${surf[id]}:" > rec.tmp
      cat rec.tmp | ${my_get_grib} $URL tmp.grib2  > /dev/null 2>&1
      cat surf.grib2 tmp.grib2 >> surf.grib2 2>/dev/null
    done
  done
  # Get all variables that are function of some special 2d variable
  rm -f misc2d.grib2
  touch misc2d.grib2
  for (( iv=0;iv<$nvar_misc2d;iv++))
  do
    echo "${t}/${HourMax} misc ${iv}/${nvar_misc2d} ${var_misc2d[iv]} ${misc2d[iv]}"
    grep ":${var_misc2d[iv]}:" my_inv | grep ":${misc2d[iv]}:" > rec.tmp
    cat rec.tmp | ${my_get_grib} $URL tmp.grib2  > /dev/null 2>&1
    cat misc2d.grib2 tmp.grib2 >> misc2d.grib2 2>/dev/null
  done

  # Now bundle all these grib2 files into a grib2 file for this timestep
  cat iso2.grib2 iso3.grib2 hag7.grib2 dbsl.grib2 surf.grib2 misc2d.grib2 > ${OUTFILE}
  ${USGSROOT}/bin/gen_GRIB_index ${OUTFILE}

  rm -f iso2.grib2 iso3.grib2 hag7.grib2 dbsl.grib2 surf.grib2 misc2d.grib2 rec.tmp tmp.grib2

  t=$(($t+${HourStep}))

done # iterate to the next forecast hour

mkdir -p $NAMDATAHOME/latest
cd $NAMDATAHOME/latest
rm -f nam.*
ln -s $NAMDATAHOME/$FC_day/nam.* .
#
#t1=`date`
#echo "download start: $t0"
#echo "download   end: $t1"
