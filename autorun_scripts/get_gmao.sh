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

# Shell script that downloads gmao data files for the date supplied on the command line.
# This script is takes four command-line arguments
#   get_gmao.sh YYYY MM DD HH

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

if [ -z $4 ]
then
  echo "Error: Insufficient command-line arguments"
  echo "Usage:  get_gmao.sh YYYY MM DD HH"
  echo "        where HH = 0, 6, 12, or 18"
  exit 1
else
  YYYY=$1
  MM=$2
  DD=$3
  HH=$4
  filename="GEOS.fp.fcst.inst3_3d_asm_Np.${YYYY}${MM}${DD}_00+${YYYY}${MM}${DD}_${HH}00.V01.nc4"
  echo "Downloading file: ${filename}"
fi

NASADATAHOME="${WINDROOT}/NASA/GEOS"
#name of directory containing current files
FC_day=${NASADATAHOME}/${yearmonthday}

#go to correct directory
mkdir -p $FC_day
cd $FC_day

wget --password="" ftp://gmao_ops@ftp.nccs.nasa.gov/fp/forecast/Y${YYYY}/M${MM}/D${DD}/H00/GEOS.fp.fcst.inst3_3d_asm_Np.${YYYY}${MM}${DD}_00+${YYYY}${MM}${DD}_${HH}00.V01.nc4
