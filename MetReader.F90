!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This file is a component of the volcanic ash transport and dispersion model Ash3d,
!      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
!      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).
!
!      The model and its source code are products of the U.S. Federal Government and therefore
!      bear no copyright.  They may be copied, redistributed and freely incorporated 
!      into derivative products.  However as a matter of scientific courtesy we ask that
!      you credit the authors and cite published documentation of this model (below) when
!      publishing or distributing derivative products.
!
!      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
!         volume, conservative numerical model for ash transport and tephra deposition,
!         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 
!
!      Although this program has been used by the USGS, no warranty, expressed or
!      implied, is made by the USGS or the United States Government as to the accuracy
!      and functioning  of the program and related program material nor shall the fact of
!      distribution constitute  any such warranty, and no responsibility is assumed by
!      the USGS in connection therewith.
!
!      We make no guarantees, expressed or implied, as to the usefulness of the software
!      and its documentation for any purpose.  We assume no responsibility to provide
!      technical support to users of this software.
!
!      
!
!      contains:
!        subroutine MR_Reset_Memory
!        subroutine MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
!        subroutine MR_Read_Met_DimVars(iy)
!        subroutine MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)
!        subroutine MR_Initialize_Met_Grids(nx,ny,nz,dumx_sp,dumy_sp,dumz_sp,periodic)
!        subroutine MR_Set_SigmaAlt_Scaling(nz,dums_sp)
!        subroutine MR_Set_Met_Times(eStartHour,Duration)
!        subroutine MR_Read_HGT_arrays(istep,reset_first_time)
!        subroutine MR_Read_3d_MetP_Variable(ivar,istep)
!        subroutine MR_Read_3d_MetH_Variable(ivar,istep)
!        subroutine MR_Read_3d_Met_Variable_to_CompP(ivar,istep,IsNext)
!        subroutine MR_Read_3d_Met_Variable_to_CompH(ivar,istep,IsNext)
!        subroutine MR_Read_2d_Met_Variable(ivar,istep)
!        subroutine MR_Read_2d_Met_Variable_to_CompGrid(ivar,istep)
!        subroutine MR_Rotate_UV_GR2ER_Met(istep,SetComp)
!        subroutine MR_Rotate_UV_ER2GR_Comp(istep)
!        subroutine MR_Regrid_MetP_to_CompH(istep)
!        subroutine MR_Regrid_MetP_to_MetH(istep)
!        subroutine MR_Regrid_Met2d_to_Comp2d()
!        subroutine MR_Set_LL_mapping()
!        subroutine MR_DelMetP_Dx()
!        subroutine MR_DelMetP_Dy()
!        subroutine MR_QC_3dvar(ivar,nx_max,ny_max,nz1_max,z_array_sp,nz2_max,&
!                               dum_array_sp,fill_val_sp,bc_low_sp, bc_high_sp)
!        subroutine MR_Check_Prerequsites(test_allocate,test_dimvars,test_compproj,&
!                                         test_initmetgrid,test_setmettimes)
!        subroutine MR_FileIO_Error_Handler(ios,linebuffer080)
!        function MR_Temp_US_StdAtm(zin)
!        function MR_Z_US_StdAtm(pin)
!        function MR_Pres_US_StdAtm(zin)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module MetReader

      ! This module requires Fortran 2003 or later
      use ieee_arithmetic, only : &
         ieee_is_nan

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public MR_Reset_Memory,MR_Allocate_FullMetFileList,MR_Set_CompProjection,&
             MR_Read_Met_DimVars,MR_Set_Met_Times,MR_Initialize_Met_Grids,MR_Read_HGT_arrays,&
             MR_Read_3d_Met_Variable_to_CompP,MR_Read_3d_Met_Variable_to_CompH,&
             MR_Set_SigmaAlt_Scaling,MR_Read_3d_MetP_Variable,MR_Read_3d_MetH_Variable,&
             MR_Read_2d_Met_Variable,MR_Read_2d_Met_Variable_to_CompGrid,&
             MR_Rotate_UV_GR2ER_Met,MR_Rotate_UV_ER2GR_Comp,&
             MR_Regrid_MetP_to_CompH,MR_Regrid_MetP_to_MetH,MR_Regrid_Met2d_to_Comp2d,&
             MR_DelMetP_Dx,MR_DelMetP_Dy,MR_Set_LL_mapping,&
             MR_Temp_US_StdAtm,MR_Z_US_StdAtm,MR_Pres_US_StdAtm,&
             MR_QC_3dvar,MR_FileIO_Error_Handler

      integer, parameter,private :: sp = selected_real_kind( 6,   37) ! single precision
      integer, parameter,private :: dp = selected_real_kind(15,  307) ! double precision
      integer, parameter,private :: qp = selected_real_kind(33, 4931) ! quad precision

        ! Publicly available variables
#include "MR_version.h"
      ! This file (MR_version.h) is the output of the simple script
!get_version.sh 
!  echo -n "      character(len=40),parameter,public :: MR_GitComID ='" > MR_version.h
!  git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> MR_version.h
!  echo -n "'" >> MR_version.h
!    which sets the variable containing the git commit ID : MR_GitComID
#ifdef USENETCDF
      integer          ,public :: NCv_datafile
      character(len=80),public :: NCv_lib
#endif
      integer      ,parameter,public :: MR_MAXVARS   = 50 ! Maximum number of variables in fixed arrays

      real(kind=sp),parameter,public :: MR_EPS_SMALL = 1.0e-7_sp  ! Small number

      real(kind=sp),parameter,public :: MR_RAD_EARTH = 6371.229_sp ! Radius of Earth in km
      real(kind=sp),parameter,public :: MR_DEG2RAD   = 1.7453292519943295e-2_sp
      real(kind=sp),parameter        :: MR_MIN_DZ    = 1.0e-4_sp   ! used in reassigning z for low-level negative GPH

      ! These output levels and stream IDs can be changed by the calling program
      ! default verbosity is 3
      ! This will write everything from verbosity_log to verbosity_error to both stdout and stdlog
      integer,dimension(2),public :: VB = (/3,3/) ! Verbosity level for stdout and logfile, respectively
      integer,public              :: MR_VERB   = 3
      integer,public              :: stdin     = 5

      integer,public :: MR_nio = 1   ! Number of output streams to use: 1 for just stdout,stderr or 2 to include a logfile
      integer,dimension(2),public :: outlog = (/6,9/) ! Assumes stdin=6, logfile=9
      integer,dimension(2),public :: errlog = (/0,9/) ! Assumes sterr=0, logfile=9

      integer,parameter,public :: verbosity_debug2       = 1  ! Additional debugging information only written to stdout
      integer,parameter,public :: verbosity_debug1       = 2  ! Debugging information only written to stdout
      integer,parameter,public :: verbosity_log          = 3  ! Time step information (this is the limit for writing to logfile)
      integer,parameter,public :: verbosity_info         = 4  ! Additional information on run set up and shutdown
      integer,parameter,public :: verbosity_statistics   = 5  ! Details on health of run (timing, mass conservation)
      integer,parameter,public :: verbosity_production   = 6  ! Major program flow info
      integer,parameter,public :: verbosity_essential    = 7  ! Only start up and shutdown messages
      integer,parameter,public :: verbosity_error        = 8  ! No logging to stdout, only stderr (and logfile)
      integer,parameter,public :: verbosity_silent       = 9  ! No logging to stdout,stderr. Logfile written as normal
      integer,parameter,public :: verbosity_dark         = 10 ! No logging to stdout,stderr or logfile

      logical          ,public :: MR_useCompGrid         = .true. ! Reset this to .false. if you only need the Met grid
      logical          ,public :: MR_useCompTime         = .true. ! Reset this to .false. if you only need the time of the file
      logical          ,public :: MR_useCompH            = .true.
      logical          ,public :: MR_useCompP            = .true.
      logical          ,public :: MR_useCompS            = .false.

      integer,public :: MR_iwind       !     MR_IWIND specifies the type of wind input to the model:
                             !   MR_IWIND=1 if a 1-D wind sounding is use, 
                             !           =2 if a 3-D grid is read from a ASCII file.
                             !           =3 if a single, multistep 3-D file is used
                             !           =4 if multiple 3-D NetCDF files are used
                             !           =5 if multiple file with multiple steps are used
      integer,public :: MR_iwindformat !      MR_iwindformat specifies the NWP product
                                       !  0 Custom format based on template
                                       !  1 ASCII profile
                                       !  2 Radiosonde data
                                       !  3 NARR3D 32 km North America files (32 km) :: ds608.0
                                       !  4 NAM Regional North America 221 (32 km)
                                       !  5 NAM AK 216  (45 km)
                                       !  6 NAM Regional 104 (90 km)
                                       !  7 NAM CONUS 212 (40 km)
                                       !  8 NAM CONUS 218 (12 km)
                                       !  9 NAM CONUS 227 (5.08 km)
                                       ! 10 NAM AK 242 (11.25 km)
                                       ! 11 NAM HI 196 (2.5 km)
                                       ! 12 NAM AK 198 (5.953 km)
                                       ! 13 NAM AK 91 (2.976 km)
                                       ! 14 NAM CONUS 1227 (3.0 km)
                                       ! 20 GFS 0.5
                                       ! 21 GFS 1.0
                                       ! 22 GFS 0.25
                                       ! 23 NCEP / DOE reanalysis 2.5 degree files  :: ds091.0
                                       ! 24 NASA-MERRA-2 reanalysis 0.625/0.5 degree files
                                       ! 25 NCEP/NCAR reanalysis 2.5 degree files   :: ds090.0  iwind=4 or 5
                                       ! 26 JRA-55                                  :: ds628.0  iwind=5
                                       ! 27 NOAA-CIRES reanalysis 2.0 degree files  :: ds131.2,3  iwind=5
                                       ! 28 ECMWF Interim Reanalysis (ERA-Interim)  :: ds627.0  requires catted GRIB files
                                       ! 29 ECMWF ERA5                              :: ds633.0  iwind=5
                                       ! 30 ECMWF 20-Century (ERA-20C)              :: ds626.0  iwind=5
                                       ! 31 NAM Caribbean 181 (0.108 deg)
                                       ! 32 Air Force Weather Agency subcenter = 0
                                       ! 33 CCSM3.0 Community Atmosphere Model (CAM)
                                       ! 34 ECMWF 0.25 degree forecast
                                       ! 40 NASA-GEOS Cp
                                       ! 41 NASA-GEOS Np
                                       ! 50 WRF - output
      logical,public :: MR_Use_RDA     = .false.  ! If reanalysis products were acquired via RDA, this flag can be changed
      integer,public :: MR_RDAcode     = -1       ! This allows the calling program to specify older versions of data
      integer,public :: MR_iversion    = -1       ! version of the product; sometimes the grid used changes with newer
                                                  ! versions. -1 indicates latest version, but can be overridden here
      integer,public :: MR_iGridCode   !   MR_iGridCode specifies the NCEP grid described in:
                                       !   http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
      integer,public :: MR_idataFormat !   Specifies the data model used
                                    !    =1 ASCII
                                    !    =2 NetCDF
                                    !    =3 GRIB1 or GRIB2

        ! These variables describe the full list of windfiles read
      integer                                  ,public :: MR_iwindfiles           ! number of files provided
      character(len=130)                       ,public :: MR_iw5_root
      character(len=42)                        ,public :: MR_iw5_prefix
      character(len=24)                        ,public :: MR_iw5_suffix1  ! e.g. 1912060100_1912063021.nc
      character(len=24)                        ,public :: MR_iw5_suffix2
      real(kind=dp)                            ,public :: MR_iw5_hours_per_file

#if USEPOINTERS      
!      character(len=130),dimension(:)  ,pointer,public :: fc_windfilename
      character(len=130),dimension(:)  ,pointer,public :: MR_windfiles                 => null() ! name of file
      real(kind=dp)     ,dimension(:)  ,pointer,public :: MR_windfile_starthour        => null()
      real(kind=dp)     ,dimension(:,:),pointer,public :: MR_windfile_stephour         => null()
      integer           ,dimension(:)  ,pointer,public :: MR_windfiles_nt_fullmet      => null() ! number of steps in files
      logical           ,dimension(:)  ,pointer,public :: MR_windfiles_Have_GRIB_index => null()
      character(len=130),dimension(:)  ,pointer,public :: MR_windfiles_GRIB_index      => null() ! name of GRIB index file`

#else
!      character(len=130), allocatable,dimension(:)  ,public :: fc_windfilename
      character(len=130),dimension(:)  ,allocatable,public :: MR_windfiles            ! name of file
      real(kind=dp)     ,dimension(:)  ,allocatable,public :: MR_windfile_starthour   ! start hour of the file
      real(kind=dp)     ,dimension(:,:),allocatable,public :: MR_windfile_stephour    ! offset hours of step
      integer           ,dimension(:)  ,allocatable,public :: MR_windfiles_nt_fullmet ! number of steps in files
      logical           ,dimension(:)  ,allocatable,public :: MR_windfiles_Have_GRIB_index
      character(len=130),dimension(:)  ,allocatable,public :: MR_windfiles_GRIB_index ! name of GRIB index file
#endif
      character(len=80)                           ,public :: MR_iwf_template         ! name of the template file
      logical                                     ,public :: MR_Reannalysis     = .false.
      logical                                     ,public :: MR_Save_Velocities = .false.

        ! These variables are a list of the same data above, but specific to the simulation
        ! duration
      integer                                     ,public :: MR_MetSteps_Total
      integer                                     ,public :: MR_iMetStep_Now
#ifdef USEPOINTERS
      character(len=130),dimension(:),pointer,public :: MR_MetStep_File                => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_findex              => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_tindex              => null()
      real(kind=dp)     ,dimension(:),pointer,public :: MR_MetStep_Hour_since_baseyear => null()
      real(kind=dp)     ,dimension(:),pointer,public :: MR_MetStep_Interval            => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_year                => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_month               => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_day                 => null()
      integer           ,dimension(:),pointer,public :: MR_MetStep_DOY                 => null()
      real(kind=dp)     ,dimension(:),pointer,public :: MR_MetStep_Hour_Of_Day         => null()
      integer           ,dimension(:),pointer,public :: MR_iwind5_year                 => null()
#else
      character(len=130),dimension(:),allocatable,public :: MR_MetStep_File
      integer           ,dimension(:),allocatable,public :: MR_MetStep_findex
      integer           ,dimension(:),allocatable,public :: MR_MetStep_tindex
      real(kind=dp)     ,dimension(:),allocatable,public :: MR_MetStep_Hour_since_baseyear
      real(kind=dp)     ,dimension(:),allocatable,public :: MR_MetStep_Interval
      integer           ,dimension(:),allocatable,public :: MR_MetStep_year
      integer           ,dimension(:),allocatable,public :: MR_MetStep_month
      integer           ,dimension(:),allocatable,public :: MR_MetStep_day
      integer           ,dimension(:),allocatable,public :: MR_MetStep_DOY
      real(kind=dp)     ,dimension(:),allocatable,public :: MR_MetStep_Hour_Of_Day
      integer           ,dimension(:),allocatable,public :: MR_iwind5_year
#endif

      !    Native grid of Met file using Pressure as vertical coordinate
#ifdef USEPOINTERS
      integer      ,dimension(:,:)  ,pointer, public :: MR_dum2d_met_int   => null() ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,pointer, public :: MR_dum2d_met       => null() ! Used for surface variables
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_metP      => null()
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d2_metP     => null()
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_geoH_metP_last  => null() ! These are needed for compH interpolation
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_geoH_metP_next  => null() 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vx_metP_last    => null() !These might need to be stored to avoid a 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vx_metP_next    => null() !second reading
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vy_metP_last    => null() 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vy_metP_next    => null() 
#else
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_met_int   ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_met       ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d2_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_last  ! These are needed for compH interpolation
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_next
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_last  !These might need to be stored to avoid a 
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_next  !second reading
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_last
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_next
#endif

      real(kind=sp), public :: MR_Max_geoH_metP_predicted   ! Highest expected height in Met file, based on pressure
      real(kind=sp)         :: Max_geoH_metP_last
      real(kind=sp)         :: Max_geoH_metP_next
      real(kind=sp)         :: Max_geoH_metP
      real(kind=sp)         :: Suppl_H
      !Parameter iHeightHandler specifies what to do if the maximum height
      !of the simulation region exceeds the maximum height in the wind files.
      !If iHeightHandler = 1, stop the program if the plume height exceeds mesoscale height
      !                    2, wind direction at levels above the highest node 
      !                       equal that of the highest node.  Temperatures in the
      !                       upper nodes don't change between 11 and 20 km; above
      !                       20 km they increase by 2 C/km, as in the Standard
      !                       atmosphere.  A warning is written to the log file.
      integer       ,public :: MR_iHeightHandler = 2

        ! The following are used for sonde data
      integer                                   ,public :: num_RadSnd_Stat
#ifdef USEPOINTERS
      character(len=4 ),dimension(:)      ,pointer,public :: MR_Snd_cd          => null() ! Station code
      integer          ,dimension(:)      ,pointer,public :: MR_Snd_id          => null() ! WMO station ID
      real(kind=sp)    ,dimension(:)      ,pointer,public :: MR_Snd_lt          => null() ! Station latitude
      real(kind=sp)    ,dimension(:)      ,pointer,public :: MR_Snd_ln          => null() ! Station longitude
      real(kind=sp)    ,dimension(:)      ,pointer,public :: MR_Snd_el          => null() ! Station elevation
      character(len=25),dimension(:)      ,pointer,public :: MR_Snd_lnm         => null() ! Station long name
      character(len=2 ),dimension(:)      ,pointer,public :: MR_Snd_st          => null() ! Station state
      character(len=2 ),dimension(:)      ,pointer,public :: MR_Snd_ct          => null() ! Station country
      integer          ,dimension(:)      ,pointer,public :: Snd_idx            => null() ! Index of radiosonde in master list
      real(kind=sp)    ,dimension(:,:,:,:),pointer,public :: MR_SndVars_metP    => null() ! (MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,300)
      integer          ,dimension(:,:)    ,pointer,public :: MR_Snd_np_fullmet  => null() ! Number of pressure values for each location/time
      integer          ,dimension(:)      ,pointer,public :: MR_SndVarsID       => null() ! Lists which vars are in which columns of MR_SndVars_metP
      real(kind=sp)    ,dimension(:,:,:)  ,pointer,public :: MR_Snd2Comp_map_wgt=> null() ! weights of nearby sondes for every comp point
      integer          ,dimension(:,:,:)  ,pointer,public :: MR_Snd2Comp_map_idx=> null() ! sonde index of weights
#else
      character(len=4 ),dimension(:)      ,allocatable,public :: MR_Snd_cd         ! Station code
      integer          ,dimension(:)      ,allocatable,public :: MR_Snd_id         ! WMO station ID
      real(kind=sp)    ,dimension(:)      ,allocatable,public :: MR_Snd_lt         ! Station latitude
      real(kind=sp)    ,dimension(:)      ,allocatable,public :: MR_Snd_ln         ! Station longitude
      real(kind=sp)    ,dimension(:)      ,allocatable,public :: MR_Snd_el         ! Station elevation
      character(len=25),dimension(:)      ,allocatable,public :: MR_Snd_lnm        ! Station long name
      character(len=2 ),dimension(:)      ,allocatable,public :: MR_Snd_st         ! Station state
      character(len=2 ),dimension(:)      ,allocatable,public :: MR_Snd_ct         ! Station country
      integer          ,dimension(:)      ,allocatable,public :: Snd_idx           ! Index of radiosonde in master list
      real(kind=sp)    ,dimension(:,:,:,:),allocatable,public :: MR_SndVars_metP   ! (MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,300)
      integer          ,dimension(:,:)    ,allocatable,public :: MR_Snd_np_fullmet ! Number of pressure values for each location/time
      integer          ,dimension(:)      ,allocatable,public :: MR_SndVarsID      ! Lists which vars are in which columns of MR_SndVars_metP
      real(kind=sp)    ,dimension(:,:,:)  ,allocatable,public :: MR_Snd2Comp_map_wgt ! weights of nearby sondes for every comp point
      integer          ,dimension(:,:,:)  ,allocatable,public :: MR_Snd2Comp_map_idx ! sonde index of weights
#endif

      integer,public :: MR_nSnd_Locs      = 1        ! Number of Sonde locations
      integer,public :: MR_Snd_nt_fullmet = 1        ! Number of times at the Sonde locations
      integer,public :: MR_Snd_nvars      = 5        ! Number of Sonde variables (P,H,U,V,T,+user-specified)
      logical,public :: Snd_Have_PT       = .false.  ! 3-col ASCII input do not have pres and temp
      logical,public :: Snd_Have_Coord    = .false.  ! If the 1-d data have the optional projection params
                                           ! then it will be used, otherwise vel will be relative
                                           ! to comp grid.
      integer                                     , public :: MR_nstat = 8      ! number of stations to consider (0 for all)
      real(kind=sp)                               , public :: MR_pexp  = 4.0_sp ! exponent for inverse distance calculation

      !    Native grid of Met file using Height as vertical coordinate
      !    (resampled onto z-gridpoints of computational grid)
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_metH => null()
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_u_ER_metP  => null()! For the cases where Met is proj and comp
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_v_ER_metP  => null()!  different we need to rotate so these
                                                                               !  store Earth-Relative velocities on MetP
#else
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metH
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_u_ER_metP ! For the cases where Met is proj and comp
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_v_ER_metP !  different we need to rotate so these
                                                                         !  store Earth-Relative velocities on MetP
#endif

      !    Computations grid
#ifdef USEPOINTERS
      integer      ,dimension(:,:)  ,pointer, public :: MR_dum2d_comp_int => null() ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,pointer, public :: MR_dum2d_comp     => null() ! Used for surface variables
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compP    => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compP_2  => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compH    => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compH_2  => null() ! Used only when a vector field
                                                                               ! rotation is needed
#else
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_comp_int  ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_comp      ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compP     !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compP_2   !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH     !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH_2   ! Used only when a vector field
                                                                               ! rotation is needed
#endif

      integer,public       :: nt_fullmet    ! length of t of file
      integer,public       :: nx_fullmet    ! length of x or lon of full met grid
      integer,public       :: ny_fullmet    ! length of y or lat of full met grid
      integer,public       :: np_fullmet    ! length of pressure of full met grid for H,U,V
      integer              :: np_fullmet_pad = 1 ! We might need to pad the top of the pressure grid.
      integer,public       :: neta_fullmet  ! Only used by WRF
      real(kind=sp),public :: Pressure_Conv_Fac = 100.0_sp  ! factor for converting from Met file units to Pa

#ifdef USEPOINTERS
      real(kind=sp),dimension(:),   pointer,public :: x_fullmet_sp    => null() ! x-coordinates of full met grid
      real(kind=sp),dimension(:),   pointer,public :: y_fullmet_sp    => null() ! y-coordinates of full met grid
      real(kind=sp),dimension(:),   pointer,public :: p_fullmet_sp    => null() ! z-coordinates of full met grid for H
      real(kind=sp),dimension(:,:), pointer,public :: levs_fullmet_sp => null() ! This hold each of the numbered level coordinates:
                                                                         !    i.e. isobaric, isobaric1, isobaric2, but also
                                                                         !    height_above_ground, depth_below_surface_layer, etc.
                                                                         !    p_fullmet_sp,p_fullmet_[Vz,RH]sp are copies of one
                                                                         !    of the slices
      integer,dimension(:), pointer,public :: nlevs_fullmet                     ! length of z coordinate
      integer,dimension(:), pointer,public :: levs_code                         ! code indicating how to map to the GPH grid
                                                                         !   0 = no mapping (not a pressure coordinate)
                                                                         !   1 = one-to-one mapping (U,V)
                                                                         !   2 = upper truncation (missing upper levels)
                                                                         !   3 = interpolation (missing mid-levels)
                                                                         !   4 = more levels than GPH grid
#else
      real(kind=sp),dimension(:),  allocatable,public :: x_fullmet_sp    ! x-coordinates of full met grid
      real(kind=sp),dimension(:),  allocatable,public :: y_fullmet_sp    ! y-coordinates of full met grid
      real(kind=sp),dimension(:),  allocatable,public :: p_fullmet_sp    ! z-coordinates of full met grid for H
      real(kind=sp),dimension(:,:),allocatable,public :: levs_fullmet_sp ! This hold each of the numbered level coordinates:
                                                                  !    i.e. isobaric, isobaric1, isobaric2, but also
                                                                  !    height_above_ground, depth_below_surface_layer, etc.
                                                                  !    p_fullmet_sp,p_fullmet_[Vz,RH]sp are copies of one
                                                                  !    of the slices
      integer,dimension(:), allocatable,public :: nlevs_fullmet
      integer,dimension(:), allocatable,public :: levs_code       ! code indicating how to map to the GPH grid
                                                                  !   0 = no mapping (not a pressure coordinate)
                                                                  !   1 = one-to-one mapping (U,V)
                                                                  !   2 = upper truncation (missing upper levels)
                                                                  !   3 = interpolation (missing mid-levels)
                                                                  !   4 = more levels than GPH grid
#endif

      logical,public :: IsLatLon_MetGrid
      logical,public :: IsGlobal_MetGrid     = .false.  ! Not all Lon/Lat grids are periodic
      logical,public :: IsLatLon_CompGrid
      logical,public :: IsPeriodic_CompGrid  = .false.
      logical,public :: UseFullMetGrid       = .false.  ! This is the special case where the comp grid
                                                        ! equals the Met grid
      logical,public :: isGridRelative = .true.  ! Most windfiles, whether Lat/Lon or projected, give
                                          ! velocities relative to the grid of the file.  Some (NARR)
                                          ! give velocities relative to earth coordinates and need to
                                          ! be rotated
      logical,public       :: MR_Have_LL_mapping = .false.
      logical,public       :: IsRegular_MetGrid                ! True if the grid-spacing is uniform
      real(kind=sp),public :: dx_met_const
      real(kind=sp),public :: dy_met_const

      integer,public       :: nx_submet ! length of x or lon of sub-grid
      integer,public       :: ny_submet ! length of y or lat of sub-grid
#ifdef USEPOINTERS
      real(kind=sp),dimension(:),  pointer,public :: x_submet_sp => null() ! x-coordinates of met sub-grid
      real(kind=sp),dimension(:),  pointer,public :: y_submet_sp => null() ! y-coordinates of met sub-grid
      real(kind=sp),dimension(:),  pointer,public :: z_approx    => null() ! approximate altidue from STD Atmos and press (in km)
      real(kind=dp),dimension(:,:),pointer,public :: theta_Met   => null() ! Earth to grid
      real(kind=dp),dimension(:,:),pointer,public :: theta_Comp  => null() ! Earth to grid
      real(kind=dp),dimension(:,:),pointer,public :: MR_xy2ll_xlon  => null() ! holds longitude of projected grid
      real(kind=dp),dimension(:,:),pointer,public :: MR_xy2ll_ylat  => null() ! holds latitude of projected grid
!      real(kind=dp),dimension(:,:),pointer,public :: MR_ll2xy_lonx  => null() ! holds x projected value of lon/lat grid
!      real(kind=dp),dimension(:,:),pointer,public :: MR_ll2xy_laty  => null() ! holds y projected value of lon/lat grid
#else
      real(kind=sp),dimension(:)  ,allocatable,public :: x_submet_sp ! x-coordinates of met sub-grid
      real(kind=sp),dimension(:)  ,allocatable,public :: y_submet_sp ! y-coordinates of met sub-grid
      real(kind=sp),dimension(:)  ,allocatable,public :: z_approx ! approximate altitude from STD Atmos and press (in km)
      real(kind=dp),dimension(:,:),allocatable,public :: theta_Met  ! Earth to grid
      real(kind=dp),dimension(:,:),allocatable,public :: theta_Comp ! Earth to grid
      real(kind=dp),dimension(:,:),allocatable,public :: MR_xy2ll_xlon  ! holds longitude of projected grid
      real(kind=dp),dimension(:,:),allocatable,public :: MR_xy2ll_ylat  ! holds latitude of projected grid
!      real(kind=dp),dimension(:,:),allocatable,public :: MR_ll2xy_lonx  ! holds x projected value of lon/lat grid
!      real(kind=dp),dimension(:,:),allocatable,public :: MR_ll2xy_laty  ! holds y projected value of lon/lat grid
#endif
      real(kind=sp),public :: MR_lonmin,MR_lonmax
      real(kind=sp),public :: MR_latmin,MR_latmax

      logical,public :: x_inverted     = .false.
      logical,public :: y_inverted     = .false. ! Some LatLon grids start at the North Pole and increment down
      logical,public :: z_inverted     = .false. ! Some grids give top pressure first
      logical,public :: y_pad_North    = .false. ! Some computational grids will require values above the top met point
      logical,public :: y_pad_South    = .false. !   

      ! Met copies of projection variables, used for proj call on Met Grid
      character(len=4),public :: Met_gridtype
      integer     ,public :: Met_iprojflag
      real(kind=8),public :: Met_Re
      real(kind=8),public :: Met_k0
      real(kind=8),public :: Met_phi0            ! latitude of projection point
      real(kind=8),public :: Met_phi1
      real(kind=8),public :: Met_phi2
      real(kind=8),public :: Met_lam0            ! longitude of projection point
      real(kind=8),public :: Met_lam1
      real(kind=8),public :: Met_lam2
      character(len=80),public :: Met_proj4

      integer     ,public :: Comp_iprojflag
      real(kind=8),public :: Comp_Re
      real(kind=8),public :: Comp_k0
      real(kind=8),public :: Comp_phi0           ! latitude of projection point
      real(kind=8),public :: Comp_phi1
      real(kind=8),public :: Comp_phi2
      real(kind=8),public :: Comp_lam0           ! longitude of projection point
      real(kind=8),public :: Comp_lam1
      real(kind=8),public :: Comp_lam2
      character(len=100),public :: Comp_proj4

      integer     ,public :: Map_Case

      ! Some geometry terms
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:)     ,pointer            :: rdphi_MetP_sp      => null()
      real(kind=sp),dimension(:,:,:)   ,pointer            :: rdlambda_MetP_sp   => null()
      real(kind=sp),dimension(:)       ,pointer    ,public :: MR_dx_met          => null()
      real(kind=sp),dimension(:)       ,pointer    ,public :: MR_dx_submet       => null()
      real(kind=sp),dimension(:)       ,pointer    ,public :: MR_dy_met          => null()
      real(kind=sp),dimension(:)       ,pointer    ,public :: MR_dy_submet       => null()
      real(kind=sp),dimension(:,:)     ,pointer    ,public :: MR_sigma_nz_submet => null()
#else
      real(kind=sp),dimension(:,:)     ,allocatable        :: rdphi_MetP_sp
      real(kind=sp),dimension(:,:,:)   ,allocatable        :: rdlambda_MetP_sp
      real(kind=sp),dimension(:)       ,allocatable,public :: MR_dx_met
      real(kind=sp),dimension(:)       ,allocatable,public :: MR_dx_submet
      real(kind=sp),dimension(:)       ,allocatable,public :: MR_dy_met
      real(kind=sp),dimension(:)       ,allocatable,public :: MR_dy_submet
      real(kind=sp),dimension(:,:)     ,allocatable,public :: MR_sigma_nz_submet
#endif
      real(kind=sp),public :: MR_minlen = 100000.0_sp    ! minimum length of the met subgrid (m)

        ! There are some computational grid variables we might need, so make local copies
      integer      ,public :: MR_BaseYear            = 1900    ! This should be reset in calling program
      logical      ,public :: MR_useLeap             = .true.  ! This too
      integer      ,public :: MR_Comp_StartYear
      integer      ,public :: MR_Comp_StartMonth
      integer      ,public :: MR_Comp_StartDay
      real(kind=dp),public :: MR_Comp_StartHour      = 0.0_dp  ! Note that these must be double-precision to
      real(kind=dp),public :: MR_Comp_Time_in_hours  = 0.0_dp  ! be passed correctly to HoursSince
      integer      ,public :: nx_comp
      integer      ,public :: ny_comp
      integer      ,public :: nz_comp
      real(kind=sp) :: dx_comp,dy_comp
      real(kind=sp) :: MaxZ_comp_sp
#ifdef USEPOINTERS
      real(kind=sp),dimension(:),    pointer,public :: x_comp_sp               => null() ! x-coordinates of computational grid
      real(kind=sp),dimension(:),    pointer,public :: y_comp_sp               => null() ! y-coordinates of computational grid
      real(kind=sp),dimension(:),    pointer,public :: z_comp_sp               => null() ! z-coordinates of computational grid
      real(kind=sp),dimension(:,:),  pointer,public :: CompPoint_X_on_Met_sp   => null() ! x-coord (on Met grid) of comp point
      real(kind=sp),dimension(:,:),  pointer,public :: CompPoint_Y_on_Met_sp   => null() ! y-coord (on Met grid) of comp point
      integer      ,dimension(:,:,:),pointer,public :: CompPoint_on_subMet_idx => null() ! index on met sub-grid of comp point
      real(kind=sp),dimension(:,:,:),pointer,public :: bilin_map_wgt           => null()
#else
      real(kind=sp),dimension(:),    allocatable,public :: x_comp_sp               ! x-coordinates of computational grid
      real(kind=sp),dimension(:),    allocatable,public :: y_comp_sp               ! y-coordinates of computational grid
      real(kind=sp),dimension(:),    allocatable,public :: z_comp_sp               ! z-coordinates of computational grid
      real(kind=sp),dimension(:,:),  allocatable,public :: CompPoint_X_on_Met_sp   ! x-coord (on Met grid) of comp point
      real(kind=sp),dimension(:,:),  allocatable,public :: CompPoint_Y_on_Met_sp   ! y-coord (on Met grid) of comp point
      integer      ,dimension(:,:,:),allocatable,public :: CompPoint_on_subMet_idx ! index on met sub-grid of comp point
      real(kind=sp),dimension(:,:,:),allocatable,public :: bilin_map_wgt
#endif

      ! Here are a few variables needed for sigma-altitude coordinates
      logical          ,public :: MR_useTopo        = .false.
      integer          ,public :: MR_ZScaling_ID    = 0  ! = 0 for no scaling (i.e. s = z)
                                                         ! = 1 for shifted-altitude (s=z-zsurf)
                                                         ! = 2 for sigma-altitude (s=Ztop*(z-zsurf)/(ztop-zsurf))
      real(kind=sp)    ,public :: MR_ztop
#ifdef USEPOINTERS
      real(kind=sp),dimension(:)   ,pointer, public :: s_comp_sp     => null() ! s-coordinates (scaled z) of comp. grid
      real(kind=sp),dimension(:,:) ,pointer, public :: MR_Topo_met   => null()
      real(kind=sp),dimension(:,:) ,pointer, public :: MR_Topo_comp  => null()
      real(kind=sp),dimension(:,:) ,pointer, public :: MR_jacob_met  => null() ! Jacobian of trans. = MR_ztop-MR_Topo_met
      real(kind=sp),dimension(:,:) ,pointer, public :: MR_jacob_comp => null() 
#else
      real(kind=sp),dimension(:)   ,allocatable, public :: s_comp_sp ! s-coordinates (scaled z) of computational grid
      real(kind=sp),dimension(:,:) ,allocatable, public :: MR_Topo_met
      real(kind=sp),dimension(:,:) ,allocatable, public :: MR_Topo_comp
      real(kind=sp),dimension(:,:) ,allocatable, public :: MR_jacob_met  ! Jacobian of trans. = MR_ztop-MR_Topo_met
      real(kind=sp),dimension(:,:) ,allocatable, public :: MR_jacob_comp
#endif

      logical                       ,public :: FoundFillVAttr = .false.
      real(kind=sp),public                  :: fill_value_sp  = -9999.0_sp  ! Initialize to the commonly used fill value
      character(len=30),dimension(9),public :: Met_dim_names      ! name of dimension
      logical          ,dimension(9),public :: Met_dim_IsAvailable
      real(kind=sp)    ,dimension(9),public :: Met_dim_fac  = 1.0_sp
      ! Here is the list of variables that can be read.  Each iwindformat will
      ! have just a sub-set available with specific names.  For now, allocate
      ! space for 50 variable names
        ! Mechanical / State variables
        !   1 = Geopotential Height              : m
        !   2 = Vx                               : m/s
        !   3 = Vy                               : m/s
        !   4 = Vz                               : m/s
        !   5 = Temperature                      : K
        !   6 = pressure (when 3d for eta grids) : Pa
        !   7 = PVV (Pressure Vertical Vel.)     : Pa/s
        ! Surface
        !  10 = Planetary Boundary Layer Height
        !  11 = U @ 10m
        !  12 = V @ 10m
        !  13 = Friction velocity
        !  15 = Snow cover
        !  16 = Soil moisture
        !  17 = Surface Roughness
        !  18 = Wind gust speed
        !  19 = surface temperature
        ! Atmospheric Structure
        !  20 = pressure at lower cloud base
        !  21 = pressure at lower cloud top
        !  22 = temperature at lower cloud top
        !  23 = Total Cloud cover
        !  24 = Cloud cover (low)
        !  25 = Cloud cover (convective)
        ! Moisture
        !  30 = Rel. Hum
        !  31 = QV (specific humidity)
        !  32 = QL (liquid)
        !  33 = QI (ice)
        ! Precipitation
        !  40 = Categorical rain
        !  41 = Categorical snow
        !  42 = Categorical frozen rain
        !  43 = Categorical ice
        !  44 = Precipitation rate large-scale (liquid)
        !  45 = Precipitation rate convective  (liquid)
        !  46 = Precipitation rate large-scale (ice)
        !  47 = Precipitation rate convective  (ice)
      logical          ,dimension(MR_MAXVARS),public   :: Met_var_IsAvailable      ! true if iwf contains the var
      character(len=80),dimension(MR_MAXVARS),public   :: Met_var_NC_names         ! name in the file
      character(len=80),dimension(MR_MAXVARS),public   :: Met_var_GRIB_names       ! name in the file
      character(len=5) ,dimension(MR_MAXVARS),public   :: Met_var_WMO_names        ! WMO version of the name
      integer          ,dimension(MR_MAXVARS),public   :: Met_var_ndim             ! number of expected dimensions for this variable
      integer          ,dimension(MR_MAXVARS),public   :: Met_var_zdim_idx         ! The index of this coordinate (used in Met_var_nlevs) 
      integer          ,dimension(MR_MAXVARS),public   :: Met_var_zdim_ncid        ! The dimID of the dimension in the nc file
      integer                                ,public   :: nlev_coords_detected = 0
      integer          ,dimension(MR_MAXVARS,4),public :: Met_var_GRIB2_DPcPnSt    ! Grib2 files have variables identified by
                                                                            ! discpln,param_cat,param_num,surf_class
      integer          ,dimension(MR_MAXVARS),public   :: Met_var_GRIB1_Param      !   indicatorOfParameter
      integer          ,dimension(MR_MAXVARS),public   :: Met_var_GRIB1_Table      !   table2Version
      character(len=3) ,dimension(MR_MAXVARS),public   :: Met_var_GRIB1_St         ! level type: pl (100), sfc (001), etc
      integer                                ,public   :: MR_GRIB_Version  = 0

      real(kind=sp)    ,dimension(MR_MAXVARS),public   :: Met_var_conversion_factor

      integer          ,dimension(MR_MAXVARS)          :: Met_var_nlevs

        ! Variables needed by netcdf reader
      real(kind=sp),public :: iwf25_scale_facs(MR_MAXVARS)
      real(kind=sp),public :: iwf25_offsets(MR_MAXVARS)

      integer,public :: istart
      integer,public :: iend
      integer,public :: jstart
      integer,public :: jend
      integer,public :: ilhalf_fm_l
      integer,public :: ilhalf_fm_r
      integer,public :: irhalf_fm_l
      integer,public :: irhalf_fm_r
      integer,public :: ilhalf_nx
      integer,public :: irhalf_nx
      logical,public :: wrapgrid

#ifdef USEPOINTERS
      character       ,dimension(:)       ,pointer,public :: temp1d_byte   => null()
      character       ,dimension(:)       ,pointer,public :: temp1d_char   => null()
      integer(kind=2) ,dimension(:)       ,pointer,public :: temp1d_intS   => null()
      integer(kind=4) ,dimension(:)       ,pointer,public :: temp1d_intL   => null()
      real(kind=sp)   ,dimension(:)       ,pointer,public :: temp1d_sp     => null()
      real(kind=dp)   ,dimension(:)       ,pointer,public :: temp1d_dp     => null()
      real(kind=sp)   ,dimension(:,:,:)   ,pointer,public :: temp2d_sp     => null()
      real(kind=sp)   ,dimension(:,:,:,:) ,pointer,public :: temp3d_sp     => null()
      integer(kind=sp),dimension(:,:,:)   ,pointer,public :: temp2d_int    => null()
      integer(kind=sp),dimension(:,:,:)   ,pointer        :: temp2d_short  => null()
      integer(kind=sp),dimension(:,:,:,:) ,pointer,public :: temp3d_short  => null()
      real(kind=dp)   ,dimension(:,:,:,:) ,pointer,public :: temp3d_dp     => null()
      real(kind=4)    ,dimension(:,:)     ,pointer,public :: Met_Proj_lat  => null()
      real(kind=4)    ,dimension(:,:)     ,pointer,public :: Met_Proj_lon  => null()
#else
      character       ,dimension(:)       ,allocatable,public :: temp1d_byte
      character       ,dimension(:)       ,allocatable,public :: temp1d_char
      integer(kind=2) ,dimension(:)       ,allocatable,public :: temp1d_intS
      integer(kind=4) ,dimension(:)       ,allocatable,public :: temp1d_intL
      real(kind=sp)   ,dimension(:)       ,allocatable,public :: temp1d_sp
      real(kind=dp)   ,dimension(:)       ,allocatable,public :: temp1d_dp
      real(kind=sp)   ,dimension(:,:,:)   ,allocatable,public :: temp2d_sp
      real(kind=sp)   ,dimension(:,:,:,:) ,allocatable,public :: temp3d_sp
      integer(kind=sp),dimension(:,:,:)   ,allocatable,public :: temp2d_int
      integer(kind=sp),dimension(:,:,:)   ,allocatable        :: temp2d_short
      real(kind=dp)   ,dimension(:,:,:,:) ,allocatable,public :: temp3d_dp
      integer(kind=sp),dimension(:,:,:,:) ,allocatable,public :: temp3d_short
      real(kind=4)    ,dimension(:,:)     ,allocatable,public :: Met_Proj_lat
      real(kind=4)    ,dimension(:,:)     ,allocatable,public :: Met_Proj_lon
#endif

      ! Status variables for error-checking
      logical :: Check_prereq_conditions            = .true.
      logical :: CALLED_MR_Allocate_FullMetFileList = .false.
      logical :: CALLED_MR_Read_Met_DimVars         = .false.
      logical :: CALLED_MR_Set_CompProjection       = .false.
      logical :: CALLED_MR_Initialize_Met_Grids     = .false.
      logical :: CALLED_MR_Set_Met_Times            = .false.

      ! Info on how to build file paths
      ! Defaults to linux-style, but this can be changed in the calling program
      integer          ,public    :: MR_OS_TYPE    = 1       ! 1=Linux, 2=MacOS, 3=MS Windows
      character (len=2),public    :: MR_DirPrefix  = ''
      character (len=1),public    :: MR_DirDelim   = '/'

      contains

!##############################################################################
!
!     MR_Reset_Memory
!
!     This subroutine reinitializes MetReader by deallocating all MR variables.
!     This is useful if a program needs to use multiple types of wind files.
!
!##############################################################################

      subroutine MR_Reset_Memory

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-------------------------------------------------------"
        write(outlog(io),*)"-------- Resetting all MetReader Memory ---------------"
        write(outlog(io),*)"-------------------------------------------------------"
      endif;enddo

#ifdef USEPOINTERS
       if(associated(MR_windfile_starthour         ))deallocate(MR_windfile_starthour)
       if(associated(MR_windfile_stephour          ))deallocate(MR_windfile_stephour)
       if(associated(MR_windfiles_nt_fullmet       ))deallocate(MR_windfiles_nt_fullmet)
       if(associated(MR_windfiles_Have_GRIB_index  ))deallocate(MR_windfiles_Have_GRIB_index)
       if(associated(MR_windfiles_GRIB_index       ))deallocate(MR_windfiles_GRIB_index)
       if(associated(MR_MetStep_File               ))deallocate(MR_MetStep_File)
       if(associated(MR_MetStep_findex             ))deallocate(MR_MetStep_findex)
       if(associated(MR_MetStep_tindex             ))deallocate(MR_MetStep_tindex)
       if(associated(MR_MetStep_Hour_since_baseyear))deallocate(MR_MetStep_Hour_since_baseyear)
       if(associated(MR_MetStep_Interval           ))deallocate(MR_MetStep_Interval)
       if(associated(MR_MetStep_year               ))deallocate(MR_MetStep_year)
       if(associated(MR_MetStep_month              ))deallocate(MR_MetStep_month)
       if(associated(MR_MetStep_day                ))deallocate(MR_MetStep_day)
       if(associated(MR_MetStep_DOY                ))deallocate(MR_MetStep_DOY)
       if(associated(MR_MetStep_Hour_Of_Day        ))deallocate(MR_MetStep_Hour_Of_Day)
       if(associated(MR_iwind5_year                ))deallocate(MR_iwind5_year)
       if(associated(MR_windfiles                  ))deallocate(MR_windfiles)
       if(associated(MR_dum2d_met_int              ))deallocate(MR_dum2d_met_int)
       if(associated(MR_dum2d_met                  ))deallocate(MR_dum2d_met)
       if(associated(MR_dum3d_metP                 ))deallocate(MR_dum3d_metP)
       if(associated(MR_dum3d2_metP                ))deallocate(MR_dum3d2_metP)
       if(associated(MR_geoH_metP_last             ))deallocate(MR_geoH_metP_last)
       if(associated(MR_geoH_metP_next             ))deallocate(MR_geoH_metP_next)
       if(associated(MR_vx_metP_last               ))deallocate(MR_vx_metP_last)
       if(associated(MR_vx_metP_next               ))deallocate(MR_vx_metP_next)
       if(associated(MR_vy_metP_last               ))deallocate(MR_vy_metP_last)
       if(associated(MR_vy_metP_next               ))deallocate(MR_vy_metP_next)
       if(associated(MR_dum3d_metH                 ))deallocate(MR_dum3d_metH)
       if(associated(MR_dum2d_comp_int             ))deallocate(MR_dum2d_comp_int)
       if(associated(MR_dum2d_comp                 ))deallocate(MR_dum2d_comp)
       if(associated(MR_dum3d_compP                ))deallocate(MR_dum3d_compP)
       if(associated(MR_dum3d_compP_2              ))deallocate(MR_dum3d_compP_2)
       if(associated(MR_dum3d_compH                ))deallocate(MR_dum3d_compH)
       if(associated(MR_dum3d_compH_2              ))deallocate(MR_dum3d_compH_2)
       if(associated(x_fullmet_sp                  ))deallocate(x_fullmet_sp)
       if(associated(y_fullmet_sp                  ))deallocate(y_fullmet_sp)
       if(associated(p_fullmet_sp                  ))deallocate(p_fullmet_sp)
       if(associated(x_submet_sp                   ))deallocate(x_submet_sp)
       if(associated(y_submet_sp                   ))deallocate(y_submet_sp)
       if(associated(levs_fullmet_sp               ))deallocate(levs_fullmet_sp)
       if(associated(nlevs_fullmet                 ))deallocate(nlevs_fullmet)
       if(associated(levs_code                     ))deallocate(levs_code)
       if(associated(z_approx                      ))deallocate(z_approx)
       if(associated(rdphi_MetP_sp                 ))deallocate(rdphi_MetP_sp)
       if(associated(rdlambda_MetP_sp              ))deallocate(rdlambda_MetP_sp)
       if(associated(MR_dx_met                     ))deallocate(MR_dx_met)
       if(associated(MR_dx_submet                  ))deallocate(MR_dx_submet)
       if(associated(MR_dy_met                     ))deallocate(MR_dy_met)
       if(associated(MR_dy_submet                  ))deallocate(MR_dy_submet)
       if(associated(MR_sigma_nz_submet            ))deallocate(MR_sigma_nz_submet)
       if(associated(x_comp_sp                     ))deallocate(x_comp_sp)
       if(associated(y_comp_sp                     ))deallocate(y_comp_sp)
       if(associated(z_comp_sp                     ))deallocate(z_comp_sp)
       if(associated(CompPoint_X_on_Met_sp         ))deallocate(CompPoint_X_on_Met_sp)
       if(associated(CompPoint_Y_on_Met_sp         ))deallocate(CompPoint_Y_on_Met_sp)
       if(associated(CompPoint_on_subMet_idx       ))deallocate(CompPoint_on_subMet_idx)
       if(associated(bilin_map_wgt                 ))deallocate(bilin_map_wgt)
       if(associated(s_comp_sp                     ))deallocate(s_comp_sp)
       if(associated(MR_Topo_met                   ))deallocate(MR_Topo_met)
       if(associated(MR_jacob_met                  ))deallocate(MR_jacob_met)
       if(associated(MR_Topo_comp                  ))deallocate(MR_Topo_comp)
       if(associated(MR_jacob_comp                 ))deallocate(MR_jacob_comp)
       if(associated(MR_u_ER_metP                  ))deallocate(MR_u_ER_metP)
       if(associated(MR_v_ER_metP                  ))deallocate(MR_v_ER_metP)
       if(associated(theta_Met                     ))deallocate(theta_Met)
       if(associated(theta_Comp                    ))deallocate(theta_Comp)
       if(associated(MR_xy2ll_xlon                 ))deallocate(MR_xy2ll_xlon)
       if(associated(MR_xy2ll_ylat                 ))deallocate(MR_xy2ll_ylat)
       if(associated(temp1d_sp                     ))deallocate(temp1d_sp)
       if(associated(temp2d_sp                     ))deallocate(temp2d_sp)
       if(associated(temp3d_sp                     ))deallocate(temp3d_sp)
       if(associated(temp2d_int                    ))deallocate(temp2d_int)
       if(associated(temp2d_short                  ))deallocate(temp2d_short)
       if(associated(temp3d_short                  ))deallocate(temp3d_short)
       if(associated(Met_Proj_lat                  ))deallocate(Met_Proj_lat)
       if(associated(Met_Proj_lon                  ))deallocate(Met_Proj_lon)
       if(associated(MR_SndVars_metP               ))deallocate(MR_SndVars_metP)
       if(associated(MR_SndVarsID                  ))deallocate(MR_SndVarsID)
       if(associated(MR_Snd_np_fullmet             ))deallocate(MR_Snd_np_fullmet)
       if(associated(MR_Snd2Comp_map_wgt           ))deallocate(MR_Snd2Comp_map_wgt)
       if(associated(MR_Snd2Comp_map_idx           ))deallocate(MR_Snd2Comp_map_idx)
       if(associated(Snd_idx                       ))deallocate(Snd_idx)
       if(associated(MR_Snd_cd                     ))deallocate(MR_Snd_cd)
       if(associated(MR_Snd_id                     ))deallocate(MR_Snd_id)
       if(associated(MR_Snd_lt                     ))deallocate(MR_Snd_lt)
       if(associated(MR_Snd_ln                     ))deallocate(MR_Snd_ln)
       if(associated(MR_Snd_el                     ))deallocate(MR_Snd_el)
       if(associated(MR_Snd_lnm                    ))deallocate(MR_Snd_lnm)
       if(associated(MR_Snd_st                     ))deallocate(MR_Snd_st)
       if(associated(MR_Snd_ct                     ))deallocate(MR_Snd_ct)
#else
       if(allocated(MR_windfile_starthour          ))deallocate(MR_windfile_starthour)
       if(allocated(MR_windfile_stephour           ))deallocate(MR_windfile_stephour)
       if(allocated(MR_windfiles_nt_fullmet        ))deallocate(MR_windfiles_nt_fullmet)
       if(allocated(MR_windfiles_Have_GRIB_index   ))deallocate(MR_windfiles_Have_GRIB_index)
       if(allocated(MR_windfiles_GRIB_index        ))deallocate(MR_windfiles_GRIB_index)
       if(allocated(MR_MetStep_File                ))deallocate(MR_MetStep_File)
       if(allocated(MR_MetStep_findex              ))deallocate(MR_MetStep_findex)
       if(allocated(MR_MetStep_tindex              ))deallocate(MR_MetStep_tindex)
       if(allocated(MR_MetStep_Hour_since_baseyear ))deallocate(MR_MetStep_Hour_since_baseyear)
       if(allocated(MR_MetStep_Interval            ))deallocate(MR_MetStep_Interval)
       if(allocated(MR_MetStep_year                ))deallocate(MR_MetStep_year)
       if(allocated(MR_MetStep_month               ))deallocate(MR_MetStep_month)
       if(allocated(MR_MetStep_day                 ))deallocate(MR_MetStep_day)
       if(allocated(MR_MetStep_DOY                 ))deallocate(MR_MetStep_DOY)
       if(allocated(MR_MetStep_Hour_Of_Day         ))deallocate(MR_MetStep_Hour_Of_Day)
       if(allocated(MR_iwind5_year                 ))deallocate(MR_iwind5_year)
       if(allocated(MR_windfiles                   ))deallocate(MR_windfiles)
       if(allocated(MR_dum2d_met_int               ))deallocate(MR_dum2d_met_int)
       if(allocated(MR_dum2d_met                   ))deallocate(MR_dum2d_met)
       if(allocated(MR_dum3d_metP                  ))deallocate(MR_dum3d_metP)
       if(allocated(MR_dum3d2_metP                 ))deallocate(MR_dum3d2_metP)
       if(allocated(MR_geoH_metP_last              ))deallocate(MR_geoH_metP_last)
       if(allocated(MR_geoH_metP_next              ))deallocate(MR_geoH_metP_next)
       if(allocated(MR_vx_metP_last                ))deallocate(MR_vx_metP_last)
       if(allocated(MR_vx_metP_next                ))deallocate(MR_vx_metP_next)
       if(allocated(MR_vy_metP_last                ))deallocate(MR_vy_metP_last)
       if(allocated(MR_vy_metP_next                ))deallocate(MR_vy_metP_next)
       if(allocated(MR_dum3d_metH                  ))deallocate(MR_dum3d_metH)
       if(allocated(MR_dum2d_comp_int              ))deallocate(MR_dum2d_comp_int)
       if(allocated(MR_dum2d_comp                  ))deallocate(MR_dum2d_comp)
       if(allocated(MR_dum3d_compP                 ))deallocate(MR_dum3d_compP)
       if(allocated(MR_dum3d_compP_2               ))deallocate(MR_dum3d_compP_2)
       if(allocated(MR_dum3d_compH                 ))deallocate(MR_dum3d_compH)
       if(allocated(MR_dum3d_compH_2               ))deallocate(MR_dum3d_compH_2)
       if(allocated(x_fullmet_sp                   ))deallocate(x_fullmet_sp)
       if(allocated(y_fullmet_sp                   ))deallocate(y_fullmet_sp)
       if(allocated(p_fullmet_sp                   ))deallocate(p_fullmet_sp)
       if(allocated(x_submet_sp                    ))deallocate(x_submet_sp)
       if(allocated(y_submet_sp                    ))deallocate(y_submet_sp)
       if(allocated(levs_fullmet_sp                ))deallocate(levs_fullmet_sp)
       if(allocated(nlevs_fullmet                  ))deallocate(nlevs_fullmet)
       if(allocated(levs_code                      ))deallocate(levs_code)
       if(allocated(z_approx                       ))deallocate(z_approx)
       if(allocated(rdphi_MetP_sp                  ))deallocate(rdphi_MetP_sp)
       if(allocated(rdlambda_MetP_sp               ))deallocate(rdlambda_MetP_sp)
       if(allocated(MR_dx_met                      ))deallocate(MR_dx_met)
       if(allocated(MR_dx_submet                   ))deallocate(MR_dx_submet)
       if(allocated(MR_dy_met                      ))deallocate(MR_dy_met)
       if(allocated(MR_dy_submet                   ))deallocate(MR_dy_submet)
       if(allocated(MR_sigma_nz_submet             ))deallocate(MR_sigma_nz_submet)
       if(allocated(x_comp_sp                      ))deallocate(x_comp_sp)
       if(allocated(y_comp_sp                      ))deallocate(y_comp_sp)
       if(allocated(z_comp_sp                      ))deallocate(z_comp_sp)
       if(allocated(CompPoint_X_on_Met_sp          ))deallocate(CompPoint_X_on_Met_sp)
       if(allocated(CompPoint_Y_on_Met_sp          ))deallocate(CompPoint_Y_on_Met_sp)
       if(allocated(CompPoint_on_subMet_idx        ))deallocate(CompPoint_on_subMet_idx)
       if(allocated(bilin_map_wgt                  ))deallocate(bilin_map_wgt)
       if(allocated(s_comp_sp                      ))deallocate(s_comp_sp)
       if(allocated(MR_Topo_met                    ))deallocate(MR_Topo_met)
       if(allocated(MR_jacob_met                   ))deallocate(MR_jacob_met)
       if(allocated(MR_Topo_comp                   ))deallocate(MR_Topo_comp)
       if(allocated(MR_jacob_comp                  ))deallocate(MR_jacob_comp)
       if(allocated(MR_u_ER_metP                   ))deallocate(MR_u_ER_metP)
       if(allocated(MR_v_ER_metP                   ))deallocate(MR_v_ER_metP)
       if(allocated(theta_Met                      ))deallocate(theta_Met)
       if(allocated(theta_Comp                     ))deallocate(theta_Comp)
       if(allocated(MR_xy2ll_xlon                  ))deallocate(MR_xy2ll_xlon)
       if(allocated(MR_xy2ll_ylat                  ))deallocate(MR_xy2ll_ylat)
       if(allocated(temp1d_sp                      ))deallocate(temp1d_sp)
       if(allocated(temp2d_sp                      ))deallocate(temp2d_sp)
       if(allocated(temp3d_sp                      ))deallocate(temp3d_sp)
       if(allocated(temp2d_int                     ))deallocate(temp2d_int)
       if(allocated(temp2d_short                   ))deallocate(temp2d_short)
       if(allocated(temp3d_short                   ))deallocate(temp3d_short)
       if(allocated(Met_Proj_lat                   ))deallocate(Met_Proj_lat)
       if(allocated(Met_Proj_lon                   ))deallocate(Met_Proj_lon)
       if(allocated(MR_SndVars_metP                ))deallocate(MR_SndVars_metP)
       if(allocated(MR_SndVarsID                   ))deallocate(MR_SndVarsID)
       if(allocated(MR_Snd_np_fullmet              ))deallocate(MR_Snd_np_fullmet)
       if(allocated(MR_Snd2Comp_map_wgt            ))deallocate(MR_Snd2Comp_map_wgt)
       if(allocated(MR_Snd2Comp_map_idx            ))deallocate(MR_Snd2Comp_map_idx)
       if(allocated(Snd_idx                        ))deallocate(Snd_idx)
       if(allocated(MR_Snd_cd                      ))deallocate(MR_Snd_cd)
       if(allocated(MR_Snd_id                      ))deallocate(MR_Snd_id)
       if(allocated(MR_Snd_lt                      ))deallocate(MR_Snd_lt)
       if(allocated(MR_Snd_ln                      ))deallocate(MR_Snd_ln)
       if(allocated(MR_Snd_el                      ))deallocate(MR_Snd_el)
       if(allocated(MR_Snd_lnm                     ))deallocate(MR_Snd_lnm)
       if(allocated(MR_Snd_st                      ))deallocate(MR_Snd_st)
       if(allocated(MR_Snd_ct                      ))deallocate(MR_Snd_ct)
#endif

       nlev_coords_detected = 0

      end subroutine MR_Reset_Memory

!##############################################################################
!
!     MR_Allocate_FullMetFileList
!
!     This subroutine allocates the list of windfiles and does some
!     error-checking based on iwind and iwindfiles.
!
!     From the calling program, this is called once the information about
!     the NWP files is available (e.g. the iwind, iwindformat, grid ID, data
!     format, and number of windfiles.
!
!     Variables allocated:
!        MR_windfiles(MR_iwindfiles)
!        MR_windfiles_nt_fullmet(MR_iwindfiles)
!       and possibly:
!        MR_windfiles_Have_GRIB_index(MR_iwindfiles)
!        MR_windfiles_GRIB_index(MR_iwindfiles)
!
!     The next step is for the calling program to populate MR_windfiles
!
!##############################################################################

      subroutine MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

      integer,intent(in)      :: iw
      integer,intent(in)      :: iwf
      integer,intent(in)      :: igrid
      integer,intent(in)      :: idf
      integer,intent(in)      :: iwfiles

      integer :: i
      integer :: io                           ! Index for output streams

      MR_iwind              = iw
      MR_iwindformat        = iwf
      MR_iGridCode          = igrid
      MR_idataFormat        = idf
      MR_iwindfiles         = iwfiles

      ! This should be the first subroutine called so reset verbosity levels
      ! in case this was changed in the calling program
      VB(1:2) = (/MR_VERB,MR_VERB/)

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Allocate_FullMetFileList                  ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo
      if ((MR_iwind.ne.1).and.(MR_iwind.ne.2).and. &
          (MR_iwind.ne.3).and.(MR_iwind.ne.4).and. &
          (MR_iwind.ne.5)) then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'MR_iwind must be between 1 and 5. Program stopped'
          write(errlog(io),*)' MR_iwind = ',MR_iwind
          write(errlog(io),*)' MR_IWIND OPTIONS:'
          write(errlog(io),*)' MR_iwind = 1 read from a 1-D wind sounding'
          write(errlog(io),*)'            2 read from 3D gridded ASCII files'
          write(errlog(io),*)'            3 read from a single, multistep file'
          write(errlog(io),*)'            4 read from multiple files'
          write(errlog(io),*)'            5 read variables from separate files'
        endif;enddo
        stop 1
      endif

      ! Initialize the dimension and variable arrays.  Select slots in these arrays will be
      ! overwritten from the calls in the case block below
      Met_dim_IsAvailable                     = .false.
      Met_var_IsAvailable                     = .false.
      Met_var_IsAvailable(1:MR_MAXVARS)       = .false.
      Met_var_zdim_idx(1:MR_MAXVARS)          = 0
      Met_var_zdim_ncid(1:MR_MAXVARS)         = 0
      Met_var_GRIB2_DPcPnSt(1:MR_MAXVARS,1:4) = 0
      Met_var_GRIB1_Param(1:MR_MAXVARS)       = 0
      Met_var_GRIB1_St(1:MR_MAXVARS)          = ""
      Met_var_conversion_factor(1:MR_MAXVARS) = 1.0_sp
      Met_var_nlevs(1:MR_MAXVARS)             = 0
      IsGridRelative = .true.

      !--------------------------------
      ! Dimensions
      !--------------------------------
      ! Assume the critical dimensions are available
      ! Dimension names for NC files are read from the first file so do not need
      ! to be set here. These could be populated by the template file.
      Met_dim_IsAvailable(1)=.true.; Met_dim_names(1)="unset name for t"
      Met_dim_IsAvailable(2)=.true.; Met_dim_names(2)="unset name for p"
      Met_dim_IsAvailable(3)=.true.; Met_dim_names(3)="unset name for y"
      Met_dim_IsAvailable(4)=.true.; Met_dim_names(4)="unset name for x"

      !--------------------------------
      ! Mechanical / State variables
      !--------------------------------
      !   Note: All default names are that assigned by netcdf-java in the grib-to-nc conversion
      !         Command used :: java -Xmx2048m -classpath ~/ncj/netcdfAll-4.5.jar ucar.nc2.dataset.NetcdfDataset \
      !                          -in ${GribFile} -out ${NCFile} -IsLargeFile
      ! 
      !         Alternatively, one could use ncl_convert2nc ${GribFile} -L
      !         This creates variable names such as HGT_P0_L100_GLL0; i.e. WMOname_P[param_class]_L[surf_class]_[gridtype]
      !         where gridtype is one of:
      !          Met_gridtype = "GLL0"  ! Latitude/Longitude
      !          Met_gridtype = "GME0"  ! Mercator
      !          Met_gridtype = "GLC0"  ! Lambert Conformal
      !          Met_gridtype = "GST0"  ! Polar stereographic
      !
      !         The NC_names listed are the names written by netcdf-java when converting GRIB forecast 
      !         files to netcdf.  Reanalysis files might be provided in netcdf format directly, in
      !         which case, the names are overwritten below. Since these name occasionally change,
      !         know variants are checked if the expected name is not found in 
        !  Geopotential Height  (m^2/s^2)
      Met_var_NC_names(1)          = "Geopotential_height_isobaric"
      Met_var_GRIB_names(1)        = "gh"
      Met_var_WMO_names(1)         = "HGT"
      Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/) ! discpln,param_cat,param_num,surf_class
      Met_var_GRIB1_Param(1)       = 7
      Met_var_GRIB1_St(1)          = "pl"
      Met_var_ndim(1)              = 4
        ! Velocity component in x (or E) direction  (m/s)
      Met_var_NC_names(2)          = "u-component_of_wind_isobaric"
      Met_var_GRIB_names(2)        = "u"
      Met_var_WMO_names(2)         = "UGRD"
      Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
      Met_var_GRIB1_Param(2)       = 33
      Met_var_GRIB1_St(2)          = "pl"
      Met_var_ndim(2)              = 4
        ! Velocity component in y (or N) direction (m/s)
      Met_var_NC_names(3)          = "v-component_of_wind_isobaric"
      Met_var_GRIB_names(3)        = "v"
      Met_var_WMO_names(3)         = "VGRD"
      Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
      Met_var_GRIB1_Param(3)       = 34
      Met_var_GRIB1_St(3)          = "pl"
      Met_var_ndim(3)              = 4
        ! Velocity component in z direction  (Pa/s)
        !  This actually returns Vz with VPP*dp/dz with a FD approx of dp/dz
      Met_var_NC_names(4)          = "Vertical_velocity_pressure_isobaric"
      Met_var_GRIB_names(4)        = "w"
      Met_var_WMO_names(4)         = "VVEL"
      Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
      Met_var_GRIB1_Param(4)       = 39
      Met_var_GRIB1_St(4)          = "pl"
      Met_var_ndim(4)              = 4
        ! Temperature  (K)
      Met_var_NC_names(5)          = "Temperature_isobaric"
      Met_var_GRIB_names(5)        = "t"
      Met_var_WMO_names(5)         = "TMP"
      Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
      Met_var_GRIB1_Param(5)       = 11
      Met_var_GRIB1_St(5)          = "pl"
      Met_var_ndim(5)              = 4
        ! Velocity component in z direction  (Pa/s)
        !  This returns the true vert. vel. pres.
      Met_var_NC_names(7)          = "Vertical_velocity_pressure_isobaric"
      Met_var_GRIB_names(7)        = "w"
      Met_var_WMO_names(7)         = "VVEL"
      Met_var_GRIB2_DPcPnSt(7,1:4) = (/0, 2, 8, 100/)
      Met_var_GRIB1_Param(7)       = 39
      Met_var_GRIB1_St(7)          = "pl"
      Met_var_ndim(7)              = 4

      !--------------------------------
      ! Surface
      !--------------------------------
        ! Height of planetary boundary layer  (m)
      Met_var_NC_names(10)          = "Planetary_Boundary_Layer_Height_surface"
      Met_var_GRIB_names(10)        = "hpbl"
      Met_var_WMO_names(10)         = "HPBL"
      Met_var_GRIB2_DPcPnSt(10,1:4) = (/0, 3, 196, 1/)
      Met_var_GRIB1_Param(10)       = 221
      Met_var_GRIB1_St(10)          = "1"
      Met_var_ndim(10)              = 3
        ! Velocity component in x (or E) direction at 10 m above ground surface  (m/s)
      Met_var_NC_names(11)          = "u-component_of_wind_height_above_ground"
      Met_var_GRIB_names(11)        = "u"
      Met_var_WMO_names(11)         = "UGRD"
      Met_var_GRIB2_DPcPnSt(11,1:4) = (/0, 2, 2, 103/)
      Met_var_GRIB1_Param(11)       = 33
      Met_var_GRIB1_St(11)          = "105"
      Met_var_ndim(11)              = 4
        ! Velocity component in y (or N) direction at 10 m above ground surface  (m/s)
      Met_var_NC_names(12)          = "v-component_of_wind_height_above_ground"
      Met_var_GRIB_names(12)        = "v"
      Met_var_WMO_names(12)         = "VGRD"
      Met_var_GRIB2_DPcPnSt(12,1:4) = (/0, 2, 3, 103/)
      Met_var_GRIB1_Param(12)       = 34
      Met_var_GRIB1_St(12)          = "105"
      Met_var_ndim(12)              = 4
        ! Friction velocity  (m/s)
      Met_var_NC_names(13)          = "Frictional_Velocity_surface"
      Met_var_GRIB_names(13)        = "fricv"
      Met_var_WMO_names(13)         = "FRICV"
      Met_var_GRIB2_DPcPnSt(13,1:4) = (/0, 2, 197, 1/)
      Met_var_GRIB1_Param(13)       = 253
      Met_var_GRIB1_St(13)          = "1"
      Met_var_ndim(13)              = 3
        ! Snow depth  (m)
      Met_var_NC_names(15)          = "Snow_depth_surface"
      Met_var_GRIB_names(15)        = "sd"
      Met_var_WMO_names(15)         = "SNOD"
      Met_var_GRIB2_DPcPnSt(15,1:4) = (/0, 1, 11, 1/)
      Met_var_GRIB1_Param(15)       = 66
      Met_var_GRIB1_St(15)          = "1"
      Met_var_ndim(15)              = 3
        ! Soil Moisture  (fraction)
      Met_var_NC_names(16)          = "Volumetric_Soil_Moisture_Content_depth_below_surface_layer"
      Met_var_GRIB_names(16)        = "soilw"
      Met_var_WMO_names(16)         = "SOILW"
      Met_var_GRIB2_DPcPnSt(16,1:4) = (/2, 0, 192, 106/)
      Met_var_GRIB1_Param(16)       = 0
      Met_var_GRIB1_St(16)          = ""
      Met_var_ndim(16)              = 4
        ! Surface roughness  (m)
      Met_var_NC_names(17)          = "Surface_roughness_surface"
      Met_var_GRIB_names(17)        = "sr"
      Met_var_WMO_names(17)         = "SFCR"
      Met_var_GRIB2_DPcPnSt(17,1:4) = (/2, 0, 1, 1/)
      Met_var_GRIB1_Param(17)       = 0
      Met_var_GRIB1_St(17)          = ""
      Met_var_ndim(17)              = 3
        ! Wind gust speed  (m/s)
      Met_var_NC_names(18)          = "Wind_speed_gust_surface"
      Met_var_GRIB_names(18)        = "gust"
      Met_var_WMO_names(18)         = "GUST"
      Met_var_GRIB2_DPcPnSt(18,1:4) = (/0, 2, 22, 1/)
      Met_var_GRIB1_Param(18)       = 0
      Met_var_GRIB1_St(18)          = ""
      Met_var_ndim(18)              = 3
        ! Surface temperature  (K)
      Met_var_NC_names(19)          = "Temperature_surface"
      Met_var_WMO_names(19)         = ""
      Met_var_GRIB2_DPcPnSt(19,1:4) = (/0, 0, 0, 0/)
      Met_var_GRIB1_Param(19)       = 11
      Met_var_GRIB1_St(19)          = "1"
      Met_var_ndim(19)              = 3
      !--------------------------------
      ! Atmospheric Structure
      !--------------------------------
        ! Pressure at base of lower cloud level  (Pa)
      Met_var_NC_names(20)          = "Pressure_cloud_base"
      Met_var_GRIB_names(20)        = "pres"
      Met_var_WMO_names(20)         = "PRES"
      Met_var_GRIB2_DPcPnSt(20,1:4) = (/0, 3, 0, 2/)
      Met_var_GRIB1_Param(20)       = 1
      Met_var_GRIB1_St(20)          = "2"
      Met_var_ndim(20)              = 3
        ! Pressure at top of lower cloud level  (Pa)
      Met_var_NC_names(21)          = "Pressure_cloud_tops"
      Met_var_GRIB_names(21)        = "pres"
      Met_var_WMO_names(21)         = "PRES"
      Met_var_GRIB2_DPcPnSt(21,1:4) = (/0, 3, 0, 3/)
      Met_var_GRIB1_Param(21)       = 1
      Met_var_GRIB1_St(21)          = "3"
      Met_var_ndim(21)              = 3
        ! Temperature at the top of lower cloud level  (K)
      Met_var_NC_names(22)          = "Temperature_cloud_tops"
      Met_var_GRIB_names(22)        = "t"
      Met_var_WMO_names(22)         = "TMP"
      Met_var_GRIB2_DPcPnSt(22,1:4) = (/0, 0, 0, 3/)
      Met_var_GRIB1_Param(22)       = 11
      Met_var_GRIB1_St(22)          = "3"
      Met_var_ndim(22)              = 3
        ! Total cloud cover  (%)
      Met_var_NC_names(23)          = "Total_cloud_cover_entire_atmosphere"
      Met_var_GRIB_names(23)        = "tcc"
      Met_var_WMO_names(23)         = "TCDC"
      Met_var_GRIB2_DPcPnSt(23,1:4) = (/0, 6, 1, 200/)
      Met_var_GRIB1_Param(23)       = 71
      Met_var_GRIB1_St(23)          = "200"
      Met_var_ndim(23)              = 3
        ! Cloud cover of lower cloud level  (%)
      Met_var_NC_names(24)          = "Low_cloud_cover_low_cloud"
      Met_var_GRIB_names(24)        = "lcc"
      Met_var_WMO_names(24)         = "LCDC"
      Met_var_GRIB2_DPcPnSt(24,1:4) = (/0, 6, 3, 214/)
      Met_var_GRIB1_Param(24)       = 0
      Met_var_GRIB1_St(24)          = ""
      Met_var_ndim(24)              = 3
      !--------------------------------
      ! Moisture
      !--------------------------------
        ! Relative humidity  (%)
      Met_var_NC_names(30)          = "Relative_humidity_isobaric"
      Met_var_GRIB_names(30)        = "r"
      Met_var_WMO_names(30)         = "RH"
      Met_var_GRIB2_DPcPnSt(30,1:4) = (/0, 1, 1, 100/)
      Met_var_GRIB1_Param(30)       = 52
      Met_var_GRIB1_St(30)          = "pl" 
      Met_var_ndim(30)              = 4
        ! Specific humidity  (kg/kg)
      Met_var_NC_names(31)          = "Specific_humidity_isobaric"
      Met_var_GRIB_names(31)        = "q"
      Met_var_WMO_names(31)         = "SPFH"
      Met_var_GRIB2_DPcPnSt(31,1:4) = (/0, 1, 0, 100/)
      Met_var_GRIB1_Param(31)       = 51
      Met_var_GRIB1_St(31)          = "pl"
      Met_var_ndim(31)              = 4
        ! Cloud water mixing ratio  (kg/kg)
      Met_var_NC_names(32)          = "Cloud_mixing_ratio_isobaric"
      Met_var_GRIB_names(32)        = "clwmr"
      Met_var_WMO_names(32)         = "CLWMR"
      Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
      Met_var_ndim(32)              = 4
        ! Snow mixing ratio  (kg/kg)
      Met_var_NC_names(33)          = "Snow_mixing_ratio_isobaric"
      Met_var_GRIB_names(33)        = "snmr"
      Met_var_WMO_names(33)         = "SNMR"
      Met_var_GRIB2_DPcPnSt(33,1:4) = (/0, 1, 25, 100/)
      Met_var_GRIB1_Param(33)       = 0
      Met_var_GRIB1_St(33)          = ""
      Met_var_ndim(33)              = 4
      !--------------------------------
      ! Precipitation
      !--------------------------------
        ! Categorical rain at ground surface  (0/1 no/yes)
      Met_var_NC_names(40)          = "Categorical_Rain_surface"
      Met_var_GRIB_names(40)        = "crain"
      Met_var_WMO_names(40)         = "CRAIN"
      Met_var_GRIB2_DPcPnSt(40,1:4) = (/0, 1, 192, 1/)
      Met_var_GRIB1_Param(40)       = 140
      Met_var_GRIB1_St(40)          = "1"
      Met_var_ndim(40)              = 3
        ! Categorical snow at ground surface  (0/1 no/yes)
      Met_var_NC_names(41)          = "Categorical_Snow_surface"
      Met_var_GRIB_names(41)        = "csnow"
      Met_var_WMO_names(41)         = "CSNOW"
      Met_var_GRIB2_DPcPnSt(41,1:4) = (/0, 1, 195, 1/)
      Met_var_GRIB1_Param(41)       = 143
      Met_var_GRIB1_St(41)          = "1"
      Met_var_ndim(41)              = 3
        ! Categorical freezing rain at ground surface  (0/1 no/yes)
      Met_var_NC_names(42)          = "Categorical_Freezing_Rain_surface"
      Met_var_GRIB_names(42)        = "cfrzr"
      Met_var_WMO_names(42)         = "CFRZR"
      Met_var_GRIB2_DPcPnSt(42,1:4) = (/0, 1, 193, 1/)
      Met_var_GRIB1_Param(42)       = 141
      Met_var_GRIB1_St(42)          = "1"
      Met_var_ndim(42)              = 3
        ! Categorical ice pellets at ground surface  (0/1 no/yes)
      Met_var_NC_names(43)          = "Categorical_Ice_Pellets_surface"
      Met_var_GRIB_names(43)        = "cicep"
      Met_var_WMO_names(43)         = "CICEP"
      Met_var_GRIB2_DPcPnSt(43,1:4) = (/0, 1, 194, 1/)
      Met_var_GRIB1_Param(43)       = 142
      Met_var_GRIB1_St(43)          = "1"
      Met_var_ndim(43)              = 3
        ! Precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(44)          = "Precipitation_rate_surface"
      Met_var_GRIB_names(44)        = "prate"
      Met_var_WMO_names(44)         = "PRATE"
      Met_var_GRIB2_DPcPnSt(44,1:4) = (/0, 1, 7, 1/)
      Met_var_GRIB1_Param(44)       = 59
      Met_var_GRIB1_St(44)          = "1"
      Met_var_ndim(44)              = 3
        ! Convective liquid precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(45)          = "Precip.rate convective  (liquid)"
      Met_var_GRIB_names(45)        = "cprat"
      Met_var_GRIB2_DPcPnSt(45,1:4) = (/0, 1, 196, 1/)
      Met_var_WMO_names(45)         = "CPRAT"
      Met_var_GRIB1_Param(45)       = 0
      Met_var_GRIB1_St(45)          = "1"
      Met_var_ndim(45)              = 3
        ! Large-scale precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(46)          = "Precip.rate large-scale (ice)"
      Met_var_WMO_names(46)         = ""
      Met_var_GRIB1_Param(46)       = 0
      Met_var_GRIB1_St(46)          = "1"
      Met_var_ndim(46)              = 3
        ! Convective frozen precipitation rate at surface for (kg/m2/s)
      Met_var_NC_names(47)          = "Precip.rate convective  (ice)"
      Met_var_WMO_names(47)         = ""
      Met_var_GRIB1_Param(47)       = 0
      Met_var_GRIB1_St(47)          = "1"
      Met_var_ndim(47)              = 3

      if (MR_iwindformat.eq.0) then
          ! Custom format based on template
        MR_Reannalysis = .false.
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "Custom format based on template"
        endif;enddo
          ! This expects that MR_iwf_template has been filled by the calling program
        call MR_Read_Met_Template

      elseif (MR_iwindformat.eq.1) then
          ! ASCII profile
        MR_Reannalysis = .false.
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "ASCII profile"
        endif;enddo
      elseif (MR_iwindformat.eq.2) then
          ! Radiosonde data
        MR_Reannalysis = .false.
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "Radiosonde data"
          write(outlog(io),*)"  Number of radiosonde stations = ",MR_iGridCode
        endif;enddo
      elseif (MR_iwindformat.eq.3) then 
        ! NARR3D NAM221 32 km North America files
          ! https://rda.ucar.edu/datasets/ds608.0/
          !   merged_AWIP32.2018062000(.nc)
          ! Note that winds are "earth-relative" and must be rotated!
          !   See  http://www.emc.ncep.noaa.gov/mmb/rreanl/faq.html#eta-winds

        if(MR_iversion.eq.-1)MR_iversion = 0 ! the lastest version of NARR data is v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NARR3D NAM221 32 km North America files with Re=6367.47"
        endif;enddo

        MR_iGridCode = 1221  ! This is almost NAM221, but uses a diff Re
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode) 
        MR_Reannalysis = .true.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.  ! Geopotential Height
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="u_wind_isobaric"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="v_wind_isobaric"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Pressure_vertical_velocity_isobaric"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="Temp_isobaric"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Pressure_vertical_velocity_isobaric"
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="Planetary_boundary_layer_height_surface"
        Met_var_IsAvailable(11)=.true.; Met_var_NC_names(11)="u_wind_height_above_ground"
        Met_var_IsAvailable(12)=.true.; Met_var_NC_names(12)="v_wind_height_above_ground"
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="Surface_friction_velocity_surface"
        Met_var_IsAvailable(15)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        ! Moisture
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="Cloud_water_isobaric"
        Met_var_IsAvailable(33)=.true.; Met_var_NC_names(33)="Ice_mixing_ratio_isobaric"
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_rain_yes1_no0_surface"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_snow_yes1_no0_surface"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_freezing_rain_yes1_no0_surface"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_ice_pellets_yes1_no0_surface"
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp
      elseif (MR_iwindformat.eq.4) then
        ! NAM Regional North America 221 32 km North America files
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !   nam.t00z.awip3200.tm00.grib2

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM221 32 km North America files with Re=6371.229"
        endif;enddo

        MR_iGridCode = 221
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(45)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.5) then
        ! NAM216 AK 45km
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !  nam.t00z.awipak00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM216 AK 45km"
        endif;enddo

        MR_iGridCode = 216
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="Planetary_Boundary_Layer_Height"
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="Frictional_Velocity"
        Met_var_IsAvailable(14)=.true.; Met_var_NC_names(16)="Volumetric_Soil_Moisture_Content"
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover"
        ! Moisture
        Met_var_IsAvailable(30)=.true.;
        Met_var_IsAvailable(31)=.true.;
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Large_scale_precipitation_non-convective"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.6) then
        ! NAM Regional 90 km grid 104
        !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod
        !    nam.t00z.grbgrd00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM Regional 90 km grid 104"
        endif;enddo

        MR_iGridCode = 104
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover_entire_atmosphere_0_Hour_Average"
        Met_var_IsAvailable(24)=.true.; Met_var_NC_names(24)="Convective_cloud_cover_entire_atmosphere_0_Hour_Average"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Large_scale_precipitation_non-convective_surface_0_Hour_Accumulation"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.7) then
        ! CONUS 212 40km
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !    nam.t00z.awp21100.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "CONUS 212 40km"
        endif;enddo

        MR_iGridCode = 212
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover"
        Met_var_IsAvailable(24)=.true.; Met_var_NC_names(24)="Convective_cloud_cover"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate_surface"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.8) then
        ! CONUS 218 (12km)
          ! ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !   nam.t00z.awphys00.grb2.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "CONUS 218 (12km)"
        endif;enddo

        MR_iGridCode = 218
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.9) then
        ! CONUS 227 (5.08 km)
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !    nam.t00z.conusnest.hiresf00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "CONUS 227 (5.1 km)"
        endif;enddo

        MR_iGridCode = 227
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.10)then
        ! NAM 242 11.25 km AK
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !    nam.t00z.awipak00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM 242 11.25 km AK"
        endif;enddo

        MR_iGridCode = 242
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Total_precipitation_surface_0_Hour_Accumulation"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.11)then
        ! NAM 196 2.5 km HI
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !    nam.t00z.hawaiinest.hiresf00.tm0

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    " NAM 196 2.5 km HI"
        endif;enddo

        MR_iGridCode = 196
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.; Met_var_NC_names(15)="Snow_Cover_surface"
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="Soil_moisture_content_depth_below_surface_layer"
        Met_var_IsAvailable(17)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.12)then
        ! NAM 198 5.953 km AK
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/ (this is actually now 91)
          !    nam.t00z.alaskanest.hiresf00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM 198 5.953 km AK"
        endif;enddo

        MR_iGridCode = 198
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="Soil_moisture_content_depth_below_surface_layer"
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.13)then  ! NAM 91 2.976 km AK
        ! NAM 91 2.976 km AK
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/ (this is actually now 91)
          !    nam.t00z.alaskanest.hiresf00.tm00
          !
          ! Note: the dimension names given below are those generated by netcdf-java 4.5
          !       acting on the truncated GRIB files generated by get_nam91.sh which 
          !       uses get_inv.pl to get just the GRIB layers needed.
          !       This is relavent because the numbering of the isobaric dimensions
          !       appears to be in the order they are processed in the GRIB file.

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM 91 2.976 km AK"
        endif;enddo

        MR_iGridCode = 91
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.14) then
        ! CONUS 1227 (3.0 km)
          !  
          !  nam.t00z.conusnest.hiresf00.tm00.grib2.nc

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "CONUS no grid ID (3.0 km)"
        endif;enddo

        MR_iGridCode = 1227
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

!      elseif (MR_iwindformat.eq.15) then
!        ! NAM Caribbean 181 (0.108 degrees)
!          !  
!          !  nam.t00z.afwaca00.tm00.grib2.nc
!
!        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
!        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
!          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
!                    "NAM Caribbean 0.108 degree"
!        endif;enddo
!
!        MR_iGridCode = 181
!        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
!        MR_Reannalysis = .false.
!
!        ! Mechanical / State variables
!        Met_var_IsAvailable(1)=.true.
!        Met_var_IsAvailable(2)=.true.
!        Met_var_IsAvailable(3)=.true.
!        Met_var_IsAvailable(4)=.true.
!        Met_var_IsAvailable(5)=.true.
!        Met_var_IsAvailable(7)=.true.
!        ! Surface
!        Met_var_IsAvailable(10)=.true.
!        Met_var_IsAvailable(11)=.true.
!        Met_var_IsAvailable(12)=.true.
!        Met_var_IsAvailable(13)=.true.
!        Met_var_IsAvailable(15)=.true.
!        Met_var_IsAvailable(16)=.true.
!        !Met_var_IsAvailable(17)=.true.
!        !Met_var_IsAvailable(18)=.true.
!        ! Atmospheric Structure
!        Met_var_IsAvailable(20)=.true.
!        Met_var_IsAvailable(21)=.true.
!        Met_var_IsAvailable(23)=.true.
!        ! Moisture
!        Met_var_IsAvailable(30)=.true.
!        Met_var_IsAvailable(31)=.true.
!        Met_var_IsAvailable(32)=.true.
!        Met_var_IsAvailable(33)=.true.
!        ! Precipitation
!        Met_var_IsAvailable(40)=.true.
!        Met_var_IsAvailable(41)=.true.
!        Met_var_IsAvailable(42)=.true.
!        Met_var_IsAvailable(43)=.true.
!        !Met_var_IsAvailable(44)=.true.
!
!        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.20)then
        ! GFS 0.5 deg
          !  http://www.nco.ncep.noaa.gov/pmb/products/gfs/
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/gfs/prod/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "GFS 0.5"
        endif;enddo

        MR_iGridCode = 4
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        !Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        !Met_var_IsAvailable(11)=.true.
        !Met_var_IsAvailable(12)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        !Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        !Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        !Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        !Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        !Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Precipitation_rate"
        !Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.21)then
        ! GFS 1.0 deg
          ! http://www.nco.ncep.noaa.gov/pmb/products/gfs/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "GFS 1.0-degree"
        endif;enddo

        MR_iGridCode = 3
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(32)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.22)then
        ! GFS 0.25 deg
          ! http://www.nco.ncep.noaa.gov/pmb/products/gfs/
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/gfs/prod/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "GFS 0.25"
        endif;enddo

        MR_iGridCode = 193
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Precipitation_rate"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.23)then
         ! NCEP / DOE reanalysis 2.5 degree files 
         ! https://rda.ucar.edu/datasets/ds091.0

        if(MR_iversion.eq.-1)MR_iversion = 0 ! latest is 0, but deprecated
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NCEP / DOE reanalysis 2.5 degree files"
        endif;enddo

        MR_iGridCode = 2
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 132

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        ! Precipitation
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.24)then
         ! NASA-MERRA reanalysis 0.625 x 0.5 degree files 

        if(MR_iversion.eq.-1)MR_iversion = 2 ! v1 finished in 2016
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NASA-MERRA-2 reanalysis 0.625/0.5 degree files"
        endif;enddo

        MR_iGridCode = 1024
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        !   Available in MERRA2_400.inst3_3d_asm_Np.YYYYMMDD.nc4
        !    from https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"
        ! Moisture
        !   Available in MERRA2_400.inst3_3d_asm_Np.YYYYMMDD.nc4
        !    from https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4
        Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="RH"          ! float percent
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="QV"          ! float cloud liquid water mixing ratio kg/kg
        Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="QL"          ! float cloud liquid water mixing ratio kg/kg
        Met_var_IsAvailable(33)=.true.; Met_var_NC_names(33)="QI"          ! float cloud ice mixing ratio kg/kg

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.25)then
         ! NCEP/NCAR reanalysis 2.5 degree files 

        if(MR_iversion.eq.-1)MR_iversion = 0 !
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                   "NCEP/NCAR reanalysis 2.5 degree files"
        endif;enddo

        MR_iGridCode = 2
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        if(MR_iwind.eq.4)then
           ! https://rda.ucar.edu/datasets/ds090.0

          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.
          Met_var_IsAvailable(2)=.true.
          Met_var_IsAvailable(3)=.true.
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Pressure_vertical_velocity_isobaric"
          Met_var_IsAvailable(5)=.true.
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Pressure_vertical_velocity_isobaric"
          ! Moisture
          Met_var_IsAvailable(30)=.true.
  
          fill_value_sp = -9999.0_sp

        elseif(MR_iwind.eq.5)then
          MR_iw5_hours_per_file = 8760.0_dp

          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="hgt"        ! short m^2/s^2 (32066.f,1.f)
          Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="uwnd"       ! short m/s (202.66f,0.01f)
          Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="vwnd"       ! short m/s (202.66f,0.01f)
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="omega"      ! short Pa/s (29.765f,0.001f)
          Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="air"        ! short K (477.66f,0.01f)
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="omega"      ! short Pa/s (29.765f,0.001f)
          ! Moisture
          !Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="rhum"      ! short  (302.66f,0.01f)
          Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="shum"      ! short SpecHum ~ mixing ratio kg/kg(0.032666f,1.e-06f)
          !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="shum"      ! short should really be QL (liquid)
  
          fill_value_sp = -9999.0_sp

        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : NCEP Reannalysis provided only as iwind 4 or iwind 5"
          endif;enddo
          stop 1
        endif

      elseif (MR_iwindformat.eq.26)then
         ! JRA-55 reanalysis 1.25 degree files 
         ! https://rda.ucar.edu/datasets/ds628.0/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! 
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "JRA-55 reanalysis 1.25 degree files"
        endif;enddo

        MR_iGridCode = 45
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.
        MR_iw5_hours_per_file = 672.0_dp

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 200

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="HGT_GDS0_ISBL"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="UGRD_GDS0_ISBL"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="VGRD_GDS0_ISBL"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="VVEL_GDS0_ISBL"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="TMP_GDS0_ISBL"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="VVEL_GDS0_ISBL"
        ! Moisture
        Met_var_IsAvailable(31)=.true.; Met_var_GRIB_names(31)= "SPFH_GDS0_ISBL"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.27)then
         ! NOAA-CIRES reanalysis 2.0 degree files 
         ! https://rda.ucar.edu/datasets/ds131.2/
         ! https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.pressure.html

        if(MR_iversion.eq.-1)then
          MR_iversion = 3 ! This is the latest
        endif
        if(MR_iversion.eq.2)then

          do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
            write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                      "NOAA-CIRES reanalysis 2.0 degree files"
          endif;enddo

          MR_iGridCode = 1027
          call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
          MR_Reannalysis = .true.
          MR_iw5_hours_per_file = 8760.0_dp
  
          Met_var_GRIB1_Table(1:MR_MAXVARS) = 2
  
          if(MR_Use_RDA)then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Version 2c (https://rda.ucar.edu/datasets/ds131.2/)
            ! Note: these must be converted from GRIB using
            !   ncl_convert2nc pgrbanl_mean_1912_VVEL_pres.grib -L
            ! netcdf-java does not seem to work on these
    
            ! Momentum / State variables
            Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="HGT_GDS0_ISBL"
            Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_GRD_GDS0_ISBL"
            Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_GRD_GDS0_ISBL"
            Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="V_VEL_GDS0_ISBL"
            Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="TMP_GDS0_ISBL"
            Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="V_VEL_GDS0_ISBL"
          else 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Version 2c (https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.pressure.html)
    
            ! Momentum / State variables
            Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="hgt"
            Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="uwnd"
            Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="vwnd"
            Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="omega"
            Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="air"
            Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="omega"
          endif
          fill_value_sp = 1.0e+20_sp
        elseif(MR_iversion.eq.3)then
          ! NOAA-CIRES-DOE Twentieth Century Reanalysis Version 3 
          ! https://rda.ucar.edu/datasets/ds131.3/
          ! Note: this is also available from
          ! https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html, but might
          ! need format verification
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
            write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                      "NOAA-CIRES-DOE reanalysis v3"
          endif;enddo

          MR_iGridCode = 1027
          call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
          MR_Reannalysis = .true.
          MR_iw5_hours_per_file = 8760.0_dp
  
          Met_var_GRIB1_Table(1:MR_MAXVARS) = 2
  
          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="gh"
          Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="u"
          Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="v"
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="w"
          Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="t"
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="w"
          fill_value_sp = 9999.0_dp
        endif

      elseif (MR_iwindformat.eq.28)then
         ! ECMWF Interim Reanalysis (ERA-Interim)
         ! https://rda.ucar.edu/datasets/ds627.0

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "ECMWF Interim Reanalysis (ERA-Interim)"
        endif;enddo

        MR_iGridCode = 170
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Geopotential_isobaric"
                                       Met_var_GRIB1_Param(1)=129
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_component_of_wind_isobaric"
                                       Met_var_GRIB1_Param(2)=131
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_component_of_wind_isobaric"
                                       Met_var_GRIB1_Param(3)=132
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"
                                       Met_var_GRIB1_Param(4)=135
        Met_var_IsAvailable(5)=.true.
                                       Met_var_GRIB1_Param(5)=130
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"
                                       Met_var_GRIB1_Param(7)=135
        Met_var_IsAvailable(30)=.true.
                                       Met_var_GRIB1_Param(30)=157
        Met_var_IsAvailable(31)=.true.
                                       Met_var_GRIB1_Param(31)=133

        fill_value_sp = -9999.0_sp
        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.29)then
         ! ECMWF ERA5
         ! This is a premier reanalysis product that is available in lots of formats
         ! From RDA, data are availabe in complete grid files with one variable/file/day
         !   Current version is at https://rda.ucar.edu/datasets/ds633.0 and has a
         !    file name format e5.oper.an.pl.128_129_z.ll025sc.YYYYMMDD00_YYYYMMMDD23.nc
         !   Previous version is at https://rda.ucar.edu/datasets/ds630.0/ and has a
         !    file name format e5.oper.an.pl.128_129_z.regn320sc.YYYYMMDD00_YYYYMMMDD23.nc
         !   These are available as GRIB or nc4 files, but variable names will change
         !   if GRIB is converted to NetCDF via netcdf-java or ncl_convert2nc
         !   Using RDA files required iwf=5 since variables on stored in different files
         ! From CDS, data are available for download as a subset of the full dataset via
         !  https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form
         ! Data are provided as a single file for all variables and timesteps (iwf=4) as
         ! both GRIB and nc. If GRIB is downloaded, then converted to nc via netcdf-java or
         ! ncl_convert2nc, variable names will be different from that of the files provided
         ! directly by CDS, though metreader will attempt to translate.

        ! Dimension names are in the direct nc4 download
        !Met_dim_names(1)="time"
        !Met_dim_names(2)="level"
        !Met_dim_names(3)="latitude"
        !Met_dim_names(4)="longitude"
        ! Dimension names are in the GRIB to NetCDF via ncj
        !Met_dim_names(1)="time"
        !Met_dim_names(2)="isobaric"
        !Met_dim_names(3)="lat"
        !Met_dim_names(4)="lon"
        ! Dimension names are in the GRIB to NetCDF via ncl_convert2nc
        !Met_dim_names(1)="initial_time0_hours"
        !Met_dim_names(2)="lv_ISBL1"
        !Met_dim_names(3)="g0_lat_2"
        !Met_dim_names(4)="g0_lon_3"

        if(MR_iversion.eq.-1)MR_iversion = 2 ! v1 finished in 2019
        if(MR_iversion.eq.1)then
          MR_RDAcode=630
        elseif(MR_iversion.eq.2)then
          MR_RDAcode=633
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(outlog(io),*)" MR ERROR: ERA5 version not recognized ",MR_iversion
          endif;enddo
          stop 1
        endif

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "ECMWF ERA5 reanalysis"
        endif;enddo

        MR_iGridCode = 1029
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128

        if(MR_iwind.eq.3.or.MR_iwind.eq.4)then
          ! These data are from CDS
          ! https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form
          !                                                       ! via ncl     via netcdf-java
          Met_var_IsAvailable( 1)=.true.; Met_var_NC_names( 1)="z" ! Z_GDS0_ISBL Geopotential_isobaric
          Met_var_IsAvailable( 2)=.true.; Met_var_NC_names( 2)="u" ! Z_GDS0_ISBL U_component_of_wind_isobaric
          Met_var_IsAvailable( 3)=.true.; Met_var_NC_names( 3)="v" ! Z_GDS0_ISBL V_component_of_wind_isobaric
          Met_var_IsAvailable( 4)=.true.; Met_var_NC_names( 4)="w" ! Z_GDS0_ISBL Vertical_velocity_isobaric
          Met_var_IsAvailable( 5)=.true.; Met_var_NC_names( 5)="t" ! Z_GDS0_ISBL Temperature_isobaric
          Met_var_IsAvailable( 7)=.true.; Met_var_NC_names( 7)="w" ! Z_GDS0_ISBL Vertical_velocity_isobaric
          Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="q" ! Z_GDS0_ISBL Specific_humidity_isobaric
        elseif(MR_iwind.eq.5)then
          ! RDA v1 or v2
          ! https://rda.ucar.edu/datasets/ds633.0
          MR_Use_RDA = .true.
          ! Note: files are provided as one variable per file and must use iwind=5
          !       Each file contains one day of data
          MR_iw5_hours_per_file = 24.0_dp
          !                                                       ! via ncl     via netcdf-java
          Met_var_IsAvailable( 1)=.true.; Met_var_NC_names( 1)="Z" ! Z_GDS0_ISBL Geopotential_isobaric
          Met_var_IsAvailable( 2)=.true.; Met_var_NC_names( 2)="U" ! U_GDS0_ISBL U_component_of_wind_isobaric
          Met_var_IsAvailable( 3)=.true.; Met_var_NC_names( 3)="V" ! V_GDS0_ISBL V_component_of_wind_isobaric
          Met_var_IsAvailable( 4)=.true.; Met_var_NC_names( 4)="W" ! W_GDS0_ISBL Vertical_velocity_isobaric
          Met_var_IsAvailable( 5)=.true.; Met_var_NC_names( 5)="T" ! T_GDS0_ISBL Temperature_isobaric
          Met_var_IsAvailable( 7)=.true.; Met_var_NC_names( 7)="W" ! W_GDS0_ISBL Vertical_velocity_isobaric
          Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="q" ! Q_GDS0_ISBL Specific_humidity_isobaric

        endif

        fill_value_sp = -9999.0_sp
        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.30)then
         ! ECMWF ERA 20C
         ! https://rda.ucar.edu/datasets/ds626.0
         ! Note: files are provided as one variable per file

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "ECMWF ERA20C reanalysis"
        endif;enddo

        MR_iGridCode = 1030
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.
        MR_iw5_hours_per_file = 672.0_dp

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Geopotential_isobaric"        ! e5.oper.an.pl.128_129_z.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_component_of_wind_isobaric" ! e5.oper.an.pl.128_129_u.ll025uv.1991061500_1991061523.nc
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_component_of_wind_isobaric" ! e5.oper.an.pl.128_129_v.ll025uv.1991061500_1991061523.nc
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"   ! e5.oper.an.pl.128_129_w.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="Temperature_isobaric"         ! e5.oper.an.pl.128_129_t.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"   ! e5.oper.an.pl.128_129_w.ll025sc.19910


        fill_value_sp = -9999.0_sp

        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.31) then
        ! NAM Caribbean 181 (0.108 degrees)
          !  
          !  nam.t00z.afwaca00.tm00.grib2.nc

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NAM Caribbean 0.108 degree"
        endif;enddo

        MR_iGridCode = 181
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        !Met_var_IsAvailable(17)=.true.
        !Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        !Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.32)then
         ! Air Force Weather Agency subcenter = 0
         ! GALWEM

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "Air Force Weather Agency subcenter = 0"
        endif;enddo

        MR_iGridCode = 1032
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.33)then
         ! CCSM3.0 Community Atmosphere Model (CAM)
         ! http://www.cesm.ucar.edu/models/atm-cam/
         ! peleoclimate monthly averages

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "CCSM3.0 Community Atmosphere Model (CAM)"
        endif;enddo

        MR_iGridCode = 1033
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_dim_fac(1) = 24.0

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Z3"      ! float m
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"
        ! Moisture
        Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="RELHUM"
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="Q"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.34)then
         ! ECMWF 0.25 degree forecast
         ! Amazon's AWS: "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"
         ! ECMWF:  https://data.ecmwf.int/forecasts

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "ECMWF "
        endif;enddo

        MR_iGridCode = 193
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.

        ! Momentum / State variables
!        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Geopotential_isobaric"        ! e5.oper.an.pl.128_129_z.ll025sc.1991061500_1991061523.nc
!        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_component_of_wind_isobaric" ! e5.oper.an.pl.128_129_u.ll025uv.1991061500_1991061523.nc
!        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_component_of_wind_isobaric" ! e5.oper.an.pl.128_129_v.ll025uv.1991061500_1991061523.nc
!        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"   ! e5.oper.an.pl.128_129_w.ll025sc.1991061500_1991061523.nc
!        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="Temperature_isobaric"         ! e5.oper.an.pl.128_129_t.ll025sc.1991061500_1991061523.nc
!        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"   ! e5.oper.an.pl.128_129_w.ll025sc.19910

        fill_value_sp = -9999.0_sp

!        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.40)then
         ! NASA-GEOS Cp

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NASA-GEOS Cp"
        endif;enddo

        MR_iGridCode = 1040
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.41)then
         ! NASA-GEOS Np

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "NASA-GEOS Np"
        endif;enddo

        MR_iGridCode = 1041
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.50)then
         ! WRF - output

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
          write(outlog(io),*)"  NWP format to be used = ",MR_iwindformat,&
                    "WRF"
        endif;enddo

        MR_iGridCode = 1050
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

!        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "Time"
!        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "bottom_top"
!        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "south_north"
!        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "west_east"
        ! for pressure, read "P"  :: perturbation pressure
        !               and  "PB" :: base pressure

        ! for geopotential, read "PH"  :: perturbation geopotential
        !                   and  "PHB" :: base-state geopotential

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="PHB"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="W"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"      ! float K perturbation potential temperature (theta-t0)
        Met_var_IsAvailable(6)=.true.; Met_var_NC_names(6)="PB"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="W"

        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="PBLH"
        Met_var_IsAvailable(11)=.true.; Met_var_NC_names(11)="U10"
                                        Met_var_ndim(11)=3
        Met_var_IsAvailable(12)=.true.; Met_var_NC_names(12)="V10"
                                        Met_var_ndim(12)=3
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="UST"
        Met_var_IsAvailable(15)=.true.; Met_var_NC_names(15)="SNOWH"
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="SMOIS" !Soil moisture m3 m-3
        ! Moisture
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="QVAPOR" !QV (specific humidity)
        ! Precipitation
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="RAINC"   ! ACCUMULATED TOTAL CUMULUS PRECIPITATION in mm
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="RAINNC"  ! ACCUMULATED TOTAL GRID SCALE PRECIPITATION in mm

        fill_value_sp = 1.0e20_sp

        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then     
          write(errlog(io),*)'MR ERROR : MR_iwindformat not supported'
          write(errlog(io),*)'           MR_iwindformat=',MR_iwindformat,&
                                       '. Program stopped in MetReader.f90'
        endif;enddo
        stop 1
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then     
        write(outlog(io),*)"                grid ID = ",MR_iGridCode
      endif;enddo
      select case (MR_idataFormat)
      case(1)
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then       
          write(outlog(io),*)"            data format = ",MR_idataFormat,"ASCII"
        endif;enddo
      case(2)
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then        
          write(outlog(io),*)"            data format = ",MR_idataFormat,"NETCDF"
        endif;enddo
#ifndef USENETCDF
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Met files are NetCDF, but MetReader was not"
          write(errlog(io),*)"          compiled with NetCDF support."
          write(errlog(io),*)"          Please recompile MetReader with NetCDF support"
        endif;enddo
        stop 1
#endif
      case(3)
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then        
          write(outlog(io),*)"            data format = ",MR_idataFormat,"GRIB"
        endif;enddo
#ifndef USEGRIB
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Met files are GRIB, but MetReader was not"
          write(errlog(io),*)"          compiled with GRIB support."
          write(errlog(io),*)"          Please recompile MetReader with GRIB support"
        endif;enddo
        stop 1
#endif
      case default
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  MR_idataFormat not 1:4."
          write(errlog(io),*)"           Exiting in MR_Allocate_FullMetFileList"
        endif;enddo
        stop 1
      end select

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then         
        write(outlog(io),*)"     Allocating space for ",MR_iwindfiles,"file(s)."
      endif;enddo

      if (MR_iwind.eq.5)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then              
          write(outlog(io),*)"For iwf=5:: one variable per file."
          write(outlog(io),*)"Only the directory should be listed."
          write(outlog(io),*)"The remaining path is hard-coded."
        endif;enddo
        ! Reset MR_iwindfiles to 2: only one "file" will be read,
        ! but we need this to be 2 to accommodate runs that might span two years/months
        ! For ERA5 files, data is provided in daily files so MR_iwindfiles might be > 2
        if(MR_Comp_Time_in_hours.gt.0.0)then
          if(MR_iwindformat.eq.25)then  ! NCEP 2.5 degree
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1  ! These are yearly files
          elseif(MR_iwindformat.eq.26)then ! 
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! Monthly files
          elseif(MR_iwindformat.eq.27)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! yearly files
          elseif(MR_iwindformat.eq.29)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! daily files
          elseif(MR_iwindformat.eq.30)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! monthly
          endif
        else
          MR_iwindfiles = 2
        endif
      endif

      allocate (MR_windfiles(MR_iwindfiles))
      do i=1,MR_iwindfiles
        write(MR_windfiles(i),'(130x)')
      enddo
      allocate (MR_windfiles_nt_fullmet(MR_iwindfiles))
      MR_windfiles_nt_fullmet(:)=0
      if(MR_idataFormat.eq.3)then
        allocate (MR_windfiles_Have_GRIB_index(MR_iwindfiles))
          ! This will be reset to true if the index files are found
        MR_windfiles_Have_GRIB_index = .false.
        allocate (MR_windfiles_GRIB_index(MR_iwindfiles))
        do i=1,MR_iwindfiles
          write(MR_windfiles_GRIB_index(i),'(130x)')
        enddo
      endif

      if(MR_iwind.eq.1)then
        ! For the 1d profile or radiosonde case, igrid is used for the number of sonde
        ! locations.  If it is not provided, set it to one
        if(igrid.eq.0)then
          MR_iGridCode = 1
          MR_nSnd_Locs = MR_iGridCode
        else
          MR_iGridCode = igrid
          MR_nSnd_Locs = MR_iGridCode
        endif
          ! Now make sure that the number of windfiles is a multiple of the number of locations
          ! since this will be the number of timesteps
        if(mod(MR_iwindfiles,MR_nSnd_Locs).eq.0)then
          MR_Snd_nt_fullmet = MR_iwindfiles / MR_nSnd_Locs
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  The grid code for 1d ASCII input is interpreted to be"
            write(errlog(io),*)"           the number of sonde locations.  Each group of sondes can"
            write(errlog(io),*)"           be repeated, correspoding to multiple timesteps.  The"
            write(errlog(io),*)"           number of windfiles must be a multiple of the number of"
            write(errlog(io),*)"           locations"
            write(errlog(io),*)"                   MR_iwind = ",MR_iwind
            write(errlog(io),*)"             MR_iwindformat = ",MR_iwindformat
            write(errlog(io),*)"               MR_iGridCode = ",MR_iGridCode
            write(errlog(io),*)"             MR_idataFormat = ",MR_idataFormat
            write(errlog(io),*)"               MR_nSnd_Locs = ",MR_iGridCode
            write(errlog(io),*)"          MR_Snd_nt_fullmet = ",real(MR_iwindfiles)/real(MR_nSnd_Locs)
          endif;enddo
          stop 1
        endif
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      CALLED_MR_Allocate_FullMetFileList = .true.

      return

      end subroutine MR_Allocate_FullMetFileList


!##############################################################################
!
!     MR_Read_Met_DimVars
!
!     This subroutine expects that MR_windfiles has been filled by the calling
!     program, and checks for their existence.
!
!     Depending on the NWP data format, a subroutine is called to read the spatial
!     structure of the NWP files, followed by one for the temporal span
!      e.g. MR_Read_Met_DimVars_netcdf
!           MR_Read_Met_Times_netcdf
!
!     From the calling program, this is called once the names of the NWP files
!     are specified.  If a custom netcdf template file is to be used (iwf=0), then
!     the variable MR_iwf_template must also be filled so that it can be read
!     from MR_Read_Met_DimVars_[].  
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!     The next step is for the calling program to specify the projection parameters
!     of the computational grid (i.e. the grid MetReader should be returning values to)
!
!##############################################################################

      subroutine MR_Read_Met_DimVars(iy)

      integer, optional,intent(in) :: iy  ! Note: this is only needed for MR_iwind=5
                                          !       since we need to know how many
                                          !       metsteps to allocate

      integer            :: i
      logical            :: IsThere      = .false.
      character(len=130) :: tmp_str      = ""
      character(len=130) :: iw5filename  = ""
      integer            :: ivar
      real(kind=8)       :: inhour
      character(len=8)   :: date
      character(len=10)  :: time2
      character(len=5)   :: zone
      integer            :: values(8)

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine MR_Set_iwind5_filenames(inhour,ivar,infile)
          real(kind=8)      ,intent(in)  :: inhour
          integer           ,intent(in)  :: ivar
          character(len=130),intent(out) :: infile
        end subroutine MR_Set_iwind5_filenames
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_Met_DimVars                          ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .false., &  ! CALLED_MR_Read_Met_DimVars
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times

#ifdef USEPOINTERS
      if(.not.associated(MR_windfiles))then
#else
      if(.not.allocated(MR_windfiles))then
#endif
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  The list of windfile names, MR_windfiles, has not been"
          write(errlog(io),*)"           allocated.  The calling program must allocate this"
          write(errlog(io),*)"           array via the subroutine MR_Allocate_FullMetFileList."
        endif;enddo
        stop 1
      endif

      if(MR_iwind.eq.5)then
        ! For iwind=5 files (NCEP 2.5 degree reanalysis, NOAA, etc. ), only the directory
        ! was read into slot 1 of MR_windfiles(:).  We need to copy to all other slots
        ! to make sure we don't throw an error
        do i=1,MR_iwindfiles
          MR_windfiles(i)   = MR_windfiles(1)
        enddo
        if(present(iy)) then
          ! This is needed at this point for allocating the number
          ! of steps per file (this depends on the year), but this
          ! is reset in MR_Set_Met_Times when the actual start time
          ! is given
          MR_Comp_StartYear = iy
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"MR WARNING: If the MR_iwind=5 (NCEP Reannalysis, NOAA Reannalysis, "
            write(outlog(io),*)"            ERA5, etc.) are used, then MR_Read_Met_DimVars"
            write(outlog(io),*)"            should be called with a start year.  This is needed"
            write(outlog(io),*)"            to allocate the correct number of time steps per file."
            write(outlog(io),*)"            Setting MR_Comp_StartYear to 2018 for a non-leap year."
            write(outlog(io),*)"            However, MR_Comp_StartYear will be checked later to"
            write(outlog(io),*)"            verify that the start year is not a leap year."
            write(outlog(io),*)"            If there is an inconsistancy, the program will stop."
            write(outlog(io),*)"            If MR_Comp_StartYear is changed to a leap year outside"
            write(outlog(io),*)"            of MetReader, then the results will be incorrect."
          endif;enddo
          ! Setting to year of program execution
          call date_and_time(date,time2,zone,values)
          MR_Comp_StartYear = values(1)
        endif

      endif

      ! Verify that the first windfile has been changed from its initialized value
      tmp_str = MR_windfiles(1)
      if(tmp_str(1:15).eq.'              ')then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: The array MR_windfiles appears to not be set."
          write(errlog(io),*)"          Please have the calling program write the names of"
          write(errlog(io),*)"          the windfiles to MR_windfiles(1:MR_iwindfiles)"
          write(errlog(io),*)" "
          write(errlog(io),*)"          Contents of the first slot is -",MR_windfiles(1),"-"
        endif;enddo
        stop 1
      endif

      ! Check the existence of the wind files
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"  Verifying existance of windfiles:"
      endif;enddo
      ! Note MR_iwind=5 cases have the number of windfiles (MR_iwindfiles)
      ! modified in MR_Allocate_FullMetFileList to be the number of
      ! anticipated files based on the length of the simulation and the number
      ! of steps per file.
      if(MR_iwind.ne.5)then
        do i=1,MR_iwindfiles
          inquire( file=trim(adjustl(MR_windfiles(i))), exist=IsThere )
          if(.not.IsThere)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: Could not find windfile ",i
              write(errlog(io),*)"          ",trim(adjustl(MR_windfiles(i)))
            endif;enddo
            stop 1
          endif
        enddo
      else
        MR_iw5_root = MR_windfiles(1)
        do i = 1,MR_iwindfiles-1
          inhour = MR_Comp_StartHour + (i-1)*MR_iw5_hours_per_file
          do ivar = 1,5
#ifdef USENETCDF
            call MR_Set_iwind5_filenames(inhour,ivar,iw5filename)
#else
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: Currently, MR_iwind=5 required netcdf."
              write(errlog(io),*)"          Recompile MetReader with netcdf enabled"
            endif;enddo
            stop 1
#endif
            inquire( file=trim(adjustl(iw5filename)), exist=IsThere )
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" ",i,trim(adjustl(iw5filename)),IsThere
            endif;enddo
            if(.not.IsThere)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then           
                write(errlog(io),*)"MR ERROR: Could not find windfile ",i
                write(errlog(io),*)"          ",trim(adjustl(iw5filename))
              endif;enddo
              stop 1
            endif
          enddo
        enddo
      endif

      ! Now set up the full spatial and temporal grids
      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_Met_DimVars_ASCII_1d
      case(2)
        !call Set_Read_Met_DimVars_ASCII_3d
      case (3:5)
         ! call routine to populate variable lists and dimension axies based on
         ! MR_iwindformat; also do a bit of error checking by reading netcdf files
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_Met_Times_netcdf
          call MR_Read_Met_DimVars_netcdf
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_Met_Times_GRIB
          call MR_Read_Met_DimVars_GRIB
#endif
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then   
            write(errlog(io),*)"MR ERROR: Unknown MR_idataFormat:",MR_idataFormat
          endif;enddo
          stop 1
        endif
      case default
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then    
          write(errlog(io),*)"MR ERROR:  MR_iwind not 1-4."
          write(errlog(io),*)"           Exiting in MR_Read_Met_DimVars"
        endif;enddo
        stop 1
      end select

      CALLED_MR_Read_Met_DimVars = .true.

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_Met_DimVars

!##############################################################################


!##############################################################################
!
!     MR_Set_CompProjection
!
!     This subroutine only sets the projection parameters for the computational
!     grid for MetReader.  The library needs to know these parameters in order
!     to know how to map the met grid onto the comp grid and whether or not to
!     rotate wind vectors.
!
!     This subroutine can be called from the calling program at any time, but
!     must be called before MR_Initialize_Met_Grids
!
!##############################################################################

      subroutine MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)

      use projection,    only : &
         PJ_proj_for,PJ_proj_inv

      logical     ,intent(in) :: LL_flag
      integer     ,intent(in) :: ipf
      real(kind=8),intent(in) :: lam0,phi0
      real(kind=8),intent(in) :: phi1
      real(kind=8),intent(in) :: phi2
      real(kind=8),intent(in) :: ko
      real(kind=8),intent(in) :: Re
!      real(kind=8) :: inx1,outx,inx2
!      real(kind=8) :: iny1,outy,iny2

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Set_CompProjection                        ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      isLatLon_CompGrid = LL_flag

      Comp_lam0 = 0.0_8
      Comp_phi0 = 0.0_8
      Comp_phi1 = 0.0_8
      Comp_phi2 = 0.0_8
      Comp_k0   = 0.0_8
      Comp_Re   = 0.0_8

      if(.not.isLatLon_CompGrid)then
        Comp_iprojflag = ipf 
        if(Comp_iprojflag.eq.1)then
          ! Polar stereographic
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_k0      = ko
          Comp_Re      = Re
        elseif(Comp_iprojflag.eq.2)then
          ! Albers Equal Area
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
        elseif(Comp_iprojflag.eq.3)then
          ! UTM
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : WARNING: UTM not yet verified"
          endif;enddo
          stop 1
        elseif(Comp_iprojflag.eq.4)then
          ! Lambert conformal conic (NARR, NAM218, NAM221)
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
          Comp_Re      = Re
        elseif(Comp_iprojflag.eq.5)then
          ! Mercator (NAM196)
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_Re      = Re
        endif
      endif

      !HFS to do: perform a sanity check on these projection parameters
      ! take a lon/lat, project, then inverse project and check
      !if(.not.isLatLon_CompGrid)then
      !  write(*,*)x_comp_sp(1),y_comp_sp(1)
      !  inx1 = x_comp_sp(1)
      !  iny1 = y_comp_sp(1)
      !  call PJ_proj_inv(inx1,iny1, &
      !                 Comp_iprojflag,&
      !                 Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re,&
      !                 outx,outy)

      !  call PJ_proj_for(outx,outy, &
      !                 Comp_iprojflag,&
      !                 Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re,&
      !                 inx2,iny2)
      !  write(*,*)inx1,outx,inx2
      !  write(*,*)iny1,outy,iny2
      !  stop 99
      !endif

      CALLED_MR_Set_CompProjection = .true.

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Set_CompProjection

!##############################################################################
!
!     MR_Initialize_Met_Grids
!
!     This subroutine sets up local copies of the computational grid and
!     evaluates the full NWP grid for the subgrid needed, through MR_Set_MetComp_Grids
!     for NWP cases and MR_Set_MetComp_Grids_1dascii for 1-d or sonde cases.
!
!     From the calling program, this must be called after all of:
!       MR_Allocate_FullMetFileList
!       MR_Read_Met_DimVars
!       MR_Set_CompProjection
!
!     In some cases, it is useful to have velocity values saved within MetReader.
!     For example, velocities might be read to calculate diffusivities, stored
!     locally, then reused later for advection.  If the calling program sets the
!     variable MR_Save_Velocities=.true. , then space is allocated for the local
!     copies of velocity.
!
!     Takes as input :: specs of computational grid defined in the calling program
!     Sets: n[x,y,z]_comp and [x,y,z]_comp_sp :: local variables holding computational grid info
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           rdphi_MetP_sp      :: length scale along y (in meters)
!           rdlambda_MetP_sp   :: length scale along x (in meters)
!
!##############################################################################

      subroutine MR_Initialize_Met_Grids(nx,ny,nz, &
                                      dumx_sp,dumy_sp,dumz_sp,periodic)

      integer      ,intent(in) :: nx,ny,nz
      real(kind=sp),intent(in) :: dumx_sp(nx)
      real(kind=sp),intent(in) :: dumy_sp(ny)
      real(kind=sp),intent(in) :: dumz_sp(nz)
      logical      ,intent(in) :: periodic

      integer :: i,j,k
      integer :: io                           ! Index for output streams

      !                allocates *_Met_P for subset of Met grid on pressure levels
      !                          *_Met_H for subset of Met grid on height levels
      !                          *_comp_H for data regridded to computational grid

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Initialize_Met_Grids                      ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times
      nx_comp = nx
      ny_comp = ny
      nz_comp = nz
      allocate(x_comp_sp(nx_comp))
      allocate(y_comp_sp(ny_comp))
      allocate(z_comp_sp(nz_comp))
      allocate(s_comp_sp(nz_comp))
      if(MR_useCompGrid.eqv..false.)then
        ! This is the case where we will not be interpolating to a computational grid, but want access
        ! to the full met grid.  All parameters to this subroutine should have been dummy values
        x_comp_sp(1)  = x_fullmet_sp(1)
        x_comp_sp(nx) = x_fullmet_sp(nx_fullmet)
        y_comp_sp(1)  = min(y_fullmet_sp(1),y_fullmet_sp(ny_fullmet))
        y_comp_sp(ny) = max(y_fullmet_sp(1),y_fullmet_sp(ny_fullmet))
        z_comp_sp(1)  = 0.0_sp
        z_comp_sp(nz) = 1.1_sp*MR_Max_geoH_metP_predicted
        s_comp_sp = z_comp_sp   ! Initialize s to z.  This might be changed in MR_Set_SigmaAlt_Scaling
      else
        x_comp_sp = dumx_sp
        y_comp_sp = dumy_sp
        z_comp_sp = dumz_sp
        s_comp_sp = z_comp_sp
      endif
      ! Do some error-checking on these grids to make sure they are
      ! strictly increasing
      do i=1,nx-1
        if(x_comp_sp(i).gt.x_comp_sp(i+1))then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  x_comp not strictly increasing"
            write(errlog(io),*)x_comp_sp
          endif;enddo
          stop 1
        endif
      enddo

      do j=1,ny-1
        if(y_comp_sp(j).gt.y_comp_sp(j+1))then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  y_comp not strictly increasing"
            write(errlog(io),*)y_comp_sp
          endif;enddo
          stop 1
        endif
      enddo

      dx_comp = x_comp_sp(2) - x_comp_sp(1)
      dy_comp = abs(y_comp_sp(2) - y_comp_sp(1))
      MaxZ_comp_sp = maxval(z_comp_sp)

      IsPeriodic_CompGrid = periodic

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Set_MetComp_Grids_ASCII_1d
        call MR_Set_Comp2Met_Map
      case(2)
        !call MR_Set_MetComp_Grids_ASCII_3d
      case (3:5)
        ! Now that we have the full grids defined in MR_Read_Met_DimVars,
        ! calculate the subgrid needed for the simulation
        call MR_Set_MetComp_Grids

        if(MR_Save_Velocities)then
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Velocities will be saved on the metP grid"
          endif;enddo
          allocate(MR_vx_metP_last(nx_submet,ny_submet,np_fullmet));MR_vx_metP_last(:,:,:)=0.0_sp
          allocate(MR_vx_metP_next(nx_submet,ny_submet,np_fullmet));MR_vx_metP_next(:,:,:)=0.0_sp
          allocate(MR_vy_metP_last(nx_submet,ny_submet,np_fullmet));MR_vy_metP_last(:,:,:)=0.0_sp
          allocate(MR_vy_metP_next(nx_submet,ny_submet,np_fullmet));MR_vy_metP_next(:,:,:)=0.0_sp
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Velocities not saved"
          endif;enddo
        endif

      case default
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then       
          write(errlog(io),*)"MR ERROR:  MR_iwind not 1:5."
          write(errlog(io),*)"           Exiting in MR_Initialize_Met_Grids"
        endif;enddo
        stop 1
      end select

      ! For LatLon Met grids, calculate some additional geometry terms
      !  These are currently only needed for calculating DelMetP_Dx and
      !  DelMetP_Dy
      if(IsLatLon_MetGrid)then
#ifdef USEPOINTERS
        if(associated(MR_dx_met))then
#else
        if(allocated(MR_dx_met))then
#endif
          allocate(rdphi_MetP_sp(ny_submet,np_fullmet))             ;   rdphi_MetP_sp(:,:)=0.0_sp
          allocate(rdlambda_MetP_sp(nx_submet,ny_submet,np_fullmet));rdlambda_MetP_sp(:,:,:)=0.0_sp

          do k=1,np_fullmet
            if(IsRegular_MetGrid)then
              ! length scale along y (in meters)
              rdphi_MetP_sp(:,k) = dy_met_const*MR_DEG2RAD * (MR_RAD_EARTH+z_approx(k))*1000.0_sp
              do j=1,ny_submet
                ! length scale along x (in meters)
                rdlambda_MetP_sp(:,j,k) =(MR_RAD_EARTH+z_approx(k))*1000.0_sp * &
                                        cos(MR_DEG2RAD*(y_submet_sp(j)-0.5_sp*dy_met_const)) * &
                                        dx_met_const*MR_DEG2RAD
              enddo
            else
              do i=1,nx_submet
                ! length scale along y (in meters)
                rdphi_MetP_sp(:,k) = MR_dy_submet(:)*MR_DEG2RAD * (MR_RAD_EARTH+z_approx(k))*1000.0_sp
                do j=1,ny_submet
                  ! length scale along x (in meters)
                  rdlambda_MetP_sp(i,j,k) =(MR_RAD_EARTH+z_approx(k))*1000.0_sp * &
                                          cos(MR_DEG2RAD*(y_submet_sp(j)-0.5_sp*MR_dy_submet(j))) * &
                                          MR_dx_submet(i)*MR_DEG2RAD
                enddo
              enddo
            endif
          enddo
          MR_minlen = maxval(rdlambda_MetP_sp(:,:,:)) / 1000.0_sp ! in km
          do i=1,nx_submet
            do j=1,ny_submet
              if(MR_minlen.gt.rdlambda_MetP_sp(i,j,1)) MR_minlen=rdlambda_MetP_sp(i,j,1)
              if(MR_minlen.gt.rdphi_MetP_sp(j,1))      MR_minlen=rdphi_MetP_sp(j,1)
              MR_sigma_nz_submet(i,j) = rdphi_MetP_sp(j,1)*rdlambda_MetP_sp(i,j,1)/1.0-6_sp ! in km2
            enddo
          enddo
        endif
      else
#ifdef USEPOINTERS
        if(associated(MR_dx_met))then
#else
        if(allocated(MR_dx_met))then
#endif
          ! This is the branch for projected NWP files
          MR_minlen = min(minval(MR_dx_met),minval(MR_dy_met)) ! (in km)
          do i=1,nx_submet
            do j=1,ny_submet
              MR_sigma_nz_submet(i,j) = MR_dx_submet(i)*MR_dy_submet(j) ! in km2
            enddo
          enddo
        else
          ! MR_dx_met and MR_dy_met might not be defined for radio sonde or ASCII 
          ! grids.  Just set MR_minlen to 10% of min domain dimension
          MR_minlen = 0.1_sp * (x_comp_sp(nx) - x_comp_sp(1))
          MR_minlen = min(MR_minlen,0.1_sp * (y_comp_sp(ny) - y_comp_sp(1))) ! in km2
          !MR_sigma_nz_submet(i,j) = 
        endif
      endif

      if(MR_useTopo)then
        if(MR_iwind.eq.1)then
          ! There is no grid for ASCII input, so just use comp grid
          allocate(MR_Topo_met(-1:nx_comp+2,-1:ny_comp+2));  MR_Topo_met(:,:)   = 0.0_sp
          allocate(MR_jacob_met(-1:nx_comp+2,-1:ny_comp+2)); MR_jacob_met(:,:)  = 1.0_sp
        else
          allocate(MR_Topo_met(nx_submet,ny_submet));  MR_Topo_met(:,:)   = 0.0_sp
          allocate(MR_jacob_met(nx_submet,ny_submet)); MR_jacob_met(:,:)  = 1.0_sp
        endif
        allocate(MR_Topo_comp(-1:nx_comp+2,-1:ny_comp+2));     MR_Topo_comp(:,:)  = 0.0_sp
        allocate(MR_jacob_comp(-1:nx_comp+2,-1:ny_comp+2));    MR_jacob_comp(:,:) = 1.0_sp
      endif

      CALLED_MR_Initialize_Met_Grids = .true.

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Initialize_Met_Grids

!##############################################################################
!
!     MR_Set_SigmaAlt_Scaling
!
!     This subroutine sets local copies of the necessary variables for vertical
!     grid modifications when including topography. No action is taken unless
!     the variable MR_useTopo is set to true. MR_ZScaling_ID must also be set:
!       MR_ZScaling = 0 :: no topography, lower boundary is s=z=0 everywhere
!       MR_ZScaling = 1 :: grid is altitude shifted; s = z-zsurf
!       MR_ZScaling = 2 :: grid is altitude scaled;  s = (z-zsurf)/(ztop-zsurf)
!     Takes as input :: vertical grid specification
!     Sets: s_comp_sp
!           MR_jacob_comp
!           MR_jacob_met
!
!##############################################################################

      subroutine MR_Set_SigmaAlt_Scaling(nz,dums_sp)

      integer      ,intent(in) :: nz
      real(kind=sp),intent(in) :: dums_sp(nz)

      integer :: io                           ! Index for output streams

      if(.not.MR_useTopo)then
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)"MR WARNING: Trying to set sigma-alt. coordinates, but MR_useTopo=F"
        endif;enddo
        return
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Set_SigmaAlt_Scaling                      ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      s_comp_sp(1:nz) = dums_sp(1:nz)
      if(MR_ZScaling_ID.eq.0)then
        ! no topo
        MR_jacob_comp(-1:nx_comp+2,-1:ny_comp+2) = 1.0_sp
        if(MR_iwind.eq.1)then
          MR_jacob_met = MR_jacob_comp
        else
          MR_jacob_met(1:nx_submet,1:ny_submet) = 1.0_sp
        endif
      elseif(MR_ZScaling_ID.eq.1)then
        ! shifted-altitude (s=z-zsurf)
        MR_jacob_comp(-1:nx_comp+2,-1:ny_comp+2) = 1.0_sp
        MR_jacob_met(1:nx_submet,1:ny_submet) = 1.0_sp
      elseif(MR_ZScaling_ID.eq.2)then
        ! sigma-altitude (s=Ztop*(z-zsurf)/(Ztop-zsurf))
        !  Note: this is not the same as in Jacobson Eq 5.89, but we want the orientation of the
        !        vertical coordinate to be the same (-> positive Jacobian). Effect of topography
        !        diminishes with altitude.

        ! Just set the s-values for now without applying topography
        MR_jacob_comp(-1:nx_comp+2,-1:ny_comp+2) = (MR_ztop - MR_Topo_comp(-1:nx_comp+2,-1:ny_comp+2))/MR_ztop
        if(MR_iwind.eq.1)then
          MR_jacob_met = MR_jacob_comp
        else
          MR_jacob_met(1:nx_submet,1:ny_submet) = (MR_ztop - MR_Topo_met(1:nx_submet,1:ny_submet))/MR_ztop
        endif
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)"MR WARNING: Topography scheme not recognized."
          write(outlog(io),*)"            Reverting to altitude."
        endif;enddo
        MR_jacob_comp(-1:nx_comp+2,-1:ny_comp+2) = 1.0_sp
        if(MR_iwind.eq.1)then
          MR_jacob_met = MR_jacob_comp
        else
          MR_jacob_met(1:nx_submet,1:ny_submet) = 1.0_sp
        endif
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Set_SigmaAlt_Scaling

!##############################################################################
!
!     MR_Set_Met_Times
!
!     This subroutine opens each file of the MR_windfile list, determines the hour
!     of each time step and logs the files/time-steps needed to bracket the
!     simulation duration.
!
!     Sets: MetStep_File(i=1:MR_MetSteps_Total)   :: contains file name for met step i
!           MetStep_findex(i=1:MR_MetSteps_Total) :: contains the index of the file in
!                                            the MR_windfiles(:) list
!           MetStep_tindex(i=1:MR_MetSteps_Total) :: contains the index within the file
!                                            for met step i
!
!##############################################################################

      subroutine MR_Set_Met_Times(eStartHour,Duration)

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer, parameter :: NT_MAXOUT = 100

      real(kind=8),intent(in) :: eStartHour
      real(kind=8),intent(in) :: Duration

      integer :: i
      integer :: iw
      integer            :: iwstep
      integer            :: istep
      real(kind=8)       :: stephour
      logical :: Found_First_Step = .false.
      logical :: Found_Last_Step  = .false.
      integer :: nMetSteps_Comp   = 0
      real(kind=8) :: StepInterval
      real(kind=8) :: met_t1,met_t2,met_dt1
      logical      :: prestep, poststep
      integer :: files_per_tstep

      integer :: io                           ! Index for output streams

      INTERFACE
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_HourOfDay
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_MonthOfEvent

        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_DayOfEvent
        integer function HS_DayOfYear(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_DayOfYear
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8) ,intent(in):: HoursSince
          integer      ,intent(in):: byear
          logical      ,intent(in):: useLeaps
        end function HS_yyyymmddhhmm_since

      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Set_Met_Times                             ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times

      if(MR_useCompTime.eqv..false.)then
        ! This is for the case where only one file is provided and steps will be accessed
        ! through the file step index.  e.g. used for evaluating windfiles
        ! Here, we assign the start hour to the first step and the duration with the last step
        MR_Comp_StartHour     = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
        MR_Comp_Time_in_hours = MR_windfile_stephour(1,nt_fullmet)
      else
        MR_Comp_StartHour        = eStartHour
        MR_Comp_Time_in_hours    = Duration
      endif
      MR_Comp_StartYear        = HS_YearOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
      MR_Comp_StartMonth       = HS_MonthOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
      MR_Comp_StartDay         = HS_DayOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)

      ! Now we want to loop through all the steps of each windfile and find the
      ! file/step that immediately preceeds the time needed (MR_Comp_StartHour)
      !  First, define MR_Comp_StartHour if we are running a forecast run

      ! Now is a good time to check to make sure the MR_Comp_StartHour is within
      !  the range of data
      ! First we need the step interval for the first step
      ! This is for checking pre- and post-steps
      if(MR_iwind.eq.1)then
        if(MR_iwindformat.eq.1)then
          ! 1-D ASCII
          if(MR_Snd_nt_fullmet.eq.1)then
            StepInterval = 1000.0_8 ! some large number
          else
            met_t1 = MR_windfile_starthour(1)
            met_t2 = MR_windfile_starthour(MR_nSnd_Locs+1)
            StepInterval = met_t2 - met_t1
          endif
        elseif(MR_iwindformat.eq.2)then
          ! 1-D Radiosonde
          StepInterval = 12.0_8
        endif
      else
        ! For all other cases, just check the arrays read in MR_Read_Met_DimVars
        met_t1  = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
        if(nt_fullmet.gt.1)then
            ! If there are multiple steps per file, use the second step
          met_t2 = MR_windfile_starthour(1)+MR_windfile_stephour(1,2)
        elseif(MR_iwindfiles.gt.1)then
            ! If only 1 step/file, use the next file (if present)
          met_t2 = MR_windfile_starthour(1)+MR_windfile_stephour(2,1)
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: Only one time step available."
            write(outlog(io),*)"         Setting interval step to 1000.0"
          endif;enddo
          met_t2 = met_t1 + 1000.0
        endif
        StepInterval = met_t2 - met_t1
      endif

      !   Checking if a pre-step is needed
      if(MR_iwind.eq.1.and.MR_iwindformat.eq.1)then
        ! For the ASCII profile cases (not the radiosonde), hours are given as offset
        ! from the start time.  So reset all hours to relative to eStartHour
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Note:  The hours value in 1d ASCII profiles are interpreted as offset"
          write(outlog(io),*)"       hours.  Shifting the file time to reflect the requested start"
          write(outlog(io),*)"       time plus offset."
        endif;enddo
        do i=1,MR_Snd_nt_fullmet
          MR_windfile_starthour(i) = MR_windfile_starthour(i) + MR_Comp_StartHour
        enddo
      endif
      met_t1  = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
      met_dt1 = StepInterval
      if(MR_Comp_StartHour.lt.met_t1)then
        ! Start time requested is before that available, check if we can extrapolate
        if(MR_Comp_StartHour.ge.met_t1-met_dt1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: Start time is before the first time step"
            write(outlog(io),*)"         However, it is within one interval so we will"
            write(outlog(io),*)"         apply the value at MR_windfile_starthour(1)."
            write(outlog(io),*)"         Using a time step interval of ",StepInterval
            if(StepInterval.gt.720.0_dp)then
              if(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)then
                write(outlog(io),*)"         Note: This interval is set high since only one"
                write(outlog(io),*)"               sonde file was provided."
              endif
            endif
          endif;enddo
          prestep = .true.
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Start time is before the first available data and"
            write(errlog(io),*)"       cannot be extrapolated."
          endif;enddo
          stop 1
        endif
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Start hour of simulation is in range of NWP data"
        endif;enddo
        prestep  = .false.
      endif
      !   Checking if a post-step is needed
      met_t1  = MR_windfile_starthour(MR_iwindfiles)+MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
      met_dt1 = StepInterval
      if(MR_Comp_StartHour+MR_Comp_Time_in_hours.ge.met_t1)then
        ! Start time requested is after that available, check if we can extrapolate
        if(MR_Comp_StartHour+MR_Comp_Time_in_hours.le.met_t1+met_dt1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: End time is at or after the last time step."
            write(outlog(io),*)"         However, it is within one interval so we will"
            write(outlog(io),*)"         apply the value at MR_windfile_starthour(MR_iwindfiles)"
            if(StepInterval.gt.720.0_dp)then
              if(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)then
                write(outlog(io),*)"         Note: This interval is set high since only one"
                write(outlog(io),*)"               sonde file was provided."
              endif
            endif
          endif;enddo
          poststep = .true.
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: End time is after the last available data and"
            write(errlog(io),*)"       cannot be extrapolated."
            write(errlog(io),*)"  MR_Comp_StartHour    = ",real(MR_Comp_StartHour,kind=4),&
                                    HS_yyyymmddhhmm_since(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
            write(errlog(io),*)"  MR_Comp_Time_in_hours= ",real(MR_Comp_Time_in_hours,kind=4)
            write(errlog(io),*)"  met_t1               = ",real(met_t1,kind=4),&
                                    HS_yyyymmddhhmm_since(met_t1,MR_BaseYear,MR_useLeap)
            write(errlog(io),*)"  met_dt1              = ",real(met_dt1,kind=4)
          endif;enddo
          stop 1
        endif
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"End hour of simulation is in range of NWP data"
        endif;enddo
        poststep  = .false.
      endif

      ! Loop through all the files and steps and count how many are needed to bracket the time
      ! needed (MR_Comp_StartHour -> MR_Comp_StartHour+MR_Comp_Time_in_hours).
      ! Note: If pre-step or post-step is needed, istep will be incremented accordingly.
      ! Once we know the number of steps needed, we will allocate space, then fill the variables
      ! with just the step info needed.
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),'(a80)')&
         '    File num  | Step in file  |  stephour   |   SimStartHour   | nMetStep | Note'
      endif;enddo
      if(prestep)then
        Found_First_Step = .true.
        istep = 1
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - StepInterval
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),150)1,1,stephour,MR_Comp_StartHour,&
                          istep,"Prestep before MR_Comp_StartHour "
        endif;enddo
      else
        Found_First_Step = .false.
        istep = 0
      endif
      Found_Last_Step  = .false.
      if(MR_iwind.eq.1.and.MR_iwindformat.eq.2)then
        ! Radiosonde files might have more than one file per time step
        files_per_tstep = MR_nSnd_Locs
      else
        files_per_tstep = 1
      endif
      do iw = 1,MR_iwindfiles,files_per_tstep
        do iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          if(stephour.lt.MR_Comp_StartHour)then
            Found_First_Step = .false.
            istep = 1
            if(MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then
              ! This suppresses outputting the windfile step info for the files that contain a year's
              ! worth of data (NCEP 2.5, etc)
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),150)                    &
                  iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                  "Before or at MR_Comp_StartHour."
              endif;enddo
            endif
          else
            !  Otherwise, increment index if we are still in the needed time bracket
            if(.not.Found_Last_Step)then
              Found_First_Step = .true.
              istep = istep + 1
              if(MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then
                do io=1,MR_nio;if(VB(io).le.verbosity_info)then 
                  write(outlog(io),150)                    &
                    iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                    "After MR_Comp_StartHour "
                endif;enddo
              endif
            endif
          endif

          ! Check if the current step hour is exclusively after the end of the simulation
          if(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)then
            if(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          endif
          if(Found_Last_Step &
             .and.MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then ! This condition suppresses output for NCEP 2.5
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then 
              write(outlog(io),150)                    &
                  iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                  "At or after END OF SIM  "
            endif;enddo
          endif
          if(Found_Last_Step)exit
        enddo
      enddo

      ! If we went through all the steps and didn't find the first step, then issue an error message
      if(.not.Found_First_Step)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  Something is wrong.  Could not find the first MetStep needed for"
          write(errlog(io),*)"           the simulation."
        endif;enddo
        stop 1
      endif
      ! If we went through all the steps and didn't find the last step, and if poststep=T, then
      ! increment istep and set nMetSteps_Comp
      if(.not.Found_Last_Step)then
        if(poststep)then
          nMetSteps_Comp = istep+1
          stephour = MR_windfile_starthour(MR_iwindfiles) + &
                       MR_windfile_stephour(MR_iwindfiles,MR_windfiles_nt_fullmet(MR_iwindfiles)) + &
                       StepInterval
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then  
            write(outlog(io),150)MR_iwindfiles,1,stephour,&
                       MR_Comp_StartHour,nMetSteps_Comp,"Poststep after END OF SIM  "
          endif;enddo
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then  
            write(errlog(io),*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
            write(errlog(io),*)"           the simulation."
          endif;enddo
          stop 1
        endif
      endif
 150  format(8x,i3,9x,i4,4x,f15.2,2x,f13.2,10x,i4,8x,a31)

      ! We now have the number of steps needed for the computation
      ! Allocate the lists
      MR_MetSteps_Total = nMetSteps_Comp
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then  
        write(outlog(io),*)"MR: Allocating space for ",MR_MetSteps_Total,"step(s)"
      endif;enddo
      if(prestep.or.poststep)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then  
          write(outlog(io),*)"         Including:"
          if(prestep) write(outlog(io),*)"             1 prestep"
          if(poststep)write(outlog(io),*)"             1 poststep"
        endif;enddo
      endif
#ifdef USEPOINTERS
      if (.not. associated(MR_MetStep_File)) then
#else
      if (.not. allocated(MR_MetStep_File)) then
#endif
        allocate(MR_MetStep_File(MR_MetSteps_Total))               ;MR_MetStep_File(:)=''
        allocate(MR_MetStep_findex(MR_MetSteps_Total))             ;MR_MetStep_findex(:)=0
        allocate(MR_MetStep_tindex(MR_MetSteps_Total))             ;MR_MetStep_tindex(:)=0
        allocate(MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total));MR_MetStep_Hour_since_baseyear(:)  =0.0_sp  
        allocate(MR_MetStep_Interval(MR_MetSteps_Total))           ;MR_MetStep_Interval(:)=0  
        allocate(MR_MetStep_year(MR_MetSteps_Total))               ;MR_MetStep_year(:)=0  
        allocate(MR_MetStep_month(MR_MetSteps_Total))              ;MR_MetStep_month(:)=  0  
        allocate(MR_MetStep_day(MR_MetSteps_Total))                ;MR_MetStep_day(:)=0  
        allocate(MR_MetStep_DOY(MR_MetSteps_Total))                ;MR_MetStep_DOY(:)=0  
        allocate(MR_MetStep_Hour_Of_Day(MR_MetSteps_Total))        ;MR_MetStep_Hour_Of_Day(:)=0.0_sp
        allocate(MR_iwind5_year(MR_MetSteps_Total))                ;MR_iwind5_year(:)=0
      endif

      ! Finally, we need to loop through the steps exactly as above, but this time populate 
      ! the lists just allocated
      if(prestep)then
        Found_First_Step = .true.
        istep = 1                                                         ! pre-step is one interval back
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - StepInterval
        MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(1)))
        MR_MetStep_findex(istep)          = 1
        MR_MetStep_tindex(istep)          = 1
        MR_MetStep_Hour_since_baseyear(istep) = stephour
        MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
        MR_iwind5_year(istep)             = MR_Comp_StartYear
      else
        Found_First_Step = .false.
        istep = 0
      endif
      Found_Last_Step  = .false.
      do iw = 1,MR_iwindfiles,files_per_tstep
        do iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          if(stephour.lt.MR_Comp_StartHour)then
            Found_First_Step = .false.
            istep = 1
          else
            !  Otherwise, increment index if we are still in the needed time bracket
            if(.not.Found_Last_Step)then
              Found_First_Step = .true.
              istep = istep + 1
            endif
          endif

          if(.not.Found_Last_Step)then
            MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(iw)))
            MR_MetStep_findex(istep)          = iw
            MR_MetStep_tindex(istep)          = iwstep
            MR_MetStep_Hour_since_baseyear(istep) = stephour
            MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
            MR_iwind5_year(istep)             = MR_MetStep_year(istep)
          endif

          ! Check if the current step hour is exclusively after the end of the simulation
          if(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)then
            if(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          endif
        enddo
      enddo
      if(.not.Found_Last_Step)then
        if(poststep)then
          istep = istep+1                                      ! post-step is one interval later
          stephour = MR_MetStep_Hour_since_baseyear(istep-1) + StepInterval
          MR_MetStep_File(istep)            = MR_MetStep_File(istep-1)
          MR_MetStep_findex(istep)          = MR_MetStep_findex(istep-1)
          MR_MetStep_tindex(istep)          = MR_MetStep_tindex(istep-1)
          MR_MetStep_Hour_since_baseyear(istep) = stephour
          MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
          MR_iwind5_year(istep)             = MR_MetStep_year(istep)
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then  
            write(errlog(io),*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
            write(errlog(io),*)"           the simulation."
          endif;enddo
          stop 1
        endif
      endif

      do istep = 1,MR_MetSteps_Total-1
        MR_MetStep_Interval(istep) = MR_MetStep_Hour_since_baseyear(istep+1) - &
                                     MR_MetStep_Hour_since_baseyear(istep)
      enddo
      MR_MetStep_Interval(MR_MetSteps_Total)=MR_MetStep_Interval(MR_MetSteps_Total-1)

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then  
        write(outlog(io),'(a12,f13.2)')"Comp start = ",MR_Comp_StartHour
        write(outlog(io),'(a12,f13.2)')"Comp end   = ",MR_Comp_StartHour+MR_Comp_Time_in_hours
      endif;enddo

      do istep = 1,MR_MetSteps_Total
        iw       = MR_MetStep_findex(istep)
        iwstep   = MR_MetStep_tindex(istep)
        stephour = MR_MetStep_Hour_since_baseyear(istep)
      enddo

      ! Make sure the wind model time window covers the entire simulation time
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then  
        write(outlog(io),99)
      endif;enddo
99    format(/,4x,'Making sure the mesoscale model time covers the simulation time . . . ')
      if (MR_MetStep_Hour_since_baseyear(1).gt.MR_Comp_StartHour) then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then  
          write(errlog(io),*)"MR ERROR:  Wind data starts after SimStartHour"
          write(errlog(io),*)"WindHour(MR_iwindfiles) = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
          write(errlog(io),*)"SimStartHour            = ",MR_Comp_StartHour
          write(errlog(io),*)"Simtime_in_hours        = ",MR_Comp_Time_in_hours
          write(errlog(io),*)"  All steps:"
          do i = 1,MR_MetSteps_Total
            write(errlog(io),*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
          enddo
        endif;enddo
        stop 1
      elseif (MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total).lt.(MR_Comp_StartHour+MR_Comp_Time_in_hours)) then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then     
          write(errlog(io),*)"MR ERROR:  Last met time is earlier than simulation time"
          write(errlog(io),*)"MR_MetSteps_Total     = ",MR_MetSteps_Total
          write(errlog(io),*)"Last time step       = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
          write(errlog(io),*)"SimStartHour         = ",MR_Comp_StartHour
          write(errlog(io),*)"Simtime_in_hours     = ",MR_Comp_Time_in_hours
          do i = 1,MR_MetSteps_Total
            write(errlog(io),*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
          enddo
          write(errlog(io) ,102)
        endif;enddo
102     format(4x,'Error.  The model ends before the end of the last eruption.  Program stopped.')
        stop 1
      endif
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),103)
      endif;enddo
103   format(4x,'Good.  It does.',/)

      CALLED_MR_Set_Met_Times = .true.

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Set_Met_Times


!##############################################################################
!
!     MR_Read_HGT_arrays
!
!     This subroutine basically does the same thing as MR_Read_3d_MetP_Variable
!     but specifically for the variable HGT (ivar=1).  This is needed because the
!     more general MR_Read_3d_MetP_Variable just returns data in dum3d_metP and
!     expects the calling program to save the results.  The HGT variables,
!     geoH_metP_last and geoH_metP_next need to persist locally, however,  since they are
!     used in QC calculations and are needed in converting PresVertVel (from Pa
!     s to m/s).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: istep :: specified the met step
!
!     Sets: MR_geoH_metP_last, MR_geoH_metP_next
!
!##############################################################################

      subroutine MR_Read_HGT_arrays(istep,reset_first_time)

      integer,intent(in)           :: istep
      logical, optional,intent(in) :: reset_first_time

      integer      :: ivar
      logical,save :: first_time = .true.
      integer      :: io                           ! Index for output streams

      if(present(reset_first_time)) then
        first_time = .true.
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_HGT_arrays                           ----------"
        write(outlog(io),*)istep,reset_first_time
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     Reading HGT array for istep = ",istep
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)
      ivar = 1 ! HGT
      if(first_time)then
        if(MR_idataFormat.eq.1)then
          call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
        elseif(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif
        MR_geoH_metP_last(:,:,:) = MR_dum3d_metP(:,:,:)
        first_time = .false.
        Max_geoH_metP_last = maxval(MR_geoH_metP_last(:,:,np_fullmet))
      else
        MR_geoH_metP_last(:,:,:) = MR_geoH_metP_next(:,:,:)
        Max_geoH_metP_last = Max_geoH_metP_next
        Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
      endif
      if(MR_idataFormat.eq.1)then
        call MR_Read_MetP_Variable_ASCII_1d(ivar,istep+1)
      elseif(MR_idataFormat.eq.2)then
#ifdef USENETCDF
        call MR_Read_MetP_Variable_netcdf(ivar,istep+1)
#endif
      elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
        call MR_Read_MetP_Variable_GRIB(ivar,istep+1)
#endif
      endif

      MR_geoH_metP_next(:,:,:) = MR_dum3d_metP(:,:,:)

      Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
        ! Now determine the maximum value provided between last and next steps
      Max_geoH_metP = max(Max_geoH_metP_last,Max_geoH_metP_next)
        ! and compare with the maximum value needed by the computational grid
      if(Max_geoH_metP.lt.z_comp_sp(nz_comp))then
        ! use highest needed point (rounded up to 5-km increment) for the
        ! padding height
        Suppl_H = real(ceiling(z_comp_sp(nz_comp)/5.0_sp),kind=sp) * 5.0_sp
      else
        ! otherwise use the height from the met files to determine padding
        Suppl_H = real(ceiling(Max_geoH_metP/5.0_sp),kind=sp) * 5.0_sp
      endif
      if(z_comp_sp(nz_comp).gt.Max_geoH_metP)then
        if(MR_iHeightHandler.eq.1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Computational grid extends higher than met grid"
            write(errlog(io),*)"           iHeightHandler = 1:  exiting"
          endif;enddo
          stop 1
        else ! MR_iHeightHandler.eq.2
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"MR WARNING : Computational grid extends higher than met grid"
            write(outlog(io),*)"             iHeightHandler = 2:  continuing values upward."
          endif;enddo
        endif
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_HGT_arrays


!##############################################################################
!
!     MR_Read_3d_MetP_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets: dum3d_metP
!
!##############################################################################

      subroutine MR_Read_3d_MetP_Variable(ivar,istep)

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_3d_MetP_Variable                     ----------"
        write(outlog(io),*)ivar,istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
      case(2)
        !call Read_3d_ASCII
      case (3:5)
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif

      case default

      end select

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_3d_MetP_Variable

!##############################################################################
!
!     MR_Read_3d_MetH_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar      :: specifies which variable to read
!                       istep     :: specified the met step
!     Reads : MR_dum3d_metP
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_Read_3d_MetH_Variable(ivar,istep)

      integer,intent(in)             :: ivar
      integer,intent(in)             :: istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext
      integer :: np_fully_padded
      real(kind=sp),dimension(:),allocatable :: dumVertCoord_sp
      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_3d_MetH_Variable                     ----------"
        write(outlog(io),*)ivar,istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      np_fully_padded = np_fullmet+1+np_fullmet_pad

      ! First get the variable on the native pressure coordinate
      call MR_Read_3d_MetP_Variable(ivar,istep)

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp
      ! Allocate 1-d arrays, padding one layer to the bottom and np_fullmet_pad
      ! layers on the top
      allocate(  z_col_metP(np_fully_padded));      z_col_metP(:) = 0.0_sp
      allocate(var_col_metP(np_fully_padded));    var_col_metP(:) = 0.0_sp
      allocate(var_col_metH(nz_comp))        ;    var_col_metH(:) = 0.0_sp
      allocate(dumVertCoord_sp(nz_comp))     ; dumVertCoord_sp(:) = 0.0_sp

      ! Create 1d arrays in P and regrid them into 1d arrays in z or s
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          if(istep.eq.MR_iMetStep_Now)then
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          else
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          endif
          ! Loop through the nodes and reset any that are negative or non-ascending in p
          ! This is usually never invoked, but would only affect potentially the bottom few nodes
          ! The intent is to maintain a monotonic GPH, while getting irrelevant nodes out of the
          ! way.
          do k=1,np_fullmet
            if(z_col_metP(k).le.real(k-1,kind=sp)*MR_MIN_DZ) z_col_metP(k)=real(k-1,kind=sp)*MR_MIN_DZ
          enddo
          var_col_metP(1)              = MR_dum3d_metP(i,j,1)
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while ((z_col_metP(knext)<z_col_metP(k-1)).and. &
                        (knext.lt.np_fullmet+1))
                 knext=knext+1
              enddo
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              enddo
            endif
          enddo
          if(np_fullmet_pad.gt.0)then
            ! Assign top node
            if(MR_iHeightHandler.eq.1)then
              z_col_metP(np_fullmet+2) = z_col_metP(np_fullmet+1) + 10.0_sp
            elseif(MR_iHeightHandler.eq.2)then
              z_col_metP(np_fullmet+2) = Suppl_H
              if(ivar.ne.5)then
                ! For all variables except temperature (ivar=5), just copy highest
                ! met node
                var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)
              else
                ! For temperature, use the 1976 US Standard Atmosphere
                var_col_metP(np_fullmet+2) = MR_Temp_US_StdAtm(Suppl_H)
              endif
            endif
          endif

          if(MR_ZScaling_ID.eq.0)then
            dumVertCoord_sp(1:nz_comp) = z_comp_sp(1:nz_comp)
          elseif(MR_ZScaling_ID.eq.1)then
            dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + s_comp_sp(1:nz_comp)
          elseif(MR_ZScaling_ID.eq.2)then
              ! this recovers the real-world z coordinate from the sigma level and topography
            if(MR_iwind.eq.1)then
              ! this is incorrect (should be MR_jacob_met(i,j)), but this isn't set for iwind=1
              dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + &
                                             s_comp_sp(1:nz_comp) * MR_jacob_comp(1,1)
            else
              dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + &
                                             s_comp_sp(1:nz_comp) * MR_jacob_met(i,j)
            endif
          else
            dumVertCoord_sp(1:nz_comp) = z_comp_sp(1:nz_comp)
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid
          do io=1,MR_nio;if(VB(io).le.verbosity_debug2)then
            write(outlog(io),*)
            write(outlog(io),*)"Callling MR_Regrid_P2H_linear from MR_Read_3d_MetH_Variable for ",i,j
          endif;enddo

          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,       var_col_metP, & 
                                    nz_comp,      dumVertCoord_sp,  var_col_metH)
          MR_dum3d_metH(i,j,:) = var_col_metH

        enddo ! j
      enddo  ! i

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Read_3d_MetH_Variable

!##############################################################################
!
!     MR_Read_3d_Met_Variable_to_CompH
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_3d_MetH_Variable                   (gets variable on Met_x, Met_y, Comp_z)
!        -> Read_3d_MetP_Variable               (gets variable on native Met subgrid)
!            -> Read_3d_MetP_Variable_[format]  (direct read of variable in
!                                                whatever format : NetCDF, GRIB1/2, ASCII)
!
!     Takes as input :: ivar      :: specifies which variable to read
!                       istep     :: specified the met step
!                       IsNext    :: T for saving velocity values
!                                    F default
!     Sets  : MR_dum3d_compH
!               MR_dum3d_metH and MR_dum3d_metP are also filled in the course of
!               generating MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Read_3d_Met_Variable_to_CompH(ivar,istep,IsNext)

      integer,intent(in)             :: ivar
      integer,intent(in)             :: istep
      logical,optional,intent(in)    :: IsNext

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp
      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_3d_Met_Variable_to_CompH             ----------"
        write(outlog(io),*)ivar,istep,IsNext
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

        ! First get the variable on the height coordinate (or s)
      call MR_Read_3d_MetH_Variable(ivar,istep)

      if(MR_Save_Velocities)then
        ! Note: this flag for saving the velocity values is useful in special
        ! cases such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then ! First check if this parameter was provided
          if(IsNext) then  ! second, check if it is true
            ! MR_dum3d_metP still contains the variable just read
            if(ivar.eq.2)then
              if(IsNext)then
                MR_vx_metP_last = MR_vx_metP_next
                MR_vx_metP_next = MR_dum3d_metP
              else
                MR_vx_metP_last = MR_dum3d_metP
              endif
            elseif(ivar.eq.3)then
              if(IsNext)then
                MR_vy_metP_last = MR_vy_metP_next
                MR_vy_metP_next = MR_dum3d_metP
              else
                MR_vy_metP_last = MR_dum3d_metP
              endif
            endif
          endif
        endif
      endif
        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp)); tmp_regrid2d_sp(:,:)=0.0_sp

      do k=1,nz_comp
        if(MR_iwindformat.eq.1.or.(MR_iwindformat.eq.2.and.MR_nSnd_Locs.eq.1))then
          ! ASCII profile data with one location
          MR_dum3d_compH(:,:,k) = MR_dum3d_metH(1,1,k)
          cycle
        elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2.and.MR_nSnd_Locs.gt.1)then
          ! ASCII profile data with multiple locations
          call MR_Regrid_MetSonde2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                       nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
        else
          ! All other cases
          call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                  nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
        endif
        do i = 1,nx_comp
          do j = 1,ny_comp
            if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          enddo
        enddo
        MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
      enddo

      deallocate(tmp_regrid2d_sp)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Read_3d_Met_Variable_to_CompH

!##############################################################################
!
!     MR_Read_3d_Met_Variable_to_CompP
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid along pressure levels.
!     This is done by calling:
!       Read_3d_MetP_Variable                   (gets variable on Met_x,Met_y, Met_p)
!        -> Read_3d_MetP_Variable               (gets variable on native Met subgrid)
!            -> Read_3d_MetP_Variable_[format]  (direct read of variable in
!                                                whatever format :
!                                                NetCDF, GRIB1/2, ASCII)
!
!     Takes as input :: ivar      :: specifies which variable to read
!                       istep     :: specified the met step
!                       IsNext    :: T for saving velocity values
!                                    F default
!     Sets  : MR_dum3d_compP
!               MR_dum3d_metP is also filled in the course of
!               generating MR_dum3d_compP
!
!##############################################################################

      subroutine MR_Read_3d_Met_Variable_to_CompP(ivar,istep,IsNext)

      integer,intent(in)           :: ivar
      integer,intent(in)           :: istep
      logical, optional,intent(in) :: IsNext

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_3d_Met_Variable_to_CompP             ----------"
        write(outlog(io),*)ivar,istep,IsNext
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times(this check is needed)

      if(.not.MR_useCompP)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Trying to run a CompP subroutine"
          write(errlog(io),*)"          with MR_useCompP = .false."
        endif;enddo
        stop 0
      endif

        ! First get the variable on the pressure levels
      call MR_Read_3d_MetP_Variable(ivar,istep)

      if(MR_Save_Velocities)then
        ! Note: this flag for saving the velocity values is useful in special
        ! cases such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then
          ! MR_dum3d_metP still contains the variable just read
          if(ivar.eq.2)then
            if(IsNext)then
              MR_vx_metP_last = MR_vx_metP_next
              MR_vx_metP_next = MR_dum3d_metP
            else
              MR_vx_metP_last = MR_dum3d_metP
            endif
          elseif(ivar.eq.3)then
            if(IsNext)then
              MR_vy_metP_last = MR_vy_metP_next
              MR_vy_metP_next = MR_dum3d_metP
            else
              MR_vy_metP_last = MR_dum3d_metP
            endif
          endif
        endif
      endif
        ! Now we have MR_dum3d_metP; interpolate onto computational grid
        !  Since MR_dum3d_metP and MR_dum3d_compP have the same
        !  pressure level, we only
        !  need to do 2d regridding on each p-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp

      do k=1,np_fullmet
        if(MR_iwindformat.eq.1.or.(MR_iwindformat.eq.2.and.MR_nSnd_Locs.eq.1))then
          ! ASCII profile data with one location
          MR_dum3d_compP(:,:,k) = MR_dum3d_metP(1,1,k)
          cycle
        elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2.and.MR_nSnd_Locs.gt.1)then
          ! ASCII profile data with multiple locations
          call MR_Regrid_MetSonde2Comp(nx_submet,ny_submet, MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                       nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          !stop 3
        else
          ! All other cases
          call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                  nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
        endif
        do i = 1,nx_comp
          do j = 1,ny_comp
            if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          enddo
        enddo
        MR_dum3d_compP(:,:,k) = tmp_regrid2d_sp(:,:)
      enddo

      deallocate(tmp_regrid2d_sp)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Read_3d_Met_Variable_to_CompP

!##############################################################################
!
!     MR_Read_2d_Met_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_met or MR_dum2d_met_int
!
!##############################################################################

      subroutine MR_Read_2d_Met_Variable(ivar,istep)

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep   

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_2d_Met_Variable                      ----------"
        write(outlog(io),*)ivar,istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        !call Read_1d_windfile
      case(2)
        !call Read_2d_ASCII
      case (3:5)
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif
      case default

      end select

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_2d_Met_Variable

!##############################################################################
!
!     MR_Read_2d_Met_Variable_to_CompGrid
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_2d_Met_Variable                   (gets variable on Met_x, Met_y)
!          -> Read_Met_Variable_[format]       (direct read of variable in
!                                                whatever format :
!                                                NetCDF, GRIB1/2, ASCII)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_comp
!               MR_dum2d_met and MR_dum2d_met are also filled in the course of
!               generating MR_dum2d_comp
!
!##############################################################################

      subroutine MR_Read_2d_Met_Variable_to_CompGrid(ivar,istep)

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Read_2d_Met_Variable_to_CompGrid          ----------"
        write(outlog(io),*)ivar,istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      call MR_Read_2d_Met_Variable(ivar,istep)

      allocate(tmp_regrid2d_sp(nx_comp,ny_comp)); tmp_regrid2d_sp(:,:)=0.0_sp

      if(MR_iwindformat.eq.1.or.(MR_iwindformat.eq.2.and.MR_nSnd_Locs.eq.1))then
        ! ASCII profile data with one location
        tmp_regrid2d_sp(:,:) = MR_dum2d_met(1,1)
      elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2.and.MR_nSnd_Locs.gt.1)then
        ! ASCII profile data with multiple locations
        call MR_Regrid_MetSonde2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                     nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
      else
        ! All other cases
        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
      endif
      do i = 1,nx_comp
        do j = 1,ny_comp
          if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        enddo
      enddo
      MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Read_2d_Met_Variable_to_CompGrid

!##############################################################################
!
!     MR_Rotate_UV_GR2ER_Met
!
!     This subroutine reads U and V on the metP grid, rotates to EarthRelative
!
!     Takes as input :: istep   :: specified the met step
!                       SetComp :: (optional) logical flag to regrid to comp grid
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_GR2ER_Met(istep,SetComp,IsNext)

      integer,intent(in)  :: istep
      logical,optional,intent(in) :: SetComp
      logical,optional,intent(in) :: IsNext

      integer             :: i,j,k

      real(kind=sp)                :: rotang
      real(kind=sp),dimension(2)   :: v_old, v_new
      real(kind=sp),dimension(2,2) :: rotmat

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Rotate_UV_GR2ER_Met                       ----------"
        write(outlog(io),*)istep,SetComp
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      ! Rotate wind vectors on the MetP grid for NARR cases with Map_Case = 2 (both projected)
      ! or rotate to LL on MetP grid for Map_Case = 4 (Met=proj, Comp=LL)

      call MR_Read_3d_MetP_Variable(2,istep)
        MR_u_ER_metP = MR_dum3d_metP
      call MR_Read_3d_MetP_Variable(3,istep)
        MR_v_ER_metP = MR_dum3d_metP
      ! Using these earth relative velocities, rotate to grid relative and
      ! copy to MR_dum3d_metP
      do i=1,nx_submet
        do j=1,ny_submet
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids
          rotang = -1.0_sp*real(theta_Met(i,j),kind=sp)
          rotmat(1,1) = cos(rotang)
          rotmat(2,1) = sin(rotang)
          rotmat(1,2) = rotmat(2,1)*(-1.0_sp)
          rotmat(2,2) = rotmat(1,1)
          do k=1,np_fullmet
            v_old(1:2) = (/MR_u_ER_metP(i,j,k),MR_v_ER_metP(i,j,k)/)
            v_new = matmul(rotmat(1:2,1:2),v_old(1:2))
            MR_u_ER_metP(i,j,k) = v_new(1)
            MR_v_ER_metP(i,j,k) = v_new(2)
          enddo
        enddo
      enddo

      if(MR_Save_Velocities)then
        ! Note: this flag for saving the velocity values is useful in special
        ! cases such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then ! First check if this parameter was provided
          if(IsNext) then  ! second, check if it is true
            ! MR_dum3d_metP still contains the variable just read
            if(IsNext)then
              MR_vx_metP_last = MR_vx_metP_next
              MR_vx_metP_next = MR_u_ER_metP
            else
              MR_vx_metP_last = MR_u_ER_metP
            endif
            if(IsNext)then
              MR_vy_metP_last = MR_vy_metP_next
              MR_vy_metP_next = MR_v_ER_metP
            else
              MR_vy_metP_last = MR_v_ER_metP
            endif
          endif
        endif
      endif

      if(present(SetComp)) then
        if(SetComp)then
          if(MR_useCompH)then
            ! copy v to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_v_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompH(istep)
            MR_dum3d_compH_2(1:nx_comp,1:ny_comp,1:nz_comp) = &
              MR_dum3d_compH(1:nx_comp,1:ny_comp,1:nz_comp)
  
            ! copy u to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_u_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompH(istep)
            !MR_dum3d_compH is set by this regrid call
          endif

          if(MR_useCompP)then
            ! copy v to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_v_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompH(istep)
            MR_dum3d_compP_2(1:nx_comp,1:ny_comp,1:np_fullmet) = &
              MR_dum3d_compP(1:nx_comp,1:ny_comp,1:np_fullmet)
            ! copy u to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_u_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompH(istep)
            !MR_dum3d_compP is set by this regrid call
          endif

        endif
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Rotate_UV_GR2ER_Met


!##############################################################################
!
!     MR_Rotate_UV_ER2GR_Comp
!
!     This subroutine is needed when we have a Earth-relative Met grid (either
!     natively or through MR_Rotate_UV_GR2ER_Met) and need to rotate that into a
!     projected Comp grid.  This corresponds to:
!       Map_Case = 3 : Met=LL,   Comp=proj
!       Map_Case = 5 : Met=proj, Comp=proj (different)
!     This subroutine reads U and V on the metP projected grid,
!     then resamples onto the projected compH grid, rotates to grid relative and
!     returns to calling program
!
!     Takes as input :: istep :: specified the met step
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_ER2GR_Comp(istep)

      integer,intent(in)  :: istep

      integer             :: i,j,k

      real(kind=sp)                :: rotang
      real(kind=sp),dimension(2)   :: v_old, v_new
      real(kind=sp),dimension(2,2) :: rotmat

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Rotate_UV_ER2GR_Met                       ----------"
        write(outlog(io),*)istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      if(Map_Case.eq.3.or.Map_Case.eq.4)then
        ! Met grid is natively LL and Comp grid is projected
        !  - or -
        ! Met grid is projected and comp grid is LL
        ! Fill the y velocities to spare array (comp_2)
        call MR_Read_3d_MetP_Variable(3,istep) ! This fills MR_dum3d_metP
        call MR_Regrid_MetP_to_CompH(istep) ! Takes MR_dum3d_metP and fills MR_dum3d_compH
          MR_dum3d_compH_2 = MR_dum3d_compH

        call MR_Read_3d_MetP_Variable(2,istep) 
        call MR_Regrid_MetP_to_CompH(istep)
        ! Now compH and compH_2 have vx and vy
      elseif(Map_Case.eq.5)then
        ! If the velocity components are grid relative, then
        ! assume MR_u_ER_metP and MR_v_ER_metP have already been set by
        ! MR_Rotate_UV_GR2ER_Met
        MR_dum3d_metP = MR_v_ER_metP
        call MR_Regrid_MetP_to_CompH(istep)
          MR_dum3d_compH_2 = MR_dum3d_compH
        MR_dum3d_metP = MR_u_ER_metP
        call MR_Regrid_MetP_to_CompH(istep)
        ! Now compH and compH_2 have vx and vy
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Should not be calling this subroutine unless Map_Case=3 or 5"
        endif;enddo
        stop 1
      endif

      ! Now loop through the comp points and rotate the vectors in place
      do i=1,nx_comp
        do j=1,ny_comp
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids
          rotang = real(theta_Comp(i,j),kind=sp)
          rotmat(1,1) = cos(rotang)
          rotmat(2,1) = sin(rotang)
          rotmat(1,2) = rotmat(2,1)*(-1.0_sp)
          rotmat(2,2) = rotmat(1,1)
          do k=1,nz_comp
            v_old(1:2) = (/MR_dum3d_compH(i,j,k),MR_dum3d_compH_2(i,j,k)/)
            v_new = matmul(rotmat(1:2,1:2),v_old(1:2))
            MR_dum3d_compH(i,j,k) = v_new(1)
            MR_dum3d_compH_2(i,j,k) = v_new(2)
          enddo
        enddo
      enddo

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Rotate_UV_ER2GR_Comp

!##############################################################################
!
!     MR_Regrid_MetP_to_CompH
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH by the subroutine Regrid_MetP_to_MetH.
!     MR_dum3d_metH is then regridded onto MR_dum3d_compH
!
!     Takes as input :: istep     :: specified the met step
!
!     Sets  : MR_dum3d_compH and MR_dum3d_compP
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_CompH(istep)

      integer,intent(in)  :: istep

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp
      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Regrid_MetP_to_CompH                      ----------"
        write(outlog(io),*)istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      if(MR_useCompH)then
        ! convert MR_dum3d_MetP to MR_dum3d_metH
        call MR_Regrid_MetP_to_MetH(istep)
      endif

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice or each p-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp
      if(MR_useCompH)then
        ! Regridding to CompH
        do k=1,nz_comp
          !if(MR_iwindformat.eq.1)then
          if(MR_iwindformat.eq.1.or.(MR_iwindformat.eq.2.and.MR_nSnd_Locs.eq.1))then
            ! ASCII profile data with one location
            MR_dum3d_compH(1:nx_comp,1:ny_comp,k) = MR_dum3d_metH(1,1,k)
            cycle
          elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2.and.MR_nSnd_Locs.gt.1)then
            ! ASCII profile data with multiple locations
            call MR_Regrid_MetSonde2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                         nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          else
            ! All other cases
            call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                    nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          endif
          do i = 1,nx_comp
            do j = 1,ny_comp
              if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
            enddo
          enddo
          MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
        enddo
      endif

      if(MR_useCompP)then
        ! Regridding to CompP
        do k=1,np_fullmet
          if(MR_iwindformat.eq.1.or.(MR_iwindformat.eq.2.and.MR_nSnd_Locs.eq.1))then
            ! ASCII profile data with one location
            MR_dum3d_compP(:,:,k) = MR_dum3d_metP(1,1,k)
            cycle
          elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2.and.MR_nSnd_Locs.gt.1)then
            ! ASCII profile data with multiple locations
            call MR_Regrid_MetSonde2Comp(nx_submet,ny_submet, MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                         nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          else
            ! All other cases
            call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                    nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          endif
          do i = 1,nx_comp
            do j = 1,ny_comp
              if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
            enddo
          enddo
          MR_dum3d_compP(:,:,k) = tmp_regrid2d_sp(:,:)
        enddo
      endif

      deallocate(tmp_regrid2d_sp)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Regrid_MetP_to_CompH

!##############################################################################
!
!     MR_Regrid_MetP_to_MetH
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH .
!
!     Takes as input :: istep     :: specified the met step
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_MetH(istep)


      integer,intent(in)  :: istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext
      real(kind=sp),dimension(:),allocatable :: dumVertCoord_sp
      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Regrid_MetP_to_MetH                       ----------"
        write(outlog(io),*)istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp

      allocate(  z_col_metP(np_fullmet+2));   z_col_metP(:)    = 0.0_sp
      allocate(var_col_metP(np_fullmet+2)); var_col_metP(:)    = 0.0_sp
      allocate(var_col_metH(nz_comp))     ; var_col_metH(:)    = 0.0_sp
      allocate(dumVertCoord_sp(nz_comp))  ; dumVertCoord_sp(:) = 0.0_sp

      ! Create 1d arrays in p, and regrid them into 1d arrays in z or s
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          if(istep.eq.MR_iMetStep_Now)then
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          else
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          endif
          ! Loop through the nodes and reset any that are negative or non-acsending in p
          ! This is usually never invoked, but would only affect potentially the bottom few nodes
          ! The intent is to maintain a monotonic GPH, while getting irrelevant nodes out of the
          ! way. 
          do k=1,np_fullmet
            if(z_col_metP(k).le.real(k,kind=sp)*MR_MIN_DZ) z_col_metP(k)=real(k,kind=sp)*MR_MIN_DZ
          enddo
          var_col_metP(1)              = 0.0_sp
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"MR_WARNING: Need to adjust value in MR_Regrid_MetP_to_MetH"
                write(outlog(io),*)"            that are not increasing in height with decreasing p"
                write(outlog(io),*)"i = ",i
                write(outlog(io),*)"j = ",j
                write(outlog(io),*)"k = ",k
                write(outlog(io),*)z_col_metP
              endif;enddo
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while (z_col_metP(knext)<z_col_metP(k-1))
                 knext=knext+1
              enddo
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              enddo
            endif
          enddo
            ! Assign top node (Should these indicies be +1?)
          z_col_metP(np_fullmet+2) = max(z_comp_sp(nz_comp)+1.0_sp,z_col_metP(np_fullmet+1)+1.0_sp)
          var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)

          if(MR_ZScaling_ID.eq.0)then
            dumVertCoord_sp(1:nz_comp) = z_comp_sp(1:nz_comp)
          elseif(MR_ZScaling_ID.eq.1)then
            dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + s_comp_sp(1:nz_comp)
          elseif(MR_ZScaling_ID.eq.2)then
              ! this recovers the real-world z coordinate from the sigma level and topography
            if(MR_iwind.eq.1)then
              ! this is incorrect (should be MR_jacob_met(i,j)), but this isn't set for iwind=1
              dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + &
                                             s_comp_sp(1:nz_comp) * MR_jacob_comp(1,1)
            else
              dumVertCoord_sp(1:nz_comp) = MR_Topo_met(i,j) + &
                                             s_comp_sp(1:nz_comp) * MR_jacob_met(i,j)
            endif
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid


          do io=1,MR_nio;if(VB(io).le.verbosity_debug2)then
            write(outlog(io),*)
            write(outlog(io),*)"Callling MR_Regrid_P2H_linear from MR_Regrid_MetP_to_MetH for ",i,j
          endif;enddo

          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,       var_col_metP, &
                                    nz_comp,      dumVertCoord_sp,  var_col_metH)
          MR_dum3d_metH(i,j,:) = var_col_metH

        enddo ! j
      enddo  ! i

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Regrid_MetP_to_MetH

!##############################################################################
!
!     MR_Regrid_Met2d_to_Comp2d
!
!     This subroutine expects the calling program to populate MR_dum2d_met.
!     MR_dum2d_met is then regridded onto MR_dum2d_comp
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum2d_comp
!
!##############################################################################

      subroutine MR_Regrid_Met2d_to_Comp2d

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------      MR_Regrid_Met2d_to_Comp2d                    ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp
      call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                              nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
      do i = 1,nx_comp
        do j = 1,ny_comp
          if(ieee_is_nan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        enddo
      enddo
      MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      return

      end subroutine MR_Regrid_Met2d_to_Comp2d

!##############################################################################
!
!     MR_Set_LL_mapping
!
!     This subroutine calculates the lon/lat coordinates for each point
!     of the Met grid.
!
!     Sets: 
!           MR_Have_LL_mapping
!           MR_xy2ll_xlon
!           MR_xy2ll_ylat
!           MR_lonmin,MR_lonmax
!           MR_latmin,MR_latmax
!
!##############################################################################

      subroutine MR_Set_LL_mapping

      use projection,      only : &
           PJ_Set_Proj_Params,&
           PJ_proj_for,&
           PJ_proj_inv

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: io
      integer :: i,j
      real(kind=dp) :: xin,yin
      real(kind=dp) :: xout,yout

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                 MR_Set_LL_mapping                 ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      ! This bit is needed to map topo values onto met grid points
      if(.not.IsLatLon_MetGrid.and. &
         .not.MR_Have_LL_mapping)then
#ifdef USEPOINTERS
        if(.not.associated(MR_xy2ll_xlon))allocate(MR_xy2ll_xlon(nx_submet,ny_submet))
        if(.not.associated(MR_xy2ll_ylat))allocate(MR_xy2ll_ylat(nx_submet,ny_submet))
#else
        if(.not.allocated(MR_xy2ll_xlon))allocate(MR_xy2ll_xlon(nx_submet,ny_submet))
        if(.not.allocated(MR_xy2ll_ylat))allocate(MR_xy2ll_ylat(nx_submet,ny_submet))
#endif
        MR_latmax =  -90.0_sp
        MR_latmin =   90.0_sp
        MR_lonmin =  360.0_sp
        MR_lonmax = -360.0_sp
        do i=1,nx_submet
          do j=1,ny_submet
            xin = real(x_submet_sp(i),kind=dp)  ! Projection routines use kind=8
            yin = real(y_submet_sp(j),kind=dp)
            call PJ_proj_inv(xin,yin, &
                           Met_iprojflag, Met_lam0,Met_phi0,Met_phi1,Met_phi2, &
                           Met_k0,Met_Re, &
                           xout,yout)
            if(xout.lt.MR_lonmin)MR_lonmin=real(xout,kind=sp)
            if(xout.gt.MR_lonmax)MR_lonmax=real(xout,kind=sp)
            if(yout.lt.MR_latmin)MR_latmin=real(yout,kind=sp)
            if(yout.gt.MR_latmax)MR_latmax=real(yout,kind=sp)
            MR_xy2ll_xlon(i,j) = real(xout,kind=sp)
            MR_xy2ll_ylat(i,j) = real(yout,kind=sp)
            if(MR_xy2ll_xlon(i,j).lt.0.0_sp) MR_xy2ll_xlon(i,j) = MR_xy2ll_xlon(i,j) + 360.0_sp
!            write(*,*)i,j,MR_xy2ll_xlon(i,j),MR_xy2ll_ylat(i,j)
          enddo
        enddo
        MR_Have_LL_mapping = .true.
      endif

      end subroutine MR_Set_LL_mapping

!##############################################################################














!##############################################################################
!
!     MR_DelMetP_Dx
!
!     Calculated the x derivative of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!     The unit of dx is always km.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dx

      integer :: i
      integer :: lside,rside
      real(kind=sp) :: dx_fac

      integer :: io                           ! Index for output streams

      do i=1,nx_submet
        if(i.eq.1)then
          if(IsPeriodic_CompGrid)then
            ! Two-sided derivative wrapping around
            lside  = nx_submet
            rside  = i+1
            dx_fac = 2.0_sp
          else
            ! One-sided derivative on left
            lside  = i+0
            rside  = i+1
            dx_fac = 1.0_sp
          endif
        elseif(i.eq.nx_submet)then
          if(IsPeriodic_CompGrid)then
            ! Two-sided derivative wrapping around
            lside  = i-1
            rside  = 1
            dx_fac = 2.0_sp
          else
            ! One-sided derivative on right
            lside  = i-1
            rside  = i+0
            dx_fac = 1.0_sp
          endif
        else
          ! Two-sided derivative
          lside  = i-1
          rside  = i+1
          dx_fac = 2.0_sp
        endif
        if(IsLatLon_MetGrid)then
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(rside,:,:)  - &
                                     MR_dum3d_metP(lside,:,:)) / &
                                    (rdlambda_MetP_sp(i,:,:)   * &
                                    dx_fac)
          else
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(rside,:,:)  - &
                                     MR_dum3d_metP(lside,:,:)) / &
                                    (0.5_sp*(rdlambda_MetP_sp(i,:,:)+rdlambda_MetP_sp(rside,:,:)) * &
                                    dx_fac)
          endif
        else
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(rside,:,:)  - &
                                     MR_dum3d_metP(lside,:,:)) / &
                                    (dx_fac*MR_dx_submet(i))
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"Need to fix DelMetP_Dx for non-regular grids."
            endif;enddo
            stop 1
          endif
        endif

      enddo

      end subroutine MR_DelMetP_Dx

!##############################################################################
!
!     MR_DelMetP_Dy
!
!     Calculated the y derivative of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!     The unit of dy is always km.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dy

      integer :: i,j
      integer :: lside,rside
      real(kind=sp) :: dy_fac

      integer :: io                           ! Index for output streams

      do j=1,ny_submet
        if(j.eq.1)then
          ! One-sided derivative on left
          lside=j+0
          rside=j+1
          dy_fac = 1.0_sp
        elseif(j.eq.ny_submet)then
          ! One-sided derivative on right
          lside=j-1
          rside=j+0
          dy_fac = 1.0_sp
        else
          ! Two-sided derivative
          lside=j-1
          rside=j+1
          dy_fac = 2.0_sp
        endif
        if(IsLatLon_MetGrid)then
          if(IsRegular_MetGrid)then
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,rside,:)  - &
                                       MR_dum3d_metP(i,lside,:)) / &
                                      (rdphi_MetP_sp(j,:)   * &
                                       dy_fac)
            enddo
          else
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,rside,:)  - &
                                       MR_dum3d_metP(i,lside,:)) / &
                                      (0.5_sp*(rdphi_MetP_sp(j,:)+rdphi_MetP_sp(rside,:))   * &
                                       dy_fac)
            enddo
          endif
        else
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(:,j,:) = (MR_dum3d_metP(:,rside,:)  - &
                                     MR_dum3d_metP(:,lside,:)) / &
                                    (dy_fac*MR_dy_submet(j))
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"Need to fix DelMetP_Dy for non-regular grids."
            endif;enddo
            stop 1
          endif
        endif
      enddo

      end subroutine MR_DelMetP_Dy

!##############################################################################
!
!  Returns the temperature (in K) given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!
!##############################################################################

      function MR_Temp_US_StdAtm(zin)

      real(kind=sp) :: MR_Temp_US_StdAtm
      real(kind=sp) :: zin

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_Tnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_temp

      integer :: k,kk
      integer :: io                           ! Index for output streams

      US_StdAtm_znodes = (/ 0.0_sp, 11.0_sp, 20.0_sp, 32.0_sp, &
                           47.0_sp, 51.0_sp, 71.0_sp, 84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.
      US_StdAtm_Tnodes = (/288.15_sp, 216.65_sp, 216.65_sp, 228.65_sp, &
                           270.65_sp, 270.65_sp, 214.65_sp, 186.946_sp, &
                           183.80_sp, 190.30_sp, 203.50_sp, 228.00_sp,  &
                           265.50_sp, 317.10_sp, 380.60_sp, 452.30_sp/)

      if(zin.le.US_StdAtm_znodes(1))then
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(1)
      elseif(zin.ge.US_StdAtm_znodes(MAXATMOSNODES))then
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(MAXATMOSNODES)
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Stopping in Temp_US_StdAtm"
        endif;enddo
        stop 1
      else
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        do k = MAXATMOSNODES-1,1,-1
          if(zin.ge.US_StdAtm_znodes(k).and.zin.lt.US_StdAtm_znodes(k+1)) kk = k
        enddo
        frac = (zin-US_StdAtm_znodes(kk)) / &
               (US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk))
        Delta_temp = US_StdAtm_Tnodes(kk+1)-US_StdAtm_Tnodes(kk)
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(kk) + Delta_temp * frac
      endif

      return

      end function MR_Temp_US_StdAtm

!##############################################################################
!
!  Returns the height in km given the pressure in hPa
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!       supplemented by doi=10.1.1.455.9475
!
!##############################################################################

      function MR_Z_US_StdAtm(pin)

      real(kind=sp) :: MR_Z_US_StdAtm  ! in km
      real(kind=sp) :: pin             ! in mb

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_pnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_z
      integer :: k,kk
      integer :: io                           ! Index for output streams

      US_StdAtm_znodes = (/ 0.0_sp,  11.0_sp,  20.0_sp,  32.0_sp,   &
                           47.0_sp,  51.0_sp,  71.0_sp,  84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.
      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 6.69389e-1_sp, 3.95642e-2_sp, 3.68501e-3_sp, &
                           1.795e-3_sp,7.345e-4_sp,3.090e-4_sp,1.422e-4_sp, &
                           7.362e-5_sp,4.236e-5_sp,2.667e-5_sp,2.666e-5_sp/)

      if(pin.ge.US_StdAtm_pnodes(1))then
        MR_Z_US_StdAtm = US_StdAtm_pnodes(1)
      elseif(pin.le.US_StdAtm_pnodes(MAXATMOSNODES))then
        MR_Z_US_StdAtm = US_StdAtm_pnodes(MAXATMOSNODES)
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Pressure provided is below lowest value of US StdAtm."
          write(errlog(io),*)"Stopping in MR_Z_US_StdAtm"
        endif;enddo
        stop 1
      else
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        do k = MAXATMOSNODES-1,1,-1
          if(pin.le.US_StdAtm_pnodes(k).and.pin.gt.US_StdAtm_pnodes(k+1)) kk = k
        enddo
        frac = (pin-US_StdAtm_pnodes(kk)) / &
               (US_StdAtm_pnodes(kk+1)-US_StdAtm_pnodes(kk))
        Delta_z = US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk)
        ! Note: this should really be an exponential fit, not linear
        MR_Z_US_StdAtm = US_StdAtm_znodes(kk) + Delta_z * frac
      endif

      return

      end function MR_Z_US_StdAtm

!##############################################################################
!
!  Returns the pressure in hPa given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!       supplemented by doi=10.1.1.455.9475
!
!##############################################################################

      function MR_Pres_US_StdAtm(zin)

      real(kind=sp) :: MR_Pres_US_StdAtm
      real(kind=sp) :: zin

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_pnodes

      real(kind=sp) :: pres0 = 1013.25_sp
      real(kind=sp) :: skinz   = 7.0_sp

      US_StdAtm_znodes = (/ 0.0_sp,  11.0_sp,  20.0_sp,  32.0_sp,   &
                           47.0_sp,  51.0_sp,  71.0_sp,  84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.

      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 6.69389e-1_sp, 3.95642e-2_sp, 3.68501e-3_sp, &
                           1.795e-3_sp,7.345e-4_sp,3.090e-4_sp,1.422e-4_sp, &
                           7.362e-5_sp,4.236e-5_sp,2.667e-5_sp,2.666e-5_sp/)

        ! Eq 1.8 of Wallace and Hobbs
        !  This starts to fall apart around 90 km or so
      MR_Pres_US_StdAtm = pres0 * exp(-zin/skinz)

      return

      end function MR_Pres_US_StdAtm

!##############################################################################

!##############################################################################
!
!    MR_QC_3dvar
!
!##############################################################################

      subroutine MR_QC_3dvar(              &
                 ivar,                     &
                 nx_max,ny_max,nz1_max,    &
                 z_array_sp,               &
                 nz2_max,                  &
                 dum_array_sp,             &
                 fill_val_sp,              &
                 bc_low_sp, bc_high_sp)

      integer               ,intent(in)    :: ivar
      integer               ,intent(in)    :: nx_max,ny_max,nz1_max
      real(kind=sp)         ,intent(in)    :: z_array_sp(nx_max,ny_max,nz1_max)
      integer               ,intent(in)    :: nz2_max
      real(kind=sp)         ,intent(inout) :: dum_array_sp(nx_max,ny_max,nz1_max)
      real(kind=sp)         ,intent(in)    :: fill_val_sp
      real(kind=sp),optional,intent(in)    :: bc_low_sp
      real(kind=sp),optional,intent(in)    :: bc_high_sp

      logical,dimension(nz1_max) :: IsFillValue
      logical,dimension(nz1_max) :: InterpolateLev
      logical,dimension(nz1_max) :: ExtrapolateLev

      integer ::  i,j,k,kk,kkk,klow,khigh
      integer :: idx
      real(kind=sp) :: pp,pv,p1,p2,fac
      integer :: io                           ! Index for output streams

      InterpolateLev(:) = .false.
      ExtrapolateLev(:) = .false.

      ! First check if this variable is on the same pressure grid as GPH
      idx = Met_var_zdim_idx(ivar)  ! index of pressure level for this variable
      if(levs_code(idx).eq.1)then
        ! one-to-one mapping: no special action needed
      else if (levs_code(idx).eq.2)then
        ! truncated grid: this will be addressed at the end of this subroutine
      else if (levs_code(idx).eq.3)then
        ! interpolated:
        ! We need the result to live on the same pressure levels as the GPH array, so
        ! for each level of the primary pressure coordinate, we need to find if a level
        ! of the variable can be copied or must be interpolated
        !
        ! Loop from the top of the variable array (because nz2_max<nz1_max)
        InterpolateLev(:) = .true.
        ExtrapolateLev(:) = .false.
        do k=nz2_max,1,-1  ! these are the p-levels the variable lives on (nz2_max<nz1_max)
          pv = levs_fullmet_sp(idx,k)  ! Pressure level in question

          do kk=nz1_max,1,-1           ! Looping over destination pressure levels
            pp = levs_fullmet_sp(1,kk)
            if(abs(pv-pp).lt.MR_EPS_SMALL)then
              ! pv matched with pp
              InterpolateLev(kk) = .false.
              ExtrapolateLev(kk) = .false.
              ! if the slices are at the same index, do nothing
              if (k.ne.kk)then
                ! Move the slice to the needed position and flag original slice with NaN/Fill_Value
                dum_array_sp(1:nx_max,1:ny_max,kk) = dum_array_sp(1:nx_max,1:ny_max,k)
                dum_array_sp(1:nx_max,1:ny_max,k)  = fill_val_sp
              endif
            elseif(levs_fullmet_sp(idx,1)-levs_fullmet_sp(1,kk).lt.MR_EPS_SMALL)then
                ! Check if this kk level of the master grid is less than the first level of the var grid
                ! pressure level is smaller (higher altitude) than top of GPH
              ExtrapolateLev(kk) = .true.
              InterpolateLev(kk) = .false.
            elseif(levs_fullmet_sp(idx,nz2_max)-levs_fullmet_sp(1,kk).gt.MR_EPS_SMALL)then
                ! This pressure level is greater (lower altitude) than bottom of GPH
              ExtrapolateLev(kk) = .true.
              InterpolateLev(kk) = .false.
            endif
          enddo

        enddo

        ! Now go to each of the layers we need to calculate, assuming we have the end layers
        do k=2,nz1_max-1
          if (InterpolateLev(k))then
            klow = 1
            do kk = k,1,-1
              if (.not.InterpolateLev(kk))then
                klow = kk
                exit
              endif
            enddo
            do kk = k,nz1_max
              if (.not.InterpolateLev(kk))then
                khigh = kk
                exit
              endif
            enddo
            pp = levs_fullmet_sp(1,k)     ! This is the pressure level we need to populate
            p1 = levs_fullmet_sp(1,klow)  ! pressure at low-altitude end of bracket
            p2 = levs_fullmet_sp(1,khigh) ! pressure at high-altitude end of bracket

            fac = (log(pp/p1))/(log(p2/p1))
            dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,klow) + &
                                                fac * (dum_array_sp(1:nx_max,1:ny_max,khigh) - &
                                                       dum_array_sp(1:nx_max,1:ny_max,klow))
          endif
        enddo
      else
        ! Should not be here!
        ! We assume that the variable lives on the GPH grid or a subset
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'MR ERROR: Variable has more levels than GPH'
        endif;enddo
        stop 1
      endif

      do i=1,nx_max
        do j=1,ny_max
          ! find all fill values in the column
          IsFillValue = .false.
          do k=1,nz2_max
            if(ieee_is_nan(dum_array_sp(i,j,k)).or.      &  ! Some windfiles use NaN's for fill values
               dum_array_sp(i,j,k).eq.fill_val_sp) &  ! Others have a specific number for fill
                   IsFillValue(k) = .true.
          enddo

          ! Set lower BC if requested and needed
          if(IsFillValue(1).and.present(bc_low_sp)) &
              dum_array_sp(i,j,1) = bc_low_sp

          if(nz2_max.gt.1)then
            ! Set upper BC if requested and needed
            if(IsFillValue(nz2_max).and.present(bc_high_sp)) &
                dum_array_sp(i,j,nz2_max) = bc_high_sp

            ! Now find the lowest non-Fill value
            do klow=1,nz2_max
              if(.not.IsFillValue(klow))exit
            enddo
            ! And the highest non-Fill value
            do khigh=nz2_max,1,-1
              if(.not.IsFillValue(khigh))exit
            enddo

            ! Set bottom values to the lowest real number
            !if(klow.gt.1.and..not.present(bc_low_sp))then
            if(klow.gt.1)then
              do k = 1,klow-1
                dum_array_sp(i,j,k) = dum_array_sp(i,j,klow)
              enddo
            endif

            ! Set top values to the highest real number
            !if(khigh.lt.nz2_max.and..not.present(bc_high_sp))then
            if(khigh.lt.nz2_max)then
              do k = nz2_max,khigh+1,-1
                dum_array_sp(i,j,k) = dum_array_sp(i,j,khigh)
              enddo
            endif

            ! Now check if there are any intermediate FillValues and
            ! interpolate
            do k=klow+1,khigh-1
              if(IsFillValue(k))then
                 ! linearly interpolate in z
                   ! Find first number bove
                 do kk = k+1,khigh,1
                   if(.not.IsFillValue(kk))exit
                 enddo
                   ! Find first number beneath
                 do kkk = max(k-1,1),klow,-1
                   if(.not.IsFillValue(kkk))exit
                 enddo
                 dum_array_sp(i,j,k) = dum_array_sp(i,j,kk) + &
                       (dum_array_sp(i,j,kk)-dum_array_sp(i,j,kkk)) * &
                       (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))/ &
                       (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))
              endif

            enddo
          endif
        enddo
      enddo

      if(levs_code(idx).eq.2)then
        ! For truncated grids, assign truncated (high altitude) values to the top bc, if
        ! present, or copy the upper-most valid value if not
        if(present(bc_high_sp))then
          do k=nz2_max+1,nz1_max
            dum_array_sp(1:nx_max,1:ny_max,k) = bc_high_sp
          enddo
        else
          do k=nz2_max+1,nz1_max
            dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,nz2_max)
          enddo
        endif
      elseif(levs_code(idx).eq.3)then
        ! Some interpolated grids are also truncated
        do k=1,nz1_max
          if(k.gt.nz1_max/2)then
            if(ExtrapolateLev(k))then
              if(present(bc_high_sp))then
                dum_array_sp(1:nx_max,1:ny_max,k) = bc_high_sp
              else
                dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,nz2_max)
              endif
            endif
          else
            if(ExtrapolateLev(k))then
              if(present(bc_low_sp))then
                dum_array_sp(1:nx_max,1:ny_max,k) = bc_low_sp
              else
                dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,1)
              endif
            endif
          endif
        enddo
      endif

      end subroutine MR_QC_3dvar

!##############################################################################
!
!    MR_Check_Prerequsites
!
!    Subroutine called at various points in MetReader to verify that the
!    prerequisite steps have been completed.  The five parameters are logicals
!    specifying whether or not a particular test is required.
!
!##############################################################################

      subroutine MR_Check_Prerequsites(test_allocate,    &
                                       test_dimvars,     &
                                       test_compproj,    &
                                       test_initmetgrid, &
                                       test_setmettimes)

      logical ,intent(in) :: test_allocate    ! Check if MR_Allocate_FullMetFileList was called
      logical ,intent(in) :: test_dimvars     ! Check if MR_Read_Met_DimVars
      logical ,intent(in) :: test_compproj    ! Check if MR_Set_CompProjection
      logical ,intent(in) :: test_initmetgrid ! Check if MR_Initialize_Met_Grids
      logical ,intent(in) :: test_setmettimes ! Check if MR_Set_Met_Times

      logical :: failtest_allocate    = .false.
      logical :: failtest_dimvars     = .false.
      logical :: failtest_compproj    = .false.
      logical :: failtest_initmetgrid = .false.
      logical :: failtest_setmettimes = .false.

      integer :: io                           ! Index for output streams

      ! Check prerequisites
      if(test_allocate.eqv..true.)then
        if(CALLED_MR_Allocate_FullMetFileList.eqv..false.)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  The subroutine MR_Allocate_FullMetFileList has not be called."
            write(errlog(io),*)"           This must first be called to allocate the space for the"
            write(errlog(io),*)"           windfiles.  The calling format is:"
            write(errlog(io),*)"             MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)"
            write(errlog(io),*)"               where iw      = windfile class (1-5)"
            write(errlog(io),*)"                     iwf     = windfile format number (1-51)"
            write(errlog(io),*)"                     igrid   = 2d grid code used (if known)"
            write(errlog(io),*)"                     idf     = data format ID (ASCII, NetCDF, GRIB, etc)"
            write(errlog(io),*)"                     iwfiles = number of windfiles to be read"
            write(errlog(io),*)" "
            write(errlog(io),*)"           The array MR_windfiles(1:iwfiles) must then be filled within"
            write(errlog(io),*)"           the calling program."
          endif;enddo
          failtest_allocate = .true.
        endif
      endif

      if(test_dimvars.eqv..true.)then
        if(CALLED_MR_Read_Met_DimVars.eqv..false.)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then     
            write(errlog(io),*)"MR ERROR:  The subroutine MR_Read_Met_DimVars has not be called."
            write(errlog(io),*)"           This must first be called to determine the size (x,y,z,t) of"
            write(errlog(io),*)"           the windfiles.  The calling format is:"
            write(errlog(io),*)"             MR_Read_Met_DimVars(optional_argument = iyear)"
            write(errlog(io),*)"               where the argument iyear determines nt for leap/non-leap years"
          endif;enddo
          failtest_dimvars = .true.
        endif
      endif

      if(test_compproj.eqv..true.)then
        if(CALLED_MR_Set_CompProjection.eqv..false.)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then     
            write(errlog(io),*)"MR ERROR:  The subroutine MR_Set_CompProjection has not be called."
            write(errlog(io),*)"           This must first be called to determine the subgrid needed"
            write(errlog(io),*)"           from the windfiles.  The calling format is:"
            write(errlog(io),*)"             MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)"
            write(errlog(io),*)"               where LL_flag = 0,1 (for projected or lon/lat, respectively)"
            write(errlog(io),*)"                     ipf     = 1-5 (projection code)"
            write(errlog(io),*)"                     lam0    = central longitude"
            write(errlog(io),*)"                     phi0    = first latitude for projection"
            write(errlog(io),*)"                     phi1    = second latitude for projection"
            write(errlog(io),*)"                     phi2    = third latitude for projection"
            write(errlog(io),*)"                     ko      = scaling term (<=1.0)"
            write(errlog(io),*)"                     Re      = radius of earth in km"
          endif;enddo
          failtest_compproj = .true.
        endif
      endif

      if(test_initmetgrid.eqv..true.)then
        if(CALLED_MR_Initialize_Met_Grids.eqv..false.)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  The subroutine MR_Initialize_Met_Grids has not be called."
            write(errlog(io),*)"           This must first be called to allocate the computational grid needed"
            write(errlog(io),*)"           and to generate the mapping from the windfiles.  The calling format is:"
            write(errlog(io),*)"             MR_Initialize_Met_Grids(nx,ny,nz,dumx_sp,dumy_sp,dumz_sp,periodic)"
            write(errlog(io),*)"               where nx,ny,nz = size of comp. grid in x,y,z"
            write(errlog(io),*)"                     dumx_sp  = cell-centered grid points in x"
            write(errlog(io),*)"                     dumy_sp  = cell-centered grid points in y"
            write(errlog(io),*)"                     dumz_sp  = cell-centered grid points in z"
            write(errlog(io),*)"                     periodic = true or false for periodic grids"
          endif;enddo
          failtest_initmetgrid = .true.
        endif
      endif

      if(test_setmettimes.eqv..true.)then
        if(CALLED_MR_Set_Met_Times.eqv..false.)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR:  The subroutine MR_Set_Met_Times has not be called."
            write(errlog(io),*)"           This must first be called to trim the list of windfiles to just"
            write(errlog(io),*)"           that needed by the simulation.  The calling format is:"
            write(errlog(io),*)"             MR_Set_Met_Times(eStartHour,Duration)"
            write(errlog(io),*)"               where eStartHour = start time needed in hours since basetime"
            write(errlog(io),*)"                     Duration   = length of time needed"
          endif;enddo
          failtest_setmettimes = .true.
        endif
      endif

      if(failtest_allocate    .or.&
         failtest_dimvars     .or.&
         failtest_compproj    .or.&
         failtest_initmetgrid .or.&
         failtest_setmettimes)then
        stop 1
      endif

      end subroutine MR_Check_Prerequsites

!##############################################################################
!
!    MR_FileIO_Error_Handler
!
!    Subroutine called whenever a line from a file is read with a non-zero
!    iostat.
!
!##############################################################################

      subroutine MR_FileIO_Error_Handler(ios,linebuffer080,iomessage)

      integer           ,intent(in) :: ios
      character(len= 80),intent(in) :: linebuffer080
      character(len=120),intent(in) :: iomessage

      integer :: io

      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        if(ios.lt.0)then
          write(errlog(io),*)'MR ERROR Reading from file:  EOF encountered',ios
        else
          write(errlog(io),*)'MR ERROR Reading line from file:  input line format error',ios
          write(errlog(io),*)linebuffer080
        endif
        write(errlog(io),*)'MR System Message: ',trim(adjustl(iomessage))
      endif;enddo
      stop 1

      end subroutine MR_FileIO_Error_Handler

!##############################################################################

      end module MetReader

