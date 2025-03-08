!##############################################################################
!##############################################################################
!##############################################################################
! Note: you can inspect the GRIB header using grib_dump tmp.grib2 > header.txt


!     MR_Read_Met_DimVars_GRIB
!     MR_Read_Met_Times_GRIB
!     MR_Read_MetP_Variable_GRIB
!     MR_GRIB_check_status

!##############################################################################
!
!     MR_Read_Met_DimVars_GRIB
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_GRIB_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_GRIB

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,verbosity_production,&
         MR_MAXVARS,x_fullmet_sp,y_fullmet_sp,nx_fullmet,ny_fullmet,&
         Met_var_zdim_idx,nlevs_fullmet,levs_code,levs_fullmet_sp,&
         p_fullmet_sp,np_fullmet,z_approx,&
         MR_dx_met,MR_dy_met,dx_met_const,dy_met_const,IsLatLon_MetGrid,IsRegular_MetGrid,&
         Met_iprojflag,Met_k0,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_Re,MR_EPS_SMALL,&
         MR_GRIB_Version,MR_iwind,MR_iwindformat,MR_Max_geoH_metP_predicted,&
         x_inverted,y_inverted,z_inverted,MR_windfiles,Met_var_IsAvailable,&
         nlev_coords_detected,nt_fullmet,Met_var_GRIB1_Param,Met_var_GRIB1_St,MR_GRIB_Version,&
         Met_var_GRIB2_DPcPnSt,Met_var_GRIB_names,&
           MR_Z_US_StdAtm        

      use projection,      only : &
           PJ_Set_Proj_Params,&
           PJ_proj_for

      use eccodes

      implicit none

      integer, parameter :: sp         = 4 ! single precision
      integer, parameter :: dp         = 8 ! double precision
      integer, parameter :: MAXGRIBREC = 10000

      integer :: i, j, k
      integer :: ir
      real(kind=sp) :: xLL_fullmet
      real(kind=sp) :: yLL_fullmet
      real(kind=sp) :: xUR_fullmet
      real(kind=sp) :: yUR_fullmet

      integer            :: ifile
      integer            :: igrib
      integer            :: iw
      integer,dimension(MAXGRIBREC) :: igribv
      integer            :: nSTAT

      integer          :: ivar,iivar
      integer          :: idx
      integer          :: maxdimlen

      character(len=130) :: grib_file_path
      integer(kind=4)  :: iv_discpl
      integer(kind=4)  :: iv_paramC
      integer(kind=4)  :: iv_paramN
      integer(kind=4)  :: iv_typeSf
      !integer(kind=4)  :: iv_Table
      character(len=3) :: iv_typeSfc

      integer(kind=4)  :: numberOfPoints
      integer(kind=4)  :: dum_int
      real(kind=8)     :: dum_dp
      character(len=20) :: dum_str
      integer          :: gg_order
      real(kind=dp) :: x_start,y_start
      real(kind=dp) :: Lon_start,Lat_start
      real(kind=dp) :: Lon_end,Lat_end
      !real(kind=dp),dimension(:),allocatable     :: lats,lons,values
      real(kind=dp),dimension(:),allocatable     :: values
      real(kind=dp), parameter :: tol = 1.0e-3_dp

      integer(kind=4)  :: typeOfFirstFixedSurface
      integer            :: count1=0
        ! Stores values of keys read from GRIB file
      character(len=6) :: grb_shortName
      character(len=36):: grb_longName
      character(len=4) :: grb_typeSfc
      integer(kind=4)  :: grb_discipline
      integer(kind=4)  :: grb_parameterCategory
      integer(kind=4)  :: grb_parameterNumber
      integer(kind=4)  :: grb_Table
      integer(kind=4)  :: grb_level
      integer(kind=4)  :: grb_scaledValueOfFirstFixedSurface

      integer(kind=4),dimension(MR_MAXVARS,100) :: zlev_dum  ! This will hold the z-levels, up to 100
      integer(kind=4),dimension(MR_MAXVARS)     :: zcount    ! This will hold the length of the z-coord
      logical :: Check
      logical :: FoundOldDim
      logical :: IsTruncatedDim
      logical :: ReadGrid
      integer(kind=4) :: kk,tmp1
      !integer :: stat
      logical :: IsNewLevel
      integer :: iz

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine MR_GRIB_check_status(nSTAT, errcode, operation)
          integer, intent(in) :: nSTAT
          integer, intent(in) :: errcode
          character(len=*), intent(in) :: operation
        end subroutine MR_GRIB_check_status
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                MR_Read_Met_DimVars_GRIB           ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      if(MR_iwind.eq.5)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR : ",&
              "GRIB reader not implemented for multi-timestep files."
          write(errlog(io),*)"         iwind=5 files are all multi-step"
        endif;enddo
        stop 1
      else
          !------------------------------------------------------------------------
          ! Start of block for all non-iwind=5 and non-iwf=50
          ! This is where the Netcdf and Grib subroutines can be compared
          !
          ! Checking for dimension length and values for x,y,t,p
          !   Assume all files have the same format
          ! Note: you can inspect the GRIB header using grib_dump tmp.grib2 > header.txt

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Opening GRIB file to find version number"
        endif;enddo
        iw = 1

        call codes_open_file(ifile,trim(adjustl(MR_windfiles(iw))),'R',nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_open_file ")
        call codes_new_from_file(ifile,igrib,CODES_PRODUCT_GRIB,nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")
        call codes_get(igrib,'editionNumber',MR_GRIB_Version,nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get editionNumber ")
        call codes_release(igrib,nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
        call codes_close_file(ifile,nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_close_file ")

      endif
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Grib version = ",MR_GRIB_Version
      endif;enddo
      !------------------------------------------------------------------------
      ! Checking for dimension length and values for x,y,t,p
      !   Assume all files have the same format
      maxdimlen = 0

      ! Loop through all the GRIB messages,
      ! If we find a message that matches a variable criterion, then log the level to 
      !  a dummy array.
      ! Finally sort the pressure values and evaluate the distinct pressure coordinates
      grib_file_path  = trim(adjustl(MR_windfiles(1)))

      call codes_open_file(ifile,grib_file_path,'R',nSTAT)
      if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_open_file ")

      count1=1
      call codes_new_from_file(ifile,igribv(count1),CODES_PRODUCT_GRIB,nSTAT)
      if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")

      do while (nSTAT.eq.CODES_SUCCESS)
        count1=count1+1
        if (count1.gt.MAXGRIBREC) then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: too many GRIB messages"
            write(errlog(io),*)"          current limit set to ",MAXGRIBREC
          endif;enddo
          stop 1
        endif
        call codes_new_from_file(ifile,igribv(count1),CODES_PRODUCT_GRIB,nSTAT)
      enddo
      count1 = count1-1
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"  Number of grib records found = ",count1
      endif;enddo
      zcount(:)     = 0
      zlev_dum(:,:) = 0
      do ir = 1,count1
        if(ir.eq.1)then
          ! For the first record, get the x,y grid info
          ! HFS: This is a bit of an assumption that the first grib record contains the grid
          !      information we need.
          ReadGrid = .false.
          call codes_get(igribv(ir),'Ni',nx_fullmet,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Ni ")
          call codes_get(igribv(ir),'Nj',ny_fullmet,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Nj ")
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))

          call codes_get(igribv(ir),'gridType',dum_str,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get gridType ")
          call codes_get(igribv(ir),'latitudeOfFirstGridPointInDegrees',dum_dp,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get latitudeOfFirstGridPointInDegrees ")
          Lat_start = dum_dp
          call codes_get(igribv(ir),'longitudeOfFirstGridPointInDegrees',dum_dp,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get longitudeOfFirstGridPointInDegrees ")
          Lon_start = dum_dp

          ! Check for end lat/lon in order to determine if either coordinate is
          ! inverted
          call codes_get(igribv(ir),'latitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)then
            !call MR_GRIB_check_status(nSTAT,0,"codes_get latitudeOfLastGridPointInDegrees ")
            ! assume it is not inverted
            y_inverted = .false.
          else
            Lat_end = dum_dp
            if(Lat_start.gt.Lat_end)then
              y_inverted = .true.
            else
              y_inverted = .false.
            endif
          endif
          call codes_get(igribv(ir),'longitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)then
            !call MR_GRIB_check_status(nSTAT,0,"codes_get longitudeOfLastGridPointInDegrees ")
            ! Note: ecmwf forecasts have longitudeOfFirstGridPointInDegrees = 180.0
            !                            longitudeOfLastGridPointInDegrees  = 179.75
            !       Grid is not actually inverted, just wraps around
            ! assume it is not inverted
            x_inverted = .false.
          else
            Lon_end = dum_dp
            if(Lon_start.gt.Lon_end)then
              x_inverted = .true.
            else
              x_inverted = .false.
            endif
          endif

          dum_int = 0
          Met_Re =  6371.229_8
          call codes_get(igribv(ir),'shapeOfTheEarth',dum_int,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get shapeOfTheEarth ")
          if (dum_int.eq.0)then
              ! 0  Earth assumed spherical with radius = 6,367,470.0 m
            Met_Re =  6367.470_8
          elseif(dum_int.eq.1)then
              ! 1  Earth assumed spherical with radius specified by data producer
              !  Try to read the radius of earth
              !  For now, just assign the default
            Met_Re =  6371.229_8
          elseif(dum_int.eq.2)then
              ! 2  Earth assumed oblate spheroid with size as determined by IAU in 1965
              !    (major axis = 6,378,160.0 m, minor axis = 6,356,775.0 m, f = 1/297.0)
            Met_Re =  6371.229_8
          elseif(dum_int.eq.3)then
              ! 3  Earth assumed oblate spheroid with major and minor axes specified by data producer
            Met_Re =  6371.229_8
          elseif(dum_int.eq.4)then
              ! 4  Earth assumed oblate spheroid as defined in IAG-GRS80 model 
              !    (major axis = 6,378,137.0 m, minor axis = 6,356,752.314 m, f = 1/298.257222101)
            Met_Re =  6371.229_8
          elseif(dum_int.eq.5)then
              ! 5  Earth assumed represented by WGS84 (as used by ICAO since 1998)
            Met_Re =  6371.229_8
          elseif(dum_int.eq.6)then
              ! 6  Earth assumed spherical with radius of 6,371,229.0 m
            Met_Re =  6371.229_8
          else
              ! 7-191 Reserved
              ! 192- 254  Reserved for local use
            Met_Re =  6371.229_8
          endif

          if(index(dum_str,'regular_ll').ne.0)then
            IsLatLon_MetGrid = .true.
            !Lat_start = y_start
            !Lon_start = x_start
            y_start = Lat_start
            x_start = Lon_start

            call codes_get(igribv(ir),'numberOfPoints',numberOfPoints,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
            allocate(values(numberOfPoints))
            call codes_get(igribv(ir),'values',values,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
            ReadGrid = .true.
            deallocate(values)
            !call codes_get(igribv(ir),'latitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
            !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get latitudeOfLastGridPointInDegrees ")
            !if(Lat_start.gt.dum_dp)then
            !  y_inverted = .true.
            !else
            !  y_inverted = .false.
            !endif
            !call codes_get(igribv(ir),&
            !      'longitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
            !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get longitudeOfLastGridPointInDegrees ")
            !if(Lon_start.gt.dum_dp)then
            !  x_inverted = .true.
            !else
            !  x_inverted = .false.
            !endif
            call codes_get(igribv(ir),'iDirectionIncrementInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get iDirectionIncrementInDegrees ")
            dx_met_const = real(dum_dp,kind=4)
            call codes_get(igribv(ir),'jDirectionIncrementInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get jDirectionIncrementInDegrees ")
            dy_met_const = real(dum_dp,kind=4)
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Setting x grid starting at",Lon_start
              write(outlog(io),*)"with a spacing of ",dx_met_const
            endif;enddo
            do i=1,nx_fullmet
              x_fullmet_sp(i) = real(Lon_start,kind=sp)+ &
                                 (i-1)* dx_met_const
            enddo
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Setting y grid starting at",Lat_start
              write(outlog(io),*)"with a spacing of ",dy_met_const
            endif;enddo
            if(y_inverted)then
              do j = 1,ny_fullmet
                y_fullmet_sp(j) = real(y_start - (j-1)*dy_met_const,kind=sp)
              enddo
            else
              do j = 1,ny_fullmet
                y_fullmet_sp(j) = real(y_start + (j-1)*dy_met_const,kind=sp)
              enddo
            endif
            !do j=1,ny_fullmet
            !  y_fullmet_sp(j) = real(Lat_start,kind=sp)+ &
            !                     (j-1)* dy_met_const
            !enddo
            x_fullmet_sp(0) = x_fullmet_sp(1)-dx_met_const
            x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(nx_fullmet)+dx_met_const
          elseif(index(dum_str,'regular_gg').ne.0)then
            IsLatLon_MetGrid = .true.
            ! For Gaussian grid, get the order
            call codes_get(igribv(ir),'N',gg_order,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get N ")
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Detected regular_gg grid of order N=",gg_order
            endif;enddo

            if(gg_order.eq.128)then
              ! this is the N128 grid (512x256) used by ERA Interim
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"Assuming this is ERA Interim; filling lon and lat"
              endif;enddo
              dx_met_const = 0.7031252_sp
              do i=1,nx_fullmet
                x_fullmet_sp(i) = 0.0_sp+(i-1)* dx_met_const
              enddo
              y_fullmet_sp(1:ny_fullmet) = &
         (/89.46282_sp, 88.76695_sp, 88.06697_sp, 87.36607_sp, 86.66480_sp, 85.96337_sp,              &
           85.26185_sp, 84.56026_sp, 83.85863_sp, 83.15699_sp, 82.45532_sp, 81.75363_sp, 81.05194_sp, &
           80.35023_sp, 79.64853_sp, 78.94681_sp, 78.24509_sp, 77.54337_sp, 76.84164_sp, 76.13991_sp, &
           75.43818_sp, 74.73644_sp, 74.03471_sp, 73.33297_sp, 72.63123_sp, 71.92949_sp, 71.22775_sp, &
           70.52601_sp, 69.82426_sp, 69.12252_sp, 68.42078_sp, 67.71903_sp, 67.01729_sp, 66.31554_sp, &
           65.61379_sp, 64.91205_sp, 64.21030_sp, 63.50855_sp, 62.80680_sp, 62.10506_sp, 61.40331_sp, &
           60.70156_sp, 59.99981_sp, 59.29806_sp, 58.59631_sp, 57.89456_sp, 57.19281_sp, 56.49106_sp, &
           55.78931_sp, 55.08756_sp, 54.38581_sp, 53.68406_sp, 52.98231_sp, 52.28056_sp, 51.57881_sp, &
           50.87706_sp, 50.17531_sp, 49.47356_sp, 48.77180_sp, 48.07005_sp, 47.36830_sp, 46.66655_sp, &
           45.96480_sp, 45.26305_sp, 44.56129_sp, 43.85954_sp, 43.15779_sp, 42.45604_sp, 41.75429_sp, &
           41.05254_sp, 40.35078_sp, 39.64903_sp, 38.94728_sp, 38.24553_sp, 37.54378_sp, 36.84202_sp, &
           36.14027_sp, 35.43852_sp, 34.73677_sp, 34.03502_sp, 33.33326_sp, 32.63151_sp, 31.92976_sp, &
           31.22800_sp, 30.52625_sp, 29.82450_sp, 29.12275_sp, 28.42099_sp, 27.71924_sp, 27.01749_sp, &
           26.31573_sp, 25.61398_sp, 24.91223_sp, 24.21048_sp, 23.50872_sp, 22.80697_sp, 22.10522_sp, &
           21.40347_sp, 20.70171_sp, 19.99996_sp, 19.29821_sp, 18.59645_sp, 17.89470_sp, 17.19295_sp, &
           16.49120_sp, 15.78944_sp, 15.08769_sp, 14.38594_sp, 13.68418_sp, 12.98243_sp, 12.28068_sp, &
           11.57893_sp, 10.87717_sp, 10.17542_sp, 9.473666_sp, 8.771913_sp,  8.07016_sp, 7.368407_sp, &
           6.666654_sp, 5.964901_sp, 5.263148_sp, 4.561395_sp, 3.859642_sp, 3.157889_sp, 2.456136_sp, &
           1.754383_sp,  1.05263_sp,0.3508765_sp,-0.350877_sp, -1.05263_sp,-1.754383_sp,-2.456136_sp, &
          -3.157889_sp,-3.859642_sp,-4.561395_sp,-5.263148_sp,-5.964901_sp,-6.666654_sp,              &
          -7.368407_sp, -8.07016_sp,-8.771913_sp,-9.473666_sp,-10.17542_sp,-10.87717_sp,              &
          -11.57893_sp,-12.28068_sp,-12.98243_sp,-13.68418_sp,-14.38594_sp,-15.08769_sp,              &
          -15.78944_sp,-16.49120_sp,-17.19295_sp, -17.8947_sp,-18.59645_sp,-19.29821_sp,              &
          -19.99996_sp,-20.70171_sp,-21.40347_sp,-22.10522_sp,-22.80697_sp,-23.50872_sp,              &
          -24.21048_sp,-24.91223_sp,-25.61398_sp,-26.31573_sp,-27.01749_sp,-27.71924_sp,              &
          -28.42099_sp,-29.12275_sp,-29.82450_sp,-30.52625_sp,-31.22800_sp,-31.92976_sp, -32.63151_sp,&
          -33.33326_sp,-34.03502_sp,-34.73677_sp,-35.43852_sp,-36.14027_sp,-36.84202_sp,              &
          -37.54378_sp,-38.24553_sp,-38.94728_sp,-39.64903_sp,-40.35078_sp,-41.05254_sp,              &
          -41.75429_sp,-42.45604_sp,-43.15779_sp,-43.85954_sp,-44.56129_sp,-45.26305_sp,              &
          -45.96480_sp,-46.66655_sp,-47.36830_sp,-48.07005_sp,-48.77180_sp,-49.47356_sp, -50.17531_sp,&
          -50.87706_sp,-51.57881_sp,-52.28056_sp,-52.98231_sp,-53.68406_sp,-54.38581_sp,              &
          -55.08756_sp,-55.78931_sp,-56.49106_sp,-57.19281_sp,-57.89456_sp,-58.59631_sp,              &
          -59.29806_sp,-59.99981_sp,-60.70156_sp,-61.40331_sp,-62.10506_sp,-62.80680_sp,              &
          -63.50855_sp,-64.21030_sp,-64.91205_sp,-65.61379_sp,-66.31554_sp,-67.01729_sp,              &
          -67.71903_sp,-68.42078_sp,-69.12252_sp,-69.82426_sp,-70.52601_sp,-71.22775_sp,              &
          -71.92949_sp,-72.63123_sp,-73.33297_sp,-74.03471_sp,-74.73644_sp,-75.43818_sp,              &
          -76.13991_sp,-76.84164_sp,-77.54337_sp,-78.24509_sp,-78.94681_sp,-79.64853_sp,              &
          -80.35023_sp,-81.05194_sp,-81.75363_sp,-82.45532_sp,-83.15699_sp,-83.85863_sp,              &
          -84.56026_sp,-85.26185_sp,-85.96337_sp, -86.6648_sp,-87.36607_sp,-88.06697_sp,              &
          -88.76695_sp,-89.46282_sp /)
            elseif(gg_order.eq.80)then
              ! this is the N80 grid (320x160) used by ERA 20C
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"Assuming this is ERA 20C; filling lon and lat"
              endif;enddo
              dx_met_const = 1.125_sp
              do i=1,nx_fullmet
                x_fullmet_sp(i) = 0.0_sp+(i-1)* dx_met_const
              enddo
              y_fullmet_sp(1:ny_fullmet) = &
         (/89.14152_sp, 88.02943_sp, 86.91077_sp, 85.79063_sp, 84.66992_sp, 83.54895_sp,              &
           82.42782_sp, 81.30659_sp, 80.18531_sp, 79.06398_sp, 77.94263_sp, 76.82124_sp, 75.69984_sp, &
           74.57843_sp, 73.45701_sp, 72.33558_sp, 71.21413_sp, 70.09269_sp, 68.97124_sp, 67.84978_sp, &
           66.72832_sp, 65.60686_sp, 64.48540_sp, 63.36393_sp, 62.24246_sp, 61.12099_sp, 59.99952_sp, &
           58.87804_sp, 57.75657_sp, 56.63509_sp, 55.51361_sp, 54.39214_sp, 53.27066_sp, 52.14917_sp, &
           51.02769_sp, 49.90621_sp, 48.78473_sp, 47.66325_sp, 46.54176_sp, 45.42028_sp, 44.29879_sp, &
           43.17731_sp, 42.05582_sp, 40.93434_sp, 39.81285_sp, 38.69136_sp, 37.56988_sp, 36.44839_sp, &
           35.32690_sp, 34.20542_sp, 33.08393_sp, 31.96244_sp, 30.84096_sp, 29.71947_sp, 28.59798_sp, &
           27.47649_sp, 26.35500_sp, 25.23351_sp, 24.11202_sp, 22.99054_sp, 21.86905_sp, 20.74756_sp, &
           19.62607_sp, 18.50458_sp, 17.38309_sp, 16.26160_sp, 15.14011_sp, 14.01862_sp, 12.89713_sp, &
           11.77564_sp, 10.65415_sp,  9.53266_sp,  8.41117_sp,  7.28968_sp, 6.168194_sp, 5.046704_sp, &
            3.92522_sp,  2.80373_sp,  1.68224_sp,  0.56074_sp, -0.56074_sp,-1.682235_sp, &
           -2.80373_sp, -3.92522_sp, -5.04670_sp, -6.16819_sp, -7.28968_sp,-8.411174_sp, &
           -9.53266_sp,-10.65415_sp,-11.77564_sp,-12.89713_sp,-14.01862_sp,-15.14011_sp, &
          -16.26160_sp,-17.38309_sp,-18.50458_sp,-19.62607_sp,-20.74756_sp,-21.86905_sp, &
          -22.99054_sp,-24.11202_sp,-25.23351_sp,-26.35500_sp,-27.47649_sp,-28.59798_sp, &
          -29.71947_sp,-30.84096_sp,-31.96244_sp,-33.08393_sp,-34.20542_sp,-35.32690_sp, &
          -36.44839_sp,-37.56988_sp,-38.69136_sp,-39.81285_sp,-40.93434_sp,-42.05582_sp, &
          -43.17731_sp,-44.29879_sp,-45.42028_sp,-46.54176_sp,-47.66325_sp,-48.78473_sp, &
          -49.90621_sp,-51.02769_sp,-52.14917_sp,-53.27066_sp,-54.39214_sp,-55.51361_sp, &
          -56.63509_sp,-57.75657_sp,-58.87804_sp,-59.99952_sp,-61.12099_sp,-62.24246_sp, &
          -63.36393_sp,-64.48540_sp,-65.60686_sp,-66.72832_sp,-67.84978_sp,-68.97124_sp, &
          -70.09269_sp,-71.21413_sp,-72.33558_sp,-73.45701_sp,-74.57843_sp,-75.69984_sp, &
          -76.82124_sp,-77.94263_sp,-79.06398_sp,-80.18531_sp,-81.30659_sp,-82.42782_sp, &
          -83.54895_sp,-84.66992_sp,-85.79063_sp,-86.91077_sp,-88.02943_sp,-89.14152_sp /)
            else
              ! We currently do not have a general gaussian grid calculator.
              ! Might incorporate this file at some point.
              ! https://github.com/NCAR/ncl/blob/develop/ni/src/lib/nfpfort/gaus.f
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"MR ERROR: Need to fix grid reading for gg grids."
              endif;enddo
              stop 1
            endif

            !Lat_start = y_start
            !Lon_start = x_start
            y_start = Lat_start 
            x_start = Lon_start 

            call codes_get(igribv(ir),'numberOfPoints',numberOfPoints,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
            allocate(values(numberOfPoints))
            call codes_get(igribv(ir),'values',values,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
            ReadGrid = .true.
            deallocate(values)
            !call codes_get(igribv(ir),'latitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
            !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get latitudeOfLastGridPointInDegrees ")
            !if(Lat_start.gt.dum_dp)then
            !  y_inverted = .true.
            !else
            !  y_inverted = .false.
            !endif
            !call codes_get(igribv(ir),&
            !      'longitudeOfLastGridPointInDegrees',dum_dp,nSTAT)
            !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get longitudeOfLastGridPointInDegrees ")
            !if(Lon_start.gt.dum_dp)then
            !  x_inverted = .true.
            !else
            !  x_inverted = .false.
            !endif
            x_fullmet_sp(0) = x_fullmet_sp(nx_fullmet)-360.0_sp
            x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(1)+360.0_sp

          elseif(index(dum_str,'polar_stereographic').ne.0)then
            IsLatLon_MetGrid = .false.
            Met_iprojflag     = 1
            call codes_get(igribv(ir),'orientationOfTheGridInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get orientationOfTheGridInDegrees ")
            Met_lam0 = dum_dp
            Met_phi0 = 90.0_8
            Met_phi1 = 90.0_8
            Met_k0   =  0.933_8
            Met_Re   =  6371.229_8

          elseif(index(dum_str,'albers').ne.0)then
            IsLatLon_MetGrid = .false.
            Met_iprojflag     = 2
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: Alber Equal Area not implemented"
            endif;enddo
            stop 1
          elseif(index(dum_str,'UTM').ne.0)then
            IsLatLon_MetGrid = .false.
            Met_iprojflag     = 3
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: UTM not implemented"
            endif;enddo
            stop 1
          elseif(index(dum_str,'lambert').ne.0)then
            IsLatLon_MetGrid = .false.
            Met_iprojflag     = 4
            call codes_get(igribv(ir),'LoVInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get LoVInDegrees ")
            Met_lam0 = dum_dp
            call codes_get(igribv(ir),'LaDInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get LaDInDegrees ")
            Met_phi0 = dum_dp
            call codes_get(igribv(ir),'Latin1InDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Latin1InDegrees ")
            Met_phi1 = dum_dp
            call codes_get(igribv(ir),'Latin2InDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Latin2InDegrees ")
            Met_phi2 = dum_dp
            Met_k0   =  0.933_8
            Met_Re   =  6371.229_8
          elseif(index(dum_str,'mercator').ne.0)then
            IsLatLon_MetGrid = .false.
            Met_iprojflag     = 5
            Met_lam0 = Lon_start
            call codes_get(igribv(ir),'LaDInDegrees',dum_dp,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get LaDInDegrees ")
            Met_phi0 = dum_dp
            Met_k0   =  0.933_8
            Met_Re   =  6371.229_8
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'MR ERROR: Cannot determine the projection from the GRIB file.'
            endif;enddo
            stop 1
          endif
          ! Override for the case of NARR
          if (MR_iwindformat.eq.3)Met_Re = 6367.470_8

          if(.not.IsLatLon_MetGrid)then
            call codes_get(igribv(ir),'DxInMetres',dum_int,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,0,"codes_get DxInMetres ")
            if(nSTAT.ne.0)then
              call codes_get(igribv(ir),'DiInMetres',dum_int,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get DiInMetres ")
            endif
            dx_met_const = real(dum_int,kind=4)/1000.0_sp
            call codes_get(igribv(ir),'DyInMetres',dum_int,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,0,"codes_get DyInMetres ")
            if(nSTAT.ne.0)then
              call codes_get(igribv(ir),'DjInMetres',dum_int,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get DjInMetres ")
            endif
            dy_met_const = real(dum_int,kind=4)/1000.0_sp

            call PJ_proj_for(Lon_start,Lat_start, Met_iprojflag, &
                     Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                     x_start,y_start)
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then   
              write(outlog(io),*)"Getting start coordinate for ",Lon_start,Lat_start
              write(outlog(io),*)" Projected coordinate = ",x_start,y_start
            endif;enddo
          endif

          if(.not.ReadGrid)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then   
              write(outlog(io),*)"Setting x grid starting at",x_start
              write(outlog(io),*)"with a spacing of ",dx_met_const
            endif;enddo
            do i = 0,nx_fullmet+1
              x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
            enddo
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then   
              write(outlog(io),*)"Setting y grid starting at",y_start
              write(outlog(io),*)"with a spacing of ",dy_met_const
            endif;enddo
            if(y_inverted)then
              do i = 1,ny_fullmet
                y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
              enddo
            else
              do i = 1,ny_fullmet
                y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
              enddo
            endif
            ReadGrid = .true.
          endif
          do i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          enddo
          do i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          enddo
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

          ! We need to check if this is a regular grid
          IsRegular_MetGrid = .true.
          do i = 1,nx_fullmet-1
            if(abs(MR_dx_met(i+1)-MR_dx_met(i)).gt.tol*MR_dx_met(i))then
              IsRegular_MetGrid = .false.
            endif
          enddo
          do i = 1,ny_fullmet-1
            if(abs(MR_dy_met(i+1)-MR_dy_met(i)).gt.tol*MR_dy_met(i))then
              IsRegular_MetGrid = .false.
            endif
          enddo

        endif ! count1.eq.1
        ! End of block reading the first record for grid info

        if(MR_GRIB_Version.eq.1)then
          call codes_get(igribv(ir),'indicatorOfTypeOfLevel',grb_typeSfc,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get indicatorOfTypeOfLevel ")
          ! for populating z-levels, we are only concerned with specific level types
          if(index(grb_typeSfc,'pl').ne.0.or. & ! Isobaric surface  (Pa)
             index(grb_typeSfc,'105').ne.0)then  ! Specified height level above ground  (m)
            call codes_get(igribv(ir),'indicatorOfParameter',grb_parameterNumber,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get indicatorOfParameter ")
            call codes_get(igribv(ir),'table2Version',grb_Table,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get table2Version ")

            ! Loop through all the variables and see if we have a match with this GRIB record
            do ivar = 1,MR_MAXVARS
              if (.not.Met_var_IsAvailable(ivar)) cycle
              iv_ParamN = Met_var_GRIB1_Param(ivar)
              iv_typeSfc   = Met_var_GRIB1_St(ivar)
              if(iv_ParamN.eq.grb_parameterNumber.and.       &
                   iv_typeSfc.eq.grb_typeSfc)then
                ! This is one we are tracking, log the level
                call codes_get(igribv(ir),'level',grb_level,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get level ")
                zcount(ivar) = zcount(ivar) + 1
                zlev_dum(ivar,zcount(ivar)) = grb_level * 100
              endif
            enddo
          endif
        elseif(MR_GRIB_Version.eq.2)then
          typeOfFirstFixedSurface            = -1
          grb_discipline                     = -1
          grb_parameterCategory              = -1
          grb_parameterNumber                = -1
          grb_scaledValueOfFirstFixedSurface = -1
          call codes_get(igribv(ir),'typeOfFirstFixedSurface', typeOfFirstFixedSurface,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get typeOfFirstFixedSurface ")
          ! for populating z-levels, we are only concerned with specific level types
          if(typeOfFirstFixedSurface.eq.100.or. & ! Isobaric surface  (Pa)
             typeOfFirstFixedSurface.eq.103.or. & ! Specified height level above ground  (m)
             typeOfFirstFixedSurface.eq.106)then  ! Depth below land surface  (m)
            call codes_get(igribv(ir),'discipline',        grb_discipline,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get discipline ")
            call codes_get(igribv(ir),'parameterCategory', grb_parameterCategory,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get parameterCategory ")
            call codes_get(igribv(ir),'parameterNumber', grb_parameterNumber,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get parameterNumber ")
            call codes_get(igribv(ir),'shortName',         grb_shortName,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get shortName ")
            call codes_get(igribv(ir),'cfNameECMF',        grb_longName,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get cfNameECMF ")
            call codes_get(igribv(ir),&
                 'scaledValueOfFirstFixedSurface',&
                 grb_scaledValueOfFirstFixedSurface,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get scaledValueOfFirstFixedSurface ")

            ! Loop through all the variables and see if we have a match with this GRIB record
            do ivar = 1,MR_MAXVARS
              if (.not.Met_var_IsAvailable(ivar)) cycle
              iv_discpl = Met_var_GRIB2_DPcPnSt(ivar,1)
              iv_paramC = Met_var_GRIB2_DPcPnSt(ivar,2)
              iv_paramN = Met_var_GRIB2_DPcPnSt(ivar,3)
              iv_typeSf = Met_var_GRIB2_DPcPnSt(ivar,4)
              if(iv_discpl.eq.grb_discipline.and.       &
                 iv_paramC.eq.grb_parameterCategory.and.&
                 iv_paramN.eq.grb_parameterNumber.and.  &
                 iv_typeSf.eq.typeOfFirstFixedSurface)then
                ! This is one we are tracking, log the level
                call codes_get(igribv(ir),&
                     'scaledValueOfFirstFixedSurface',&
                     grb_scaledValueOfFirstFixedSurface,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get scaledValueOfFirstFixedSurface ")

                ! Double-check that this isn't a duplicate record
                IsNewLevel = .true.
                do iz = 1,zcount(ivar)
                  if (zlev_dum(ivar,iz).eq.grb_scaledValueOfFirstFixedSurface) then
                    IsNewLevel = .false.
                    exit
                  endif
                enddo
                if (IsNewLevel) then
                  zcount(ivar) = zcount(ivar) + 1
                  zlev_dum(ivar,zcount(ivar)) = grb_scaledValueOfFirstFixedSurface
                endif
              endif
            enddo
          endif
        endif ! MR_GRIB_Version.eq.1 or .eq.2
      enddo

      do ir = 1,count1
        call codes_release(igribv(ir),nSTAT)
        if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
      enddo
      call codes_close_file(ifile,nSTAT)
      if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_close_file ")

      maxdimlen = maxval(zcount(:))
      nlev_coords_detected = 1
      Met_var_zdim_idx(:) = 0
      Met_var_zdim_idx(1) = 1
      do ivar = 2,MR_MAXVARS
        if (.not.Met_var_IsAvailable(ivar)) cycle
        FoundOldDim = .false.
        do iivar = 1,ivar-1
          if (zcount(iivar).eq.zcount(ivar))then  ! This check for a different coordinate is
                                                  ! solely on the size of the dimension.
            FoundOldDim = .true.
            Met_var_zdim_idx(ivar)  = Met_var_zdim_idx(iivar)
            exit
          endif
        enddo
        if(.not.FoundOldDim)then
          nlev_coords_detected = nlev_coords_detected + 1
          Met_var_zdim_idx(ivar)  = nlev_coords_detected
        endif
      enddo
      ! The V part of velocity is typically the second of a multi-component record and the above code
      ! does not catch it.  Assume V has the same characteristics as U and copy U
      ! V @ isobaric
      zcount(3) = zcount(2)
      zlev_dum(3,1:zcount(3)) = zlev_dum(2,1:zcount(2))
      Met_var_zdim_idx(3) = Met_var_zdim_idx(2)
      ! V @ 10 m
      zcount(12) = zcount(11)
      zlev_dum(12,1:zcount(12)) = zlev_dum(11,1:zcount(11))
      Met_var_zdim_idx(12) = Met_var_zdim_idx(11)

      ! We have all the level dimension names and dim_ids; now we need to get the sizes
      allocate(nlevs_fullmet(nlev_coords_detected))
      allocate(levs_code(nlev_coords_detected))
      allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
      do ivar = 1,MR_MAXVARS
        if (.not.Met_var_IsAvailable(ivar)) cycle
        ! Check if this variable has a z-dimension (pressure, height, depth, etc.)
        if(Met_var_zdim_idx(ivar).gt.0)then
          idx = Met_var_zdim_idx(ivar)
          nlevs_fullmet(idx) = zcount(ivar)
          ! Check that the pressure level is monotonically increasing
          Check=.true.
          do k=1,zcount(ivar)-1
            if(zlev_dum(ivar,k).gt.zlev_dum(ivar,k+1)) Check=.false.
          enddo
          if (Check)then
            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = &
              real(zlev_dum(ivar,1:nlevs_fullmet(idx)),kind=sp)
            z_inverted = .true.
            ! Pressure is expected to be from the bottom up so invert
            allocate(p_fullmet_sp(nlevs_fullmet(idx)))
            do k=1,zcount(ivar)
              p_fullmet_sp(k) = levs_fullmet_sp(idx,nlevs_fullmet(idx)+1-k)
            enddo
            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = p_fullmet_sp(1:nlevs_fullmet(idx))
            deallocate(p_fullmet_sp)
          else
            ! Need to first sort pressure variable here
            !   Using insertion sort on p321 of Numerical Recipes
            do k=2,nlevs_fullmet(idx)
              tmp1 = zlev_dum(ivar,k)
              do kk=k-1,1,-1
                if (zlev_dum(ivar,kk).ge.tmp1) goto 101
                zlev_dum(ivar,kk+1) = zlev_dum(ivar,kk)
              enddo
              kk=0
 101          zlev_dum(ivar,kk+1) = tmp1
            enddo

            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = &
              real(zlev_dum(ivar,1:nlevs_fullmet(idx)),kind=sp)
            z_inverted = .false.
          endif
        endif
      enddo

      ! Now log all pressure coordinates as one-to-one, truncated, or interrupted
      levs_code(1:nlev_coords_detected) = 0
      levs_code(1) = 1                       ! The first var checked (GPH) should have a one-to-one mapping
      ! Check how each of the pressure coordinates map onto the GPH grid
      if (nlev_coords_detected.gt.1)then
        ! Only bother if there are multiple pressure coordinates
        do idx = 2,nlev_coords_detected
          if (nlevs_fullmet(idx).gt.nlevs_fullmet(1))then
            ! This coordinate has more values than the GPH pressure coordinate
            levs_code(idx) = 4
          elseif (nlevs_fullmet(idx).lt.nlevs_fullmet(1))then
            ! It there are fewer levels, check if this is a truncated coordinate (code = 2)
            ! or one with missing levels that requires interpolation (code = 3)
            IsTruncatedDim = .true.
            do i=1,nlevs_fullmet(idx)
              if(abs(levs_fullmet_sp(idx,i)-levs_fullmet_sp(1,i)).gt.MR_EPS_SMALL)then
                IsTruncatedDim = .false.
                exit
              endif
            enddo
            if(IsTruncatedDim)then
              levs_code(idx) = 2
            else
              levs_code(idx) = 3
            endif
          else
            ! This coordinate has the same dimension as the GPH pressure coordinate.
            ! They are probably the same
            levs_code(idx) = 1
          endif
        enddo
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then            
        write(outlog(io),*)" Found these levels"
        write(outlog(io),*)&
          "  VaribleID    LevelIdx       dimID      length"
        do ivar = 1,MR_MAXVARS
          if (Met_var_IsAvailable(ivar))then 
            if(Met_var_zdim_idx(ivar).eq.0)then
              write(outlog(io),*)ivar,Met_var_zdim_idx(ivar),0,0,&
                                           trim(adjustl(Met_var_GRIB_names(ivar)))
            else
              write(outlog(io),*)ivar,Met_var_zdim_idx(ivar),0,&
                                           nlevs_fullmet(Met_var_zdim_idx(ivar)),&
                                           trim(adjustl(Met_var_GRIB_names(ivar)))
            endif
          endif
        enddo
      endif;enddo

      ! Now assign these levels to the working arrays
      ! Geopotential is the first variable checked, use this for np_fullmet
      nt_fullmet = 1
      np_fullmet = nlevs_fullmet(Met_var_zdim_idx(1))  ! Assign fullmet the length of H,U,V
      allocate(p_fullmet_sp(np_fullmet))
      idx = Met_var_zdim_idx(1)
      p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))
      MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp)

      allocate(z_approx(np_fullmet))
      do k=1,np_fullmet
        ! Calculate heights for US Std Atmos while pressures are still in mbars
        ! or hPa
        z_approx(k) = MR_Z_US_StdAtm(p_fullmet_sp(k))
      enddo

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then            
        write(outlog(io),*)"Dimension info:"
        write(outlog(io),*)"  record (time): ",nt_fullmet
        write(outlog(io),*)"  level  (z)   : ",np_fullmet
        write(outlog(io),*)"  y            : ",ny_fullmet
        write(outlog(io),*)"  x            : ",nx_fullmet
      endif;enddo

      !************************************************************************
      ! assign boundaries of mesoscale model
      if(x_inverted)then
          ! I know of no windfiles with x-coordinate reversed
        xLL_fullmet = x_fullmet_sp(nx_fullmet)
        xUR_fullmet = x_fullmet_sp(1)
      else
        xLL_fullmet = x_fullmet_sp(1)
        xUR_fullmet = x_fullmet_sp(nx_fullmet)
      endif

      if(y_inverted)then
          ! Most lon/lat grids have y reversed
        yLL_fullmet = y_fullmet_sp(ny_fullmet)
        yUR_fullmet = y_fullmet_sp(1)
      else
          ! Projected grids have y not reversed
        yLL_fullmet = y_fullmet_sp(1)
        yUR_fullmet = y_fullmet_sp(ny_fullmet)
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_Met_DimVars_GRIB

!##############################################################################

!##############################################################################
!
!     MR_Read_Met_Times_GRIB
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine opens each GRIB file and determines the time of each
!     time step of each file in the number of hours since MR_BaseYear.
!     In most cases, the length of the time variable (nt_fullmet) will be 
!     read directly from the file and overwritten (is was set in MR_Read_Met_DimVars_GRIB
!     above).
!
!     After this subroutine completes, the following variables will be set:
!       MR_windfile_starthour(MR_iwindfiles)
!       MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
!
!##############################################################################

      subroutine MR_Read_Met_Times_GRIB

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,verbosity_production,&
         MR_windfile_starthour,MR_iwindfiles,MR_windfile_stephour,&
         MR_BaseYear,MR_useLeap,MR_Comp_StartYear,MR_iwindformat,nt_fullmet,&
         Met_dim_IsAvailable,MR_windfiles,MR_windfiles_nt_fullmet,MR_GRIB_Version

      use eccodes

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      !integer, parameter :: dp        = 8 ! double precision
      integer, parameter :: MAXGRIBREC = 10000

      integer :: iw,iws
      integer :: itstart_year,itstart_month
      integer :: itstart_day
      real(kind=sp) :: filestart_hour

      integer :: itstart_hour,itstart_min,itstart_sec

      !real(kind=8)       :: HS_hours_since_baseyear
      character(len=130) :: dumstr
      integer            :: iwstep

      integer            :: dataDate
      integer            :: dataTime
      integer            :: forecastTime
      integer            :: ifile
      integer            :: count1
      integer,dimension(MAXGRIBREC) :: igribv
      integer(kind=4)    :: typeOfFirstFixedSurface
      integer            :: nSTAT

      integer :: io                           ! Index for output streams

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
        subroutine MR_GRIB_check_status(nSTAT, errcode, operation)
          integer, intent(in) :: nSTAT
          integer, intent(in) :: errcode
          character(len=*), intent(in) :: operation
        end subroutine MR_GRIB_check_status
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                MR_Read_Met_Times_GRIB             ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      if(.not.Met_dim_IsAvailable(1))then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Time dimension is required and not listed"
          write(errlog(io),*)"          in custom windfile specification file."
        endif;enddo
        stop 1
      endif

      allocate(MR_windfile_starthour(MR_iwindfiles))
      if(MR_iwindformat.eq.27)then
        ! GRIB1 reader not yet working!!
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: iwf=27 is a GRIB1 format."
          write(errlog(io),*)"       The GRIB1 reader is not yet working"
        endif;enddo
        stop 1

        ! Here the branch for when MR_iwindformat = 27
        ! First copy path read in to slot 2
        !if(MR_runAsForecast)then
        !  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        !    write(errlog(io),*)"MR ERROR: iwf=27 cannot be used for forecast runs."
        !    write(errlog(io),*)"          These are reanalysis files."
        !  endif;enddo
        !  stop 1
        !endif
        dumstr = MR_windfiles(1)
 110    format(a50,a1,i4,a1)
        write(MR_windfiles(1),110)trim(adjustl(dumstr)),'/', &
                                   MR_Comp_StartYear,'/'
        write(MR_windfiles(2),110)trim(adjustl(dumstr)),'/', &
                                   MR_Comp_StartYear+1,'/'
        MR_windfile_starthour(1) = real(HS_hours_since_baseyear( &
                                    MR_Comp_StartYear,1,1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
        MR_windfile_starthour(2) = real(HS_hours_since_baseyear( &
                                    MR_Comp_StartYear+1,1,1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
        if  ((mod(MR_Comp_StartYear,4).eq.0)     .and.                     &
             (mod(MR_Comp_StartYear,100).ne.0).or.(mod(MR_Comp_StartYear,400).eq.0))then
          nt_fullmet = 1464     ! Leap year
        else
          nt_fullmet = 1460     ! Not a leap year
        endif
        MR_windfiles_nt_fullmet(1)=nt_fullmet
        MR_windfiles_nt_fullmet(2)=nt_fullmet  ! Note: we don't care if the next
                                               !       year is a leap year since
                                               !       the simulation will never
                                               !       be long enough to matter.

        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

          ! the interval for iwf27 is 6 hours
        do iwstep = 1,nt_fullmet
          MR_windfile_stephour(:,iwstep) = (iwstep-1)*6.0_4
        enddo
      else
        ! For all other formats, try to read the first GRIB message and get
        ! dataDate, dataTime and forecastTime
        ! Loop through all the windfiles
        do iw = 1,MR_iwindfiles

          ! Each wind file needs a ref-time which in almost all cases is given
          ! in the 'units' attribute of the time variable
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)iw,trim(adjustl(MR_windfiles(iw)))
          endif;enddo

          if(iw.eq.1)then
            ! For now, assume one time step per file
            nt_fullmet = 1
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"  Assuming all NWP files have the same number of steps."
              write(outlog(io),*)"   For GRIB, assume one time step per file."
              write(outlog(io),*)"   Allocating time arrays for ",MR_iwindfiles,"file(s)"
              write(outlog(io),*)"                              ",nt_fullmet,"step(s) each"
            endif;enddo
            allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
          endif

          call codes_open_file(ifile,trim(adjustl(MR_windfiles(iw))),'R',nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_open_file ")

          count1=1
          call codes_new_from_file(ifile,igribv(count1),CODES_PRODUCT_GRIB,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")

          typeOfFirstFixedSurface = -1

          do while (count1.lt.100)  ! Read up to the first 100 grib records looking for a
                                    ! pressure level
            call codes_get(igribv(count1),'typeOfFirstFixedSurface', typeOfFirstFixedSurface,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS) exit
            ! for populating z-levels, we are only concerned with specific level types
            if(typeOfFirstFixedSurface.eq.100) exit ! Isobaric surface  (Pa)
            ! if we haven't found a pressure level, keep looking, some surface variables
            ! are not for the expected forcast hour
            count1 = count1+1
            call codes_new_from_file(ifile,igribv(count1),CODES_PRODUCT_GRIB,nSTAT)
          enddo


          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")
          if(iw.eq.1) then
            call codes_get(igribv(count1),'editionNumber',MR_GRIB_Version,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get editionNumber ")
          endif
          call codes_get(igribv(count1),'dataDate',dataDate,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get dataDate ")
          call codes_get(igribv(count1),'dataTime',dataTime,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get dataTime ")

          if(MR_GRIB_Version.eq.1)then
            ! The only grib1 files we deal with are reanalysis files with no FC time
            forecastTime = 0
          else
            call codes_get(igribv(count1),'forecastTime',forecastTime,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get forecastTime ")
          endif

          itstart_year  = int(dataDate/10000)
          itstart_month = int((dataDate-10000*itstart_year)/100)
          itstart_day   = mod(dataDate,100)
          itstart_hour  = int(dataTime/100)
          itstart_min   = mod(dataTime,100)
          itstart_sec   = 0

          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),2100)&
               "Ref time = ",itstart_year,itstart_month,itstart_day, &
               itstart_hour,itstart_min,itstart_sec
          endif;enddo

          call codes_release(igribv(count1),nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
          call codes_close_file(ifile,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_close_file ")

          filestart_hour = real(itstart_hour,kind=sp) + &
                           real(itstart_min,kind=sp)/60.0_sp      + &
                           real(itstart_sec,kind=sp)/3600.0_sp

          MR_windfiles_nt_fullmet(iw)=nt_fullmet
          MR_windfile_starthour(iw) =  real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                         itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
          MR_windfile_stephour(iw,1) = real(forecastTime,kind=4)

        enddo
      endif

2100  format(20x,a11,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)

      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
        write(outlog(io),*)"  File,  step,        Ref,     Offset,  HoursSince"
        do iw = 1,MR_iwindfiles
          do iws = 1,nt_fullmet
            write(outlog(io),800)iw,iws,real(MR_windfile_starthour(iw),kind=4),&
                             real(MR_windfile_stephour(iw,iws),kind=4),&
                             real(MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws),kind=4)
          enddo
        enddo
      endif;enddo
 800  format(i7,i7,3f12.2)

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_Met_Times_GRIB
!##############################################################################

!##############################################################################
!
!     MR_Read_MetP_Variable_GRIB
!
!     Called from Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_GRIB(ivar,istep)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,verbosity_debug1,&
         np_fullmet,levs_fullmet_sp,nlevs_fullmet,MR_MetStep_File,temp3d_sp,&
         MR_dum3d_metP,temp2d_sp,MR_dum2d_met,MR_EPS_SMALL,MR_geoH_metP_last,&
         fill_value_sp,ilhalf_fm_l,ilhalf_nx,irhalf_fm_l,irhalf_nx,istart,jstart,&
         MR_GRIB_Version,MR_iMetStep_Now,MR_iwindformat,nx_fullmet,ny_fullmet,&
         nx_submet,ny_submet,wrapgrid,Met_var_IsAvailable,Met_var_GRIB_names,&
         Met_var_GRIB1_Param,Met_var_GRIB1_Table,Met_var_GRIB1_St,Met_var_GRIB2_DPcPnSt,&
         y_inverted,MR_MetStep_findex,MR_MetStep_tindex,MR_iwind5_year,Met_var_zdim_idx,&
         MR_geoH_metP_next,p_fullmet_sp,Met_var_conversion_factor,&
           MR_Temp_US_StdAtm,&
           MR_Z_US_StdAtm,&
           MR_QC_3dvar

      use eccodes

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      !integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: ivar
      integer,intent(in) :: istep

      integer :: iw,iwstep
      integer :: np_met_loc
      character(len=80)  :: invar
      character(len=130) :: index_file
      character(len=130) :: grib_file_path
      character(len=130) :: grib_file

      real(kind=sp) :: del_H,del_P,dpdz

      integer :: i,j,k
      integer :: kk
      integer :: kkk,itmp
      integer :: ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      integer :: Dimension_of_Variable
      logical :: IsCategorical

      integer,dimension(np_fullmet) :: p_met_loc

      integer            :: ifile
      integer            :: igrib
      integer            :: idx
!      integer            :: iret
      integer            :: nSTAT
      integer            :: l,m,t
      integer            :: count1=0
      integer            :: rstrt, rend
      real(kind=8),dimension(:),allocatable     :: values
      real(kind=8),dimension(:,:),allocatable   :: slice
      integer(kind=4)  :: numberOfPoints
      integer(kind=4)  :: Ni
      integer(kind=4)  :: Nj
      integer(kind=4)  :: typeOfFirstFixedSurface
        ! Stores values of keys read from GRIB file
      !character(len=7) :: grb_marsParam
      character(len=4) :: grb_typeSfc
      integer(kind=4)  :: grb_discipline
      integer(kind=4)  :: grb_parameterCategory
      integer(kind=4)  :: grb_parameterNumber
      integer(kind=4)  :: grb_level
      integer(kind=4)  :: grb_scaledValueOfFirstFixedSurface
      !character(len=9)   :: sName
      !integer(kind=4)  :: forecastTime
      !integer(kind=4),dimension(:),allocatable :: marsParam_idx
      integer(kind=4),dimension(:),allocatable :: discipline_idx
      integer(kind=4),dimension(:),allocatable :: parameterCategory_idx
      integer(kind=4),dimension(:),allocatable :: parameterNumber_idx
      integer(kind=4),dimension(:),allocatable :: level_idx
      integer(kind=4),dimension(:),allocatable :: forecastTime_idx

      !integer(kind=4)  :: marsParamSize
      integer(kind=4)  :: disciplineSize
      integer(kind=4)  :: parameterCategorySize
      integer(kind=4)  :: parameterNumberSize
      integer(kind=4)  :: levelSize
      integer(kind=4)  :: forecastTimeSize
        ! temporary versions of those stored in Met_var_GRIB[...](ivar)
      integer(kind=4)  :: iv_discpl
      integer(kind=4)  :: iv_paramC
      integer(kind=4)  :: iv_paramN
      integer(kind=4)  :: iv_typeSf
      integer(kind=4)  :: iv_Table
      !character(len=7) :: iv_marsParam
      character(len=3) :: iv_typeSfc

      real(kind=sp) :: Z_top, T_top
      real(kind=sp),dimension(:,:,:),allocatable :: full_values

      logical :: Use_GRIB_Index = .false.
      integer :: fn_idx
      character(len=40)  :: fileposstr

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine MR_GRIB_check_status(nSTAT, errcode, operation)
          integer, intent(in) :: nSTAT
          integer, intent(in) :: errcode
          character(len=*), intent(in) :: operation
        end subroutine MR_GRIB_check_status
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                MR_Read_MetP_Variable_GRIB         ----------"
        write(outlog(io),*)ivar,istep
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      if(.not.Met_var_IsAvailable(ivar))then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  Variable not available for this windfile"
          write(errlog(io),*)"             ivar = ",ivar
          write(errlog(io),*)"            vname = ",Met_var_GRIB_names(ivar)
          write(errlog(io),*)"             iwf  = ",MR_iwindformat
        endif;enddo
        stop 1
      endif

      if(ivar.eq.3 .or. &   ! V_isobaric
         ivar.eq.12)then    ! V_height_above_ground
        ! The v components might be not reachable with the indexing if they are the second
        ! component of a multi-variable message
        Use_GRIB_Index = .false.
!      else
!        Use_GRIB_Index = .true.
      endif

      if(MR_GRIB_Version.eq.1)then
        ! GRIB1 uses an index based on Param
        iv_paramN  = Met_var_GRIB1_Param(ivar)
        iv_Table   = Met_var_GRIB1_Table(ivar)
        iv_typeSfc = Met_var_GRIB1_St(ivar)
      elseif(MR_GRIB_Version.eq.2)then
        ! Get the variable discipline, Parameter Category, Parameter Number, and
        ! level type for this variable
        iv_discpl = Met_var_GRIB2_DPcPnSt(ivar,1)
        iv_paramC = Met_var_GRIB2_DPcPnSt(ivar,2)
        iv_paramN = Met_var_GRIB2_DPcPnSt(ivar,3)
        iv_typeSf = Met_var_GRIB2_DPcPnSt(ivar,4)
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  GRIB type not determined"
        endif;enddo
        stop 1
      endif

      iw     = MR_MetStep_findex(istep)
      iwstep = MR_MetStep_tindex(istep)

      if(Met_var_GRIB_names(ivar).eq."")then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Variable ",ivar,&
                    " not available for MR_iwindformat = ",&
                    MR_iwindformat
        endif;enddo
        stop 1
      endif

      ! Get the dimension of the variable requested (either 2 or 3-D)
      if(ivar.eq.1 ) Dimension_of_Variable = 3 ! Geopotential Height
      if(ivar.eq.2 ) Dimension_of_Variable = 3 ! Vx
      if(ivar.eq.3 ) Dimension_of_Variable = 3 ! Vy
      if(ivar.eq.4 ) Dimension_of_Variable = 3 ! Vz
      if(ivar.eq.5 ) Dimension_of_Variable = 3 ! Temperature
      if(ivar.eq.6 ) Dimension_of_Variable = 3 ! Pressure (only for WRF or other eta-level files)
      if(ivar.eq.7 ) Dimension_of_Variable = 3 ! PVV

      if(ivar.eq.10) Dimension_of_Variable = 2 ! Planetary Boundary Layer Height
      if(ivar.eq.11) Dimension_of_Variable = 2 ! U @ 10m
      if(ivar.eq.12) Dimension_of_Variable = 2 ! V @ 10m
      if(ivar.eq.13) Dimension_of_Variable = 2 ! Friction velocity
      if(ivar.eq.14) Dimension_of_Variable = 2 ! Displacement Height
      if(ivar.eq.15) Dimension_of_Variable = 2 ! Snow cover
      if(ivar.eq.16) Dimension_of_Variable = 2 ! Soil moisture
      if(ivar.eq.17) Dimension_of_Variable = 2 ! Surface roughness
      if(ivar.eq.18) Dimension_of_Variable = 2 ! Wind_speed_gust_surface

      if(ivar.eq.20) Dimension_of_Variable = 2 ! pressure at lower cloud base
      if(ivar.eq.21) Dimension_of_Variable = 2 ! pressure at lower cloud top
      if(ivar.eq.22) Dimension_of_Variable = 2 ! temperature at lower cloud top
      if(ivar.eq.23) Dimension_of_Variable = 2 ! Total Cloud cover
      if(ivar.eq.24) Dimension_of_Variable = 2 ! Cloud cover (low)
      if(ivar.eq.25) Dimension_of_Variable = 2 ! Cloud cover (convective)

      if(ivar.eq.30) Dimension_of_Variable = 3 ! Rel. Hum
      if(ivar.eq.31) Dimension_of_Variable = 3 ! QV (specific humidity)
      if(ivar.eq.32) Dimension_of_Variable = 3 ! QL (liquid)
      if(ivar.eq.33) Dimension_of_Variable = 3 ! QI (ice)

      if(ivar.eq.40) Dimension_of_Variable = 2 ! Categorical rain
      if(ivar.eq.41) Dimension_of_Variable = 2 ! Categorical snow
      if(ivar.eq.42) Dimension_of_Variable = 2 ! Categorical frozen rain
      if(ivar.eq.43) Dimension_of_Variable = 2 ! Categorical ice
      if(ivar.eq.44) Dimension_of_Variable = 2 ! Precipitation rate large-scale (liquid)
      if(ivar.eq.45) Dimension_of_Variable = 2 ! Precipitation rate convective (liquid)
      if(ivar.eq.46) Dimension_of_Variable = 3 ! Precipitation rate large-scale (ice)
      if(ivar.eq.47) Dimension_of_Variable = 3 ! Precipitation rate convective (ice)

      if(ivar.eq.40.or.&
         ivar.eq.41.or.&
         ivar.eq.42.or.&
         ivar.eq.43)then
          ! Categorical variables are integers and need special interpolation
        IsCategorical = .true.
      else
          ! The default is to read floating point values
        IsCategorical = .false.
      endif

      if(MR_iwindformat.eq.27)then
        ! Get correct GRIB1 file
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: iwf27 not working for GRIB1"
        endif;enddo
        stop 1
        if(ivar.eq.1)then
          write(index_file,125)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_HGT_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.2)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_UGRD_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.3)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_VGRD_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.4)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_VVEL_pres.nc"
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"+++++++++++++++++++++++++++++++++++++++++++"
            write(errlog(io),*)"  NEED TO FIX THIS Vz"
          endif;enddo
          stop 1
          np_met_loc = np_fullmet
        elseif(ivar.eq.5)then
          write(index_file,125)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_TMP_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.10)then
          write(index_file,128)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_HPBL_sfc.nc"
        elseif(ivar.eq.22)then
          write(index_file,130)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_TMP_low-cldtop.nc"
        elseif(ivar.eq.23)then
          write(index_file,131)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_TCDC_low-cldlay.nc"
        elseif(ivar.eq.30)then
          write(index_file,127)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_RH_pres.nc"
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"+++++++++++++++++++++++++++++++++++++++++++"
            write(errlog(io),*)"  NEED TO FIX THIS : RH"
          endif;enddo
          stop 1
          np_met_loc = np_fullmet
        elseif(ivar.eq.44)then
          write(index_file,129)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_PRATE_sfc.nc"
        elseif(ivar.eq.45)then
          write(index_file,129)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_CPRAT_sfc.nc"
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then 
            write(errlog(io),*)"MR ERROR : Requested variable not available."
          endif;enddo
          stop 1
        endif
        index_file = trim(adjustl(index_file))

 125      format(a50,a13,i4,a12)
 126      format(a50,a13,i4,a13)
 127      format(a50,a13,i4,a11)
 128      format(a50,a15,i4,a12)
 129      format(a50,a15,i4,a13)
 130      format(a50,a15,i4,a18)
 131      format(a50,a15,i4,a19)
      else  ! all other cases besides iwf27
        ! Set up pressure level index that we will search for
        p_met_loc = 0
        idx = Met_var_zdim_idx(ivar)
        if(ivar.eq.4)then      ! Vertical_velocity_pressure_isobaric
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        elseif(ivar.eq.10)then ! Planetary_Boundary_Layer_Height_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.11)then ! u-component_of_wind_height_above_ground
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 10
        elseif(ivar.eq.12)then ! v-component_of_wind_height_above_ground
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 10
        elseif(ivar.eq.13)then ! Frictional_Velocity_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.15)then ! Snow_depth_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.16)then ! Volumetric_Soil_Moisture_Content_depth_below_surface_layer
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.17)then ! Surface_roughness_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.18)then ! Wind_speed_gust_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.20)then ! Pressure_cloud_base
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.21)then ! Pressure_cloud_topw
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.23)then ! Total_cloud_cover_entire_atmosphere
           ! Something is wrong reading this
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.30)then ! Relative_humidity_isobaric
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        elseif(ivar.eq.40.or.ivar.eq.41.or.ivar.eq.42.or.ivar.eq.43)then ! categorical precip
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.44)then ! Precipitation_rate_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        else
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = &
              int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        endif
        allocate(full_values(nx_fullmet,ny_fullmet,np_met_loc))
          ! Files are listed directly, not through directories (as in MR_iwindformat=25,27)
        grib_file_path  = trim(adjustl(MR_MetStep_File(istep)))
        fn_idx = index(MR_MetStep_File(istep), '/' , BACK=.true.)
        grib_file  = MR_MetStep_File(istep)(fn_idx+1:)
        index_file = trim(adjustl(MR_MetStep_File(istep))) // ".index"
      endif
      invar = Met_var_GRIB_names(ivar)

      ! Load data variables for just the subgrid defined above
      if (wrapgrid) then
        ict        = 2
          ! index on the sub-met
        ileft(1)   = 1;         ileft(2)   = ilhalf_nx+1
        iright(1)  = ilhalf_nx; iright(2)  = nx_submet
          ! indes on the full-met
        iistart(1) = ilhalf_fm_l; iistart(2) = irhalf_fm_l
        iicount(1) = ilhalf_nx  ; iicount(2) = irhalf_nx
      else
        ict        = 1
        ileft(1)   = 1
        iright(1)  = nx_submet
        iistart(1) = istart
        iicount(1) = nx_submet
      endif

      !----------------------------------------------------
      ! This is the part where we actually read the files
      ! 4 branches total
      !  Grib1 with index file
      !  Grib1 without index file
      !  Grib2 with index file
      !  Grib2 without index file
      if(MR_GRIB_Version.eq.1)then
        if(Use_GRIB_Index)then
          write(fileposstr,'(a9,i4,a9,i4,a10,i4)')"  step = ",istep,&
                         ", file = ",iw,&
                         ", slice = ",iwstep
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Reading ",trim(adjustl(invar)),&
               " from file : ",trim(adjustl(index_file)),fileposstr
          endif;enddo

          call codes_index_read(idx,index_file)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_read ")
          !call codes_grib_multi_support_on()
          !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_grib_multi_support_on")

            ! get the number of distinct values of all the keys in the index
          call codes_index_get_size(idx,'indicatorOfParameter',parameterNumberSize)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")
          call codes_index_get_size(idx,'level',levelSize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")

            ! allocate the array to contain the list of distinct values
          allocate(parameterNumber_idx(parameterNumberSize))
          allocate(level_idx(levelSize))

            ! get the list of distinct key values from the index
          call codes_index_get(idx,'indicatorOfParameter',parameterNumber_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")
          call codes_index_get(idx,'level',level_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")

          ! Start marching through the index file and look for the match with the 
          ! keys
          count1=0
          do l=1,parameterNumberSize
            call codes_index_select(idx,'indicatorOfParameter',parameterNumber_idx(l),nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")

            do i=1,levelSize
              call codes_index_select(idx,'level',level_idx(i),nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
              call codes_new_from_index(idx,igrib,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_index ")

              do while (nSTAT /= GRIB_END_OF_INDEX)
                count1=count1+1
                call codes_get(igrib,'indicatorOfTypeOfLevel',grb_typeSfc,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get indicatorOfTypeOfLevel ")

                if( parameterNumber_idx(l)           .eq. iv_paramN .and. &
                    grb_typeSfc(1:3)        .eq. iv_typeSfc) then
                  call codes_get(igrib,'numberOfPoints',numberOfPoints,nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
                  call codes_get(igrib,'Ni',Ni,nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Ni ")
                  call codes_get(igrib,'Nj',Nj,nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Nj ")
                  if(nx_fullmet.ne.Ni)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then 
                     write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                      write(errlog(io),*)"nx_fullmet = ",nx_fullmet
                      write(errlog(io),*)"Ni         = ",Ni
                    endif;enddo
                    stop 1
                  endif
                  if(ny_fullmet.ne.Nj)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                      write(errlog(io),*)"ny_fullmet = ",ny_fullmet
                      write(errlog(io),*)"Nj         = ",Nj
                    endif;enddo
                    stop 1
                  endif
                  allocate(values(numberOfPoints))
                  allocate(slice(Ni,Nj))
                  call codes_get(igrib,'values',values,nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
                  do m = 1,Nj
                    rstrt = (m-1)*Ni + 1
                    rend  = m*Ni
                    slice(1:Ni,m) = values(rstrt:rend)
                  enddo
                  deallocate(values)

                  ! There is no guarantee that GRIB levels are in order so...
                  ! Now loop through the pressure values for this variable and put
                  ! this slice at the correct level.
                  call codes_get(igrib,'level',grb_level,nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get level ")
                  do kk = 1,np_met_loc
                    if(p_met_loc(kk).eq.grb_level)then
                      full_values(:,:,kk) = real(slice(:,:),kind=sp)
                      exit
                    endif
                  enddo
                  deallocate(slice)
                endif

                call codes_release(igrib,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
                call codes_new_from_index(idx,igrib,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_index ")
              enddo ! while
              call codes_release(igrib,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
            enddo ! loop on level
          enddo ! loop on marsParam

          !call grib_index_release(idx)

        else ! Non-index Grib1 case

          ! We don't have/(can't make) the index file so scan all messages of the
          ! GRIB1 file
          write(fileposstr,'(a9,i4,a9,i4,a10,i4)')"  step = ",istep,&
                         ", file = ",iw,&
                         ", slice = ",iwstep
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Reading ",trim(adjustl(invar)),&
                  " from file : ",&
                  trim(adjustl(grib_file_path)),fileposstr
          endif;enddo

          ifile=5
          call codes_open_file(ifile,grib_file_path,'R',nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_open_file ")

          !     turn on support for multi fields messages */
          call codes_grib_multi_support_on()
          !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_grib_multi_support_on ")

          ! Loop on all the messages in a file.
          call codes_new_from_file(ifile,igrib,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")
          count1=0

          do while (nSTAT.eq.CODES_SUCCESS.and.igrib.gt.0)
            count1=count1+1

            call codes_get(igrib,'indicatorOfParameter',grb_parameterNumber,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get indicatorOfParameter ")
            call codes_get(igrib,'indicatorOfTypeOfLevel',grb_typeSfc,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get indicatorOfTypeOfLevel ")

            if ( grb_parameterNumber  .eq. iv_paramN .and. &
                 grb_typeSfc(1:3)     .eq. iv_typeSfc) then
              call codes_get(igrib,'numberOfPoints',numberOfPoints,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
              call codes_get(igrib,'Ni',Ni,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Ni ")
              call codes_get(igrib,'Nj',Nj,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Nj ")

              if(nx_fullmet.ne.Ni)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"nx_fullmet = ",nx_fullmet
                  write(errlog(io),*)"Ni         = ",Ni
                endif;enddo
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*) "MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"ny_fullmet = ",ny_fullmet
                  write(errlog(io),*)"Nj         = ",Nj
                endif;enddo
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
              call codes_get(igrib,'values',values,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
              do m = 1,Nj
                rstrt = (m-1)*Ni + 1
                rend  = m*Ni
                slice(1:Ni,m) = values(rstrt:rend)
              enddo
              deallocate(values)

               ! There is no guarantee that GRIB levels are in order so...
               ! Now loop through the pressure values for this variable and put
               ! this slice at the correct level.
               call codes_get(igrib,'level',grb_level,nSTAT)
               if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get level ")
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.grb_level)then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
            endif
            call codes_release(igrib,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
            call codes_new_from_file(ifile,igrib,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_file ")

          enddo
          call codes_close_file(ifile,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_close_file ")

        endif

      elseif(MR_GRIB_Version.eq.2)then
        if(Use_GRIB_Index)then
          write(fileposstr,'(a9,i4,a9,i4,a10,i4)')"  step = ",istep,&
                         ", file = ",iw,&
                         ", slice = ",iwstep
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)&
               "Reading ",trim(adjustl(invar))," from file : ",&
               trim(adjustl(index_file)),fileposstr
          endif;enddo

          call codes_index_read(idx,index_file,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_read ")
          call codes_grib_multi_support_on()
          !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_grib_multi_support_on")

            ! get the number of distinct values of all the keys in the index
          call codes_index_get_size(idx,'discipline',disciplineSize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")
          call codes_index_get_size(idx,'parameterCategory',parameterCategorySize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")
          call codes_index_get_size(idx,'parameterNumber',parameterNumberSize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")
          call codes_index_get_size(idx,'scaledValueOfFirstFixedSurface',levelSize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")
          call codes_index_get_size(idx,'forecastTime',forecastTimeSize,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get_size ")

            ! allocate the array to contain the list of distinct values
          allocate(discipline_idx(disciplineSize))
          allocate(parameterCategory_idx(parameterCategorySize))
          allocate(parameterNumber_idx(parameterNumberSize))
          allocate(level_idx(levelSize))
          allocate(forecastTime_idx(forecastTimeSize))
          
            ! get the list of distinct key values from the index
          call codes_index_get(idx,'discipline',discipline_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")
          call codes_index_get(idx,'parameterCategory',parameterCategory_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")
          call codes_index_get(idx,'parameterNumber',parameterNumber_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")
          call codes_index_get(idx,'scaledValueOfFirstFixedSurface',level_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")
          call codes_index_get(idx,'forecastTime',forecastTime_idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_get ")

          ! Start marching through the index file and look for the match with the 
          ! keys
          count1=0
          do l=1,disciplineSize
            call codes_index_select(idx,'discipline',discipline_idx(l),nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
            do j=1,parameterCategorySize
              call codes_index_select(idx,'parameterCategory',parameterCategory_idx(j),nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
              do k=1,parameterNumberSize
                call codes_index_select(idx,'parameterNumber',parameterNumber_idx(k),nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
                do i=1,levelSize
                  call codes_index_select(idx,'level',level_idx(i),nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
                  call codes_index_select(idx,&
                         'scaledValueOfFirstFixedSurface',level_idx(i),nSTAT)
                  if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")

                  do t=1,forecastTimeSize
                    call codes_index_select(idx,'forecastTime',forecastTime_idx(t),nSTAT)
                    if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_select ")
                    call codes_new_from_index(idx,igrib,nSTAT)
                    if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_new_from_index ")

                    do while (nSTAT /= GRIB_END_OF_INDEX)
                      count1=count1+1
    
            call codes_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get typeOfFirstFixedSurface ")

            if ( discipline_idx(l)       .eq. iv_discpl .and. &
                 parameterCategory_idx(j).eq. iv_paramC .and. &
                 parameterNumber_idx(k)  .eq. iv_paramN .and. &
                 typeOfFirstFixedSurface .eq. iv_typeSf) then
              call codes_get(igrib,'numberOfPoints',numberOfPoints,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
              call codes_get(igrib,'Ni',Ni,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Ni ")
              call codes_get(igrib,'Nj',Nj,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Nj ")

              if(nx_fullmet.ne.Ni)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"nx_fullmet = ",nx_fullmet
                  write(errlog(io),*)"Ni         = ",Ni
                endif;enddo
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"ny_fullmet = ",ny_fullmet
                  write(errlog(io),*)"Nj         = ",Nj
                endif;enddo
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
                call codes_get(igrib,'values',values,nSTAT)
                if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
                do m = 1,Nj
                  rstrt = (m-1)*Ni + 1
                  rend  = m*Ni
                  slice(1:Ni,m) = values(rstrt:rend)
                enddo
                deallocate(values)
        
               ! There is no guarantee that GRIB levels are in order so...
               ! Now loop through the pressure values for this variable and put this
               ! slice at the correct level.
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.level_idx(i))then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
             endif

                      call grib_release(igrib)
                      call grib_new_from_index(idx,igrib,nSTAT)
                    enddo
                    call codes_release(igrib,nSTAT)
                    if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
                  enddo ! loop on forecastTime
                enddo ! loop on level
              enddo ! loop on parameterNumber
            enddo ! loop on parameterCategory
          enddo ! loop on discipline

          call codes_index_release(idx,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_index_release ")
        else
          ! We don't have/(can't make) the index file so scan all messages of the
          ! GRIB2 file
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)&
               istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
               trim(adjustl(grib_file_path))
          endif;enddo
          ifile=5
          call codes_open_file(ifile,grib_file_path,'R',nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_open_file ")

          !     turn on support for multi fields messages */
          call codes_grib_multi_support_on()
          !if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"")

          ! Loop on all the messages in a file.
          call codes_new_from_file(ifile,igrib,CODES_PRODUCT_GRIB,nSTAT)
          count1=0
          do while (nSTAT.eq.CODES_SUCCESS.and.igrib.gt.0)
            count1=count1+1
            call codes_get(igrib,'discipline',              grb_discipline,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get discipline ")
            call codes_get(igrib,'parameterCategory',       grb_parameterCategory,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get parameterCategory ")
            call codes_get(igrib,'parameterNumber', grb_parameterNumber,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get parameterNumber ")
            call codes_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get typeOfFirstFixedSurface ")
            call codes_get(igrib,'scaledValueOfFirstFixedSurface',grb_level,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get scaledValueOfFirstFixedSurface ")
            call codes_get(igrib,'scaledValueOfFirstFixedSurface',grb_scaledValueOfFirstFixedSurface,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get scaledValueOfFirstFixedSurface ")

            if ( grb_discipline              .eq. iv_discpl .and. &
                 grb_parameterCategory       .eq. iv_paramC .and. &
                 grb_parameterNumber         .eq. iv_paramN .and. &
                 typeOfFirstFixedSurface     .eq. iv_typeSf) then
              call codes_get(igrib,'numberOfPoints',numberOfPoints,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get numberOfPoints ")
              call codes_get(igrib,'Ni',Ni,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Ni ")
              call codes_get(igrib,'Nj',Nj,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get Nj ")

              if(nx_fullmet.ne.Ni)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"nx_fullmet = ",nx_fullmet
                  write(errlog(io),*)"Ni         = ",Ni
                endif;enddo
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*) "MR ERROR:  Grid is not the expected size"
                  write(errlog(io),*)"ny_fullmet = ",ny_fullmet
                  write(errlog(io),*)"Nj         = ",Nj
                endif;enddo
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
              call codes_get(igrib,'values',values,nSTAT)
              if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_get values ")
              do m = 1,Nj
                rstrt = (m-1)*Ni + 1
                rend  = m*Ni
                slice(1:Ni,m) = values(rstrt:rend)
              enddo
              deallocate(values)
  
               ! There is no guarantee that GRIB levels are in order so...
               ! Now loop through the pressure values for this variable and put
               ! this slice at the correct level.
               if(ivar.eq.16)grb_level = grb_scaledValueOfFirstFixedSurface
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.grb_level)then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
            endif
            call codes_release(igrib,nSTAT)
            if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_release ")
            call codes_new_from_file(ifile,igrib,CODES_PRODUCT_GRIB,nSTAT)
          enddo
          call codes_close_file(ifile,nSTAT)
          if(nSTAT.ne.CODES_SUCCESS)call MR_GRIB_check_status(nSTAT,1,"codes_close_file ")
        endif ! Use_GRIB_index
      endif ! MR_GRIB_Version eq 1 or 2

      if(Dimension_of_Variable.eq.3)then
        MR_dum3d_metP = 0.0_sp
        allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))

        do i=1,ict        !read subgrid at current time step
            ! for any other 3d variable (non-WRF, non-NCEP)
          temp3d_sp(ileft(i):iright(i)              ,1:ny_submet            ,1:np_met_loc,1) = &
          full_values(iistart(i):iistart(i)+iicount(i)-1,jstart:jstart+ny_submet-1,1:np_met_loc)
        enddo

          do j=1,ny_submet
            itmp = ny_submet-j+1
            !reverse the j indices (since they increment from N to S)
            if(y_inverted)then
              MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,itmp,1:np_met_loc,1)
            else
              MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,j,1:np_met_loc,1)
            endif
          enddo

        deallocate(temp3d_sp)

      elseif(Dimension_of_Variable.eq.2)then
!        if(IsCategorical)then
!        else
          allocate(temp2d_sp(nx_submet,ny_submet,1))
          if(ivar.eq.11.or.ivar.eq.12)then
              ! Surface winds usually have a z coordinate as well
            allocate(temp3d_sp(nx_submet,ny_submet,1,1))
          endif
  
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then

            else
              ! 2d variables for iwf .ne. 25
!              if(ivar.eq.11.or.ivar.eq.12)then
!                ! Surface velocities do have a z dimension
!                nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
!                         start = (/iistart(i),jstart,1,iwstep/),       &
!                         count = (/iicount(i),ny_submet,1,1/))
!                if(nSTAT.ne.0)then
!                   do io=1,MR_nio;if(VB(io).le.verbosity_error)then
!                     write(errlog(io),*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   endif;enddo
!                   stop 1
!                endif
!                do j=1,ny_submet
!                  itmp = ny_submet-j+1
!                  if(y_inverted)then
!                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,itmp,1,1)
!                  else
!                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,j,1,1)
!                  endif
!                enddo
!              else
!                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(ileft(i):iright(i),:,:), &
!                         start = (/iistart(i),jstart,iwstep/),       &
!                         count = (/iicount(i),ny_submet,1/))
!                if(nSTAT.ne.0)then
!                   do io=1,MR_nio;if(VB(io).le.verbosity_error)then
!                     write(errlog(io),*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   endif;enddo
!                   stop 1
!                endif
          temp2d_sp(ileft(i):iright(i)              ,1:ny_submet,1) = &
        full_values(iistart(i):iistart(i)+iicount(i)-1,jstart:jstart+ny_submet-1,1)

                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,itmp,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,j,1)
                  endif
                enddo
!              endif
            endif
          enddo
          deallocate(temp2d_sp)
          if(ivar.eq.11.or.ivar.eq.12) deallocate(temp3d_sp)
!        endif ! IsCategorical
      endif ! Dimension_of_Variable.eq.2

      if(ivar.eq.1)then
        ! If this is filling HGT, then we need to do a special QC check
        do i=1,nx_submet
          do j=1,ny_submet
            do k=1,np_met_loc
              if(abs(MR_dum3d_metP(i,j,k)-fill_value_sp).lt.MR_EPS_SMALL.or.&
                     MR_dum3d_metP(i,j,k).lt.0.0_sp)then  ! also flag value as to be reset if it
                                                          ! maps below sea level
                 ! linearly interpolate in z
                 ! find the first non NaN above k
                 do kk = k+1,np_met_loc,1
                   if(abs(MR_dum3d_metP(i,j,kk)-fill_value_sp).gt.MR_EPS_SMALL.and.&
                          MR_dum3d_metP(i,j,kk).ge.0.0_sp)exit
                 enddo
                 if(kk.eq.np_met_loc+1)then
                   kk=np_met_loc
                   MR_dum3d_metP(i,j,kk) = 0.0_sp
                 endif
                 ! find the first non NaN below k if k!=1
                 do kkk = max(k-1,1),1,-1
                   if(abs(MR_dum3d_metP(i,j,kkk)-fill_value_sp).gt.MR_EPS_SMALL.and.&
                          MR_dum3d_metP(i,j,kkk).ge.0.0_sp)exit
                 enddo
                 if(kkk.eq.0)then
                   kkk=1
                   MR_dum3d_metP(i,j,kkk) = 0.0_sp
                 endif
                 MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,kkk) + &
                       (MR_dum3d_metP(i,j,kk)-MR_dum3d_metP(i,j,kkk)) * &
                       real(k-kkk,kind=sp)/real(kk-kkk,kind=sp)
              endif
            enddo
          enddo
        enddo
        ! convert m to km
        MR_dum3d_metP = MR_dum3d_metP / 1000.0_sp
      elseif(Dimension_of_Variable.eq.3)then
        ! Do QC checking of all other 3d variables
        if(ivar.eq.2.or.ivar.eq.3.or.ivar.eq.4)then
          ! taper winds (vx,vy,vz) to zero at ground surface
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=0.0_sp)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=0.0_sp)
          endif
        elseif(ivar.eq.5)then
          ! set ground and top-level conditions for temperature
          Z_top = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/real(100.0,kind=sp))
          T_top = MR_Temp_US_StdAtm(Z_top)
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          endif
        else
          ! For other variables, use the top and bottom non-fill values
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp)
          endif
        endif
      endif

      if(ivar.eq.4)then
          ! For pressure vertical velocity, convert from Pa s to m/s by dividing
          ! by pressure gradient
        idx = Met_var_zdim_idx(ivar)
        do k=1,np_met_loc
          do i=1,nx_submet
            do j=1,ny_submet
              if(k.eq.1)then
                ! Use one-sided gradients for bottom
                !del_P = p_fullmet_Vz_sp(2)-p_fullmet_Vz_sp(1)
                del_p = levs_fullmet_sp(idx,2) - levs_fullmet_sp(idx,1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,2) - MR_geoH_metP_last(i,j,1)
                else
                  del_H = MR_geoH_metP_next(i,j,2) - MR_geoH_metP_next(i,j,1)
                endif
              elseif(k.eq.np_met_loc)then
                ! Use one-sided gradients for top
                !del_P = p_fullmet_Vz_sp(np_met_loc) - &
                !         p_fullmet_Vz_sp(np_met_loc-1)
                del_p = levs_fullmet_sp(idx,nlevs_fullmet(idx)) - &
                        levs_fullmet_sp(idx,nlevs_fullmet(idx)-1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,np_met_loc) - &
                           MR_geoH_metP_last(i,j,np_met_loc-1)
                else
                  del_H = MR_geoH_metP_next(i,j,np_met_loc) - &
                           MR_geoH_metP_next(i,j,np_met_loc-1)
                endif
              else
                ! otherwise, two-sided calculation
                !del_P = p_fullmet_Vz_sp(k+1)-p_fullmet_Vz_sp(k-1)
                del_p = levs_fullmet_sp(idx,k+1) - &
                        levs_fullmet_sp(idx,k)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,k+1) - MR_geoH_metP_last(i,j,k-1)
                else
                  del_H = MR_geoH_metP_next(i,j,k+1) - MR_geoH_metP_next(i,j,k-1)
                endif
              endif
              del_h = del_H * 1000.0_sp ! convert to m
              if(abs(del_H).gt.MR_EPS_SMALL)then
                dpdz  = del_P/del_H
              else
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)'MR ERROR: failed to calculate dpdz'
                  write(errlog(io),*)i,j,k,del_P,del_H
                  write(errlog(io),*)MR_geoH_metP_last(i,j,:)
                endif;enddo
                stop 1
              endif
              MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,k) / dpdz
            enddo
          enddo
        enddo
      endif
      MR_dum3d_metP = MR_dum3d_metP * Met_var_conversion_factor(ivar)

      do io=1,MR_nio;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Read_MetP_Variable_GRIB

!##############################################################################
!
!     MR_GRIB_check_status
!
!     iSTAT   = error code returned from eccodes call
!     errcode = user-supplied return value on stopping of code
!     operation = string descriptor of function call causing error
!
!     Error-checking routine for eccodes function calls.
!     Modeled after a subroutine posted at:
!     https://climate-cms.org/2018/10/12/create-netcdf.html
!
!##############################################################################

      subroutine MR_GRIB_check_status(nSTAT, errcode, operation)

      use MetReader,       only : &
         MR_nio,VB,errlog,verbosity_error

      use eccodes

      implicit none

      integer, intent(in) :: nSTAT
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: operation

      character(len=50) :: err_message
      character(len=12) :: severity

      integer :: io                           ! Index for output streams

      if (errcode.eq.0)then
        severity = "MR WARNING: "
       else
        severity = "MR ERROR:   "
      endif

      if (nSTAT == CODES_SUCCESS) return
      call codes_get_error_string(nSTAT,err_message)
      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io) ,*)severity,errcode,operation,' ',trim(adjustl(err_message))
      endif;enddo

      ! If user-supplied error code is 0, then consider this a warning,
      ! otherwise do a hard stop
      if (errcode.ne.0) stop

      end subroutine MR_GRIB_check_status

