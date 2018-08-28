      subroutine MR_Set_Gen_Index_GRIB(grib_file)

      use grib_api

      implicit none

      integer            :: MR_global_info   = 6
      integer            :: MR_global_error  = 0

      character(len=130)  :: grib_file

      character(len=130)  :: index_file

      integer            :: idx
      integer            :: iret
      integer            :: i,j,k,l,t
      integer            :: count1=0
        ! Used for keys
      character(len=7),dimension(:),allocatable :: marsParam_idx
      integer(kind=4) ,dimension(:),allocatable :: level_idx
      integer(kind=4)  :: marsParamSize
      integer(kind=4)  :: levelSize

      integer(kind=4) ,dimension(:),allocatable :: discipline_idx
      integer(kind=4) ,dimension(:),allocatable :: parameterCategory_idx
      integer(kind=4) ,dimension(:),allocatable :: parameterNumber_idx
      integer(kind=4) ,dimension(:),allocatable :: forecastTime_idx
      integer(kind=4)  :: disciplineSize
      integer(kind=4)  :: parameterCategorySize
      integer(kind=4)  :: parameterNumberSize
      integer(kind=4)  :: forecastTimeSize

      ! Message identifier.
      integer          :: igrib
      integer          :: gribver
      logical          :: Got_Version = .false.

      ! First, open grib file and determine if we have a grib1 or grib2 file
      call grib_open_file(idx,adjustl(trim(grib_file)),'r')
      call grib_new_from_file(idx,igrib,iret)
      call grib_get(igrib,'editionNumber',gribver)
      call grib_close_file(idx)
      if(gribver.eq.1.or.gribver.eq.2)then
        write(*,*)"Grib version detected = ",gribver
        Got_Version = .true.
      else
        write(MR_global_error,*)'MR ERROR : Could not detect grib version.'
        write(MR_global_error,*)'           version must be 1 or 2.'
        stop 0
      endif

      index_file = adjustl(trim(grib_file)) // ".index"
      write(MR_global_info,*)"Generating index file: ",index_file

      if(gribver.eq.1)then
        ! For grib1 files, we use the MARS parameter numbers
        !!!  grib_ls -p marsParam,date,time,level
          ! create an index from a grib file using some keys
        call grib_index_create(idx,adjustl(trim(grib_file)),&
              'marsParam,level')
        call grib_multi_support_on()
          ! get the number of distinct values of all the keys in the index
        call grib_index_get_size(idx,'marsParam',marsParamSize)
        call grib_index_get_size(idx,'level',levelSize)

          ! allocate the arry to contain the list of distinct values
        allocate(marsParam_idx(marsParamSize))
        allocate(level_idx(levelSize))

          ! get the list of distinct key values from the index
        call grib_index_get(idx,'marsParam',marsParam_idx)
        call grib_index_get(idx,'level',level_idx)

        count1=0
        do l=1,marsParamSize
          call grib_index_select(idx,'marsParam',marsParam_idx(l))

          do i=1,levelSize
            call grib_index_select(idx,'level',level_idx(i))

            call grib_new_from_index(idx,igrib, iret)
            do while (iret /= GRIB_END_OF_INDEX)
              count1=count1+1
              call grib_release(igrib)
              call grib_new_from_index(idx,igrib, iret)
            enddo
            call grib_release(igrib)

          enddo ! loop on level
        enddo ! loop on marsParam

        call grib_index_write(idx,adjustl(trim(index_file)))

        call grib_index_release(idx)

      elseif(gribver.eq.2)then
        ! For grib2 files, we use the message discipline, parameter category, parameter number,
        ! and layer type (pressure level, surface, tropopause, etc.) as keys for builing the index

        !----------------------------------------------------------------------
        ! discipline
        ! 0 Meteorological products
        ! 1 Hydrological products
        ! 2 Land surface products
        !10 Oceanographic products
        
        ! parameterCategory
        ! 0 Temperature
        ! 1 Moisture
        ! 2 Momentum
        ! 3 Mass
        ! 4 Short-wave Radiation
        ! 5 Long-wave Radiation
        ! 6 Cloud
        ! 7 Thermodynamic Stability indices
        
        ! parameterNumber
        !  This is a sub-category of the discipline and parameterCategory
        
        ! typeOfFirstFixedSurface
        !  1 Ground or water surface
        !  2 Cloud base level
        !  3 Level of cloud tops
        !100 Isobaric surface  (Pa)
        !103 Specified height level above ground  (m)
        !106 Depth below land surface  (m)
        !200 Unknown code table entry
        ! http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
        !    shortName 
        !            discipline 
        !              parameterCategory 
        !                  parameterNumber
        ! 1  gh      0 3   5 100   HGT   Geopotential_height_isobaric
        ! 2  u       0 2   2 100  UGRD   u-component_of_wind_isobaric
        ! 3  v       0 2   3 100  VGRD   v-component_of_wind_isobaric
        ! 4  w       0 2   8 100  VVEL   Vertical_velocity_pressure_isobaric
        ! 5  t       0 0   0 100   TMP   Temperature_isobaric
        ! 10 hpbl    0 3 196   1  HPBL   Planetary_Boundary_Layer_Height_surface
        ! 11 u       0 2   2 103  UGRD   u-component_of_wind_height_above_ground
        ! 12 v       0 2   3 103  VGRD   v-component_of_wind_height_above_ground
        ! 13 fricv   0 2 197   1 FRICV   Frictional_Velocity_surface
        ! 15 sd      0 1  11   1  SNOD   Snow_depth_surface
        ! 16 soilw   2 0 192 106 SOILW
        ! Volumetric_Soil_Moisture_Content_depth_below_surface_layer
        ! 17 sr      2 0   1   1  SFCR   Surface_roughness_surface
        ! 18 gust    0 2  22   1  GUST   Wind_speed_gust_surface
        ! 20 pres    0 3   0   2  PRES   Pressure_cloud_base
        ! 21 pres    0 3   0   3  PRES   Pressure_cloud_tops
        ! 23 tcc     0 6   1 200  TCDC   Total_cloud_cover_entire_atmosphere
        ! 30 r       0 1   1 100    RH   Relative_humidity_isobaric
        ! 31 q       0 1   0 100  SPFH   Specific_humidity_isobaric
        ! 32 clwmr   0 1  22 100 CLWMR   Cloud_mixing_ratio_isobaric
        ! 33 snmr    0 1  25 100  SNMR   Snow_mixing_ratio_isobaric
        ! 40 crain   0 1 192   1 CRAIN   Categorical_Rain_surface
        ! 41 csnow   0 1 195   1 CSNOW   Categorical_Snow_surface
        ! 42 cfrzr   0 1 193   1 CFRZR   Categorical_Freezing_Rain_surface
        ! 43 cicep   0 1 194   1 CICEP   Categorical_Ice_Pellets_surface
        ! 44 prate   0 1   7   1 PRATE   Precipitation_rate_surface
        ! 45                   1 CPRAT   Convective_Precipitation_Rate_surface
        ! 45 tp      0 1   8   APCP   Total_precipitation_surface_0_Hour_Accumulation
        ! 46                  ACPCP
        ! Convective_precipitation_surface_0_Hour_Accumulation
        ! 47 ncpcp   0 1   9  NCPCP
        ! Large-scale_precipitation_non-convective_surface_0_Hour_Accumulation
        !----------------------------------------------------------------------

          ! create an index from a grib file using some keys
        call grib_index_create(idx,adjustl(trim(grib_file)),&
              'discipline,parameterCategory,parameterNumber,scaledValueOfFirstFixedSurface,forecastTime')

        call grib_multi_support_on()
      
          ! get the number of distinct values of all the keys in the index
        call grib_index_get_size(idx,'discipline',disciplineSize)
        call grib_index_get_size(idx,'parameterCategory',parameterCategorySize)
        call grib_index_get_size(idx,'parameterNumber',parameterNumberSize)
        call grib_index_get_size(idx,'scaledValueOfFirstFixedSurface',levelSize)
        call grib_index_get_size(idx,'forecastTime',forecastTimeSize)
      
          ! allocate the arry to contain the list of distinct values
        allocate(discipline_idx(disciplineSize))
        allocate(parameterCategory_idx(parameterCategorySize))
        allocate(parameterNumber_idx(parameterNumberSize))
        allocate(level_idx(levelSize))
        allocate(forecastTime_idx(forecastTimeSize))
      
          ! get the list of distinct key values from the index
        call grib_index_get(idx,'discipline',discipline_idx)
        call grib_index_get(idx,'parameterCategory',parameterCategory_idx)
        call grib_index_get(idx,'parameterNumber',parameterNumber_idx)
        call grib_index_get(idx,'scaledValueOfFirstFixedSurface',level_idx)
        call grib_index_get(idx,'forecastTime',forecastTime_idx)
      
        count1=0
        do l=1,disciplineSize
          call grib_index_select(idx,'discipline',discipline_idx(l))
      
          do j=1,parameterCategorySize
            call grib_index_select(idx,'parameterCategory',parameterCategory_idx(j))
      
            do k=1,parameterNumberSize
              call grib_index_select(idx,'parameterNumber',parameterNumber_idx(k))
      
              do i=1,levelSize
                call grib_index_select(idx,'scaledValueOfFirstFixedSurface',level_idx(i))
      
                do t=1,forecastTimeSize
                  call grib_index_select(idx,'forecastTime',forecastTime_idx(t))
      
      
                call grib_new_from_index(idx,igrib, iret)
                do while (iret /= GRIB_END_OF_INDEX)
                   count1=count1+1
                   call grib_release(igrib)
                   call grib_new_from_index(idx,igrib, iret)
                enddo
                call grib_release(igrib)
      
                enddo ! loop on forecastTime
              enddo ! loop on level
            enddo ! loop on parameterNumber
          enddo ! loop on parameterCategory
        enddo ! loop on discipline
      
        call grib_index_write(idx,adjustl(trim(index_file)))
      
        call grib_index_release(idx)
      endif  ! grib version number

      end subroutine MR_Set_Gen_Index_GRIB

