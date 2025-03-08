!##############################################################################
!##############################################################################
!
! MetTraj_[FB]
!
! This is a stand-alone program that uses the MetReader interface to calculate
! trajectories from a point (at various altitudes) using any of the NWP products
! available via MetReader.  This is designed to run with a series of
! command-line arguments which will default to GFS forecast packages (if the
! requested time is within 14 days of the execution time) or the NCEP 50-year
! reanalysis data.  This program assumes that the wind files are stored in
! /data/WindFiles/ as expected by the download scripts in autorun_scripts/.
! The minimal set of command-line arguments are:  lon lat YYYY MM DD h.h
!
!     MetTraj_[F,B] -169.9468 52.8217 2022 8 29 5.5
!
! Optionally, the simulation time (in hours) can be provided (default is 24).
! Additionally, the trajectory levels can be provided, first with the number of
! levels, followed by the level altitudes in km.  To specify 6 hours
! integrations at 3 levels (5, 10, and 15 km), use:
!
!     MetTraj_[F,B] -169.9468 52.8217 2022 8 29 5.5 6.0 3 5.0 10.0 15.0
!
! If levels are not provided, then 6 levels are assumed at heights of 1.52,
! 3.05, 6.10, 9.14, 12.19, 15.24 km; corresponding to 5000, 10000, 20000, 30000,
! 40000, and 50000 ft.
!
! The full set of options, including streamline vs streakline, projected grids,
! etc. are available using a control file.
!
! This program was designed for short-term trajectory plots, so it includes some
! simplifying assumptions that will break down at longer integration times
! Caveats:
!  (1) Currently, the GPH values at the probe point are used to map from
!      pressure levels to altitude throughout the simulation. This obviously
!      falls apart at long integration times.
!  (2) No vertical advection is calculated
!  (3) First-order forward Euler integration is used
!  (4) The time-step is fixed to 1 minute
!  (5) The sub-grid of the windfile is assumed from the total integration time.
!      Anything over 24.0 hours loads the full global grid.
!
!
!##############################################################################

      program MetTraj

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,verbosity_production,&
         dx_met_const,dy_met_const,IsLatLon_CompGrid,&
           MR_Initialize_Met_Grids,&
           MR_Reset_Memory,&
           MR_FileIO_Error_Handler

      implicit none

      ! These are the variables that must be set in the input file or command line
      real(kind=8)        :: inlon
      real(kind=8)        :: inlat
      integer             :: inyear,inmonth,inday
      real(kind=8)        :: inhour
      real(kind=8)        :: Simtime_in_hours
      real(kind=8)        :: Met_hours_needed
      integer             :: StreamFlag
      integer             :: OutStepInc_Minutes
      integer             :: ntraj
      real(kind=8), dimension(9) :: OutputLevels
      integer             :: iw
      integer             :: iwf
      integer             :: igrid
      integer             :: idf
      integer             :: iwfiles
      integer             :: autoflag
      integer             :: FC_freq = 12 ! This is the frequency the GFS packages are downloaded
      integer             :: GFS_Archive_Days = 14
      integer             :: GFS_FC_TotHours = 198  ! This is the max anticipated

      integer             :: TrajFlag  ! .ge. 0 for forward, .lt.0 for backward

      integer             :: nxmax,nymax,nzmax
      real(kind=8)        :: dx,dy
      real(kind=8)        :: xwidth,ywidth
      real(kind=4),dimension(:)    ,allocatable :: lon_grid
      real(kind=4),dimension(:)    ,allocatable :: lat_grid
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      logical             :: IsPeriodic
      integer             :: i

      real(kind=8)        :: starty

      logical             :: IsGlobal

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine Read_ComdLine_InpFile(inlon,inlat, &
                      inyear,inmonth,inday,inhour,Simtime_in_hours,&
                      StreamFlag,OutStepInc_Minutes,ntraj,OutputLevels,&
                      iw,iwf,igrid,idf,iwfiles,&
                      autoflag,FC_freq,GFS_Archive_Days)
          integer,parameter   :: dp        = 8 ! double precision
          real(kind=dp)              ,intent(out) :: inlon
          real(kind=dp)              ,intent(out) :: inlat
          integer                    ,intent(out) :: inyear,inmonth,inday
          real(kind=dp)              ,intent(out) :: inhour
          real(kind=dp)              ,intent(out) :: Simtime_in_hours
          integer                    ,intent(out) :: StreamFlag
          integer                    ,intent(out) :: OutStepInc_Minutes
          integer                    ,intent(out) :: ntraj
          real(kind=dp), dimension(9),intent(out) :: OutputLevels
          integer                    ,intent(out) :: iw
          integer                    ,intent(out) :: iwf
          integer                    ,intent(out) :: igrid
          integer                    ,intent(out) :: idf
          integer                    ,intent(out) :: iwfiles
          integer                    ,intent(out) :: autoflag
          integer                    ,intent(out) :: FC_freq
          integer                    ,intent(out) :: GFS_Archive_Days
        end subroutine Read_ComdLine_InpFile
        subroutine GetWindFile(inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,&
                                iw,iwf,igrid,idf,iwfiles,&
                                autoflag,FC_freq,GFS_Archive_Days,GFS_FC_TotHours)
          integer,parameter   :: dp        = 8 ! double precision
          integer         ,intent(in) :: inyear
          integer         ,intent(in) :: inmonth
          integer         ,intent(in) :: inday
          real(kind=dp)   ,intent(in) :: inhour
          real(kind=dp)   ,intent(in) :: Simtime_in_hours
          integer         ,intent(in) :: TrajFlag
          integer         ,intent(inout) :: iw
          integer         ,intent(inout) :: iwf
          integer         ,intent(inout) :: igrid
          integer         ,intent(inout) :: idf
          integer         ,intent(inout) :: iwfiles
          integer         ,intent(in) :: autoflag
          integer         ,intent(in) :: FC_freq
          integer         ,intent(in) :: GFS_Archive_Days
          integer         ,intent(in) :: GFS_FC_TotHours
        end subroutine GetWindFile
        subroutine Integrate_ConstH_Traj(IsGlobal,inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,output_interv)
          integer,parameter   :: dp        = 8 ! double precision
          logical      , intent(in)      :: IsGlobal
          real(kind=dp), intent(in)      :: inlon
          real(kind=dp), intent(in)      :: inlat
          integer      , intent(in)      :: inyear
          integer      , intent(in)      :: inmonth
          integer      , intent(in)      :: inday
          real(kind=dp), intent(in)      :: inhour
          real(kind=dp), intent(in)      :: Simtime_in_hours
          integer      , intent(in)      :: TrajFlag
          integer      , intent(in)      :: ntraj
          integer      , intent(in)      :: output_interv
        end subroutine Integrate_ConstH_Traj
      END INTERFACE

      call Read_ComdLine_InpFile(inlon,inlat, &
                      inyear,inmonth,inday,inhour,Simtime_in_hours,&
                      StreamFlag,OutStepInc_Minutes,ntraj,OutputLevels,&
                      iw,iwf,igrid,idf,iwfiles,&
                      autoflag,FC_freq,GFS_Archive_Days)

      ! Now set up the computational grid
      TrajFlag = 0
#ifdef FORWARD
      TrajFlag = 1
#endif
#ifdef BACKWARD
      TrajFlag = -1
#endif

      if(TrajFlag.eq.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Forward/Backward not specified."
          write(errlog(io),*)"          Recompile with preprocessor flags"
          write(errlog(io),*)"           FORWARD for forward trajectories"
          write(errlog(io),*)"           BACKWARD for backward trajectories"
        endif;enddo
        stop 1
      elseif(TrajFlag.lt.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)"Calculating Backward trajectories"
        endif;enddo
      elseif(TrajFlag.gt.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)"Calculating Forward trajectories"
        endif;enddo
      endif
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Calling GetWindFile"
      endif;enddo
      if(StreamFlag.eq.0)then
        Met_hours_needed = 0.0_8
      else
        Met_hours_needed = Simtime_in_hours
      endif
      call GetWindFile(inyear,inmonth,inday,inhour, &
                          Met_hours_needed,TrajFlag,&
                          iw,iwf,igrid,idf,iwfiles,&
                          autoflag,FC_freq,GFS_Archive_Days,GFS_FC_TotHours)

      ! Now set up computational grid.
      nzmax = ntraj
      IsGlobal = .false.
      ! Define grid padding based on the integration time
      if(Simtime_in_hours.le.8.0_8)then
        if(IsLatLon_CompGrid)then
          ! +-15 in lon ; +-10 in lat
          xwidth = 30.0_4
          ywidth = 20.0_4
        else
          ! +-300km in x ; +-300 in y
          xwidth = 300.0_4
          ywidth = 300.0_4
        endif
      elseif(Simtime_in_hours.le.16.0_8)then
        if(IsLatLon_CompGrid)then
          ! +-25 in lon ; +-15 in lat
          xwidth = 50.0_4
          ywidth = 30.0_4
        else
          ! +-600km in x ; +-600 in y
          xwidth = 600.0_4
          ywidth = 600.0_4
        endif
      elseif(Simtime_in_hours.le.24.0_8)then
        if(IsLatLon_CompGrid)then
          ! +-35 in lon ; +-20 in lat
          xwidth = 70.0_4
          ywidth = 40.0_4
        else
          ! +-1000km in x ; +-1000 in y
          xwidth = 1000.0_4
          ywidth = 1000.0_4
        endif
      else
        if(IsLatLon_CompGrid)then
          ! Full globe
          xwidth = 360.0_4
          ywidth = 180.0_4
          IsGlobal = .true.
        else
          ! +-1500km in x ; +-1500 in y
          xwidth = 1500.0_4
          ywidth = 1500.0_4
        endif
      endif

      if(iw.eq.1)then
        ! For the ASCII case, the met grid is not a 2-d grid, but maybe
        ! scattered points.  We need to set up comp grid independent of
        ! the met grid.
        dx = 0.5_4
        dy = 0.5_4
      else
        if(dx_met_const.gt.0.0.and.dy_met_const.gt.0.0)then
          dx = dx_met_const
          dy = dy_met_const
        else
          dx = 0.5_4
          dy = 0.5_4
        endif
      endif
      nxmax = ceiling(xwidth/dx) + 1
      nymax = ceiling(ywidth/dy) + 1

      allocate(lon_grid(0:nxmax+1))
      allocate(lat_grid(0:nymax+1))
      allocate(z_cc(nzmax))
      if(IsGlobal)then
        do i=0,nxmax+1
          lon_grid(i) = real((i-1) * dx,kind=4)
        enddo
        IsPeriodic = .true.
      else
        do i=0,nxmax+1
          lon_grid(i) = real(inlon - 0.5_4*(nxmax-1) * dx + (i-1) * dx,kind=4)
        enddo
        if(lon_grid(1).lt.0.0_4)then
          lon_grid(:) = lon_grid(:) + 360.0_4
        endif
      endif
      ! We can specify the longitude grid with no problems, but the latitude
      ! grid might be padded up across the pole.  We need to set the cap at
      ! an extreme point (89 degrees N or S)
      if(IsGlobal)then
        do i=0,nymax+1
          lat_grid(i) = real(-90.0_4 + (i-1) * dy,kind=4)
        enddo
      else
        if((inlat + (nymax-1)/2 * dy).gt.89.0_8)then
          ! Start from 89.0 N and count down nymax
          starty = 89.0_8 - (nymax-1) * dy
          do i=0,nymax+1
            lat_grid(i) = real(starty + (i-1) * dy,kind=4)
          enddo
        elseif((inlat - (nymax-1)/2 * dy).lt.-89.0_8)then
          ! Start from 89.0 N and count down nymax
          starty = -89.0_8
          do i=0,nymax+1
            lat_grid(i) = real(starty + (i-1) * dy,kind=4)
          enddo
        else
          ! lat grid doesn't involve poles; center grid over inlat
          do i=0,nymax+1
            lat_grid(i) = real(inlat - 0.5_8*(nymax-1) * dy + (i-1) * dy,kind=4)
          enddo
        endif
        IsPeriodic = .false.
      endif
      do i = 1,ntraj
        z_cc(i) = real(OutputLevels(i),4)
      enddo

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Setting up wind grids"
      endif;enddo
      ! Again, since we only need the metH grid, lon_grid and lat_grid are dummy
      ! arrays
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now integrating from start point"
      endif;enddo
      call Integrate_ConstH_Traj(IsGlobal,inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,OutStepInc_Minutes)

      call MR_Reset_Memory
      if(allocated(        z_cc)) deallocate(z_cc)
      if(allocated(    lon_grid)) deallocate(lon_grid)
      if(allocated(    lat_grid)) deallocate(lat_grid)

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Program ended normally."
      endif;enddo

      end program MetTraj

!##############################################################################
!
!     Subroutines
!
!##############################################################################

!##############################################################################
!  Read_ComdLine_InpFile
!
!  This subroutine will parse the command-line options and set up the run
!  specifications.  This will either be determined through a command file if
!  only one command-line argument is given, or through a set of 6 or more
!  command-line arguments.  All aspects of the run will be set, except the
!  wind-file names, which are set in a separate subroutine (GetWindFile).
!
! Example control file:
! -169.9468 52.8217                     ! lon lat
! 1875 6 20 5.5                         ! YYYY MM DD HH.H
! 24.0                                  ! simtime
! 1                                     ! streamflag (0 for streak lines (static
! windfield), 1 for streamlines)
! 60                                    ! output time step (minutes)
! 6                                     ! ntraj (<10)
! 1.524 3.048 6.096 9.144 12.192 15.240 ! level values in km
! 1 4 -107.0 50.0 50.0 50.0 6367.470    ! Output projection
! 5 27 1027 2                           ! iwind iwindformat igrid iformat
! 0 12 14                               ! autoflag (0 for auto, 1 for specified)
! FC_freq GFS_Archive_Days
! 1                                     ! number of windfiles
! NOAA20CRv3
!
! Example command-line argument:
!
!     MetTraj_[F,B] -169.9468 52.8217 2022 8 29 5.5 6.0 3 5.0 10.0 15.0
!
!##############################################################################

      subroutine Read_ComdLine_InpFile(inlon,inlat, &
                      inyear,inmonth,inday,inhour,Simtime_in_hours,&
                      StreamFlag,OutStepInc_Minutes,ntraj,OutputLevels,&
                      iw,iwf,igrid,idf,iwfiles,&
                      autoflag,FC_freq,GFS_Archive_Days)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,&
         MR_BaseYear,MR_useLeap,MR_useCompH,Comp_lam1,Comp_lam2,&
         Comp_iprojflag,Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re,&
         IsLatLon_CompGrid,&
           MR_Set_CompProjection,&
           MR_FileIO_Error_Handler

      use projection,      only : &
         PJ_ilatlonflag,PJ_iprojflag,PJ_k0,PJ_lam0,PJ_lam1,PJ_lam2,PJ_phi0,PJ_phi1,PJ_phi2,PJ_Re,&
           PJ_Set_Proj_Params

      implicit none

      integer,parameter :: fid_ctrlfile = 10

      ! These are the variables that must be set in the input file or command line
      real(kind=8)              ,intent(out) :: inlon
      real(kind=8)              ,intent(out) :: inlat
      integer                   ,intent(out) :: inyear,inmonth,inday
      real(kind=8)              ,intent(out) :: inhour
      real(kind=8)              ,intent(out) :: Simtime_in_hours
      integer                   ,intent(out) :: StreamFlag
      integer                   ,intent(out) :: OutStepInc_Minutes
      integer                   ,intent(out) :: ntraj
      real(kind=8), dimension(9),intent(out) :: OutputLevels
      integer                   ,intent(out) :: iw
      integer                   ,intent(out) :: iwf
      integer                   ,intent(out) :: igrid
      integer                   ,intent(out) :: idf
      integer                   ,intent(out) :: iwfiles
      integer                   ,intent(out) :: autoflag
      integer                   ,intent(out) :: FC_freq
      integer                   ,intent(out) :: GFS_Archive_Days

      !logical             :: IsLatLon
      integer             :: nargs
      integer             :: iostatus
      integer             :: inlen
      character(len=120)  :: iomessage
      character(len=130)  :: arg

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.
      integer :: i
      real(kind=4) :: tmp_4

      character(len=100):: infile
      logical           :: IsThere
      character(len=80) :: linebuffer080
      character(len=80) :: Comp_projection_line

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine Print_Usage
        end subroutine Print_Usage
      END INTERFACE

      ! Test read command line arguments
      nargs = command_argument_count()
      if (nargs.gt.1.and.nargs.lt.6) then
        ! We need either one command-line argument (input file name) or at least
        ! 6 parameter to define the run.
        ! Write usage to stdout and exit
        call Print_Usage
      elseif(nargs.ge.6)then
        ! Here, everything is set from the command-line with lots of assumed
        ! values.  Only GFS and NCEP 50-year are used here.
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Reading command-line"
        endif;enddo
        !  And here is what we assume:
        StreamFlag = 1  ! This means we are doing streamlines, NOT streaklines
        Simtime_in_hours = 24.0_8   ! Length of time to integrate (can be changed on command-line)
        OutStepInc_Minutes = 60     ! Minutes between output points
        ntraj              = 0      ! Number of trajectories (can be changed on command-line)
        ! OutputLevels : this is allocated once ntraj is locked in
        IsLatLon_CompGrid  = .true. ! Assume LonLat output coordinates
        autoflag           = 1      ! This command-line branch necesarily means auto windfile selection
                                    !  with all the hard-wired paths to GFS and NCEP
        FC_freq            = 12     ! Number of hours between GFS package downloads
        GFS_Archive_Days   = 14     ! Number of days GFS data are archived on local machine

        iw      = 0 ! These are all set for autoruns in GetWindFile
        iwf     = 0 
        igrid   = 0 
        idf     = 0
        iwfiles = 0 

        ! Make user MetReader is using the same calendar
        MR_BaseYear = BaseYear  ! This defaults to 1900 for autoruns, but can be something else
                                ! for command-file runs
        MR_useLeap  = useLeap
        MR_useCompH = .false.

        ! Minimum required is lon, lat, YYYY, MM, DD, hours
        call get_command_argument(1, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read first command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inlon
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inlon
        if(inlon.lt.-360.0_8.or.&
           inlon.gt.360.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"MR ERROR: longitude not in range -180->360"
            write(errlog(io),*)" inlon = ",inlon
          endif;enddo
          stop 1
        endif
        if(inlon.lt.0.0_8.or.inlon.gt.360.0_8)inlon=mod(inlon+360.0_8,360.0_8)

        call get_command_argument(2, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read second command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inlat
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inlat
        if(inlat.lt.-90.0_4.or.&
           inlat.gt.90.0_4)then
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"MR ERROR: latitude not in range -90->90"
            write(errlog(io),*)" inlat = ",inlat
          endif;enddo
          stop 1
        endif

        call get_command_argument(3, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read third command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inyear
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)

        call get_command_argument(4, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read fourth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inmonth
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inmonth
        if(inmonth.lt.1.or.&
           inmonth.gt.12)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: month must be in range 1-12"
            write(errlog(io),*)" inmonth = ",inmonth
          endif;enddo
          stop 1
        endif

        call get_command_argument(5, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read fifth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inday
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inday
        if(inday.lt.1.or.&
           inday.gt.31)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: day must be in range 1-31"
            write(errlog(io),*)" inday = ",inday
          endif;enddo
          stop 1
        endif

        call get_command_argument(6, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read sixth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inhour
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inhour
        if(inhour.lt.0.0_8.or.&
           inhour.gt.48.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: hour must be in range 0.0-48.0"
            write(errlog(io),*)" inhour = ",inhour
          endif;enddo
          stop 1
        endif

        if(nargs.gt.6)then
          ! First optional parameter is the simulation time
          call get_command_argument(7, arg, length=inlen, status=iostatus)
          if(iostatus.ne.0)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read seventh command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)Simtime_in_hours
          linebuffer080 = arg(1:80)
          if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          ! Error-checkSimtime_in_hours
          if(Simtime_in_hours.lt.0.0_8)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: Sim Time must be positive."
              write(errlog(io),*)" inhour = ",inhour
            endif;enddo
            stop 1
          endif

          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Calculating trajectories for ",&
                       Simtime_in_hours," hours."
            endif;enddo

          if(nargs.gt.7)then
            ! Next optional is the number of trajectories followed by the
            ! trajectory levels (in km)
            call get_command_argument(8, arg, length=inlen, status=iostatus)
            if(iostatus.ne.0)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"MR ERROR : Could not read eighth command-line argument"
                write(errlog(io),*)" arg = ",arg
              endif;enddo
              call Print_Usage
            endif
            read(arg,*,iostat=iostatus,iomsg=iomessage)ntraj
            linebuffer080 = arg(1:80)
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
            if(ntraj.le.0)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"MR ERROR reading ntraj."
                write(errlog(io),*)"ntraj = ",ntraj
                write(errlog(io),*)"ntraj must be positive."
              endif;enddo
              stop 1
            elseif(ntraj.gt.9)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then 
                write(errlog(io),*)"MR ERROR: ntraj is currently limited to 9"
              endif;enddo
              stop 1
            endif
            !allocate(OutputLevels(ntraj))
            if(nargs-8.lt.ntraj)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then 
                write(errlog(io),*)"MR ERROR:  There are not enough arguments for ",&
                          ntraj," levels"
                endif;enddo
            elseif(nargs-8.gt.ntraj)then
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"WARNING:  There are more trajectory levels given than needed"
                write(outlog(io),*)"  Expected ntraj = ",ntraj
                write(outlog(io),*)"  Extra command line arguments = ",nargs-8
              endif;enddo
            endif
            do i=1,ntraj
              call get_command_argument(8+i, arg, length=inlen, status=iostatus)
              if(iostatus.ne.0)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR : Could not read 8+i command-line argument"
                  write(errlog(io),*)" arg = ",arg
                endif;enddo
                call Print_Usage
              endif
              read(arg,*,iostat=iostatus,iomsg=iomessage)OutputLevels(i)
              linebuffer080 = arg(1:80)
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              if(OutputLevels(i).le.0.0.or.OutputLevels(i).gt.30.0)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then 
                  write(errlog(io),*)"MR ERROR: trajectory levels must be in range 0-30 km"
                  write(errlog(io),*)"          Failing on trajectory ",i,OutputLevels(i)
                endif;enddo
                stop 1
              endif
            enddo
          endif
        else
          ! Default simulation time is 24 hours
          Simtime_in_hours = 24.0_8
        endif
        if(ntraj.eq.0)then
          ! These are the default trajectory levels if none are specified
          ntraj = 6
          OutputLevels(1) =  1.524_4 !  5000 ft
          OutputLevels(2) =  3.048_4 ! 10000 ft
          OutputLevels(3) =  6.096_4 ! 20000 ft
          OutputLevels(4) =  9.144_4 ! 30000 ft
          OutputLevels(5) = 12.192_4 ! 40000 ft
          OutputLevels(6) = 15.240_4 ! 50000 ft
        endif

        ! Now we need to set the projection for the computational grid, which
        ! for the command-line runs will always be lon/lat
        PJ_iprojflag = 1
        PJ_lam0      = -105.0_8
        PJ_phi0      = 90.0_8
        PJ_phi1      = 90.0_8
        PJ_phi2      = 90.0_8
        PJ_k0        = 0.933_8
        PJ_Re        = 6371.229_8

      elseif(nargs.eq.1)then
        ! we're using a control file.  This is the most general case where non-
        ! GFS and NCEP wind files can be used
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Reading control file"
        endif;enddo
        call get_command_argument(1, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read first command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage) infile
        linebuffer080 = arg(1:80)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        inquire( file=infile, exist=IsThere )
        if(.not.IsThere)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Cannot find input file"
          endif;enddo
          stop 1
        endif
        open(unit=fid_ctrlfile,file=infile,status='old',err=1900)
        ! Line 1: lon, lat
        !  Note that input coordinates are always lon,lat.
        !  If the windfile is projected, travectories will be calcualted on the projected
        !  grid and reported back as lon,lat or whatever is the projection on line 8
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) inlon, inlat
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(inlon.lt.-360.0_8.or.inlon.gt.360.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Longitude must be in range -360->360"
          endif;enddo
          stop 1
        endif
        if(inlon.lt.0.0_8.or.inlon.gt.360.0_8)inlon=mod(inlon+360.0_8,360.0_8)
        if(inlat.lt.-90.0_8.or.inlat.gt.90.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Latitude must be in range -90->90"
          endif;enddo
          stop 1
        endif

        ! Line 2: YYYY MM DD HH.H
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) inyear,inmonth,inday,inhour
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inmonth
        if(inmonth.lt.1.or.&
           inmonth.gt.12)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: month must be in range 1-12"
            write(errlog(io),*)" inmonth = ",inmonth
          endif;enddo
          stop 1
        endif
        ! Error-check inday
        if(inday.lt.1.or.&
           inday.gt.31)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: day must be in range 1-31"
            write(errlog(io),*)" inday = ",inday
          endif;enddo
          stop 1
        endif
        ! Error-check inhour
        if(inhour.lt.0.0_8.or.&
           inhour.gt.24.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: hour must be in range 0.0-24.0"
            write(errlog(io),*)" inhour = ",inhour
          endif;enddo
          stop 1
        endif

        ! Line 3: Length of integration in hours
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) Simtime_in_hours
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check Simtime_in_hours
        if(Simtime_in_hours.lt.0.0_8)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: hour must be positive"
            write(errlog(io),*)" Simtime_in_hours = ",Simtime_in_hours
          endif;enddo
          stop 1
        endif

        ! Line 4: Streamline v.s. Streakline
        !   This is where we could put optional parameters on 2d vs 3d or on
        !     Euler vs something higher-order
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) StreamFlag
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(StreamFlag.ne.0.and.StreamFlag.ne.1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: StreamFlag must be 0 or 1"
            write(errlog(io),*)" StreamFlag = ",StreamFlag
          endif;enddo
          stop 1
        endif

        ! Line 5: Output interval in minutes
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) OutStepInc_Minutes
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(OutStepInc_Minutes.lt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: OutStepInc_Minutes must be positive"
            write(errlog(io),*)" OutStepInc_Minutes = ",OutStepInc_Minutes
          endif;enddo
          stop 1
        endif

        ! Line 6: number of trajectories (must be < 10)
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ntraj
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(ntraj.lt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: ntraj must be positive"
            write(errlog(io),*)" ntraj = ",ntraj
          endif;enddo
          stop 1
        endif
        !allocate(OutputLevels(ntraj))

        ! Line 7: Trajectory levels in km
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) OutputLevels(1:ntraj)
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        do i=1,ntraj
          if(OutputLevels(i).lt.0.0_8)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR: OutputLevels must be positive"
              write(errlog(io),*)" i OutputLevels(i) = ",i,OutputLevels(i)
            endif;enddo
            stop 1
          endif
        enddo

        ! Line 8: Projection of computational grid
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        Comp_projection_line = linebuffer080
        call PJ_Set_Proj_Params(Comp_projection_line)
        Comp_iprojflag  = PJ_iprojflag
        Comp_k0         = PJ_k0
        Comp_Re         = PJ_Re
        Comp_lam0       = PJ_lam0
        Comp_lam1       = PJ_lam1
        Comp_lam2       = PJ_lam2
        Comp_phi0       = PJ_phi0
        Comp_phi1       = PJ_phi1
        Comp_phi2       = PJ_phi2
        if (PJ_ilatlonflag.eq.0)then
          IsLatLon_CompGrid = .false.
        else
          IsLatLon_CompGrid = .true.
        endif

        ! Line 9: iwind iwindformat iformat
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iw,iwf,igrid,idf
        if(iostatus.ne.0.or.&
           iw.lt.1.or.iw.gt.5.or.&
           iwf.lt.0.or.iwf.gt.50.or.&
           idf.lt.1.or.idf.gt.5)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)  'error reading iw,iwf,igrid,idf.'
            write(errlog(io),*)  'You gave: ', linebuffer080
            write(errlog(io),*)  'wind format must be one of 1,2,3,4, or 5'
            write(errlog(io),*)  'wind product must be in range 1-50'
            write(errlog(io),*)  'wind data format must be one of 1,2, or 3'
            write(errlog(io),*)  'Program stopped.'
          endif;enddo
          call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          stop 1
        endif

        ! Line 10: autoflag (0 for auto, 1 for specified) [FC_freq] [GFS_Archive_Days]
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) autoflag
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(autoflag.ne.0)then
          read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) autoflag,FC_freq,GFS_Archive_Days
          if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        else
          FC_freq            = 12     ! Number of hours between GFS package downloads
          GFS_Archive_Days   = 14     ! Number of days GFS data are archived
        endif

        ! Line 11: number of windfiles
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iwfiles
        if(iostatus.ne.0.or.&
           iwfiles.lt.1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)  'error reading iwfiles.'
            write(errlog(io),*)  'You gave: ', linebuffer080
            write(errlog(io),*)  'Program stopped.'
          endif;enddo
          call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          stop 1
        endif

        if(inyear.lt.BaseYear.or.inyear-BaseYear.gt.200)then
          ! Reset BaseYear to the start of the century containing the starttime
          BaseYear = inyear - mod(inyear,100)
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: Resetting BaseYear to ",BaseYear
          endif;enddo
        endif

        MR_BaseYear = BaseYear
        MR_useLeap  = useLeap
        MR_useCompH = .false.

        close(fid_ctrlfile)
      endif

      call MR_Set_CompProjection(IsLatLon_CompGrid,PJ_iprojflag,PJ_lam0,&
                                 PJ_phi0,PJ_phi1,PJ_phi2,&
                                 PJ_k0,PJ_Re)

      ! write out values of parameters defining the run
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"inlon              = ",real(inlon,kind=4)
        write(outlog(io),*)"inlat              = ",real(inlat,kind=4)
        write(outlog(io),*)"inyear             = ",inyear
        write(outlog(io),*)"inmonth            = ",inmonth
        write(outlog(io),*)"inday              = ",inday
        write(outlog(io),*)"inhour             = ",real(inhour,kind=4)
        write(outlog(io),*)"Simtime_in_hours   = ",real(Simtime_in_hours,kind=4)
        write(outlog(io),*)"StreamFlag         = ",StreamFlag
        write(outlog(io),*)"OutStepInc_Minutes = ",OutStepInc_Minutes
        write(outlog(io),*)"ntraj              = ",ntraj
        write(outlog(io),*)"OutputLevels       = "
        do i=1,ntraj
          tmp_4 = real(OutputLevels(i),kind=4)
          write(outlog(io),*)"                  ",i," at ",tmp_4,"km (",tmp_4*3280.8_4," ft)."
        enddo
        write(outlog(io),*)"IsLatLon           = ",IsLatLon_CompGrid
        write(outlog(io),*)"iw                 = ",iw
        write(outlog(io),*)"iwf                = ",iwf
        write(outlog(io),*)"igrid              = ",igrid
        write(outlog(io),*)"idf                = ",idf
        write(outlog(io),*)"autoflag           = ",autoflag
        write(outlog(io),*)"FC_freq            = ",FC_freq
        write(outlog(io),*)"GFS_Archive_Days   = ",GFS_Archive_Days
        write(outlog(io),*)"iwfiles            = ",iwfiles
        write(outlog(io),*)"--------------------------------------------------------------"
      endif;enddo

      return

!******************************************************************************
!     ERROR TRAPS

1900  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine Read_ComdLine_InpFile

!##############################################################################
!##############################################################################
!
!  Print_Usage
!
!  This subroutine is called if there is an error reading the command-line.
!  Expected usage is written to stdout and the program exits.
!
!##############################################################################

      subroutine Print_Usage

      use MetReader,       only : &
         MR_nio,VB,errlog,verbosity_error

      implicit none

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)"Too few command-line arguments:"
        write(errlog(io),*)"  Usage: MetTraj_[F,B] lon lat YYYY MM DD HH.H (FC_hours nlev lev1 lev2 ...)"
        write(errlog(io),*)"           lon       = longitude of start point"
        write(errlog(io),*)"           lat       = latitude of start point"
        write(errlog(io),*)"           YYYY      = start year"
        write(errlog(io),*)"           MM        = start month"
        write(errlog(io),*)"           DD        = start day"
        write(errlog(io),*)"           HH.H      = start hour"
        write(errlog(io),*)"           FC_hours  = [Opt] number of hours to calculate"
        write(errlog(io),*)"           nlev      = [Opt] number of levels"
        write(errlog(io),*)"           lev1 lev2 ... = [Opt] list of nlev level heights in km"
      endif;enddo
      stop 1

      end subroutine Print_Usage


!##############################################################################
!##############################################################################
!  GetWindFile
!
!  This subroutine sets the list of windfiles to be used in the calculation.
!  These will either be an explicit list provided by the control file, or 
!  through an assessment of the current GFS and NCEP files on the system.
!
!##############################################################################

      subroutine GetWindFile(inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,&
                                iw,iwf,igrid,idf,iwfiles,&
                                autoflag,FC_freq,GFS_Archive_Days,GFS_FC_TotHours)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,&
         MR_windfiles,MR_BaseYear,MR_useLeap,&
         MR_Comp_StartHour,MR_Comp_Time_in_hours,MR_iwindfiles,MR_iwind,&
           MR_Allocate_FullMetFileList,&
           MR_Read_Met_DimVars,&
           MR_Set_Met_Times,&
           MR_FileIO_Error_Handler

      implicit none

      integer,parameter :: fid_ctrlfile = 10

      integer        ,intent(in) :: inyear
      integer        ,intent(in) :: inmonth
      integer        ,intent(in) :: inday
      real(kind=8)   ,intent(in) :: inhour
      real(kind=8)   ,intent(in) :: Simtime_in_hours
      integer        ,intent(in) :: TrajFlag
      integer        ,intent(inout) :: iw
      integer        ,intent(inout) :: iwf
      integer        ,intent(inout) :: igrid
      integer        ,intent(inout) :: idf
      integer        ,intent(inout) :: iwfiles
      integer        ,intent(in) :: autoflag
      integer        ,intent(in) :: FC_freq
      integer        ,intent(in) :: GFS_Archive_Days
      integer        ,intent(in) :: GFS_FC_TotHours

      character(len=8)   :: date
      character(len=10)  :: time2       ! time argument used to get current
                                        ! date and time.
      character(len=5)   :: zone        ! variables used by the date_and_time subroutine
      integer            :: values(8)   ! return values from date_and_time
      integer            :: timezone    ! timezone of grid relative to UTC

      real(kind=8)      :: StartHour
      real(kind=8)      :: RunStartHour    ! Current UTC time, in hours since MR_BaseYear
      character(len=17) :: RunStartHour_ch
      real(kind=8)      :: Probe_StartHour
      real(kind=8)      :: Probe_EndHour
      real(kind=8)      :: Met_needed_StartHour
      real(kind=8)      :: Met_needed_EndHour

      integer      :: RunStartYear
      integer      :: RunStartMonth
      integer      :: RunStartDay
      integer      :: RunStartHr

      integer      :: FC_Package_year
      integer      :: FC_Package_month
      integer      :: FC_Package_day
      integer      :: FC_Package_hour
      integer      :: FC_hour_int
      real(kind=8) :: FC_Package_StartHour
      real(kind=8) :: FC_Archive_StartHour

      character(len=100) :: WINDROOT = '/data/WindFiles'
      character(len=47) :: string1,string2

      integer :: i,ii
      integer :: FC_year,FC_mon,FC_day
      real(kind=8) :: FC_hour,FC_Package_hour_dp,FC_intvl
      integer :: NumFCpackages

      character (len=130):: testfile
      real(kind=8)       :: FCStartHour
      real(kind=8)       :: FCEndHour
      logical            :: IsThere

      logical,dimension(:),allocatable :: GFS_candidate
      integer,dimension(:),allocatable :: GFS_FC_step_avail
      integer :: OptimalPackageNum
      integer             :: iostatus
      integer             :: inlen
      character(len=120)  :: iomessage
      character(len=100)  :: arg
      character(len=100)  :: infile
      character(len=80 )  :: linebuffer080
      character(len=130)  :: linebuffer130

      integer :: io                           ! Index for output streams

      INTERFACE
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_yyyymmddhhmm_since
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_MonthOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_DayOfEvent
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_HourOfDay
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE

       ! Get the UTC time for program execution
       !   This will be used to determine if gfs or NCEP winds are to be used
      call date_and_time(date,time2,zone,values)
      read(zone,'(i3)') timezone
        ! Find time in UTC
      StartHour = real(values(5)-timezone,kind=8) + &
                  real(values(6)/60.0_8,kind=8)
        ! find time in hours since BaseYear
      RunStartHour = HS_hours_since_baseyear(values(1),values(2),values(3),&
                                             StartHour,MR_BaseYear,MR_useLeap)
        ! get character string
      RunStartHour_ch = HS_yyyymmddhhmm_since(RunStartHour,MR_BaseYear,MR_useLeap)
      read(RunStartHour_ch,'(i4)')    RunStartYear
      read(RunStartHour_ch,'(4x,i2)') RunStartMonth
      read(RunStartHour_ch,'(6x,i2)') RunStartDay
      read(RunStartHour_ch,'(8x,i2)') RunStartHr

      if(autoflag.eq.1)then
        ! We are using the automatic selection for wind files.  We need to first
        ! determine the current date/time and check if the requested start time is
        ! more recent than the length of the GFS archive.  If so, use GFS;
        ! otherwise, use the NCEP 50-year Reanalysis.

        ! Find the start time given on the command line
        Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                  MR_BaseYear,MR_useLeap)
        MR_Comp_StartHour     = Probe_StartHour
        MR_Comp_Time_in_hours = Simtime_in_hours
  
        ! Calculate the earliest Met data needed
        ! We want this to be as close to the beginning of a forecast package as
        ! possible
        if(TrajFlag.gt.0)then
          Probe_EndHour = Probe_StartHour + Simtime_in_hours
        else
          Probe_EndHour = Probe_StartHour - Simtime_in_hours
        endif
        Met_needed_StartHour = min(Probe_StartHour,Probe_EndHour)
        Met_needed_EndHour   = max(Probe_StartHour,Probe_EndHour)

        if(RunStartHour-Met_needed_StartHour.gt.24.0_8*GFS_Archive_Days)then
          ! NCEP case
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Requested start time is too old for GFS archive."
            write(outlog(io),*)"Start time is older than the hardwired threshold of ",&
                      GFS_Archive_Days," days"
            write(outlog(io),*)"Using NCEP 50-year Reanalysis"
          endif;enddo
          iw  = 5
          iwf = 25
          igrid   = 0
          idf     = 2
          iwfiles = 1

          call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

          do i=1,MR_iwindfiles
            write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)),'/NCEP'
          enddo

          FC_year    = HS_YearOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_mon     = HS_MonthOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_day     = HS_DayOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_hour_dp = HS_HourOfDay(RunStartHour,MR_BaseYear,MR_useLeap)
          ! now round down the hour to the forecast package increment
          FC_Package_hour    = floor(FC_Package_hour_dp/FC_freq) * FC_freq

        elseif(inyear.lt.1948)then
          ! Too old for NCEP runs, must use control file
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Requested start time is too old for NCEP Reanalysis."
            write(errlog(io),*)"Please use a control file and specify an older"
            write(errlog(io),*)"product such as NOAA20CRv3."
          endif;enddo
          stop 1
        else
          ! GFS case
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Requested start time is within the GFS archive."
            write(outlog(io),*)"Using archived global forecast data (GFS 0.5-degree)"
          endif;enddo
          if(RunStartHour-Met_needed_StartHour.lt.0.0_8)then
            ! GFS case for future run
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Requested start time is later than current time,"
              write(outlog(io),*)"but it might fit in the current forecast package."
            endif;enddo
            if (Met_needed_EndHour-RunStartHour.ge.real(GFS_FC_TotHours,kind=8))then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)" Run cannot complete with the current FC package."
              endif;enddo
              stop 1
            endif
          endif
          FC_intvl = 3.0_8

          ! calculate the number of forecast packages stored on system
          NumFCpackages = GFS_Archive_Days * (24/FC_freq)
          allocate(GFS_candidate(NumFCpackages))
          allocate(GFS_FC_step_avail(NumFCpackages))
          GFS_candidate(:) = .true. ! every package is a candidate until proven otherwise
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Approximate Number of FC packages on system:",NumFCpackages
            write(outlog(io),*)"Checking to see which packages might work for the"
            write(outlog(io),*)"requested start-time"
          endif;enddo

          ! First, get the start hour of the forecast package immediately before
          ! the execution time
          FC_Package_year    = HS_YearOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_month   = HS_MonthOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_day     = HS_DayOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_hour_dp = HS_HourOfDay(RunStartHour,MR_BaseYear,MR_useLeap)
          ! now round down the hour to the forecast package increment
          FC_Package_hour    = floor(FC_Package_hour_dp/FC_freq) * FC_freq
          ! and get the hours since base time for that FC package
          FC_Package_StartHour = HS_hours_since_baseyear(FC_Package_year,&
                                   FC_Package_month,FC_Package_day,real(FC_Package_hour,kind=8),&
                                   MR_BaseYear,MR_useLeap)
          ! estimate the start time of the oldest forecast package on the system
          FC_Archive_StartHour = FC_Package_StartHour - GFS_Archive_Days*24.0_8

          ! Loop through all the packages and check which ones might span the needed
          ! time range
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'Looping backward through packages'
          endif;enddo
          do i = NumFCpackages,1,-1
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Package # ",i
            endif;enddo
            ! Get the start hour of this forecast package
            FCStartHour = FC_Archive_StartHour + real((i-1)*FC_freq,kind=8)
            if (FCStartHour.gt.Met_needed_StartHour)then
              ! This package starts after the requested time so dismiss it
              GFS_candidate(i)     = .false.
              GFS_FC_step_avail(i) = 0
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"   Package starts too late"
              endif;enddo
              cycle
            endif

            FCEndHour = FCStartHour + real(GFS_FC_TotHours,kind=8)
            if (FCEndHour.lt.Met_needed_EndHour)then
              ! This package ends before the needed time so dismiss it
              GFS_candidate(i)     = .false.
              GFS_FC_step_avail(i) = 0
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"   Package ends too early"
              endif;enddo
              cycle
            endif

            FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)'   This one could work. Testing each file. ', &
                         FC_year,FC_mon,FC_day,real(FC_hour,kind=4)
            endif;enddo

            FC_hour_int = 0

            ! Now we loop through all the files in this package assuming 198 hours
            ! of forecast data and 3 hour steps (67 files overall)
            do ii = 1,67
              if(.not.GFS_candidate(i)) cycle  ! This will be true until a
                                               ! needed file is missing
              FC_hour_int = nint((ii-1)*FC_intvl)
              ! skip all the files that are not needed to bracket the simulation
              if (FCStartHour+FC_hour_int.lt.Met_needed_StartHour-FC_intvl) cycle ! to early
              if (FCStartHour+FC_hour_int.gt.Met_needed_EndHour+FC_intvl)   cycle ! to late

              ! if we are here, then we are inspecting a file that would be needed for
              ! the requested time span, starting with the most recent forecast package
              ! See if the file actually exists
              write(string1,'(a9,I4.4,I2.2,I2.2,I2.2,a1)')'/gfs/gfs.', &
                            FC_year,FC_mon,FC_day,FC_Package_hour,'/'
              write(string2,'(I4.4,I2.2,I2.2,I2.2,a2,I3.3,a3)')&
                            FC_year,FC_mon,FC_day,FC_Package_hour, &
                            '.f',FC_hour_int,'.nc'
    
              write(testfile,*)trim(ADJUSTL(WINDROOT)), &
                                   trim(ADJUSTL(string1)), &
                                   trim(ADJUSTL(string2))
              inquire( file=trim(adjustl(testfile)), exist=IsThere )
              if (IsThere)then
                GFS_FC_step_avail(i) = ii
                FCEndHour = FCStartHour + real(FC_hour_int,kind=8)
                do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                  write(outlog(io),*)"     Found: ",trim(adjustl(testfile))
                endif;enddo
              else
                ! If a needed file is not available, mark the whole package as not a
                ! candidate
                GFS_candidate(i) = .false.
                cycle
              endif
            enddo

            ! If after all that testing, if we still have a valid candidate forecast package,
            ! set this package as the one we will use and exit the do loop
            if (GFS_candidate(i)) then
              OptimalPackageNum = i
              exit
            endif

          enddo ! 1,NumFCpackages

          if (OptimalPackageNum.eq.0)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"No GFS package available on this system will span the"
              write(errlog(io),*)"requested simulation time.  Exiting"
              write(errlog(io),*)"  Sim start time       = ",Met_needed_StartHour
              write(errlog(io),*)"  Sim end time         = ",Met_needed_EndHour
              write(errlog(io),*)"  FC_Archive_StartHour =",FC_Archive_StartHour
              write(errlog(io),*)"  FCStartHour          = ",FCStartHour
              write(errlog(io),*)"  FCEndHour            = ",FCEndHour
            endif;enddo
            stop 1
          endif

          iw      = 4
          iwf     = 20
          igrid   = 0
          idf     = 2
          iwfiles = GFS_FC_step_avail(OptimalPackageNum)
          call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

          !  Now list the windfiles we will use and copy to MR_windfiles
          FCStartHour = FC_Archive_StartHour + real((OptimalPackageNum-1)*FC_freq,kind=8)
          FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq

          do i=1,MR_iwindfiles
            
            FCStartHour = FC_Archive_StartHour + real((OptimalPackageNum-1)*FC_freq,kind=8)

            FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq
            FC_hour_int = nint((i-1)*FC_intvl)
            write(string1,'(a9,I4.4,I2.2,I2.2,I2.2,a1)')'/gfs/gfs.', &
                          FC_year,FC_mon,FC_day,FC_Package_hour,'/'
            write(string2,'(I4.4,I2.2,I2.2,I2.2,a2,I3.3,a3)')&
                          FC_year,FC_mon,FC_day,FC_Package_hour, &
                          '.f',FC_hour_int,'.nc'

            write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)), &
                                 trim(ADJUSTL(string1)), &
                                 trim(ADJUSTL(string2))
          enddo
        endif  ! automatic cases: NCEP, then GFS
      else

        Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                  MR_BaseYear,MR_useLeap)
        MR_Comp_StartHour     = Probe_StartHour
        MR_Comp_Time_in_hours = Simtime_in_hours

        ! Calculate the earliest Met data needed
        if(TrajFlag.gt.0)then
          Probe_EndHour = Probe_StartHour + Simtime_in_hours
        else
          Probe_EndHour = Probe_StartHour - Simtime_in_hours
        endif
        Met_needed_StartHour = min(Probe_StartHour,Probe_EndHour)
        Met_needed_EndHour   = max(Probe_StartHour,Probe_EndHour)

        ! This is a run controlled by an input file.
        call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
        ! Reread the input file to get the windfile names

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Reading control file" 
        endif;enddo
        call get_command_argument(1, arg, length=inlen, status=iostatus)
        read(arg,*,iostat=iostatus,iomsg=iomessage) infile
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,arg(1:80),iomessage)
        inquire( file=infile, exist=IsThere )
        if(.not.IsThere)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Cannot find input file"
          endif;enddo
          stop 1
        endif
        open(unit=fid_ctrlfile,file=infile,status='old')
        ! Line 1: lon, lat
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 2: YYYY MM DD HH.H
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 3: Length of integration in hours
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 4: Streamline v.s. Streakline
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 5: Output interval in minutes
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 6: number of trajectories (must be < 10)
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 7: Trajectory levels in km
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 8: Projection of computational grid
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 9: iwind iwindformat iformat
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 10: autoflag (0 for auto, 1 for specified) [FC_freq] [GFS_Archive_Days]
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Line 11 -> 11 + #of windfiles
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        if(MR_iwind.eq.5)then
          ! For NCEP 2.5 degree (25), NOAA product (27), ERA5 (29), or ERA-20C (30)
          ! just read the path to the files
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer130
          if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer130(1:80),iomessage)
          read(linebuffer130,'(a130)',iostat=iostatus,iomsg=iomessage) MR_windfiles(1)
          if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer130(1:80),iomessage)
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Read windfile name: ",adjustl(trim(MR_windfiles(1)))
          endif;enddo
        else
          ! For all other iwf (MR_iwindformats), read the full list
          do i=1,iwfiles
            read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer130
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer130(1:80),iomessage)
            read(linebuffer130,'(a130)',iostat=iostatus,iomsg=iomessage) MR_windfiles(i)
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer130(1:80),iomessage)
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Read windfile name: ",adjustl(trim(MR_windfiles(i)))
          endif;enddo
          enddo
        endif
        close(fid_ctrlfile)

      endif

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(FC_year)

      call MR_Set_Met_Times(Met_needed_StartHour, Simtime_in_hours)

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Traj time: ",inyear,inmonth,inday,inhour
        write(outlog(io),*)"Now      : ",RunStartYear,RunStartMonth,RunStartDay,RunStartHr
      endif;enddo

      end subroutine GetWindFile

!##############################################################################
!##############################################################################
!  Integrate_ConstH_Traj
!
!  This subroutine actually performs the integration (forward or backward) and
!  writes the trajectory data to [f,b]traj[1-9].dat
!
!##############################################################################

      subroutine Integrate_ConstH_Traj(IsGlobal,inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,output_interv)

      use MetReader,       only : &
         MR_nio,VB,outlog,verbosity_info,&
         MR_dum3d_CompH,nx_comp,ny_comp,dx_met_const,dy_met_const,&
         IsRegular_MetGrid,MR_BaseYear,MR_useLeap,MR_iMetStep_Now,&
         MR_MetSteps_Total,MR_MetStep_Interval,MR_MetStep_Hour_since_baseyear,&
         x_comp_sp,y_comp_sp,&
           MR_Allocate_FullMetFileList,&
           MR_Read_HGT_arrays,&
           MR_Read_3d_Met_Variable_to_CompH

      implicit none

      logical     , intent(in)      :: IsGlobal
      real(kind=8), intent(in)      :: inlon
      real(kind=8), intent(in)      :: inlat
      integer     , intent(in)      :: inyear
      integer     , intent(in)      :: inmonth
      integer     , intent(in)      :: inday
      real(kind=8), intent(in)      :: inhour
      real(kind=8), intent(in)      :: Simtime_in_hours
      integer     , intent(in)      :: TrajFlag
      integer     , intent(in)      :: ntraj
      integer     , intent(in)      :: output_interv

      real(kind=8), parameter :: PI        = 3.141592653589793
      real(kind=8), parameter :: DEG2RAD   = 1.7453292519943295e-2
      real(kind=8), parameter :: KM_2_M    = 1.0e3
      real(kind=8), parameter :: RAD_EARTH = 6371.229 ! Radius of Earth in km

      real(kind=8)       :: Probe_StartHour
      integer            :: ivar
      integer            :: kk
      integer,dimension(9) :: ixold,iyold
      integer            :: ix,iy

      real(kind=4),dimension(:,:,:),allocatable :: Vx_meso_last_step_MetH_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vx_meso_next_step_MetH_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vy_meso_last_step_MetH_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vy_meso_next_step_MetH_sp

      real(kind=8) :: tfrac,tc
      real(kind=8) :: xfrac,xc,yfrac,yc
      real(kind=4) :: a1,a2,a3,a4
      real(kind=8), dimension(ntraj) :: x1,y1

      real(kind=8),dimension(:,:,:,:),allocatable :: Vx_full
      real(kind=8),dimension(:,:,:,:),allocatable :: Vy_full
      real(kind=8),dimension(:)      ,allocatable :: Step_Time_since1900
      real(kind=8),dimension(:,:,:)  ,allocatable :: dvxdt
      real(kind=8),dimension(:,:,:)  ,allocatable :: dvydt

      integer      :: istep,stepindx
      integer      :: ti,iit,it
      real(kind=8) :: vx1,vx2,vx3,vx4
      real(kind=8) :: vy1,vy2,vy3,vy4
      real(kind=8),dimension(2)  :: vel_1
      real(kind=8) :: dt
      real(kind=8) :: mstodeghr
      real(kind=8) :: t1
      real(kind=8) :: x_fin,y_fin
      real(kind=8) :: xstep,ystep
      real(kind=8) :: lonmin,lonmax,latmin,latmax

      integer :: io                           ! Index for output streams

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE

      allocate(Vx_meso_last_step_MetH_sp(nx_comp,ny_comp,ntraj))
      allocate(Vx_meso_next_step_MetH_sp(nx_comp,ny_comp,ntraj))
      allocate(Vy_meso_last_step_MetH_sp(nx_comp,ny_comp,ntraj))
      allocate(Vy_meso_next_step_MetH_sp(nx_comp,ny_comp,ntraj))

      ! These store just the layers relevant, but for all time
      allocate(Vx_full(0:nx_comp+1,0:ny_comp+1,ntraj,MR_MetSteps_Total))
      allocate(Vy_full(0:nx_comp+1,0:ny_comp+1,ntraj,MR_MetSteps_Total))
      allocate(Step_Time_since1900(MR_MetSteps_Total))
      ! These are needed for each integration point
      allocate(dvxdt(0:nx_comp+1,0:ny_comp+1,ntraj))
      allocate(dvydt(0:nx_comp+1,0:ny_comp+1,ntraj))

      lonmin = 360.0_8
      lonmax =   0.0_8
      latmin =  90.0_8
      latmax = -90.0_8

       ! Load the full sub-grid for all times
        ! First load the Met grids for Geopotential
      if(TrajFlag.gt.0)then
        ! Foreward trajectory
        MR_iMetStep_Now = 1
      else
        ! Backward trajectory
        MR_iMetStep_Now = MR_MetSteps_Total-1 
      endif
      !call MR_Read_HGT_arrays(MR_iMetStep_Now)

      ! Get the fractional time between forecast steps
      Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      tfrac = (Probe_StartHour-MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now))/ &
               MR_MetStep_Interval(MR_iMetStep_Now)
      tc    = 1.0_8-tfrac

      ! Loop through all the steps in proper chronological order, but store in
      ! Vx_full, and Vy_full in order of integration (forward or backward)
      do istep = 1,MR_MetSteps_Total
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Reading step of ",istep,MR_MetSteps_Total
        endif;enddo
        MR_iMetStep_Now = istep
        if(TrajFlag.gt.0)then
          stepindx=istep
        else
          stepindx=MR_MetSteps_Total-istep+1
        endif
        Step_Time_since1900(stepindx) = MR_MetStep_Hour_since_baseyear(istep)
        if(istep.lt.MR_MetSteps_Total)call MR_Read_HGT_arrays(istep)
        ivar = 2 ! Vx

        call MR_Read_3d_Met_Variable_to_CompH(ivar,istep)
        Vx_full(1:nx_comp,1:ny_comp,:,stepindx) = &
          MR_dum3d_CompH(1:nx_comp,1:ny_comp,:)
        ivar = 3 ! Vy
        call MR_Read_3d_Met_Variable_to_CompH(ivar,istep)
        Vy_full(1:nx_comp,1:ny_comp,:,stepindx) = &
          MR_dum3d_CompH(1:nx_comp,1:ny_comp,:)

      enddo

      ! Finally B.C.'s
      Vx_full(:,          0,:,:)=Vx_full(:,        1,:,:)
      Vx_full(:,ny_comp+1,:,:)=Vx_full(:,ny_comp,:,:)
      Vy_full(:,          0,:,:)=Vy_full(:,        1,:,:)
      Vy_full(:,ny_comp+1,:,:)=Vy_full(:,ny_comp,:,:)

      if(IsGlobal)then
        Vx_full(          0,:,:,:)=Vx_full(nx_comp,:,:,:)
        Vx_full(nx_comp+1,:,:,:)=Vx_full(        1,:,:,:)
        Vy_full(          0,:,:,:)=Vy_full(nx_comp,:,:,:)
        Vy_full(nx_comp+1,:,:,:)=Vy_full(        1,:,:,:)
      else
        Vx_full(          0,:,:,:)=Vx_full(        1,:,:,:)
        Vx_full(nx_comp+1,:,:,:)=Vx_full(nx_comp,:,:,:)
        Vy_full(          0,:,:,:)=Vy_full(        1,:,:,:)
        Vy_full(nx_comp+1,:,:,:)=Vy_full(nx_comp,:,:,:)
      endif

      ! We now have the full x,y,z,vx,vy data needed from the Met file
      ! for the full forward/backward simulation

      ! Initialize the start coordinates
      x1(:) = inlon
      y1(:) = inlat
      t1    = Probe_StartHour
      it    = 1

      ! Assume an integration step of 1 min and a max v of around 100m/s
      if(TrajFlag.ge.0)then
        dt = 1.0_8/60.0_8
      else
        dt = -1.0_8/60.0_8
      endif

      mstodeghr = 3600.0_8*360.0_8/(2.0_8*PI*RAD_EARTH*KM_2_M)

      if(TrajFlag.ge.0)then
        if(ntraj.ge.1)open(unit=21,file='ftraj1.dat')
        if(ntraj.ge.2)open(unit=22,file='ftraj2.dat')
        if(ntraj.ge.3)open(unit=23,file='ftraj3.dat')
        if(ntraj.ge.4)open(unit=24,file='ftraj4.dat')
        if(ntraj.ge.5)open(unit=25,file='ftraj5.dat')
        if(ntraj.ge.6)open(unit=26,file='ftraj6.dat')
        if(ntraj.ge.7)open(unit=27,file='ftraj7.dat')
        if(ntraj.ge.8)open(unit=28,file='ftraj8.dat')
        if(ntraj.ge.9)open(unit=29,file='ftraj9.dat')
      else
        if(ntraj.ge.1)open(unit=21,file='btraj1.dat')
        if(ntraj.ge.2)open(unit=22,file='btraj2.dat')
        if(ntraj.ge.3)open(unit=23,file='btraj3.dat')
        if(ntraj.ge.4)open(unit=24,file='btraj4.dat')
        if(ntraj.ge.5)open(unit=25,file='btraj5.dat')
        if(ntraj.ge.6)open(unit=26,file='btraj6.dat')
        if(ntraj.ge.7)open(unit=27,file='btraj7.dat')
        if(ntraj.ge.8)open(unit=28,file='btraj8.dat')
        if(ntraj.ge.9)open(unit=29,file='btraj9.dat')
      endif

      if(ntraj.ge.1)write(21,*)real(x1(1),kind=4),real(y1(1),kind=4)
      if(ntraj.ge.2)write(22,*)real(x1(2),kind=4),real(y1(2),kind=4)
      if(ntraj.ge.3)write(23,*)real(x1(3),kind=4),real(y1(3),kind=4)
      if(ntraj.ge.4)write(24,*)real(x1(4),kind=4),real(y1(4),kind=4)
      if(ntraj.ge.5)write(25,*)real(x1(5),kind=4),real(y1(5),kind=4)
      if(ntraj.ge.6)write(26,*)real(x1(6),kind=4),real(y1(6),kind=4)
      if(ntraj.ge.7)write(27,*)real(x1(7),kind=4),real(y1(7),kind=4)
      if(ntraj.ge.8)write(28,*)real(x1(8),kind=4),real(y1(8),kind=4)
      if(ntraj.ge.9)write(29,*)real(x1(9),kind=4),real(y1(9),kind=4)

      ! Find location of initial trajectory heights
      ! Get interpolation coefficients
      it = 1
      if(TrajFlag.gt.0)then
        tfrac = (t1-Step_Time_since1900(it))/MR_MetStep_Interval(it)
      else
        tfrac = -(t1-Step_Time_since1900(it))/MR_MetStep_Interval(it)
      endif

      ! integrate out Simtime_in_hours hours
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now integrating out ",abs(int(Simtime_in_hours/dt))," steps"
      endif;enddo
      it = 0
      ixold(:)=1
      iyold(:)=1
      do ti = 1,abs(int(Simtime_in_hours/dt))
        if(TrajFlag.gt.0)then
          ! Get the interval by assuming all MetStep_Intervals are the same
          iit = floor((t1-Step_Time_since1900(1))/MR_MetStep_Interval(1)) + 1
          tfrac = (t1-Step_Time_since1900(iit))/MR_MetStep_Interval(1)
        else
          ! Get the interval by assuming all MetStep_Intervals are the same
          iit = floor(-(t1-Step_Time_since1900(1))/MR_MetStep_Interval(1)) + 1
          tfrac = -(t1-Step_Time_since1900(iit))/MR_MetStep_Interval(1)
        endif

        if(iit.ne.it)then
          it = iit
           ! Get the change in velocity in m/s/hr
          dvxdt(:,:,:) = (Vx_full(:,:,:,it+1)-Vx_full(:,:,:,it))/MR_MetStep_Interval(it)
          dvydt(:,:,:) = (Vy_full(:,:,:,it+1)-Vy_full(:,:,:,it))/MR_MetStep_Interval(it)
        endif

        do kk = 1,ntraj
          ! Get current time and position indecies
          if(IsRegular_MetGrid)then
            ! For regular grids, finding the indecies is trivial
            ix = floor((x1(kk)-x_comp_sp(1))/abs(dx_met_const)) + 1
            iy = floor((y1(kk)-y_comp_sp(1))/abs(dy_met_const)) + 1

          else
            ! For non-regular grids (e.g. Gaussian), we need to march over the
            ! subgrid to find the index of the current point.  This could be
            ! faster.
            do ix=max(ixold(kk),1),nx_comp-1
              if (x1(kk).ge.x_comp_sp(ix).and.x1(kk).lt.x_comp_sp(ix+1))then
                exit
              endif
            enddo
            do iy=max(iyold(kk),1),ny_comp-1
              if (y1(kk).ge.y_comp_sp(iy).and.y1(kk).lt.y_comp_sp(iy+1))then
                exit
              endif
            enddo

          endif
          if(.not.IsGlobal)then
            ! Skip over points that leave the domain
            if(ix.le.0.or.ix.ge.nx_comp)cycle
            if(iy.le.0.or.iy.ge.ny_comp)cycle
          endif

          ! Get the fractional position within the cell
          xfrac = (x1(kk)-x_comp_sp(ix))/(x_comp_sp(ix+1)-x_comp_sp(ix))
          yfrac = (y1(kk)-y_comp_sp(iy))/(y_comp_sp(iy+1)-y_comp_sp(iy))
          xc = 1.0_4-xfrac
          yc = 1.0_4-yfrac
          ! Build interpolation coefficients
          a1 = real(xc*yc,kind=4)
          a2 = real(xfrac*yc,kind=4)
          a3 = real(xfrac*yfrac,kind=4)
          a4 = real(yfrac*xc,kind=4)

          ! Corner velocities for current time
          vx1 = Vx_full(ix  ,iy  ,kk,it) + tfrac*dvxdt(ix  ,iy  ,kk)
          vx2 = Vx_full(ix+1,iy  ,kk,it) + tfrac*dvxdt(ix+1,iy  ,kk)
          vx3 = Vx_full(ix+1,iy+1,kk,it) + tfrac*dvxdt(ix+1,iy+1,kk)
          vx4 = Vx_full(ix  ,iy+1,kk,it) + tfrac*dvxdt(ix  ,iy+1,kk)
          vy1 = Vy_full(ix  ,iy  ,kk,it) + tfrac*dvydt(ix  ,iy  ,kk)
          vy2 = Vy_full(ix+1,iy  ,kk,it) + tfrac*dvydt(ix+1,iy  ,kk)
          vy3 = Vy_full(ix+1,iy+1,kk,it) + tfrac*dvydt(ix+1,iy+1,kk)
          vy4 = Vy_full(ix  ,iy+1,kk,it) + tfrac*dvydt(ix  ,iy+1,kk)

          ! Interpolate velocity onto current position and time (in m/s)
          vel_1(1) = (a1*vx1+a2*vx2+a3*vx3+a4*vx4)
          vel_1(2) = (a1*vy1+a2*vy2+a3*vy3+a4*vy4)
          !  now convert to deg/hr
          vel_1(1) = vel_1(1)*mstodeghr/sin((90.0_4-y1(kk))*DEG2RAD)
          vel_1(2) = vel_1(2)*mstodeghr

          ! Now advect via Forward Euler
          xstep = vel_1(1) * dt
          ystep = vel_1(2) * dt
          x_fin = x1(kk) + xstep
          y_fin = y1(kk) + ystep
          if (x_fin.ge.360.0_8)x_fin=x_fin - 360.0_8
          if (x_fin.lt.0.0_8)x_fin=x_fin + 360.0_8

          x1(kk) = x_fin
          y1(kk) = y_fin
        enddo

        t1 = t1 + dt
        if(mod(ti,output_interv).eq.0)then
        !if(mod(ti,60).eq.0)then
          do kk = 1,ntraj
            if(kk.eq.1)write(21,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.2)write(22,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.3)write(23,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.4)write(24,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.5)write(25,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.6)write(26,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.7)write(27,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.8)write(28,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.9)write(29,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            ! Check min/max of trajectories

            if(.not.IsGlobal)then
              if(x1(kk).lt.lonmin)lonmin=x1(kk)
              if(x1(kk).gt.lonmax)lonmax=x1(kk)
              if(y1(kk).lt.latmin)latmin=y1(kk)
              if(y1(kk).gt.latmax)latmax=y1(kk)
            endif
          enddo
        endif
      enddo
      open(unit=40,file='map_range_traj.txt')
      write(40,*)real(lonmin,kind=4),&
                 real(lonmax,kind=4),&
                 real(latmin,kind=4),&
                 real(latmax,kind=4),&
                  real(inlon,kind=4),&
                  real(inlat,kind=4)
      close(40)

      if(ntraj.ge.1)close(21)
      if(ntraj.ge.2)close(22)
      if(ntraj.ge.3)close(23)
      if(ntraj.ge.4)close(24)
      if(ntraj.ge.5)close(25)
      if(ntraj.ge.6)close(26)
      if(ntraj.ge.7)close(27)
      if(ntraj.ge.8)close(28)
      if(ntraj.ge.9)close(29)

      end subroutine Integrate_ConstH_Traj

!##############################################################################
!##############################################################################
