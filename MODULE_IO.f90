MODULE MODULE_IO

USE VARIABLES

contains

!####################################################################
! Subroutine to read in and pass meteorological data (e.g. rainfall).
!####################################################################

      subroutine rdatmo(nrow,ncol,ipixnum,iyear,iday,ihour,i,ATMOS)

      implicit none
      type (GRID_MET_template) :: ATMOS(nrow*ncol)
      type (OPTIONS_template) :: OPTIONS
      integer :: nrow,ncol,ipixnum(nrow,ncol)
      integer :: iyear,iday,ihour
      integer :: forcingnvars,i
      real,dimension(:,:,:),allocatable :: TempArray
      forcingnvars = 7
      OPTIONS%ncol = ncol
      OPTIONS%nrow = nrow
      allocate(TempArray(OPTIONS%ncol,OPTIONS%nrow,forcingnvars))

! ====================================================================
! Read year, day and hour
! ====================================================================

      read(61,*) iyear,iday,ihour

! ####################################################################
! Read all variables in at once for each time step
! ####################################################################

      read(1004,rec=i) TempArray(:,:,:)

! Longwave Radiation

      call rdforc(ATMOS%rld,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,1))

! Air Pressure

      call rdforc(ATMOS%press,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,2))

! Relative Humidity

      call rdforc(ATMOS%rh,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,3))

! Shortwave Radiation

      call rdforc(ATMOS%rsd,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,4))

! Air Temperature

      call rdforc(ATMOS%tdry,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,5))

! Wind Speed

      call rdforc(ATMOS%uzw,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,6))

! Precipitation

      call rdforc(ATMOS%pptms,OPTIONS%nrow,OPTIONS%ncol,ipixnum,TempArray(:,:,7))
       
      contains

              subroutine rdforc(pixdat,nrow,ncol,ipixnum,data_in)

              implicit none
              integer,intent(in) :: nrow,ncol
              integer :: l,m
              real,dimension(:,:),intent(in) :: data_in(:,:)
              integer,dimension(:,:),intent(in) :: ipixnum
              real*8,dimension(:),intent(inout) :: pixdat

                do l=1,ncol
                        do m=1,nrow
                                if (ipixnum(m,l) .gt. 0)pixdat(ipixnum(m,l)) = data_in(l,m)
                        enddo
                enddo

              end subroutine

      end subroutine

! ####################################################################
! Subroutine to open input/output files, read and initialize time 
! in-variant data.
! ####################################################################

      subroutine rddata(OPTIONS,STORM_PARAM,TOPMODEL_PARAM,&
                SOIL_MOISTURE,&
                INF_PARAM,SNOW_VARS,GRID,REG,GLOBAL,CAT)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/rddata.h"
      type (OPTIONS_template) :: OPTIONS
      type (STORM_PARAM_template) :: STORM_PARAM
      type (TOPMODEL_PARAM_template) :: TOPMODEL_PARAM
      type (SOIL_PARAM_template) :: SOIL_PARAM
      type (SOIL_MOISTURE_template) :: SOIL_MOISTURE
      type (INF_PARAM_template) :: INF_PARAM
      type (SNOW_VARS_template) :: SNOW_VARS
      type (GRID_template),dimension(:),allocatable :: GRID
      type (REGIONAL_template) :: REG
      type (GLOBAL_template) :: GLOBAL
      type (CATCHMENT_template),dimension(:),allocatable :: CAT

! ====================================================================
! Open input/output files and set variable to control output file
! printing.
! ====================================================================

      call filopn(OPTIONS%iprn,OPTIONS%ioutst,OPTIONS%ioutsp,OPTIONS%iouten,OPTIONS%nseries,OPTIONS%icurser,OPTIONS%fnimg)

! ====================================================================
! Read in the way the input images are named.
! ====================================================================

      read(61,*) OPTIONS%img_opt

      if (OPTIONS%img_opt.eq.0) then

         print*, 'rddata : Images are listed aINF_PARAM%ccording to timestep'

      endif

      if (OPTIONS%img_opt.eq.1) then

         print*, 'rddata : Images are listed aINF_PARAM%ccording to year_day_hour'

      endif

! ====================================================================
! Read in simulation time constants and control variables.
! ====================================================================

      read(1000,*) OPTIONS%ndata
      read(1000,*) STORM_PARAM%dt
      read(1000,*) STORM_PARAM%endstm
      read(1000,*) iophd

      print*, 'rddata:  Done reading time parameters'
      print*, 'rddata:  Total time steps = ',OPTIONS%ndata

! ====================================================================
! Write output file headers if requested.
! ====================================================================

      if (iophd.eq.0) call hdprnt(OPTIONS%iprn)

      print*, 'rddata:  Done printing headers'

! ====================================================================
! Read and initialize topmodel parameters, atb distribution and
! initial water table depth.
! ====================================================================

      call rdtpmd(OPTIONS%iopbf,OPTIONS%iopwt0,OPTIONS%ncatch,OPTIONS%nrow,&
       OPTIONS%ncol,STORM_PARAM%pixsiz,OPTIONS%ipixnum,OPTIONS%iprn,OPTIONS%ixpix,OPTIONS%iypix,&
       OPTIONS%npix,TOPMODEL_PARAM%q0,TOPMODEL_PARAM%ff,INF_PARAM%qb0,&
       TOPMODEL_PARAM%dd,TOPMODEL_PARAM%xlength,TOPMODEL_PARAM%basink,TOPMODEL_PARAM%xlamda,OPTIONS%icatch,&
       TOPMODEL_PARAM%area,TOPMODEL_PARAM%atanb,TOPMODEL_PARAM%dtil,TOPMODEL_PARAM%zbar1,&
       TOPMODEL_PARAM%iwel,TOPMODEL_PARAM%wslp,&
       OPTIONS%lat_deg,OPTIONS%lat_min,OPTIONS%lng_deg,&
       OPTIONS%lng_min,OPTIONS%lng_mer,OPTIONS%rlatitude,&
       OPTIONS%rlongitude,OPTIONS%rlng_merid,GRID,CAT)

      print*,'rddata:  Done reading TOPMODEL parameters'

! ====================================================================
! Read in and initialize vegatation parameters.
! ====================================================================

      call rdveg(OPTIONS%npix,OPTIONS%nrow,OPTIONS%ncol,OPTIONS%ilandc,&
       OPTIONS%ipixnum,OPTIONS%nlandc,OPTIONS%iopveg,OPTIONS%ivgtyp,OPTIONS%iprn,&
       GRID%VEG%xlai,GRID%VEG%xlai_wsc,GRID%VEG%albd,&
       GRID%VEG%albw,GRID%VEG%emiss,GRID%VEG%za,&
       GRID%VEG%zww,GRID%VEG%z0m,GRID%VEG%z0h,&
       GRID%VEG%zpd,GRID%VEG%rsmin,GRID%VEG%rsmax,&
       GRID%VEG%Rpl,GRID%VEG%f3vpdpar,GRID%VEG%f4temppar,&
       GRID%VEG%trefk,GRID%VEG%tcbeta,GRID%VEG%extinct,&
       GRID%VEG%canclos,GRID%VEG%Tslope1,GRID%VEG%Tint1,&
       GRID%VEG%Tslope2,GRID%VEG%Tint2,GRID%VEG%Twslope1,&
       GRID%VEG%Twint1,GRID%VEG%Twslope2,GRID%VEG%Twint2,&
       GRID%VEG%Tsep,GRID%VEG%Twsep,GRID%VEG%eps,&
       GRID%VEG%rtact,GRID%VEG%rtdens,GRID%VEG%rtres,&
       GRID%VEG%psicri,GRID%VEG%rescan,GRID%VEG%respla,&
       GRID%VEG%wsc,GRID%VEG%wcip1,&
       STORM_PARAM%pixsiz,TOPMODEL_PARAM%area,CAT%fbs,&
       REG%fbsrg,OPTIONS%ncatch)

      print*,'rddata:  Done reading vegetation parameters'

! ====================================================================
! Read in soil parameters and root and transmission zone information.
! ====================================================================

      call rdsoil(OPTIONS%nsoil,OPTIONS%irestype,OPTIONS%ikopt,GLOBAL%zrzmax,OPTIONS%iopsmini,SOIL_PARAM%smpet0,&
       OPTIONS%isoil,OPTIONS%nrow,OPTIONS%ncol,OPTIONS%ipixnum,SOIL_PARAM%bcbeta,&
       SOIL_PARAM%psic,SOIL_PARAM%thetas,SOIL_PARAM%thetar,SOIL_PARAM%xk0,SOIL_PARAM%zdeep,SOIL_PARAM%tdeep,SOIL_PARAM%zmid,&
       SOIL_PARAM%tmid0,SOIL_PARAM%rocpsoil,SOIL_PARAM%quartz,OPTIONS%ifcoarse,&
       SOIL_PARAM%srespar1,SOIL_PARAM%srespar2,SOIL_PARAM%srespar3,SOIL_PARAM%a_ice,SOIL_PARAM%b_ice,&
       SOIL_PARAM%bulk_dens,SOIL_PARAM%amp,SOIL_PARAM%phase,SOIL_PARAM%shift,&
       OPTIONS%inc_frozen,SOIL_PARAM%bcgamm,SOIL_PARAM%par,SOIL_PARAM%corr,OPTIONS%idifind,&
       OPTIONS%ncatch,OPTIONS%icatch,STORM_PARAM%pixsiz,TOPMODEL_PARAM%area,&
       OPTIONS%npix,SOIL_PARAM%psicav,OPTIONS%iprn,GRID%VEG%tc,GRID%VEG%tw)

      print*,'rddata:  Done reading soil parameters'

! ====================================================================
! Read in options for energy balance and calculate soil thermal
! conductivity. 
! ====================================================================

      call rdebal(OPTIONS%ioppet,OPTIONS%iopwv,OPTIONS%iopstab,OPTIONS%iopgveg,&
                  OPTIONS%iopthermc,OPTIONS%iopthermc_v,OPTIONS%maxnri,STORM_PARAM%toleb)

      print*,'rddata:  Done reading energy balance parameters'

! ====================================================================
! Read in the mode in which to run the program.
! ====================================================================

      OPTIONS%MODE = 1
      OPTIONS%FRCOV = 0
      OPTIONS%frcbeta = 999

      if (OPTIONS%MODE.eq.1) print*,'rddata:  Running the model in dist. mode'
      if (OPTIONS%FRCOV.eq.0) print*,'rddata:  Frac. rain.l cover not included'

! ====================================================================
! Initialize the simulation sum variables and storm.and.&
! interstorm flags and times.
! ====================================================================

      call inisim(OPTIONS%iopsmini,OPTIONS%nrow,OPTIONS%ncol,OPTIONS%ipixnum,OPTIONS%ilandc,&
       OPTIONS%npix,OPTIONS%inc_frozen,STORM_PARAM%istorm,&
       STORM_PARAM%intstm,STORM_PARAM%istmst,intstp,STORM_PARAM%istorm_moss,&
       STORM_PARAM%intstm_moss,STORM_PARAM%istmst_moss,STORM_PARAM%intstp_moss,&
       OPTIONS%isoil,OPTIONS%idifind,SOIL_PARAM%smpet0,r_mossmpet0,STORM_PARAM%endstm,&
       SOIL_MOISTURE%rzsm1,SOIL_MOISTURE%tzsm1,SOIL_MOISTURE%r_mossm1,&
       SOIL_MOISTURE%r_mossm,SOIL_MOISTURE%rzsm1_u,SOIL_MOISTURE%tzsm1_u,&
       SOIL_MOISTURE%rzsm1_f,SOIL_MOISTURE%tzsm1_f,SOIL_MOISTURE%r_mossm1_u,&
       SOIL_MOISTURE%r_mossm_u,&
       SOIL_MOISTURE%r_mossm1_f,SOIL_MOISTURE%r_mossm_f,SOIL_MOISTURE%rzdthetaidt,&
       SOIL_MOISTURE%tzdthetaidt,SOIL_MOISTURE%zmoss,r_moss_depth,&
       thetas_moss,INF_PARAM%xintst,INF_PARAM%xintst_moss,INF_PARAM%cuminf,SOIL_PARAM%xk0,SOIL_PARAM%psic,&
       SOIL_PARAM%thetas,SOIL_PARAM%thetar,SOIL_PARAM%bcgamm,&
       bcbeta,INF_PARAM%sorp,INF_PARAM%cc,STORM_PARAM%dt,INF_PARAM%sesq,SOIL_PARAM%corr,SOIL_PARAM%par,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,r_MeltEnergy_us,Outflow_us,&
       PackWater,SurfWater,Swq,VaporMassFlux,r_MeltEnergy,Outflow)

      read (1000,*) OPTIONS%dtveg

      print*,'rddata:  Done initializing simulation'

      !TEMPORARY
      !SOIL
      GRID%SOIL%bcbeta = SOIL_PARAM%bcbeta
      GRID%SOIL%psic = SOIL_PARAM%psic
      GRID%SOIL%thetas = SOIL_PARAM%thetas
      GRID%SOIL%thetar = SOIL_PARAM%thetar
      GRID%SOIL%xk0 = SOIL_PARAM%xk0
      GRID%SOIL%zdeep = SOIL_PARAM%zdeep
      GRID%SOIL%tdeep = SOIL_PARAM%tdeep
      GRID%SOIL%zmid = SOIL_PARAM%zmid
      GRID%SOIL%tmid0 = SOIL_PARAM%tmid0
      GRID%SOIL%rocpsoil = SOIL_PARAM%rocpsoil
      GRID%SOIL%quartz = SOIL_PARAM%quartz
      GRID%SOIL%srespar1 = SOIL_PARAM%srespar1
      GRID%SOIL%srespar2 = SOIL_PARAM%srespar2
      GRID%SOIL%srespar3 = SOIL_PARAM%srespar3
      GRID%SOIL%a_ice = SOIL_PARAM%a_ice
      GRID%SOIL%b_ice = SOIL_PARAM%b_ice
      GRID%SOIL%bulk_dens = SOIL_PARAM%bulk_dens
      GRID%SOIL%amp = SOIL_PARAM%amp
      GRID%SOIL%phase = SOIL_PARAM%phase
      GRID%SOIL%shift = SOIL_PARAM%shift
      GRID%SOIL%bcgamm = SOIL_PARAM%bcgamm
      GRID%SOIL%par = SOIL_PARAM%par
      GRID%SOIL%corr = SOIL_PARAM%corr
      CAT%psicav = SOIL_PARAM%psicav
      GLOBAL%smpet0 = SOIL_PARAM%smpet0

      !Water Balance Variables
      GRID%VARS%rzsm1 = SOIL_MOISTURE%rzsm1
      GRID%VARS%tzsm1 = SOIL_MOISTURE%tzsm1
      GRID%VARS%rzsm1_f = SOIL_MOISTURE%rzsm1_f
      GRID%VARS%tzsm1_f = SOIL_MOISTURE%tzsm1_f
      GRID%VARS%rzdthetaidt = SOIL_MOISTURE%rzdthetaidt
      GRID%VARS%tzdthetaidt = SOIL_MOISTURE%tzdthetaidt
      
      return

      end subroutine


! ====================================================================
!
!			subroutine rdveg_update
!
! ====================================================================
!
! Subroutine to update the vegetation parameters.
!
! ====================================================================

      subroutine rdveg_update (&
       OPTIONS,ntdveg,GRID)

      implicit none
      include "wgtpar.h"
      include "help/rdveg_update.h"
      type VegDataTemplate
        real*8,dimension(:,:),allocatable :: xlai
        real*8,dimension(:,:),allocatable :: albd
      end type VegDataTemplate
      type (VegDataTemplate) VegData
      integer :: dvegnvars,ipos,jpos
      integer :: ntdveg
      real,dimension(:,:,:),allocatable :: TempArray
      type (OPTIONS_template) :: OPTIONS
      type (GRID_template),dimension(:),allocatable :: GRID
      dvegnvars = 2

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      ntdveg = ntdveg + 1
      allocate(TempArray(OPTIONS%ncol,OPTIONS%nrow,dvegnvars))
      allocate(VegData%xlai(OPTIONS%ncol,OPTIONS%nrow))
      allocate(VegData%albd(OPTIONS%ncol,OPTIONS%nrow))
      read(1003,rec=ntdveg)TempArray(:,:,:)
      VegData%xlai(:,:) = dble(TempArray(:,:,1))
      VegData%albd(:,:) = dble(TempArray(:,:,2))

! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

        do kk=1,OPTIONS%nlandc

                !Map the kk position to the i,j position
                if(mod(kk,OPTIONS%nrow) .ne. 0)then
                        ipos = kk/OPTIONS%nrow+1
                        jpos = mod(kk,OPTIONS%nrow)
                else
                        ipos = kk/OPTIONS%nrow
                        jpos = OPTIONS%nrow
                endif

                GRID(kk)%VEG%xlai = VegData%xlai(ipos,jpos) !dveg
                GRID(kk)%VEG%albd = VegData%albd(ipos,jpos) !dveg
                GRID(kk)%VEG%tcbeta = exp(-0.5*GRID(kk)%VEG%xlai)
                GRID(kk)%VEG%xlai_wsc = VegData%xlai(ipos,jpos)
                GRID(kk)%VEG%albw = VegData%albd(ipos,jpos)

        enddo


! ====================================================================
! Calculate parameters for each land cover type.
! ====================================================================

      do 300 kk=1,OPTIONS%nlandc

! --------------------------------------------------------------------&
! If not bare soil then calculate the canopy resistance.
! --------------------------------------------------------------------&

         if (OPTIONS%ivgtyp(kk).ne.0) then

            GRID(kk)%VEG%rescan = GRID(kk)%VEG%rsmin/GRID(kk)%VEG%xlai

! --------------------------------------------------------------------&
! If bare soil then set canopy resistance to zero.
! --------------------------------------------------------------------&

         else

            GRID(kk)%VEG%rescan = 0.

         endif

! --------------------------------------------------------------------&
! Calculate canopy storage capacity and initial canopy storage.
! --------------------------------------------------------------------&

         GRID(kk)%VEG%wsc = 0.0002*GRID(kk)%VEG%xlai_wsc

300   continue

      return

      end subroutine rdveg_update

! ====================================================================
!
!			subroutine rdtpmd
!
! ====================================================================
!
! Subroutine to read in and initialize topmodel parameters and the
!  soils-topographi! index map.
!
! ====================================================================

      subroutine rdtpmd(iopbf,iopwt0,ncatch,nrow,ncol,pixsiz,ipixnum,iprn,&
       ixpix,iypix,npix,q0,ff,qb0,dd,xlength,basink,xlamda,icatch,area,&
       atanb,dtil,zbar1,&
       iwel,wslp,&
       lat_deg,lat_min,lng_deg,lng_min,lng_mer,rlatitude,rlongitude,rlng_merid,&
       GRID,CAT)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/rdtpmd.h"
      type (GRID_template),dimension(:),allocatable :: GRID
      type (CATCHMENT_template),dimension(:),allocatable :: CAT

! ====================================================================
! Read the option for baseflow calculation, the option .or.&
! initial water table entry, the number of basin in the area
! of interest, and the dimensions of the soils-topographi!
! index map.
! ====================================================================

      read(1000,*) iopbf
      read(1000,*) iopwt0
      read(1000,*) ncatch

      print*,'rdtpmd:  Number of Catchments:  ',ncatch

      read(1000,*) nrow
      read(1000,*) ncol
      read(1000,*) pixsiz

      print*,'Rows: ',nrow
      print*,'Columns: ',ncol
      print*,'Pixel Size: ',pixsiz

!#####################################################################
! Allocate memory
!#####################################################################

allocate(CAT(ncatch))
allocate(GRID(nrow*ncol))

! ====================================================================
! Read in the soils-topographi! index map and set up the 
! transformation from pixel number to row/column.
! ====================================================================

      call rdatb(atb,nrow,ncol,ipixnum,ixpix,iypix,npix)

      print*,'rdtpmd:  Read topographi! index image '

! ====================================================================
! Read in the catchment look-up table - read different values based
! on what is necessary for baseflow and initial condition calculations.
! ====================================================================

      if (iopwt0.eq.1) then

         do 100 kk=1,ncatch

            read(71,*) jj,q0(kk),ff(kk),qb0(kk),dd(kk),&
                       xlength(kk),basink(kk) 

100      continue

      else if (iopbf.eq.1) then

         do 200 kk=1,ncatch

            read(71,*) jj,q0(kk),ff(kk),zbar0(kk),dd(kk),&
                       xlength(kk),basink(kk) 

200      continue

      else 

         do 300 kk=1,ncatch

            read(71,*) jj,q0(kk),ff(kk),zbar0(kk)

300      continue

      endif

      print*,'rdtpmd:  Read catchment lookup table'

! ====================================================================
! Read image of transmissivities for use in calculating the 
! soils-topographi! index.
! ====================================================================

      call rdimgr(ti,8,nrow,ncol,ipixnum)
! NWC 11/06/11 Read in K0 instead of T0, then divide by ff (exponential decay parameter) to obtain ti (Only valid for 1 catchment)
		do kk=1,nrow*ncol
			ti(kk) = ti(kk)/ff(1);
		enddo
      print*,'rdtpmd:  Read transmissivity image '

! ====================================================================
! Read the catchment image.
! ====================================================================

      call rdimgi(icatch,10,nrow,ncol,ipixnum)

      print*,'rdtpmd:  Read catchment image'

! --------------------------------------------------------------------&
! Check the dimensioning of the number of catchments.
! --------------------------------------------------------------------&

      do kk=1,npix

         if (icatch(kk).gt.MAX_CAT) then

            write (*,*) 'rdtpmd : highest catchment number '
            write (*,*) 'higher than MAX_CAT ',icatch(kk),MAX_CAT
            stop

         endif

      enddo

! ====================================================================
! Calculate the average topographi! index value, the average of 
! the natural log of transmissivity and area for each basin.
! ====================================================================

      do 400 kk=1,ncatch

         sumatb(kk) = 0.0
         sumlti(kk) = 0.0
         icount(kk) = 0

400   continue

      do 500 kk=1,npix

         sumatb(icatch(kk)) = sumatb(icatch(kk)) + atb(kk)
         sumlti(icatch(kk)) = sumlti(icatch(kk)) + dlog(ti(kk))
         icount(icatch(kk)) = icount(icatch(kk)) + 1

500   continue

      do 600 kk=1,ncatch

         if (icount(kk).eq.0) then

            xlamda(kk) = 0
            lte(kk) = 0

         else

            xlamda(kk) = sumatb(kk)/icount(kk)
            lte(kk) = sumlti(kk)/icount(kk)

         endif

         area(kk) = icount(kk)*pixsiz*pixsiz

600   continue

      print*,'Area',area(1)
      print*,'ln Te',lte(1)
      print*,'Lambda',xlamda(1)

! ====================================================================
! Calculate soils-topographi! index for each pixel.
! ====================================================================

      do 50 kk=1,npix

         atanb(kk) = atb(kk) + lte(icatch(kk)) - dlog(ti(kk))

50    continue

! ====================================================================
! Calculate the initial water table depth for each catchment.
! ====================================================================

      if ((iopwt0.eq.1).or.(iopbf.eq.1)) then

         do 700 kk=1,ncatch

            dtil(kk) = sqrt(q0(kk)/(3.45*basink(kk)*dd(kk)*xlength(kk)))

            if (iopwt0.eq.1) then

               hbar0 = sqrt(qb0(kk)/(5.772*basink(kk)*dd(kk)*xlength(kk)))
               zbar0(kk) = dtil(kk) - hbar0

            endif

700      continue

      endif

! ====================================================================
! Set initial average water table depth to the water table depth 
! at the end of the previous time step since program updates 
! this depth.
! ====================================================================

      do 800 kk=1,ncatch

         zbar1(kk) = zbar0(kk)

800   continue

! ====================================================================
! Output the soils-topographi! index image.
! ====================================================================

      if (iprn(7).eq.1) then

         call wrimgr(atanb,7,1.,nrow,ncol,ipixnum)
         close(7)

      endif

! ====================================================================
! Print table of results for each catchment.
! ====================================================================

      if (iprn(75).eq.1) then

         do 900 kk=1,ncatch

            write(75,1000) kk,area(kk)/1000000.,xlamda(kk),lte(kk),zbar0(kk),&
            ff(kk),q0(kk)

900      continue

         close(75)

      endif

! ====================================================================
! Format statement.
! ====================================================================

1000  format(i5,f15.3,5f10.3)

      return

      end subroutine rdtpmd

END MODULE MODULE_IO
