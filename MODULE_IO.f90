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

!##################################################################
!Subroutine to handle opening and closing the input files
!##################################################################

subroutine FILE_OPEN()

implicit none

character(len=200) :: filename
integer :: soilnvars,vegnvars,dvegnvars,nforcingvars,noutvars
integer :: nrow,ncol
soilnvars = 23
vegnvars = 20
noutvars = 1
dvegnvars = 2
nforcingvars = 7;
nrow = 60
ncol = 66

!Global Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GLOBAL_PARAMETER.txt"
open(1000,file=trim(filename))

!Soil Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/soil.bin"
open(1001,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*soilnvars*4)

!Vegetation Static Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/veg.bin"
open(1002,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*vegnvars*4)

!Vegetation Dynamic Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/dveg.bin"
open(1003,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*dvegnvars*4)

!Forcing Data Set
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/FORCING_01202005_10202005.bin"
open(1004,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*nforcingvars*4)

!Output Data Set
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/OUTPUT.bin"
open(1005,file=trim(filename),status='unknown',access='direct',form='unformatted',recl=ncol*nrow*noutvars*4)

!Regional Actual Energy Fluxes (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/EF_fruit.txt"
open(2091,file=trim(filename))

!Regional Canopy water balance (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/CWB_fruit.txt"
open(2092,file=trim(filename))

end subroutine FILE_OPEN

subroutine FILE_CLOSE()

implicit none
integer :: FILEUNIT

!Global Parameter File
close(1000)

!Soil Parameter File
close(1001)

!Vegetation Static Parameter File
close(1002)

!Vegetation Dynamic Parameter File
close(1003)

!Forcing Data Set
close(1004)

!Output Data Set
close(1005)

!Regional Actual Energy Fluxes (VALIDATION FILE)
close(2091)

!Regional Canopy water balance (VALIDATION FILE)
close(2092)

end subroutine FILE_CLOSE

! ====================================================================
! 
!			subroutine imgctl
!
! ====================================================================
!
! Subroutine to control output of time varying images by deciding
!   what images to print each time step, openning the output
!   file and writing the image
!
! ====================================================================

      subroutine imgctl(tkact,gact,hact,xleact,rnact,etpix,&
       runtot,xinact,irntyp,zw,tzsm1,rzsm1,tkpet,wcip1,gpet,hpet,xlepet,rnpet,&
       r_mossm,pptms,ievcon,ipixnum,icurser,nseries,ioutst,iouten,iprn,ioutsp,&
       pptsumrg,i,ncol,nrow,fnimg,img_opt,iyear,iday,ihour,Swq,Swq_us)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/imgctl.h"

! ====================================================================
! Loop through each of the possible output image unit numbers.
! ====================================================================

      do 100 iu=35,55

! ====================================================================
! Set the default value of the flag to signify image output
! for this time step to not print (zero).
! ====================================================================

         imgprn(iu) = 0

! ====================================================================
! Check if any images of this type are to be printed.
! ====================================================================

         if (iprn(iu).eq.1) then

! ====================================================================
! Check if this time step has image output.  If so, then 
! open the image file and reset the flag to print the 
! image.
! ====================================================================

            if ( (iouten(iu,icurser(iu)).lt.i).and. &
                 (nseries(iu).gt.icurser(iu)) ) icurser(iu) = icurser(iu) + 1 

            if ( (mod(i-ioutst(iu,icurser(iu)),ioutsp(iu,icurser(iu))).eq.0)&
               .and.&
                 (i.ge.ioutst(iu,icurser(iu)))&
               .and.&
                 (i.le.iouten(iu,icurser(iu))) ) then

                   call imgopn(fnimg(iu),i,iu,img_opt,iyear,iday,ihour)
                   imgprn(iu) = 1

            endif

         endif

100   continue

      do 200 iu=151,199

! ====================================================================
! Set the default value of the flag to signify image output
! for this time step to not print (zero).
! ====================================================================

         imgprn(iu) = 0

! ====================================================================
! Check if any images of this type are to be printed.
! ====================================================================

         if (iprn(iu).eq.1) then

! ====================================================================
! Check if this time step has image output.  If so, then 
! open the image file and reset the flag to print the 
! image.
! ====================================================================

            if ( (iouten(iu,icurser(iu)).lt.i).and. &
                 (nseries(iu).gt.icurser(iu)) ) icurser(iu) = icurser(iu) + 1 

            if ( (mod(i-ioutst(iu,icurser(iu)),ioutsp(iu,icurser(iu))).eq.0)&
               .and.&
                 (i.ge.ioutst(iu,icurser(iu)))&
               .and.&
                 (i.le.iouten(iu,icurser(iu))) ) then

                   call imgopn(fnimg(iu),i,iu,img_opt,iyear,iday,ihour)
                   imgprn(iu) = 1

            endif

         endif

200   continue

! ====================================================================
! Print the requested images.
! ====================================================================

      if (imgprn(35).eq.1) then

         if (pptsumrg.gt.(0.d0))&
             call wrimgr(pptms,35,3600000.,nrow,ncol,ipixnum)

      endif

      !if (imgprn(36).eq.1) call wrimgr(r_mossm,36,100.,nrow,ncol,ipixnum)
      if (imgprn(37).eq.1) call wrimgr(rnpet,37,1.,nrow,ncol,ipixnum)
      if (imgprn(38).eq.1) call wrimgr(xlepet,38,1.,nrow,ncol,ipixnum)
      if (imgprn(39).eq.1) call wrimgr(hpet,39,1.,nrow,ncol,ipixnum)
      if (imgprn(40).eq.1) call wrimgr(gpet,40,1.,nrow,ncol,ipixnum)
      if (imgprn(41).eq.1) call wrimgr(tkpet,41,1.,nrow,ncol,ipixnum)
      if (imgprn(42).eq.1) call wrimgr(wcip1,42,1000000.,nrow,ncol,ipixnum)
      if (imgprn(43).eq.1) call wrimgr(rzsm1,43,100.,nrow,ncol,ipixnum)
      if (imgprn(43).eq.1) call WRITE_BINARY(rzsm1,100.,nrow,ncol,ipixnum,i)
      if (imgprn(44).eq.1) call wrimgr(tzsm1,44,100.,nrow,ncol,ipixnum)
      if (imgprn(45).eq.1) call wrimgr(zw,45,100.,nrow,ncol,ipixnum)

! ====================================================================
! Close any open image files.
! ====================================================================

      do 203 iu=35,45

         if (imgprn(iu).eq.1) close(iu)

203   continue

      if (imgprn(46).eq.1) call wrimgr(xinact,46,3600000.,nrow,ncol,ipixnum)
      if (imgprn(47).eq.1) call wrimgr(runtot,47,3600000.,nrow,ncol,ipixnum)
      if (imgprn(48).eq.1) call wrimgi(irntyp,48,1,nrow,ncol,ipixnum)
      if (imgprn(49).eq.1) call wrimgr(etpix,49,3600000.,nrow,ncol,ipixnum)
      if (imgprn(50).eq.1) call wrimgi(ievcon,50,1,nrow,ncol,ipixnum)
      if (imgprn(51).eq.1) call wrimgr(rnact,51,1.,nrow,ncol,ipixnum)
      if (imgprn(52).eq.1) call wrimgr(xleact,52,1.,nrow,ncol,ipixnum)
      if (imgprn(53).eq.1) call wrimgr(hact,53,1.,nrow,ncol,ipixnum)
      if (imgprn(54).eq.1) call wrimgr(gact,54,1.,nrow,ncol,ipixnum)
      if (imgprn(55).eq.1) call wrimgr(tkact,55,1.,nrow,ncol,ipixnum)
      if (imgprn(151).eq.1) call wrimgr(Swq,151,1.,nrow,ncol,ipixnum)
      !if (imgprn(152).eq.1) call wrimgr(Swq_us,152,1.,nrow,ncol,ipixnum)

! ====================================================================
! Close any open image files.
! ====================================================================

      do 201 iu=46,55

         if (imgprn(iu).eq.1) close(iu)

201   continue

      do 202 iu=151,199

         if (imgprn(iu).eq.1) close(iu)

202   continue

      return

      end subroutine imgctl

! ====================================================================
!
!			subroutine rdveg
!
! ====================================================================
!
! Subroutine to read and initiailize simulation constant vegetation
!   and land cover parameters.
!
! ====================================================================

      subroutine rdveg(npix,nrow,ncol,ilandc,ipixnum,nlandc,iopveg,ivgtyp,iprn,&
       xlai,xlai_wsc,albd,albw,emiss,za,zww,z0m,z0h,zpd,rsmin,rsmax,Rpl,&
       f3vpdpar,f4temppar,trefk,tcbeta,extinct,canclos,Tslope1,Tint1,&
       Tslope2,Tint2,Twslope1,Twint1,Twslope2,Twint2,Tsep,Twsep,&
       eps,&
       rtact,rtdens,rtres,psicri,&
       rescan,respla,wsc,wcip1,&
       pixsiz,area,fbs,fbsrg,ncatch)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      include "help/rdveg.h"
      type VegDataTemplate
        real*8,dimension(:,:),allocatable :: ivgtyp
        real*8,dimension(:,:),allocatable :: xlai
        real*8,dimension(:,:),allocatable :: xlai_wsc
        real*8,dimension(:,:),allocatable :: albd
        real*8,dimension(:,:),allocatable :: albw
        real*8,dimension(:,:),allocatable :: emiss
        real*8,dimension(:,:),allocatable :: za
        real*8,dimension(:,:),allocatable :: zww
        real*8,dimension(:,:),allocatable :: z0m
        real*8,dimension(:,:),allocatable :: z0h
        real*8,dimension(:,:),allocatable :: zpd
        real*8,dimension(:,:),allocatable :: rsmin
        real*8,dimension(:,:),allocatable :: rsmax
        real*8,dimension(:,:),allocatable :: Rpl
        real*8,dimension(:,:),allocatable :: f3vpdpar
        real*8,dimension(:,:),allocatable :: f4temppar
        real*8,dimension(:,:),allocatable :: trefk
        real*8,dimension(:,:),allocatable :: tcbeta
        real*8,dimension(:,:),allocatable :: extinct
        real*8,dimension(:,:),allocatable :: canclos
      end type VegDataTemplate
      type (VegDataTemplate) VegData
      character(len=200) :: filename
      integer :: vegnvars,dvegnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      vegnvars = 20
      dvegnvars = 2
      allocate(TempArray(ncol,nrow,vegnvars))

! ====================================================================
! Read the image with the land cover clasifications.
! ====================================================================

      call rdimgi(ilandc,11,nrow,ncol,ipixnum)

! ====================================================================
! Read spatially constant vegetation parameters.
! ====================================================================

      nlandc = nrow*ncol
      iopveg = 0
      read(1000,*) iopwc0

      print*,"rdveg:  Read spatially constant veg pars"

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      print*,"rdveg:  Read lookup table"

! ====================================================================
! Read the binary vegetation binary file
! ====================================================================

      allocate(VegData%ivgtyp(ncol,nrow))
      allocate(VegData%xlai(ncol,nrow))
      allocate(VegData%xlai_wsc(ncol,nrow))
      allocate(VegData%albd(ncol,nrow))
      allocate(VegData%albw(ncol,nrow))
      allocate(VegData%emiss(ncol,nrow))
      allocate(VegData%za(ncol,nrow))
      allocate(VegData%zww(ncol,nrow))
      allocate(VegData%z0m(ncol,nrow))
      allocate(VegData%z0h(ncol,nrow))
      allocate(VegData%zpd(ncol,nrow))
      allocate(VegData%rsmin(ncol,nrow))
      allocate(VegData%rsmax(ncol,nrow))
      allocate(VegData%Rpl(ncol,nrow))
      allocate(VegData%f3vpdpar(ncol,nrow))
      allocate(VegData%f4temppar(ncol,nrow))
      allocate(VegData%trefk(ncol,nrow))
      allocate(VegData%tcbeta(ncol,nrow))
      allocate(VegData%extinct(ncol,nrow))
      allocate(VegData%canclos(ncol,nrow))
        
      print*,"rdveg:  Reading in all the vegetation properties at once"

      read(1002,rec=1)TempArray(:,:,:)

      VegData%ivgtyp(:,:) = dble(TempArray(:,:,1))
      VegData%xlai(:,:) = dble(TempArray(:,:,2))
      VegData%xlai_wsc(:,:) = dble(TempArray(:,:,3))
      VegData%albd(:,:) = dble(TempArray(:,:,4))
      VegData%albw(:,:) = dble(TempArray(:,:,5))
      VegData%emiss(:,:) = dble(TempArray(:,:,6))
      VegData%za(:,:) = dble(TempArray(:,:,7))
      VegData%zww(:,:) = dble(TempArray(:,:,8))
      VegData%z0m(:,:) = dble(TempArray(:,:,9))
      VegData%z0h(:,:) = dble(TempArray(:,:,10))
      VegData%zpd(:,:) = dble(TempArray(:,:,11))
      VegData%rsmin(:,:) = dble(TempArray(:,:,12))
      VegData%rsmax(:,:) = dble(TempArray(:,:,13))
      VegData%Rpl(:,:) = dble(TempArray(:,:,14))
      VegData%f3vpdpar(:,:) = dble(TempArray(:,:,15))
      VegData%f4temppar(:,:) = dble(TempArray(:,:,16))
      VegData%trefk(:,:) = dble(TempArray(:,:,17))
      VegData%tcbeta(:,:) = dble(TempArray(:,:,18))
      VegData%extinct(:,:) = dble(TempArray(:,:,19))
      VegData%canclos(:,:) = dble(TempArray(:,:,20))

! ####################################################################
! Read the vegetation dynamic parameter file
! ####################################################################

      deallocate(TempArray)
      allocate(TempArray(ncol,nrow,dvegnvars))

      print*,"rdveg:  Reading in the dynamic vegetation properties"

      read(1003,rec=1)TempArray(:,:,:)

      VegData%xlai(:,:) = dble(TempArray(:,:,1))
      VegData%albd(:,:) = dble(TempArray(:,:,2))


! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

        do kk=1,nlandc

                !Map the kk position to the i,j position
                if(mod(kk,nrow) .ne. 0)then
                        ipos = kk/nrow+1
                        jpos = mod(kk,nrow)
                else
                        ipos = kk/nrow
                        jpos = nrow
                endif
                ivgtyp(kk) = VegData%ivgtyp(ipos,jpos)
                emiss(kk) = VegData%emiss(ipos,jpos)
                za(kk) = VegData%za(ipos,jpos)
                zww(kk) = VegData%zww(ipos,jpos)
                z0m(kk) = VegData%z0m(ipos,jpos)
                z0h(kk) = VegData%z0h(ipos,jpos)
                zpd(kk) = VegData%zpd(ipos,jpos)
                rsmin(kk) = VegData%rsmin(ipos,jpos)
                rsmax(kk) = VegData%rsmax(ipos,jpos)
                Rpl(kk) = VegData%Rpl(ipos,jpos)
                f3vpdpar(kk) = VegData%f3vpdpar(ipos,jpos)
                f4temppar(kk) = VegData%f4temppar(ipos,jpos)
                trefk(kk) = VegData%trefk(ipos,jpos)
                xlai(kk) = VegData%xlai(ipos,jpos) !dveg
                albd(kk) = VegData%albd(ipos,jpos) !dveg
                tcbeta(kk) = exp(-0.5*xlai(kk)) !dveg
                xlai_wsc(kk) = xlai(kk)
                albw(kk) = albd(kk) !Move to its own file
                extinct(kk) = 0.00!VegData%extinct(ipos,jpos)
                canclos(kk) = 1.00!VegData%canclos(ipos,jpos)
                Tslope1(kk) = 0.00!VegData%Tslope1(ipos,jpos)
                Tint1(kk) = 0.00!VegData%Tint1(ipos,jpos)
                Tslope2(kk) = 0.00!VegData%Tslope2(ipos,jpos)
                Tint2(kk) = 0.00!VegData%Tint2(ipos,jpos)
                Twslope1(kk) = 0.00!VegData%Twslope1(ipos,jpos)
                Twint1(kk) = 0.00!VegData%Twint1(ipos,jpos)
                Twslope2(kk) = 0.00!VegData%Twslope2(ipos,jpos)
                Twint2(kk) = 0.00!VegData%Twint2(ipos,jpos)
                Tsep(kk) = 0.00!VegData%Tsep(ipos,jpos)
                Twsep(kk) = 0.00!VegData%Twsep(ipos,jpos)

        enddo

! ====================================================================
! Calculate parameters for each land cover type.
! ====================================================================

      do kk=1,nlandc

! --------------------------------------------------------------------&
! If not bare soil then calculate the canopy resistance.
! --------------------------------------------------------------------&

         if (ivgtyp(kk).ne.0) then

            rescan(kk) = rsmin(kk)/xlai(kk)

! --------------------------------------------------------------------&
! If bare soil then set canopy resistance to zero.
! --------------------------------------------------------------------&

         else

            rescan(kk) = 0.

         endif

! --------------------------------------------------------------------&
! Calculate canopy storage capacity and initial canopy storage.
! --------------------------------------------------------------------&

         wsc(kk) = 0.0002*xlai_wsc(kk)

      enddo

      print*,"rdveg:  Set minimum st. resist. and and wet can.stor.cap."

! ====================================================================
! Read the initial canopy storage from an image or as a constant
! as requested.
! ====================================================================

      if (iopwc0.eq.0) then

         read(1000,*) wc0

         do kk=1,npix

            wcip1(kk) = wc0

         enddo  

         print*,"rdveg:  Read initial wet canopy storages"

      else 
		 
         call rdimgr(wcip1,13,nrow,ncol,ipixnum)

      endif

! ====================================================================
! Calculate the fraction of different cover types in each catchment
! and in total area.  First array index for icount and frcov is
! the land cover type, second array index is the catchment number.
! Catchment ncatch+1 is total area.
! ====================================================================

      do 550 kk=1,nlandc

         do 540 jj=1,ncatch

            frcov(kk,jj) = pixsiz*pixsiz/area(jj)

540      continue

         frcov(kk,ncatch+1) = 1/real(npix)

550   continue

      print*,"rdveg:  Calculated fractional covers for vegetation"

! ====================================================================
! Find fraction of bare soil in each catchment.
! ====================================================================

      fbsrg = zero

      do 570 jj=1,ncatch+1

         fbs(jj)  = zero

         do 560 kk=1,nlandc

            if (ivgtyp(kk).eq.0) then
               
               if (jj .eq. ncatch+1) then
               
                  fbsrg = fbsrg + frcov(kk,jj)
        
               else

                  fbs(jj) = fbs(jj) + frcov(kk,jj)

               endif

            endif

560      continue

570   continue

      print*,"rdveg:  Calculated fractional covers for bare soil"

! ====================================================================
! Print land cover table summary.
! ====================================================================

      if (iprn(76).eq.1) then

         do 600 kk = 1,nlandc

            write(76,1000) kk,xlai(kk),rescan(kk),wsc(kk)

600      continue

      endif

! ====================================================================
! Print land cover fractions in each catchment.
! ====================================================================

      if (iprn(77).eq.1) then

         do 700 jj=1,ncatch

            write (77,1100) jj,area(jj)/1000000.

            do 650 kk=1,nlandc

               write(77,1110) kk,frcov(kk,jj)*100

650         continue

            write (77,*)

700      continue

         write(77,1120)

         do 800 kk=1,nlandc

            write(77,1110) kk,frcov(kk,ncatch+1)*100

800      continue

      endif        

! ====================================================================
! Format statements.
! ====================================================================

1000  format (i5,2f10.3,f11.7)
1100  format ('Catchment Number',i4,' (Area = ',f8.2,' km^2):')
1110  format ('   Land Cover Type',i4,':',f7.2,'%')
1120  format ('Total Area:')

      return

      end subroutine rdveg

! ====================================================================
!
!			subroutine hdrprnt
!
! ====================================================================
!
! Subroutine to print the output file headings.
!
! ====================================================================

      subroutine hdprnt(iprn)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      integer iprn(MAX_FIL)

! ====================================================================
! Write output file headers.
! ====================================================================

      if (iprn(75).eq.1) then

         write(75,*)'  Catchment Table'
         write(75,*)
         write(75,*)&
       '                     Average              Initial            Qo'
         write(75,*)'Basin              Topographic              W.T.'
         write(75,*)' No.       Area       Index      ln Te     Depth'
         write(75,*)&
           '(km^2)                            (m)      (m-1) (m^3/s)'
      endif


      if (iprn(76).eq.1) then

         write(76,*)'  Land Cover Table'
         write(76,*)
         write(76,*)'Land    Leaf               Canopy ' 
         write(76,*)'Cover   Area     Canopy    Storage'
         write(76,*)' No.    Index    Resist   Capacity.'
         write(76,*)'                 (s/m)       (m)'

      endif

      if (iprn(77).eq.1) then

         write(77,*)'  Land Cover Fractions by Catchment'
         write(77,*)

      endif

      if (iprn(78).eq.1) then

         write(78,*)'  Soil Type Table'
         write(78,*)
         write(78,*)&
        'Soil Pore Size           Saturated  Residual Saturated' 
         write(78,*)&
        'Type   Dist.    Bubbling    Soil      Soil   Hydraulic'
         write(78,*)&
        ' No.   Index    Pressure  Moisture  Moisture    Cond'
         write(78,*)'                  (cm)                         (mm/h)'

      endif      

      if (iprn(79).eq.1) then

         write(79,*)'  Soil Type Fractions by Catchment'
         write(79,*)

      endif

      if (iprn(80).eq.1) then

         write(80,*)' Catchment Evaporation Results'
         write(80,*)
         write(80,*)&
                '                       Bare       Dry       Wet    Frac   Frac'
         write(80,*)&
                'Time Catch   Total     Soil     Canopy    Canopy   Wet    Bare'
         write(80,*)&
                'Step  No.    Evap      Evap      Evap      Evap   Canopy  Soil'
         write(80,*)'            (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(81).eq.1) then

         write(81,*)' Catchment Precipitation/Infiltration Results'
         write(81,*)
         write(81,*)'                                                                Sat       Inf'
         write(81,*)'Time Catch   Total      Net                          Total    Excess    Excess'
         write(81,*)'Step  No.   Precip    Precip    Condens   Runoff    Infilt    Infilt    Infilt'
         write(81,*)'            (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(82).eq.1) then
! NWC 13/06/11
!         write(82,*)' Catchment Water Table Balance Results'
!         write(82,*)
!         write(82,*)'            New    Old                                            RZ     TZ'
!         write(82,*)'Time Catch  Ave    Ave   Capillary Drainage    Evap              Avail  Avail'
!         write(82,*)'Step  No.   W.T.   W.T.    Rise     to W.T.  from WT   Baseflow  Pore   Pore'
!         write(82,*)
!                '            (m)    (m)    (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(90).eq.1) then

         write(90,*)' Regional Energy Fluxes at PET'
         write(90,*)
         write(90,*)' Time    Net     Latent   Sensible   Ground     Energy   G + DS   Surface   Surface   Surface'
         write(90,*)' Step Radiation   Heat      Heat      Heat     Balance              Temp     Temp      Temp'
         write(90,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)    (K)       (K)       (K)'

      endif
         
      if (iprn(91).eq.1) then

         write(91,*)' Regional Actual Energy Fluxes'
         write(91,*)
         write(91,*)' Time    Net     Latent   Sensible   Ground     Energy   G + DS   Surface   Surface   Surface'
         write(91,*)' Step Radiation   Heat      Heat      Heat     Balance              Temp     Temp      Temp'
         write(91,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)    (K)       (K)       (K)'
!cdp      write(91,*)' Time    Net     Latent   Sensible   Ground     Energy   Surface'
!cdp      write(91,*)' Step Radiation   Heat      Heat      Heat     Balance    Temp'
!cdp      write(91,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)'

      endif

      if (iprn(92).eq.1) then

         write(92,*)' Regional Canopy Water Balance'
         write(92,*)
         write(92,*)&
              '         New       Old                           Wet    Fraction'
         write(92,*)' Time   Canopy    Canopy   Total      Net      Canopy      Wet    Change in  Sum of     Water'
         write(92,*)' Step  Storage   Storage  Rainfall  Rainfall    Evap     Canopy    Storage   Fluxes    Balance'
         write(92,*)'         (mm)      (mm)    (mm/h)    (mm/h)    (mm/h)                (m)       (m)       (m)'

      endif
         
      if (iprn(93).eq.1) then

         write(93,*)' Regional Precipitation/Infiltration/Runoff'
         write(93,*)
         write(93,*)'                                                       Saturation  Infilt'
         write(93,*)' Time   Total      Net               Surface    Total    Excess    Excess'
         write(93,*)' Step Rainfall  Rainfall   Condens   Runoff    Infilt    Runoff    Runoff'
         write(93,*)'       (mm/h)    (mm/h)     (mm/h)   (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(94).eq.1) then

         write(94,*)' Regional Evapotranspiration Rates'
         write(94,*)
         write(94,*)'                  Bare       Dry       Wet     Frac   Frac'
         write(94,*)' Time   Total     Soil     Canopy    Canopy    Wet    Bare'
         write(94,*)' Step   Evap      Evap      Trans     Evap    Canopy  Soil'
         write(94,*)'       (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif
     
      if (iprn(95).eq.1) then

         write(95,*)' Regional Root and Transmission Zone Water Balance'
         write(95,*)
         write(95,*)&
              '            Root Zone                      Transmision Zone'
         write(95,*)&
               ' Time  Soil  Change in  Sum of         Soil  Change in   Sum of'
         write(95,*)&
               ' Step Moist   Storage   Fluxes        Moist   Storage    Fluxes'
         write(95,*)&
               '               (mm)      (mm)                   (mm)       (mm)'

      endif
         
      if (iprn(96).eq.1) then

! NWC 06/13/11
!         write(96,*)' Regional Water Table Balance and Vertical Fluxes'
!         write(96,*)
!         write(96,*)'        New       Old'
!         write(96,*)' Time   Avg       Avg     Capillary       Drainage       Evap
!                       Drainage      Drainage      Diffusion  Diffusion  '
!         write(96,*)' Step  Depth     Depth    Rise fr WT       to WT        from WT
!              Baseflow  fr RZ        fr TZ          to RZ       to TZ   '
!         write(96,*)'        (mm)      (mm)     (mm/h)         (mm/h)         (mm/h)   
!             (mm/h)     (mm/h)      (mm/h)         (mm/h)      (mm/h)  '

      endif

      if (iprn(97).eq.1) then

         write(97,*)' Regional Fractional Saturation States'
         write(97,*)
         write(97,*)&
              ' Time  Region 3         Region 2                        Region 1 '
         write(97,*)' Step                  Saturated   Unsat         Saturated  TZ Sat  RZ Sat   Unsat'

      endif

      if (iprn(98).eq.1) then

         write(98,*)' Regional Evapotranspiration Controls and'
         write(98,*)'    Infiltration Mechanisms'
         write(98,*)
         write(98,*)&
               ' Time    Total  Atmos  Atmos   Land         Net    Sat   Infilt'
         write(98,*)&
                ' Step     Evap   Sat   Unsat  Unsat      Rainfall Exces  Exces'
         write(98,*)&
                              '         (mm/h)                           (mm/h)'

      endif

      if (iprn(110).eq.1) then

         write (110,*) 'Results for the moss layer, actual fluxes '
         write (110,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin      Ttopsoil        Tmid   theta_mois'

      endif

      if (iprn(111).eq.1) then

         write (111,*) 'Results for the understory layer, actual fluxes'
         write (111,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin        Tmid'

      endif

      if (iprn(112).eq.1) then

         write (112,*) 'Results for the overstory layer, actual fluxes'
          write (112,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin        Tmid'

      endif

      return

      end subroutine hdprnt

! ====================================================================
!
!			subroutine inisim
!
! ====================================================================
!
! Subroutine to initialize simulation total water balance variables.
!
! ====================================================================

      subroutine inisim(iopsmini,nrow,ncol,ipixnum,ilandc,npix,inc_frozen,&
       istorm,intstm,istmst,intstp,istorm_moss,intstm_moss,istmst_moss,&
       intstp_moss,isoil,idifind,smpet0,r_mossmpet0,endstm,rzsm1,&
       tzsm1,r_mossm1,r_mossm,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,r_mossm1_u,&
       r_mossm_u,r_mossm1_f,r_mossm_f,rzdthetaidt,tzdthetaidt,zmoss,&
       r_moss_depth,thetas_moss,xintst,xintst_moss,cuminf,xk0,psic,thetas,&
       thetar,bcgamm,bcbeta,sorp,cc,dt,sesq,corr,par,&
       PackWater_us,SurfWater_us,Swq_us,VaporMassFlux_us,r_MeltEnergy_us,&
       Outflow_us,PackWater,SurfWater,Swq,VaporMassFlux,r_MeltEnergy,Outflow)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      include "help/inisim.h"

      if (iopsmini.eq.1) then

! ====================================================================
! If initial root and transmission zone soil moistures are    
! input as images then read the images and assign to the
! initial condition.
! ====================================================================

! --------------------------------------------------------------------
! Alert the user.
! --------------------------------------------------------------------

	 print *, "Reading initial soil moisture images"

         call rdimgr(rzsm1,14,nrow,ncol,ipixnum)
         call rdimgr(tzsm1,15,nrow,ncol,ipixnum)

	 close(14)
	 close(15)

      else

! ====================================================================
! If initial root zone is not entered, then set the root
! zone soil moisture for calculation of thermal conductivity.
! This is changed later (initsm) to a new initial condition      
! based on Brooks-Corey and local water table depth.
! ====================================================================

         do 50 kk=1,npix

            if (MOS_FLG.eq.1) then

               m_kk=kk
               v_kk=ilandc(kk)

            endif

            if (MOS_FLG.eq.0) then

               m_kk=1
               v_kk=1

            endif

            rzsm1(kk) = smpet0
            tzsm1(kk) = smpet0
            r_mossm1(m_kk) = r_mossmpet0(v_kk)
            r_mossm(m_kk) = r_mossmpet0(v_kk)
            rzsm1_u(kk) = smpet0
            tzsm1_u(kk) = smpet0
            rzsm1_f(kk) = 0.d0
            tzsm1_f(kk) = 0.d0
            r_mossm1_u(m_kk) = r_mossmpet0(v_kk)
            r_mossm_u(m_kk) = r_mossmpet0(v_kk)
            r_mossm1_f(m_kk) = 0.d0
            r_mossm_f(m_kk) = 0.d0
            rzdthetaidt(kk)=0.d0
            tzdthetaidt(kk)=0.d0

50       continue

      endif

      do kk=1,npix

         if (MOS_FLG.eq.1) then

            m_kk=kk
            v_kk=ilandc(kk)

         endif

         if (MOS_FLG.eq.0) then

            m_kk=1
            v_kk=1

         endif

         if (inc_frozen.eq.0) then

          if(thetas_moss(v_kk).eq.0.)then
            zmoss(m_kk)=0.
          else
            zmoss(m_kk)=r_moss_depth(v_kk)*r_mossm(m_kk)/&
                        thetas_moss(v_kk)
          endif
         else
          if(thetas_moss(v_kk).eq.0.)then
            zmoss(m_kk)=0.
          else

            zmoss(m_kk)=r_moss_depth(v_kk)*r_mossm_u(m_kk)/&
                        thetas_moss(v_kk)
          endif

         endif

      enddo

! ====================================================================
! Read data to tell how program will initialize the storm
! and interstorm event flags and times.
! ====================================================================

      read(1000,*) iopflg

      if (iopflg.eq.0) then

! --------------------------------------------------------------------
! If one event flag value is used then set all flags and times
! accordingly.
! --------------------------------------------------------------------

         read(1000,*) istflg

         if (istflg.eq.1) then

! ....................................................................
! If the event is a storm event.
! ....................................................................

            do 100 kk=1,npix

               if (MOS_FLG.eq.1) m_kk=kk
               if (MOS_FLG.eq.0) m_kk=1

               istorm(kk) = 1
               intstm(kk) = 0
               istmst(kk) = 0
               intstp(kk) = 0
               xintst(kk) = 0.0
               istorm_moss(m_kk) = 1
               intstm_moss(m_kk) = 0
               istmst_moss(m_kk) = 0
               intstp_moss(m_kk) = 0
               xintst_moss(m_kk) = 0.0

100         continue

         else

! ....................................................................
! If the event is an interstorm event.
! ....................................................................

            do 200 kk=1,npix

               if (MOS_FLG.eq.1) m_kk=kk
               if (MOS_FLG.eq.0) m_kk=1

               istorm(kk) = 0
               intstm(kk) = 1
               istmst(kk) = 0
               intstp(kk) = 0
               xintst(kk) = endstm
               istorm_moss(m_kk) = 0
               intstm_moss(m_kk) = 1
               istmst_moss(m_kk) = 0
               intstp_moss(m_kk) = 0
               xintst_moss(m_kk) = endstm

200         continue

         endif

      else

! --------------------------------------------------------------------
! If spatially varying data is used read image data.
! --------------------------------------------------------------------

         call rdimgi(istorm,16,nrow,ncol,ipixnum)
         call rdimgi(istep,17,nrow,ncol,ipixnum)
         call rdimgr(cumdep,18,nrow,ncol,ipixnum)
         call rdimgr(smbeg,19,nrow,ncol,ipixnum)

! --------------------------------------------------------------------
! Set event flags, time of events, cumulative values.
! --------------------------------------------------------------------

         do 300 kk=1,npix

! ....................................................................
! For pixels under storm event.
! ....................................................................

            if (istorm(kk).eq.1) then

               intstm(kk) = 0
               istmst(kk) = istep(kk)
               intstp(kk) = 0
               xintst(kk) = 0.0
               cuminf(kk) = cumdep(kk)

! ....................................................................
! Find philip's equation parameters.
! ....................................................................

               sorp(kk) = (((two*xk0(isoil(kk))*&
                          ((thetas(isoil(kk))-smbeg(kk))**two)&
                          *psic(isoil(kk)))/&
                          (thetas(isoil(kk))-thetar(isoil(kk))))*&
                          ((one/(bcgamm(isoil(kk))+&
                               0.5d0*bcbeta(isoil(kk))-one)) +&
                          ((thetas(isoil(kk))-thetar(isoil(kk))) /&
                           (thetas(isoil(kk))-smbeg(kk))))) ** 0.5d0
               deltrz = smbeg(kk) - thetar(isoil(kk))
               if (deltrz.le.zero) deltrz=zero
               cc(kk) = 0.5d0 *&
                       (one+((deltrz/&
                              (thetas(isoil(kk))-thetar(isoil(kk)))) **&
                             (bcgamm(isoil(kk))/bcbeta(isoil(kk)))))

! ....................................................................
! For pixels under interstorm event.
! ....................................................................

            else

               intstm(kk) = 1
               istmst(kk) = 0
               intstp(kk) = istep(kk)
               xintst(kk) = intstp(kk)*dt + endstm
               
               relrze = (smbeg(kk) - thetar(isoil(kk))) /&
                        (thetas(isoil(kk)) - thetar(isoil(kk)))

               if (relrze.le.zero) relrze=zero
               if (relrze.ge.one) relrze=one


            endif

300      continue

      endif

! ====================================================================
! Initialize snow pack variables.
! ====================================================================

      do kk=1,npix

         if (SNW_FLG.eq.1) s_kk=kk
         if (SNW_FLG.eq.0) s_kk=1
         if (SNOW_RUN.eq.1) sw_kk=kk
         if (SNOW_RUN.eq.0) sw_kk=1

         PackWater_us(s_kk)=zero
         SurfWater_us(s_kk)=zero
         Swq_us(s_kk)=zero
         VaporMassFlux_us(s_kk)=zero
         r_MeltEnergy_us(s_kk)=zero
         Outflow_us(s_kk)=zero
         PackWater(sw_kk)=zero
         SurfWater(sw_kk)=zero
         Swq(sw_kk)=zero
         VaporMassFlux(sw_kk)=zero
         r_MeltEnergy(sw_kk)=zero
         Outflow(sw_kk)=zero

      enddo

      return

      end subroutine inisim

! ====================================================================
!
!			subroutine filopn
!
! ====================================================================
!
! Subroutine to open input/output files and set up variable to 
! control the printing of output files.
!
! ====================================================================

      subroutine filopn(iprn,ioutst,ioutsp,iouten,nseries,icurser,fnimg)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      !include "sun_sgi.h"
      include "help/filopn.h"

! ====================================================================
! Initialize variable which controls output file printing.
! ====================================================================
     

      do 40 jj = 1,MAX_FIL

        iprn(jj) = 0

40    continue 
     
! ====================================================================
! Open the file which lists the input/output file names.  This
! file should be given as the second word in the command line.
! ====================================================================

       call getarg(1,topfil)
     
       open(unit=4,file=topfil)

! ====================================================================
! Read the filenames from file containing i/o filenames.
! ====================================================================

50    read(4,1000) iu,fname
    
      if (iu.gt.MAX_FIL) then

         write (*,*) 'filopn.f : File unit number is greater than '
         write (*,*) 'MAX_FIL ',iu,MAX_FIL
         stop

      endif

! ====================================================================
! If unit number is between 35 and 55, then read the 
! image output start time and spacing and the image
! name prefix for each series of image output.
! ====================================================================

      if ( ((iu.ge.35).and.(iu.le.55)).or.&
           ((iu.ge.151).and.(iu.le.199)) ) then

         iseries = 0
60       iseries = iseries + 1 

         if (iseries.gt.MAX_SER) then

            write (*,*) 'filopn.f : Number of input/output series'
            write (*,*) 'is greater MAX_SER ',iseries,MAX_SER
            stop

         endif
 
! --------------------------------------------------------------------
! As long as there are parts of the timeseries for which this variable
! is to be printed out keep on reading these parts in.
! --------------------------------------------------------------------

         read(4,*) ioutst(iu,iseries),ioutsp(iu,iseries),&
                   iouten(iu,iseries),icont

         if (icont.eq.1) goto 60

         nseries(iu) = iseries
         icurser(iu) = 1
         fnimg(iu) = fname

! ====================================================================
! If unit number is between 20 and 34, then set the
! input image prefix.
! ====================================================================

      else if ((iu.ge.20).and.(iu.le.34)) then

         fnimg(iu) = fname      

! ====================================================================
! If the unit number is between 7 and 19 then open file for
! binary image.
! ====================================================================

      else if ((iu.ge.7).and.(iu.le.19)) then

         open(unit=iu,file=fname,form='unformatted',&
                 access='direct',recl=4)

! ====================================================================
! If the unit number is 55 or greater then
! open that file for ASCII sequential input/output.
! ====================================================================

      else if ((iu.gt.55)) then

         open(unit=iu,file=fname)

! ====================================================================
! If the unit number is negative then done reading and return
! to calling routine.
! ====================================================================

      else

         return

      endif

! ====================================================================
! Set 'iprn' so that the output file will print.
! ====================================================================

      iprn(iu) = 1

! ====================================================================
! Read the next filename.
! ====================================================================
   
      goto 50

! ====================================================================
! Format statement
! ====================================================================

1000  format(i5,a100)        

      end subroutine filopn

! ====================================================================
!
!			subroutine rdebal
!
! ====================================================================
!
! Subroutine to read in energy balance calculation specification.
!
! ====================================================================

      subroutine rdebal(ioppet,iopwv,iopstab,iopgveg,iopthermc,&
                        iopthermc_v,maxnri,toleb)

      implicit none
      include "help/rdebal.h"

! ====================================================================
! Read in options for energy balance and calculation specifications. 
! ====================================================================

      ioppet = 0 !Always run in full water and energy balance
      iopwv = 1 !Always read in water vapor using relative humidity
      iopstab = 1 !Always perform stability correction on aero. resis.
      read(1000,*) iopgveg
      read(1000,*) iopthermc
      read(1000,*) iopthermc_v
      read(1000,*) maxnri
      read(1000,*) toleb

      return

      end subroutine rdebal

! ====================================================================
!
!			subroutine imgopn
!
! ====================================================================
!
! Subroutine to set the name of and open a file time varying
!   image input or output.
!
! ====================================================================

      subroutine imgopn(prefix,i,iu,img_opt,iyear,iday,ihour)

      implicit none
      include "help/imgopn.h"

! ====================================================================
! Initialize the file name to null and the file name suffix to 
! the time step, or if requested to the year_day_hour string.
! ====================================================================

      fname=''
      write (suffix,1001) iyear,iday,ihour

! ====================================================================
! Find the length of the prefix.
! ====================================================================

      lgth = index(prefix,' ')

! ====================================================================
! Build the file name one character at a time until a blank space
! in the prefix is encountered then add a period and the suffix.
! ====================================================================

      do 100 k=1,lgth-1

        fname(k:k)=prefix(k:k)

100   continue

      fname(lgth:lgth)='_'
      fname(lgth+1:lgth+13)=suffix

! ====================================================================
! Open the file for binary I/O.
! ====================================================================

         open(unit=iu,file=fname,form='unformatted',&
              access='direct',recl=4)

! ====================================================================
! Format statment.
! ====================================================================

1000  format (I5.5)
1001  format (I2.2,'_',I3.3,'_',I2.2,'.bin')

      return

      end subroutine imgopn

! ====================================================================
!
!  Parameter definitions:
!
!    icol:      loop index for image column
!    dummy:     dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    rmult:     all image output to be multiplied by this amount
! ====================================================================

      subroutine WRITE_BINARY(datain,rmult,nrow,ncol,ipixnum,i)

      implicit none
      real :: dummy 
      real :: rmult
      real*8 :: datain(nrow*ncol)
      real :: dataout(ncol,nrow)
      integer :: ipixnum(nrow,ncol)
      integer :: irow,icol,nrow,ncol,i,x,y

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      dummy = -999.0
       x = 1
       y = 0
       do irow = 1,nrow

         do icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

               if (y .eq. nrow)then
                        y = 0
                        x = x + 1
               endif

               y = y + 1

               if (ipixnum(irow,icol).gt.0) then

                  dataout(x,y) = rmult*datain(ipixnum(irow,icol))

               else
                  
                  dataout(x,y) = -999

               endif

        enddo

      enddo

      write (1005,rec=i) dataout

      return

      end subroutine WRITE_BINARY

! ====================================================================
!
!			subroutine wrimgi
!
! ====================================================================
!
! Subroutine to write an image of 4 byte integers from an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!  Parameter definitions:
!
!    ia:        values to write to the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    idummy:    dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    itmpvl:    temporary value used to hold the output value
!    iu:        unit number to read data
!    mult:      all image output to be multiplied by this amount
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
! ====================================================================

      subroutine wrimgi(ia,iu,mult,nrow,ncol,ipixnum)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      !include "sun_sgi.h"
      include "help/wrimgi.h"

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      do 200 irow = 1,nrow

         do 100 icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

            if (ipixnum(irow,icol).gt.0) then

                  itmpvl = mult*ia(ipixnum(irow,icol))
                  write(iu,rec=((irow-1)*ncol) + icol) itmpvl

            else
                  write(iu,rec=((irow-1)*ncol) + icol) idummy
         
            endif

100      continue

200   continue

      return

      end subroutine wrimgi

! ====================================================================
!
!			subroutine wrimgr
!
! ====================================================================
!
! Subroutine to write an image of 4 byte reals from an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!  Parameter definitions:
!
!    a:         values to write to the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    dummy:     dummy value to be writen to row-column location
!                 with no pixel number
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    rmult:     all image output to be multiplied by this amount
!    tmpval:    temporary value used to write the real value
!                 in four byte form when the input array
!                 is eight byte
! ====================================================================

      subroutine wrimgr(a,iu,rmult,nrow,ncol,ipixnum)

      implicit real*8 (a-h,o-z)
      !include "SNOW.h"
      !include "wgtpar.h"
      !include "sun_sgi.h"
      include "help/wrimgr.h"

! ====================================================================
! Loop through the image and write each value in proper location.
! ====================================================================

      dummy = 0.0

      do 200 irow = 1,nrow

         do 100 icol = 1,ncol

! --------------------------------------------------------------------&
! If the location is within the area of interest then
! write the correct value to the image, otherwise 
! write the dummy value.
! --------------------------------------------------------------------&

            if (ipixnum(irow,icol).gt.0) then

                  tmpval = rmult*a(ipixnum(irow,icol))
!	          call swap_r(tmpval,1)
                  write(iu,rec=((irow-1)*ncol) + icol) tmpval

            else
!	          call swap_r(dummy,1)	
                  write(iu,rec=((irow-1)*ncol) + icol) dummy
         
            endif

100      continue

200   continue

      return

      end subroutine wrimgr

! ====================================================================
!
!			subroutine rdsoil
!
! ====================================================================
!
! Subroutine to read and initialize time in-variant soil parameters
!
! ====================================================================

      subroutine rdsoil(nsoil,irestype,ikopt,zrzmax,iopsmini,&
       smpet0,isoil,nrow,ncol,ipixnum,bcbeta,psic,thetas,thetar,xk0,zdeep,&
       tdeep,zmid,tmid0,rocpsoil,quartz,ifcoarse,srespar1,srespar2,srespar3,&
       a_ice,b_ice,bulk_dens,amp,phase,shift,inc_frozen,bcgamm,par,corr,&
       idifind,ncatch,icatch,pixsiz,area,npix,psicav,iprn,tc,tw)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      include "help/rdsoil.h"
      type SoilDataTemplate
        real*8,dimension(:,:),allocatable :: bcbeta
        real*8,dimension(:,:),allocatable :: psic
        real*8,dimension(:,:),allocatable :: thetas
        real*8,dimension(:,:),allocatable :: thetar
        real*8,dimension(:,:),allocatable :: xk0
        real*8,dimension(:,:),allocatable :: zdeep
        real*8,dimension(:,:),allocatable :: tdeep
        real*8,dimension(:,:),allocatable :: zmid
        real*8,dimension(:,:),allocatable :: tmid0
        real*8,dimension(:,:),allocatable :: rocpsoil
        real*8,dimension(:,:),allocatable :: quartz
        integer,dimension(:,:),allocatable :: ifcoarse
        real*8,dimension(:,:),allocatable :: srespar1
        real*8,dimension(:,:),allocatable :: srespar2
        real*8,dimension(:,:),allocatable :: srespar3
        real*8,dimension(:,:),allocatable :: a_ice
        real*8,dimension(:,:),allocatable :: b_ice
        real*8,dimension(:,:),allocatable :: bulk_dens
        real*8,dimension(:,:),allocatable :: amp
        real*8,dimension(:,:),allocatable :: phase
        real*8,dimension(:,:),allocatable :: shift
        real*8,dimension(:,:),allocatable :: thetaw
        real*8,dimension(:,:),allocatable :: thetac
      end type SoilDataTemplate
      type (SoilDataTemplate) SoilData
      integer :: soilnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      soilnvars = 23
      allocate(TempArray(ncol,nrow,soilnvars))

! ====================================================================
! Read spatially constant bare soil parameters and options.
! Then read root and transmission zone data.
! ====================================================================

      !read(1000,*)nsoil
      nsoil = nrow*ncol

      read(1000,*)irestype
      read(1000,*)ikopt
      read(1000,*)zrzmax
      read(1000,*)iopsmini

      if (iopsmini.eq.0) read(1000,*)smpet0

      print*,"rdsoil:  Read spatially constant soil pars"

      if (iopsmini.eq.1)&
         print*,"rdsoil:  Will read initial soil moisture images"

! ====================================================================
! Read the binary soil file
! ====================================================================

      allocate(SoilData%bcbeta(ncol,nrow))
      allocate(SoilData%psic(ncol,nrow))
      allocate(SoilData%thetas(ncol,nrow))
      allocate(SoilData%thetar(ncol,nrow))
      allocate(SoilData%xk0(ncol,nrow))
      allocate(SoilData%zdeep(ncol,nrow))
      allocate(SoilData%tdeep(ncol,nrow))
      allocate(SoilData%zmid(ncol,nrow))
      allocate(SoilData%tmid0(ncol,nrow))
      allocate(SoilData%rocpsoil(ncol,nrow))
      allocate(SoilData%quartz(ncol,nrow))
      allocate(SoilData%ifcoarse(ncol,nrow))
      allocate(SoilData%srespar1(ncol,nrow))
      allocate(SoilData%srespar2(ncol,nrow))
      allocate(SoilData%srespar3(ncol,nrow))
      allocate(SoilData%a_ice(ncol,nrow))
      allocate(SoilData%b_ice(ncol,nrow))
      allocate(SoilData%bulk_dens(ncol,nrow))
      allocate(SoilData%amp(ncol,nrow))
      allocate(SoilData%phase(ncol,nrow))
      allocate(SoilData%shift(ncol,nrow))
      allocate(SoilData%thetaw(ncol,nrow))
      allocate(SoilData%thetac(ncol,nrow))

      print*,"rdsoil:  Reading in all soil properties at once"
      read(1001,rec=1)TempArray(:,:,:)
      SoilData%bcbeta(:,:) = dble(TempArray(:,:,1))
      SoilData%psic(:,:) = dble(TempArray(:,:,2))
      SoilData%thetas(:,:) = dble(TempArray(:,:,3))
      SoilData%thetar(:,:) = dble(TempArray(:,:,4))
      SoilData%xk0(:,:) = dble(TempArray(:,:,5))
      SoilData%zdeep(:,:) = dble(TempArray(:,:,6))
      SoilData%tdeep(:,:) = dble(TempArray(:,:,7))
      SoilData%zmid(:,:) = dble(TempArray(:,:,8))
      SoilData%tmid0(:,:) = dble(TempArray(:,:,9))
      SoilData%rocpsoil(:,:) = dble(TempArray(:,:,10))
      SoilData%quartz(:,:) = dble(TempArray(:,:,11))
      SoilData%ifcoarse(:,:) = TempArray(:,:,12)
      SoilData%srespar1(:,:) = dble(TempArray(:,:,13))
      SoilData%srespar2(:,:) = dble(TempArray(:,:,14))
      SoilData%srespar3(:,:) = dble(TempArray(:,:,15))
      SoilData%a_ice(:,:) = dble(TempArray(:,:,16))
      SoilData%b_ice(:,:) = dble(TempArray(:,:,17))
      SoilData%bulk_dens(:,:) = dble(TempArray(:,:,18))
      SoilData%amp(:,:) = dble(TempArray(:,:,19))
      SoilData%phase(:,:) = dble(TempArray(:,:,20))
      SoilData%shift(:,:) = dble(TempArray(:,:,21))
      SoilData%thetaw(:,:) = dble(TempArray(:,:,22))
      SoilData%thetac(:,:) = dble(TempArray(:,:,23))

! ====================================================================
! Read the soil classification image.
! ====================================================================

      call rdimgi(isoil,12,nrow,ncol,ipixnum)

      print*,"rdsoil:  Read soil texture image"

! ====================================================================
! Pass the soil properties from the original i,j pos. to the kk pos.
!  ====================================================================
      do kk=1,nsoil

         !Map the kk position to the i,j position
         if(mod(kk,nrow) .ne. 0)then
                ipos = kk/nrow+1
                jpos = mod(kk,nrow)
         else
                ipos = kk/nrow
                jpos = nrow
         endif

         bcbeta(kk) = SoilData%bcbeta(ipos,jpos)
         psic(kk) = SoilData%psic(ipos,jpos)
         thetas(kk) = SoilData%thetas(ipos,jpos)
         thetar(kk) = SoilData%thetar(ipos,jpos)
         xk0(kk) = SoilData%xk0(ipos,jpos)
         zdeep(kk) = SoilData%zdeep(ipos,jpos)
         tdeep(kk) = SoilData%tdeep(ipos,jpos)
         zmid(kk) = SoilData%zmid(ipos,jpos)
         tmid0(kk) = SoilData%tmid0(ipos,jpos)
         rocpsoil(kk) = SoilData%rocpsoil(ipos,jpos)
         quartz(kk) = SoilData%quartz(ipos,jpos)
         ifcoarse(kk) = SoilData%ifcoarse(ipos,jpos)
         srespar1(kk) = SoilData%srespar1(ipos,jpos)
         srespar2(kk) = SoilData%srespar2(ipos,jpos)
         srespar3(kk) = SoilData%srespar3(ipos,jpos)
         a_ice(kk) = SoilData%a_ice(ipos,jpos)
         b_ice(kk) = SoilData%b_ice(ipos,jpos)
         bulk_dens(kk) = SoilData%bulk_dens(ipos,jpos)
         amp(kk) = SoilData%amp(ipos,jpos)
         phase(kk) = SoilData%phase(ipos,jpos)
         shift(kk) = SoilData%shift(ipos,jpos)
         tc(kk) = SoilData%thetac(ipos,jpos)
         tw(kk) = SoilData%thetaw(ipos,jpos)

      enddo
 
      inc_frozen = 1 !THIS MEANS THAT THE FROZEN ALGORITHM IS ALWAYS RUN

! ====================================================================
! Calculate time in-variant soil parameters for each soil class.
! ====================================================================

      do 400 kk=1,nsoil

! --------------------------------------------------------------------&
! Calculate soil parameters based on Brooks-Corey parameters.
! --------------------------------------------------------------------&

         bcgamm(kk) = two + three * bcbeta(kk)

! --------------------------------------------------------------------&
! Calculate constants for bare soil evaporation desorptivity 
! equation used in Famiglietti PhD Thesis, Princetion Univ, 1992.
! --------------------------------------------------------------------&

         par(kk) = one + ((bcgamm(kk)-one)/bcbeta(kk))
         corr(kk)=((one/(one+three*bcbeta(kk)))-&
                   (0.85d0/(one+four*bcbeta(kk)))-&
                   (0.85d0*0.15d0*0.5d0/(one+five*bcbeta(kk)))+&
                   (0.85d0*0.15d0*1.15d0/&
                   (six*(one+six*bcbeta(kk)))))

! --------------------------------------------------------------------&
! Calculate diffusivity index and dimensionless exfiltration
! diffusivity from Eagleson, WRR, 1972.
! --------------------------------------------------------------------&

         idifind(kk) = ((1.0+2.0*bcbeta(kk))/bcbeta(kk))+0.5
         tempsum=0 

         do 300 nn=1,idifind(kk)

            dtaken = exp(factln(idifind(kk))-factln(nn)-&
                       factln(idifind(kk)-nn))
            tempsum = tempsum+(((-1)**nn)*dtaken/(1.85+nn))

300      continue

400   continue

      print*,"rdsoil:  Calculated time-invariant soil pars"

! ====================================================================
! Calculate the fraction of different soil types in each catchment
! and in total area.  First array index for icount and frcov is
! the soil type, second array index is the catchment number.
! Catchment ncatch+1 is total area.
! ====================================================================

      do 450 kk=1,nsoil

         do 440 jj=1,ncatch+1

            icount(kk,jj) = 0

440      continue

450   continue

      do 500 kk=1,npix

        icount(isoil(kk),icatch(kk))=icount(isoil(kk),icatch(kk))+1
        icount(isoil(kk),ncatch+1) = icount(isoil(kk),ncatch+1) + 1

500   continue

      do 550 kk=1,nsoil

         do 540 jj=1,ncatch

            frsoil(kk,jj) = icount(kk,jj)*pixsiz*pixsiz/area(jj)

540      continue

         frsoil(kk,ncatch+1) = icount(kk,ncatch+1)/real(npix)

550   continue

      print*,"rdsoil:  Calculated fractional coverage for soil types"

! ====================================================================
! Calculate average bubbling pressure in each catchment (used
! in updating average water table depths.
! ====================================================================

      do 570 jj=1,ncatch

         psicav(jj) = zero

         do 560 kk=1,nsoil

            psicav(jj) = psicav(jj) + frsoil(kk,jj)*psic(kk)

560      continue

570   continue

      print*,"rdsoil:  Calculated average psi! for each catchment"

! ====================================================================
! Print summary table.
! ====================================================================

      if (iprn(78).eq.1) then

         do 580 kk=1,nsoil

            write(78,1000) kk,bcbeta(kk),psic(kk)*100,thetas(kk),&
                           thetar(kk),xk0(kk)*3600000.

580      continue

      endif

! ====================================================================
! Print fractional coverage of soils types in each catchment.
! ====================================================================

      if (iprn(79).eq.1) then

         do 700 jj=1,ncatch

            write (79,1100) jj,area(jj)/1000000.

            do 650 kk=1,nsoil

               write(79,1110) kk,frsoil(kk,jj)*100

650         continue

            write(79,*)

700      continue

         write(79,1120)

         do 800 kk=1,nsoil

            write(79,1110) kk,frsoil(kk,ncatch+1)*100

800      continue

      endif        

      print*,"rdsoil:  Printed summary tables"

! ====================================================================
! Format statement.
! ====================================================================

1000  format (i5,5f10.3)
1100  format ('Catchment Number',i4,' (Area = ',f8.2,' km^2):')
1110  format ('   Soil Type',i4,':',f7.2,'%')
1120  format ('Total Area:')

      return

      end subroutine rdsoil

! ====================================================================
!
!			subroutine rdpet
!
! ====================================================================
!
! Subroutine to read in potential evapotranspiration (m/s) from a file.
! It also sets evaporation from bare soil equal to evaporation
! from wet canopy equal to unstressed transpiration .
!
! ====================================================================

      subroutine rdpet(epetd,epetw,ebspot,xled,xlew,xlhv,row)

      implicit none
      include "help/rdpet.h"

! ====================================================================
! Set potential evapotranspiration for wet and dry canopy.and.&
! bare soil to the input value.
! ====================================================================
      
      epetd = ebspot/(3600*1000)
      epetw = ebspot/(3600*1000)

! ====================================================================
! Calculate mass flux.
! ====================================================================

      epetdmf = epetd*row
      epetwmf = epetw*row

! ====================================================================
! Calculate latent heat flux.
! ====================================================================

      xled = epetdmf*xlhv
      xlew = epetwmf*xlhv

      return

      end subroutine rdpet

! ====================================================================
!
!			subroutine rdimgi
!
! ====================================================================
!
! Subroutine to read in an image of integers and return an array
!   of these values indexed by the soils-topographi! index
!   pixel numbers.
!
! ====================================================================
!
!
!  Parameter definitions:
!
!    ia:        values read from the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    itmpval:   temporary 4 byte integer value read from the image
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
! ====================================================================

      subroutine rdimgi(ia,iu,nrow,ncol,ipixnum)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      !include "sun_sgi.h"
      include "help/rdimgi.h"

! ====================================================================
! Loop through the image and read each value.
! ====================================================================

      do irow = 1,nrow

         do icol = 1,ncol
                
             if (iu .eq. 12 .or. iu .eq. 11)then
                itmpval = (irow-1)*ncol + icol
             else
                read(iu,rec=((irow-1)*ncol) + icol) itmpval
             endif
! --------------------------------------------------------------------
! If the location is within the area of interest then
! assign to array 'ia', otherwise read next value.
! --------------------------------------------------------------------

            if (ipixnum(irow,icol).gt.0) then

               ia(ipixnum(irow,icol)) = itmpval

            endif

          enddo
        
      enddo

      return

      end subroutine rdimgi

! ====================================================================
!
!			subroutine rdimgr
!
! ====================================================================
!
! Subroutine to read in an image of 4 byte reals and return an array
!   of these values indexed by the soils-topographic index
!   pixel numbers.
!
! ====================================================================
!
!
!  Parameter definitions:
!
!    a:         values read from the image in an array indexed by 
!                 pixel number
!    icol:      loop index for image column
!    ipixnum:   pixel number of image row/column location
!    irow:      loop index for image row
!    iu:        unit number to read data
!    ncol:      number of columns in the image
!    nrow:      number of rows in the image
!    tmpval:    temporary 4 byte real value read from the image
! ====================================================================

      subroutine rdimgr(a,iu,nrow,ncol,ipixnum)

      implicit real*8 (a-h,o-z)
      !include "SNOW.h"
      !include "wgtpar.h"
      include "help/rdimgr.h"

! ====================================================================
! Loop through the image and read each value.
! ====================================================================
      do 200 irow = 1,nrow
         do 100 icol = 1,ncol
               read(iu,rec=((irow-1)*ncol) + icol) tmpval

! --------------------------------------------------------------------
! If the location is within the area of interest then
! assign to array 'a', otherwise read next value.
! --------------------------------------------------------------------

            if (ipixnum(irow,icol).gt.0) then

               a(ipixnum(irow,icol)) = tmpval
            
            endif

100      continue

200   continue

      return
      end subroutine rdimgr

! ====================================================================
!
!			subroutine rdatb
!
! ====================================================================
!
! Subroutine to read in the topographi! index image and set up
! at translation between the pixel number and the image row.and.&
! column.
!
! ====================================================================

      subroutine rdatb(atb,nrow,ncol,ipixnum,ixpix,iypix,npix)

      implicit none
      !include "SNOW.h"
      !include "wgtpar.h"
      !include "sun_sgi.h"
      include "help/rdatb.h"

! ====================================================================
! Initialize pixel number counter.
! ====================================================================

      ip = 1

! ====================================================================
! Loop through image row by row.
! ====================================================================

      do 500 irow=1,nrow

         do 400 icol=1,ncol

! --------------------------------------------------------------------
! Read image value at irow,icol.
! --------------------------------------------------------------------


               read(9,rec=((irow-1)*ncol) + icol) tempatb

! --------------------------------------------------------------------
! Check for atb value at current location.  If there is
! a value then record the value and location.  If not
! then mark location as outside the area of interest.
! --------------------------------------------------------------------

            if (tempatb.gt.0.0) then

               atb(ip) = tempatb
               ipixnum(irow,icol) = ip
               ixpix(ip) = icol
               iypix(ip) = irow
               ip = ip + 1

            else

               ipixnum(irow,icol) = 0

            endif

400      continue

500   continue

! ====================================================================
! Set the total number of pixels.
! ====================================================================

      npix = ip - 1

      if (npix.gt.MAX_PIX) then

         write (*,*) 'rdatb : npix greater then MAX_PIX ',npix,MAX_PIX
         stop

      endif

      return

      end subroutine rdatb

! ====================================================================
!
!             subroutine clc_ind
!
! ====================================================================
!
! Calculate the matrix indexes.
!
! ====================================================================

      subroutine clc_ind (ilandc,ipix,SNOW_RUN,sw_lc,sw_px,SNW_FLG,s_lc,s_px)

      implicit none

      include 'help/clc_ind.h'

! --------------------------------------------------------------------
! Snow indexes, overstory.
! --------------------------------------------------------------------

      if (SNOW_RUN.eq.1) then

         sw_lc=ilandc
         sw_px=ipix

      endif

      if (SNOW_RUN.eq.0) then

         sw_lc=1
         sw_px=1

      endif

! --------------------------------------------------------------------
! Snow indexes, understory.
! --------------------------------------------------------------------

      if (SNW_FLG.eq.1) then

         s_lc=ilandc
         s_px=ipix

      endif

      if (SNW_FLG.eq.0) then

         s_lc=1
         s_px=1

      endif

      return

      end subroutine clc_ind

END MODULE MODULE_IO
