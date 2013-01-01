MODULE MODULE_IO

USE MODULE_VARIABLES

!Add the variables that are reused throughout the subroutines

implicit none

contains

!####################################################################
! Subroutine to read in and pass meteorological data (e.g. rainfall).
!####################################################################

      subroutine rdatmo(nrow,ncol,ipixnum,iyear,iday,ihour,i,ATMOS)

      implicit none
      type (GRID_MET_template) :: ATMOS(nrow*ncol)
      type (GLOBAL_template) :: GLOBAL
      integer :: nrow,ncol,ipixnum(nrow,ncol)
      integer :: iyear,iday,ihour
      integer :: forcingnvars,i
      real,dimension(:,:,:),allocatable :: TempArray
      forcingnvars = 7
      GLOBAL%ncol = ncol
      GLOBAL%nrow = nrow
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,forcingnvars))

! ====================================================================
! Read year, day and hour
! ====================================================================

      read(61,*) iyear,iday,ihour

! ####################################################################
! Read all variables in at once for each time step
! ####################################################################

      read(1004,rec=i) TempArray(:,:,:)

! Longwave Radiation

      call rdforc(ATMOS%rld,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,1))

! Air Pressure

      call rdforc(ATMOS%press,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,2))

! Relative Humidity

      call rdforc(ATMOS%rh,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,3))

! Shortwave Radiation

      call rdforc(ATMOS%rsd,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,4))

! Air Temperature

      call rdforc(ATMOS%tdry,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,5))

! Wind Speed

      call rdforc(ATMOS%uzw,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,6))

! Precipitation

      call rdforc(ATMOS%pptms,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,7))
       
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

      subroutine rddata(GLOBAL,GRID,REG,CAT)

      implicit none
      integer iophd
      type (GLOBAL_template) :: GLOBAL
      !type (GRID%SOIL_template) :: GRID%SOIL
      type (GRID_template),dimension(:),allocatable :: GRID
      type (REGIONAL_template) :: REG
      type (CATCHMENT_template),dimension(:),allocatable :: CAT
      character(len=200) :: filename

! ====================================================================
! Open input/output files and set variable to control output file
! printing.
! ====================================================================

      filename = '/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/TI'
      open(unit=9,file=filename,form='unformatted',access='direct',recl=4)
      filename = '/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/subbasin'
      open(unit=10,file=filename,form='unformatted',access='direct',recl=4)
      filename = '/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/K_0.img'
      open(unit=8,file=filename,form='unformatted',access='direct',recl=4)
      filename = '/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/dat.61.input'
      open(unit=61,file=filename)
      filename = '/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/CL_Table'
      open(unit=71,file=filename)

! ====================================================================
! Read in simulation time constants and control variables.
! ====================================================================

      read(1000,*) GLOBAL%ndata
      read(1000,*) GLOBAL%dt
      read(1000,*) GLOBAL%endstm
      read(1000,*) iophd

      print*, 'rddata:  Done reading time parameters'
      print*, 'rddata:  Total time steps = ',GLOBAL%ndata

! ====================================================================
! Read and initialize topmodel parameters, atb distribution and
! initial water table depth.
! ====================================================================

      call rdtpmd(GLOBAL%iopbf,GLOBAL%iopwt0,GLOBAL%ncatch,GLOBAL%nrow,&
       GLOBAL%ncol,GLOBAL%pixsiz,GLOBAL%ipixnum,GLOBAL%iprn,GLOBAL%ixpix,GLOBAL%iypix,&
       GLOBAL%npix,GLOBAL%q0,GLOBAL%ff,GLOBAL%qb0,&
       GLOBAL%dd,GLOBAL%xlength,GLOBAL%basink,GLOBAL%xlamda,GLOBAL%icatch,&
       GLOBAL%area,GLOBAL%atanb,GLOBAL%dtil,GLOBAL%zbar1,&
       GLOBAL%iwel,GLOBAL%wslp,&
       GLOBAL%lat_deg,GLOBAL%lat_min,GLOBAL%lng_deg,&
       GLOBAL%lng_min,GLOBAL%lng_mer,GLOBAL%rlatitude,&
       GLOBAL%rlongitude,GLOBAL%rlng_merid,GRID,CAT)

      print*,'rddata:  Done reading TOPMODEL parameters'

! ====================================================================
! Read in and initialize vegatation parameters.
! ====================================================================

      call rdveg(GLOBAL%npix,GLOBAL%nrow,GLOBAL%ncol,GLOBAL%ilandc,&
       GLOBAL%ipixnum,GLOBAL%nlandc,GLOBAL%iopveg,GRID%VEG%ivgtyp,GLOBAL%iprn,&
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
       GLOBAL%pixsiz,GLOBAL%area,CAT%fbs,&
       REG%fbsrg,GLOBAL%ncatch)

      print*,'rddata:  Done reading vegetation parameters'

! ====================================================================
! Read in soil parameters and root and transmission zone information.
! ====================================================================

      call rdsoil(GLOBAL%nsoil,GLOBAL%irestype,GLOBAL%ikopt,GLOBAL%zrzmax,GLOBAL%iopsmini,GLOBAL%smpet0,&
       GLOBAL%isoil,GLOBAL%nrow,GLOBAL%ncol,GLOBAL%ipixnum,GRID%SOIL%bcbeta,&
       GRID%SOIL%psic,GRID%SOIL%thetas,GRID%SOIL%thetar,GRID%SOIL%xk0,GRID%SOIL%zdeep,GRID%SOIL%tdeep,GRID%SOIL%zmid,&
       GRID%SOIL%tmid0,GRID%SOIL%rocpsoil,GRID%SOIL%quartz,GLOBAL%ifcoarse,&
       GRID%SOIL%srespar1,GRID%SOIL%srespar2,GRID%SOIL%srespar3,GRID%SOIL%a_ice,GRID%SOIL%b_ice,&
       GRID%SOIL%bulk_dens,GRID%SOIL%amp,GRID%SOIL%phase,GRID%SOIL%shift,&
       GLOBAL%inc_frozen,GRID%SOIL%bcgamm,GRID%SOIL%par,GRID%SOIL%corr,GLOBAL%idifind,&
       GLOBAL%ncatch,GLOBAL%icatch,GLOBAL%pixsiz,GLOBAL%area,&
       GLOBAL%npix,CAT%psicav,GLOBAL%iprn,GRID%VEG%tc,GRID%VEG%tw)

      print*,'rddata:  Done reading soil parameters'

! ====================================================================
! Read in GLOBAL for energy balance and calculate soil thermal
! conductivity. 
! ====================================================================

      GLOBAL%ioppet = 0 !Always run in full water and energy balance
      GLOBAL%iopwv = 1 !Always read in water vapor using relative humidity
      GLOBAL%iopstab = 1 !Always perform stability correction on aero. resis.
      read(1000,*) GLOBAL%iopgveg
      read(1000,*) GLOBAL%iopthermc
      read(1000,*) GLOBAL%iopthermc_v
      read(1000,*) GLOBAL%maxnri
      read(1000,*) GLOBAL%toleb

      print*,'rddata:  Done reading energy balance parameters'

! ====================================================================
! Read in the mode in which to run the program.
! ====================================================================

      GLOBAL%MODE = 1
      GLOBAL%FRCOV = 0
      GLOBAL%frcbeta = 999

      if (GLOBAL%MODE.eq.1) print*,'rddata:  Running the model in dist. mode'
      if (GLOBAL%FRCOV.eq.0) print*,'rddata:  Frac. rain.l cover not included'

! ====================================================================
! Initialize the simulation sum variables and storm.and.&
! interstorm flags and times.
! ====================================================================

      call inisim(GLOBAL%iopsmini,GLOBAL%nrow,GLOBAL%ncol,GLOBAL%ipixnum,GLOBAL%ilandc,&
       GLOBAL%npix,GLOBAL%inc_frozen,GLOBAL%istorm,&
       GLOBAL%intstm,GLOBAL%istmst,intstp,GLOBAL%istorm_moss,&
       GLOBAL%intstm_moss,GLOBAL%istmst_moss,GLOBAL%intstp_moss,&
       GLOBAL%isoil,GLOBAL%idifind,GLOBAL%smpet0,r_mossmpet0,GLOBAL%endstm,&
       GRID%VARS%rzsm1,GRID%VARS%tzsm1,GRID%VARS%r_mossm1,&
       GRID%VARS%r_mossm,GRID%VARS%rzsm1_u,GRID%VARS%tzsm1_u,&
       GRID%VARS%rzsm1_f,GRID%VARS%tzsm1_f,GRID%VARS%r_mossm1_u,&
       GRID%VARS%r_mossm_u,&
       GRID%VARS%r_mossm1_f,GRID%VARS%r_mossm_f,GRID%VARS%rzdthetaidt,&
       GRID%VARS%tzdthetaidt,GRID%VARS%zmoss,r_moss_depth,&
       thetas_moss,GLOBAL%xintst,GLOBAL%xintst_moss,GLOBAL%cuminf,GRID%SOIL%xk0,GRID%SOIL%psic,&
       GRID%SOIL%thetas,GRID%SOIL%thetar,GRID%SOIL%bcgamm,&
       bcbeta,GLOBAL%sorp,GLOBAL%cc,GLOBAL%dt,GLOBAL%sesq,GRID%SOIL%corr,GRID%SOIL%par,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,r_MeltEnergy_us,Outflow_us,&
       PackWater,SurfWater,Swq,VaporMassFlux,r_MeltEnergy,Outflow)

      read (1000,*) GLOBAL%dtveg

      print*,'rddata:  Done initializing simulation'

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
       GLOBAL,ntdveg,GRID)

      implicit none
      !include "wgtpar.h"
      include "help/rdveg_update.h"
      type VegDataTemplate
        real*8,dimension(:,:),allocatable :: xlai
        real*8,dimension(:,:),allocatable :: albd
      end type VegDataTemplate
      type (VegDataTemplate) VegData
      integer :: dvegnvars,ipos,jpos
      integer :: ntdveg
      real,dimension(:,:,:),allocatable :: TempArray
      type (GLOBAL_template) :: GLOBAL
      type (GRID_template),dimension(:),allocatable :: GRID
      dvegnvars = 2

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      ntdveg = ntdveg + 1
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,dvegnvars))
      allocate(VegData%xlai(GLOBAL%ncol,GLOBAL%nrow))
      allocate(VegData%albd(GLOBAL%ncol,GLOBAL%nrow))
      read(1003,rec=ntdveg)TempArray(:,:,:)
      VegData%xlai(:,:) = dble(TempArray(:,:,1))
      VegData%albd(:,:) = dble(TempArray(:,:,2))

! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

        do kk=1,GLOBAL%nlandc

                !Map the kk position to the i,j position
                if(mod(kk,GLOBAL%nrow) .ne. 0)then
                        ipos = kk/GLOBAL%nrow+1
                        jpos = mod(kk,GLOBAL%nrow)
                else
                        ipos = kk/GLOBAL%nrow
                        jpos = GLOBAL%nrow
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

      do 300 kk=1,GLOBAL%nlandc

! --------------------------------------------------------------------&
! If not bare soil then calculate the canopy resistance.
! --------------------------------------------------------------------&

         if (GRID(kk)%VEG%ivgtyp.ne.0) then

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
      !include "SNOW.h"
      !include "wgtpar.h"
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

!Regional Precipitation/Infiltration/Runoff (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/PRI_fruit.txt"
open(2093,file=trim(filename))

!Regional Evapotranspiration Rates (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/ET_fruit.txt"
open(2094,file=trim(filename))

!Regional Root and Transmission Zone Water Balance (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/STZB_fruit.txt"
open(2095,file=trim(filename))

!Regional Water Table Balance (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/WTB_fruit.txt"
open(2096,file=trim(filename))

!Regional Fractional Saturation States (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/SSR_fruit.txt"
open(2097,file=trim(filename))

!Regional Evapotranspiration Controls and Infiltration Mechanisms (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/RMEC_fruit.txt"
open(2098,file=trim(filename))

!Regional Snow Cover (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/SWE_fruit.txt"
open(2099,file=trim(filename))

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

!Regional Precipitation/Infiltration/Runoff (VALIDATION FILE)
close(2093)

!Regional Evapotranspiration Rates (VALIDATION FILE)
close(2094)

!Regional Root and Transmission Zone Water Balance (VALIDATION FILE)
close(2095)

!Regional Water Table Balance (VALIDATION FILE)
close(2096)

!Regional Fractional Saturation States (VALIDATION FILE)
close(2097)

!Regional Evapotranspiration Controls and Infiltration Mechanisms (VALIDATION FILE)
close(2098)

!Regional Snow Cover (VALIDATION FILE)
close(2099)

end subroutine FILE_CLOSE

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

  read(1000,*) wc0

  do kk=1,npix

    wcip1(kk) = wc0

  enddo  

  print*,"rdveg:  Read initial wet canopy storages"


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

      return

      end subroutine rdveg

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
! Read spatially constant bare soil parameters and GLOBAL.
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

! ====================================================================
!
!                       function factln
!
! ====================================================================
!
! Returns the natural log of the i factorial.
!
! ====================================================================
!
!
!  Variable definitions:
!  
!    factln:     The return value (ln (i factorial))
!    i:          Argument for which factorial is calculated
!    jj:         Loop index for series of adding logs
!    x:          Real*8 value used to send to dlog function
! ====================================================================

      function factln(i)
      implicit none
      include "help/factln.h"

! ====================================================================
! Calculate the log of factorial by adding all logs up
! to the argument value.
! ====================================================================

      factln = 0.0

      do 100 jj=1,i

        x=jj
        factln = factln + dlog(x)

100   continue

      return

      end function factln


! ====================================================================
!
!                   subroutine bilinear_interpolation
!
! ====================================================================
!
! Subroutine to downscale data using bilinear interpolation
!
! ====================================================================

!subroutine bilinear_interpolation()

!end subroutine biliner_interpolation

END MODULE MODULE_IO
