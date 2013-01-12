MODULE MODULE_IO

USE MODULE_VARIABLES

!Add the variables that are reused throughout the subroutines

implicit none

type IO_template
  integer,allocatable,dimension(:,:) :: ipixnum
  integer,allocatable,dimension(:) :: ixpix,iypix
end type IO_template

contains

!####################################################################
! Subroutine to read in and pass meteorological data (e.g. rainfall).
!####################################################################

      subroutine rdatmo(i,MET,GLOBAL,IO)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_MET_template),intent(inout) :: MET(GLOBAL%nrow*GLOBAL%ncol)
      type (IO_template),intent(in) :: IO
      integer :: ipixnum(GLOBAL%nrow,GLOBAL%ncol)
      integer :: forcingnvars,i
      real,dimension(:,:,:),allocatable :: TempArray
      ipixnum = IO%ipixnum
      forcingnvars = 7
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,forcingnvars))

! ====================================================================
! Read year, day and hour
! ====================================================================

      read(61,*) GLOBAL%iyear,GLOBAL%iday,GLOBAL%ihour

! ####################################################################
! Read all variables in at once for each time step
! ####################################################################

      read(1004,rec=i) TempArray(:,:,:)

! Longwave Radiation

      call rdforc(MET%rld,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,1))

! Air Pressure

      call rdforc(MET%press,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,2))

! Relative Humidity

      call rdforc(MET%rh,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,3))

! Shortwave Radiation

      call rdforc(MET%rsd,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,4))

! Air Temperature

      call rdforc(MET%tdry,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,5))

! Wind Speed

      call rdforc(MET%uzw,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,6))

! Precipitation

      call rdforc(MET%pptms,GLOBAL%nrow,GLOBAL%ncol,ipixnum,TempArray(:,:,7))
       
      end subroutine

! ####################################################################
! Subroutine to convert data to original setup
! ####################################################################

      subroutine rdforc(pixdat,nrow,ncol,ipixnum,data_in)

      implicit none
      integer,intent(in) :: nrow,ncol
      integer :: l,m
      real,dimension(:,:),intent(in) :: data_in(ncol,nrow)
      integer,dimension(:,:),intent(in) :: ipixnum
      real*8,dimension(:),intent(inout) :: pixdat

      do l=1,ncol
        do m=1,nrow
          if (ipixnum(m,l) .gt. 0)pixdat(ipixnum(m,l)) = data_in(l,m)
        enddo
      enddo

      end subroutine rdforc

! ####################################################################
! Subroutine to open input/output files, read and initialize time 
! in-variant data.
! ####################################################################

      subroutine rddata(GLOBAL,GRID,REG,CAT,IO)

      implicit none
      integer iophd
      type (GLOBAL_template) :: GLOBAL
      type (GRID_template),dimension(:),allocatable :: GRID
      type (REGIONAL_template) :: REG
      type (CATCHMENT_template),dimension(:),allocatable :: CAT
      type (IO_template),intent(inout) :: IO
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
      read(1000,*) GLOBAL%iophd

      print*, 'rddata:  Done reading time parameters'
      print*, 'rddata:  Total time steps = ',GLOBAL%ndata

! ====================================================================
! Read and initialize topmodel parameters, atb distribution and
! initial water table depth.
! ====================================================================

      call rdtpmd(GRID,CAT,IO,GLOBAL)

      print*,'rddata:  Done reading TOPMODEL parameters'

! ====================================================================
! Read in and initialize vegatation parameters.
! ====================================================================

       call rdveg(GLOBAL,CAT,GRID,REG,IO)

      print*,'rddata:  Done reading vegetation parameters'

! ====================================================================
! Read in soil parameters and root and transmission zone information.
! ====================================================================

      call rdsoil(GLOBAL,CAT,GRID,IO)

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

      GLOBAL%frcbeta = 999

! ====================================================================
! Initialize the simulation sum variables and storm.and.&
! interstorm flags and times.
! ====================================================================

      call inisim(GLOBAL,IO,GRID)

      read (1000,*) GLOBAL%dtveg

      print*,'rddata:  Done initializing simulation'

! Initialize the vegetation time step

      GLOBAL%ntdveg = 1

! Set constants
  GRID%VARS%row = 997.d0
  GRID%VARS%cph2o = 4186.d0
  GRID%VARS%cp = 1005.d0
  GRID%VARS%roi = 850.d0
  GRID%VARS%rzsmold = 0.d0

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

      subroutine rdveg_update (GLOBAL,GRID)

      implicit none
      integer :: kk,jj
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),intent(inout) :: GRID
      integer :: dvegnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%ncol,GLOBAL%nrow)
      dvegnvars = 2

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      GLOBAL%ntdveg = GLOBAL%ntdveg + 1
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,dvegnvars))
      read(1003,rec=GLOBAL%ntdveg)TempArray
      GRID_VEG_2D%xlai = dble(TempArray(:,:,1))
      GRID_VEG_2D%albd = dble(TempArray(:,:,2))

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

                GRID(kk)%VEG%xlai = GRID_VEG_2D(ipos,jpos)%xlai !dveg
                GRID(kk)%VEG%albd = GRID_VEG_2D(ipos,jpos)%albd !dveg
                GRID(kk)%VEG%tcbeta = exp(-0.5*GRID(kk)%VEG%xlai)
                GRID(kk)%VEG%xlai_wsc = GRID_VEG_2D(ipos,jpos)%xlai
                GRID(kk)%VEG%albw = GRID_VEG_2D(ipos,jpos)%albd

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

!#####################################################################
!			subroutine rdtpmd
!#####################################################################
! Subroutine to read in and initialize topmodel parameters and the
! soils-topographic index map.
!#####################################################################

subroutine rdtpmd(GRID,CAT,IO,GLOBAL)

  implicit none
  type (GLOBAL_template) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),allocatable :: CAT
  type (IO_template),intent(inout) :: IO
  integer kk,jj
  real*8 hbar0
  integer,dimension(:),allocatable :: icount
  real*8,dimension(:),allocatable :: atb,ti,zbar0,sumatb,sumlti,qb0,lte
  real*8,dimension(:),allocatable :: ki

  !###################################################################
  ! Read the option for baseflow calculation, the option .or.&
  ! initial water table entry, the number of basin in the area
  ! of interest, and the dimensions of the soils-topographi!
  ! index map.
  !###################################################################

  read(1000,*) GLOBAL%iopbf
  read(1000,*) GLOBAL%iopwt0
  read(1000,*) GLOBAL%ncatch
  read(1000,*) GLOBAL%nrow
  read(1000,*) GLOBAL%ncol
  read(1000,*) GLOBAL%pixsiz

  !#####################################################################
  ! Allocate memory
  !#####################################################################

  allocate(CAT(GLOBAL%ncatch))
  allocate(GRID(GLOBAL%nrow*GLOBAL%ncol))
  allocate(IO%ipixnum(GLOBAL%nrow,GLOBAL%ncol))
  allocate(IO%ixpix(GLOBAL%nrow*GLOBAL%ncol))
  allocate(IO%iypix(GLOBAL%nrow*GLOBAL%ncol))
  allocate(icount(GLOBAL%ncatch))
  allocate(atb(GLOBAL%nrow*GLOBAL%ncol))
  allocate(ti(GLOBAL%nrow*GLOBAL%ncol))
  allocate(zbar0(GLOBAL%ncatch))
  allocate(sumatb(GLOBAL%ncatch))
  allocate(sumlti(GLOBAL%ncatch))
  allocate(qb0(GLOBAL%ncatch))
  allocate(lte(GLOBAL%ncatch))
  allocate(ki(GLOBAL%nrow*GLOBAL%ncol))

  !####################################################################
  ! Read in the soils-topographic index map and set up the 
  ! transformation from pixel number to row/column.
  !####################################################################

  call rdatb(atb,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum,IO%ixpix,IO%iypix,GLOBAL%npix)

! ====================================================================
! Read in the catchment look-up table - read different values based
! on what is necessary for baseflow and initial condition calculations.
! ====================================================================

      if (GLOBAL%iopwt0.eq.1) then

         do 100 kk=1,GLOBAL%ncatch

            read(71,*) jj,CAT(kk)%q0,CAT(kk)%ff,CAT(kk)%qb0,CAT(kk)%dd,&
                       CAT(kk)%xlength,CAT(kk)%basink

100      continue

      else if (GLOBAL%iopbf.eq.1) then

         do 200 kk=1,GLOBAL%ncatch

            read(71,*) jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk),CAT(kk)%dd,&
                       CAT(kk)%xlength,CAT(kk)%basink 

200      continue

      else 

         do 300 kk=1,GLOBAL%ncatch

            read(71,*) jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk)

300      continue

      endif

! ====================================================================
! Read the catchment image.
! ====================================================================

      call rdimgi(GRID%VARS%icatch,10,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)

! ====================================================================
! Read image of transmissivities for use in calculating the 
! soils-topographi! index.
! ====================================================================

      call rdimgr(ki,8,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)
      do kk=1,GLOBAL%nrow*GLOBAL%ncol
        if (GRID(kk)%VARS%icatch .gt. 0)then
          ti(kk) = ki(kk)/CAT(GRID(kk)%VARS%icatch)%ff
        else
          ti(kk) = 0.0
        endif
      enddo

! ====================================================================
! Calculate the average topographi! index value, the average of 
! the natural log of transmissivity and area for each basin.
! ====================================================================

      do 400 kk=1,GLOBAL%ncatch

         sumatb(kk) = 0.0
         sumlti(kk) = 0.0
         icount(kk) = 0

400   continue

      do 500 kk=1,GLOBAL%npix

         sumatb(GRID(kk)%VARS%icatch) = sumatb(GRID(kk)%VARS%icatch) + atb(kk)
         sumlti(GRID(kk)%VARS%icatch) = sumlti(GRID(kk)%VARS%icatch) + dlog(ti(kk))
         icount(GRID(kk)%VARS%icatch) = icount(GRID(kk)%VARS%icatch) + 1

500   continue

      do 600 kk=1,GLOBAL%ncatch

         if (icount(kk).eq.0) then

            CAT(kk)%xlamda = 0
            lte(kk) = 0

         else

            CAT(kk)%xlamda = sumatb(kk)/icount(kk)
            lte(kk) = sumlti(kk)/icount(kk)

         endif

         CAT(kk)%area = icount(kk)*GLOBAL%pixsiz*GLOBAL%pixsiz

600   continue

      print*,'Area',CAT(1)%area
      print*,'ln Te',lte(1)
      print*,'Lambda',CAT(1)%xlamda

! ====================================================================
! Calculate soils-topographi! index for each pixel.
! ====================================================================

      do 50 kk=1,GLOBAL%npix

         GRID(kk)%VARS%atanb = atb(kk) + lte(GRID(kk)%VARS%icatch) - dlog(ti(kk))

50    continue

! ====================================================================
! Calculate the initial water table depth for each catchment.
! ====================================================================

      if ((GLOBAL%iopwt0.eq.1).or.(GLOBAL%iopbf.eq.1)) then

         do 700 kk=1,GLOBAL%ncatch

            CAT(kk)%dtil = sqrt(CAT(kk)%q0/(3.45*CAT(kk)%basink*CAT(kk)%dd*CAT(kk)%xlength))

            if (GLOBAL%iopwt0.eq.1) then

               hbar0 = sqrt(CAT(kk)%qb0/(5.772*CAT(kk)%basink*CAT(kk)%dd*CAT(kk)%xlength))
               zbar0(kk) = CAT(kk)%dtil - hbar0

            endif

700      continue

      endif

! ====================================================================
! Set initial average water table depth to the water table depth 
! at the end of the previous time step since program updates 
! this depth.
! ====================================================================

      do 800 kk=1,GLOBAL%ncatch

         CAT(kk)%zbar1 = zbar0(kk)

800   continue

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

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GLOBAL_PARAMETER_1.txt"
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

!Regional Variables Output
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/NEW/Regional_Variables.txt"
open(2000,file=trim(filename))

!Old Regional Variables Output
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/Regional_Variables.txt"
open(2001,file=trim(filename))

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

!Regional Variables Output
close(2000)

!Regional Variables Output (OLD)
close(2001)

end subroutine FILE_CLOSE

subroutine Write_Regional(i,REG)

  type(REGIONAL_template),intent(in) :: REG
  integer,intent(in) :: i

  !Write regional variables to file
  write(2000,*)i,REG

end subroutine Write_Regional

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

      subroutine rdveg(GLOBAL,CAT,GRID,REG,IO)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
      type (REGIONAL_template),intent(inout) :: REG
      type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
      type (IO_template),intent(inout) :: IO
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%ncol,GLOBAL%nrow)
      character(len=200) :: filename
      integer :: vegnvars,dvegnvars,ipos,jpos
      real,dimension(:,:,:),allocatable :: TempArray
      real*8 :: frcov(GLOBAL%nrow*GLOBAL%ncol,GLOBAL%ncatch+1),wc0
      integer :: jj,kk
      vegnvars = 20
      dvegnvars = 2
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,vegnvars))

! ====================================================================
! Read the image with the land cover clasifications.
! ====================================================================

      call rdimgi(GRID%VEG%ilandc,11,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)

! ====================================================================
! Read spatially constant vegetation parameters.
! ====================================================================

      GLOBAL%nlandc = GLOBAL%nrow*GLOBAL%ncol
      GLOBAL%iopveg = 0

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

      print*,"rdveg:  Reading in all the vegetation properties at once"

      read(1002,rec=1)TempArray(:,:,:)

      GRID_VEG_2D%ivgtyp = dble(TempArray(:,:,1))
      GRID_VEG_2D%xlai = dble(TempArray(:,:,2))
      GRID_VEG_2D%xlai_wsc = dble(TempArray(:,:,3))
      GRID_VEG_2D%albd = dble(TempArray(:,:,4))
      GRID_VEG_2D%albw = dble(TempArray(:,:,5))
      GRID_VEG_2D%emiss = dble(TempArray(:,:,6))
      GRID_VEG_2D%za = dble(TempArray(:,:,7))
      GRID_VEG_2D%zww = dble(TempArray(:,:,8))
      GRID_VEG_2D%z0m = dble(TempArray(:,:,9))
      GRID_VEG_2D%z0h = dble(TempArray(:,:,10))
      GRID_VEG_2D%zpd = dble(TempArray(:,:,11))
      GRID_VEG_2D%rsmin = dble(TempArray(:,:,12))
      GRID_VEG_2D%rsmax = dble(TempArray(:,:,13))
      GRID_VEG_2D%Rpl = dble(TempArray(:,:,14))
      GRID_VEG_2D%f3vpdpar = dble(TempArray(:,:,15))
      GRID_VEG_2D%f4temppar = dble(TempArray(:,:,16))
      GRID_VEG_2D%trefk = dble(TempArray(:,:,17))
      GRID_VEG_2D%tcbeta = dble(TempArray(:,:,18))
      GRID_VEG_2D%extinct = dble(TempArray(:,:,19))
      GRID_VEG_2D%canclos = dble(TempArray(:,:,20))

! ####################################################################
! Read the vegetation dynamic parameter file
! ####################################################################

      deallocate(TempArray)
      allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,dvegnvars))

      print*,"rdveg:  Reading in the dynamic vegetation properties"

      read(1003,rec=1)TempArray(:,:,:)

      GRID_VEG_2D%xlai = dble(TempArray(:,:,1))
      GRID_VEG_2D%albd = dble(TempArray(:,:,2))


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
                GRID(kk)%VEG%ivgtyp = GRID_VEG_2D(ipos,jpos)%ivgtyp
                GRID(kk)%VEG%emiss = GRID_VEG_2D(ipos,jpos)%emiss
                GRID(kk)%VEG%za = GRID_VEG_2D(ipos,jpos)%za
                GRID(kk)%VEG%zww = GRID_VEG_2D(ipos,jpos)%zww
                GRID(kk)%VEG%z0m = GRID_VEG_2D(ipos,jpos)%z0m
                GRID(kk)%VEG%z0h = GRID_VEG_2D(ipos,jpos)%z0h
                GRID(kk)%VEG%zpd = GRID_VEG_2D(ipos,jpos)%zpd
                GRID(kk)%VEG%rsmin = GRID_VEG_2D(ipos,jpos)%rsmin
                GRID(kk)%VEG%rsmax = GRID_VEG_2D(ipos,jpos)%rsmax
                GRID(kk)%VEG%Rpl = GRID_VEG_2D(ipos,jpos)%Rpl
                GRID(kk)%VEG%f3vpdpar = GRID_VEG_2D(ipos,jpos)%f3vpdpar
                GRID(kk)%VEG%f4temppar = GRID_VEG_2D(ipos,jpos)%f4temppar
                GRID(kk)%VEG%trefk = GRID_VEG_2D(ipos,jpos)%trefk
                GRID(kk)%VEG%xlai = GRID_VEG_2D(ipos,jpos)%xlai
                GRID(kk)%VEG%albd = GRID_VEG_2D(ipos,jpos)%albd
                GRID(kk)%VEG%tcbeta = exp(-0.5*GRID(kk)%VEG%xlai)
                GRID(kk)%VEG%xlai_wsc = GRID(kk)%VEG%xlai
                GRID(kk)%VEG%albw = GRID(kk)%VEG%albd !Move to its own file
                GRID(kk)%VEG%extinct = 0.00!VegData%extinct(ipos,jpos)
                GRID(kk)%VEG%canclos = 1.00!VegData%canclos(ipos,jpos)
                GRID(kk)%VEG%Tslope1 = 0.00!VegData%Tslope1(ipos,jpos)
                GRID(kk)%VEG%Tint1 = 0.00!VegData%Tint1(ipos,jpos)
                GRID(kk)%VEG%Tslope2 = 0.00!VegData%Tslope2(ipos,jpos)
                GRID(kk)%VEG%Tint2 = 0.00!VegData%Tint2(ipos,jpos)
                GRID(kk)%VEG%Twslope1 = 0.00!VegData%Twslope1(ipos,jpos)
                GRID(kk)%VEG%Twint1 = 0.00!VegData%Twint1(ipos,jpos)
                GRID(kk)%VEG%Twslope2 = 0.00!VegData%Twslope2(ipos,jpos)
                GRID(kk)%VEG%Twint2 = 0.00!VegData%Twint2(ipos,jpos)
                GRID(kk)%VEG%Tsep = 0.00!VegData%Tsep(ipos,jpos)
                GRID(kk)%VEG%Twsep = 0.00!VegData%Twsep(ipos,jpos)

        enddo

! ====================================================================
! Calculate parameters for each land cover type.
! ====================================================================

      do kk=1,GLOBAL%nlandc

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

      enddo

      print*,"rdveg:  Set minimum st. resist. and and wet can.stor.cap."

! ====================================================================
! Read the initial canopy storage from an image or as a constant
! as requested.
! ====================================================================

  read(1000,*) wc0

  do kk=1,GLOBAL%npix

    GRID(kk)%VARS%wcip1 = wc0

  enddo  

  print*,"rdveg:  Read initial wet canopy storages"


! ====================================================================
! Calculate the fraction of different cover types in each catchment
! and in total area.  First array index for icount and frcov is
! the land cover type, second array index is the catchment number.
! Catchment ncatch+1 is total area.
! ====================================================================

      do 550 kk=1,GLOBAL%nlandc

         do 540 jj=1,GLOBAL%ncatch

            frcov(kk,jj) = GLOBAL%pixsiz**2/CAT(jj)%area

540      continue

         frcov(kk,GLOBAL%ncatch+1) = 1/real(GLOBAL%npix)

550   continue

      print*,"rdveg:  Calculated fractional covers for vegetation"

! ====================================================================
! Find fraction of bare soil in each catchment.
! ====================================================================

      REG%fbsrg = zero

      do jj=1,GLOBAL%ncatch
        CAT(jj)%fbs = zero
      enddo

      do 570 jj=1,GLOBAL%ncatch+1

         !fbs(jj)  = zero

         do 560 kk=1,GLOBAL%nlandc

            if (GRID(kk)%VEG%ivgtyp.eq.0) then
               
               if (jj .eq. GLOBAL%ncatch+1) then
               
                  REG%fbsrg = REG%fbsrg + frcov(kk,jj)
        
               else

                  CAT(jj)%fbs = CAT(jj)%fbs + frcov(kk,jj)

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

      subroutine inisim(GLOBAL,IO,GRID)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
      type (IO_template),intent(inout) :: IO
      integer :: kk

! ====================================================================
! If initial root zone is not entered, then set the root
! zone soil moisture for calculation of thermal conductivity.
! This is changed later (initsm) to a new initial condition      
! based on Brooks-Corey and local water table depth.
! ====================================================================

         do 50 kk=1,GLOBAL%npix

            GRID(kk)%VARS%rzsm1 = GLOBAL%smpet0
            GRID(kk)%VARS%tzsm1 = GLOBAL%smpet0
            GRID(kk)%VARS%rzsm1_u = GLOBAL%smpet0
            GRID(kk)%VARS%tzsm1_u = GLOBAL%smpet0
            GRID(kk)%VARS%rzsm1_f = 0.d0
            GRID(kk)%VARS%tzsm1_f = 0.d0
            GRID(kk)%VARS%rzdthetaidt=0.d0
            GRID(kk)%VARS%tzdthetaidt=0.d0

50       continue

! ====================================================================
! Read data to tell how program will initialize the storm
! and interstorm event flags and times.
! ====================================================================

      read(1000,*) GLOBAL%iopflg

      if (GLOBAL%iopflg.eq.0) then

! --------------------------------------------------------------------
! If one event flag value is used then set all flags and times
! accordingly.
! --------------------------------------------------------------------

         read(1000,*) GLOBAL%istflg

         if (GLOBAL%istflg.eq.1) then

! ....................................................................
! If the event is a storm event.
! ....................................................................

            do 100 kk=1,GLOBAL%npix

               GRID(kk)%VARS%istorm = 1
               GRID(kk)%VARS%intstm = 0
               GRID(kk)%VARS%istmst = 0
               GRID(kk)%VARS%intstp = 0
               GRID(kk)%VARS%xintst = 0.0

100         continue

         else

! ....................................................................
! If the event is an interstorm event.
! ....................................................................

            do 200 kk=1,GLOBAL%npix

               GRID(kk)%VARS%istorm = 0
               GRID(kk)%VARS%intstm = 1
               GRID(kk)%VARS%istmst = 0
               GRID(kk)%VARS%intstp = 0
               GRID(kk)%VARS%xintst = GLOBAL%endstm

200         continue

         endif

      else

! --------------------------------------------------------------------
! Set event flags, time of events, cumulative values.
! --------------------------------------------------------------------

         do 300 kk=1,GLOBAL%npix

! ....................................................................
! For pixels under storm event.
! ....................................................................

            if (GRID(kk)%VARS%istorm.eq.1) then

               GRID(kk)%VARS%intstm = 0
               GRID(kk)%VARS%istmst = 0
               GRID(kk)%VARS%intstp = 0
               GRID(kk)%VARS%xintst = zero
               GRID(kk)%VARS%cuminf = zero

! ....................................................................
! For pixels under interstorm event.
! ....................................................................

            else

               GRID(kk)%VARS%intstm = 1
               GRID(kk)%VARS%istmst = 0
               GRID(kk)%VARS%intstp = 0
               GRID(kk)%VARS%xintst = GRID(kk)%VARS%intstp*GLOBAL%dt + GLOBAL%endstm
               
            endif

300      continue

      endif

! ====================================================================
! Initialize snow pack variables.
! ====================================================================

      do kk=1,GLOBAL%npix

         GRID(kk)%VARS%PackWater_us=zero
         GRID(kk)%VARS%SurfWater_us=zero
         GRID(kk)%VARS%Swq_us=zero
         GRID(kk)%VARS%VaporMassFlux_us=zero
         GRID(kk)%VARS%r_MeltEnergy_us=zero
         GRID(kk)%VARS%Outflow_us=zero
         GRID(kk)%VARS%PackWater=zero
         GRID(kk)%VARS%SurfWater=zero
         GRID(kk)%VARS%Swq=zero
         GRID(kk)%VARS%VaporMassFlux=zero
         GRID(kk)%VARS%r_MeltEnergy=zero
         GRID(kk)%VARS%Outflow=zero

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

      subroutine rdsoil(GLOBAL,CAT,GRID,IO)

      implicit none

      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
      type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
      type (IO_template),intent(inout) :: IO
      type (GRID_SOIL_template) :: GRID_SOIL_2D(GLOBAL%ncol,GLOBAL%nrow)
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%ncol,GLOBAL%nrow)
      integer :: soilnvars,ipos,jpos,jj,kk,nn
      integer :: icount(GLOBAL%ncol*GLOBAL%nrow,GLOBAL%ncatch+1)
      real,dimension(:,:,:),allocatable :: TempArray
      real*8 :: psic(GLOBAL%ncol*GLOBAL%nrow),tempsum,dtaken
      real*8 :: frsoil(GLOBAL%ncol*GLOBAL%nrow,GLOBAL%ncatch+1)
      soilnvars = 23

 allocate(TempArray(GLOBAL%ncol,GLOBAL%nrow,soilnvars))

! ====================================================================
! Read spatially constant bare soil parameters and GLOBAL.
! Then read root and transmission zone data.
! ====================================================================

      GLOBAL%nsoil = GLOBAL%nrow*GLOBAL%ncol

      read(1000,*)GLOBAL%irestype
      read(1000,*)GLOBAL%ikopt
      read(1000,*)GLOBAL%zrzmax
      read(1000,*)GLOBAL%iopsmini

      if (GLOBAL%iopsmini.eq.0) read(1000,*)GLOBAL%smpet0

      print*,"rdsoil:  Read spatially constant soil pars"

      if (GLOBAL%iopsmini.eq.1)&
         print*,"rdsoil:  Will read initial soil moisture images"

! ====================================================================
! Read the binary soil file
! ====================================================================

      print*,"rdsoil:  Reading in all soil properties at once"
      read(1001,rec=1)TempArray(:,:,:)
      GRID_SOIL_2D%bcbeta = dble(TempArray(:,:,1))
      GRID_SOIL_2D%psic = dble(TempArray(:,:,2))
      GRID_SOIL_2D%thetas = dble(TempArray(:,:,3))
      GRID_SOIL_2D%thetar = dble(TempArray(:,:,4))
      GRID_SOIL_2D%xk0 = dble(TempArray(:,:,5))
      GRID_SOIL_2D%zdeep = dble(TempArray(:,:,6))
      GRID_SOIL_2D%tdeep = dble(TempArray(:,:,7))
      GRID_SOIL_2D%zmid = dble(TempArray(:,:,8))
      GRID_SOIL_2D%tmid0 = dble(TempArray(:,:,9))
      GRID_SOIL_2D%rocpsoil = dble(TempArray(:,:,10))
      GRID_SOIL_2D%quartz = dble(TempArray(:,:,11))
      GRID_SOIL_2D%ifcoarse = TempArray(:,:,12)
      GRID_SOIL_2D%srespar1 = dble(TempArray(:,:,13))
      GRID_SOIL_2D%srespar2 = dble(TempArray(:,:,14))
      GRID_SOIL_2D%srespar3 = dble(TempArray(:,:,15))
      GRID_SOIL_2D%a_ice = dble(TempArray(:,:,16))
      GRID_SOIL_2D%b_ice = dble(TempArray(:,:,17))
      GRID_SOIL_2D%bulk_dens = dble(TempArray(:,:,18))
      GRID_SOIL_2D%amp = dble(TempArray(:,:,19))
      GRID_SOIL_2D%phase = dble(TempArray(:,:,20))
      GRID_SOIL_2D%shift = dble(TempArray(:,:,21))
      GRID_VEG_2D%tw = dble(TempArray(:,:,22))
      GRID_VEG_2D%tc = dble(TempArray(:,:,23))

! ====================================================================
! Read the soil classification image.
! ====================================================================

      call rdimgi(GRID%SOIL%isoil,12,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)

      print*,"rdsoil:  Read soil texture image"

! ====================================================================
! Pass the soil properties from the original i,j pos. to the kk pos.
!  ====================================================================
      do kk=1,GLOBAL%nsoil

         !Map the kk position to the i,j position
         if(mod(kk,GLOBAL%nrow) .ne. 0)then
                ipos = kk/GLOBAL%nrow+1
                jpos = mod(kk,GLOBAL%nrow)
         else
                ipos = kk/GLOBAL%nrow
                jpos = GLOBAL%nrow
         endif

         GRID(kk)%SOIL%bcbeta = GRID_SOIL_2D(ipos,jpos)%bcbeta
         GRID(kk)%SOIL%psic = GRID_SOIL_2D(ipos,jpos)%psic
         GRID(kk)%SOIL%thetas = GRID_SOIL_2D(ipos,jpos)%thetas
         GRID(kk)%SOIL%thetar = GRID_SOIL_2D(ipos,jpos)%thetar
         GRID(kk)%SOIL%xk0 = GRID_SOIL_2D(ipos,jpos)%xk0
         GRID(kk)%SOIL%zdeep = GRID_SOIL_2D(ipos,jpos)%zdeep
         GRID(kk)%SOIL%tdeep = GRID_SOIL_2D(ipos,jpos)%tdeep
         GRID(kk)%SOIL%zmid = GRID_SOIL_2D(ipos,jpos)%zmid
         GRID(kk)%SOIL%tmid0 = GRID_SOIL_2D(ipos,jpos)%tmid0
         GRID(kk)%SOIL%rocpsoil = GRID_SOIL_2D(ipos,jpos)%rocpsoil
         GRID(kk)%SOIL%quartz = GRID_SOIL_2D(ipos,jpos)%quartz
         GRID(kk)%SOIL%ifcoarse = GRID_SOIL_2D(ipos,jpos)%ifcoarse
         GRID(kk)%SOIL%srespar1 = GRID_SOIL_2D(ipos,jpos)%srespar1
         GRID(kk)%SOIL%srespar2 = GRID_SOIL_2D(ipos,jpos)%srespar2
         GRID(kk)%SOIL%srespar3 = GRID_SOIL_2D(ipos,jpos)%srespar3
         GRID(kk)%SOIL%a_ice = GRID_SOIL_2D(ipos,jpos)%a_ice
         GRID(kk)%SOIL%b_ice = GRID_SOIL_2D(ipos,jpos)%b_ice
         GRID(kk)%SOIL%bulk_dens = GRID_SOIL_2D(ipos,jpos)%bulk_dens
         GRID(kk)%SOIL%amp = GRID_SOIL_2D(ipos,jpos)%amp
         GRID(kk)%SOIL%phase = GRID_SOIL_2D(ipos,jpos)%phase
         GRID(kk)%SOIL%shift = GRID_SOIL_2D(ipos,jpos)%shift
         GRID(kk)%VEG%tc = GRID_VEG_2D(ipos,jpos)%tc
         GRID(kk)%VEG%tw = GRID_VEG_2D(ipos,jpos)%tw

      enddo
 
      GLOBAL%inc_frozen = 1 !THIS MEANS THAT THE FROZEN ALGORITHM IS ALWAYS RUN

! ====================================================================
! Calculate time in-variant soil parameters for each soil class.
! ====================================================================

      do 400 kk=1,GLOBAL%nsoil

! --------------------------------------------------------------------&
! Calculate soil parameters based on Brooks-Corey parameters.
! --------------------------------------------------------------------&

         GRID(kk)%SOIL%bcgamm = two + three * GRID(kk)%SOIL%bcbeta

! --------------------------------------------------------------------&
! Calculate constants for bare soil evaporation desorptivity 
! equation used in Famiglietti PhD Thesis, Princetion Univ, 1992.
! --------------------------------------------------------------------&

         GRID(kk)%SOIL%par = one + ((GRID(kk)%SOIL%bcgamm-one)/GRID(kk)%SOIL%bcbeta)
         GRID(kk)%SOIL%corr =((one/(one+three*GRID(kk)%SOIL%bcbeta))-&
                   (0.85d0/(one+four*GRID(kk)%SOIL%bcbeta))-&
                   (0.85d0*0.15d0*0.5d0/(one+five*GRID(kk)%SOIL%bcbeta))+&
                   (0.85d0*0.15d0*1.15d0/&
                   (six*(one+six*GRID(kk)%SOIL%bcbeta))))

! --------------------------------------------------------------------&
! Calculate diffusivity index and dimensionless exfiltration
! diffusivity from Eagleson, WRR, 1972.
! --------------------------------------------------------------------&

         GRID(kk)%SOIL%idifind = ((1.0+2.0*GRID(kk)%SOIL%bcbeta)/GRID(kk)%SOIL%bcbeta)+0.5
         tempsum=0 

         do 300 nn=1,GRID(kk)%SOIL%idifind

            dtaken = exp(factln(GRID(kk)%SOIL%idifind)-factln(nn)-&
                       factln(GRID(kk)%SOIL%idifind-nn))
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

      do 450 kk=1,GLOBAL%nsoil

         do 440 jj=1,GLOBAL%ncatch+1

            icount(kk,jj) = 0

440      continue

450   continue

      do 500 kk=1,GLOBAL%npix

        icount(GRID(kk)%SOIL%isoil,GRID(kk)%VARS%icatch)=icount(GRID(kk)%SOIL%isoil,GRID(kk)%VARS%icatch)+1
        icount(GRID(kk)%SOIL%isoil,GLOBAL%ncatch+1) = icount(GRID(kk)%SOIL%isoil,GLOBAL%ncatch+1) + 1

500   continue

      do 550 kk=1,GLOBAL%nsoil

         do 540 jj=1,GLOBAL%ncatch

            frsoil(kk,jj) = icount(kk,jj)*GLOBAL%pixsiz*GLOBAL%pixsiz/CAT(jj)%area

540      continue

         frsoil(kk,GLOBAL%ncatch+1) = icount(kk,GLOBAL%ncatch+1)/real(GLOBAL%npix)

550   continue

      print*,"rdsoil:  Calculated fractional coverage for soil types"

! ====================================================================
! Calculate average bubbling pressure in each catchment (used
! in updating average water table depths.
! ====================================================================

      do 570 jj=1,GLOBAL%ncatch

         CAT(jj)%psicav = zero

         do 560 kk=1,GLOBAL%nsoil

            CAT(jj)%psicav = CAT(jj)%psicav + frsoil(kk,jj)*psic(kk)

560      continue

570   continue

      print*,"rdsoil:  Calculated average psi! for each catchment"

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
      integer :: ia(nrow*ncol),ipixnum(nrow,ncol)
      integer :: nrow,ncol,iu,irow,icol,itmpval

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

      integer ipixnum(nrow,ncol)
      integer nrow,ncol,iu,irow,icol
      real*8 :: a(nrow*ncol)
      real*4 :: tmpval

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
      integer :: ipixnum(nrow,ncol),ixpix(nrow*ncol),iypix(nrow*ncol)
      integer :: nrow,ncol,npix
      integer :: ip,irow,icol
      real*8 :: atb(nrow*ncol)
      real*4 :: tempatb

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

      return

      end subroutine rdatb

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
      integer :: i,jj
      real*8 :: x,factln

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
