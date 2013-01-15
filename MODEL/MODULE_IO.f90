MODULE MODULE_IO

USE MODULE_VARIABLES

!Add the variables that are reused throughout the subroutines

implicit none

type IO_template
  integer,allocatable,dimension(:,:) :: ipixnum
  integer,allocatable,dimension(:) :: ixpix,iypix
end type IO_template

interface Extract_Info_General_File
  !====== begin of generated interface ======
  module procedure Extract_Info_General_File_Int
  module procedure Extract_Info_General_File_String
  module procedure Extract_Info_General_File_Double
end interface Extract_Info_General_File

contains

!###################################################################
! Subroutine to read in model parameters and initialize the variables
!###################################################################

subroutine Initialize_Model(GLOBAL,GRID,REG,CAT,IO)
  
  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (REGIONAL_template),intent(inout) :: REG
  type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
  type (IO_template),intent(inout) :: IO

!####################################################################
! Read the general filename
!####################################################################

  call Read_General_File(GLOBAL)

!####################################################################
! Open all files
!####################################################################

  call FILE_OPEN(GLOBAL)

!####################################################################
! Call rddata to open files, read in time in-variant parameters,&
! and initialize simulation sums.
!####################################################################

  call rddata(GLOBAL,GRID,REG,CAT,IO)

GRID%VARS%rzsm_f = 0.d0
GRID%VARS%rzsm1_f = 0.d0
GRID%VARS%tzsm_f = 0.d0
GRID%VARS%tzsm1_f = 0.d0
GRID%VARS%rzsm1 = zero
GRID%VARS%tzsm1 = zero
GRID%VARS%zw = 0.0d0
GRID%VARS%Sdepth_us = 0.0d0
GRID%VARS%Swq_us = 0.d0
  GRID%VARS%rnpet = 0.d0
  GRID%VARS%xlepet = 0.d0
  GRID%VARS%hpet = 0.d0
  GRID%VARS%gpet = 0.d0
  GRID%VARS%fw = 0.d0
  REG%zbar1rg = 0.d0
  REG%wcip1sum = 0.d0
  GRID%VEG%i_und = 0
GRID%VARS%PackWater_us = zero
GRID%VARS%SurfWater_us = zero
GRID%VARS%VaporMassFlux_us = zero
GRID%VARS%r_MeltEnergy_us = zero
GRID%VARS%Outflow_us = zero
GRID%VARS%alb_snow = zero
GRID%VARS%TPack = zero
GRID%VARS%TSurf = zero
GRID%VARS%xleact_snow = zero
GRID%VARS%hact_snow = zero
GRID%VARS%rn_snow = zero
GRID%VARS%TPack_us = zero
GRID%VARS%TSurf_us = zero
GRID%VARS%xleact_snow_us = zero
GRID%VARS%hact_snow_us = zero
GRID%VARS%rn_snow_us = zero
GRID%VARS%dens = zero
GRID%VARS%dens_us = zero
GRID%VARS%dsty_us = zero
GRID%VARS%Sdepth = zero
GRID%VARS%PackWater = zero
GRID%VARS%SurfWater = zero
GRID%VARS%Swq = zero
GRID%VARS%VaporMassFlux = zero
GRID%VARS%r_MeltEnergy = zero
GRID%VARS%Outflow = zero

end subroutine Initialize_Model

subroutine Read_Data(GLOBAL,GRID,CAT,IO,i)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
  type (IO_template),intent(inout) :: IO
  integer,intent(in) :: i

!#####################################################################
! Define the new water table
!#####################################################################

  CAT%zbar = CAT%zbar1

!#####################################################################
! Update the vegetation parameters if required.
!#####################################################################

  if (mod(i,GLOBAL%dtveg).eq.0) call rdveg_update(GLOBAL,GRID)

!#####################################################################
! Update the decimal Julian day.
!#####################################################################

  GLOBAL%djday = GLOBAL%djday + 0.0416666667d0*2.777777777d-4*GLOBAL%dt

!#####################################################################
! Read meteorological data.
!#####################################################################

  call rdatmo(i,GRID%MET,GLOBAL,IO)

end subroutine Read_Data

!#####################################################################
! Subroutine to output data to file
!#####################################################################

subroutine Write_Data(GLOBAL,GRID,IO,REG,i)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),intent(inout) :: GRID
  type (REGIONAL_template),intent(inout) :: REG
  type (IO_template),intent(inout) :: IO
  integer,intent(inout) :: i

!#####################################################################
! Output regional variables
!#####################################################################

  call Write_Regional(i,REG,GLOBAL)

!#####################################################################
! Output spatial field
!#####################################################################

  call Write_Binary(GRID%VARS%rzsm,1.0,GLOBAL%nrow,GLOBAL%ncol,&
                    IO%ipixnum,i,GLOBAL)

end subroutine Write_Data

!#####################################################################
! Finalize model and close files
!#####################################################################

subroutine Finalize_Model(GLOBAL)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL

!#####################################################################
! Close all files
!#####################################################################

  call FILE_CLOSE(GLOBAL)

end subroutine

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

! ####################################################################
! Read all variables in at once for each time step
! ####################################################################

      read(GLOBAL%FORCING_FILE%fp,rec=i) TempArray(:,:,:)

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
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable :: GRID
      type (REGIONAL_template) :: REG
      type (CATCHMENT_template),dimension(:),allocatable :: CAT
      type (IO_template),intent(inout) :: IO
      character(len=200) :: filename

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
      read(GLOBAL%DVEG_FILE%fp,rec=GLOBAL%ntdveg)TempArray
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

  call rdatb(atb,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum,IO%ixpix,IO%iypix,GLOBAL%npix,GLOBAL%TI_FILE%fp)

! ====================================================================
! Read in the catchment look-up table - read different values based
! on what is necessary for baseflow and initial condition calculations.
! ====================================================================

      if (GLOBAL%iopwt0.eq.1) then

         do 100 kk=1,GLOBAL%ncatch

            read(GLOBAL%CL_table_FILE%fp,*) jj,CAT(kk)%q0,CAT(kk)%ff,CAT(kk)%qb0,CAT(kk)%dd,&
                       CAT(kk)%xlength,CAT(kk)%basink

100      continue

      else if (GLOBAL%iopbf.eq.1) then

         do 200 kk=1,GLOBAL%ncatch

            read(GLOBAL%CL_table_FILE%fp,*) jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk),CAT(kk)%dd,&
                       CAT(kk)%xlength,CAT(kk)%basink 

200      continue

      else 

         do 300 kk=1,GLOBAL%ncatch

            read(GLOBAL%CL_table_FILE%fp,*) jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk)

300      continue

      endif

! ====================================================================
! Read the catchment image.
! ====================================================================

      GRID%VARS%icatch = 0
      call rdimgi(GRID%VARS%icatch,GLOBAL%Subbasin_FILE%fp,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)

! ====================================================================
! Read image of transmissivities for use in calculating the 
! soils-topographi! index.
! ====================================================================

      call rdimgr(ki,GLOBAL%K0_FILE%fp,GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum)
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

subroutine FILE_OPEN(GLOBAL)

implicit none
type (GLOBAL_template),intent(inout) :: GLOBAL
integer :: soilnvars,vegnvars,dvegnvars,nforcingvars,noutvars
integer :: nrow,ncol
soilnvars = 23
vegnvars = 20
noutvars = 1
dvegnvars = 2
nforcingvars = 7;
nrow = 60
ncol = 66

!Define the pointer number of each file
GLOBAL%TI_FILE%fp = 101
GLOBAL%Subbasin_FILE%fp = 102
GLOBAL%K0_FILE%fp = 103
GLOBAL%CL_Table_FILE%fp = 104
GLOBAL%VEG_FILE%fp = 105
GLOBAL%DVEG_FILE%fp = 106
GLOBAL%FORCING_FILE%fp = 107
GLOBAL%OUTPUT_FILE%fp = 108
GLOBAL%SOIL_FILE%fp = 109
GLOBAL%REGIONAL_FILE%fp = 110

!Open the files

!Topographic index file
open(unit=GLOBAL%TI_FILE%fp,file=GLOBAL%TI_FILE%fname,form='unformatted',access='direct',recl=4)

!Subbasin distribution file
open(unit=GLOBAL%Subbasin_FILE%fp,file=GLOBAL%Subbasin_FILE%fname,form='unformatted',access='direct',recl=4)

!Saturated hydraulic conductivity file
open(unit=GLOBAL%K0_FILE%fp,file=GLOBAL%K0_FILE%fname,form='unformatted',access='direct',recl=4)

!Catchment table parameters file
open(unit=GLOBAL%CL_table_FILE%fp,file=GLOBAL%CL_table_FILE%fname)

!Soil Parameter File
open(GLOBAL%SOIL_FILE%fp,file=trim(GLOBAL%SOIL_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=ncol*nrow*soilnvars*4)

!Vegetation Static Parameter File
open(GLOBAL%VEG_FILE%fp,file=trim(GLOBAL%VEG_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=ncol*nrow*vegnvars*4)

!Vegetation Dynamic Parameter File
open(GLOBAL%DVEG_FILE%fp,file=trim(GLOBAL%DVEG_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=ncol*nrow*dvegnvars*4)

!Forcing Data Set
open(GLOBAL%FORCING_FILE%fp,file=trim(GLOBAL%FORCING_FILE%fname),status='old',access='direct',&
     form='unformatted',recl=ncol*nrow*nforcingvars*4)

!Output Data Set
open(GLOBAL%OUTPUT_FILE%fp,file=trim(GLOBAL%OUTPUT_FILE%fname),status='unknown',access='direct',&
     form='unformatted',recl=ncol*nrow*noutvars*4)

!Regional Variables Output
open(GLOBAL%REGIONAL_FILE%fp,file=trim(GLOBAL%REGIONAL_FILE%fname))

end subroutine FILE_OPEN

subroutine FILE_CLOSE(GLOBAL)

implicit none
type(GLOBAL_template),intent(in) :: GLOBAL

!Soil Parameter File
close(GLOBAL%SOIL_FILE%fp)

!Vegetation Static Parameter File
close(GLOBAL%VEG_FILE%fp)

!Vegetation Dynamic Parameter File
close(GLOBAL%DVEG_FILE%fp)

!Forcing Data Set
close(GLOBAL%FORCING_FILE%fp)

!Output Data Set
close(GLOBAL%OUTPUT_FILE%fp)

!Regional Variables Output
close(GLOBAL%REGIONAL_FILE%fp)

!Topographic index
close(GLOBAL%TI_FILE%fp)

!Subbasin distribution file
close(GLOBAL%Subbasin_FILE%fp)

!Saturated hydraulic conductivity file
open(GLOBAL%K0_FILE%fp)

!Catchment table parameters file
open(GLOBAL%CL_table_FILE%fp)


end subroutine FILE_CLOSE

subroutine Write_Regional(i,REG,GLOBAL)

  type(REGIONAL_template),intent(in) :: REG
  type(GLOBAL_template),intent(in) :: GLOBAL
  integer,intent(in) :: i

  !print*,REG%fbsrg
  !print*,REG%Swqsum
  !print*,REG%Swq_ussum
  !print*,REG%Sdepthsum
  !print*,REG%Sdepth_ussum
  !print*,REG%fwreg
  !print*,REG%rzsmav
  !print*,REG%tzsmav
  !print*,REG%wcsum
  !print*,REG%wcip1sum
  !print*,REG%ettotrg
  !print*,REG%etstsumrg
  !print*,REG%etwtsumrg
  !print*,REG%etbssumrg
  !print*,REG%etdcsumrg
  !print*,REG%etwcsumrg
  !print*,REG%etlakesumrg
  !print*,REG%pptsumrg
  !print*,REG%pnetsumrg
  !print*,REG%contotrg
  !print*,REG%sxrtotrg
  !print*,REG%xixtotrg
  !print*,REG%qsurfrg
  !print*,REG%ranrunrg
  !print*,REG%conrunrg
  !print*,REG%qbreg
  !print*,REG%capsumrg
  !print*,REG%difrzsumrg
  !print*,REG%gwtsumrg
  !print*,REG%grzsumrg
  !print*,REG%gtzsumrg
  !print*,REG%zbarrg
  !print*,REG%zbar1rg
  !print*,REG%dswcsum
  !print*,REG%dsrzsum
  !print*,REG%dstzsum
  !print*,REG%dssum
  !print*,REG%wcrhssum
  !print*,REG%rzrhssum
  !print*,REG%tzrhssum
  !print*,REG%svarhssum
  !print*,REG%rnsum
  !print*,REG%xlesum
  !print*,REG%hsum
  !print*,REG%gsum
  !print*,REG%tksum
  !print*,REG%dshsum
  !print*,REG%tkmidsum
  !print*,REG%tkdeepsum
  !print*,REG%rnpetsum
  !print*,REG%xlepetsum
  !print*,REG%hpetsum
  !print*,REG%gpetsum
  !print*,REG%tkpetsum
  !print*,REG%tkmidpetsum
  !print*,REG%dshpetsum
  !print*,REG%perrg1
  !print*,REG%perrg2
  !print*,REG%pr3sat
  !print*,REG%pr2sat
  !print*,REG%pr2uns
  !print*,REG%pr1sat
  !print*,REG%pr1rzs
  !print*,REG%pr1tzs
  !print*,REG%pr1uns
  !print*,REG%persac
  !print*,REG%peruac
  !print*,REG%perusc
  !print*,REG%persxr
  !print*,REG%perixr
  write(GLOBAL%REGIONAL_FILE%fp,*)i,REG

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

      read(GLOBAL%VEG_FILE%fp,rec=1)TempArray(:,:,:)

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

      read(GLOBAL%DVEG_FILE%fp,rec=1)TempArray(:,:,:)

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

  do kk=1,GLOBAL%npix

    GRID(kk)%VARS%wcip1 = GLOBAL%wc0

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

      if (GLOBAL%iopflg.eq.0) then

! --------------------------------------------------------------------
! If one event flag value is used then set all flags and times
! accordingly.
! --------------------------------------------------------------------

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

      subroutine WRITE_BINARY(datain,rmult,nrow,ncol,ipixnum,i,GLOBAL)

      implicit none
      type(GLOBAL_template),intent(in) :: GLOBAL
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

      write(GLOBAL%OUTPUT_FILE%fp,rec=i) dataout

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

      print*,"rdsoil:  Read spatially constant soil pars"

      if (GLOBAL%iopsmini.eq.1)&
         print*,"rdsoil:  Will read initial soil moisture images"

! ====================================================================
! Read the binary soil file
! ====================================================================

      print*,"rdsoil:  Reading in all soil properties at once"
      read(GLOBAL%SOIL_FILE%fp,rec=1)TempArray(:,:,:)
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

      subroutine rdatb(atb,nrow,ncol,ipixnum,ixpix,iypix,npix,fp)

      implicit none
      integer :: ipixnum(nrow,ncol),ixpix(nrow*ncol),iypix(nrow*ncol)
      integer :: nrow,ncol,npix
      integer :: ip,irow,icol,fp
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


               read(fp,rec=((irow-1)*ncol) + icol) tempatb

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

subroutine Read_General_File(GLOBAL)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL

  !Read the general filename from the command line argument
  call get_command_argument(1,GLOBAL%GENERAL_FILE%fname)
  if (GLOBAL%GENERAL_FILE%fname .eq. '')then
    print*,'Mising general parameter file. Exiting model.'
    stop
  endif

  !Define the pointer to the file
  GLOBAL%GENERAL_FILE%fp = 100

  !Open the file
  open(GLOBAL%GENERAL_FILE%fp,file=trim(GLOBAL%GENERAL_FILE%fname),status='old')

  !Read through the file extracting the necessary information

  !Number of time steps
  call Extract_Info_General_File('ndata',GLOBAL,GLOBAL%ndata) 
  !Time step (seconds)
  call Extract_Info_General_File('dt',GLOBAL,GLOBAL%dt)
  !Time between end of precip event and end of storm event
  call Extract_Info_General_File('endstm',GLOBAL,GLOBAL%endstm)
  !Option for calculation of baseflow 
  call Extract_Info_General_File('iopbf',GLOBAL,GLOBAL%iopbf)
  !Option for initial condition specification 
  call Extract_Info_General_File('iopwt0',GLOBAL,GLOBAL%iopwt0)
  !Number of catchments
  call Extract_Info_General_File('ncatch',GLOBAL,GLOBAL%ncatch)
  !Number of rows in soils-TI image 
  call Extract_Info_General_File('nrow',GLOBAL,GLOBAL%nrow)
  !Number of columns in soils-TI image 
  call Extract_Info_General_File('ncol',GLOBAL,GLOBAL%ncol)
  !Pixel resolution for soils-TI imagex (m) 
  call Extract_Info_General_File('pixsiz',GLOBAL,GLOBAL%pixsiz)
  !Initial canopy water storage (m)
  call Extract_Info_General_File('wc0',GLOBAL,GLOBAL%wc0)
  !Type of soil resistance parameterization 
  call Extract_Info_General_File('irestype',GLOBAL,GLOBAL%irestype)
  !Option for vertical Ks change with depth 
  call Extract_Info_General_File('ikopt',GLOBAL,GLOBAL%ikopt)
  !Maximum root depth 
  call Extract_Info_General_File('zrzmax',GLOBAL,GLOBAL%zrzmax)
  !Initial conditions (yes/no)
  call Extract_Info_General_File('iopsmini',GLOBAL,GLOBAL%iopsmini)
  !Soil moisture used to calculate the thermal cond for t1 
  call Extract_Info_General_File('smpet0',GLOBAL,GLOBAL%smpet0)
  !Option for ground heat flux under vegetation 
  call Extract_Info_General_File('iopgveg',GLOBAL,GLOBAL%iopgveg)
  !Option for thermal conductivity calculation 
  call Extract_Info_General_File('iopthermc',GLOBAL,GLOBAL%iopthermc)
  !Adjust soil thermal conductivity based on LAI 
  call Extract_Info_General_File('iopthermc_v',GLOBAL,GLOBAL%iopthermc_v)
  !Maximum number of iterations for energy balance solution 
  call Extract_Info_General_File('maxnri',GLOBAL,GLOBAL%maxnri)
  !Tolerance for skin temperature (K) 
  call Extract_Info_General_File('toleb',GLOBAL,GLOBAL%toleb)
  !Option for the initial input of storm/interstorm event status  
  call Extract_Info_General_File('iopflg',GLOBAL,GLOBAL%iopflg)
  !Initial storm/interstorm event flag 
  call Extract_Info_General_File('istflg',GLOBAL,GLOBAL%istflg)
  !Number of time steps between updating vegetation parameters 
  call Extract_Info_General_File('dtveg',GLOBAL,GLOBAL%dtveg)
  !Topographic index filename
  call Extract_Info_General_File('TI_fname',GLOBAL,GLOBAL%TI_FILE%fname)
  !Basin distribution filename
  call Extract_Info_General_File('Subbasin_fname',GLOBAL,GLOBAL%Subbasin_FILE%fname)
  !Saturated Hydraulic Conductivity filename
  call Extract_Info_General_File('K_0_fname',GLOBAL,GLOBAL%K0_FILE%fname)
  !Catchment parameter table
  call Extract_Info_General_File('CL_Table_fname',GLOBAL,GLOBAL%CL_table_FILE%fname)
  !Soil file
  call Extract_Info_General_File('SOIL_fname',GLOBAL,GLOBAL%SOIL_FILE%fname)
  !Vegetation file
  call Extract_Info_General_File('VEG_fname',GLOBAL,GLOBAL%VEG_FILE%fname)
  !Dynamic vegetation file
  call Extract_Info_General_File('DVEG_fname',GLOBAL,GLOBAL%DVEG_FILE%fname)
  !Forcing file
  call Extract_Info_General_File('FORCING_fname',GLOBAL,GLOBAL%FORCING_FILE%fname)
  !Output file
  call Extract_Info_General_File('OUTPUT_fname',GLOBAL,GLOBAL%OUTPUT_FILE%fname)
  !Output file
  call Extract_Info_General_File('REGIONAL_fname',GLOBAL,GLOBAL%REGIONAL_FILE%fname)
  !Number of threads used in openmp
  call Extract_Info_General_File('nthreads',GLOBAL,GLOBAL%nthreads)
  
  !Close the file
  close(GLOBAL%GENERAL_FILE%fp)
  
end subroutine

subroutine Extract_Info_General_File_Int(strid,GLOBAL,read_arg)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  integer :: read_arg,flag
  character(*) :: strid
  character(len=200) :: string
  !Read until we find what we need
  do 
    read(GLOBAL%GENERAL_FILE%fp,*,iostat=flag)string
    if (string .eq. strid)then
      !go back one record
      backspace(GLOBAL%GENERAL_FILE%fp)
      !read the desired parameter/file
      read(GLOBAL%GENERAL_FILE%fp,*)string,read_arg
      !Go back to the beginning of the file and exit
      rewind(GLOBAL%GENERAL_FILE%fp)
      exit
    endif
    if (flag .ne. 0)then
      print*,'Missing the ',trim(strid),' input parameter. Check the',&
             ' General Parameter file' 
      stop
    endif
  enddo

end subroutine Extract_Info_General_File_Int

subroutine Extract_Info_General_File_String(strid,GLOBAL,read_arg)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  integer :: flag
  character(len=*) :: strid
  character(len=400) :: read_arg
  character(len=200) :: string
  !Read until we find what we need
  do
    read(GLOBAL%GENERAL_FILE%fp,*,iostat=flag)string
    if (string .eq. strid)then
      !go back one record
      backspace(GLOBAL%GENERAL_FILE%fp)
      !read the desired parameter/file
      read(GLOBAL%GENERAL_FILE%fp,*)string,read_arg
      !Go back to the beginning of the file and exit
      rewind(GLOBAL%GENERAL_FILE%fp)
      exit
    endif
    if (flag .ne. 0)then
      print*,'Missing the ',trim(strid),' input parameter. Check the',& 
             ' General Parameter file'
      stop
    endif
  enddo

end subroutine Extract_Info_General_File_String

subroutine Extract_Info_General_File_Double(strid,GLOBAL,read_arg)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  integer :: flag
  real*8 :: read_arg
  character(len=*) :: strid
  character(len=200) :: string
  !Read until we find what we need
  do
    read(GLOBAL%GENERAL_FILE%fp,*,iostat=flag)string
    if (string .eq. strid)then
      !go back one record
      backspace(GLOBAL%GENERAL_FILE%fp)
      !read the desired parameter/file
      read(GLOBAL%GENERAL_FILE%fp,*)string,read_arg
      !Go back to the beginning of the file and exit
      rewind(GLOBAL%GENERAL_FILE%fp)
      exit
    endif
    if (flag .ne. 0)then
      print*,'Missing the ',trim(strid),' input parameter. Check the',&
             ' General Parameter file'
      stop
    endif
  enddo

end subroutine Extract_Info_General_File_Double

END MODULE MODULE_IO
