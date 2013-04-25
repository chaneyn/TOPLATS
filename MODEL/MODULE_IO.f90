MODULE MODULE_IO

USE MODULE_VARIABLES

USE MODULE_TOPMODEL

USE NETCDF

!Add the variables that are reused throughout the subroutines

implicit none

type IO_template
  integer,allocatable,dimension(:,:) :: ipixnum
end type IO_template

interface Extract_Info_General_File
  !====== begin of generated interface ======
  module procedure Extract_Info_General_File_Int
  module procedure Extract_Info_General_File_String
  module procedure Extract_Info_General_File_Double
end interface Extract_Info_General_File

interface Convert_Grads2Model
  !====== begin of generated interface ======
  module procedure convert_grads2model_double
  module procedure convert_grads2model_int
end interface Convert_Grads2Model

contains

!###################################################################
!> Subroutine to read in model parameters and initialize the variables
!###################################################################

subroutine Initialize_Model(GLOBAL,GRID,REG,CAT,IO)
  
  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (REGIONAL_template),intent(inout) :: REG
  type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
  type (IO_template),intent(inout) :: IO
  integer :: i

!####################################################################
! Read the general filename
!####################################################################

  call Read_General_File(GLOBAL)
  GLOBAL%iopflg = 0
  GLOBAL%istflg = 0

!####################################################################
! Open all files
!####################################################################

  call FILE_OPEN(GLOBAL)

!####################################################################
! Call rddata to open files, read in time in-variant parameters,&
! and initialize simulation sums.
!####################################################################

  call rddata(GLOBAL,GRID,REG,CAT,IO)

!####################################################################
! Initialize necessary grid values to 0
!####################################################################

  do i = 1,GLOBAL%npix
    GRID(i)%VARS%sm_f = zero
    GRID(i)%VARS%sm1_f = zero
    GRID(i)%VARS%sm1 = zero
  enddo
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
  GRID%VARS%precip_u = zero

end subroutine Initialize_Model

!###################################################################
!> Subroutine to read in the model parameters
!###################################################################

subroutine Read_Data(GLOBAL,GRID,CAT,IO,i)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
  type (IO_template),intent(inout) :: IO
  integer,intent(in) :: i
  integer :: date(9)
  integer :: itime,time

!#####################################################################
! Define the new water table
!#####################################################################

  CAT%zbar = CAT%zbar1

!#####################################################################
! Find the date information for the current time step
!#####################################################################

  GLOBAL%time = GLOBAL%itime + GLOBAL%dt*(i-1)
  GLOBAL%old_date = GLOBAL%new_date
  call gmtime(GLOBAL%time,GLOBAL%new_date)

!#####################################################################
! Update the vegetation parameters when required.
!#####################################################################

  if (GLOBAL%dynamic_vegetation .eq. 0)then
    !Monthly Climatology
    if (GLOBAL%new_date(5) .ne. GLOBAL%old_date(5))then
      GLOBAL%ntdveg = GLOBAL%ntdveg + 1
      if (GLOBAL%ntdveg .gt. 12)GLOBAL%ntdveg = 1
      call rdveg_update(GLOBAL,GRID,IO)
    endif

  elseif (GLOBAL%dynamic_vegetation .eq. 1)then
    !Dynamic update
    if (mod(i,GLOBAL%dtveg).eq.0)then
      GLOBAL%ntdveg = GLOBAL%ntdveg + 1
      call rdveg_update(GLOBAL,GRID,IO)
    endif
  endif

!#####################################################################
! Read meteorological data.
!#####################################################################

  call rdatmo(i,GRID%MET,GLOBAL,IO)

end subroutine Read_Data

!#####################################################################
!> Subroutine to output data to file
!#####################################################################

subroutine Write_Data(GLOBAL,GRID,IO,REG,i,CAT)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (GRID_template),dimension(:),intent(inout) :: GRID
  type (REGIONAL_template),intent(inout) :: REG
  type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
  type (IO_template),intent(inout) :: IO
  integer,intent(in) :: i

!#####################################################################
! Output regional variables
!#####################################################################

  call Write_Regional(i,REG,GLOBAL)

!#####################################################################
! Output catchment variables
!#####################################################################

  call Write_Catchment(i,CAT,GLOBAL)

!#####################################################################
! Output spatial field
!#####################################################################

  !Root Zone Soil Moisture
  call WRITE_NETCDF(GRID%VARS%rzsm_zrzmax,IO%ipixnum,i,GLOBAL,GLOBAL%NETCDF_OUTPUT_FILE%varid(1))
  !Evaportranspiration
  call WRITE_NETCDF(GRID%VARS%etpix,IO%ipixnum,i,GLOBAL,GLOBAL%NETCDF_OUTPUT_FILE%varid(2))
  !Surface Runoff
  call WRITE_NETCDF(GRID%VARS%runtot,IO%ipixnum,i,GLOBAL,GLOBAL%NETCDF_OUTPUT_FILE%varid(3))

end subroutine Write_Data

!#####################################################################
!> Finalize model and close files
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
!> Subroutine to read in and pass meteorological data (e.g. rainfall).
!####################################################################

      subroutine rdatmo(i,MET,GLOBAL,IO)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_MET_template),intent(inout) :: MET(GLOBAL%nrow*GLOBAL%ncol)
      type (IO_template),intent(in) :: IO
      integer :: ipixnum(GLOBAL%nrow,GLOBAL%ncol)
      integer :: forcingnvars,i,ilat,ilon,ip
      real,dimension(:,:,:),allocatable :: TempArray
      ipixnum = IO%ipixnum
      forcingnvars = 7
      allocate(TempArray(GLOBAL%FORCING_FILE%nlon,GLOBAL%FORCING_FILE%nlat,forcingnvars))
     

! ####################################################################
! Read all variables in at once for each time step
! ####################################################################

  read(GLOBAL%FORCING_FILE%fp,rec=i) TempArray(:,:,:)

  !Set all missing values to the areal average
  call Replace_Undefined(TempArray,GLOBAL%Forcing_File%undef,GLOBAL%Forcing_File%nlon,&
                   GLOBAL%Forcing_File%nlat,forcingnvars)

  do ilat = 1,GLOBAL%nrow
    do ilon = 1,GLOBAL%ncol 
      ip = ipixnum(ilat,ilon)
      if (ip .gt. 0)then
        !Air temperature
        MET(ip)%tdry = TempArray(MET(ip)%tdry_MAP%ilon,MET(ip)%tdry_MAP%ilat,5)
        !Longwave radiation
        MET(ip)%rld = TempArray(MET(ip)%rld_MAP%ilon,MET(ip)%rld_MAP%ilat,1)
        !Air pressure
        MET(ip)%press = TempArray(MET(ip)%press_MAP%ilon,MET(ip)%press_MAP%ilat,2)
        !Relative humidity
        MET(ip)%rh = TempArray(MET(ip)%rh_MAP%ilon,MET(ip)%rh_MAP%ilat,3)
        !Downward shortwave radiation
        MET(ip)%rsd = TempArray(MET(ip)%rsd_MAP%ilon,MET(ip)%rsd_MAP%ilat,4)
        !Wind speed
        MET(ip)%uzw = TempArray(MET(ip)%uzw_MAP%ilon,MET(ip)%uzw_MAP%ilat,6)
        !Precipitation
        MET(ip)%pptms = TempArray(MET(ip)%pptms_MAP%ilon,MET(ip)%pptms_MAP%ilat,7)
      endif
    enddo
  enddo

      end subroutine

!>Subroutine to replace all undefined values with the areal average for an array
subroutine Replace_Undefined(array,undef,nlon,nlat,nvars)

  implicit none
  integer,intent(in) :: nlat,nlon,nvars
  real*4,intent(inout) :: array(nlon,nlat,nvars)
  real*8,intent(in) :: undef
  integer :: ivar,val_count,ilon,ilat
  real*8 :: val

  !Replace undefined values with the areal average for every variable
  do ivar = 1,nvars
    !Find average
    val = 0.0d0
    val_count = 0
    do ilon=1,nlon
      do ilat=1,nlat
        if (int(array(ilon,ilat,ivar)) .ne. int(undef))then
          val = val + array(ilon,ilat,ivar)
          val_count = val_count + 1
        endif
      enddo
    enddo
    val = val/val_count
    !Replace all missing values with the average value
    where (int(array(:,:,ivar)) .eq. int(undef))
      array(:,:,ivar) = val
    endwhere
  enddo

end subroutine Replace_Undefined


! ####################################################################
!> Subroutine to open input/output files, read and initialize time in-variant data.
! ####################################################################

      subroutine rddata(GLOBAL,GRID,REG,CAT,IO)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable :: GRID
      type (REGIONAL_template) :: REG
      type (CATCHMENT_template),dimension(:),allocatable :: CAT
      type (IO_template),intent(inout) :: IO
      character(len=200) :: filename
      integer :: i

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

      return

      end subroutine


! ====================================================================
!
!			subroutine rdveg_update
!
! ====================================================================
!
!>Subroutine to update the vegetation parameters.
!
! ====================================================================

      subroutine rdveg_update (GLOBAL,GRID,IO)

      implicit none
      integer :: kk,jj
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),intent(inout) :: GRID
      type (IO_template),intent(in) :: IO
      integer :: dvegnvars,ipos,jpos,ip,ilat,ilon
      real,dimension(:,:,:),allocatable :: TempArray
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%DVEG_FILE%nlon,GLOBAL%DVEG_FILE%nlat)
      dvegnvars = 2

! ====================================================================
! Read lookup table with parameters for vegetation/land cover.
! Read either critical and wilting soil moistures or
! plant/root resistance parameter depending on actual
! transpiration option.
! ====================================================================

      allocate(TempArray(GLOBAL%DVEG_FILE%nlon,GLOBAL%DVEG_FILE%nlat,dvegnvars))
      read(GLOBAL%DVEG_FILE%fp,rec=GLOBAL%ntdveg)TempArray
  !Set all missing values to the areal average
  call Replace_Undefined(TempArray,GLOBAL%DVEG_File%undef,GLOBAL%DVEG_File%nlon,&
                   GLOBAL%DVEG_File%nlat,dvegnvars)

      GRID_VEG_2D%xlai = dble(TempArray(:,:,1))
      GRID_VEG_2D%albd = dble(TempArray(:,:,2))

! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

  do ilat = 1,GLOBAL%nrow
    do ilon = 1,GLOBAL%ncol
      ip = IO%ipixnum(ilat,ilon)
      if (ip .gt. 0)then
        !Dynamic Vegetation
        GRID(ip)%VEG%xlai = GRID_VEG_2D(GRID(ip)%VEG%dynamic_MAP%ilon,GRID(ip)%VEG%dynamic_MAP%ilat)%xlai
        GRID(ip)%VEG%albd = GRID_VEG_2D(GRID(ip)%VEG%dynamic_MAP%ilon,GRID(ip)%VEG%dynamic_MAP%ilat)%albd
      endif
    enddo
  enddo
  GRID%VEG%tcbeta = exp(-0.5*GRID%VEG%xlai)
  GRID%VEG%xlai_wsc = GRID%VEG%xlai
  GRID%VEG%albw = GRID%VEG%albd

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
!> Subroutine to read in and initialize topmodel parameters and the soils-topographic index map.
!#####################################################################

subroutine rdtpmd(GRID,CAT,IO,GLOBAL)

  implicit none
  type (GLOBAL_template) :: GLOBAL
  type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
  type (CATCHMENT_template),dimension(:),allocatable :: CAT
  type (IO_template),intent(inout) :: IO
  integer kk,jj,ilat,ilon,ip
  real*8 hbar0
  integer,dimension(:),allocatable :: icount
  real*8,dimension(:),allocatable :: atb,ti,zbar0,sumatb,sumlti,qb0,lte
  real*8,dimension(:),allocatable :: ki
  real*8,dimension(:,:),allocatable :: array_2d
  real*8,dimension(:),allocatable :: array_1d
  real*4,dimension(:,:),allocatable :: temp

  ! Set the global grid information to the same as the topographic index
  GLOBAL%minlat = GLOBAL%TI_FILE%minlat
  GLOBAL%minlon = GLOBAL%TI_FILE%minlon
  GLOBAL%spatial_res = GLOBAL%TI_FILE%spatial_res
  GLOBAL%nrow = GLOBAL%TI_FILE%nlat
  GLOBAL%ncol = GLOBAL%TI_FILE%nlon

  !#####################################################################
  ! Allocate memory
  !#####################################################################

  allocate(CAT(GLOBAL%ncatch))
  allocate(GRID(GLOBAL%nrow*GLOBAL%ncol))
  allocate(IO%ipixnum(GLOBAL%nrow,GLOBAL%ncol))
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

  call rdatb(atb,GLOBAL,IO)
  GRID%VARS%TI = atb

! ====================================================================
! Read in the catchment look-up table - read different values based
! on what is necessary for baseflow and initial condition calculations.
! ====================================================================

    if (GLOBAL%KS_TYPE.eq.0)then

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

            read(GLOBAL%CL_table_FILE%fp,*)jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk)

300      continue

      endif
  
    else if (GLOBAL%KS_TYPE.eq.1)then

      do kk=1,GLOBAL%ncatch

        read(GLOBAL%CL_table_FILE%fp,*)jj,CAT(kk)%q0,CAT(kk)%ff,zbar0(kk),CAT(kk)%n

      enddo
 
    endif
  ! Map other variables using nni to the topographic index grid
  !Air temperature
  call spatial_mapping(GLOBAL,GRID%MET%tdry_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Precipitation
  call spatial_mapping(GLOBAL,GRID%MET%pptms_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Downward Longwave Radiation
  call spatial_mapping(GLOBAL,GRID%MET%rld_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Downward Shortwave Radiation
  call spatial_mapping(GLOBAL,GRID%MET%rsd_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Pressure
  call spatial_mapping(GLOBAL,GRID%MET%press_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Mean Wind Speed
  call spatial_mapping(GLOBAL,GRID%MET%uzw_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Relative Humidity
  call spatial_mapping(GLOBAL,GRID%MET%rh_MAP,GLOBAL%FORCING_FILE,IO%ipixnum)
  !Static Vegetation Properties
  call spatial_mapping(GLOBAL,GRID%VEG%static_MAP,GLOBAL%VEG_FILE,IO%ipixnum)
  !Dynamic Vegetation Properties
  call spatial_mapping(GLOBAL,GRID%VEG%dynamic_MAP,GLOBAL%DVEG_FILE,IO%ipixnum)
  !Soil Properties
  call spatial_mapping(GLOBAL,GRID%SOIL%MAP,GLOBAL%SOIL_FILE,IO%ipixnum)
  !Saturated Hydraulic Conductivity
  call spatial_mapping(GLOBAL,GRID%SOIL%K0_MAP,GLOBAL%K0_FILE,IO%ipixnum)

! ====================================================================
! Read the catchment image.
! ====================================================================

  GRID%VARS%icatch = 0
  allocate(array_2d(GLOBAL%Subbasin_FILE%nlon,GLOBAL%Subbasin_FILE%nlat))
  allocate(temp(GLOBAL%Subbasin_FILE%nlon,GLOBAL%Subbasin_FILE%nlat))
  allocate(array_1d(GLOBAL%Subbasin_FILE%nlon*GLOBAL%Subbasin_FILE%nlat))

  read(GLOBAL%Subbasin_FILE%fp,rec=1)temp

  ! Convert from single to double point precision
  array_2d = temp

  ! Convert to model format
  call convert_grads2model(array_1d,array_2d,IO%ipixnum,GLOBAL%nrow,GLOBAL%ncol,GLOBAL%Subbasin_FILE%undef)
  GRID%VARS%icatch = int(array_1d)

  deallocate(array_2d)
  deallocate(temp)
  deallocate(array_1d)
      
! ====================================================================
! Read image of transmissivities for use in calculating the 
! soils-topographi! index.
! ====================================================================

  allocate(array_2d(GLOBAL%K0_FILE%nlon,GLOBAL%K0_FILE%nlat))
  allocate(temp(GLOBAL%K0_FILE%nlon,GLOBAL%K0_FILE%nlat))
  allocate(array_1d(GLOBAL%K0_FILE%nlon*GLOBAL%K0_FILE%nlat))

  read(GLOBAL%K0_FILE%fp,rec=1)temp
  !Set all missing values to the areal average
  call Replace_Undefined(temp,GLOBAL%K0_File%undef,GLOBAL%K0_File%nlon,&
                   GLOBAL%K0_File%nlat,1)

  ! Convert from single to double point precision
  array_2d = temp
  ki = 0.d0

  do ilat = 1,GLOBAL%nrow
    do ilon = 1,GLOBAL%ncol
      ip = IO%ipixnum(ilat,ilon)
      if (ip .gt. 0)then
        !Static Soil
        ki(ip)=array_2d(GRID(ip)%SOIL%K0_MAP%ilon,GRID(ip)%SOIL%K0_MAP%ilat)
      endif
    enddo
  enddo

  do kk=1,GLOBAL%nrow*GLOBAL%ncol
    if (GRID(kk)%VARS%icatch .gt. 0)then
      ti(kk) = ki(kk)/CAT(GRID(kk)%VARS%icatch)%ff
    else
      ti(kk) = 0.0
    endif
  enddo 
  GRID%VARS%T0 = ti

  deallocate(array_2d)
  deallocate(temp)
  deallocate(array_1d)


! ====================================================================
! Calculate the generalized soils-topographic index for each catchment
! ====================================================================

  call Calculate_GSTI(GLOBAL,CAT,GRID)

! ====================================================================
! Calculate the average topographi! index value, the average of 
! the natural log of transmissivity and area for each basin.
! ====================================================================

      do 400 kk=1,GLOBAL%ncatch

         sumatb(kk) = 0.0
         sumlti(kk) = 0.0
         icount(kk) = 0
         CAT(GRID(kk)%VARS%icatch)%lambda = zero

400   continue

      do 500 kk=1,GLOBAL%npix

         sumatb(GRID(kk)%VARS%icatch) = sumatb(GRID(kk)%VARS%icatch) + atb(kk)
         sumlti(GRID(kk)%VARS%icatch) = sumlti(GRID(kk)%VARS%icatch) + dlog(ti(kk))
         CAT(GRID(kk)%VARS%icatch)%lambda = CAT(GRID(kk)%VARS%icatch)%lambda + GRID(kk)%VARS%GSTI
         icount(GRID(kk)%VARS%icatch) = icount(GRID(kk)%VARS%icatch) + 1

500   continue

      do 600 kk=1,GLOBAL%ncatch

         if (icount(kk).eq.0) then

            CAT(kk)%xlamda = 0
            lte(kk) = 0
            CAT(kk)%lambda = zero

         else

            CAT(kk)%xlamda = sumatb(kk)/icount(kk)
            lte(kk) = sumlti(kk)/icount(kk)
            CAT(kk)%lambda = CAT(kk)%lambda/icount(kk)

         endif

         CAT(kk)%area = icount(kk)*GLOBAL%pixsiz*GLOBAL%pixsiz

600   continue

  !if (GLOBAL%KS_TYPE .eq. 1)then

!###################################################################################
! Calculate the initial baseflow using the generalized topmodel (Chaney et al.,2013)
!###################################################################################

   ! do kk = 1,GLOBAL%ncatch

   !   CAT(kk)%q0 = CAT(kk)%area/(CAT(kk)%lambda**CAT(kk)%n)

   ! enddo

    print*,'Area',CAT%area
    print*,'Lambda',CAT%lambda
    print*,'Q0',CAT%q0
    print*,'te',exp(lte)
    print*,'xlamda',CAT%xlamda
    print*,'Q0old',CAT%area*dexp(lte)*dexp(-CAT%xlamda)

  !endif

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
!> Subroutine to handle opening and closing the input files
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
nrow = GLOBAL%nrow
ncol = GLOBAL%ncol

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
GLOBAL%CATCHMENT_FILE%fp = 111
GLOBAL%GSTI_FILE%fp = 112
GLOBAL%QS_FILE%fp = 113
GLOBAL%ET_FILE%fp = 113

!Open the files

!Topographic index file
open(unit=GLOBAL%TI_FILE%fp,file=GLOBAL%TI_FILE%fname,form='unformatted',&
     access='direct',recl=4*GLOBAL%TI_FILE%nlat*GLOBAL%TI_FILE%nlon)

!Subbasin distribution file
open(unit=GLOBAL%Subbasin_FILE%fp,file=GLOBAL%Subbasin_FILE%fname,form='unformatted',&
     access='direct',recl=4*GLOBAL%Subbasin_FILE%nlat*GLOBAL%Subbasin_FILE%nlon)

!Saturated hydraulic conductivity file
open(unit=GLOBAL%K0_FILE%fp,file=GLOBAL%K0_FILE%fname,form='unformatted',&
     access='direct',recl=4*GLOBAL%K0_FILE%nlat*GLOBAL%K0_FILE%nlon)

!Catchment table parameters file
open(unit=GLOBAL%CL_table_FILE%fp,file=GLOBAL%CL_table_FILE%fname)

!Soil Parameter File
open(GLOBAL%SOIL_FILE%fp,file=trim(GLOBAL%SOIL_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=GLOBAL%SOIL_FILE%nlon*GLOBAL%SOIL_FILE%nlat*soilnvars*4)

!Vegetation Static Parameter File
open(GLOBAL%VEG_FILE%fp,file=trim(GLOBAL%VEG_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=GLOBAL%VEG_FILE%nlon*GLOBAL%VEG_FILE%nlat*vegnvars*4)

!Vegetation Dynamic Parameter File
open(GLOBAL%DVEG_FILE%fp,file=trim(GLOBAL%DVEG_FILE%fname),status='old',access='direct',form='unformatted',&
     recl=GLOBAL%DVEG_FILE%nlon*GLOBAL%DVEG_FILE%nlat*dvegnvars*4)

!Forcing Data Set
open(GLOBAL%FORCING_FILE%fp,file=trim(GLOBAL%FORCING_FILE%fname),status='old',access='direct',&
     form='unformatted',recl=GLOBAL%FORCING_FILE%nlon*GLOBAL%FORCING_FILE%nlat*nforcingvars*4)

!Output Data Set
open(GLOBAL%OUTPUT_FILE%fp,file=trim(GLOBAL%OUTPUT_FILE%fname),status='unknown',access='direct',&
     form='unformatted',recl=GLOBAL%TI_FILE%nlon*GLOBAL%TI_FILE%nlat*noutvars*4)

!Output Data Set
open(GLOBAL%ET_FILE%fp,file=trim(GLOBAL%ET_FILE%fname),status='unknown',access='direct',&
     form='unformatted',recl=GLOBAL%TI_FILE%nlon*GLOBAL%TI_FILE%nlat*noutvars*4)

!Output Data Set
open(GLOBAL%QS_FILE%fp,file=trim(GLOBAL%QS_FILE%fname),status='unknown',access='direct',&
     form='unformatted',recl=GLOBAL%TI_FILE%nlon*GLOBAL%TI_FILE%nlat*noutvars*4)

!Regional Variables Output
open(GLOBAL%REGIONAL_FILE%fp,file=trim(GLOBAL%REGIONAL_FILE%fname))

!Catchment Variables Output
open(GLOBAL%CATCHMENT_FILE%fp,file=trim(GLOBAL%CATCHMENT_FILE%fname))

!GSTI Output
open(GLOBAL%GSTI_FILE%fp,file=trim(GLOBAL%GSTI_FILE%fname),status='unknown',form='unformatted',&
     access='direct',recl=4*GLOBAL%TI_FILE%nlat*GLOBAL%TI_FILE%nlon)

!NETCDF Output
GLOBAL%NETCDF_OUTPUT_FILE%var_name(1) = 'RZSM'
GLOBAL%NETCDF_OUTPUT_FILE%var_name(2) = 'ET'
GLOBAL%NETCDF_OUTPUT_FILE%var_name(3) = 'QSURF'
call Create_Netcdf_Output(GLOBAL%NETCDF_OUTPUT_FILE,GLOBAL,3)


end subroutine FILE_OPEN

!>Subroutine to close all open files
subroutine FILE_CLOSE(GLOBAL)

implicit none
type(GLOBAL_template),intent(in) :: GLOBAL
integer :: status

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
close(GLOBAL%K0_FILE%fp)

!Catchment table parameters file
close(GLOBAL%CL_table_FILE%fp)

!Netcdf OUTPUT
status = nf90_close(GLOBAL%NETCDF_OUTPUT_FILE%fp)

end subroutine FILE_CLOSE

!>subroutine to write regional variables to file
subroutine Write_Regional(i,REG,GLOBAL)

  type(REGIONAL_template),intent(in) :: REG
  type(GLOBAL_template),intent(in) :: GLOBAL
  integer,intent(in) :: i

  write(GLOBAL%REGIONAL_FILE%fp,*)i,REG

end subroutine Write_Regional

!>subroutine to write catchment variables to file
subroutine Write_Catchment(i,CAT,GLOBAL)

  type(CATCHMENT_template),dimension(:),intent(in) :: CAT
  type(GLOBAL_template),intent(in) :: GLOBAL
  integer,intent(in) :: i
  write(GLOBAL%CATCHMENT_FILE%fp,*)i,CAT%qb,CAT%qsurf*CAT%area,CAT%qsurf,CAT%ettot,CAT%pptsum,CAT%zbar,CAT%qb/CAT%area

end subroutine Write_Catchment

! ====================================================================
!
!			subroutine rdveg
!
! ====================================================================
!
!>Subroutine to read and initiailize simulation constant vegetation and land cover parameters.
!
! ====================================================================

      subroutine rdveg(GLOBAL,CAT,GRID,REG,IO)

      implicit none
      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
      type (REGIONAL_template),intent(inout) :: REG
      type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
      type (IO_template),intent(inout) :: IO
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%VEG_FILE%nlon,GLOBAL%VEG_FILE%nlat)
      type (GRID_VEG_template) :: GRID_DVEG_2D(GLOBAL%DVEG_FILE%nlon,GLOBAL%DVEG_FILE%nlat)
      character(len=200) :: filename
      integer :: vegnvars,dvegnvars,ipos,jpos,ilat,ilon,ip
      real,dimension(:,:,:),allocatable :: TempArray
      real*8 :: frcov(GLOBAL%nrow*GLOBAL%ncol,GLOBAL%ncatch+1),wc0
      integer :: jj,kk
      vegnvars = 20
      dvegnvars = 2
      allocate(TempArray(GLOBAL%VEG_FILE%nlon,GLOBAL%VEG_FILE%nlat,vegnvars))

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
  !Set all missing values to the areal average
  call Replace_Undefined(TempArray,GLOBAL%Veg_File%undef,GLOBAL%Veg_File%nlon,&
                   GLOBAL%Veg_File%nlat,vegnvars)


      GRID_VEG_2D%ivgtyp = int(TempArray(:,:,1))
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
      allocate(TempArray(GLOBAL%DVEG_FILE%nlon,GLOBAL%DVEG_FILE%nlat,dvegnvars))

      print*,"rdveg:  Reading in the dynamic vegetation properties"

      read(GLOBAL%DVEG_FILE%fp,rec=1)TempArray(:,:,:)
  !Set all missing values to the areal average
  call Replace_Undefined(TempArray,GLOBAL%DVEG_File%undef,GLOBAL%DVEG_File%nlon,&
                   GLOBAL%DVEG_File%nlat,dvegnvars)

      GRID_DVEG_2D%xlai = dble(TempArray(:,:,1))
      GRID_DVEG_2D%albd = dble(TempArray(:,:,2))

  GRID%VEG%ivgtyp = 0
  GRID%VEG%ilandc = 0
  do ilat = 1,GLOBAL%nrow
    do ilon = 1,GLOBAL%ncol
      ip = IO%ipixnum(ilat,ilon)
      if (ip .gt. 0)then
        !Static Vegetation
        GRID(ip)%VEG%ivgtyp = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%ivgtyp
        GRID(ip)%VEG%emiss = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%emiss
        GRID(ip)%VEG%za = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%za
        GRID(ip)%VEG%zww = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%zww
        GRID(ip)%VEG%z0m = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%z0m
        GRID(ip)%VEG%z0h = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%z0h
        GRID(ip)%VEG%zpd = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%zpd
        GRID(ip)%VEG%rsmin = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%rsmin
        GRID(ip)%VEG%rsmax = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%rsmax
        GRID(ip)%VEG%Rpl = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%Rpl
        GRID(ip)%VEG%f3vpdpar = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%f3vpdpar
        GRID(ip)%VEG%f4temppar = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%f4temppar
        GRID(ip)%VEG%trefk = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%trefk
        GRID(ip)%VEG%tcbeta = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%tcbeta
        GRID(ip)%VEG%extinct = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%extinct
        GRID(ip)%VEG%canclos = GRID_VEG_2D(GRID(ip)%VEG%static_MAP%ilon,GRID(ip)%VEG%static_MAP%ilat)%canclos
        !Dynamic Vegetation
        GRID(ip)%VEG%xlai = GRID_DVEG_2D(GRID(ip)%VEG%dynamic_MAP%ilon,GRID(ip)%VEG%dynamic_MAP%ilat)%xlai
        GRID(ip)%VEG%albd = GRID_DVEG_2D(GRID(ip)%VEG%dynamic_MAP%ilon,GRID(ip)%VEG%dynamic_MAP%ilat)%albd
      endif
    enddo
  enddo

! ####################################################################
! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

  do kk=1,maxval(IO%ipixnum)
    GRID(kk)%VEG%ilandc = kk
  enddo

! ####################################################################
!! Convert the 2-d arrays to the model's 1-d arrays
! ####################################################################

        do kk=1,GLOBAL%nlandc

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
  GRID%VEG%tcbeta = exp(-0.5*GRID%VEG%xlai)
  GRID%VEG%xlai_wsc = GRID%VEG%xlai
  GRID%VEG%albw = GRID%VEG%albd !Move to its own file
!
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

         do 560 kk=1,GLOBAL%nlandc

            if (GRID(kk)%VEG%ivgtyp.eq.0 .and. GRID(kk)%VEG%ilandc .gt. 0) then
               
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
!>Subroutine to initialize simulation total water balance variables.
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

            GRID(kk)%VARS%sm1 = GLOBAL%smpet0
            GRID(kk)%VARS%sm1_u = GLOBAL%smpet0
            GRID(kk)%VARS%sm1_f = zero
            GRID(kk)%VARS%smdthetaidt = zero

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
!> Subroutine to write raw binary data to the output file
! ====================================================================

subroutine WRITE_NETCDF(datain,ipixnum,i,GLOBAL,varid)

  implicit none
  type(GLOBAL_template),intent(in) :: GLOBAL
  real*8,intent(in) :: datain(GLOBAL%nrow*GLOBAL%ncol)
  integer,intent(in) :: ipixnum(GLOBAL%nrow,GLOBAL%ncol)
  integer,intent(in) :: varid
  real :: dataout(GLOBAL%ncol,GLOBAL%nrow)
  integer :: i
  integer :: start(3),count(3),status

  call MODEL2GRID(datain,GLOBAL%nrow,GLOBAL%ncol,ipixnum,GLOBAL%NETCDF_OUTPUT_FILE%undef,dataout)

  !Write the netcdf data
  count = [GLOBAL%ncol,GLOBAL%nrow,1]
  start = [1,1,i]
  status = nf90_put_var(GLOBAL%NETCDF_OUTPUT_FILE%fp,varid,dataout,start,count)

end subroutine WRITE_NETCDF

! ====================================================================
!> Subroutine to convert from model array to grid
! ====================================================================

subroutine MODEL2GRID(datain,nrow,ncol,ipixnum,undef,dataout)

  implicit none
  integer,intent(in) :: nrow,ncol
  real*8,intent(in) :: datain(nrow*ncol)
  real,intent(inout) :: dataout(ncol,nrow)
  real*8,intent(in) :: undef
  integer,intent(in) :: ipixnum(nrow,ncol)
  integer :: irow,icol,x,y

  ! ====================================================================
  ! Loop through the image and write each value in proper location.
  ! ====================================================================

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

          dataout(x,y) = datain(ipixnum(irow,icol))

        else

          dataout(x,y) = real(undef)

        endif

      enddo

    enddo

  end subroutine MODEL2GRID

! ====================================================================
!
!			subroutine rdsoil
!
! ====================================================================
!
!>Subroutine to read and initialize time in-variant soil parameters
!
! ====================================================================

      subroutine rdsoil(GLOBAL,CAT,GRID,IO)

      implicit none

      type (GLOBAL_template),intent(inout) :: GLOBAL
      type (GRID_template),dimension(:),allocatable,intent(inout) :: GRID
      type (CATCHMENT_template),dimension(:),allocatable,intent(inout) :: CAT
      type (IO_template),intent(inout) :: IO
      type (GRID_SOIL_template) :: GRID_SOIL_2D(GLOBAL%SOIL_FILE%nlon,GLOBAL%SOIL_FILE%nlat)
      type (GRID_VEG_template) :: GRID_VEG_2D(GLOBAL%SOIL_FILE%nlon,GLOBAL%SOIL_FILE%nlat)
      integer :: soilnvars,ipos,jpos,jj,kk,nn,ip,ilat,ilon
      integer :: icount(GLOBAL%ncol*GLOBAL%nrow,GLOBAL%ncatch+1)
      real,dimension(:,:,:),allocatable :: TempArray
      real*8 :: psic(GLOBAL%ncol*GLOBAL%nrow),tempsum,dtaken
      real*8 :: frsoil(GLOBAL%ncol*GLOBAL%nrow,GLOBAL%ncatch+1)
      soilnvars = 23

 allocate(TempArray(GLOBAL%SOIL_FILE%nlon,GLOBAL%SOIL_FILE%nlat,soilnvars))

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
  !Set all missing values to the areal average
  call Replace_Undefined(TempArray,GLOBAL%Soil_File%undef,GLOBAL%Soil_File%nlon,&
                   GLOBAL%Soil_File%nlat,soilnvars)
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
! Read the soil properties
! ====================================================================

  do kk=1,maxval(IO%ipixnum)
    GRID(kk)%SOIL%isoil = kk
  enddo

  do ilat = 1,GLOBAL%nrow
    do ilon = 1,GLOBAL%ncol
      ip = IO%ipixnum(ilat,ilon)
      if (ip .gt. 0)then
        !Static Soil
        GRID(ip)%SOIL%thetas =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%thetas
        GRID(ip)%SOIL%thetar =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%thetar
        GRID(ip)%SOIL%bcbeta =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%bcbeta
        GRID(ip)%SOIL%psic =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%psic
        GRID(ip)%SOIL%xk0 =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%xk0
        GRID(ip)%SOIL%zdeep =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%zdeep
        GRID(ip)%SOIL%tdeep =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%tdeep
        GRID(ip)%SOIL%zmid =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%zmid
        GRID(ip)%SOIL%tmid0 =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%tmid0
        GRID(ip)%SOIL%rocpsoil =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%rocpsoil
        GRID(ip)%SOIL%quartz =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%quartz
        GRID(ip)%SOIL%ifcoarse =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%ifcoarse
        GRID(ip)%SOIL%srespar1 = GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%srespar1
        GRID(ip)%SOIL%srespar2 = GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%srespar2
        GRID(ip)%SOIL%srespar3 = GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%srespar3
        GRID(ip)%SOIL%a_ice =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%a_ice
        GRID(ip)%SOIL%b_ice =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%b_ice
        GRID(ip)%SOIL%bulk_dens =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%bulk_dens
        GRID(ip)%SOIL%amp =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%amp
        GRID(ip)%SOIL%phase =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%phase
        GRID(ip)%SOIL%shift =GRID_SOIL_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%shift
        !Static Vegetation
        GRID(ip)%VEG%tc=GRID_VEG_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%tc
        GRID(ip)%VEG%tw=GRID_VEG_2D(GRID(ip)%SOIL%MAP%ilon,GRID(ip)%SOIL%MAP%ilat)%tw
      endif
    enddo
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

            CAT(jj)%psicav = CAT(jj)%psicav + frsoil(kk,jj)*GRID(kk)%SOIL%psic

560      continue

570   continue

      print*,"rdsoil:  Calculated average psi! for each catchment"

      return

      end subroutine rdsoil

! ====================================================================
!
!			subroutine rdatb
!
! ====================================================================
!
!> Subroutine to read in the topographi! index image and set up
!> at translation between the pixel number and the image row.and.&
!> column.
!
! ====================================================================

subroutine rdatb(atb,GLOBAL,IO)

  implicit none
  type (GLOBAL_template) :: GLOBAL
  type (IO_template),intent(inout) :: IO
  integer :: ip,irow,icol,fp,x,y
  real*8 :: atb(GLOBAL%nrow*GLOBAL%ncol),atb_2d(GLOBAL%ncol,GLOBAL%nrow)
  real*4 :: temp(GLOBAL%ncol,GLOBAL%nrow)
  real*8 :: undef

  !Read in the topographic index data
  read(GLOBAL%TI_FILE%fp,rec=1) temp

  ! Convert from single to double point precision
  atb_2d = temp

  ! Extract all the positioning information of the original model
  x = 1
  y = 0
  IO%ipixnum = 0
  ip = 0
  do irow = 1,GLOBAL%nrow
    do icol =1,GLOBAL%ncol
      if (y.eq.GLOBAL%nrow)then
        y=0
        x=x+1
      endif
      y = y +1
      if (nint(atb_2d(x,y)).ne.nint(GLOBAL%TI_FILE%undef))then
        ip = ip + 1
        IO%ipixnum(irow,icol) = ip
      endif
    enddo
  enddo

  ! Convert to model format
  call convert_grads2model(atb,atb_2d,IO%ipixnum,GLOBAL%nrow,&
                           GLOBAL%ncol,GLOBAL%TI_FILE%undef)

  ! Set the total number of pixels.
  GLOBAL%npix = ip

end subroutine rdatb

! ====================================================================
!
!                       function factln
!
! ====================================================================
!
!> Returns the natural log of the i factorial.
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

!>Subroutine to read the general parameters file
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
  !Number of time steps between updating vegetation parameters 
  call Extract_Info_General_File('dtveg',GLOBAL,GLOBAL%dtveg)
  !Topographic index filename
  call Extract_Info_General_File_File_Info('TI_fname',GLOBAL,GLOBAL%TI_FILE)
  !Basin distribution filename
  call Extract_Info_General_File_File_Info('Subbasin_fname',GLOBAL,GLOBAL%Subbasin_FILE)
  !Saturated Hydraulic Conductivity filename
  call Extract_Info_General_File_File_Info('K_0_fname',GLOBAL,GLOBAL%K0_FILE)
  !Catchment parameter table
  call Extract_Info_General_File('CL_Table_fname',GLOBAL,GLOBAL%CL_table_FILE%fname)
  !Soil file
  call Extract_Info_General_File_File_Info('SOIL_fname',GLOBAL,GLOBAL%SOIL_FILE)
  !Vegetation file
  call Extract_Info_General_File_File_Info('VEG_fname',GLOBAL,GLOBAL%VEG_FILE)
  !Dynamic vegetation file
  call Extract_Info_General_File_File_Info('DVEG_fname',GLOBAL,GLOBAL%DVEG_FILE)
  !Forcing file
  call Extract_Info_General_File_File_Info('FORCING_fname',GLOBAL,GLOBAL%FORCING_FILE)
  !Output file
  call Extract_Info_General_File('OUTPUT_fname',GLOBAL,GLOBAL%OUTPUT_FILE%fname)
  !ET Output file
  call Extract_Info_General_File('ET_fname',GLOBAL,GLOBAL%ET_FILE%fname)
  !QSURF Output file
  call Extract_Info_General_File('QS_fname',GLOBAL,GLOBAL%QS_FILE%fname)
  !Output file
  call Extract_Info_General_File('REGIONAL_fname',GLOBAL,GLOBAL%REGIONAL_FILE%fname)
  !Output Catchment file
  call Extract_Info_General_File('CATCHMENT_fname',GLOBAL,GLOBAL%CATCHMENT_FILE%fname)
  !Output Catchment file
  call Extract_Info_General_File('GSTI_fname',GLOBAL,GLOBAL%GSTI_FILE%fname)
  !Number of threads used in openmp
  call Extract_Info_General_File('nthreads',GLOBAL,GLOBAL%nthreads)
  !Initial time stamp in epoch time
  call Extract_Info_General_File('itime',GLOBAL,GLOBAL%itime)
  !Flag defining how to read in the vegetation
  call Extract_Info_General_File('dynamic_vegetation',GLOBAL,GLOBAL%dynamic_vegetation)
  !Flag defining the saturate hydraulic conductivity profile (0-Sivapalan,1987,1-Chaney,2013)
  call Extract_Info_General_File('KS_PROFILE_TYPE',GLOBAL,GLOBAL%KS_TYPE)
  !Output file
  call Extract_Info_General_File_File_Info('NETCDF_OUTPUT_fname',GLOBAL,GLOBAL%NETCDF_OUTPUT_FILE)
  !Option for vertical Ks change with depth 
  GLOBAL%ikopt = GLOBAL%KS_TYPE
  
  !Close the file
  close(GLOBAL%GENERAL_FILE%fp)
  
end subroutine

!>Subroutine to extract all integer variables from the general parameters file
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

!>Subroutine to extract all string variables from the general parameters file
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

!>Subroutine to extract all double precision variables from the general parameters file
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

!>Subroutine to extract all string variables from the general parameters file
subroutine Extract_Info_General_File_File_Info(strid,GLOBAL,FILE_INFO)

  implicit none
  type (GLOBAL_template),intent(inout) :: GLOBAL
  type (FILE_template),intent(inout) :: FILE_INFO
  integer :: flag
  character(len=*) :: strid
  character(len=400) :: filename
  character(len=200) :: string
  !Read until we find what we need
  do
    read(GLOBAL%GENERAL_FILE%fp,*,iostat=flag)string
    if (string .eq. strid)then
      !go back one record
      backspace(GLOBAL%GENERAL_FILE%fp)
      !read the desired parameter/file
      read(GLOBAL%GENERAL_FILE%fp,*)string,FILE_INFO%fname,FILE_INFO%spatial_res,&
        FILE_INFO%undef,FILE_INFO%minlat,FILE_INFO%minlon,FILE_INFO%nlat,FILE_INFO%nlon
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

end subroutine Extract_Info_General_File_File_Info

!>Subroutine to convert the 1-d format back to 2-d grads format
subroutine convert_model2grads(array_1d,array_2d,ipixnum,nrow,ncol,undef)

  implicit none
  integer,intent(in) :: nrow,ncol,ipixnum(nrow,ncol)
  real*8,intent(in) :: array_1d(nrow*ncol),undef
  real*8,intent(inout) :: array_2d(ncol,nrow)
  integer :: x,y,irow,icol

  !Map the kk position ot hte i,j position
  x = 1
  y = 0
  do irow = 1,nrow
    do icol =1,ncol
      if (y.eq.nrow)then
        y=0
        x=x+1
      endif
      y = y +1
      if (ipixnum(irow,icol).gt.0)then
        array_2d(x,y) = array_1d(ipixnum(irow,icol))
      else
        array_2d(x,y) = undef
      endif
    enddo
  enddo

end subroutine

!>Subroutine to convert the 1-d format back to grads format
subroutine convert_grads2model_double(array_1d,array_2d,ipixnum,nrow,ncol,undef)

  implicit none
  integer,intent(in) :: nrow,ncol
  integer,intent(in) :: ipixnum(nrow,ncol)
  real*8,intent(in) :: array_2d(ncol,nrow),undef
  real*8,intent(out) :: array_1d(ncol*nrow)
  integer :: x,y,irow,icol

  !Map the kk position ot hte i,j position
  x = 1
  y = 0
  array_1d = 0
  do irow = 1,nrow
    do icol =1,ncol
      if (y.eq.nrow)then
        y=0
        x=x+1
      endif
      y = y +1
      if (ipixnum(irow,icol).ne.0)then
        array_1d(ipixnum(irow,icol)) = array_2d(x,y)
      endif
    enddo
  enddo

end subroutine

subroutine convert_grads2model_int(array_1d,array_2d,ipixnum,nrow,ncol,undef)

  implicit none
  integer,intent(in) :: nrow,ncol
  integer,intent(in) :: ipixnum(nrow,ncol)
  integer,intent(in) :: array_2d(ncol,nrow)
  real*8,intent(in) :: undef
  integer,intent(out) :: array_1d(ncol*nrow)
  integer :: x,y,irow,icol

  !Map the kk position ot hte i,j position
  x = 1
  y = 0
  do irow = 1,nrow
    do icol =1,ncol
      if (y.eq.nrow)then
        y=0
        x=x+1
      endif
      y = y +1
      if (ipixnum(irow,icol).ne.0)then
        array_1d(ipixnum(irow,icol)) = array_2d(x,y)
      endif
    enddo
  enddo

end subroutine

!>Subroutine to create the mapping between two grids
subroutine spatial_mapping(GLOBAL,MAP,FILE_INFO,ipixnum)

  implicit none
  type(GLOBAL_template),intent(in) :: GLOBAL
  type(FILE_template),intent(in) :: FILE_INFO
  type(MAP_template),intent(inout) :: MAP(GLOBAL%nrow*GLOBAL%ncol)
  integer :: ipixnum(GLOBAL%nrow,GLOBAL%ncol)
  integer :: irow,icol,irow_alt,icol_alt,x,y
  real*8 :: lat_ref,lon_ref,lat_res,lon_res

  !Find the latitude and longitude of the point of the reference grid
  x = 1
  y = 0
  irow = 0
  icol = 0
  MAP%ilat = 0
  MAP%ilon = 0
  do irow = 1,GLOBAL%nrow
    do icol = 1,GLOBAL%ncol
      if (y.eq.GLOBAL%nrow)then
        y=0
        x=x+1
      endif
      y = y +1
      if (ipixnum(irow,icol) .gt. 0)then
        lat_ref = GLOBAL%minlat + (y-1)*GLOBAL%spatial_res 
        lon_ref = GLOBAL%minlon + (x-1)*GLOBAL%spatial_res
        ! Find the closest latitude and longitude on the other grid
        irow_alt = nint((lat_ref - FILE_INFO%minlat)/FILE_INFO%spatial_res) + 1
        icol_alt = nint((lon_ref - FILE_INFO%minlon)/FILE_INFO%spatial_res) + 1
        ! Ensure it is within the bounds
        if (irow_alt .gt. FILE_INFO%nlat)irow_alt = FILE_INFO%nlat
        if (irow_alt .lt. 1)irow_alt = 1
        if (icol_alt .gt. FILE_INFO%nlon)icol_alt = FILE_INFO%nlon
        if (icol_alt .lt. 1)icol_alt = 1
        ! Place mapping in array
        MAP(ipixnum(irow,icol))%ilat = irow_alt
        MAP(ipixnum(irow,icol))%ilon = icol_alt
      endif
    enddo
  enddo

end subroutine

!>Subroutine to create output netcdf file
subroutine Create_Netcdf_Output(FILE_INFO,GLOBAL,nvars)

  implicit none
  type(FILE_template),intent(inout) :: FILE_INFO
  type(GLOBAL_template),intent(in) :: GLOBAL 
  integer,intent(in) :: nvars
  integer :: status,i
  integer :: LonDimId,LatDimId,TimeDimId,lon_varid,lat_varid,time_varid
  real*8 :: time(GLOBAL%ndata),lats(FILE_INFO%nlat),lons(FILE_INFO%nlon)

  !Create and open the new file
  status = nf90_create(FILE_INFO%fname,nf90_hdf5,FILE_INFO%fp)

  !Define the dimensions of the new file
  status = nf90_def_dim(FILE_INFO%fp,'lon',FILE_INFO%nlon, LonDimId)
  status = nf90_def_dim(FILE_INFO%fp,'lat',FILE_INFO%nlat, LatDimId)
  status = nf90_def_dim(FILE_INFO%fp,'t',GLOBAL%ndata, TimeDimId)

  !Define the coordinate variables
  status = nf90_def_var(FILE_INFO%fp,'lon',NF90_DOUBLE,LonDimId,lon_varid)
  status = nf90_def_var(FILE_INFO%fp,'lat',NF90_DOUBLE,LatDimId,lat_varid)
  status = nf90_def_var(FILE_INFO%fp,'t',NF90_DOUBLE,TimeDimId,time_varid)

  !Assign units attributes to coordinate var data.
  status = nf90_put_att(FILE_INFO%fp,lat_varid,'units','degrees_north')
  status = nf90_put_att(FILE_INFO%fp,lon_varid,'units','degrees_east')
  status = nf90_put_att(FILE_INFO%fp,time_varid,'units','hours')
  status = nf90_put_att(FILE_INFO%fp,lat_varid,'long_name','Latitude')
  status = nf90_put_att(FILE_INFO%fp,lon_varid,'long_name','Longitude')
  status = nf90_put_att(FILE_INFO%fp,time_varid,'units','hours since 1970-01-01 00:00:00.0')
  status = nf90_put_att(FILE_INFO%fp,time_varid,'long_name','Time')

  !Set the data attributes
  do i=1,nvars
   FILE_INFO%varid(i) = 0
   status = nf90_def_var(FILE_INFO%fp,FILE_INFO%var_name(i),NF90_REAL,&
            [LonDimId,LatDimId,TimeDimId],FILE_INFO%varid(i))
   status = nf90_put_att(FILE_INFO%fp,FILE_INFO%varid(i),'long_name',FILE_INFO%var_name(i))
   status = nf90_put_att(FILE_INFO%fp,FILE_INFO%varid(i),'_FillValue',real(FILE_INFO%undef))
  enddo

  !End define mode
  status = nf90_enddef(FILE_INFO%fp)

  !Create the coordainte variable data
  do i=1,FILE_INFO%nlat
   lats(i) = FILE_INFO%minlat + (i-1)*FILE_INFO%spatial_res
  enddo
  do i=1,FILE_INFO%nlon
   lons(i) = FILE_INFO%minlon + (i-1)*FILE_INFO%spatial_res
  enddo
  do i=1,GLOBAL%ndata
   time(i) = GLOBAL%itime + (i-1)*GLOBAL%dt
  enddo
  lats = lats/3600.0d0
  lons = lons/3600.0d0
  print*,lats
  time = time/3600.0d0

  !Write the coordinate variable data
  status = nf90_put_var(FILE_INFO%fp,lat_varid,lats/3600.0)
  status = nf90_put_var(FILE_INFO%fp,lon_varid,lons/3600.0)
  status = nf90_put_var(FILE_INFO%fp,time_varid,time)
!  print*,FILE_INFO,lat_varid,lon_varid,time_varid
!  stop

end subroutine Create_Netcdf_Output


END MODULE MODULE_IO
