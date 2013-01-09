PROGRAM TOPLATS_DRIVER

! ====================================================================
!
!			TOPLATS Version 3.0
!
! ====================================================================
!
!               October, 1992
!               Last revised:  January, 2013
!
! TOPographically-based Land-Atmosphere Transfer Scheme
!
! ====================================================================

!Module containing the unit tests
USE FRUIT

!Module containing all the variables used in the model
USE MODULE_VARIABLES

!Module containing all the tests
USE MODULE_TESTS,ONLY: lswb

!Module containing all the I/O for the interface
USE MODULE_IO,ONLY: IO_template,FILE_OPEN,rddata,rdveg_update,rdatmo,&
                    file_close

!Module containing topmodel
USE MODULE_TOPMODEL,ONLY: instep_catchment,catflx,upzbar,Update_Catchments

!Module containing regional subroutines
USE MODULE_REGIONAL,ONLY: Update_Regional,instep_regional

!Module containing the cell model
USE MODULE_CELL,ONLY: Update_Cell

implicit none
type (GLOBAL_template) :: GLOBAL
type (GRID_template),dimension(:),allocatable :: GRID
type (REGIONAL_template) :: REG
type (CATCHMENT_template),dimension(:),allocatable :: CAT
type (IO_template) :: IO
integer :: i,ic,ipix,isoil,ilandc,icatch
GLOBAL%nthreads = 8

!####################################################################
! Initialize unit testing
!####################################################################

call init_fruit

!####################################################################
! Open all files
!####################################################################

call FILE_OPEN()

!####################################################################
! Call rddata to open files, read in time in-variant parameters,&
! and initialize simulation sums.
!####################################################################

call rddata(GLOBAL,GRID,REG,CAT,IO)

!####################################################################
! Loop through the simulation time.
!####################################################################

do i=1,GLOBAL%ndata

  print*, "Time Step: ",i," Year: ",GLOBAL%iyear," Julian Day: ",&
             GLOBAL%iday," Hour: ",GLOBAL%ihour

!#####################################################################
! Update the vegetation parameters if required.
!#####################################################################

  if (mod(i,GLOBAL%dtveg).eq.0) call rdveg_update(GLOBAL,GRID)

!#####################################################################
! Update the decimal Julian day.
!#####################################################################

  GLOBAL%djday = GLOBAL%djday + 0.0416666667d0*2.777777777d-4*GLOBAL%dt

!#####################################################################
! Initialize water balance variables for the time step.
!#####################################################################

  call instep_catchment(GLOBAL%ncatch,CAT)
  call instep_regional(REG)

!#####################################################################
! Read meteorological data.
!#####################################################################

  call rdatmo(i,GRID%MET,GLOBAL,IO)

!#####################################################################
! Loop through each grid cell
!#####################################################################

  call OMP_SET_NUM_THREADS(GLOBAL%nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipix,isoil,ilandc,icatch) 
!$OMP DO SCHEDULE(DYNAMIC) ORDERED

  do ipix=1,GLOBAL%npix

    isoil = GRID(ipix)%SOIL%isoil
    ilandc = GRID(ipix)%VEG%ilandc
    icatch = GRID(ipix)%VARS%icatch

!#####################################################################
! Update the current grid cell
!#####################################################################

    call Update_Cell(ipix,i,GRID(ipix)%MET,GRID(isoil)%SOIL,&
       GRID(ilandc)%VEG,GRID(ipix)%VARS,GRID(ipix)%VARS%wcip1,&
       CAT(icatch),GLOBAL)

  enddo

!$OMP END DO
!$OMP END PARALLEL

!#####################################################################
! Update the catchments
!#####################################################################

  call Update_Catchments(GLOBAL,CAT,GRID)

!#####################################################################
! Update regional variables
!#####################################################################

  call Update_Regional(REG,GRID,GLOBAL)

!#####################################################################
! Run Tests to compare to previous model
!#####################################################################

  call lswb(i,REG,GLOBAL,GRID,CAT)

enddo

! ####################################################################
! Close all files
! ####################################################################

call FILE_CLOSE()

! ####################################################################
! Finalize unit testing and print summary
! ####################################################################

call fruit_summary
call fruit_finalize

END PROGRAM
