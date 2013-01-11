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

!Module containing all the variables used in the model
USE MODULE_VARIABLES

!Module containing all the I/O for the interface
USE MODULE_IO,ONLY: IO_template,FILE_OPEN,rddata,rdveg_update,rdatmo,&
                    file_close,Write_Regional

!Module containing topmodel
USE MODULE_TOPMODEL,ONLY: instep_catchment,Update_Catchments

!Module containing regional subroutines
USE MODULE_REGIONAL,ONLY: Update_Regional

!Module containing the cell model
USE MODULE_CELL,ONLY: Update_Cells

implicit none
type (GLOBAL_template) :: GLOBAL
type (GRID_template),dimension(:),allocatable :: GRID
type (REGIONAL_template) :: REG
type (CATCHMENT_template),dimension(:),allocatable :: CAT
type (IO_template) :: IO
integer :: i,ic,ipix,isoil,ilandc,icatch
GLOBAL%nthreads = 8

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

!#####################################################################
! Read meteorological data.
!#####################################################################

  call rdatmo(i,GRID%MET,GLOBAL,IO)

!#####################################################################
! Update each grid cell
!#####################################################################

  call Update_Cells(GRID,CAT,GLOBAL,i)

!#####################################################################
! Update the catchments
!#####################################################################

  call Update_Catchments(GLOBAL,CAT,GRID)

!#####################################################################
! Update regional variables
!#####################################################################

  call Update_Regional(REG,GRID,GLOBAL,CAT)

!#####################################################################
! Output regional variables
!#####################################################################

  call Write_Regional(i,REG)

enddo

!#####################################################################
! Close all files
!#####################################################################

call FILE_CLOSE()

END PROGRAM
