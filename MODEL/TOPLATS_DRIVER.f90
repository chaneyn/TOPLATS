!>TOPLATS interface
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
USE MODULE_IO,ONLY: IO_template,Read_Data,Finalize_Model,Write_Data,Initialize_Model

!Module containing topmodel
USE MODULE_TOPMODEL,ONLY: Update_Catchments

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

!####################################################################
! Read the general filename
!####################################################################

call Initialize_Model(GLOBAL,GRID,REG,CAT,IO)

!####################################################################
! Loop through the simulation time.
!####################################################################

do i=1,GLOBAL%ndata

  print*, "Time Step: ",i

!#####################################################################
! Read in the input data for this time step
!#####################################################################

  call Read_Data(GLOBAL,GRID,CAT,IO,i)

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
! Write out the data for this time step
!#####################################################################

  call Write_Data(GLOBAL,GRID,IO,REG,i)

enddo

!#####################################################################
! Finalize model and close files
!#####################################################################

call Finalize_Model(GLOBAL)

END PROGRAM
