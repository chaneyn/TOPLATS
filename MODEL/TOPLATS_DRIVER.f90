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
real*8 :: start_time,end_time,omp_get_wtime

!####################################################################
! Read the general filename
!####################################################################

call Initialize_Model(GLOBAL,GRID,REG,CAT,IO)

GRID%VARS%rzsm1 = GRID%VARS%sm1(1)
GRID%VARS%tzsm1 = GRID%VARS%sm1(2)
GRID%VARS%rzsm1_u = GRID%VARS%sm1_u(1)
GRID%VARS%tzsm1_u = GRID%VARS%sm1_u(2)
GRID%VARS%rzsm1_f = GRID%VARS%sm1_f(1)
GRID%VARS%tzsm1_f = GRID%VARS%sm1_f(2)
GRID%VARS%rzdthetaidt = GRID%VARS%smdthetaidt(1)
GRID%VARS%tzdthetaidt = GRID%VARS%smdthetaidt(2)
GLOBAL%zmax_layer(1) = GLOBAL%zrzmax
GLOBAL%nlayer = 2


!####################################################################
! Loop through the simulation time.
!####################################################################

do i=1,GLOBAL%ndata

  print*, "Time Step: ",i

!#####################################################################
! Read in the input data for this time step
!#####################################################################

  start_time = omp_get_wtime()
  call Read_Data(GLOBAL,GRID,CAT,IO,i)
  end_time = omp_get_wtime()
  print*,end_time-start_time

!#####################################################################
! Update each grid cell
!#####################################################################
 
  print*,'Updating the cells'

  call Update_Cells(GRID,CAT,GLOBAL,i)

!#####################################################################
! Update the catchments
!#####################################################################

  !GRID%VARS%z_layer(1) = GRID%VARS%zrz
  !GRID%VARS%z_layer(2) = GRID%VARS%ztz
  !GRID%VARS%sm(1) = GRID%VARS%rzsm
  !GRID%VARS%sm(2) = GRID%VARS%tzsm
  !GRID%VARS%sm1(1) = GRID%VARS%rzsm1
  !GRID%VARS%sm1(2) = GRID%VARS%tzsm1

  print*,'Updating the catchments'
  call Update_Catchments(GLOBAL,CAT,GRID)

  !GRID%VARS%rzsm1_u=GRID%VARS%sm1_u(1)
  !GRID%VARS%tzsm1_u=GRID%VARS%sm1_u(2)



!#####################################################################
! Update regional variables
!#####################################################################

  !print*,'Updating the regional variables'
  !call Update_Regional(REG,GRID,GLOBAL,CAT)

!#####################################################################
! Write out the data for this time step
!#####################################################################
  
  GRID%VARS%sm(1) = GRID%VARS%rzsm
  print*,'Writing the data'
  call Write_Data(GLOBAL,GRID,IO,REG,i,CAT)

enddo

!#####################################################################
! Output the GSTI
!#####################################################################

  !call Write_Binary(GRID%VARS%GSTI,1.0,GLOBAL%nrow,GLOBAL%ncol,&
  !                  IO%ipixnum,1,GLOBAL,GLOBAL%GSTI_file%fp)

!#####################################################################
! Finalize model and close files
!#####################################################################

call Finalize_Model(GLOBAL)

END PROGRAM
