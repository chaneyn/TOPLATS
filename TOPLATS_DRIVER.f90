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
USE MODULE_TOPMODEL,ONLY: instep,catflx,upzbar,sumflx,Update_Catchment

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
! Initialize water balance variables for the time step.
!#####################################################################

  call instep(i,GLOBAL%ncatch,GLOBAL%djday,GLOBAL%dt,REG,CAT)

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
       REG,CAT(icatch),GLOBAL)

!$OMP ORDERED
!$OMP CRITICAL


!$OMP END CRITICAL
!$OMP END ORDERED


!#####################################################################
! Sum the local water and energy balance fluxes.
!#####################################################################

    call sumflx(REG,CAT(icatch),GRID(ipix)%VARS,GLOBAL,& 
       GRID(ilandc)%VEG,GRID(isoil)%SOIL,GRID(ipix)%MET,&
       ilandc)

  enddo

!$OMP END DO
!$OMP END PARALLEL

!#####################################################################
! Loop through each catchment and update with the catchment module
!#####################################################################

  do ic=1,GLOBAL%ncatch

    call Update_Catchment(GLOBAL,CAT(ic),GRID,REG,ic)

  enddo

!#####################################################################
! Run Tests to compare to previous model
!#####################################################################

  call lswb(i,REG,GLOBAL,GRID)

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
