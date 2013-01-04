! ==================================================================
!
!			TOPLATS Version 4.0
!
! ====================================================================
!
!               October, 1992
!               Last revised:  July, 2012
!
! TOPographically-based Land-Atmosphere Transfer Scheme
!    for regional and global atmosREG%peric models and
!    studies of macroscale water and energy balance --
!    Distributed Version.
!
! Model developed at Princeton University under direction of:
!
!  Eric F. Wood
!  Department of Civil and Environmental Engineering and Water Resources
!  Program of Environmental Engineering and Water Resources
!  Princeton University
!  Princeton, NJ 08544
!  Tel. (609) 258-4675 
!  Fax. (609) 258-1270
!  Email : efwood@.princeton.edu
!
! TOPLATS 4.0 takes the basis of TOPLATS 3.1. The code has been
! rewritten in Fortran 2003. For any questions, please contact:
!
!  Nathaniel W. Chaney
!  Department of Civil and Environmental Engineering and Water Resources
!  Program of Environmental Enginerring and Water Resources
!  Princeton University
!  Princeton, NJ 08544
!  Email : nchaney@princeton.edu
!
! ====================================================================

!Module containing the unit tests
USE FRUIT

!Module containing all the variables used in the model
USE MODULE_VARIABLES

USE MODULE_VARIABLES_OLD

!Module containing all the tests
USE MODULE_TESTS

!Module containing all the I/O for the interface
USE MODULE_IO

!Module containing topmodel
USE MODULE_TOPMODEL

!Module containing the cell model
USE MODULE_CELL

implicit none
type (GLOBAL_template) :: GLOBAL
type (GRID_template),dimension(:),allocatable :: GRID
type (REGIONAL_template) :: REG
type (CATCHMENT_template),dimension(:),allocatable :: CAT
type (IO_template) :: IO
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

  call instep(i,GLOBAL%ncatch,djday,GLOBAL%dt,REG,CAT)

!#####################################################################
! Read meteorological data.
!#####################################################################

  call rdatmo(i,GRID%MET,GLOBAL,IO)

!#####################################################################
! Loop through each pixel in atanb map.
!#####################################################################

  call OMP_SET_NUM_THREADS(GLOBAL%nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipix) 
!$OMP DO SCHEDULE(DYNAMIC) ORDERED

  do ipix=1,GLOBAL%npix

!#####################################################################
! Update the current grid cell
!#####################################################################

    call Update_Cell(ipix,i,GRID(ipix)%MET,GRID(GRID(ipix)%SOIL%isoil)%SOIL,&
       GRID(GRID(ipix)%VEG%ilandc)%VEG,GRID(ipix)%VARS,GRID(ipix)%VARS%wcip1,&
       REG,CAT(GRID(ipix)%VARS%icatch),GLOBAL)

!#####################################################################
! Sum the local water and energy balance fluxes.
!#####################################################################

    call sumflx(REG,CAT(GRID(ipix)%VARS%icatch),&
       GRID(ipix)%VARS,&        

! Factor to rescale all the local fluxes with

       1.d0,&

! General vegetation parameters

       GRID(GRID(ipix)%VEG%ilandc)%VEG%ivgtyp,&
       i,&
       canclos(GRID(ipix)%VEG%ilandc),GRID(ipix)%VEG%ilandc,GLOBAL%dt,&

! Grid data

       GRID(ipix)%MET%tdry,GRID(ipix)%MET%pptms,GRID(ipix)%VARS%wcip1,&

! Soil moisture variables

       GLOBAL%inc_frozen,&
       GRID(GRID(ipix)%SOIL%isoil)%SOIL%thetas,&
       GRID(ipix)%VARS%Swq,GRID(ipix)%VARS%Swq_us,&
       GRID(ipix)%VARS%Sdepth,GRID(ipix)%VARS%Sdepth_us,&

! GRID Variables

       GRID(GRID(ipix)%SOIL%isoil)%SOIL%Tdeepstep)

            enddo

!$OMP END DO
!$OMP END PARALLEL

!GRID variables
etpix = GRID%VARS%etpix

! --------------------------------------------------------------------
! Loop through each catchment to calculate catchment total fluxes
! (catflx) and update average water table depths (upzbar).
! --------------------------------------------------------------------

         do ic=1,GLOBAL%ncatch

            call catflx(i,ic,CAT(ic)%area,GLOBAL%pixsiz,&
                r_lakearea(ic),CAT(ic)%ettot,&
       CAT(ic)%etstsum,CAT(ic)%etwtsum,CAT(ic)%etlakesum,&
       CAT(ic)%etbssum,CAT(ic)%fbs,CAT(ic)%etdcsum,&
       CAT(ic)%etwcsum,CAT(ic)%pptsum,CAT(ic)%pnetsum,CAT(ic)%contot,&
       CAT(ic)%qsurf,CAT(ic)%sxrtot,CAT(ic)%xixtot,CAT(ic)%ranrun,&
       CAT(ic)%conrun,CAT(ic)%gwtsum,CAT(ic)%capsum,CAT(ic)%tzpsum,&
       CAT(ic)%rzpsum,CAT(ic)%fwcat)

               call upzbar(i,ic,GLOBAL%iopbf,CAT(ic)%q0,&
       CAT(ic)%ff,CAT(ic)%zbar,CAT(ic)%dtil,&
       CAT(ic)%basink,CAT(ic)%dd,CAT(ic)%xlength,CAT(ic)%gwtsum,CAT(ic)%capsum,CAT(ic)%area,&
       r_lakearea(ic),GLOBAL%dt,CAT(ic)%etwtsum,CAT(ic)%rzpsum,CAT(ic)%tzpsum,CAT(ic)%psicav,&
       GRID%VEG%ivgtyp,GRID%VEG%ilandc,GLOBAL%npix,GRID%VARS%icatch,zw,&
       GRID%SOIL%psic,GRID%SOIL%isoil,GLOBAL%zrzmax,GRID%VARS%tzsm1,GRID%SOIL%thetas,&
       GRID%VARS%rzsm1,CAT(ic)%zbar1,REG%qbreg,REG%zbar1rg,GLOBAL%pixsiz)

         enddo

! --------------------------------------------------------------------
! Call lswb to ouput areal average flux rates for the time step
! and sum simulation totals.  Then goto next time step.
! --------------------------------------------------------------------

         call lswb(i,r_lakearea,f_lake,veg_pdf,nlcs,veg,REG,GLOBAL,GRID)

      enddo

! ####################################################################
! Close all files
! ####################################################################

call FILE_CLOSE()


      write (*,*)
      write (*,*) 'Simulation terminated'
      write (*,*)

! ####################################################################
! Finalize unit testing and print summary
! ####################################################################

call fruit_summary !Summarize the fruit output for this time step
call fruit_finalize !Finalize the fruit l
      stop

      end
