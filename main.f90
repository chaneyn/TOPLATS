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
      integer :: nthreads,chunksize
      integer :: ntdveg
      ntdveg = 1 !remember the dynamic vegetation parameter time step
                !rdveg_update.f90
      chunksize = 1
      nthreads = 1

! ####################################################################
! Initialize unit testing
! ####################################################################

call init_fruit

! ####################################################################
! Open all files
! ####################################################################

      call FILE_OPEN()

! ####################################################################
! Call rddata to open files, read in time in-variant parameters,&
! and initialize simulation sums.
! ####################################################################

      call rddata(GLOBAL,GRID,REG,CAT,IO)

! ####################################################################
! Loop through the simulation time.
! ####################################################################

      do i=1,GLOBAL%ndata

          print*, "Time Step: ",i," Year: ",iyear," Julian Day: ",&
                    iday," Hour: ",ihour

! ####################################################################
! Update the vegetation parameters if required.
! ####################################################################

         if (mod(i,GLOBAL%dtveg).eq.0) then

            call rdveg_update(GLOBAL,ntdveg,GRID)

         endif

! ####################################################################
! Initialize water balance variables for the time step.
! ####################################################################

        call instep(i,GLOBAL%ncatch,djday,GLOBAL%dt,&
                REG,CAT)

! ####################################################################
! Read meteorological data.
! ####################################################################

        call rdatmo(GLOBAL%nrow,GLOBAL%ncol,IO%ipixnum,iyear,&
                    iday,ihour,i,GRID%MET)

! ####################################################################
! Loop through each pixel in atanb map.
! ####################################################################

            lakpix=0
            m_lc = 1
            m_px = 1
            u_lc = 1
            u_px = 1
            l_lc = 1
            l_px = 1

GRID%VARS%row = 997.d0
GRID%VARS%cph2o = 4186.d0
GRID%VARS%cp = 1005.d0
GRID%VARS%roi = 850.d0
GRID%VARS%rzsmold = 0.d0

call OMP_SET_NUM_THREADS(8)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(sw_lc,sw_px,s_lc,s_px) 
!$OMP DO SCHEDULE(DYNAMIC) ORDERED

            do ipix=1,GLOBAL%npix

! ####################################################################
! Calculate the matrix indices depending on the land cover geometry.
! ####################################################################

               call clc_ind (GRID(ipix)%VEG%ilandc,ipix,SNOW_RUN,sw_lc,sw_px,SNW_FLG,s_lc,s_px)

! ####################################################################
! Run the point model
! ####################################################################

               call land_lake(NEWSTORM(1,1),ipix,i,GLOBAL%dt,GLOBAL%inc_frozen,GLOBAL%i_2l,l_px,&

! Factor to multiply the regional parameters with

       1.d0,&

! General vegetation parameters

       GRID(GRID(ipix)%VEG%ilandc)%VEG%ivgtyp,&

! Snow pack variables

       GRID(sw_px)%VARS%PackWater,GRID(sw_px)%VARS%SurfWater,Swq(sw_px),GRID(sw_px)%VARS%VaporMassFlux,&
       TPack(sw_px),TSurf(sw_px),GRID(sw_px)%VARS%r_MeltEnergy,GRID(sw_px)%VARS%Outflow,&
       xleact_snow(sw_px),hact_snow(sw_px),rn_snow(sw_px),GRID(s_px)%VARS%PackWater_us,&
       GRID(s_px)%VARS%SurfWater_us,Swq_us(s_px),GRID(s_px)%VARS%VaporMassFlux_us,TPack_us(s_px),&
       TSurf_us(s_px),GRID(s_px)%VARS%r_MeltEnergy_us,&
       GRID(s_px)%VARS%Outflow_us,xleact_snow_us(s_px),hact_snow_us(s_px),rn_snow_us(s_px),dens(sw_px),dens_us(s_px),&
       dsty(sw_px),dsty_us(s_px),Sdepth(sw_px),Sdepth_us(s_px),&

! Albedos of the over story, under story,&
! and moss layer

       alb_snow(ipix),&

! Meteorological data

       GRID(ipix)%MET,&

! Temperature variables

       Tdeepstep(GRID(ipix)%SOIL%isoil),&

! Soil parameters
        
       GRID(GRID(ipix)%SOIL%isoil)%SOIL,&
       GRID(GRID(ipix)%SOIL%isoil)%SOIL%ifcoarse,&
       GLOBAL%zrzmax,&

! Vegetation parameters

       GRID(GRID(ipix)%VEG%ilandc)%VEG,&

! Constants
       GLOBAL%toleb,GLOBAL%maxnri,&

! Energy balance variables

       rib(ipix),&

! Water balance variables
       
       GRID(ipix)%VARS,&
       GRID(ipix)%VARS%cuminf,&
       GRID(ipix)%VARS%sorp,GRID(ipix)%VARS%cc,&
       GRID(ipix)%VARS%sesq,GRID(GRID(ipix)%SOIL%isoil)%SOIL%corr,&
       GRID(GRID(ipix)%SOIL%isoil)%SOIL%idifind,&
       GRID(ipix)%VARS%wcip1,GRID(GRID(ipix)%SOIL%isoil)%SOIL%par,&
       GLOBAL%smpet0,&

! Storm parameters

       istmst(ipix),intstm(ipix),&
       intstp(ipix),GLOBAL%endstm,istorm(ipix),&
       GRID(ipix)%VARS%xintst,&

! Topmodel parameters

       CAT(GRID(ipix)%VARS%icatch)%ff,GRID(ipix)%VARS%atanb,CAT(GRID(ipix)%VARS%icatch)%xlamda,&

! Regional saturation parameters

       REG,&

! Different option parameters

       GLOBAL%iopthermc,GLOBAL%iopgveg,GLOBAL%iopthermc_v,GLOBAL%iopsmini,GLOBAL%ikopt,&
       GLOBAL%irestype,GLOBAL%ioppet,GLOBAL%iopveg,GLOBAL%iopstab,GLOBAL%iopwv,&

! Catchment data
       CAT(GRID(ipix)%VARS%icatch))

!TEMPORARY PASS BACK TO ORIGINAL VARIABLES
       !GRID
       rzsm(ipix) = GRID(ipix)%VARS%rzsm
       tzsm(ipix) = GRID(ipix)%VARS%tzsm
       rzsm_f(ipix) = GRID(ipix)%VARS%rzsm_f
       tzsm_f(ipix) = GRID(ipix)%VARS%tzsm_f
       pnet(ipix) = GRID(ipix)%VARS%pnet
       xinact(ipix) = GRID(ipix)%VARS%xinact
       runtot(ipix) = GRID(ipix)%VARS%runtot
       irntyp(ipix) = GRID(ipix)%VARS%irntyp
       !Meteorological Variables
       Tincan(ipix) = GRID(ipix)%VARS%Tincan
       rh_ic(ipix) = GRID(ipix)%VARS%rh_ic
       precip_o(ipix) = GRID(ipix)%VARS%precip_o
       precip_u(ipix) = GRID(ipix)%VARS%precip_u
       !Energy Fluxes
       rnetpn(ipix) = GRID(ipix)%VARS%rnetpn
       gbspen(ipix) = GRID(ipix)%VARS%gbspen
       evtact(ipix) = GRID(ipix)%VARS%evtact
       ievcon(ipix) = GRID(ipix)%VARS%ievcon
       ebspot(ipix) = GRID(ipix)%VARS%ebspot

! ....................................................................
! Sum the local water and energy balance fluxes.
! ....................................................................

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
       Swq(sw_px),Swq_us(s_px),&
       Sdepth(sw_px),Sdepth_us(s_px),&

! GRID Variables

       Tdeepstep(GRID(ipix)%SOIL%isoil))

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
