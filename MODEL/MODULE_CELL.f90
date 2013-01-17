MODULE MODULE_CELL

USE MODULE_VARIABLES

USE MODULE_LAND

USE MODULE_ATMOS

USE MODULE_CANOPY

USE MODULE_SNOW

contains

!#####################################################################
!
!                        subroutine Update_Cells
!
!#####################################################################
!
! Solve the water and energy budget for all land surface area
!
!#####################################################################

  subroutine Update_Cells(GRID,CAT,GLOBAL,i)
  
    implicit none
    type (GRID_template),dimension(:),intent(inout) :: GRID
    type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
    type (GLOBAL_template),intent(in) :: GLOBAL
    integer,intent(in) :: i
    integer :: ipix,isoil,icatch,ilandc

!#####################################################################
! Update each grid cell
!#####################################################################

    call OMP_SET_NUM_THREADS(GLOBAL%nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipix,isoil,ilandc,icatch) 
!$OMP DO SCHEDULE(DYNAMIC) ORDERED

    do ipix=1,GLOBAL%npix

      isoil = GRID(ipix)%SOIL%isoil
      ilandc = GRID(ipix)%VEG%ilandc
      icatch = GRID(ipix)%VARS%icatch

      call Update_Cell(ipix,i,GRID(ipix)%MET,GRID(isoil)%SOIL,&
         GRID(ilandc)%VEG,GRID(ipix)%VARS,&
         CAT(icatch),GLOBAL)

    enddo

!$OMP END DO
!$OMP END PARALLEL

  end subroutine Update_Cells


! ====================================================================
!
!                        subroutine Update_Cell
!
! ====================================================================
!
! Solve the water and energy budget for a land surface area
!
! ====================================================================


      subroutine Update_Cell(ipix,i,GRID_MET,GRID_SOIL,GRID_VEG,&
               GRID_VARS,CAT,GLOBAL)

      implicit none
      include 'help/land_lake.h'
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_SOIL_template) :: GRID_SOIL
      type (GRID_MET_template) :: GRID_MET
      type (GRID_VARS_template) :: GRID_VARS
      type (CATCHMENT_template) :: CAT
      type (GLOBAL_template) :: GLOBAL
      type (CELL_VARS_template) :: CELL_VARS
      GLOBAL%mul_fac = 1.0d0

! TEMPORARY LOCATION TO PASS STRUCTURE INFORMATION TO OLD FORMAT

!Removal Causes Failure
sesq = GRID_VARS%sesq
xintst = GRID_VARS%xintst
smpet0 = GLOBAL%smpet0
mul_fac = GLOBAL%mul_fac
row = GRID_VARS%row!row
roi = GRID_VARS%roi!roi

!Point Data Initializations
!Water Balance
GRID_VARS%zrz = 0.d0
GRID_VARS%ztz = 0.d0
GRID_VARS%smold = 0.d0
GRID_VARS%rzsmold = 0.d0
GRID_VARS%tzsmold = 0.d0
GRID_VARS%capflx = 0.d0
GRID_VARS%difrz = 0.d0
GRID_VARS%diftz = 0.d0
GRID_VARS%grz = 0.d0
GRID_VARS%gtz = 0.d0
GRID_VARS%satxr = 0.d0
GRID_VARS%xinfxr = 0.d0
GRID_VARS%dc = 0.d0!d
GRID_VARS%fw = 0.d0!fw
GRID_VARS%dsrz = 0.d0!dsrz
GRID_VARS%rzrhs = 0.d0!rzrhs
GRID_VARS%dstz = 0.d0!dstz
GRID_VARS%tzrhs = 0.d0!tzrhs
GRID_VARS%dswc = 0.d0!dswc
GRID_VARS%wcrhs = 0.d0!wcrhs
!Energy Fluxes
GRID_VARS%epwms = 0.d0!epwms

!Catchment
zbar = CAT%zbar
ff = CAT%ff
xlamda = CAT%xlamda

! ====================================================================
! If the vegetation type is greater than or equal to zero then
! solve the water and energy balance for a land area.
! ====================================================================

      if (ivgtyp.ge.0) then

! ....................................................................
! Calculate the local energy fluxes and set
! up the storm/interstorm event times and flags.
! ..................................................................

       call atmos(ipix,i,GRID_VEG,&

! Meteorological data
       GRID_MET,tcel,vppa,psychr,xlhv,tkel,uzw,&
       appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&
       twet_ic,twet,qv,qv_ic,ra,ra_ic,&

! Temperature variables

       GRID_VARS,tkmid,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,&

! Energy fluxes and states

       epetd,epetd_us,dshact_moss,xle_act_moss,rnetd,&
       tkd,tkmidd,&
       tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss,&
       dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,epetw_us,&
       rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,&
       tkmidd_us,rnet_pot_moss,xle_p_moss,&
       h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss,&
       tsoilold,tkmidpet_us,tkmidpet_moss,&
       dspet_us,dspet_moss,&

       GRID_SOIL,&
       
! Moss parameters

       rib_moss,&
       epet_moss,&

! Vegetation parameters

       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,&
       f1,f2,f3,&
       f3vpdpar_us,f4temppar_us,&

! Constants

       roa,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,RaSnow,rib_us,&
       ravw,ravw_us,rahw,rahw_us,&

       GLOBAL)

! ....................................................................
! Calculate local wet canopy water balance.
! ....................................................................

 
         call canopy(ipix,GRID_VARS,GRID_VEG,GRID_MET,GLOBAL)

! ....................................................................
! Calculate the local land surface water/energy balance.
! ....................................................................

! ....................................................................
! Option 2 : the incoming long wave radiation for both under and over
! story is equal and is the atmospheri! incoming long wave radiation.
! The uncouples the radiation balances for both layers from each
! other.  This option is also used when under story is not represented.
! ....................................................................
       
       CELL_VARS%tcel = tcel 
       CELL_VARS%vppa = vppa
       CELL_VARS%psychr = psychr
       CELL_VARS%xlhv = xlhv
       CELL_VARS%tkel = tkel
       CELL_VARS%appa = appa
       CELL_VARS%vpsat = vpsat
       CELL_VARS%tcel_ic = tcel_ic
       CELL_VARS%vppa_ic = vppa_ic
       CELL_VARS%psychr_ic = psychr_ic
       CELL_VARS%xlhv_ic = xlhv_ic
       CELL_VARS%tkel_ic = tkel_ic
       CELL_VARS%vpsat_ic = vpsat_ic
       CELL_VARS%twet_ic = twet_ic
       CELL_VARS%twet = twet
       CELL_VARS%qv = qv
       CELL_VARS%qv_ic = qv_ic
       CELL_VARS%ra = ra
       CELL_VARS%ra_ic = ra_ic
       CELL_VARS%tkmid_us = tkmid_us
       CELL_VARS%tkact_us = tkact_us
       CELL_VARS%tskinact_moss = tskinact_moss
       CELL_VARS%tkact_moss = tkact_moss
       CELL_VARS%tkmid_moss = tkmid_moss
       CELL_VARS%epetd = epetd
       CELL_VARS%epetd_us = epetd_us
       CELL_VARS%dshact_moss = dshact_moss
       CELL_VARS%xle_act_moss = xle_act_moss
       CELL_VARS%tskinact_moss = tskinact_moss
       CELL_VARS%tkactd_moss = tkactd_moss
       CELL_VARS%tkmidactd_moss = tkmidactd_moss
       CELL_VARS%ds_p_moss = ds_p_moss
       CELL_VARS%dshact_us = dshact_us
       CELL_VARS%rnetw_us = rnetw_us
       CELL_VARS%xlew_us = xlew_us
       CELL_VARS%hw_us = hw_us
       CELL_VARS%gw_us = gw_us
       CELL_VARS%dshw_us = dshw_us
       CELL_VARS%tkw_us = tkw_us
       CELL_VARS%tkmidw_us = tkmidw_us
       CELL_VARS%epetw_us = epetw_us
       CELL_VARS%rnetd_us = rnetd_us
       CELL_VARS%xled_us = xled_us
       CELL_VARS%hd_us = hd_us
       CELL_VARS%gd_us = gd_us
       CELL_VARS%dshd_us = dshd_us
       CELL_VARS%tkd_us = tkd_us
       CELL_VARS%tkmidd_us = tkmidd_us
       CELL_VARS%rnet_pot_moss = rnet_pot_moss
       CELL_VARS%xle_p_moss = xle_p_moss
       CELL_VARS%h_p_moss = h_p_moss
       CELL_VARS%g_p_moss = g_p_moss
       CELL_VARS%tk_p_moss = tk_p_moss
       CELL_VARS%tkmid_p_moss = tkmid_p_moss
       CELL_VARS%tskin_p_moss = tskin_p_moss
       CELL_VARS%eact_moss = eact_moss
       CELL_VARS%tsoilold = tsoilold
       CELL_VARS%tkmidpet_us = tkmidpet_us
       CELL_VARS%tkmidpet_moss = tkmidpet_moss
       CELL_VARS%dspet_us = dspet_us
       CELL_VARS%dspet_moss = dspet_moss
       CELL_VARS%rib_moss = rib_moss
       CELL_VARS%epet_moss = epet_moss
       CELL_VARS%xleact_us = xleact_us
       CELL_VARS%hact_us = hact_us
       CELL_VARS%gact_us = gact_us
       CELL_VARS%evtact_us = evtact_us
       CELL_VARS%ievcon_moss = ievcon_moss
       CELL_VARS%bsdew_moss = bsdew_moss
       CELL_VARS%evtact_moss = evtact_moss
       CELL_VARS%rnact_moss = rnact_moss
       CELL_VARS%xleact_moss = xleact_moss
       CELL_VARS%hact_moss = hact_moss
       CELL_VARS%gact_moss = gact_moss
       CELL_VARS%rnact_us = rnact_us
       CELL_VARS%ievcon_us = ievcon_us
       CELL_VARS%f1par = f1par
       CELL_VARS%f3vpd = f3vpd
       CELL_VARS%f4temp = f4temp
       CELL_VARS%f1par_us = f1par_us
       CELL_VARS%f3vpd_us = f3vpd_us
       CELL_VARS%f4temp_us = f4temp_us
       CELL_VARS%f1 = f1
       CELL_VARS%f2 = f2
       CELL_VARS%f3 = f3
       CELL_VARS%f3vpdpar_us = f3vpdpar_us
       CELL_VARS%f4temppar_us = f4temppar_us
       CELL_VARS%roa = roa
       CELL_VARS%roa_ic = roa_ic
       CELL_VARS%ravd = ravd
       CELL_VARS%rahd = rahd
       CELL_VARS%ravd_us = ravd_us
       CELL_VARS%rahd_us = rahd_us
       CELL_VARS%rav_moss = rav_moss
       CELL_VARS%rah_moss = rah_moss
       CELL_VARS%RaSnow = RaSnow
       CELL_VARS%rib_us = rib_us
       CELL_VARS%ravw = ravw
       CELL_VARS%ravw_us = ravw_us
       CELL_VARS%rahw = rahw
       CELL_VARS%rahw_us = rahw_us
       CELL_VARS%rzsm_u = rzsm_u
       CELL_VARS%tzsm_u = tzsm_u
       CELL_VARS%r_mossmold = r_mossmold
       CELL_VARS%deltrz = deltrz
       CELL_VARS%dc_us = dc_us
       CELL_VARS%fw_us = fw_us
       CELL_VARS%dewrun = dewrun



       call land(ipix,i,&
! Energy fluxes

      bsdew,&

       CELL_VARS,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

       tcel = CELL_VARS%tcel
       vppa = CELL_VARS%vppa
       psychr = CELL_VARS%psychr
       xlhv = CELL_VARS%xlhv
       tkel = CELL_VARS%tkel
       appa = CELL_VARS%appa
       vpsat = CELL_VARS%vpsat
       tcel_ic = CELL_VARS%tcel_ic
       vppa_ic = CELL_VARS%vppa_ic
       psychr_ic = CELL_VARS%psychr_ic
       xlhv_ic = CELL_VARS%xlhv_ic
       tkel_ic = CELL_VARS%tkel_ic
       vpsat_ic = CELL_VARS%vpsat_ic
       twet_ic = CELL_VARS%twet_ic
       twet = CELL_VARS%twet
       qv = CELL_VARS%qv
       qv_ic = CELL_VARS%qv_ic
       ra = CELL_VARS%ra
       ra_ic = CELL_VARS%ra_ic
       tkmid_us = CELL_VARS%tkmid_us
       tkact_us = CELL_VARS%tkact_us
       tskinact_moss = CELL_VARS%tskinact_moss
       tkact_moss = CELL_VARS%tkact_moss
       tkmid_moss = CELL_VARS%tkmid_moss
       epetd = CELL_VARS%epetd
       epetd_us = CELL_VARS%epetd_us
       dshact_moss = CELL_VARS%dshact_moss
       xle_act_moss = CELL_VARS%xle_act_moss
       tskinactd_moss = CELL_VARS%tskinact_moss
       tkactd_moss = CELL_VARS%tkactd_moss
       tkmidactd_moss = CELL_VARS%tkmidactd_moss
       ds_p_moss = CELL_VARS%ds_p_moss
       dshact_us = CELL_VARS%dshact_us
       rnetw_us = CELL_VARS%rnetw_us
       xlew_us = CELL_VARS%xlew_us
       hw_us = CELL_VARS%hw_us
       gw_us = CELL_VARS%gw_us
       dshw_us = CELL_VARS%dshw_us
       tkw_us = CELL_VARS%tkw_us
       tkmidw_us = CELL_VARS%tkmidw_us
       epetw_us = CELL_VARS%epetw_us
       rnetd_us = CELL_VARS%rnetd_us
       xled_us = CELL_VARS%xled_us
       hd_us = CELL_VARS%hd_us
       gd_us = CELL_VARS%gd_us
       dshd_us = CELL_VARS%dshd_us
       tkd_us = CELL_VARS%tkd_us
       tkmidd_us = CELL_VARS%tkmidd_us
       rnet_pot_moss = CELL_VARS%rnet_pot_moss
       xle_p_moss = CELL_VARS%xle_p_moss
       h_p_moss = CELL_VARS%h_p_moss
       g_p_moss = CELL_VARS%g_p_moss
       tk_p_moss = CELL_VARS%tk_p_moss
       tkmid_p_moss = CELL_VARS%tkmid_p_moss
       tskin_p_moss = CELL_VARS%tskin_p_moss
       eact_moss = CELL_VARS%eact_moss
       tsoilold = CELL_VARS%tsoilold
       tkmidpet_us = CELL_VARS%tkmidpet_us
       tkmidpet_moss = CELL_VARS%tkmidpet_moss
       dspet_us = CELL_VARS%dspet_us
       dspet_moss = CELL_VARS%dspet_moss
       rib_moss = CELL_VARS%rib_moss
       epet_moss = CELL_VARS%epet_moss
       xleact_us = CELL_VARS%xleact_us
       hact_us = CELL_VARS%hact_us 
       gact_us = CELL_VARS%gact_us  
       evtact_us = CELL_VARS%evtact_us
       ievcon_moss = CELL_VARS%ievcon_moss
       bsdew_moss = CELL_VARS%bsdew_moss
       evtact_moss = CELL_VARS%evtact_moss
       rnact_moss = CELL_VARS%rnact_moss
       xleact_moss = CELL_VARS%xleact_moss
       hact_moss = CELL_VARS%hact_moss
       gact_moss = CELL_VARS%gact_moss
       rnact_us = CELL_VARS%rnact_us
       ievcon_us = CELL_VARS%ievcon_us
       f1par = CELL_VARS%f1par
       f3vpd = CELL_VARS%f3vpd
       f4temp = CELL_VARS%f4temp
       f1par_us = CELL_VARS%f1par_us
       f3vpd_us = CELL_VARS%f3vpd_us
       f4temp_us = CELL_VARS%f4temp_us
       f1 = CELL_VARS%f1
       f2 = CELL_VARS%f2
       f3 = CELL_VARS%f3
       f3vpdpar_us = CELL_VARS%f3vpdpar_us
       f4temppar_us = CELL_VARS%f4temppar_us
       roa = CELL_VARS%roa
       roa_ic = CELL_VARS%roa_ic
       ravd = CELL_VARS%ravd
       rahd = CELL_VARS%rahd
       ravd_us = CELL_VARS%ravd_us
       rahd_us = CELL_VARS%rahd_us
       rav_moss = CELL_VARS%rav_moss
       rah_moss = CELL_VARS%rah_moss
       RaSnow = CELL_VARS%RaSnow
       rib_us = CELL_VARS%rib_us
       ravw = CELL_VARS%ravw
       ravw_us = CELL_VARS%ravw_us
       rahw = CELL_VARS%rahw
       rahw_us = CELL_VARS%rahw_us
       rzsm_u = CELL_VARS%rzsm_u
       tzsm_u = CELL_VARS%tzsm_u
       r_mossmold = CELL_VARS%r_mossmold
       deltrz = CELL_VARS%deltrz
       dc_us = CELL_VARS%dc_us
       fw_us = CELL_VARS%fw_us
       dewrun = CELL_VARS%dewrun
       
! ====================================================================
! Calculate the density and depth of the snow layers.
! ====================================================================

         if ( (GRID_VARS%Swq.gt.0.d0) ) then

            call calcrain (tcel,snow,rain,GRID_VARS%precip_o,GLOBAL%dt)
            call snow_density(GRID_VARS%dsty,snow,tcel,GRID_VARS%Swq,GRID_VARS%Sdepth,GRID_VARS%TSurf,GLOBAL%dt)

         else

           GRID_VARS%Sdepth=0.d0
           GRID_VARS%dsty=0.d0

         endif

      endif

! ====================================================================
! In the vegetation type is lower than zero then solve the open 
! water energy and water balance.
! ====================================================================

      if (GRID_VEG%ivgtyp.eq.(-1)) then

! ....................................................................
! Calculate the deep soil temperature.
! ....................................................................

         if ( (GRID_SOIL%amp.eq.(0.d0)).and.&
              (GRID_SOIL%phase.eq.(0.d0)).and.&
              (GRID_SOIL%shift.eq.(0.d0)) ) then

           GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep

         else

            rrr=real(i)

            GRID_SOIL%Tdeepstep=GRID_SOIL%tdeep + GRID_SOIL%amp*cos ( rrr*GRID_SOIL%phase - GRID_SOIL%shift )

         endif

      endif

      return

      end subroutine Update_Cell

END MODULE MODULE_CELL
