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
      GLOBAL%mul_fac = 1.0d0

! TEMPORARY LOCATION TO PASS STRUCTURE INFORMATION TO OLD FORMAT

!Removal Causes Failure
f3vpdpar = GRID_VEG%f3vpdpar
f4temppar = GRID_VEG%f4temppar
sesq = GRID_VARS%sesq
xintst = GRID_VARS%xintst
smpet0 = GLOBAL%smpet0

!Global variables
iopgveg = GLOBAL%iopgveg
iopsmini = GLOBAL%iopsmini
ikopt = GLOBAL%ikopt
irestype = GLOBAL%irestype
ioppet = GLOBAL%ioppet
iopveg = GLOBAL%iopveg
iopstab = GLOBAL%iopstab
iopwv = GLOBAL%iopwv
i_2l = GLOBAL%i_2l
newstorm = GLOBAL%newstorm
mul_fac = GLOBAL%mul_fac

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
!Constants
row = GRID_VARS%row!row
roi = GRID_VARS%roi!roi

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

       call atmos(ipix,i,i_2l,GRID_VEG,&

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
       f3vpdpar,f3vpdpar_us,f4temppar,f4temppar_us,&

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
       

       call land(newstorm,ipix,i,i_2l,&
! Meteorological data

       tcel,vppa,psychr,xlhv,tkel,appa,&
       vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&

! Temperature variables

       tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,tkmidpet_us,tkmidpet_moss,tsoilold,&

! Energy fluxes

       epetd,bsdew,&
       epetd_us,dshact_moss,xle_act_moss,rnetd,&
       tkd,tkmidd,ievcon_us,rnact_us,&
       xleact_us,hact_us,gact_us,dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,evtact_us,rnetd_us,xled_us,hd_us,gd_us,dshd_us,&
       tkd_us,tkmidd_us,ievcon_moss,bsdew_moss,evtact_moss,rnet_pot_moss,&
       xle_p_moss,h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,&
       tskin_p_moss,eact_moss,rnact_moss,xleact_moss,hact_moss,gact_moss,&
       ds_p_moss,&

! Moss parameters

       r_moss_depth,thetas_moss,srespar1_moss,srespar2_moss,srespar3_moss,&
       eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,&
       a_ice_moss,b_ice_moss,bulk_dens_moss,&

! Vegetation parameters

       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,&
       f1,f2,f3,&

! Constants

       roa,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,&

! Water balance variables

       rzsm_u,tzsm_u,r_mossmold,&
       deltrz,dc_us,fw_us,dewrun,&

! Different option paramters

       iopgveg,iopsmini,ikopt,&
       irestype,ioppet,iopveg,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

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
