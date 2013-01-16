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
         GRID(ilandc)%VEG,GRID(ipix)%VARS,GRID(ipix)%VARS%wcip1,&
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
               GRID_VARS,wcip1,CAT,GLOBAL)

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

!VEGETATION
xlai = GRID_VEG%xlai
emiss = GRID_VEG%emiss
zpd = GRID_VEG%zpd
z0m = GRID_VEG%z0m
z0h = GRID_VEG%z0h
rescan = GRID_VEG%rescan
tc = GRID_VEG%tc
tw = GRID_VEG%tw
rtact = GRID_VEG%rtact
rtdens = GRID_VEG%rtdens
psicri = GRID_VEG%psicri
respla = GRID_VEG%respla
rsmax = GRID_VEG%rsmax
Rpl = GRID_VEG%Rpl
f3vpdpar = GRID_VEG%f3vpdpar
rsmin = GRID_VEG%rsmin
trefk = GRID_VEG%trefk
f4temppar = GRID_VEG%f4temppar
canclos = GRID_VEG%canclos
extinct = GRID_VEG%extinct
albd = GRID_VEG%albd
albw = GRID_VEG%albw
zww = GRID_VEG%zww
za = GRID_VEG%za
Tslope1 = GRID_VEG%Tslope1
Tint1 = GRID_VEG%Tint1
Tslope2 = GRID_VEG%Tslope2
Tint2 = GRID_VEG%Tint2
Tsep = GRID_VEG%Tsep
Twslope1 = GRID_VEG%Twslope1
Twint1 = GRID_VEG%Twint1
Twslope2 = GRID_VEG%Twslope2
Twint2 = GRID_VEG%Twint2
Twsep = GRID_VEG%Twsep
wsc = GRID_VEG%wsc
tcbeta = GRID_VEG%tcbeta
ivgtyp = GRID_VEG%ivgtyp
i_und = GRID_VEG%i_und

!SOIL PROPERTIES
bcbeta = GRID_SOIL%bcbeta
psic = GRID_SOIL%psic
thetas = GRID_SOIL%thetas
thetar = GRID_SOIL%thetar
xk0 = GRID_SOIL%xk0
zdeep = GRID_SOIL%zdeep
tdeep = GRID_SOIL%tdeep
zmid = GRID_SOIL%zmid
tmid0 = GRID_SOIL%tmid0
rocpsoil = GRID_SOIL%rocpsoil
quartz = GRID_SOIL%quartz
srespar1 = GRID_SOIL%srespar1
srespar2 = GRID_SOIL%srespar2
srespar3 = GRID_SOIL%srespar3
a_ice = GRID_SOIL%a_ice
b_ice = GRID_SOIL%b_ice
bulk_dens = GRID_SOIL%bulk_dens
amp = GRID_SOIL%amp
phase = GRID_SOIL%phase
shift = GRID_SOIL%shift
bcgamm = GRID_SOIL%bcgamm
par = GRID_SOIL%par
corr = GRID_SOIL%corr
idifind = GRID_SOIL%idifind
par = GRID_SOIL%par
ifcoarse = GRID_SOIL%ifcoarse

!METEOROLOGY
rsd = GRID_MET%rsd
rld = GRID_MET%rld
tdry = GRID_MET%tdry
rh = GRID_MET%rh
uzw = GRID_MET%uzw
press = GRID_MET%press
pptms = GRID_MET%pptms

!GRID VARIABLES
!Water Balance Variables
rzsm = GRID_VARS%rzsm
tzsm = GRID_VARS%tzsm
rzsm1 = GRID_VARS%rzsm1
tzsm1 = GRID_VARS%tzsm1
rzsm_f = GRID_VARS%rzsm_f
tzsm_f = GRID_VARS%tzsm_f
rzsm1_f = GRID_VARS%rzsm1_f
tzsm1_f = GRID_VARS%tzsm1_f
rzdthetaidt = GRID_VARS%rzdthetaidt
tzdthetaidt = GRID_VARS%tzdthetaidt
rzdthetaudtemp = GRID_VARS%rzdthetaudtemp
pnet = GRID_VARS%pnet
xinact = GRID_VARS%xinact
runtot = GRID_VARS%runtot
irntyp = GRID_VARS%irntyp
!Meteorological Variables
Tincan = GRID_VARS%Tincan
rh_ic = GRID_VARS%rh_ic
precip_o = GRID_VARS%precip_o
precip_u = GRID_VARS%precip_u
!Temperature Variables
tkmid = GRID_VARS%tkmid
tkact = GRID_VARS%tkact
tkmidpet = GRID_VARS%tkmidpet
tkpet = GRID_VARS%tkpet
!Energy Fluxes
dshact = GRID_VARS%dshact
rnetpn = GRID_VARS%rnetpn
gbspen = GRID_VARS%gbspen
evtact = GRID_VARS%evtact
ievcon = GRID_VARS%ievcon
gact = GRID_VARS%gact
rnact = GRID_VARS%rnact
xleact = GRID_VARS%xleact
hact = GRID_VARS%hact
ebspot = GRID_VARS%ebspot
dspet = GRID_VARS%dspet
rnpet = GRID_VARS%rnpet
xlepet = GRID_VARS%xlepet
hpet = GRID_VARS%hpet
gpet = GRID_VARS%gpet
rib = GRID_VARS%rib
!Snow
PackWater = GRID_VARS%Packwater
SurfWater = GRID_VARS%SurfWater
VaporMassFlux = GRID_VARS%VaporMassFlux
r_MeltEnergy = GRID_VARS%r_MeltEnergy
Outflow = GRID_VARS%Outflow
PackWater_us = GRID_VARS%PackWater_us
Swq = GRID_VARS%Swq
TPack = GRID_VARS%Tpack
TSurf = GRID_VARS%TSurf
xleact_snow = GRID_VARS%xleact_snow
hact_snow = GRID_VARS%hact_snow
rn_snow = GRID_VARS%rn_snow
Swq_us = GRID_VARS%Swq_us
TPack_us = GRID_VARS%TPack_us
TSurf_us = GRID_VARS%TSurf_us
xleact_snow_us = GRID_VARS%xleact_snow_us
hact_snow_us = GRID_VARS%hact_snow_us
rn_snow_us = GRID_VARS%rn_snow_us
dens = GRID_VARS%dens
dens_us = GRID_VARS%dens_us
dsty = GRID_VARS%dsty
dsty_us = GRID_VARS%dsty_us
Sdepth = GRID_VARS%Sdepth
Sdepth_us = GRID_VARS%Sdepth_us

!wcip1 = GRID_VARS%wcip1
cuminf = GRID_VARS%cuminf
sorp = GRID_VARS%sorp
cc = GRID_VARS%cc
sesq = GRID_VARS%sesq
alb_snow = GRID_VARS%alb_snow
istmst = GRID_VARS%istmst
intstm = GRID_VARS%intstm
intstp = GRID_VARS%intstp
istorm = GRID_VARS%istorm
xintst = GRID_VARS%xintst
atanb = GRID_VARS%atanb

!Global variables
dt = GLOBAL%dt
inc_frozen = GLOBAL%inc_frozen
zrzmax = GLOBAL%zrzmax
toleb = GLOBAL%toleb
maxnri = GLOBAL%maxnri
smpet0 = GLOBAL%smpet0
endstm = GLOBAL%endstm
iopthermc = GLOBAL%iopthermc
iopgveg = GLOBAL%iopgveg
iopthermc_v = GLOBAL%iopthermc_v
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

!Point Data
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
cph2o = GRID_VARS%cph2o!cph2o
cp = GRID_VARS%cp!cp
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
       

       call land(newstorm,ipix,i,dt,inc_frozen,i_2l,&
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

! Soil parameters
       zrzmax,&

! Moss parameters

       r_moss_depth,thetas_moss,srespar1_moss,srespar2_moss,srespar3_moss,&
       eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,&
       a_ice_moss,b_ice_moss,bulk_dens_moss,&

! Vegetation parameters

       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,&
       f1,f2,f3,&

! Constants

       roa,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,&

! Water balance variables

       rzsm_u,tzsm_u,r_mossmold,&
       deltrz,dc_us,fw_us,dewrun,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopsmini,ikopt,&
       irestype,ioppet,iopveg,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

! ====================================================================
! Calculate the density and depth of the snow layers.
! ====================================================================

         if ( (GRID_VARS%Swq.gt.0.d0) ) then

            call calcrain (tcel,snow,rain,GRID_VARS%precip_o,dt)
            call snow_density(GRID_VARS%dsty,snow,tcel,GRID_VARS%Swq,GRID_VARS%Sdepth,GRID_VARS%TSurf,dt)

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
