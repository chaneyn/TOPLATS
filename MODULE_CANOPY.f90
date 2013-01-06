MODULE MODULE_CANOPY

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES

USE MODULE_SHARED

USE MODULE_SNOW

implicit none

type CANOPY_template
  real*8 wc
end type CANOPY_template

contains

! ====================================================================
!
! 			subroutine canopy
!
! ====================================================================
!
!Subroutine to calculate interception by vegetation and calculate
! canopy water balance.
!
! ====================================================================

  subroutine canopy(ipix,wc,fw,dc,epetw,epwms,pnet,intstm,&
rnetd,xled,hd,gd,rnetw,xlew,hw,gw,&
tkd,tkmidd,dshd,tkw,tkmidw,dshw,&
GRID_VARS, GRID_VEG, GRID_MET, GLOBAL )

  implicit none
  include "help/canopy.h"

  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (GRID_MET_template) :: GRID_MET
  type (GLOBAL_template) :: GLOBAL
  type (CANOPY_template) :: GRID_CANOPY

!dt = GLOBAL%dt
wc = GRID_CANOPY%wc
!wcip1 = GRID_VARS%wcip1
!Swq = GRID_VARS%Swq
fw = GRID_VARS%fw
!wsc = GRID_VEG%wsc
dc = GRID_VARS%dc
!epetw-Shared with MODULE_ATMOS
epwms = GRID_VARS%epwms
pnet = GRID_VARS%pnet
!pptms = GRID_MET%pptms
!precip_o = GRID_VARS%precip_o
!dswc = GRID_VARS%dswc
!wcrhs = GRID_VARS%wcrhs
!endstm = GLOBAL%endstm
!xintst = GRID_VARS%xintst
!intstp = GRID_VARS%intstp
!istmst = GRID_VARS%istmst
!istorm = GRID_VARS%istorm
intstm = GRID_VARS%intstm
!Outflow = GRID_VARS%Outflow
!PackWater = GRID_VARS%PackWater
!SurfWater = GRID_VARS%SurfWater
!rnpet = GRID_VARS%rnpet
!xlepet = GRID_VARS%xlepet
!hpet = GRID_VARS%hpet
!gpet = GRID_VARS%gpet
!rnetd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!xled-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!hd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!gd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!rnetw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!xlew-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!hw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!gw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!ioppet = GLOBAL%ioppet
!tkpet = GRID_VARS%tkpet
!tkmidpet = GRID_VARS%tkmidpet
!dspet = GRID_VARS%dspet
!tkd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkmidd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!dshd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkmidw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!dshw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL

! ====================================================================
! Initialize interception storage depth to value from previous time
! step.
! ====================================================================

  wc = GRID_VARS%wcip1

! ====================================================================
! Set fraction of wet canopy which is considered 100% if canopy
! is wet and 0% if canopy is dry.
! ====================================================================

  call calcfw(GRID_VARS%Swq,wc,zero,fw,GRID_VEG%wsc)

! ====================================================================
! If potential evaporation is negative, dew forms over the
! whole canopy.  Set dc to zero for this case.
! ====================================================================

  call calcdc(dc,one,epetw,zero)

! ====================================================================
! Calculate evaporation from the wet canopy
! ====================================================================

  call calcepw(epwms,epetw,one,dc,fw,GLOBAL%dt,wc)

! ====================================================================
! Calculate through fall of rainfall.  This is the part of the
! rainfall that can get through the canopy to the underlying soil
! ====================================================================
 
  pnet = zero 

  if ((GRID_MET%pptms-epwms)*GLOBAL%dt.gt.(GRID_VEG%wsc-wc)) then

    pnet = (GRID_MET%pptms-epwms)-((GRID_VEG%wsc-wc)/GLOBAL%dt)

  endif

! ====================================================================
! Perform water balance on canopy storage, calculate the new
! interception storage.
! ====================================================================

  GRID_VARS%wcip1 = wc + GLOBAL%dt*(GRID_MET%pptms-epwms-pnet)

! --------------------------------------------------------------------
! Don't allow canopy storage to go below zero.
! --------------------------------------------------------------------

  if (GRID_VARS%wcip1.lt.zero) then
 
    epwms = epwms + (GRID_VARS%wcip1)/(GLOBAL%dt)
    GRID_VARS%wcip1 = zero

  endif

! ====================================================================
! Calculate the precipitation that will go to the overstory
! layer and that will not fall through.
! This is the precipitation input for the snow melt model for the
! over story.
! ====================================================================

  GRID_VARS%precip_o=GRID_MET%pptms

! ====================================================================
! Check canopy water balance, calculate the change in water storage.
! ====================================================================

  GRID_VARS%dswc = GRID_VARS%wcip1-wc
  GRID_VARS%wcrhs=(GRID_MET%pptms-epwms-pnet)*(GLOBAL%dt)

! --------------------------------------------------------------------
! Double check : if no rain there is no precipitation input to the
! under story.
! --------------------------------------------------------------------

  if (GRID_MET%pptms.eq.(0.d0)) pnet =0.d0

! ====================================================================
! Check if the present time step can be considered an interstorm
! period or not in the calculation of the soil water balance.
! ====================================================================

  call interstorm(ipix,pnet,GRID_VARS%Outflow,GRID_VARS%PackWater+GRID_VARS%SurfWater+GRID_VARS%Swq,&
  GRID_VARS%xintst,GLOBAL%dt,GRID_VARS%intstp,GLOBAL%endstm,GRID_VARS%istmst,GRID_VARS%istorm,intstm)

! ====================================================================
! Add up pet terms of the over story to get average values.
! ====================================================================

  GRID_VARS%rnpet = rnetd*dc*(one-fw) + rnetw*(one-dc*(one-fw))
  GRID_VARS%xlepet = xled*dc*(one-fw) + xlew*(one-dc*(one-fw))
  GRID_VARS%hpet = hd*dc*(one-fw) + hw*(one-dc*(one-fw))
  GRID_VARS%gpet = gd*dc*(one-fw) + gw*(one-dc*(one-fw))

! --------------------------------------------------------------------
! Solve for temperature and heat storage only when energy balance
! method is used.
! --------------------------------------------------------------------

  if (GLOBAL%ioppet.eq.0) then

    GRID_VARS%tkpet = tkd*dc*(one-fw) + tkw*(one-dc*(one-fw))
    GRID_VARS%tkmidpet = tkmidd*dc*(one-fw) + tkmidw*(one-dc*(one-fw))
    GRID_VARS%dspet = dshd*dc*(one-fw) + dshw*(one-dc*(one-fw))

  else

    GRID_VARS%tkpet = zero
    GRID_VARS%tkmidpet = zero
    GRID_VARS%dspet = zero

  endif

!GLOBAL%dt = dt
GRID_CANOPY%wc = wc
!GRID_VARS%wcip1 = wcip1
!GRID_VARS%Swq = Swq
GRID_VARS%fw = fw
!GRID_VEG%wsc = wsc
GRID_VARS%dc = dc
!epetw-Shared with MODULE_ATMOS
GRID_VARS%epwms = epwms
GRID_VARS%pnet = pnet
!GRID_MET%pptms = pptms
!GRID_VARS%precip_o = precip_o
!GRID_VARS%dswc = dswc
!GRID_VARS%wcrhs = wcrhs
!GLOBAL%endstm = endstm
!GRID_VARS%xintst = xintst
!GRID_VARS%intstp = intstp
!GRID_VARS%istmst = istmst
!GRID_VARS%istorm = istorm
GRID_VARS%intstm = intstm
!GRID_VARS%Outflow = Outflow
!GRID_VARS%PackWater = PackWater
!GRID_VARS%SurfWater = SurfWater
!GRID_VARS%rnpet = rnpet
!GRID_VARS%xlepet = xlepet
!GRID_VARS%hpet = hpet
!GRID_VARS%gpet = gpet
!rnetd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!xled-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!hd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!gd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!rnetw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!xlew-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!hw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!gw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!GLOBAL%ioppet = ioppet
!GRID_VARS%tkpet = tkpet
!GRID_VARS%tkmidpet = tkmidpet
!GRID_VARS%dspet = dspet 
!tkd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkmidd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!dshd-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!tkmidw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL
!dshw-Shared with MODULE_ATMOS, MODULE_LAND, MODULE_CELL

  return

  end subroutine canopy

!====================================================================
!
!                   subroutine calcfw
!
! ====================================================================
!
! Set fraction of wet canopy which is considered 100% if canopy
! is wet and 0% if canopy is dry.
!
! ====================================================================

  subroutine calcfw(Swq,wc,zero,fw,wsc)

  implicit none
  !include "help/calcfw.h"
  real*8 Swq,wc,zero,fw,wsc

  if (Swq.le.zero) then

  if (wc.gt.zero) then

  fw = (wc/wsc)**(0.667d0)

  else

  fw = zero

  endif

  if (wsc.eq.zero) fw=zero

  else

  fw=1.d0

  endif

  if (fw.ge.1.d0) fw=1.d0

  if ( (fw.ge.0.d0).and.(fw.le.1.d0) ) then

  fw=fw

  else

  write (*,*) 'CALCFW : fw : ',fw
  write (*,*) Swq,wc,zero,fw,wsc
  stop

  endif

  return

  end subroutine calcfw

! ====================================================================
!
!                  subroutine calcdc
!
! ====================================================================
!
! If potential evaporation is negative, dew forms over the
! whole canopy.
!
! ====================================================================

  subroutine calcdc(dc,one,epetw,zero)

  implicit none
  !include "help/calcdc.h"
  real*8 dc,one,epetw,zero

  dc = one
  if (epetw.lt.zero) dc=zero

  if ( (dc.ge.0.d0).and.(dc.le.1.d0) ) then

  dc=dc

  else

  write (*,*) 'CALCD! : d! out of bounds ',dc
  if (dc.lt.0.d0) dc=zero
  if (dc.gt.1.d0) dc=one

  endif

  return

  end subroutine calcdc

! ====================================================================
!
!                   subroutine calcepw
!
! ====================================================================
!
! Calculate evaporation from the wet canopy.
!
! ====================================================================

  subroutine calcepw(epwms,epetw,one,dc,fw,dt,wc)

  implicit none
  !include "help/calcepw.h"
  real*8 epwms,epetw,one,dc,fw,dt,wc

  epwms = epetw * (one-dc*(one-fw))

  if ((epwms*dt).gt.wc) then
 
  fw = fw*wc/(epwms*dt)
  epwms=epetw*(one-dc*(one-fw))

  endif

  if ( (fw.ge.0.d0).and.(fw.le.1.d0) ) then

  fw=fw

  else

  write (*,*) 'CALEPW : fw : ',fw
  write (*,*) epwms,epetw,one,dc,fw,dt,wc
  stop

  endif

  return

  end subroutine calcepw

! ====================================================================
!
!                   subroutine interstorm
!
! ====================================================================
!
! This subroutine checks if the soil under snow is treated as   *
! interstorm or storm period
!
! ====================================================================

  subroutine interstorm(ipix,precipi,outf,snowp,xintst,&
dt,intstp,endstm,istmst,istorm,intstm)

  implicit none
  !include "help/interstorm.h"
  integer :: ipix,intstp,istmst,istorm,intstm
  real*8 :: precipi,outf,snowp,xintst,dt,endstm,r_input

! ====================================================================
! Calculate the water input to the ground.
! ====================================================================

  if (snowp.gt.(0.001d0)) then

  r_input=outf

  else

  r_input=precipi

  endif

! ====================================================================
! Define storm and interstorm events
! ====================================================================

! --------------------------------------------------------------------
! First if there is no precipitation this time step then add
! to the time since ppt ended.  Then check if time since end
! of ppt is greater than threshold which defines beginning of
! interstorm period.
! --------------------------------------------------------------------

  if (r_input.le.(0.d0))then

  xintst=xintst+dt

! --------------------------------------------------------------------
! Now check if time since end of ppt is past threshold
! then add one to number of time steps into the
! interstorm period.
! --------------------------------------------------------------------

  if (xintst.gt.endstm)then

  intstp=intstp+1

! --------------------------------------------------------------------
! Now, if this is the first step in the interstorm period then
! reset storm flags (istmst,istorm) and reset number of
! steps into storm period to zero.
! --------------------------------------------------------------------

  if (intstp.eq.1) then

  istmst=0
  istorm=0
  intstm=1

  endif

! --------------------------------------------------------------------
! If time since end of ppt is within threshold then continue
! to add step to the storm period.
! --------------------------------------------------------------------

  else

  istmst=istmst+1

  endif

! --------------------------------------------------------------------
! If there is precipitation then storm event is in progress --
! increment the number of steps in the storm period and reset
! the time from the end of precipitation to zero.
! --------------------------------------------------------------------

  else

  istmst=istmst+1
  xintst=0.d0

! --------------------------------------------------------------------
! If this is the first time step in the storm period
! then reset the storm flags and the number of time
! steps in the interstorm period.
! --------------------------------------------------------------------

  if (istmst.eq.1)then

  intstp=0
  intstm=0
  istorm=1

  endif

  endif

  return

  end subroutine interstorm

END MODULE MODULE_CANOPY
