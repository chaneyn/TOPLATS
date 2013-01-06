MODULE MODULE_CANOPY

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES

USE MODULE_SHARED

USE MODULE_SNOW

implicit none

type CANOPY_template
  real*8 wc
  real*8 zero,one
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

  subroutine canopy(ipix,epetw,&
rnetd,xled,hd,gd,rnetw,xlew,hw,gw,&
tkd,tkmidd,dshd,tkw,tkmidw,dshw,&
GRID_VARS, GRID_VEG, GRID_MET, GLOBAL)

  implicit none
  !include "help/canopy.h"

  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (GRID_MET_template) :: GRID_MET
  type (GLOBAL_template) :: GLOBAL
  type (CANOPY_template) :: GRID_CANOPY

  integer ipix
  real*8 epetw
  real*8 rnetd,xled,hd,gd,rnetw,xlew,hw,gw
  real*8 tkd,tkmidd,dshd,tkw,tkmidw,dshw

  GRID_CANOPY%zero = 0.0d0
  GRID_CANOPY%one = 1.0d0

!epetw-Shared with MODULE_ATMOS, used in MODULE_CELL
!rnetd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!xled-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!hd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!gd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!rnetw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!xlew-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!hw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!gw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!tkd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!tkmidd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!dshd-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!tkw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!tkmidw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL
!dshw-Shared with MODULE_ATMOS, MODULE_LAND, used in MODULE_CELL

! ====================================================================
! Initialize interception storage depth to value from previous time
! step.
! ====================================================================

  GRID_CANOPY%wc = GRID_VARS%wcip1

! ====================================================================
! Set fraction of wet canopy which is considered 100% if canopy
! is wet and 0% if canopy is dry.
! ====================================================================

  call calcfw(GRID_VARS%Swq,GRID_CANOPY%wc,GRID_CANOPY%zero,GRID_VARS%fw,GRID_VEG%wsc)

! ====================================================================
! If potential evaporation is negative, dew forms over the
! whole canopy.  Set dc to zero for this case.
! ====================================================================

  call calcdc(GRID_VARS%dc,GRID_CANOPY%one,epetw,GRID_CANOPY%zero)

! ====================================================================
! Calculate evaporation from the wet canopy
! ====================================================================

  call calcepw(GRID_VARS%epwms,epetw,GRID_CANOPY%one,&
  GRID_VARS%dc,GRID_VARS%fw,GLOBAL%dt,GRID_CANOPY%wc)

! ====================================================================
! Calculate through fall of rainfall.  This is the part of the
! rainfall that can get through the canopy to the underlying soil
! ====================================================================
 
  GRID_VARS%pnet = GRID_CANOPY%zero 

  if ((GRID_MET%pptms-GRID_VARS%epwms)*GLOBAL%dt.gt.(GRID_VEG%wsc-GRID_CANOPY%wc)) then

    GRID_VARS%pnet = (GRID_MET%pptms-GRID_VARS%epwms)-((GRID_VEG%wsc-GRID_CANOPY%wc)/GLOBAL%dt)

  endif

! ====================================================================
! Perform water balance on canopy storage, calculate the new
! interception storage.
! ====================================================================

  GRID_VARS%wcip1 = GRID_CANOPY%wc + GLOBAL%dt*(GRID_MET%pptms-GRID_VARS%epwms-GRID_VARS%pnet)

! --------------------------------------------------------------------
! Don't allow canopy storage to go below zero.
! --------------------------------------------------------------------

  if (GRID_VARS%wcip1.lt.GRID_CANOPY%zero) then
 
    GRID_VARS%epwms = GRID_VARS%epwms + (GRID_VARS%wcip1)/(GLOBAL%dt)
    GRID_VARS%wcip1 = GRID_CANOPY%zero

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

  GRID_VARS%dswc = GRID_VARS%wcip1-GRID_CANOPY%wc
  GRID_VARS%wcrhs=(GRID_MET%pptms-GRID_VARS%epwms-GRID_VARS%pnet)*(GLOBAL%dt)

! --------------------------------------------------------------------
! Double check : if no rain there is no precipitation input to the
! under story.
! --------------------------------------------------------------------

  if (GRID_MET%pptms.eq.(0.d0)) GRID_VARS%pnet =0.d0

! ====================================================================
! Check if the present time step can be considered an interstorm
! period or not in the calculation of the soil water balance.
! ====================================================================

  call interstorm(ipix,GRID_VARS,GLOBAL)

! ====================================================================
! Add up pet terms of the over story to get average values.
! ====================================================================

  GRID_VARS%rnpet = rnetd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                     rnetw*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))
  GRID_VARS%xlepet = xled*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                      xlew*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))
  GRID_VARS%hpet = hd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                    hw*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))
  GRID_VARS%gpet = gd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                    gw*(GRID_CANOPY%one-GRID_VARS%dc*(one-GRID_VARS%fw))

! --------------------------------------------------------------------
! Solve for temperature and heat storage only when energy balance
! method is used.
! --------------------------------------------------------------------

  if (GLOBAL%ioppet.eq.0) then

    GRID_VARS%tkpet = tkd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                       tkw*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))
    GRID_VARS%tkmidpet = tkmidd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                          tkmidw*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))
    GRID_VARS%dspet = dshd*GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw)+&
                       dshw*(GRID_CANOPY%one-GRID_VARS%dc*(GRID_CANOPY%one-GRID_VARS%fw))

  else

    GRID_VARS%tkpet = GRID_CANOPY%zero
    GRID_VARS%tkmidpet = GRID_CANOPY%zero
    GRID_VARS%dspet = GRID_CANOPY%zero

  endif

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

  subroutine interstorm(ipix,GRID_VARS,GLOBAL)

  implicit none

  integer :: ipix
  real*8 :: r_input

  type (GRID_VARS_template) :: GRID_VARS
  type (GLOBAL_template) :: GLOBAL
  
! ====================================================================
! Calculate the water input to the ground.
! ====================================================================

  if ((GRID_VARS%PackWater+GRID_VARS%SurfWater+GRID_VARS%Swq).gt.(0.001d0)) then

  r_input=GRID_VARS%Outflow

  else

  r_input=GRID_VARS%pnet

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

  GRID_VARS%xintst=GRID_VARS%xintst+GLOBAL%dt

! --------------------------------------------------------------------
! Now check if time since end of ppt is past threshold
! then add one to number of time steps into the
! interstorm period.
! --------------------------------------------------------------------

  if (GRID_VARS%xintst.gt.GLOBAL%endstm)then

  GRID_VARS%intstp=GRID_VARS%intstp+1

! --------------------------------------------------------------------
! Now, if this is the first step in the interstorm period then
! reset storm flags (istmst,istorm) and reset number of
! steps into storm period to zero.
! --------------------------------------------------------------------

  if (GRID_VARS%intstp.eq.1) then

  GRID_VARS%istmst=0
  GRID_VARS%istorm=0
  GRID_VARS%intstm=1

  endif

! --------------------------------------------------------------------
! If time since end of ppt is within threshold then continue
! to add step to the storm period.
! --------------------------------------------------------------------

  else

  GRID_VARS%istmst=GRID_VARS%istmst+1

  endif

! --------------------------------------------------------------------
! If there is precipitation then storm event is in progress --
! increment the number of steps in the storm period and reset
! the time from the end of precipitation to zero.
! --------------------------------------------------------------------

  else

  GRID_VARS%istmst=GRID_VARS%istmst+1
  GRID_VARS%xintst=0.d0

! --------------------------------------------------------------------
! If this is the first time step in the storm period
! then reset the storm flags and the number of time
! steps in the interstorm period.
! --------------------------------------------------------------------

  if (GRID_VARS%istmst.eq.1)then

  GRID_VARS%intstp=0
  GRID_VARS%intstm=0
  GRID_VARS%istorm=1

  endif

  endif

  return

  end subroutine interstorm

END MODULE MODULE_CANOPY
