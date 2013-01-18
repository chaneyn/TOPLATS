MODULE MODULE_CANOPY

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
! Subroutine to calculate interception by vegetation and calculate
! canopy water balance.
!
! ====================================================================

  subroutine canopy(ipix,GRID_VARS, GRID_VEG, GRID_MET, GLOBAL)

  implicit none

  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (GRID_MET_template) :: GRID_MET
  type (GLOBAL_template) :: GLOBAL
  type (CANOPY_template) :: GRID_CANOPY

  integer :: ipix
  
  real*8 :: zero,one
  zero = 0.0d0
  one = 1.0d0

! ====================================================================
! Initialize interception storage depth to value from previous time
! step.
! ====================================================================

  GRID_CANOPY%wc = GRID_VARS%wcip1

! ====================================================================
! Set fraction of wet canopy which is considered 100% if canopy
! is wet and 0% if canopy is dry.
! ====================================================================

  call calcfw(GRID_VARS,GRID_VEG,GRID_CANOPY)

! ====================================================================
! If potential evaporation is negative, dew forms over the
! whole canopy.  Set dc to zero for this case.
! ====================================================================

  call calcdc(GRID_VARS)

! ====================================================================
! Calculate evaporation from the wet canopy
! ====================================================================

  call calcepw(GRID_VARS,GRID_CANOPY,GLOBAL)

! ====================================================================
! Calculate canopy rainfall and water balance
! ====================================================================
 
  call calcwt(GRID_VARS,GRID_VEG,GRID_MET,GLOBAL,GRID_CANOPY)

! ====================================================================
! Check if the present time step can be considered an interstorm
! period or not in the calculation of the soil water balance.
! ====================================================================

  call interstorm(ipix,GRID_VARS,GLOBAL)

! ====================================================================
! Add up pet terms of the over story to get average values.
! ====================================================================

  call calcnet(GRID_VARS,GRID_VEG,GLOBAL)

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

  subroutine calcfw(GRID_VARS,GRID_VEG,GRID_CANOPY)

  implicit none

  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (CANOPY_template) :: GRID_CANOPY
  
  real*8 :: zero,one
  zero = 0.0d0
  one = 1.0d0

  if (GRID_VARS%Swq.le.zero) then

  if (GRID_CANOPY%wc.gt.zero) then

  GRID_VARS%fw = (GRID_CANOPY%wc/GRID_VEG%wsc)**(0.667d0)

  else

  GRID_VARS%fw = zero

  endif

  if (GRID_VEG%wsc.eq.zero) GRID_VARS%fw=zero

  else

  GRID_VARS%fw=one

  endif

  if (GRID_VARS%fw.ge.one) GRID_VARS%fw=one

  if ( (GRID_VARS%fw.ge.zero).and.(GRID_VARS%fw.le.one) ) then

  GRID_VARS%fw=GRID_VARS%fw

  else

  write (*,*) 'CALCFW : fw : ',GRID_VARS%fw
  write (*,*) GRID_VARS%Swq,GRID_CANOPY%wc,GRID_VARS%fw,GRID_VEG%wsc
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

  subroutine calcdc(GRID_VARS)

  implicit none
  
  type (GRID_VARS_template) :: GRID_VARS

  real*8 :: zero,one
  zero = 0.0d0
  one = 1.0d0

  GRID_VARS%dc = one
  if (GRID_VARS%epetw.lt.zero) GRID_VARS%dc=zero

  if ( (GRID_VARS%dc.ge.zero).and.(GRID_VARS%dc.le.one) ) then

  GRID_VARS%dc=GRID_VARS%dc

  else

  write (*,*) 'CALCD! : d! out of bounds ',GRID_VARS%dc
  if (GRID_VARS%dc.lt.zero) GRID_VARS%dc=zero
  if (GRID_VARS%dc.gt.one) GRID_VARS%dc=one

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

  subroutine calcepw(GRID_VARS,GRID_CANOPY,GLOBAL)

  implicit none
  
  type (GRID_VARS_template) :: GRID_VARS
  type (CANOPY_template) :: GRID_CANOPY
  type (GLOBAL_template) :: GLOBAL

  real*8 :: zero,one
  zero = 0.0d0
  one = 1.0d0

  GRID_VARS%epwms = GRID_VARS%epetw * (one-GRID_VARS%dc*(one-GRID_VARS%fw))

  if ((GRID_VARS%epwms*GLOBAL%dt).gt.GRID_CANOPY%wc) then
 
  GRID_VARS%fw = GRID_VARS%fw*GRID_CANOPY%wc/(GRID_VARS%epwms*GLOBAL%dt)
  GRID_VARS%epwms=GRID_VARS%epetw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
 
  endif

  if ( (GRID_VARS%fw.ge.zero).and.(GRID_VARS%fw.le.one) ) then

  GRID_VARS%fw=GRID_VARS%fw
 
  else

  write (*,*) 'CALEPW : fw : ',GRID_VARS%fw
  write (*,*) GRID_VARS%epwms,GRID_VARS%epetw,GRID_VARS%dc,GRID_VARS%fw,GLOBAL%dt,GRID_CANOPY%wc
  stop

  endif

  return

  end subroutine calcepw

! ====================================================================
!
!                   subroutine calcwt
!
! ====================================================================
!
! Calculate canopy rainfall and water balance
!
! ====================================================================

  subroutine calcwt(GRID_VARS,GRID_VEG,GRID_MET,GLOBAL,GRID_CANOPY)

  implicit none
 
  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (GRID_MET_template) :: GRID_MET
  type (GLOBAL_template) :: GLOBAL 
  type (CANOPY_template) :: GRID_CANOPY

  real*8 :: zero
  zero = 0.0d0

! --------------------------------------------------------------------
! Calculate through fall of rainfall.  This is the part of the
! rainfall that can get through the canopy to the underlying soil
! --------------------------------------------------------------------
  GRID_VARS%pnet = zero
  if (((GRID_MET%pptms-GRID_VARS%epwms)*GLOBAL%dt).gt.(GRID_VEG%wsc-GRID_CANOPY%wc)) then
    GRID_VARS%pnet=(GRID_MET%pptms-GRID_VARS%epwms)-((GRID_VEG%wsc-GRID_CANOPY%wc)/GLOBAL%dt)
  endif

! --------------------------------------------------------------------
! Perform water balance on canopy storage, calculate the new
! interception storage.
! --------------------------------------------------------------------
  GRID_VARS%wcip1=GRID_CANOPY%wc+GLOBAL%dt*(GRID_MET%pptms-GRID_VARS%epwms-GRID_VARS%pnet)

! --------------------------------------------------------------------
! Don't allow canopy storage to go below zero.
! --------------------------------------------------------------------
  if (GRID_VARS%wcip1.lt.zero) then
    GRID_VARS%epwms=GRID_VARS%epwms+GRID_VARS%wcip1/GLOBAL%dt
    GRID_VARS%wcip1=zero
  endif

! --------------------------------------------------------------------
! Calculate the precipitation that will go to the overstory
! layer and that will not fall through.
! This is the precipitation input for the snow melt model for the
! over story. 
! --------------------------------------------------------------------
  GRID_VARS%precip_o=GRID_MET%pptms

! --------------------------------------------------------------------
! Check canopy water balance, calculate the change in water storage.
! --------------------------------------------------------------------
  GRID_VARS%dswc=GRID_VARS%wcip1-GRID_CANOPY%wc
  GRID_VARS%wcrhs=(GRID_MET%pptms-GRID_VARS%epwms-GRID_VARS%pnet)*GLOBAL%dt

! --------------------------------------------------------------------
! Double check : if no rain there is no precipitation input to the
! under story.
! --------------------------------------------------------------------
  if (GRID_MET%pptms.eq.zero) GRID_VARS%pnet=zero

  return

  end subroutine calcwt

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

! ====================================================================
!
!                   subroutine calcnet
!
! ====================================================================
!
! Add up pet terms of the over story to get average values.
!
! ====================================================================
  subroutine calcnet(GRID_VARS,GRID_VEG,GLOBAL)

  implicit none

  type (GRID_VARS_template) :: GRID_VARS
  type (GRID_VEG_template) :: GRID_VEG
  type (GLOBAL_template) :: GLOBAL

  real*8 :: zero,one
  zero=0.0d0
  one=1.0d0

  GRID_VARS%rnpet = GRID_VEG%rnetd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                     GRID_VEG%rnetw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
  GRID_VARS%xlepet = GRID_VEG%xled*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                      GRID_VEG%xlew*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
  GRID_VARS%hpet = GRID_VEG%hd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                    GRID_VEG%hw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
  GRID_VARS%gpet = GRID_VEG%gd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                    GRID_VEG%gw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))

! --------------------------------------------------------------------
! Solve for temperature and heat storage only when energy balance
! method is used.
! --------------------------------------------------------------------

  if (GLOBAL%ioppet.eq.0) then

    GRID_VARS%tkpet = GRID_VEG%tkd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                       GRID_VEG%tkw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
    GRID_VARS%tkmidpet = GRID_VEG%tkmidd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                          GRID_VEG%tkmidw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))
    GRID_VARS%dspet = GRID_VEG%dshd*GRID_VARS%dc*(one-GRID_VARS%fw)+&
                       GRID_VEG%dshw*(one-GRID_VARS%dc*(one-GRID_VARS%fw))

  else

    GRID_VARS%tkpet = zero
    GRID_VARS%tkmidpet = zero
    GRID_VARS%dspet = zero

  endif

  return 

  end subroutine calcnet  

end module MODULE_CANOPY
