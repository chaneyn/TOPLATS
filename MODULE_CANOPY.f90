MODULE MODULE_CANOPY

USE VARIABLES

implicit none

contains
!====================================================================
!
!                   subroutine Canopy
!
! ====================================================================
  subroutine canopy(ipix,dt,wc,wcip1,Swq,fw,wsc,dc,epetw,epwms,pnet,&
        pptms,precip_o,dswc,wcrhs,endstm,xintst,intstp,istmst,istorm,&
        intstm,Outflow,PackWater,SurfWater,rnpet,xlepet,hpet,gpet,&
        rnetd,xled,hd,gd,rnetw,xlew,hw,gw,ioppet,tkpet,tkmidpet,dspet,&
        tkd,tkmidd,dshd,tkw,tkmidw,dshw)
  implicit none
 
  real*8 Swq,wc,wsc,zero
  real*8 fw
  integer ipix
  integer intstp,istmst,istorm,intstm,ioppet
  real*8 wcip1
  real*8 dc,epetw
  real*8 epwms,pnet,pptms,precip_o
  real*8 dswc,wcrhs
  real*8 endstm
  real*8 xintst,Outflow,PackWater,SurfWater,rnpet,xlepet,hpet
  real*8 gpet,rnetd,xled,hd,gd,rnetw,xlew,hw,gw,tkpet,tkmidpet
  real*8 dspet,tkd,tkmidd,dshd,tkw,tkmidw,dshw
  real*8 one,two,three,four,five,six,dt,dummy

!    include "help/canopy.h"

! --------------------------------------------------------------------
! Declare variables and classes
! --------------------------------------------------------------------
    !data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0/
    !print*,'Calculates canopy'

! --------------------------------------------------------------------
! Initialize interception storage depth to value from previous time
! step.
! --------------------------------------------------------------------
    wc = wcip1

! --------------------------------------------------------------------
! Set fraction of wet canopy which is considered 100% if canopy
! is wet and 0% if canopy is dry.
! --------------------------------------------------------------------
    !Subroutine calcfw
    call calcfw(Swq, wc, wsc, zero, fw)
    ! Flow_in: Swq, wc, wsc, zero
    ! Flow_out: fw

! --------------------------------------------------------------------
! If potential evaporation is negative, dew forms over the
! whole canopy.  Set dc to zero for this case.
! --------------------------------------------------------------------
    !Subroutine calcdc
    call calcdc(epetw, zero, one, dc)
    ! Dew_in: epetw, zero, one
    ! Dew_out: dc

! --------------------------------------------------------------------
! Calculate evaporation from the wet canopy
! --------------------------------------------------------------------
    !Subroutine calcepw
    call calcepw(epetw, dc, dt, wc, fw, epwms)
    ! Evapo_in: epetw, dc, dt, wc, fw
    ! Evapo_out: epwms, fw

! --------------------------------------------------------------------
! Calculate canopy rainfall and water balance
! --------------------------------------------------------------------
    !Subroutine calcewt
    call calcwt(pptms, epwms, dt, wsc, wc, zero,&
                pnet, wcip1, precip_o, wcrhs)
    ! Water_in: pptms, epwms, dt, wsc, wc, zero
    ! Water_out: pnet, wcip1, epwms, precip_o, precip_o, wcrhs

! --------------------------------------------------------------------
! Check if the present time step can be considered an interstorm
! period or not in the calculation of the soil water balance.
! --------------------------------------------------------------------
    !Subroutine interstorm
    call interstorm(PackWater+SurfWater+Swq, Outflow, pnet, dt, endstm,&
                    xintst, intstp, istmst, istorm, intstm)
    ! interstorm_in: Outflow, pnet, dt, endstm, xintst, intstp, istmst
    ! interstorm_out: xintst, intstp, istmst, istorm, intstm

! --------------------------------------------------------------------
! Add up pet terms of the over story to get average values.
! --------------------------------------------------------------------
    !Subroutine calcnet
    call calcnet(ioppet, dc, fw, zero, one, rnetd, rnetw, xled, xlew, hd, hw,&
                 gd, gw, tkd, tkw, tkmidd, tkmidw, dshd, dshw,&
                 rnpet, xlepet, hpet, gpet, tkpet, tkmidpet, dspet)
    ! Average_in: ioppet, dc, fw, zero, one, rnetd, rnetw, xled, xlew, hd, hw,
    !             gd, gw, tkd, tkw, tkmidd, tkmidw, dshd, dshw
    ! Average_out: rnpet, xlepet, hpet, gpet, tkpet, tkmidpet, dspet
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
  subroutine calcfw(Swq, wc, wsc, zero, fw)
  implicit none
  !Sets fraction of wet canopy

    real*8 Swq, wc, wsc, zero, fw

    !print*,'    --Sets fraction of wet canopy'

    if ((Swq).le.(zero)) then

      if ((wc).gt.(zero)) then
        fw = (wc/wsc)**(0.667d0)
      else
        fw = zero
      endif

      if ((wsc).eq.(zero)) fw = zero
    else
      fw = 1.d0
    endif

    if (fw.ge.1.d0) fw = 1.d0

    if ( ((fw).ge.(0.d0)).and.((fw).le.(1.d0)) ) then
      fw = fw
    else
      write (*,*) 'CALCFW : fw = ',fw
      write (*,*) Swq,wc,wsc,fw
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
  subroutine calcdc(epetw, zero, one, dc)
  implicit none
  !Calculates dew formed over the whole canopy

    real*8 epetw, zero, one
    real*8 dc

    !print*,'    --Calculates dew formed over the whole canopy'

    dc = one
    if ((epetw).lt.(zero)) dc = zero

    if ( ((dc).ge.(0.d0)).and.((dc).le.(1.d0)) ) then
      dc = dc
    else
      write (*,*) 'CALCD! : d! out of bounds ',dc
      if ((dc).lt.(0.d0)) dc = zero
      if ((dc).gt.(1.d0)) dc = one
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
  subroutine calcepw(epetw, dc, dt, wc, fw, epwms)
  implicit none
  !Calculates evaporation from the wet canopy

    real*8 epetw, dc, dt, wc, fw, epwms

    !print*,'    --Calculates evaporation from the wet canopy'

    epwms = (epetw)*(one-(dc)*(one-fw))

    if (((epwms)*(dt)).gt.(wc)) then
      fw = (fw)*(wc)/((epwms)*(dt))
      epwms = (epetw)*(one-(dc)*(one-fw))
    endif

    if ( ((fw).ge.(0.d0)).and.((fw).le.(1.d0)) ) then
      fw = fw
    else
      write (*,*) 'CALEPW : fw = ',fw
      write (*,*) epwms,fw
      write (*,*) epetw,dc,dt,wc,one
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
! Input variable(ipix) deleted.
!
! ===================================================================
  subroutine interstorm(snowp, outf, pnet, dt, endstm,&
                        xintst, intstp, istmst, istorm, intstm)
  implicit none
  !Checks if the soil under snow is treated as interstorm

    real*8 snowp, outf, pnet, dt, endstm, xintst
    integer intstp, istmst, istorm, intstm

    real*8 r_input

    !print*,'    --Checks if the soil under snow is treated as interstorm'

    ! Calculate the water input to the ground.
    if ((snowp).gt.(0.001d0)) then
      r_input = outf
    else
      r_input = pnet
    endif

    ! Define storm and interstorm events
    if ((r_input).le.(0.d0))then

      xintst = xintst+dt
      if ((xintst).gt.(endstm))then
        intstp = intstp+1
        if ((intstp).eq.(1)) then
          istmst = 0
          istorm = 0
          intstm = 1
        endif
      else
        istmst = istmst+1
      endif

    else
      istmst = istmst+1
      xintst = 0.d0

      if ((istmst).eq.(1))then
        intstp = 0
        intstm = 0
        istorm = 1
      endif

    endif

  return
  end subroutine interstorm

! ====================================================================
!
!                   subroutine calcwt
!
! ====================================================================
!
! Calculate canopy rainfall and water balance
!
! ====================================================================
subroutine calcwt(pptms, epwms, dt, wsc, wc, zero,&
                  pnet, wcip1, precip_o, wcrhs)
  implicit none
  !Calculates precipitation and canopy water storage

    real*8 pptms, epwms, dt, wsc, wc, zero
    real*8 pnet, wcip1, precip_o, wcrhs

    !print*,'    --Calculates precipitation and canopy water storage'

    ! Calculate through fall of rainfall.  This is the part of the
    ! rainfall that can get through the canopy to the underlying soil
    pnet = zero
    if (((pptms-epwms)*(dt)).gt.(wsc-wc)) then
      pnet = (pptms-epwms)-((wsc-wc)/dt)
    endif

    ! Perform water balance on canopy storage, calculate the new
    ! interception storage.
    wcip1 = wc + (dt)*(pptms-epwms-pnet)

    ! Don't allow canopy storage to go below zero.
    if ((wcip1).lt.(zero)) then
      epwms = epwms + wcip1/dt
      wcip1 = zero
    endif

    ! Calculate the precipitation that will go to the overstory
    ! layer and that will not fall through.
    ! This is the precipitation input for the snow melt model for the
    ! over story.
    precip_o = pptms

    ! Check canopy water balance, calculate the change in water storage.
    dswc = wcip1-wc
    wcrhs=(pptms-epwms-pnet)*(dt)

    ! Double check : if no rain there is no precipitation input to the
    ! under story.
    if ((pptms).eq.(0.d0)) pnet = 0.d0

  return
  end subroutine calcwt

! ====================================================================
!
!                   subroutine calcnet
!
! ====================================================================
!
! Add up pet terms of the over story to get average values.
!
! ====================================================================
subroutine calcnet(ioppet, dc, fw, zero, one, rnetd, rnetw, xled, xlew,&
                   hd, hw, gd, gw, tkd, tkw, tkmidd, tkmidw, dshd, dshw,&
                   rnpet, xlepet, hpet, gpet, tkpet, tkmidpet, dspet)
  implicit none
  !Adds up overstory PET and calculates averages

    integer ioppet
    real*8 dc, fw, zero, one, rnetd, rnetw, xlew, xled
    real*8 hw, hd, gw, gd, tkw, tkd, tkmidd, tkmidw, dshd, dshw
    real*8 rnpet, xlepet, hpet, gpet, tkpet, tkmidpet, dspet

    real*8 tmpd,tmpw

    !print*,'    --Adds up overstory PET and calculates averages'

    tmpd = (dc)*(one-fw)
    tmpw = one-(dc)*(one-fw)

    rnpet = (rnetd)*tmpd + (rnetw)*tmpw
    xlepet = (xled)*tmpd + (xlew)*tmpw
    hpet = (hd)*tmpd + (hw)*tmpw
    gpet = (gd)*tmpd + (gw)*tmpw

    ! Solve for temperature and heat storage only when energy balance
    ! method is used.
    if ((ioppet).eq.(0)) then
      tkpet = (tkd)*tmpd + (tkw)*tmpw
      tkmidpet = (tkmidd)*tmpd + (tkmidw)*tmpw
      dspet = (dshd)*tmpd + (dshw)*tmpw
    else
      tkpet = zero
      tkmidpet = zero
      dspet = zero
    endif

  return
  end subroutine calcnet

end module MODULE_CANOPY

