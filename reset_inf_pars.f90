! ====================================================================
!
!                       subroutine reset_inf_pars
!
! ====================================================================
!
! Reset cumulative infiltration, initial root zone soil moisture,&
! sorptivity and dimensionless gravity parameter, cc.
!
! ====================================================================

      subroutine reset_inf_pars(cuminf,zero,rzsmst,rzsm,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

      implicit none
      include "help/reset_inf_pars.h"

! ====================================================================
! Calculate the root zone soil moisture
! ====================================================================

      cuminf = zero

      rzsmst = rzsm

      if (rzsmst.ge.thetas)  then

         rzsmst=thetas-tolinf

      endif

! ====================================================================
! Calculate sorptivity and gravity parameters at the first step
! of the storm event.
! ====================================================================

      sorp = (((two*xk0*((thetas-rzsmst)**two)*psic)&
                   /(thetas-thetar))*((one/(bcgamm+0.5d0*bcbeta-one))+&
                    ((thetas-thetar)/ (thetas-rzsmst))))**0.5d0

      deltrz = rzsmst-thetar
      if(deltrz.le.zero) deltrz=zero
      cc = 0.5d0*(one+ ((deltrz/(thetas-thetar))**(bcgamm/bcbeta)))

      return

      end
