! ====================================================================
!
!			subroutine sursat
!
! ====================================================================
!
! Define land surface saturation states for the region.
!
! ====================================================================

      subroutine sursat(fwcat,fwreg,zw,psic,fw,mul_fac,pr3sat,zrzmax,perrg2,&
       rzsm1,thetas,pr2sat,pr2uns,perrg1,tzsm1,pr1sat,pr1rzs,pr1tzs,pr1uns,&
       satxr,persxr,xinfxr,perixr,ievcon,persac,peruac,perusc,i,ipix)

      implicit none
      include "help/sursat.h"
      integer :: i,ipix

      data tolsat / 0.0001d0 /

! ====================================================================
! Calculate regional fraction of wet canopy.
! ====================================================================

      !fwcat = fwcat + fw*mul_fac
      fwcat = fw*mul_fac
      !fwreg = fwreg + fw*mul_fac
      fwreg = fw*mul_fac

! ====================================================================
! Define land surface saturation states for the region:
!
! Region 3:  Water Table at surface.
! Region 2:  Water Table in root zone.
! Region 1:  Water Table below root zone.
! ====================================================================

      if ((zw-psic).le.zero) then

! --------------------------------------------------------------------&
! First find areas in region 3.
! --------------------------------------------------------------------&

         !pr3sat = pr3sat + one*mul_fac
         pr3sat = one*mul_fac

      else if (((zw-psic).lt.zrzmax).and.((zw-psic).gt.zero)) then

! --------------------------------------------------------------------&
! For all pixels not in area 3 : first see if the water table is
! in the root zone and if the root zone is not saturated (region 2).
! --------------------------------------------------------------------&

         !perrg2 = perrg2 + one*mul_fac
         perrg2 = one*mul_fac

         if (rzsm1.ge.thetas-tolsat) then

            !pr2sat = pr2sat + one*mul_fac
            pr2sat = one*mul_fac

         else

            !pr2uns = pr2uns + one*mul_fac
            pr2uns = one*mul_fac

         endif

      else

! --------------------------------------------------------------------&
! If a pixel is not in in region 3 or 2 it has to be in region 1.
! Ssplit into four possibilities for root and transmission zone
! saturation:
! --------------------------------------------------------------------&

         !perrg1 = perrg1 + one
         perrg1 = one

         if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then

! ....................................................................
! 1) Boot root and transmsission zone are saturated.
! ....................................................................

            !pr1sat = pr1sat + one*mul_fac
            pr1sat = one*mul_fac

         else if ((rzsm1.ge.thetas-tolsat).and.(tzsm1.lt.thetas-tolsat)) then

! ....................................................................
! 2) Root zone is saturated and transmsission zone is not
!    saturated.
! ....................................................................

            !pr1rzs = pr1rzs + one*mul_fac
            pr1rzs = one*mul_fac

         else if ((rzsm1.lt.thetas-tolsat).and.(tzsm1.ge.thetas-tolsat)) then 

! ....................................................................
! 3) Root zone is not saturated and transmsission zone is
!    saturated.
! ....................................................................

            !pr1tzs = pr1tzs + one*mul_fac
            pr1tzs = one*mul_fac

         else

! ....................................................................
! 4) Both root and transmsission zone are not saturated.
! ....................................................................

            !pr1uns = pr1uns + one*mul_fac
            pr1uns = one*mul_fac

         endif

      endif

! ====================================================================
! Determine fractions of land surface contribtuting saturation
! or infiltration excess runoff.
! ====================================================================

      if (satxr.gt.zero) then

         !persxr = persxr + one*mul_fac
         persxr = one*mul_fac

      else if (xinfxr.gt.zero) then

         !perixr = perixr + one*mul_fac
         perixr = one*mul_fac

      endif

! ====================================================================
! Determine areal fractions of bare soil evaporation
! controls - check for atmospheri! contolled (saturated),&
! atmospheri! contolled (unsaturated) and soil controlled.
! ====================================================================

      if (ievcon.eq.3) then

         !persac = persac + one*mul_fac
         persac = one*mul_fac

      else if (ievcon.eq.2) then

         !peruac = peruac + one*mul_fac
         peruac = one*mul_fac

      else 

         !perusc = perusc + one*mul_fac
         perusc = one*mul_fac

      endif

      return

      end
