! ====================================================================
! Calculate the soil temperatures given that a snow pack is on top
! of the soil.
! ====================================================================

      subroutine nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
                          tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp)

      implicit none
      include "help/nreb_snow.h"

! --------------------------------------------------------------------
! The skin temperature is assumed to be equal to the temperature
! at the bottom of the snow pack.
! --------------------------------------------------------------------

      tskink=tairc
      tmidk=tmidknew

! --------------------------------------------------------------------
! Calculate the distance between the bottom 2 nodes.
! --------------------------------------------------------------------

      dzdeep = (zdeep - zmid)

! --------------------------------------------------------------------
! Calculate the ground heat flux.
! --------------------------------------------------------------------

      gdenom = thermc1*dzdeep*dt*2.d0 + thermc2*zmid*dt*2.d0 +&
               heatcap2*dzdeep*dzdeep*zmid

      gtmp = (thermc1*thermc2*dt*2.d0*(tskink-tdeep)+&
             thermc1*heatcap2*dzdeep*dzdeep*(tskink-tmidk))/&
             gdenom

! --------------------------------------------------------------------
! Calculate the mid soil temperature.
! --------------------------------------------------------------------
      
      tmidknew = (heatcap2*tmidk/(2.d0*dt) +&
               0.5*gtmp/dzdeep + thermc2*tdeep/(dzdeep*dzdeep))/&
               (heatcap2/(2.d0*dt) + thermc2/(dzdeep*dzdeep))

      if ( (tskink.le.(0.)).or.(tskink.ge.(343.)).or.&
           (tmidknew.le.(10.)).or.(tmidknew.ge.(303.)) ) then

         write (*,*) 'Error in the energy balance under snow '
         write (*,*) thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
                     tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp
         stop

      endif

      return

      end
