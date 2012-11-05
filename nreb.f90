! ====================================================================
!
!		subroutine nreb
!
! ====================================================================
!
! Subroutine to calculate the skin temperature at either the
! potential or actual (depending on value of ipetopt)
! evapotranspiration rate using a newton-raphson interation scheme.
!
! ====================================================================
!
!       Rn - G - H - LE = DS
! NOTE: sign convention: all radiative fluxes directed toward the 
!         surface are positive (e.g. net radiation).  All 
!         non-radiative fluxes directed away from surface are 
!         positive (e.g. latent, sensible and soil heat fluxes)
!
! ====================================================================

      subroutine nreb(ipetopt,alb,emiss,thermc1,thermc2,heatcap1,heatcap2,&
       heatcapold,rav,rah,tskink,tmidk,rn,xle,epot,h,g,ds,tairc,vppa,roa,&
       psychr,xlhv,zdeep,tdeep,zmid,rsd,rld,toleb,maxnri,dt,itime)

      implicit none
      include "help/nreb.h"
     
      data cp,row,vkc,sbc,astab,cph2o,GRAV/1005.d0,&
           997.d0,0.4d0,5.6696d-8,1.d0,4186.d0,9.81d0/
      data HMIN/-250.d0/,DTMAX/4.d0/
 

! ====================================================================
! Initialize variable conv to 0 so non-convergence will be recognized.
! ====================================================================

      itmpflg=0
      ihflg=0
      ileflg=0
      iconv=0
      T0new=tskink
      T1new=tmidk

! ====================================================================
! Initialize temperatures and flux derivatives.
! ====================================================================

! --------------------------------------------------------------------
! Save old skin temperature--set to air temp if skin temp unrealistic.
! --------------------------------------------------------------------

      if (tskink.ge.200.and.tskink.le.400) then	

         toldk = tskink
         toldmidk=tmidk

      else

         toldk = tairc+273.15d0

      endif

! --------------------------------------------------------------------
! New Skin temperature initialized as air temperature.
! --------------------------------------------------------------------

      tskinc=tairc
      tskink=tairc+273.15d0

! --------------------------------------------------------------------
! Save old storage term.
! --------------------------------------------------------------------

      if (ds.ge.-500.and.ds.le.500) then	

         dsold=ds

      else

         dsold=0.d0

      endif

! --------------------------------------------------------------------
! Sensible Heat Flux Derivative.
! --------------------------------------------------------------------

      if (ihflg.eq.0) then

         dhdt=((roa*cp)/rah)

      else

         dhdt=0.d0

      endif

! --------------------------------------------------------------------
! Ground Heat Flux Derivative.
! --------------------------------------------------------------------

! ....................................................................
! Calculate denominator for gradient
! ....................................................................

      dzdeep = (zdeep - zmid)
	 
! ....................................................................
! Xu Liang Formula with input thermodynami! parameters
! ....................................................................

      gdenom= thermc1*dzdeep*dt*2.d0 + thermc2*zmid*dt*2.d0 +&
              heatcap2*dzdeep*dzdeep*zmid
      dgdt = (thermc1*thermc2*dt*2.d0 &
              + thermc1*heatcap2*dzdeep*dzdeep)/&
              gdenom

! --------------------------------------------------------------------
! Ground Heat Storage Derivative.
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Set ds to zero--save expressions below.  Appears still missing some
! storage effect.
! --------------------------------------------------------------------

      ddsdt = 0.d0

! ====================================================================
! Iterate up to maxnri times.
! ====================================================================

      do 100 iter=1,maxnri

! --------------------------------------------------------------------
! Calculate the vapor pressure deficit and temperature 
! difference at beginning of each iteration.
! --------------------------------------------------------------------

         vpsat=611.d0*dexp((17.27d0*tskinc)/(237.3d0+tskinc))
         vpdef=vpsat-vppa
         dftfac=(237.d0*17.27d0)/((237.3d0+tskinc)**two)

! --------------------------------------------------------------------
! Calculate fluxes with current skin temperature.
! --------------------------------------------------------------------

! ....................................................................
! Net radiation.
! ....................................................................

          rntmp=rsd*(one-alb)+rld-emiss*sbc*(tskink**four)

! ....................................................................
! Net radiation Derivative.
! ....................................................................

         drndt=-four*emiss*sbc*(tskink**three)

! ....................................................................
! Latent Heat Flux.
! Compute potential energy balance if ipetopt==1, otherwise, use input
! soil- or vegetation- controlled" value.
! ....................................................................

	 if (ipetopt.eq.1) then

            xletmp=((roa*cp)/(psychr*rav))*(vpdef)

	 else

            xletmp=xle

	 endif

! ....................................................................
! Latent Heat Flux Derivative.
! ....................................................................

	 if (ipetopt.eq.1) then

            dxledt=((roa*cp)/(psychr*rav))*vpsat*dftfac

	 else

            dxledt=0.d0

	 endif

! ....................................................................
! Sensible Heat Flux.
! ....................................................................

         if (ihflg.eq.0) htmp=((roa*cp)/rah)*(tskinc-tairc)

! ....................................................................
! Ground Heat Flux.
! Xu Liang Formula with input thermodynami! parameters.
! ....................................................................

         gtmp = (thermc1*thermc2*dt*2.d0*(tskink-tdeep)+&
                thermc1*heatcap2*dzdeep*dzdeep*(tskink-tmidk))/&
                gdenom

! ....................................................................
! Ground Heat Storage.
! ....................................................................

         tmidknew = (heatcap2*tmidk/(2.d0*dt) +&
                    0.5*gtmp/dzdeep + thermc2*tdeep/(dzdeep*dzdeep))/&
                    (heatcap2/(2.d0*dt) + thermc2/(dzdeep*dzdeep))
         dstmp= -zmid*heatcap1*(tmidknew-tmidk)/dt

! --------------------------------------------------------------------
! Calculate energy balance with current skin temperature
! and latent heat flux.
! --------------------------------------------------------------------

         ft= rntmp - xletmp - htmp - gtmp - dstmp

! --------------------------------------------------------------------
! Calculate the derivative of the energy balance w.r.t. 
! skin temperature.
! --------------------------------------------------------------------

         dftdt= drndt - dxledt - dhdt - dgdt - ddsdt

! --------------------------------------------------------------------
! Calculate skin temperature adjustment and update
! for next interation.
! --------------------------------------------------------------------

         deltnr=-ft/dftdt

	 if (ihflg.eq.0.and.ileflg.eq.0.and.itmpflg.eq.0) then

            tskink=tskink+deltnr
            tskinc=tskink-273.15d0

	 endif

! --------------------------------------------------------------------
! Check for convergence 
! --------------------------------------------------------------------
!
         if (abs(deltnr).lt.toleb.and.iconv.eq.1) then

           goto 200

         endif

         if (abs(deltnr).lt.toleb.and.iconv.eq.0) then

           iconv=1

         endif

100   continue

! ====================================================================
! No convergence.
! ====================================================================

      if (iconv.eq.0)  then

         print*,'nreb: no convergence for time step',itime
         print*,'thermc1 =             ',thermc1
         print*,'thermc2 =             ',thermc2
         print*,'heatcap1 =             ',heatcap1
         print*,'heatcap2 =             ',heatcap2
         print*,'dt =                  ',dt
         print*,'tairk =             ',tairc + 273.15
         print*,'tskink =             ',tskink
         print*,'tmidk =             ',tmidk
         print*,'tdeep =             ',tdeep
         print*,'gdenom =             ',gdenom
         print*,'rah =             ',rah
         print*,'rav =             ',rav
         print*,'roa =             ',roa
         print*,'error =             ',deltnr
         print*,'tolerance =         ',toleb
         print*,
         print*,'dhdt =              ',dhdt
         print*,'dgdt =              ',dgdt
         print*,'ddsdt =             ',ddsdt
         print*,'dzdeep =             ',dzdeep
         print*,'gdenom =             ',gdenom
         print*,'gtmp =              ',gtmp
         print*,'rsd =               ',rsd
         print*,'rld =               ',rld
         print*,

	 stop

      endif

! ====================================================================
! Reset energy balance fluxes 
! ====================================================================

200   rn = rntmp

      if (ipetopt.eq.1) then

         xle = xletmp
         epot=xle/(xlhv*row)

       endif

       h = htmp
       g = gtmp 
       ds= rn - xle - h - g 

! ====================================================================
! Update mid soil temperature given ground heat flux.
! ====================================================================

      if (dabs(thermc1).gt.1e-10) then

         tmidk = (heatcap2*tmidk/(2.d0*dt) +&
               g/dzdeep + thermc2*tdeep/(dzdeep*dzdeep))/&
               (heatcap2/(2.d0*dt) + thermc2/(dzdeep*dzdeep))

      else

         tmidk = toldmidk

      endif

! ====================================================================
! If necessary, do convergence time backstepping.
! ====================================================================

      if ( (tskink.ge.0.).and.(tskink.le.533.).and.&
           (tmidk.ge.0.).and.(tmidk.le.533.) ) then

         h=h

      else

         write (*,*) 'NREB convergence error '
         write (*,*) g,h,xle,tskink,tmidknew,ds
         write (*,*)
         write (*,*) ipetopt,alb,emiss,thermc1,thermc2,heatcap1,heatcap2,&
       heatcapold,rav,rah,tskink,tmidk,rn,xle,epot,h,g,ds,tairc,vppa,roa,&
       psychr,xlhv,zdeep,tdeep,zmid,rsd,rld,toleb,maxnri,dt,itime

         deltnew=dt/2.d0
         if (deltnew.le.0.5) stop
!cw!         write (*,*) 'Solution time backstepping algorithm ',deltnew

!cw!         call nreb(ipetopt,alb,emiss,thermc1,thermc2,heatcap1,heatcap2,&
!cw!       heatcapold,rav,rah,T0new,T1new,rn,xle,epot,h,g,ds,tairc,vppa,roa,&
!cw!       psychr,xlhv,zdeep,tdeep,zmid,rsd,rld,toleb,maxnri,deltnew,itime)

!cw!         call nreb(ipetopt,alb,emiss,thermc1,thermc2,heatcap1,heatcap2,&
!cw!       heatcapold,rav,rah,T0new,T1new,rn,xle,epot,h,g,ds,tairc,vppa,roa,&
!cw!       psychr,xlhv,zdeep,tdeep,zmid,rsd,rld,toleb,maxnri,deltnew,itime)

         tskink=T0new
         tmidknew=T1new
         tmidk=T1new

      endif

      return

      end
