! ====================================================================
!
!			subroutine infilt
!
! ====================================================================
!
! Subroutine infilt calculates the actual rate of infiltration 
! and the rate of infiltration excess runoff.
!
! ====================================================================   

      subroutine infilt(pnet,i_moss,i_und,PackWater_us,SurfWater_us,Swq_us,&
       Outflow_us,dt,PackWater,SurfWater,Swq,Outflow,&
       istmst,cuminf,inc_frozen,rzsmst,rzsm,rzsm_u,thetas,thetar,&
       tolinf,sorp,xk0,psic,bcgamm,bcbeta,deltrz,cc,&
       zw,xinact,satxr,xinfxr,intstm,xinfcp,runtot,irntyp)

      implicit none
      include "help/infilt.h"

! ====================================================================
! Calculate the precipitation input to the soil column.
! ====================================================================

      precipitation=pnet

      if ( (i_moss+i_und) .gt.0) then

         if ((PackWater_us+SurfWater_us+Swq_us).gt.(0.d0))then

! --------------------------------------------------------------------
! In case of snow on top of the under story : the precipitation is
! the liquid water outflow from the snow pack.
! --------------------------------------------------------------------

            precipitation=Outflow_us/dt

         endif

      else

         if ((PackWater+SurfWater+Swq).gt.(0.d0))then

! --------------------------------------------------------------------
! In case of snow on top of the over story : the precipitation is
! the liquid water outflow from the snow pack.
! --------------------------------------------------------------------

            precipitation=Outflow/dt

         endif

      endif

      if(istmst.eq.1)then

! ====================================================================
! If this is the first step of storm event: reset cumulative    
! infiltration, initial root zone soil moisture, sorptivity
! and dimensionless gravity parameter, cc.
! ====================================================================

         if (inc_frozen.eq.0) then

! ....................................................................
! Option 1 : Treat frozen water as liquid water.
! ....................................................................

            call reset_inf_pars(cuminf,zero,rzsmst,rzsm,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

         else

! ....................................................................
! Option 2 : Treat frozen water as soil particles.
! ....................................................................

            call reset_inf_pars(cuminf,zero,rzsmst,rzsm_u,thetas,tolinf,&
       sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)

         endif

      endif

! ====================================================================
! If this is the first step of the storm event all the infiltration
! parameters are calculated now, if this is not the first step
! the parameters are calculated from the previous timestep.
! ====================================================================

      if ((zw-psic).le.zero) then

! ====================================================================
! If surface is saturated then set infiltration to zero.  Also    
! calclulate saturation excess runoff and set infiltration
! excess to zero.
! ====================================================================

         xinact = zero
         satxr = precipitation
         xinfxr = zero

      else

! ====================================================================
! If surface is unsaturated then calculate infiltration.  Set
! saturation excess runoff to zero and calculate infiltration
! excess runoff.
! ====================================================================

         if (intstm.eq.1) then

! --------------------------------------------------------------------
! In interstorm period let bare soil infiltration rate
! equal zero, and let vegetated soil infiltration rate
! equal pnetms to allow for surface wetting due to throughfall
! of condensation.
! --------------------------------------------------------------------

            xinact = precipitation

         else if(cuminf.le.zero)then

! --------------------------------------------------------------------
! If first time step with infiltration then all precipitation is
! infiltrated to avoid division by zero.
! --------------------------------------------------------------------

            xinact = precipitation

         else

! --------------------------------------------------------------------
! Calculate infitration capacity as a function of cumulative      
! infiltration for all other time steps of storm.
! --------------------------------------------------------------------

            xinfcp = cc*xk0*(one+(one/(((one+((four*cc*xk0*cuminf)/&
                                     (sorp**two)))**0.5d0)-one)))

! --------------------------------------------------------------------
! Take the actual infiltration rate as the minimum of the 
! precipitation rate or the infiltration capacity.
! --------------------------------------------------------------------

            if (precipitation.gt.xinfcp) then

               xinact = xinfcp

            else

               xinact = precipitation

            endif

         endif

! --------------------------------------------------------------------
! Calculate infiltration excess runoff; set saturation 
! excess runoff to zero.
! --------------------------------------------------------------------

         xinfxr = precipitation - xinact
         satxr = zero

      endif

! ====================================================================
! Set the value of flag used in output image for type of
! runoff and find total runoff.
! ====================================================================

      runtot = precipitation - xinact
      irntyp = 0

      if (xinfxr.gt.zero) irntyp=1
      if (satxr.gt.zero) irntyp=2

      return

      end
