      subroutine lakemain(lakpix,mixmax,tp_in,hice_in,hsnw_in,xlat_a,xlong_a,&
       tax,avpx,qax,psurfx,uax,rhx,rlwdx,swx,Qe,Qh,eta_a,dt,surface,temp_a,&
       tempi_a,hice_a,hsnow_a,preca_a,fraci_a,precacc,numnod,rnet,i)

! ====================================================================
! Sets up the boundary conditions for simulations with the lake
! module.
!
! Parameters :
!
! lakpix	Lake pixel number (-).
! tp_in		Initial lake skin temperature (C).
! hice_in	Initial ice depth (m).
! hsnw_in	Initial snow depth (m).
! xtime		Number of minutes into the simulation (-).
! xlat_a	Latitiude of the lake (degrees).
! xlong_a	Longitude of the lake (degrees).
! tax		Air temperature (K).
! avpx		Air vapor pressure (Pa).
! Qax		Air humidity (-).
! psurfx	Air pressure (Pa).
! uax		Wind stpeed (m/s).
! rhx		Relative humidity (-).
! rlwdx		Downwelling long wave radiation (W/m2).
! swx		Incoming short wave radiation (W/m2).
! Qe		Latent heat flux (W/m2).
! Qh		Sensible heat flux (W/m2).
! eta_a		Decline of solar radiation input with depth (m-1).
! dt		Time step size (s).
! surface	Area of the lake (m2).
! temp_a	Temperature of the lake water (C).
! tempi_a	Temperature of the lake ice (C).
! hice_a	Height of the lake ice (m).
! hsnow_a	Height of the snow on top of the lake (m).
! preca_a	Precipitation (m/s).
! precacc	Accumulation precipitation (m).
! fraci_a	Fractional coverage of ice (-).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/lakemain.h"

! ====================================================================
! Initialize the lake location parameters and the water freezing
! temperature.
! ====================================================================

      xlat=xlat_a
      xlong=xlong_a

      Tcutoff=0.
      Tcutoff=Tcutoff+273.15

! ====================================================================
! Calculate the net long wave radiation.
! ====================================================================

      ts=tp_in-273.15
      rlu=0.97*delta*tp_in**4.
      rlnet=rlwdx-rlu
      precl=0.d0

! ====================================================================
! Solve the energy budget for the lake.
! ====================================================================

      call lake(lakpix,tax,uax,qax,psurfx,&
                swx,rlwdx,rlnet,mixmax,rhx,Qe,Qh,&
                xlat,xlong,eta_a,dt,surface,temp_a,&
                tempi_a,hice_a,hsnow_a,preca_a,&
                fraci_a,precacc,numnod,rnet)

      return

      end

      subroutine adjflux (sw, tp, tcutoff, hice, hsnow, ta,&
                           qa, ua, psurf, rhosurf, delq, evap, hsen, rlwd )

! ====================================================================
! Adjust the sensible and latent heat fluxes calculated for liquid
! water for an ice cover.
! 
! Parameters :
!
! sw		Incoming solar radiation (W/m2).
! tp		Temperature of the ice (C).
! tcutoff	Temperature for wich below there is no liquid water (C).
! hice		Ice height (m).
! hsnow		Snow height (m).
! ta		Air temperature (C).
! qa		Absolute air humididity (-).
! ua		Wind speed (m/s).
! psurf		Air pressure (Pa).
! rhosurf	Density of the air at the surface (kg/m3).
! delq		Difference in water vapor mixing ratio at the surface
!		of the lake and in the air (-).
! evap		Corrected evaporation over ice (mg/m2 sec).
! hsen		Sensible heat flux over ice (W/m2).
! rlwd		Downwelling long wave radiation (W/m2).
! ====================================================================

      implicit none
      integer :: temporary
      include "help/adjflux.h"

      t4(x)=(x+273.15)**4.
      tcutc=tcutoff-273.15
      ntimes=1
      a=-0.1
      b=10.0

! --------------------------------------------------------------------
! Start the iteration for the correct ice temperature.
! --------------------------------------------------------------------

99    continue

      do 100 temporary=1,601,1!tposs=b,-50.0,a
        tposs = b + (temporary-1)*a
! --------------------------------------------------------------------
! Calculate the energy balance of the lake.
! --------------------------------------------------------------------

         t=tposs
         tlat=t+273.15

         call latsens (tlat, tcutoff, hice, ta,&
                      qa, ua, psurf, rhosurf, delq, evap, qsen )

         hsen=qsen
         qlat=-evap*Lei

! --------------------------------------------------------------------
! Calculate the closure of the energy balance.
! --------------------------------------------------------------------

         qmet=rlwd-emice*stefbol*t4(t)+qsen+qlat

! --------------------------------------------------------------------
! Calculate the radiation balance of the ice.
! --------------------------------------------------------------------

         call icerad (sw, hice, hsnow, condbar, val, val2 )
         q0=-qmet
         t0=condbar*(sw-q0)+tcutc-val

         if ( (t0.ge.t) .and. (ntimes.eq.1) ) then

! --------------------------------------------------------------------
! If the newly calculated ice temperature is larger than the original
! (which means convergence) and if this is the first iteration
! series, make the iteration for the ice temperature vary between
! the new ice temperature and -50 with smaller temperature steps.
! --------------------------------------------------------------------

             b=t-a
             a=-.001
             ntimes=2
             go to 99

         else if ( (t0.ge.t) .and. (ntimes.eq.2) ) then

! --------------------------------------------------------------------
! If the newly calculated ice temperature is larger than the original
! (which means convergence) and if this is the second iteration
! series, make the iteration for the ice temperature vary between
! the new ice temperature and -50 with even smaller temperature steps.
! --------------------------------------------------------------------

             b=t-a
             a=-.00001
             ntimes=3
             go to 99

         else if ( (t0.ge.t) .and. (ntimes.eq.3 )) then

! --------------------------------------------------------------------
! If the newly calculated ice temperature is larger than the original
! (which means convergence) and if this is the third iteration
! series, than we have found the ice temperature up to a precision
! of 0.00001 degrees.  Return this value.
! --------------------------------------------------------------------

            tp=t+273.15

            return

         endif

100   continue

      return

      end

      subroutine alblake ( tair, tcutoff, albs, albi, albw )

! ====================================================================
! Calculate the albedo of snow, ice and water of the lake.
!
! Parameters :
!
! tair		Air temperature (C).
! tcutoff	Temperature at which water freezes (C).
! albs		Snow albedo (-).
! albi		Ice albedo (-).
! albw		Water albedo (-).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/alblake.h"

      albs=0.2
      albi=0.2
      albw=0.1
      albw=0.05

      return

      end

      subroutine colavg (t,ti,fracprv,dnsty,numnod)

! ====================================================================
! Calculate the temperature of the lake water at the different node
! levels and calculate the water density.
!
! Parameters :
!
! t		Water temperature (C).
! ti		Ice temperature (C).
! fracprv	Fraction of lake covered with ice (-).
! densty	Water density at each node (kg/m3).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/colavg.h"

      do 100 j=1,numnod

! --------------------------------------------------------------------
! Calculate the densities of the ice and water fractions.
! --------------------------------------------------------------------

         call density (t(j,1),dnstyw)
         call density (ti(j,1),dnstyi)
         dnstyw=dnstyw+1000.
         dnstyi=dnstyi+1000.

! --------------------------------------------------------------------
! Calculate the specifi! heats of the ice and water fractions.
! --------------------------------------------------------------------

         call specheat (t(j,1),cpw)
         call specheat (ti(j,1),cpi)

! --------------------------------------------------------------------
! Calculate the depth differences.
! --------------------------------------------------------------------

         z=dz_lake
         if (j.eq.1) z=surf

! --------------------------------------------------------------------
! Calculate the lake temperature as a weight of ice and water
! temperatures.
! --------------------------------------------------------------------

         temp=((1.-fracprv)*t(j,1)*z*dnstyw*cpw+&
              fracprv*ti(j,1)*z*dnstyi*cpi)/&
              ((z*dnstyw*cpw+z*dnstyi*cpi)*0.5)

         t(j,1)=temp

! --------------------------------------------------------------------
! Recalculate the water density at the different nodes.
! --------------------------------------------------------------------

         call density (t(j,1),dnsty(j))

100      continue

      return

      end

      subroutine density (ts,rhostps)

! ====================================================================
! Calculate the water density.
!
! Parameters :
!
! ts		Temperature at the node (C).
! rhostps	Water density at the node (kg/m3).
! ====================================================================

      implicit none
      include "help/density.h"

! --------------------------------------------------------------------
! Convert the single precision parameter into double precision.
! --------------------------------------------------------------------

      t=dble(ts)

! --------------------------------------------------------------------
! Calculate the water density as a function of temperature.
! --------------------------------------------------------------------

      rhow=999.842594D0 + 6.793952D-2*t - 9.095290D-3*t**2. +&
          1.001685D-4*t**3. - 1.120083D-6*t**4. + 6.536332D-9*t**5.

      rhostp=rhow

! --------------------------------------------------------------------
! Return the difference between 1000 kg/m3 and the calculated
! density, not the actual value.
! --------------------------------------------------------------------

      rhostps=rhostp-1.D3

      return

      end


      subroutine eddy (iwater,u2,T, dnsty, de, deglat,numnod)

! ====================================================================
! Calculate the eddy diffusivity.
!
! Parameters :
!
! iwater	1 if liquid water, 0 if ice.
! u2		Wind speed at 2 meter (m/s).
! T		Water temperature at the different nodes (C).
! dnsty		Water density at the different nodes (C).
! de		Eddy diffusivity (m2/d).
! deglat	Latitude of the pixel (degrees).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/eddy.h"

      do 10 k=1,numnod

! ====================================================================
! Calculate the density at all nodes.
! ====================================================================

         call density (T(k,1),dnsty(k))

! ====================================================================
! Calculate the distance between the nodes.
! ====================================================================

         zhalf(k) = dz_lake

10    continue

! ====================================================================
! Calculate the distance between the first node and the surface.
! ====================================================================

      zhalf(1) = (surf + dz_lake) * 0.5

! ====================================================================
! If there is ice only molecular diffusivity is taken in account,
! no eddy diffusivity.
! ====================================================================

      if (iwater.ne.1) then

         do 22 k=1,numnod-1

            de(k)=dm

22       continue

         de(numnod)=dm
         return

      endif

! ====================================================================
! Avoid to low wind speeds for computational stability.
! ====================================================================

      if(u2.lt.1.0) u2 = 1.0

! ====================================================================
! Calculate the latitude in radians.
! ====================================================================

      zlat=deglat*pi/180.

! ====================================================================
! Determine the latitudinalyy dependent parameter of the Ekman
! profile, the surface value of the friction velocity and the neutral
! value of the Prandtl number.
! ====================================================================

      ks=6.6*sin(zlat)**0.5*u2**(-1.84)
      ws=0.0012*u2
      Po=1.0

! ====================================================================
! Determine the maximum value of the second term in the calculation
! of the Richardson number for computational stability.
! ====================================================================

      radmax=4.e4

      do 20 k= 1,numnod-1

! ====================================================================
! Calculate the eddy diffusivity for each node.
! ====================================================================

! --------------------------------------------------------------------
! Calculate the derivative of density with depth, the Brunt-Vaisala
! frequency and the depth of the node.
! --------------------------------------------------------------------

         dpdz=(dnsty(k+1)-dnsty(k))/zhalf(k)
         N2=(dpdz/(1.e3+dnsty(k)))*9.8   
         z=surf+float(k-1)*dz_lake 

! --------------------------------------------------------------------
! Calculate the second term in the calculation of the Richardson
! number, make sure this number does not get too large for
! computational stability.  Also make sure this term does not go
! below 1 to avoid negative Richardson numbers.
! --------------------------------------------------------------------
 
         if ((z*exp(ks*z)/ws).gt.1.e8) then

            rad = radmax

         else

            rad=1.+40.*N2*(kv*z)**2./(ws**2.*exp(-ks*z))

            if (rad.gt.radmax) rad=radmax

         endif

         if (rad.lt.1.0) rad=1.0

! --------------------------------------------------------------------
! Calculate the Richardson number and the eddy diffusivity.
! --------------------------------------------------------------------

         Ri=(-1.0+sqrt(rad))/20.0 
         de(k)=dm+kv*ws*z*Po*exp(-ks*z)&
               /(1.0+37.0*Ri**2)

20    continue

! --------------------------------------------------------------------
! The eddy diffusivity of the last node is assumed to equal the
! eddy diffusivity of the second last node.
! --------------------------------------------------------------------

      de(numnod)=de(numnod-1)

      return

      end

      subroutine iceform (qnetice,T,Tcutoff,&
          fracprv,fracadd, fracice,hi,dt,numnod)

! ====================================================================
! Calculate the form of new ice in the lake as long as the fractional
! coverage of the ice is not 1 yet.
!
! Parameters :
!
! qnetice	Heat flux absorbed into the ice (W/m2).
! T		Water temperatures for all the nodes (C).
! Tcutoff	Temperature at which water freezes (C).
! fracprv	Fractional coverage of ice before this calculation (-).
! fracadd	Added fractional ice coverage (-).
! fracice	New fractional ice coverage (-).
! hice		Ice height (m).
! dt		Time step (s).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/iceform.h"

! ====================================================================
! Calculate the newly added ice to the ice layer.
! ====================================================================

      sum=0.

      do 100 j=1,numnod

         if (t(j,1).lt.Tcutoff) then

! --------------------------------------------------------------------
! Recalculate density and specifi! heat for frozen water.
! --------------------------------------------------------------------

            call density (t(j,1),dnsty)
            call specheat (t(j,1),cp)

! --------------------------------------------------------------------
! Calculate the ice growth.
! --------------------------------------------------------------------

            extra=(Tcutoff-t(j,1))*dz_lake*(dnsty+1.e3)*cp

            if (j.eq.1) then

               extra=(Tcutoff-t(j,1))*surf*(dnsty+1.e3)*cp
 
            endif

! --------------------------------------------------------------------
! If ice, the temperature of the water is the freezing temperature.
! --------------------------------------------------------------------

            t(j,1)=Tcutoff

! --------------------------------------------------------------------
! Calculate the total ice growth.
! --------------------------------------------------------------------

            sum=sum+extra

         endif

100   continue

! ====================================================================
! Calculate the heat flux absorbed into the ice and the thickness of the
! new ice.
! ====================================================================

      qnetice=(sum/dt)*(1.0-fracprv)
      if (fracprv.le.0.0) hi=fracmin
      di=sum/(fusion*rhoice)

! ====================================================================
! Calculate the added fractional coverage of ice, make sure the
! total fractional coverage does not exceed 1.
! ====================================================================

      fracadd=(di/fracmin)*(1.0-fracprv)

      if ((fracadd+fracice).gt.1.0) then

         xfrac=(fracice+fracadd)-1.0
         di=xfrac*fracmin/1.0
         hi=hi+di
         fracadd=-999.

      endif

      return

      end

      subroutine icerad ( sw, hi, hs, condbar, val, val2 )

! ====================================================================
! Calculate the radiation balance over ice.
!
! Paramterers :
!
! sw		Incoming solar radiation (W/m2).
! hi		Ice depth (m).
! hs		Snow depth (m).
! condbar	Thermal conductivity of the ice and snow pack
!		combined (W/mK).
! val1		Net short wave radiation at the top of the lake.
! val2		Incoming short wave radiation at the bottom of
!		the snow-ice layer.
! ====================================================================
 
      implicit none
      include "help/icerad.h"

! ====================================================================
! Calculate the thermal conductivity of the combined snow and ice
! layer.
! ====================================================================

      condbar=(hs*condi+hi*conds)/(condi*conds)

! ====================================================================
! Calculate the incoming radiation at different levels in the
! ice-snow layer.
! ====================================================================

! --------------------------------------------------------------------
! Calculation of constants.
! --------------------------------------------------------------------

      a=(1.-exp(-lamssw*hs))/(conds*lamssw)
      b=exp(-lamssw*hs)*(1-exp(-lamisw*hi))/(condi*lamisw)
      c=(1.-exp(-lamslw*hs))/(conds*lamslw)
      d=exp(-lamslw*hs)*(1-exp(-lamilw*hi))/(condi*lamilw)

! --------------------------------------------------------------------
! Solar radiation at bottom of snow pack.
! --------------------------------------------------------------------

      val=sw*a1*(a+b)+sw*a2*(c+d)

! --------------------------------------------------------------------
! Solar radiation at bottom of ice pack.
! --------------------------------------------------------------------

      val2=-a1*sw*(1-exp(-(lamssw*hs+lamisw*hi)))&
           -a2*sw*(1-exp(-(lamslw*hs+lamilw*hi)))

      return

      end

      subroutine initlk(lakpix,xlat_a,xlong_a,&
                        eta_a,surface,temp_a,&
                        tempi_a,hice_a,hsnow_a,preca_a,&
                        fraci_a,&
                        hice_in,hsnw_in,mixmax,numnod)

! ====================================================================
! Initialize a list of lake variables and constants.
!
! Parameters :
!
! lapkix	Lake pixel number (-).
! xlat_a	Latitude of the pixel (degrees).
! xlong_a	Longitude of the pixel (degrees).
! eta_a		Decline of solar radiation input with depth (m-1).
! tempi_a	Temperature of the ice layer (C).
! hice_a	Depth of the ice layer (m).
! hsnow_a	Snow height (m).
! preca_a	Accumulation precipitation (m).
! fraci_a	Fraction of lake covered with ice (-).
! hice_in	Initial ice height (m).
! hsnw_in	Initial snow height (m).
! mixmax	Maximal depth of convective mixing (nodes number).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/initlk.h"

! ====================================================================
! Initialize the unit number containing information about the lake
! location, depth, ...
! ====================================================================

      idata = 59

! ====================================================================
! Read in the lake locations, depth, solar radiation versus depth
! decay, surface per node number (if required).
! ====================================================================

      if (ar_flag) then

! --------------------------------------------------------------------
! Area changes with depth.
! --------------------------------------------------------------------

         read (idata,*) numnod,xlat_a,xlong_a,eta_a,&
                        (surface(lakpix,k),k=1,numnod)  

         if (numnod.gt.MAX_NOD) then

            write (*,*) 'numnod in unit ',idata,' larger than MAX_NOD'
            write (*,*) numnod,MAX_NOD
            stop

         endif

      else

! --------------------------------------------------------------------
! Area does not change with depth.
! --------------------------------------------------------------------

         read (idata,*) numnod,xlat_a,xlong_a,eta_a

         if (numnod.gt.MAX_NOD) then

            write (*,*) 'numnod in unit ',idata,' larger than MAX_NOD'
            write (*,*) numnod,MAX_NOD
            stop

         endif

         do k=1,numnod

! ....................................................................
! Initialize the surface of all the lake pixel node numbers in
! case the surface of all nodes is the same.
! ....................................................................

            surface (lakpix,k) = 1. 

         enddo

      endif         

      if (numnod.gt.MAX_NOD) then

         write (*,*) 'initlk.f (in lake_mod.f) : numnod is greater '
         write (*,*) 'than MAX_NOD ',numnod,MAX_NOD
         stop

      endif

! ====================================================================
! 2. Initialize lake temperatures.
! ====================================================================

      do 10 j=1,numnod

! --------------------------------------------------------------------
! Surface temperature is 25 C, deep temperature 10 C, the rest a
! linear interpolation.
! --------------------------------------------------------------------

         rrr=real(j-numnod)/real(1-numnod)
         temp_a(lakpix,j)=15.+10.*rrr
!CVAL for nsa bp
         temp_a(lakpix,j)=25.
!CVAL for distributed simulations over study areas
         temp_a(lakpix,j)=-15.
!CVAL for regional runs
         temp_a(lakpix,j)=5.-30.*rrr

10    continue

! ====================================================================
! 3. initialize other lake variables.
! ====================================================================
 
      preca_a = 0. 

      tempi_a = tice  

      hice_a = hice_in
      hsnow_a = hsnw_in

      fraci_a=0.
      if (hice_a .gt. (0.)) fraci_a=1.

      mixmax=0

!CVAL      write (*,*) lakpix,xlat_a,xlong_a,eta_a,tempi_a,hice_a,hsnow_a,
!CVAL                  preca_a,fraci_a,hice_in,hsnw_in,mixmax,numnod
!CVAL      write (*,*)

      return

      end


      subroutine lake(lakpix,ta,ua,qa,psurf,&
                      sw,rlwd,lnet,mixmax,rh,Qe,Qh,&
                      xlat,xlong,eta_a,dt,surface,temp_a,&
                      tempi_a,hice_a,hsnow_a,preca_a,&
                      fraci_a,precacc,numnod,rnet)

! ====================================================================
! This subroutine solves the energy budget for open water bodies.
!
! Paremeter description :
!
! lakpix	Current lake pixel number, starting at top left
!		corner of the image, going right till the end of the
!		top line, than starting again at the second top line
!		at the left.
! ta		Air temperature (K).
! ua		Wind speed (m/s).
! qa		Specifi! humidity (kg/kg).
! psurf		Air pressure (Pa).
! prec		Precipitation (m/s).
! sw		Incoming short wave radiation (W/m2).
! rlwd		Downwelling long wave radiation (W/m2).
! lnet 		Downwelling - upwelling long wave radiation (W/m2).
! mixmax	Top depth of local instability (node number).
! rh		Relative humidity (%).
! Qe		Latent heat flux (W/m2).
! Qh		Sensible heat flux (W/m2).
! xlat		Latitude of the pixel (Degrees).
! xlong		Longitude of the pixel (Degrees).
! eta_a		Decline of solar radiation input with depth (m-1).
! surface	Area of the lake pixel at different depths, node at
!		surface being node 1 (m2).
! dt		Time step (seconds).
! temp_a	Temperature of each of the nodes (C).
! tempi_a	Ice temperature (C).
! hice_a	Height of the ice layer (m).
! hsnow_a	Snow height (m).
! preca_a	Precipitation (m/s).
! precacc	Accumulated precipitation (m/s).
! fraci_a	Fraction of the lake covered with ice (-).
! numnod	Number of nodes in the lake (-).
! rnet		Net radiation (W/m2).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/lake.h"

      real*8 xfac,Qg

! ====================================================================
! 1. Initialize and read in info from previous dt.
! ====================================================================

! --------------------------------------------------------------------
! Save the incoming solar radiation in a separate variable.
! --------------------------------------------------------------------

      sworig=sw

! --------------------------------------------------------------------
! Take the lake temperature used in the rest of TOPLATS and put the into
! lake temperature used in the lake module alone.
! --------------------------------------------------------------------

      i_shuf = 1

      call shuffle (lakpix, i_shuf, t,temp_a, numnod)

! --------------------------------------------------------------------
! Calculate the water density depending on lake temperature.
! --------------------------------------------------------------------
 
      call rhoinit (Tcutoff,rhostp  )

! --------------------------------------------------------------------
! Convert air temperature from K into C.
! --------------------------------------------------------------------

      taC= ta - 273.15

! --------------------------------------------------------------------
! Initialize liquid water and ice temperature profiles to be used
! later in the module.
! --------------------------------------------------------------------

      do k=1,numnod

         t(k,2) = temp_a(lakpix,k)
         ti(k,1) = temp_a(lakpix,k)
         ti(k,2) = temp_a(lakpix,k)

      enddo

! ====================================================================
! 2. Calculate added precipitation and total snow height.
! ====================================================================

! --------------------------------------------------------------------
! Calculate the precipitation in meter and the total precipitation
! over the simulation period.
! --------------------------------------------------------------------

      prec=preca_a

      addprec = prec*dt
      precacc=precacc+addprec

! --------------------------------------------------------------------
! If the air temperature is colder than 2.2 C and there is ice on
! the lake than the precipitation fallen is treated as snow and
! added to the snow layer.
! --------------------------------------------------------------------

      if (hice_a.gt.(0.)) then

         if ( (taC.le.(2.2)) .and. (hice_a.gt.(0.)) ) then

            hsnow_a=hsnow_a+addprec*(rhowat/rhosnow)

         endif

      endif

! ====================================================================
! 3. Calculate incoming solar radiation over water and ice.
! ====================================================================

! --------------------------------------------------------------------
! Calculate the albedo of the lake depending on the solar zenith
! angle for ice, snow and liquid water.
! --------------------------------------------------------------------

      call alblake ( taC, tcutoff,&
                     albs, albi, albw )

! --------------------------------------------------------------------
! Calculate the incoming solar radiaton for both the ice fraction
! and the liquid water fraction of the lake.
! --------------------------------------------------------------------

      if (hsnow_a.gt.snocrit) then

         swi=sw*(1.-albs)

      else if ( (hsnow_a.gt.(0.)) .and. (hsnow_a.le.snocrit) ) then

         albi=(albi+albs)/2.
         swi=sw*(1.-albi)

      else if ( (hice_a.gt.(0.)) .and. (hsnow_a.le.(0.)) ) then

         swi=sw*(1.-albi)

      endif

      sww=sw*(1.-albw)

! ====================================================================
! 4. Calculate evaporation, adjust fluxes over ice.
! ====================================================================

! --------------------------------------------------------------------
! Calculate fluxes for open water.
! --------------------------------------------------------------------

! ....................................................................
! Pass the skin temperature of the lake in Kelvin since the
! Stefan-Boltzmann formula uses K.
! ....................................................................

      tin=T(1,1)+273.15
      Tcutk=Tcutoff+273.15

! ....................................................................
! Send an ice height of 0 meters to latsens for the calculation
! of the energy balance over liquid water.
! ....................................................................

      hicedum=0.0

      call latsens ( tin, Tcutk, hicedum, ta,&
                     qa, ua, psurf, rhosurf, delq, evapw, Qhw )

      Qew = -evapw*Le

! --------------------------------------------------------------------
! Adjust fluxes for ice cover.
! --------------------------------------------------------------------

      if (hice_a.gt.(0.)) then

! ....................................................................
! Pass the skin temperature of the ice in Kelvin since the
! Stefan-Boltzmann formula uses K, the initial estimate of the
! evaporation over ice is the eveporation over water, same for
! the sensible heat flux.
! ....................................................................

         tempi_a=tempi_a+273.15
         evapi=evapw
         Qhi=Qhw

         call adjflux ( swi, tempi_a , Tcutk, hice_a, hsnow_a, ta,&
                        qa, ua, psurf, rhosurf, delq, evapi, Qhi, rlwd )

         Qei = -evapi*Lei
         tempi_a=tempi_a-273.15

      endif

! ====================================================================
! 5. Calculate long wave fluxes over ice and water.
! ====================================================================

! --------------------------------------------------------------------
! Convert the ice and water temperatures into Kelvin.
! --------------------------------------------------------------------

      tkw=273.15+T(1,1)
      tki=273.15+tempi_a

! --------------------------------------------------------------------
! Calculate the outgoing long wave fluxes, positive downwards.
! --------------------------------------------------------------------

      luw= -0.97*delta*tkw**4.
      lui= -0.97*delta*tki**4.

! --------------------------------------------------------------------
! Calculate the net long wave balance.
! --------------------------------------------------------------------

      lnetw=rlwd+Luw
      lneti=rlwd+Lui

! ====================================================================
! Calculate the initial total energy balance of the lake of
! sensible and latent heat fluxes and total short wave and long wave
! radiation.
! ====================================================================

      fracprv=fraci_a
      evap=(1.-fracprv)*evapw+fracprv*evapi
      Qe=(1.-fracprv)*Qew+fracprv*Qei
      Qh=(1.-fracprv)*Qhw+fracprv*Qhi
      lu=(1.-fracprv)*luw+fracprv*lui
      lnet=(1.-fracprv)*lnetw+fracprv*lneti
      sw=(1.-fracprv)*sww+fracprv*swi

      eflux=(sw+lnet+Qe+Qh)
      eadd=eflux*dt

! ====================================================================
! 7. Calculate change in ice thickness and fraction
!    within fraction that already has ice.
! ====================================================================
 
      if (fraci_a.gt.(0.)) then

         call lakeice ( rlwd, tempi_a, Qhi, Qei, Tcutoff,&
                        swi, hice_a, hsnow_a, ti(1,1), qbot, qw, snowmlt,&
                        evapi, qnetice , fraci_a , evaps, dt)

      endif

! --------------------------------------------------------------------
! Qnetice was only used for ice formation, set again to zero.
! --------------------------------------------------------------------

      qnetice=0.

! ====================================================================
! 8. Adjust temperatures of water column in open water fraction.
! ====================================================================

      if (fracprv.lt.(1.0)) then

! --------------------------------------------------------------------
! Iwater 1 means we are calculation for water, not ice.
! --------------------------------------------------------------------

         iwater=1

! --------------------------------------------------------------------
! Calculate the eddy diffusivity.
! --------------------------------------------------------------------

         call eddy ( iwater, u2, Ti, dnsty, de, xlat, numnod )

! --------------------------------------------------------------------
! Calculate the lake temperatures at different levels for the
! new timestep.
! --------------------------------------------------------------------

         call temp_area (iwater, qbot, qw, T, lakpix,&
                         sww, lnetw, Qew, Qhw, dnsty, de, eta_a,&
                         dt, surface, numnod)

! --------------------------------------------------------------------
! Do the convective mixing of the lake water.
! --------------------------------------------------------------------

         mixdep  = 0

         call tracer_mixer (T, dnsty, lakpix,&
                            mixdep,iwater,surface,numnod)

         if (mixdep .gt. mixmax) mixmax=mixdep

      endif

! ====================================================================
! 9. Adjust temperatures of water column in ice fraction.
! ====================================================================

      if (fracprv.gt.(0.)) then

! --------------------------------------------------------------------
! Iwater 0 means we are calculating for ice.
! --------------------------------------------------------------------

         iwater=0

! --------------------------------------------------------------------
! Calculate the eddy diffusivity.
! --------------------------------------------------------------------

         call eddy (iwater, u2, T, dnsty, de, xlat, numnod )

! --------------------------------------------------------------------
! Calculate the lake temperatures at different levels for the
! new timestep.
! --------------------------------------------------------------------

         call temp_area (iwater, qbot, qw, Ti, lakpix,&
                         swi, lneti, Qei, Qhi, dnsty, de, eta_a,&
                         dt, surface, numnod)

! --------------------------------------------------------------------
! Do the convective mixing of the lake water.
! --------------------------------------------------------------------

         call tracer_mixer (T, dnsty, lakpix,&
                            mixdep, iwater,surface,numnod)

         if (mixdep .gt. mixmax) mixmax=mixdep
   
      endif

! ====================================================================
! 10. Calculate ice formation in open water fraction.
! ====================================================================

      if ( (fracprv.lt.(1.)) .and. (t(1,1).lt.Tcutoff) ) then

! --------------------------------------------------------------------
! If the fraction of the lake covered by ice is lower than 1 and
! if the water surface temperature is freezing than ice will
! form over the lake.
! --------------------------------------------------------------------

         fracadd=0.

         if (iceflag) then

            call iceform (qnetice,T,Tcutoff,fracprv,fracadd,&
                          fraci_a,hice_a,dt,numnod)

!CVAL added
            call adjflux ( swi, tempi_a , Tcutk, hice_a, hsnow_a, ta,&
                           qa, ua, psurf, rhosurf, delq, evapi, Qhi, rlwd )

            Qei = -evapi*Lei

            Qe=(1.-fraci_a)*Qew+fraci_a*Qei
            Qh=(1.-fraci_a)*Qhw+fraci_a*Qhi

!CVAL added end
 
            fraci_a=fraci_a+fracadd 
            if (fracadd.eq.-999.) fraci_a=1.0

         endif

      endif

! ====================================================================
! 11. Average ice and water columns.
! ====================================================================

      call colavg (t,ti,fracprv,dnsty,numnod)

! ====================================================================
! 12. Calculate the net radiation.
! ====================================================================

!CVAL added

      if ( (t(1,1).lt.-100.).or.(t(1,1).gt.100.) ) then

         write (*,*) 'LAKE_MAIN : check 1'
         t(1,1)=ta

      endif

      if ( (tempi_a.lt.-100.).or.(tempi_a.gt.100.) ) then

         tempi_a=0.
         write (*,*) 'LAKE_MAIN : check 2'

      endif

      if ( (Qe.lt.-2500.).or.(Qe.gt.2500.) ) then

         Qe=0.
         write (*,*) 'LAKE_MAIN : check 3'

      endif

      if ( (Qh.lt.-2500.).or.(Qh.gt.2500.) ) then

         Qh=0.
         write (*,*) 'LAKE_MAIN : check 4'

      endif

!CVALEND

      tw1=t(1,1)+273.15
      tw2=tempi_a+273.15

      r1=0.97*delta*tw1**4.
      r2=0.97*delta*tw2**4.

      rtu=fraci_a*r2+(1.-fraci_a)*r1

      rnet=rlwd-rtu+sw

! ====================================================================
! 13. Place the lake temperatures that are used only for the lake module
! back into the lake temperatures for the rest of TOPLATS.
! ====================================================================

      i_shuf = 2
      call shuffle (lakpix, i_shuf, t,temp_a, numnod)

! ====================================================================
! 14. In the lake module latent and sensible heat fluxes are positive
! downwards, to pass them back to TOPLATS the fluxes have to be
! changed sign.
! ====================================================================

      Qe=-Qe
      Qh=-Qh

      Qg=(de(1)+dm)*(temp_a(lakpix,1)-temp_a(lakpix,2))/surf

      xfac=rnet-(Qe+Qh)
!CVAL      xfac=rnet-(Qe+Qh+Qg)

!CVAL      rnet=rnet-xfac
      Qe=Qe+0.5*xfac-0.5*Qg
      Qh=Qh+0.5*xfac-0.5*Qg

! ====================================================================
! 15. Reset the incoming solar radiation to its original value.
! ====================================================================

      sw=sworig

      return

      end

      subroutine lakeice ( radlwd, tempice, qsen, qlat,&
          tcutoff, sw, hi, hs, twater,  qbot, qw, ds, evapi,&
          qnetice, fracice , evaps, dt )

! ====================================================================
! Calculate the growth and decrease in the lake ice cover.
!
! Parameters :
!
! radlwd	Downwelling long wave radiation (W/m2).
! tempice	Temperature of the ice (K).
! qsen		Sensible heat flux (W/m2).
! qlat		Latent heat flux (W/m2).
! tcutoff	Temperature at which water freezes (K).
! sw		Incoming short wave radiation (W/m2).
! hi		Ice height (m).
! hs		Snow height (m).
! twater	Temperature of the top water layer (K).
! qbot		Incoming short wave radiation at bottom of the ice (W/m2).
! qw		Heat storage in the lake ice (J/m3).
! ds		Amount of melt of snow pack (m).
! evapi		Evaporation from the ice (m/s).
! qnetice	Heat flux absorbed into the ice (W/m2).
! fracice	Fraction of the lake covered by ice (-).
! evaps		Evaporation from the snow (m).
! dt		Time step (s).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/lakeice.h"

! ====================================================================
! Initialize some functions that will be used later in the program.
! ====================================================================

      t4(x)=(x+273.15)**4.
      q0t0(x)=sw+(1./condbar)*(tcutoff-(x)-val)

! ====================================================================
! Calculate the thermal conductivity of water.
! ====================================================================

      condqw = rhowat * cpw_ice * surf / qwtau

! ====================================================================
! Convert evaporation from mm/s to m, keep track of original incoming
! variables, calculate the net radiation for the ice layer.
! ====================================================================

      evapl=evapi*dt/1000.
      tprev=tempice
      qmet=radlwd-emice*delta*t4(tprev)+qsen+qlat

! ====================================================================
! Calculate the radiation input for the top of the ice, calculate the
! temperature of the ice.
! ====================================================================

      call icerad (sw, hi, hs, condbar, val, val2 )
      q0=-qmet
      tempice=condbar*(sw-q0)+tcutoff-val 
      qbot=sw+val2

      if (tempice.gt.tmelt) then

! --------------------------------------------------------------------
! If the ice temperature is above the melt temperature, adjust the
! ice temperature and calculate the heat flux used to melt the ice.
! --------------------------------------------------------------------

         q0=q0t0(tmelt)
         tempice=tmelt
         qmelts=q0+qmet

      else

! --------------------------------------------------------------------
! If the ice temperature is lower than the melt temperature there is
! no energy used to melt the ice.
! --------------------------------------------------------------------

         qmelts=0.0

      endif

! ====================================================================
! Calculate fluxes at the base of the ice.
! ====================================================================

! --------------------------------------------------------------------
! Flux of heat into ice from base of ice.
! --------------------------------------------------------------------

      qf=q0+val2

! --------------------------------------------------------------------
! Amount of heat needed to melt the ice.
! --------------------------------------------------------------------

      qw=-condqw*(tcutoff-twater)

! --------------------------------------------------------------------
! Amount of heat used to melt the ice.
! --------------------------------------------------------------------

      qmeltb=qf-qw

! ====================================================================
! Total of amount of heat used to melt and freeze the ice from bottom
! and top of the ice.
! ====================================================================

      qnetice=qmeltb-qmelts

! ====================================================================
! If qnetice > 0 there is freezing occuring, the amount of heat
! released from this has to be included in the energy balance.
! ====================================================================

! --------------------------------------------------------------------
! Initialize the amount of snow melt.
! --------------------------------------------------------------------

      qmeltsx=0.0

! --------------------------------------------------------------------
! Calculate the evaporation from the snow pack and adjust the snow
! height.
! --------------------------------------------------------------------

      if (hs.gt.(0.)) then

! ....................................................................
! In presence of a snow pack adjust the snow height for evaporation
! and condensation.
! ....................................................................

         if ( (evapl*(rhowat/rhosnow)) .le.hs) then

! ....................................................................
! Condensation adds to the snow layer, evaporation removes from the
! snow layer (if the snow layer is large enough) and all the ice
! evaporation is used to evaporate from the snow.
! ....................................................................

            hs=hs-evapl*(rhowat/rhosnow)
            evapl=0.0
            evaps=evapi

         else

! ....................................................................
! All the snow layer is evaporated (if the evporative demand is
! larger than the snow height), the remaining evaporation comes from
! the ice.
! ....................................................................

            evapl=evapl-hs*(rhosnow/rhowat)
            evaps=hs*(rhosnow/rhowat)
            hs=0.0

         endif

      endif

! ====================================================================
! Calculate the amount of snow melt.
! ====================================================================

      if (hs.gt.0.0) then

! --------------------------------------------------------------------
! Calculate the growth of the snow pack, negative : melting.
! --------------------------------------------------------------------

         ds=(-qmelts/(rhosnow*fusion))*dt

         if (-ds.gt.hs) then

! ....................................................................
! If the amount of melt is larger than the snow pack ice will be
! melted as well.  Calculate the melt energy remain after melting
! the snow pack and put the snow height to zero.
! ....................................................................

            qmeltsx=qmelts-(hs*rhosnow*fusion/dt)
            hs=0.0 

         else

! ....................................................................
! If the snow pack grows calculate its new height.
! ....................................................................

            hs=hs+ds

         endif

      endif

! ====================================================================
! Calculate the growth of the ice pack.
! ====================================================================

      if (hs.le.(0.)) then

! --------------------------------------------------------------------
! If now snow present : calculate the growth of the ice pack at the
! top.  First term : energy from the ice layer itself, second term : extra
! energy from the snow pack, third term : evaporation from the ice layer.
! --------------------------------------------------------------------

         disurf=((-qmelts/(rhoice*fusion))+&
                (-qmeltsx/(rhoice*fusion)))*dt+&
                (-evapl*(rhowat/rhoice))

      else

! --------------------------------------------------------------------
! If there is snow the ice layer does not change at the top.
! --------------------------------------------------------------------

         disurf=0.0

      endif

! --------------------------------------------------------------------
! Calculate the growth of the ice pack at the bottom.
! --------------------------------------------------------------------

      dibot=(qmeltb/(rhoice*fusion))*dt

! ====================================================================
! Calculate the new height of the ice pack and the fractional ice
! cover.
! ====================================================================

      if (fracice.ge.(1.)) then

! --------------------------------------------------------------------
! If the lake is fully covered with ice change only the thickness and
! not the fractional cover.
! --------------------------------------------------------------------

         hiprv=hi
         hi=hi+disurf+dibot

         if (hi.lt.fracmin) then

! ....................................................................
! If the ice height is lower than the minimum ice height increase the
! height and decrease the fractional cover in order to conserve
! numerical stability.
! ....................................................................

            extradi=fracmin-hi
            df=(extradi/fracmin)*fracice
            fracice=1.0-df
            hi=fracmin

         endif

      else

! --------------------------------------------------------------------
! If the lake is not fully covered with ice change both the fractional
! cover and the ice height.
! --------------------------------------------------------------------

         df=fracice*(disurf+dibot)/fracmin
         fracice=fracice+df

         if (fracice.gt.1.0) then

! ....................................................................
! If the fractioncal cover is larger than 100 % decrease the
! fractional cover to 100 % and increase the ice depth.
! ....................................................................

            xfrac=fracice-1.0
            di=xfrac*fracmin/1.0
            hi=hi+di
            fracice=1.0

         endif

         if ( (fracice.lt.fraclim) .and. (df.lt.(0.)) ) then

! ....................................................................
! If the ice thickness is lower than the minimal ice thickness
! and melting has occured set the ice and snow layer thicknesses to zero.
! ....................................................................

            xfrac=fracice
            diextra=xfrac*fracmin/1.0
            extraf=-diextra*rhoice*fusion*(1./dt)

            qw=qw-extraf
            qnetice=qnetice+extraf
            fracice=0.0
            hi=0.0
            hs=0.0

         endif

      endif

      return

      end

      subroutine latsens ( tsurf, Tcutoff, hice, ta,&
           qa, ua, psurf, rhosurf, delq, evap, qsen )

! ====================================================================
! Calculate the partitioning of the energy balance into latent and
! sensible heat fluxes.
!
! Parameters :
!
! tsurf		Lake surface temperature (C).
! Tcutoff	Temperature at which water freezes (K).
! hice		Ice height (m).
! ta		Air temperature (K).
! qa		Absolute air humididity (-).
! ua		Wind speed (m/s).
! psurf		Air pressure (Pa).
! rhosurf	Air density (kg/m3).
! delq		Difference in absolute humidity between the lake
!		surface and higher up (-).
! evap		Evaporation rate (m/s).
! qsen		Sensible heat flux (W/m2).
! eog		Vapor pressure at lake level (Pa).
! ====================================================================

      implicit none
      include "help/latsens.h"

! ====================================================================
! Calculate the drag coefficient.
! ====================================================================

      call lkdrag (tsurf, ta, ua, cdrx )

! ====================================================================
! Determine the coefficients to be used in the calculation of
! the vapor pressure depending on whether the lake is covered with
! ice or not.
! ====================================================================

      if ( (hice.le.(0.)) .and. (tsurf.gt.Tcutoff) ) then

! --------------------------------------------------------------------
! Lake is not covered with ice.
! --------------------------------------------------------------------

         a=c72
         b=c73

      else

! --------------------------------------------------------------------
! Lake is covered with ice.
! --------------------------------------------------------------------

         a=c70
         b=c71

      endif

! ====================================================================
! Calculate the vapor pressure at lake level.
! ====================================================================

      eog=ca*exp(a*(tsurf-cb)/(tsurf-b))

! ====================================================================
! Calculate the absolute air humidity at lake level and the difference
! between this humidity and the air humidity at measurement height.
! ====================================================================

      qg=0.622*(eog/(psurf-0.378*eog))
      delq=qa-qg

! ====================================================================
! Calculate the evaporation rate.
! ====================================================================

      rai=cdrx*ua*rhosurf
      evap=-rai*delq

! ====================================================================
! Calculate the difference in lake surface temperature and the
! air temperature at measurement height and consequently the
! sensible heat flux.
! ====================================================================

      delt=ta-tsurf
      rash=cdrx*ua*rhosurf*cpair
      qsen=rash*delt

      return

      end


      subroutine lkdrag (tlakek, t1k, u1, cdrx )

! ====================================================================
! Calculate the lake drag coefficient.
!
! Parameter :
!
! tlakek 	Lake surface temperature (K).
! t1k		Air temperature (K).
! u1		Wind speed (m/s).
! cdrx		Drag coefficient.
! ====================================================================

      implicit none
      include "help/lkdrag.h"

! ====================================================================
! Calculate the Richardson number.
! ====================================================================

      cdrn=(vk/log(z1/zwater))**2.

      ribn=z1*g*(1.-tlakek/t1k)

      if ((tlakek/t1k).le.1.0) then

         ribd=u1**2.+0.1**2.

      else

         ribd=u1**2.+1.0**2.

      endif

      rib=ribn/ribd

! ====================================================================
! Calculate the drag coefficient using the Richardson number.
! ====================================================================

      if (rib.lt.0.) then

         cdr=cdrn*(1.0+24.5*(-cdrn*rib)**0.5)

      else

         cdr=cdrn/(1.0+11.5*rib)

      endif

      cdrmin=.25*cdrn
      if (cdrmin.lt.6.e-4) cdrmin=6.e-4

      if (cdr.lt.cdrmin) cdr=cdrmin

      cdrx=cdr

      return

      end


      subroutine readt ( tp, hice, hsnow, temp_a,&
                          tempi_a, hice_a, hsnow_a, lakpix)

! ====================================================================
! Read in parameters from the previous time step.
!
! tp		Lake temperature (C).
! hice		Ice height (m).
! hsnow		Snow height (m).
! temp_a	Lake temperature (C).
! tempi_a	Ice temperature (C).
! hice_a	Ice height (m).
! hsnow_a	Snow height (m).
! lakpix	Lake pixel number (-).
! ====================================================================
   
      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/readt.h"

      if (hice.le.(0.)) then

         tp= temp_a (lakpix,1) + 273.15

      else

         tp= tempi_a + 273.15

      endif

      hice = hice_a
      hsnow = hsnow_a

      return

      end


      subroutine rhoinit (tfsp,rhostp)

! ====================================================================
! Calculate the temperature at which water freezes depending on
! salinity and air pressure.
!
! Paramters :
!
! tfsp		Lake freezing point (C).
! rhostp	Water density (kg/m3).
! ====================================================================

      implicit none
      include "help/rhoinit.h"

! ====================================================================
! Air pressure is assumed to be 1.013 Pa, salinity is assumed to
! be zero, initial water temperature is assumed to be 10 C.
! ====================================================================

      s=0.
      p=1.013
      s=0.
      t=10.

! ====================================================================
! Calculate the lake freezing temperature.
! ====================================================================

      tfsp=-0.0575*s + 1.710523e-3*s**(3./2.) - 2.154996e-4*s**2. -&
           7.53e-3*p

! ====================================================================
! Calculate the lake water density.
! ====================================================================

      rhow=999.842594 + 6.793952e-2*t - 9.095290e-3*t**2.+&
           1.001685e-4*t**3. - 1.120083e-6*t**4. + 6.536332e-9*t**5.

      sbmw=19652.21 + 148.4206*t - 2.327105*t**2. + 1.360477e-2*t**3. -&
           5.155288e-5*t**4.

      rhostp=rhow/(1.-p/sbmw)

      return

      end

      subroutine shuffle (lakpix, i_shuf, t_shuf,temp_a, numnod)

! ====================================================================
! This subroutine switches lake temperatures used by TOPLATS with
! lake temperatures used by the lake module and vice versa.
!
! Parameters :
!
! lakpix		Lake pixel number.
! i_shuf		1 : Swictch variables used by TOPLATS to variables
!			    used by the lake module.
!			2: Swictch variables used by the lake module to
!			   variables used by TOPLATS.c
! t_shuf, temp_a	Lake water temperatures (C).
! numnod		Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/shuffle.h"

      if (i_shuf.eq.1) then

! ====================================================================
! Put variables used by TOPLATS into variables used by the lake
! module.
! ====================================================================

         do k = 1, numnod 

            t_shuf(k,1) = temp_a (lakpix,k) 

         enddo 

      else if (i_shuf.eq.2) then

! ====================================================================
! Put variables used by the lake module back into variables used by
! TOPLATS.
! ====================================================================

         do k = 1, numnod 

!CVAL t_shuf(k,2) rather than t_shuf(k,1)
            temp_a (lakpix,k) = t_shuf(k,2)  

         enddo 

      endif 

      return 

      end 
 
      subroutine specheat (t,cpt)

! ====================================================================
! Calculate the specifi! heat of the water depending on water
! temperature.  Salinity is assumed to be zero.
!
! Paramterers :
!
! t	Water temperature (C).
! cpt	Specifi! heat (J/Kg K).
! ====================================================================

      implicit none
      include "help/specheat.h"

      cpt=4217.4 - 3.720283*t + 0.1412855*t**2. - 2.654387e-3*t**3. +&
           2.093236e-5*t**4.

      return

      end

      subroutine temp_area( iwater, qbot, qw,&
                     T, lakpix, sw, lnet, Qe, Qh, dnsty, de,&
                     eta, dt, surface, numnod)

! ====================================================================
! Calculate the water temperature for different levels in the lake.
!
! Parameters :
!
! iwater	0 : calculation for ice.
!		1 : calculation for liquid water.
! qbot		Incoming short wave radiation at bottom of the ice (W/m2).
! qw		Heat storage in the lake ice (J/m3).
! T		Lake water temperature at different levels (K).
! lakpix	Number of lake pixel.
! sw		Incoming short wave radiation (W/m2).
! lnet		Net long wave radiation (W/m2).
! Qe		Latent heat flux (W/m2).
! Qh		Sensible heat flux (W/m2).
! dnsty		Water density at different levels (kg/m3).
! de		Diffusivity of water (or ice) (m2/d).
! eta		Decline of solar radiation input with depth (m-1).
! dt		Time step size (s).
! surface	Area of the lake at different levels (m2).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/temp_area.h"

! ====================================================================
! Calculate the distance between the centers of the surface and first
! lake layers.
! ====================================================================

      dist12=surf*0.5+dz_lake*0.5

! ====================================================================
! Initialize the water density at all the nodes and the depth of all
! and distance between all nodes.
! ====================================================================

      do 10 k= 1,numnod

         T(k,2)=T(k,1)
         told(k)=T(k,1)
         call density(T(k,1),dnsty(k))
         call specheat (T(k,1),cpz(k))
         z(k)=dz_lake
         z(k)=surf+(k-1)*dz_lake
         zhalf(k)=dz_lake

10    continue

      z(1)=surf
      zhalf(1)=0.5*(z(1)+z(2))

! ====================================================================
! Calculate the right hand side vector in the tridiagonal matrix system
! of equations.
! ====================================================================

! --------------------------------------------------------------------
! First calculate d for the surface layer of the lake.
! --------------------------------------------------------------------

      k  = 1 

      surface_1 =( surface(lakpix,k) + surface(lakpix,k)) / 2.
      surface_2 =( surface(lakpix,k) + surface(lakpix,k+1)) / 2.

      if (iwater.eq.1) then

         T1 = sw*beta +&
                (1.-beta)*&
                sw*(1.-exp(-eta*surf))*surface_2/surface(lakpix,k) +&
                (Lnet+Qe+Qh) * surface_1/surface(lakpix,k)

      else

         T1 = qbot*beta +&       
                (1.-beta)*&
                qbot*(1.-exp(-eta*surf))*surface_2/surface(lakpix,k) -&
                qw * surface_1/surface(lakpix,k) 

      endif

      cnextra = 0.5 * surface_2/surface(lakpix,k) *&
                ((de(1) / zhalf(1)) * (dt / z(1)) * ( T(2,1) - T(1,1) ))

      d(1) = T(1,1)+T1*dt/((1.e3+dnsty(1))*cpz(1)*z(1))+cnextra

! --------------------------------------------------------------------
! Calculate d for the remainder of the column.
! --------------------------------------------------------------------

! ....................................................................
! All nodes bot the deepest node.
! ....................................................................

      swtop = qbot
      if (iwater.eq.1) swtop = sw

      do 11 k=2,numnod-1

         top = (surf+(k-2)*dz_lake)
         bot = (surf+(k-1)*dz_lake)


         surface_1 =( surface(lakpix,k-1) + surface(lakpix,k)) / 2. 
         surface_2 =( surface(lakpix,k)  + surface(lakpix,k+1)) / 2. 

         T1 = (1.-beta)*swtop*&
              ( (surface_1*exp(-eta*top)-surface_2*exp(-eta*bot))/&
               surface(lakpix,k))

         cnextra=0.5 *1./surface(lakpix,k)*&
                 (((de(k)/zhalf(k))*(dt/z(k))*(T(k+1,1)-T(k,1)))*&
                 surface_2- &
                 ((de(k-1)/zhalf(k-1))*(dt/z(k))*(T(k,1)-T(k-1,1)))*&
                 surface_1)

         d(k) = T(k,1)+T1*dt/((1.e3+dnsty(k))*cpz(k)*z(k))+cnextra

11    continue

! ....................................................................
! Calculation for the deepest node.
! ....................................................................

      surface_1 =( surface(lakpix,k-1) + surface(lakpix,k)) / 2.

      swtop = qbot
      if (iwater.eq.1) swtop = sw
        k = numnod
      top = (surf+(k-2)*dz_lake)

      T1 = (1.-beta)*swtop*(exp(-eta*top))*surface_1/surface(lakpix,k)

      cnextra=0.5 * surface_1/surface(lakpix,k) *&
              ((de(k-1)/zhalf(k-1))*(dt/z(k))*(T(k,1)-T(k-1,1)))

      d(k) = T(k,1)+T1*dt/((1.e3+dnsty(k))*cpz(k)*z(k))+cnextra


! ====================================================================
! Calculate arrays for tridiagonal matrix.
! ====================================================================

! --------------------------------------------------------------------
! Top node of the column.
! --------------------------------------------------------------------

      k = 1
      surface_2 =( surface(lakpix,k)  + surface(lakpix,k+1)) / 2.

      b(1) = -0.5 * ( de(1) / zhalf(1) ) * &
                ( dt / z(1) ) * surface_2/surface(lakpix,k)
      a(1) = 1. - b(1)

! --------------------------------------------------------------------
! Second to second last node of the column.
! --------------------------------------------------------------------

      do 13 k = 2,numnod-1

         surface_1 =( surface(lakpix,k-1) + surface(lakpix,k)) / 2.
         surface_2 =( surface(lakpix,k)  + surface(lakpix,k+1)) / 2.

         b(k) = -0.5 * ( de(k) / zhalf(k) ) *&
              (  dt / z(k) )*surface_2/surface(lakpix,k)
         c(k) = -0.5 * ( de(k-1) / zhalf(k-1) ) *&
              (  dt / z(k) )*surface_1/surface(lakpix,k)
         a(k) = 1. - b(k) - c(k)

13    continue

! --------------------------------------------------------------------
! Deepest node of the column.
! --------------------------------------------------------------------

      surface_1 =( surface(lakpix,k-1) + surface(lakpix,k)) / 2.

      c(numnod) = -0.5 * ( de(numnod) / zhalf(numnod) ) *&
                 (  dt / z(numnod) ) * surface_1/surface(lakpix,k)
      a(numnod) = 1. - c(numnod)

! ====================================================================
! Solve the tridiagonal matrix.
! ====================================================================

      call tridia(1,1,numnod,c,a,b,d,tnew,work1,work2)

! ====================================================================
! Put the solutions of the tridiagonal system into the temperature
! matrix of the present time step.  Recalculate the water
! densities.
! ====================================================================

      do 40 k = 1, numnod

         t(k,2) = tnew(k)
         t(k,1) = t(k,2)
         call density (T(k,2),dnsty(k))

40    continue

      return

      end

      subroutine tracer_mixer ( T, dnsty, lakpix,&
                       mixdep,iwater,surface,numnod)

! ====================================================================
! Simulate the convective mixing in the lake.
!
! Paramters :
!
! T		Water temperatures (K).
! dnsty		Water densities (kg/m3).
! lakpix	Lake pixel number (-).
! mixdep	Top depth of local instability (node number, -).
! iwater	0 for ice, 1 for liquid water.
! surface	Area of the lake per node number (m2).
! numnod	Number of nodes in the lake (-).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "LAKE.h"
      include "help/tracer_mixer.h"

! ====================================================================
! Initialize the top depth of local instability.
! ====================================================================

      mixprev = 1

9     continue

      do 10 k=1,numnod-1

! ====================================================================
! Check for instability at each slice in water column.
! ====================================================================

         avet=0.0
         avev=0.0

         if (dnsty(k).gt.dnsty(k+1)) then

! --------------------------------------------------------------------
! If there is instability apply the mixing scheme.
! --------------------------------------------------------------------

            if (mixprev.eq.1.and.(k+1).gt.mixdep) then

! ....................................................................
! Correct the top depth of local instability.
! ....................................................................

                mixdep=k+1

            endif

!-----------------------------------------------------------------------
! Mix from mixprev to k+1  
!------------------------------------------------------------------------

            do 20 m = mixprev, k+1

! --------------------------------------------------------------------
! Apply the mixing scheme from the previous depth to the instability
! up to this node.
! --------------------------------------------------------------------

! ....................................................................
! Calculate the specifi! heat of the lake water.
! ....................................................................

               call specheat (T(m,2),cp)

               if (m .eq. 1) then

! ....................................................................
! Calculate the heat content and volume of the surface layer.
! ....................................................................

                  vol = surf*(1.e3+dnsty(m))*cp*surface(lakpix,m)
                  vol_tr = surf * surface(lakpix,m)

               else

! ....................................................................
! Calculate the heat content and volume of all layers but the surface
! layer.
! ....................................................................

                  vol = dz_lake*(1.e3+dnsty(m))*cp*surface(lakpix,m)
                  vol_tr = dz_lake * surface(lakpix,m)

               endif

! ....................................................................
! Calculate the volumetri! weighed average lake temperature.
! ....................................................................

               avet=avet+T(m,2)*vol
               avev=avev+vol

20          continue

            Tav=avet/avev

! --------------------------------------------------------------------
! Calculate the density of the surface layer.
! --------------------------------------------------------------------

            call density (Tav,densnew)

! --------------------------------------------------------------------
! Calculate the maximum density up to the local level.
! --------------------------------------------------------------------

            rho_max = 0.e3

            do 23 kk=1,mixprev-1

               if ( (1.e3+dnsty(kk)) .gt. rho_max ) then

                  rho_max=1.e3+dnsty(kk)

               endif

23          continue

! --------------------------------------------------------------------
! Adjust temperatures and density in the mixed part of column.
! --------------------------------------------------------------------

            do 30 k2 = mixprev, k+1

               T(k2,2)=Tav 
               dnsty(k2)=densnew

30          continue

! --------------------------------------------------------------------
! Check to make sure that the mixing has not generated new instabilities
! above the previous depth to the instabilities.
! --------------------------------------------------------------------

            if (rho_max.gt.(1.e3+densnew)) then

! ....................................................................
! If there are still instabilities iterate again.
! ....................................................................

               mixprev = 1
               go to 9

            endif
 
         else

! ====================================================================
! If there are no instabilities up to now then the depth to the
! instability has to be increased by 1 node.
! ====================================================================

            mixprev=k+1

         endif

10    continue

! ====================================================================
! Recalculate the water density.
! ====================================================================

      do 40 k = 1, numnod

         T(k,1)=T(k,2) 
         call density (T(k,1),dnsty(k))

40    continue

      return

      end

      subroutine tridia (ns, nd, ne, a, b, c, y, x, alpha, gamma)

! ====================================================================
! Solve a tridiagonal system of equations.
!
! Parameters :
!
! ns		The number of systems to be solved.
! nd		First dimension of arrays (larger than or equal to ns).
! ne		The number of unknowns in each system. This must 
!		be larger than or equal to 2.
! a		The sub diagonals of the matrices are stored in locations
!		a(j,2) through a(j,ne).
! b		The main diagonals of the matrices are stored in
!		locations b(j,1) through b(j,ne).
! c		The super diagonals of the matrices are stored in
!		locations c(j,1) through c(j,ne-1).
! y		The right hand side of the equations is stored in
!		y(j,1) through y(j,ne).
! x		The solutions of the systems are returned in
!		locations x(j,1) through x(j,ne).
! alpha		Work array dimensioned alpha(nd,ne).
! gamma		Work array dimensioned gamma(nd,ne).
!
! History : Based on a streamlined version of the old NCAR ULIB subroutine
! TRDI used in the PHOENIX climate model of Schneider and Thompson
! (J.G.R., 1981). Revised by Starley Thompson to solve multiple systems
! and vectorize well on the CRAY-1. Later revised to include a PARAMETER
! statement to define loop limits and thus enable Cray short vector
! loops.
!
! Algorithm:  LU decomposition followed by solution.  NOTE: This
! subroutine executes satisfactorily if the input matrix is diagonally
! dominant and non-singular.  The diagonal elements are used to pivot, and
! no tests are made to determine singularity. If a singular or numerically
! singular matrix is used as input a divide by zero or floating point
! overflow will result.
! ====================================================================

      implicit none
      include "help/tridia.h"

      nm1 = ne-1

! ====================================================================
! Obtain the LU decompositions.
! ====================================================================

      do 10 j=1,ns

         alpha(j,1) = 1./b(j,1)
         gamma(j,1) = c(j,1)*alpha(j,1)

10    continue

      do 11 i=2,nm1

         do 12 j=1,ns

            alpha(j,i) = 1./(b(j,i)-a(j,i)*gamma(j,i-1))
            gamma(j,i) = C(j,i)*alpha(j,i)

12       continue

11    continue

! ====================================================================
! Solve the system.
! ====================================================================

      do 20 j=1,ns

         x(j,1) = y(j,1)*alpha(j,1)

20    continue

      do 21 i=2,nm1

         do 22 j=1,ns

            x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)

22       continue

21    continue

      do 23 j=1,ns

         x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1))/&
                   (b(j,ne)-a(j,ne)*gamma(j,nm1))

23    continue

      do 24 i=1,nm1

         ib = ne-i

         do 25 j=1,ns

            x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)

25       continue

24    continue

      return

      end
  
      subroutine zerolake(tp, hice, hsnw)

! ====================================================================
! initialize the ice and snow heights and the lake surface temperature.
!
! Parameters :
!
! tp		Lake surface temperature (C).
! hice		ice height (m).
! hsnow		Snow height (m).
! ====================================================================

      implicit none
      include "LAKE.h"
      include "help/zerolake.h"

      tp=0.0
      hice=0.0
      hsnw=0.0

      return

      end

      subroutine convert (tdry,twet,pres,rld,rsd,pptms,uza,&
                          avpp,tap,rhp,psurfp,qap,rlwdp,swd,preca,ua)

! ====================================================================
! This subroutine converts some meteorological variables with the
! units used in TOPLATS into the units used by the lake module.
!
! Parameters : 
!
! tdry          Air temperature (C).
! twet          Dew point temperature (C).
! pres          Air pressure (mbar).
! avpp          Vapor pressure (Pa).
! tap           Air temperature (K).
! rhp           Relative humidity (-).
! psurfp        Air pressure (Pa).
! qap           Absolute humidity (-).
! rld		Down welling long wave radiation (W/m2).
! rlwdp		Down welling long wave radiation (W/m2).
! rsd		Incoming short wave radiation (W/m2).
! swd		Incoming short wave radiation (W/m2).
! pptms		Precipitation rate (m/s).
! preca		Precipitation rate (m/s).
! uza		Wind speed (m/s).
! ua		Wind speed (m/s).
! ====================================================================

      implicit none
      include "help/convert.h"

! ====================================================================
! Calculate the vapor pressure and relative humididty.
! ====================================================================

      avpp=610.8*exp( (17.27*twet)/(237.3+twet) )
      es=610.8*exp( (17.27*tdry)/(237.3+tdry) )
      rhp=avpp/es

! ====================================================================
! Convert air temperature from Centigrades into Kelvin.
! ====================================================================

      tap=tdry+273.15

! ====================================================================
! Convert air pressure from mbar into Pa.
! ====================================================================

      psurfp=pres*100.

! ====================================================================
! Calculate the absolute air humdity.
! ====================================================================

      qap=0.622*( avpp/(psurfp-0.378*avpp) )

! ====================================================================
! Initialize the long and short wave radiation.
! ====================================================================

      swd=rsd
      rlwdp=rld

! ====================================================================
! Initialize the precipitation and wind speed.
! ====================================================================

      preca=pptms
      ua=uza

      return

      end

      subroutine lakevars(evt,pn,pt,xs,xi,rt,ta,row,xle)

! ====================================================================
! Calculate lake variables for the lake pixel to be used in
! calculation of regional averages and totals in TOPLATS.
!
! Parameter definitions :
!
! evt		Evaporation rate (m/s).
! pn		Net precipitation (m/s).
! pt		Total precipitation (m/s).
! xs		Saturation excess runoff (m/s).
! xi		Infiltration excess runoff (m/s).
! rt		Total runoff (m/s).
! ta		Air temperature (C).
! row		Water density (kg/m3).
! rlv		Latent heat of vaporization (J/kg).
! xle		Latent heat flux (W/m2).
! ====================================================================

      implicit none
      include "help/lakevars.h"

! --------------------------------------------------------------------
! Calculate the latent heat of vaporization of water.
! --------------------------------------------------------------------

      rlv=2501000.-2361.*ta

! --------------------------------------------------------------------
! Calculate the evaporation rate.
! --------------------------------------------------------------------

      evt=xle/(rlv*row)

! --------------------------------------------------------------------
! Calculate the net precipitation.
! --------------------------------------------------------------------

      pn=pt

! --------------------------------------------------------------------
! Calculate the runoff and its components.
! --------------------------------------------------------------------

      xs=0.
      xi=0.
      rt=0.

      return

      end

      subroutine catchlake(ncatch,r_lakearea,npix,ivgtyp,ilandc,&
                           pixsiz,icatch)

! ====================================================================
! Calculate, for each catchment, the area of the lakes in the
! catchment.
!
! Parameter definition :
!
! ncatch	Number of catchments.
! r_lakearea	Area of lakes in each catchment (m2).
! ilandc	Land cover of each pixel.
! pixsiz	Pixel resolution (m).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/catchlake.h"

      do i=1,ncatch

         r_lakearea(i)=0.d0

      enddo

      do i=1,npix

         if (ivgtyp(ilandc(i)).eq.(-1)) then

            r_lakearea(icatch(i))=r_lakearea(icatch(i))+pixsiz*pixsiz

         endif

      enddo

      return

      end

      subroutine catchlake_s(ncatch,r_lakearea,npix,pixsiz,f_lake,icatch)

! ====================================================================
! Calculate, for each catchment, the area of the lakes in the
! catchment.
!
! Parameter definition :
!
! ncatch        Number of catchments.
! r_lakearea    Area of lakes in each catchment (m2).
! iland!        Land cover of each pixel.
! pixsiz        Pixel resolution (m).
! ====================================================================

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/catchlake_s.h"

      do i=1,ncatch

         r_lakearea(i)=0.d0

      enddo

      do i=1,npix

         r_lakearea(icatch(i))=r_lakearea(icatch(i))+f_lake(i)*pixsiz*pixsiz

      enddo

      return

      end
