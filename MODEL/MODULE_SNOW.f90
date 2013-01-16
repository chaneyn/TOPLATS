MODULE MODULE_SNOW

  USE MODULE_VARIABLES

  !implicit none

  contains

!**********************************************************************
!* SUMMARY:      SnowMelt.! - Calculate snow accumulation and melt
!* USAGE:        
!*
!* AUTHOR:       Mark Wigmosta and Pascal Storck
!* ORG:          University of Washington, Department of Civil Engineering
!* E-MAIL:       nijssen@u.washington.edu
!* ORIG-DATE:    8-Oct-1996 at 08:50:06
!* LAST-MOD:     Mon Apr 14 17:39:40 1997 by Bart Nijssen
!*               <nijssen@meter.ce.washington.edu>
!* DESCRIPTION:  Calculate snow accumulation and melt using an energy balance
!*               approach for a two layer snow model
!* DESCRIP-END.
!* FUNCTIONS:    SnowMelt()
!* COMMENTS:     This subroutine has been translated into fortran
!*               by Valentijn Pauwels, Princeton University, May 1997.
!*               Contact : vpauwels@earth.Princeton.EDU
!**********************************************************************

!*****************************************************************************
! Function name: SnowMelt()

! Purpose      : Calculate snow accumulation and melt using an energy balance
!                approach for a two layer snow model

! Required     :
!   integer y               - Row counter
!   integer x               - Column counter
!   real*8 Tstep            - Model timestep (hours)
!   real*8 Z                - Reference height (m) 
!   real*8 Displacement     - Displacement height (m)
!   real*8 Z0               - Surface roughness (m)
!   real*8 BaseRa           - Aerodynami! resistance (uncorrected for
!                              stability) (s/m)
!   real*8 AirDens          - Density of air (kg/m3)
!   real*8 EactAir          - Actual vapor pressure of air (Pa) 
!   real*8 Lv               - Latent heat of vaporization (J/kg)
!   real*8 ShortRad         - Net exchange of shortwave radiation (W/m2)
!   real*8 LongRadIn        - Incoming long wave radiation (W/m2)
!   real*8 Press            - Air pressure (Pa)
!   real*8 RainFall         - Amount of rain (m)
!   real*8 SnowFall         - Amount of snow (m)
!   real*8 Tair             - Air temperature (C)
!   real*8 Vpd              - Vapor pressure deficit (Pa)
!   real*8 Wind             - Wind speed (m/s)
!   real*8 PackWater        - Liquid water content of snow pack 
!   real*8 SurfWater        - Liquid water content of surface layer 
!   real*8 Swq              - Snow water equivalent at current pixel (m)
!   real*8 VaporMassFlux    - Mass flux of water vapor to or from the
!                             intercepted snow (m)
!   real*8 TPack            - Temperature of snow pack (C)
!   real*8 TSurf            - Temperature of snow pack surface layer (C)
!   real*8 MeltEnergy       - Energy used for melting and heating of snow pack
!                              (W/m2)
!   real*8 Density	    - Density of the snow pack (g/cm3)

! Returns      :
!   real*8 Outflow          - Amount of snowpack outflow (m)
!   real*8 LatentHeat	    - Latent heat flux (W/m2)
!   real*8 SensibleHeat	    - Sensible heat flux (W/m2)
!   real*8 NetRad	    - Net Radiation (W/m2)
!   real*8 Emiss	    - Emissivity (-)

! Modifies     :
!   real*8 PackWater        - Liquid water content of snow pack 
!   real*8 SurfWater        - Liquid water content of surface layer 
!   real*8 Swq              - Snow water equivalent at current pixel (m)
!   real*8 VaporMassFlux    - Mass flux of water vapor to or from the
!                             intercepted snow (m)
!   real*8 TPack            - Temperature of snow pack (C)
!   real*8 TSurf            - Temperature of snow pack surface layer (C)
!   real*8 MeltEnergy       - Energy used for melting and heating of snow pack
!                             (W/m2)
!   real*8 Density	    - Density of the snow pack (g/cm3)
! ********************************************************************

      subroutine calcsnowmelt(y, x, Tstep, Z, Displacement, Z0,&
                               BaseRa, AirDens, EactAir, Lv,& 
                               ShortRad, LongRadIn, Press, RainFall,&
                               SnowFall, Tair, Vpd, Wind,& 
                               PackWater, SurfWater, Swq,& 
                               VaporMassFlux, TPack, TSurf,&
                               MeltEnergy, Outflow,&
                               LatentHeat, SensibleHeat, NetRad, Emiss,&
                               Density, Gact)

      !implicit none
      parameter (CH_ICE=2100.0d3)
!         Volumetri! heat capacity (J/(m3*C) of ice (0C)
      parameter (CH_WATER=4186.8d3)
!         Volumetri! heat capacity (J/(m3*C) of water
      parameter (CP=1013.d0)
!         Specifi! heat of moist air at constant pressure (J/(kg*C))
      parameter (DELTAT=25.d0)
!         Used in SensibleHeatFlux to bracket the
!         effective surface temperature (C)
      parameter (EPS=0.622d0)
!         ratio of molecular weight of water vapor to that for dry air
      parameter (G=9.81d0)
!         gravitational accelleration (m/(s^2))
      parameter (GRAMSPKG=1000.d0)
!         grams per kilogram
      parameter (R_JOULESPCAL=4.1868d0)
!         Joules per calorie
      parameter (LF=333.7d3)
!         latent heat of fusion (J/kg)
      parameter (SECPHOUR=3600.d0)
      parameter (STEFAN=5.6696d-8)
!         Stefan-Boltzmann constant (W/(M^2*C^4)
      parameter (WATER_DENSITY=1000.d0)
!         Density of water in kg/m3
      parameter (DHSVM_HUGE=1e20)
      parameter (R_LIQUID_WATER_CAPACITY=0.03d0)
!         water holding capacity of snow as a fraction of
!         snow-water-equivalent
!CVAL      parameter (R_MAX_SURFACE_SWE=0.005d0)
      parameter (R_MAX_SURFACE_SWE=0.125d0)
!         maximum depth of the surface layer in water equivalent (m)
      parameter (Rate=0.002)
!	  rate of snow density increase (g/cm3/day)

      integer y,x

      real*8 Tstep, Z, Displacement, Z0, Gact
      real*8 BaseRa, AirDens, EactAir, Lv
      real*8 ShortRad, LongRadIn, Press, RainFall
      real*8 SnowFall, Tair, Vpd, Wind
      real*8 PackWater, SurfWater, Swq
      real*8 VaporMassFlux, TPack, TSurf
      real*8 MeltEnergy, aaa
      real*8 var1,var2,var3
      real*8 LatentHeat, SensibleHeat, NetRad, Emiss
      real*8 zero, Density

      real*8 DeltaPackCC
!          Change in cold content of the pack
      real*8 DeltaPackSwq           
!          Change in snow water equivalent of the pack (m)
      real*8 Ice                    
!          Ice content of snow pack (m)*/
      real*8 InitialSwq             
!          Initial snow water equivalent (m)
      real*8 MassBalanceError       
!          Mass balance error (m)
      real*8 MaxLiquidWater         
!          Maximum liquid water content of pack (m)
      real*8 OldTSurf               
!          Old snow surface temperature (C)
      real*8 Outflow                
!          Amount water flowing out of the snow pack
!          during the time interval (m)
      real*8 PackCC                 
!          Cold content of snow pack (J)
      real*8 PackSwq                
!          Snow pack snow water equivalent (m)
      real*8 Qnet                   
!          Net energy exchange at the surface (W/m2)
      real*8 RefreezeEnergy         
!          refreeze energy (W/m2)
      real*8 RefrozenWater          
!          Amount of refrozen water (m)
      real*8 SnowFallCC
!          Cold content of new snowfall (J)
      real*8 SnowMelt               
!          Amount of snow melt during time interval (m water equivalent)
      real*8 SurfaceCC              
!          Cold content of snow pack (J)
      real*8 SurfaceSwq             
!          Surface layer snow water equivalent (m)
      real*8 Dnew
! 	   Density of the new snow (g/cm3)
      real*8 Tfar
!	   Air temperature (F)
      real*8 P1,P2
!	   Amount of snow already present and amount of new snow

      zero=0.0d0

      Tfar=Tair*1.8d0+32.d0

!CVAL      Dnew=0.000015625d0*Tfar*Tfar+0.001375d0*Tfar+0.02d0
      Dnew=0.000015625d0*Tfar*Tfar+0.001375d0*Tfar+0.005d0

      if (Dnew.le.(0.02d0)) then

         Dnew=0.02d0

      endif

      if (Dnew.ge.(0.1d0)) then

         Dnew=0.1d0

      endif

      if (Density.eq.(0.0d0)) then

         Density=Dnew

      else

         Density = Density + 0.004*Tstep/24.d0

         P1=Swq
         P2=SnowFall
         Density = Density*P1/(P1+P2)+Dnew*P2/(P1+P2)

         if (Density.gt.(0.7d0)) then

            Density=0.7d0

         endif

      endif

      InitialSwq = Swq
      OldTSurf = TSurf

! Initialize snowpack variables 
  
      Ice  = Swq - PackWater - SurfWater
  
! Reconstruct snow pack

      if (Ice.gt.R_MAX_SURFACE_SWE) then

         SurfaceSwq = R_MAX_SURFACE_SWE

      else

         SurfaceSwq = Ice

      endif

      PackSwq = Ice - SurfaceSwq

! Calculate cold contents
      SurfaceCc = CH_ICE * SurfaceSwq * TSurf
      PackCc = CH_ICE * PackSwq * TPack

      if (Tair.gt.(zero)) then

         SnowFallCC = zero

      else

         SnowFallCC = CH_ICE * SnowFall * Tair

      endif
  
! Distribute fresh snowfall
      if (SnowFall.gt.(R_MAX_SURFACE_SWE - SurfaceSwq)) then

         DeltaPackSwq = SurfaceSwq + SnowFall - R_MAX_SURFACE_SWE

         if (DeltaPackSwq.gt.SurfaceSwq) then

            DeltaPackCC = (SnowFall - R_MAX_SURFACE_SWE)/SnowFall *&
                          SnowFallCC
            DeltaPackCC = DeltaPackCC + SurfaceCC

         else

            DeltaPackCC = DeltaPackSwq/SurfaceSwq * SurfaceCC

         endif

!CVAL In case the snow pack is divided in 2 layers add the ground
!CVAL heat flux to the snow pack.

         if (Swq.gt.R_MAX_SURFACE_SWE) then

            DeltaPackCC=DeltaPackCC+Gact

         endif

         SurfaceSwq = R_MAX_SURFACE_SWE
         SurfaceCC = SurfaceCC + SnowFallCC - DeltaPackCC
         PackSwq = PackSwq + DeltaPackSwq
         PackCC = PackCC + DeltaPackCC

      else

         SurfaceSwq = SurfaceSwq + SnowFall
         SurfaceCC = SurfaceCC + SnowFallCC

      endif

      if (SurfaceSwq.gt.(zero)) then

         TSurf = SurfaceCC/(CH_ICE * SurfaceSwq)

      else 

         TSurf = zero

      endif

      if (PackSwq.gt.(zero)) then

         TPack = PackCC/(CH_ICE * PackSwq)

      else

         TPack = zero

      endif

! Adjust ice and *SurfWater
      Ice = Ice + SnowFall
      SurfWater = SurfWater + RainFall
  
! Calculate the surface energy balance for snow_temp = zero

      var1=BaseRa
      var2=VaporMassFlux
      var3=SurfaceSwq

      call snowpackenergybalance(zero, Tstep, var1, Z,&
                                  Displacement, Z0, Wind, ShortRad,&
                                  LongRadIn, AirDens,& 
                                 Lv, Tair, Press, Vpd, EactAir,&
                                  RainFall, SurfaceSwq, SurfWater,&
                                  OldTSurf, RefreezeEnergy,&
                                  VaporMassFlux, Qnet,&
                                  LatentHeat, SensibleHeat, NetRad, Emiss,&
                                  Gact, R_MAX_SURFACE_SWE, Swq)

! If Qnet == zero, then set the surface temperature to zero
      if (Qnet.eq.(zero)) then

         TSurf = zero

         if (RefreezeEnergy.ge.(zero)) then

            RefrozenWater=RefreezeEnergy/(LF*WATER_DENSITY)*Tstep*SECPHOUR

            if (RefrozenWater.gt.SurfWater) then

               RefrozenWater = SurfWater
               RefreezeEnergy = RefrozenWater*LF*WATER_DENSITY/&
                                (Tstep*SECPHOUR)

            endif

            MeltEnergy = MeltEnergy + RefreezeEnergy
            SurfaceSwq = SurfaceSwq + RefrozenWater
            Ice = Ice + RefrozenWater
            SurfWater = SurfWater - RefrozenWater

            if (SurfWater.lt.zero) then

               write (*,*) 'Surfwater lt. zero ',y,x

            endif

            SnowMelt = zero

         else
      
! Calculate snow melt
            SnowMelt=abs(RefreezeEnergy)/(LF*WATER_DENSITY)*Tstep*SECPHOUR
            MeltEnergy = MeltEnergy + RefreezeEnergy

         endif

! Convert vapor mass flux to a depth per timestep and adjust SurfWater

         VaporMassFlux = VaporMassFlux * Tstep * SECPHOUR
    
         if (SurfWater.lt.(-VaporMassFlux)) then

            VaporMassFlux = -(SurfWater)
            SurfWater    = zero

         else

            SurfWater = SurfWater + VaporMassFlux

         endif
    
! If SnowMelt < Ice, there was incomplete melting of the pack
    
         if (SnowMelt.lt.Ice) then

             if (SnowMelt.le.PackSwq) then

                SurfWater = SurfWater + SnowMelt
                PackSwq = PackSwq - SnowMelt
                Ice = Ice - SnowMelt

             else

                SurfWater = SurfWater + SnowMelt
                PackWater = zero
                PackSwq = zero
                Ice = Ice - SnowMelt
                SurfaceSwq  = Ice

             endif

! Else, SnowMelt > Ice and there was complete melting of the pack
         else

            SnowMelt = Ice
            SurfWater = SurfWater + Ice
            SurfaceSwq = zero
            TSurf = zero
            PackSwq = zero
            TPack = zero
            Ice = zero

         endif
  
! Else, SnowPackEnergyBalance(T=zero) <= zero
      else
! Calculate surface layer temperature using "Brent method"

         var1=VaporMassFlux
         var2=BaseRa
         aaa=zero
         call rootbrent(y, x, TSurf-DELTAT, zero,&
                         aaa, Tstep, BaseRa, Z,&
                         Displacement, Z0, Wind, ShortRad,&
                         LongRadIn, AirDens, Lv, Tair, &
                        Press, Vpd, EactAir, RainFall, SurfaceSwq,&
                         SurfWater, OldTSurf, RefreezeEnergy,&
                         VaporMassFlux, TSurf,&
                         LatentHeat, SensibleHeat, NetRad, Emiss,&
                         Gact, R_MAX_SURFACE_SWE, Swq)

! Since we iterated, the surface layer is below freezing and no snowmelt
    
         SnowMelt = zero
    
! Since updated snow_temp < zero, all of the liquid water
! in the surface layer has been frozen

         SurfaceSwq = SurfaceSwq + SurfWater
         Ice = Ice + SurfWater
         SurfWater = zero
         MeltEnergy=MeltEnergy+SurfWater*LF*WATER_DENSITY/(Tstep*SECPHOUR)

! Convert mass flux to a depth per timestep and adjust SurfaceSwq
    
         VaporMassFlux = VaporMassFlux * Tstep * SECPHOUR
         if (SurfaceSwq.lt.(-VaporMassFlux)) then

            VaporMassFlux = -SurfaceSwq
            SurfaceSwq = zero
            Ice = PackSwq

         else

            SurfaceSwq = SurfaceSwq + VaporMassFlux
            Ice = Ice + VaporMassFlux

         endif

      endif

! Done with iteration etc, now Update the liquid water content of the
! surface layer
  
      MaxLiquidWater = R_LIQUID_WATER_CAPACITY * SurfaceSwq

      if  (SurfWater.gt.MaxLiquidWater) then

          Outflow = SurfWater - MaxLiquidWater
          SurfWater = MaxLiquidWater
    
      else
    
          Outflow = zero

      endif
  
! Refreeze liquid water in the pack.                                   
! variable 'RefreezeEnergy' is the heat released to the snow pack            
! if all liquid water were refrozen.                                   
! if RefreezeEnergy < PackCC then all water IS refrozen           
! PackCC always <=zero 

! WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does 
! not involve energy transported to the pixel.  Instead heat from the snow 
! pack is used to refreeze water.
  
      PackWater = PackWater + Outflow
! Add surface layer outflow to pack liquid water
      RefreezeEnergy = PackWater * LF * WATER_DENSITY

! calculate energy released to freeze
  
      if (PackCC.lt.(-RefreezeEnergy)) then
! cold content not fully depleted

         PackSwq = PackSwq + PackWater 
! refreeze all water and update
         Ice = Ice + PackWater
         PackWater = zero

         if (PackSwq.gt.(zero)) then

            TPack = (PackCC + RefreezeEnergy) / (CH_ICE * PackSwq)

         else 

            TPack = zero

         endif

      else

! cold content has been either exactly satisfied or exceeded. If
! PackCC = refreeze then pack is ripe and all pack water is
! refrozen, else if energy released in refreezing exceeds PackCC 
! then exactly the right amount of water is refrozen to satify PackCC.
! The refrozen water is added to PackSwq and Ice

         TPack = zero
         DeltaPackSwq = -PackCC/(LF * WATER_DENSITY) 
         PackWater = PackWater - DeltaPackSwq
         PackSwq = PackSwq + DeltaPackSwq
         Ice = Ice + DeltaPackSwq

      endif
  
! Update the liquid water content of the pack
  
      MaxLiquidWater = R_LIQUID_WATER_CAPACITY * PackSwq

      if (PackWater.gt.MaxLiquidWater) then

         Outflow = PackWater - MaxLiquidWater
         PackWater = MaxLiquidWater

      else

         Outflow = zero

      endif

! Update snow properties

      Ice = PackSwq + SurfaceSwq

      if (Ice.gt.R_MAX_SURFACE_SWE) then

         SurfaceCC = CH_ICE * TSurf * SurfaceSwq
         PackCC = CH_ICE * TPack * PackSwq

         if (SurfaceSwq.gt.R_MAX_SURFACE_SWE) then

            PackCC = PackCC + SurfaceCC *&
                     (SurfaceSwq - R_MAX_SURFACE_SWE) / SurfaceSwq
            SurfaceCC = SurfaceCC - SurfaceCC *&
                        (SurfaceSwq - R_MAX_SURFACE_SWE) / SurfaceSwq
            PackSwq = PackSwq + SurfaceSwq - R_MAX_SURFACE_SWE
            SurfaceSwq = SurfaceSwq - SurfaceSwq + R_MAX_SURFACE_SWE

         else

            if (SurfaceSwq.lt.R_MAX_SURFACE_SWE) then

               PackCC = PackCC - PackCC *&
                        (R_MAX_SURFACE_SWE - SurfaceSwq) / PackSwq
               SurfaceCC = SurfaceCC + PackCC *&
                           (R_MAX_SURFACE_SWE - SurfaceSwq) / PackSwq
               PackSwq = PackSwq - R_MAX_SURFACE_SWE + SurfaceSwq
               SurfaceSwq = SurfaceSwq + R_MAX_SURFACE_SWE - SurfaceSwq
    
            endif

         endif

         TPack = PackCC / (CH_ICE * PackSwq)
         TSurf = SurfaceCC / (CH_ICE * SurfaceSwq)

      else

         PackSwq = zero
         PackCC = zero
         TPack = zero

      endif

      Swq = Ice + PackWater + SurfWater

      if (Swq.le.(zero)) then

         TSurf = zero
         TPack = zero

      endif

      if (Swq.le.R_MAX_SURFACE_SWE) then

         TPack = zero

      endif

! Mass balance test
  
      MassBalanceError = (InitialSwq - Swq) + (RainFall + SnowFall) -&
                         Outflow + VaporMassFlux 

      if ( (TSurf.le.(Tair-273.d0)).or.&
           (TPack.le.(Tair-273.d0)) ) then

         if (Swq.le.(0.1d0*R_MAX_SURFACE_SWE)) then

            Outflow=Swq
            Swq=0.d0

            TSurf=Tair
            TPack=Tair

         else

            if ( ((Swq-SurfaceSwq).le.(0.1d0*R_MAX_SURFACE_SWE)).and.&
                 (Swq.gt.R_MAX_SURFACE_SWE) ) then

               Outflow=Swq-SurfaceSwq
               Swq=R_MAX_SURFACE_SWE
               PackSwq=0.d0

               TPack=0.d0

            else

               write (*,*) y, x, Tstep, Z, Displacement, Z0,&
                            BaseRa, AirDens, EactAir, Lv,&
                            ShortRad, LongRadIn, Press, RainFall,&
                            SnowFall, Tair, Vpd, Wind,&
                            PackWater, SurfWater, ' SWQ ',Swq,&
                            VaporMassFlux, TPack, TSurf,&
                            MeltEnergy, Outflow,&
                            LatentHeat, SensibleHeat, NetRad, Emiss,&
                            Density, Gact, MassBalanceError, ' PACKSWQ ',&
                            PackSwq, SurfaceSwq, PackCC, SurfaceCC,&
                            Qnet, RefreezeEnergy

!CVAL               stop

            endif

         endif

      endif

      return

      end subroutine calcsnowmelt


!*************************************************************************
!* SUMMARY:      RootBrent.! - Determine surface temperature iteratively
!* USAGE:        Part of DHSVM
!*
!* AUTHOR:       Bart Nijssen
!* ORG:          University of Washington, Department of Civil Engineering
!* E-MAIL:       nijssen@u.washington.edu
!* ORIG-DATE:    Apr-96
!* LAST-MOD:      2-Oct-1996 at 16:25:36 by Bart Nijssen
!* DESCRIPTION:  Determine surface temperature iteratively using the Brent
!*               method.  
!* DESCRIP-END.
!* FUNCTIONS:    RootBrent()
!* COMMENTS:     
!*************************************************************************

!*****************************************************************************
! GENERAL DOCUMENTATION FOR THIS MODULE
! -------------------------------------

! Source: Brent, R. P., 1973, Algorithms for minimization without derivatives,&
 !                       Prentice Hall, Inc., Englewood Cliffs, New Jersey
!			Chapter 4
! This source includes an implementation of the algorithm in ALGOL-60, which
! was translated into ! for this application.

! The method is also discussed in:
! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, 1992,&
 !               Numerical Recipes in FORTRAN, The art of scientifi! computing,&
 !		Second edition, Cambridge University Press
! (Be aware that this book discusses a Brent method for minimization (brent), 
! and one for root finding (zbrent).  The latter one is similar to the one 
! implemented here and is also copied from Brent [1973].)

! The function returns the surface temperature, TSurf, for which the sum
! of the energy balance terms is zero, with TSurf in the interval 
! [MinTSurf, MaxTSurf].  The surface temperature is calculated to within
! a tolerance (6 * MACHEPS * |TSurf| + 2 * T), where MACHEPS is the relative
! machine precision and T is a positive tolerance, as specified in brent.h.

! The function assures that f(MinTSurf) and f(MaxTSurf) have opposite signs.
! If this is not the case the program will abort.  In addition the program
! will perform not more than a certain number of iterations, as specified
! in brent.h, and will abort if more iterations are needed.
!*****************************************************************************
  
!*****************************************************************************
! Function name: RootBrent()

! Purpose      : Calculate the surface temperature in the absence of snow

! Required     :
!   int y                   - Row number of current pixel
!   int x                   - Column number of current pixel 
!   float LowerBound        - Lower bound for root
!   float UpperBound        - Upper bound for root
!   float (*Function)(float Estimate, va_list ap)
!   ...                     - Variable arguments 
!                             The number and order of arguments has to be
!                             appropriate for the Function pointed to, since
!                             the list of arguments after Nargs will be passed
!                             on to Function.
!                             See the appropriate Function for the correct
!                             arguments. 
!   real*8 Tstep              - Model time step (hours) 
!   real*8 Ra                 - Aerodynami! resistance (s/m) 
!   real*8 Z                  - Reference height (m) 
!   real*8 Displacement       - Displacement height (m) 
!   real*8 Z0                 - Roughness length (m) 
!   real*8 Wind               - Wind speed (m/s) 
!   real*8 ShortRad           - Net incident shortwave radiation (W/m2) 
!   real*8 LongRadIn          - Incoming longwave radiation (W/m2) 
!   real*8 AirDens            - Density of air (kg/m3) 
!   real*8 Lv                 - Latent heat of vaporization (J/kg) 
!   real*8 Tair               - Air temperature (C) 
!   real*8 Press              - Air pressure (Pa) 
!   real*8 Vpd                - Vapor pressure deficit (Pa) 
!   real*8 EactAir            - Actual vapor pressure of air (Pa) 
!   real*8 Rain               - Rain fall (m/timestep) 
!   real*8 SweSurfaceLayer    - Snow water equivalent in surface layer (m)
                           
!   real*8 SurfaceLiquidWater - Liquid water in the surface layer (m) 
!   real*8 OldTSurf           - Surface temperature during previous time step 
!   real*8 RefreezeEnergy     - Refreeze energy (W/m2) 
!   real*8 VaporMassFlux      - Mass flux of water vapor to or from
!                             the intercepted snow 


! Returns      :
!   real*8 b                  - Effective surface temperature (C)
!   real*8 LatentHeat	    - Latent heat flux (W/m2)
!   real*8 SensibleHeat	    - Sensible heat flux (W/m2)

! Modifies     : none

!****************************************************************************

      subroutine rootbrent(y, x, LowerBound, UpperBound, aaa,&
                            Tstep, BaseRa, Z, Displacement, Z0, Wind,&
                            ShortRad, LongRadIn, AirDens, Lv, Tair,& 
                           Press, Vpd, EactAir, RainFall, SurfaceSwq,&
                            SurfWater, OldTSurf, RefreezeEnergy,&
                            VaporMassFlux, b, LatentHeat, SensibleHeat,&
                            NetRad, Emiss, Gact, R_MAX_SURFACE_SWE, Swq)

      parameter (MACHEPS=0.00000008d0)
!         machine floating point precision (float)
      parameter (T=1d-5)
!         tolerance
      parameter (MAXITER=1000)
!         maximum number of allowed iterations
      parameter (MAXTRIES=5)
!         maximum number of tries to bracket the root
      parameter (T_STEP=10)
!         step to take in both directions if attempting to
!         bracket the root 

      integer y,x

      real*8 LowerBound, UpperBound, aaa, Gact
      real*8 Tstep, BaseRa, Z, Displacement, Z0, Wind
      real*8 ShortRad, LongRadIn, AirDens, Lv, Tair
      real*8 Press, Vpd, EactAir, RainFall, SurfaceSwq
      real*8 SurfWater, OldTSurf, RefreezeEnergy
      real*8 VaporMassFlux
      real*8 var1,var2,var3
      real*8 LatentHeat, SensibleHeat, zero, NetRad, Emiss, Swq

      real*8 a
      real*8 b
      real*8 c
      real*8 d
      real*8 e
      real*8 fa
      real*8 fb
      real*8 fc
      real*8 m
      real*8 p
      real*8 q
      real*8 r
      real*8 s
      real*8 tol
      integer i
      integer j
      integer eval

! initialize variable argument list

      zero=0.0d0

      eval=0

      a = LowerBound
      b = UpperBound
      var1=BaseRa
      var2=VaporMassFlux
      var3=RefreezeEnergy
      call snowpackenergybalance(a,Tstep,var1,Z,Displacement,Z0,Wind,&
                                  ShortRad,LongRadIn,AirDens,Lv,Tair,&
                                  Press,Vpd,EactAir,RainFall,&
                                  SurfaceSwq,&
                                  SurfWater,OldTSurf,&
                                  RefreezeEnergy,VaporMassFlux,fa,&
                                  LatentHeat, SensibleHeat, NetRad,&
                                  Emiss, Gact, R_MAX_SURFACE_SWE, Swq)
      eval=eval+1
      var1=BaseRa
      var2=VaporMassFlux
      var3=RefreezeEnergy
      call snowpackenergybalance(b,Tstep,var1,Z,Displacement,Z0,Wind,&
                                  ShortRad,LongRadIn,AirDens,Lv,Tair,&
                                  Press,Vpd,EactAir,RainFall,&
                                  SurfaceSwq,&
                                  SurfWater,OldTSurf,&
                                  RefreezeEnergy,VaporMassFlux,fb,&
                                  LatentHeat, SensibleHeat, NetRad,&
                                  Emiss, Gact, R_MAX_SURFACE_SWE, Swq)
      eval=eval+1
      
! if root not bracketed attempt to bracket the root

      j=0

      if ( ((fa*fb).ge.zero).and.(j.lt.MAXTRIES) ) then

1000     a = a - T_STEP
         b = b + T_STEP
         var1=BaseRa
         var2=VaporMassFlux
         var3=RefreezeEnergy
         call snowpackenergybalance(a,Tstep,var1,Z,Displacement,Z0,Wind,&
                                     ShortRad,LongRadIn,AirDens,Lv,Tair,&
                                     Press,Vpd,EactAir,RainFall,&
                                     SurfaceSwq,&
                                     SurfWater,OldTSurf,&
                                     RefreezeEnergy,VaporMassFlux,fa,&
                                     LatentHeat, SensibleHeat, NetRad,&
                                     Emiss, Gact, R_MAX_SURFACE_SWE, Swq)
         eval=eval + 1
         var1=BaseRa
         var2=VaporMassFlux
         var3=RefreezeEnergy
         call snowpackenergybalance(b,Tstep,var1,Z,Displacement,Z0,Wind,&
                                     ShortRad,LongRadIn,AirDens,Lv,Tair,&
                                     Press,Vpd,EactAir,RainFall,&
                                     SurfaceSwq,&
                                     SurfWater,OldTSurf,&
                                     RefreezeEnergy,VaporMassFlux,fb,&
                                     LatentHeat, SensibleHeat, NetRad,&
                                     Emiss, Gact, R_MAX_SURFACE_SWE, Swq)
         eval=eval + 1
         j = j + 1

         if ( ((fa*fb).ge.zero).and.(j.lt.MAXTRIES) ) then

            goto 1000

         endif

      endif

      if ((fa * fb).ge.zero) then

!CVAL         write (*,*) 'fa x fb equals ',fa * fb

      endif
  
      fc = fb

      do i=1,MAXITER

         if ((fb*fc).gt.zero) then

           c  = a
           fc = fa
           d = b - a
           e = d

         endif

         if (abs(fc).lt.abs(fb)) then

            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
     
         endif

         tol = 2.d0 * MACHEPS * abs(b) + T
         m = 0.5d0 * (c - b)

         if ( (abs(m).le.tol).or.(fb.eq.zero) ) then

            aaa=fb
            return

         else

            if ( (abs(e).lt.tol).or.(abs(fa).le.abs(fb)) ) then

	       d = m
	       e = d

            else

	       s = fb/fa
	
	       if (a.eq.c) then
	  
! linear interpolation
          
	          p = 2.d0 * m * s
	          q = 1.d0 - s

	       else
	  
! inverse quadrati! interpolation
	  
	          q = fa/fc
	          r = fb/fc
	          p = s * (2.d0 * m * q * (q - r) - (b - a) * (r - 1.d0))
	          q = (q - 1.d0) * (r - 1) * (s - 1)

	       endif
	
	       if (p.gt.zero) then

	           q = -q

	       else

	           p = -p

               endif

               s = e
               e = d

	       if ( (2.d0 * p).lt.( 3.d0 * m * q - abs(tol * q)).and.&
                   ( p.lt.abs(0.5d0 * s * q)) ) then

	          d = p/q

	       else

	          d = m
	          e = d

	       endif

            endif

            a = b
            fa = fb

            if (abs(d).gt.tol) then
 
               b=b+d

            else

               if (m.gt.zero) then

                  b=b+tol

               else

                  b=b-tol

               endif
                 
            endif

            var1=BaseRa
            var2=VaporMassFlux
            var3=RefreezeEnergy
            call snowpackenergybalance(b,Tstep,var1,Z,Displacement,Z0,&
                                        Wind,&
                                        ShortRad,LongRadIn,AirDens,Lv,&
                                        Tair,&
                                        Press,Vpd,EactAir,RainFall,&
                                        SurfaceSwq, SurfWater,&
                                        OldTSurf,&
                                        RefreezeEnergy,VaporMassFlux,fb,&
                                        LatentHeat, SensibleHeat, NetRad,&
                                        Emiss, Gact, R_MAX_SURFACE_SWE, Swq)
            eval=eval+1

         endif
  
      enddo

      write (*,*) 'Error 33 '

      return

      end subroutine rootbrent

! * SUMMARY:      SnowPackEnergyBalance.! - Calculate snow pack energy balance
! * USAGE:        Part of DHSVM
! *
! * AUTHOR:       Bart Nijssen
! * ORG:          University of Washington, Department of Civil Engineering
! * E-MAIL:              nijssen@u.washington.edu
! * ORIG-DATE:     8-Oct-1996 at 09:09:29
! * LAST-MOD: Wed Apr  9 15:08:01 1997 by Bart Nijssen
!                 <nijssen@meter.ce.washington.edu>
! * DESCRIPTION:  Calculate snow pack energy balance
! * DESCRIP-END.
! * FUNCTIONS:    SnowPackEnergyBalance()
! * COMMENTS:     

!*****************************************************************************
! Function name: SnowPackEnergyBalance()

! Purpose      : Calculate the surface energy balance for the snow pack

! Required     :
!   real*8 TSurf              - new estimate of effective surface temperature
!   real*8 Tstep              - Model time step (hours) 
!   real*8 Ra                 - Aerodynami! resistance (s/m) 
!   real*8 Z                  - Reference height (m) 
!   real*8 Displacement       - Displacement height (m) 
!   real*8 Z0                 - Roughness length (m) 
!   real*8 Wind               - Wind speed (m/s) 
!   real*8 ShortRad           - Net incident shortwave radiation (W/m2) 
!   real*8 LongRadIn          - Incoming longwave radiation (W/m2) 
!   real*8 AirDens            - Density of air (kg/m3) 
!   real*8 Lv                 - Latent heat of vaporization (J/kg) 
!   real*8 Tair               - Air temperature (C) 
!   real*8 Press              - Air pressure (Pa) 
!   real*8 Vpd                - Vapor pressure deficit (Pa) 
!   real*8 EactAir            - Actual vapor pressure of air (Pa) 
!   real*8 Rain               - Rain fall (m/timestep) 
!   real*8 SweSurfaceLayer    - Snow water equivalent in surface layer (m)
                           
!   real*8 SurfaceLiquidWater - Liquid water in the surface layer (m) 
!   real*8 OldTSurf           - Surface temperature during previous time step 
!   real*8 RefreezeEnergy     - Refreeze energy (W/m2) 
!   real*8 VaporMassFlux      - Mass flux of water vapor to or from
!                             the intercepted snow 


! Returns      :
!   real*8 RestTerm         - Rest term in the energy balance
!   real*8 LatentHeat	  - Latent heat flux (W/m2)
!   real*8 SensibleHeat	  - Sensible heat flux (W/m2)

! Modifies     : 
!   real*8 RefreezeEnergy   - Refreeze energy (W/m2) 
!   real*8 VaporMassFlux    - Mass flux of water vapor to or from the
!                           intercepted snow 

! Comments     :
!   Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
!               science, Addisson Wesley, Inc., Reading, etc., 1990.
!****************************************************************************

      subroutine snowpackenergybalance(TSurf,Tstep,Ra,Z,Displacement,Z0,&
                                        Wind, ShortRad,LongRadIn,AirDens,&
                                        Lv,Tair, Press,Vpd,EactAir,Rain,&
                                        SweSurfaceLayer,&
                                        SurfaceLiquidWater,OldTSurf,&
                                        RefreezeEnergy, VaporMassFlux,&
                                        RestTerm,&
                                        LatentHeat, SensibleHeat, NetRad,&
                                        Emiss, Gact, R_MAX_SURFACE_SWE, Swq)

      parameter (DHSVM_HUGE=1e20)
      parameter (CH_ICE=2100.0d3)
!         Volumetri! heat capacity (J/(m3*C) of ice (0C)
      parameter (CH_WATER=4186.8d3)
!         Volumetri! heat capacity (J/(m3*C) of water
      parameter (CP=1013.0d0)
!         Specifi! heat of moist air at constant pressure (J/(kg*C))
      parameter (EPS=0.622d0)
!         ratio of molecular weight of water vapor to that for dry air
      parameter (GRAMSPKG=1000.d0)
!         grams per kilogram
      parameter (R_JOULESPCAL=4.1868d0)
!         Joules per calorie
      parameter (LF=333.7d3)
!         latent heat of fusion (J/kg)
      parameter (SECPHOUR=3600.d0)
      parameter (STEFAN=5.6696d-8)
!         Stefan-Boltzmann constant (W/(M^2*C^4)
      parameter (WATER_DENSITY=1000.d0)
!         Density of water in kg/m3

      real*8 TSurf,Tstep,Ra,Z,Displacement,Z0,Wind,Swq,Gact
      real*8 ShortRad,LongRadIn,AirDens,Lv,Tair
      real*8 Press,Vpd,EactAir,Rain,SweSurfaceLayer
      real*8 SurfaceLiquidWater,OldTSurf
      real*8 RefreezeEnergy, VaporMassFlux,Tmean, zero, Emiss

      real*8 AdvectedEnergy         
!          Energy advected by precipitation (W/m2) 
      real*8 DeltaColdContent       
!          Change in cold content (W/m2) 
      real*8 EsSnow                 
!          saturated vapor pressure in the snow pack (Pa)  
      real*8 LatentHeat		
!          Latent heat exchange at surface (W/m2) 
      real*8 LongRadOut		
!          long wave radiation emitted by surface (W/m2) 
      real*8 Ls                     
!          Latent heat of sublimation (J/kg) 
      real*8 NetRad			
!          Net radiation exchange at surface (W/m2) 
      real*8 RestTerm		
!          Rest term in surface energy balance (W/m2) 
      real*8 SensibleHeat		
!          Sensible heat exchange at surface (W/m2) 
!          Mean temperature during interval (C) 
      real*8 Tmp                   
!            temporary variable 

! Calculate active temp for energy balance as average of old and new

      zero=0.0d0
  
      TMean = 0.5d0 * (OldTSurf + TSurf)

! Correct aerodynami! conductance for stable conditions
!    Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
!    velocity (vel_2m) is expected to be in m/sec
! 
! Apply the stability correction to the aerodynami! resistance 
!    NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
!    think that it is more correct to calculate ALL fluxes at the same
!    reference level

      if (Wind.gt.(zero)) then

         Ra = Ra / stabilitycorrection(2.d00, zero, TMean, Tair, Wind, Z0)

      else

         Ra = DHSVM_HUGE

      endif

! Calculate longwave exchange and net radiation

      Tmp = TMean+273.15d0
      LongRadOut = Emiss * STEFAN * (Tmp*Tmp*Tmp*Tmp)
      NetRad = ShortRad + LongRadIn - LongRadOut
  
! Calculate the sensible heat flux

      SensibleHeat=AirDens*CP*(Tair-Tmean)/Ra

! Calculate the mass flux of ice to or from the surface layer
 
! Calculate the saturated vapor pressure in the snow pack, 
!    (Equation 3.32, Bras 1990)

      call satvaporpressure(Tmean,EsSnow)

      VaporMassFlux = AirDens * (EPS/Press) * (EactAir - EsSnow)/Ra
      VaporMassFlux = VaporMassFlux / WATER_DENSITY
      if ( (Vpd.eq.(zero)).and.(VaporMassFlux.lt.(zero)) ) then

         VaporMassFlux = zero

      endif

      if ( (VaporMassFlux.ge.-5000.d0).and.(VaporMassFlux.le.5000.d0) ) then

         VaporMassFlux=VaporMassFlux

      else

         write (*,*) 'CALCSNOWMELT : VaporMassFlux out of bounds ',VaporMassFlux
         write (*,*) AirDens,EPS,Press,EactAir,EsSnow,Ra,WATER_DENSITY,Tmean

      endif
  
! Calculate latent heat flux
   
      if (TMean.ge.(zero)) then
! Melt conditions: use latent heat of vaporization

         LatentHeat = Lv * VaporMassFlux * WATER_DENSITY

      endif

      if ( (LatentHeat.ge.-5000.d0).and.(LatentHeat.le.5000.d0) ) then

        ! Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990) 

         Ls = (677.d0 - 0.07d0 * TMean) * R_JOULESPCAL * GRAMSPKG
         LatentHeat=LatentHeat
         LatentHeat = Ls * VaporMassFlux * WATER_DENSITY

      else

         write (*,*) 'CALCSNOWMELT : LatentHeat out of bounds ',LatentHeat
         write (*,*) TMean,VaporMassFlux,WATER_DENSITY
         write (*,*) R_JOULESPCAL,GRAMSPKG
         write (*,*) Ls,Lv
         write (*,*)

      endif
  
! Calculate advected heat flux from rain 
!     WORK IN PROGRESS:  Should the following read (Tair - Tsurf) ??
  
      AdvectedEnergy = (CH_WATER * Tair * Rain) / (Tstep*SECPHOUR)
  
! Calculate change in cold content

      DeltaColdContent = CH_ICE * SweSurfaceLayer * (TSurf - OldTSurf)/&
                         (Tstep * SECPHOUR)

!CVAL Add the ground heat flux to the change in cold content
!CVAL if the top layer is smaller than the maximum top layer depth.

      if (Swq.le.R_MAX_SURFACE_SWE) then

         DeltaColdContent = DeltaColdContent+Gact

      endif

! Calculate net energy exchange at the snow surface

      RestTerm = NetRad + SensibleHeat + LatentHeat + AdvectedEnergy -&
                 DeltaColdContent

      RefreezeEnergy = (SurfaceLiquidWater * LF * WATER_DENSITY)/&
                       (Tstep * SECPHOUR)

      if ( (TSurf.eq.(zero)).and.(RestTerm.gt.(-RefreezeEnergy)) ) then

         RefreezeEnergy = -RestTerm 
! available energy input over cold content
! used to melt, i.e. Qrf is negative value
! (energy out of pack)
         RestTerm = zero

      else

         RestTerm = RestTerm + RefreezeEnergy
! add this positive value to the pack
  
      endif

      return

      end subroutine snowpackenergybalance


!*
!* SUMMARY:      StabilityCorrection.! - Calculate the stability correction
!* USAGE:        Part of DHSVM
!*
!* AUTHOR:       Bart Nijssen and Pascal Storck
!* ORG:          University of Washington, Department of Civil Engineering
!* E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
!* ORIG-DATE:    Apr-1996
!* LAST-MOD:     Tue Mar 25 16:54:52 1997 by Bart Nijssen
!*               <nijssen@meter.ce.washington.edu>
!* DESCRIPTION:  Calculate the stability correction for exchange of sensible
!*               heat between the surface and the atmosphere 
!* DESCRIP-END.
!* FUNCTIONS:    StabilityCorrection()
!* COMMENTS:     
!*

!*****************************************************************************
! Function name: StabilityCorrection()

! Purpose      : Calculate atmospheri! stability correction for non-neutral
!                conditions

! Required     :
!   real*8 Z          - Reference height (m)
!   real*8 d          - Displacement height (m)
!   real*8 TSurf      - Surface temperature (C)
!   real*8 Tair       - Air temperature (C)
!   real*8 Wind       - Wind speed (m/s)
!   real*8 Z0         - Roughness length (m)

! Returns      :
!   real*8 Correction - Multiplier for aerodynami! resistance

! Modifies     : None
    
!****************************************************************************
      function  stabilitycorrection(Z,d,TSurf,Tair,Wind,Z0)

      parameter (G=9.81)
!         gravitational accelleration (m/(s^2))

      real*8 Z,d,TSurf,Tair,Wind,Z0

      real*8 Correction		
!          Correction to aerodynami! resistance
      real*8 Ri                     
!          Richardson's Number
      real*8 RiCr             
!          Critical Richardson's Number
      real*8 RiLimit                
!          Upper limit for Richardson's Number
      real*8 zero

      Correction = 1.0d0
      RiCr = 0.2d0
      zero=0.0d0

! Calculate the effect of the atmospheri! stability using a Richardson 
! Number approach

      if (TSurf.ne.Tair) then

! Non-neutral conditions
    
         Ri = G * (Tair - TSurf) * (Z - d)/&
           (((Tair+273.15d0)+(TSurf+273.15d0))/2.d0 * Wind * Wind)
    
         RiLimit = (Tair+273.15d0)/&
           (((Tair+273.15d0)+(TSurf+273.15d0))/2.d0 * (log((Z-d)/Z0) + 5.d0))
    
         if (Ri.gt.RiLimit) then

            Ri = RiLimit

         endif
    
         if (Ri.gt.(zero)) then

            Correction = (1.d0 - Ri/RiCr) * (1.d0 - Ri/RiCr)
      
         else 

            if (Ri.lt.(-0.5d0)) then

               Ri = -0.5d0

            endif

            Correction = sqrt(1.d0 - 16.d0 * Ri)

         endif

      endif

      StabilityCorrection=Correction

      return

      end function stabilitycorrection

      subroutine satvaporpressure(tair,satvap)

      real*8 tair,satvap,zero

      zero=0.0d0

      if (tair.gt.(-10.d0)) then

         satvap= 610.78d0 * exp((17.269d0 * tair) / (237.3d0 + tair))

         if (tair.lt.zero) then

            satvap=satvap*(1.0d0+0.00972d0*tair+.000042d0*tair*tair)

         endif

      else

         satvap= 610.78d0 * exp((17.269d0 * (-10.d0)) / (237.3d0 + (-10.d0)))

         if ((-10.d0).lt.zero) then

            satvap=satvap*(1.0d0+0.00972d0*(-10.d0)+.000042d0*(-10.d0)*(-10.d0))

         endif

      endif

      return

      end subroutine satvaporpressure

! ====================================================================
!
!             subroutine snow_density
!
! ====================================================================
!
! Calculate the density of the snow pack.
!
! ====================================================================
!
! Input variables:
!
! snow_m        New snow                        m
! tcel          Air temperature                 C
! Swq_m         Snow water equivalent           m
! snow_depth    Snow depth                      m
! Tsurf         Snow surface temperature        C
! tstep         Time step                       s
!
! Output variables:
!
! dens          Snow pack density               kg/m3
!
! ====================================================================

      subroutine snow_density(dens,snow_m,tcel,Swq_m,snow_depth,Tsurf,tstep)

      implicit none

      real*8 tstep
      real*8 snow_mm,Tfar,dens_new,delta_depth,depth_new,dt
      real*8 dens,snow_m,tcel,Swq_m,snow_depth,Tsurf
      real*8 overburden,viscosity,deltadepth
      real*8 NEW_SNOW_DENSITY,MAX_CHANGE,ETA0,G,C5,C6,RHO_W,SECPHOUR

      parameter (NEW_SNOW_DENSITY=50.d0)
      parameter (MAX_CHANGE=0.9d0)
      parameter (ETA0=3600000.d0)
      parameter (G=9.81)
      parameter (SECPHOUR=3600.d0)
      parameter (C5=0.08d0)
      parameter (C6=0.021d0)
      parameter (RHO_W=997.d0)

      dt=tstep/SECPHOUR
      snow_mm=1000.*snow_m

      if (snow_mm.gt.0.d0) then

! ====================================================================
! Estimate density of new snow based on air temperature.
! ====================================================================

         Tfar=32.+tcel*9./5.

         if (Tfar.gt.0.d0) then

            dens_new=NEW_SNOW_DENSITY+1000.d0*(Tfar/100.)*(Tfar/100.)

         else

            dens_new=NEW_SNOW_DENSITY

         endif

! ====================================================================
! Compact current snowpack by weight of new snowfall.
! ====================================================================

         if (snow_depth.gt.0.d0) then

            delta_depth=( ((snow_mm/25.4d0)*(snow_depth/0.0254d0))/&
                          (Swq_m/0.0254d0)*&
                          ( ((snow_depth/0.0254d0)/10.d0)**(0.35d0) )&
                        )*0.0254d0

            if (delta_depth.ge.snow_depth) then

               delta_depth=MAX_CHANGE*snow_depth

            endif

            depth_new=snow_mm/dens_new
            snow_depth=snow_depth-delta_depth+depth_new
            dens=1000.d0*Swq_m/snow_depth
         else

! ====================================================================
! No snowpack present, so snow density equals that of new snow.
! ====================================================================

            dens=dens_new
            snow_depth=1000.d0*Swq_m/dens

         endif

      else

          dens=1000.d0*Swq_m/snow_depth

      endif

! ====================================================================
! Densification of the snow pack due to aging.
! ====================================================================

      snow_depth=1000.d0*Swq_m/dens
      overburden=0.5d0*G*RHO_W*Swq_m
      viscosity=ETA0*exp(-C5*Tsurf + C6*dens)
      deltadepth=-overburden/viscosity*snow_depth*dt*SECPHOUR
      snow_depth=snow_depth+deltadepth
      dens=1000.d0*Swq_m/snow_depth

      return

      end subroutine snow_density

END MODULE MODULE_SNOW
