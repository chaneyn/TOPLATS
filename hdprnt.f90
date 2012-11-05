! ====================================================================
!
!			subroutine hdrprnt
!
! ====================================================================
!
! Subroutine to print the output file headings.
!
! ====================================================================

      subroutine hdprnt(iprn)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      integer iprn(MAX_FIL)

! ====================================================================
! Write output file headers.
! ====================================================================

      if (iprn(75).eq.1) then

         write(75,*)'  Catchment Table'
         write(75,*)
         write(75,*)&
       '                     Average              Initial            Qo'
         write(75,*)'Basin              Topographic              W.T.'
         write(75,*)' No.       Area       Index      ln Te     Depth'
         write(75,*)&
           '(km^2)                            (m)      (m-1) (m^3/s)'
      endif


      if (iprn(76).eq.1) then

         write(76,*)'  Land Cover Table'
         write(76,*)
         write(76,*)'Land    Leaf               Canopy ' 
         write(76,*)'Cover   Area     Canopy    Storage'
         write(76,*)' No.    Index    Resist   Capacity.'
         write(76,*)'                 (s/m)       (m)'

      endif

      if (iprn(77).eq.1) then

         write(77,*)'  Land Cover Fractions by Catchment'
         write(77,*)

      endif

      if (iprn(78).eq.1) then

         write(78,*)'  Soil Type Table'
         write(78,*)
         write(78,*)&
        'Soil Pore Size           Saturated  Residual Saturated' 
         write(78,*)&
        'Type   Dist.    Bubbling    Soil      Soil   Hydraulic'
         write(78,*)&
        ' No.   Index    Pressure  Moisture  Moisture    Cond'
         write(78,*)'                  (cm)                         (mm/h)'

      endif      

      if (iprn(79).eq.1) then

         write(79,*)'  Soil Type Fractions by Catchment'
         write(79,*)

      endif

      if (iprn(80).eq.1) then

         write(80,*)' Catchment Evaporation Results'
         write(80,*)
         write(80,*)&
                '                       Bare       Dry       Wet    Frac   Frac'
         write(80,*)&
                'Time Catch   Total     Soil     Canopy    Canopy   Wet    Bare'
         write(80,*)&
                'Step  No.    Evap      Evap      Evap      Evap   Canopy  Soil'
         write(80,*)'            (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(81).eq.1) then

         write(81,*)' Catchment Precipitation/Infiltration Results'
         write(81,*)
         write(81,*)'                                                                Sat       Inf'
         write(81,*)'Time Catch   Total      Net                          Total    Excess    Excess'
         write(81,*)'Step  No.   Precip    Precip    Condens   Runoff    Infilt    Infilt    Infilt'
         write(81,*)'            (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(82).eq.1) then
! NWC 13/06/11
!         write(82,*)' Catchment Water Table Balance Results'
!         write(82,*)
!         write(82,*)'            New    Old                                            RZ     TZ'
!         write(82,*)'Time Catch  Ave    Ave   Capillary Drainage    Evap              Avail  Avail'
!         write(82,*)'Step  No.   W.T.   W.T.    Rise     to W.T.  from WT   Baseflow  Pore   Pore'
!         write(82,*)
!                '            (m)    (m)    (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(90).eq.1) then

         write(90,*)' Regional Energy Fluxes at PET'
         write(90,*)
         write(90,*)' Time    Net     Latent   Sensible   Ground     Energy   G + DS   Surface   Surface   Surface'
         write(90,*)' Step Radiation   Heat      Heat      Heat     Balance              Temp     Temp      Temp'
         write(90,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)    (K)       (K)       (K)'

      endif
         
      if (iprn(91).eq.1) then

         write(91,*)' Regional Actual Energy Fluxes'
         write(91,*)
         write(91,*)' Time    Net     Latent   Sensible   Ground     Energy   G + DS   Surface   Surface   Surface'
         write(91,*)' Step Radiation   Heat      Heat      Heat     Balance              Temp     Temp      Temp'
         write(91,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)    (K)       (K)       (K)'
!cdp      write(91,*)' Time    Net     Latent   Sensible   Ground     Energy   Surface'
!cdp      write(91,*)' Step Radiation   Heat      Heat      Heat     Balance    Temp'
!cdp      write(91,*)'       (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)   (W/m^2)'

      endif

      if (iprn(92).eq.1) then

         write(92,*)' Regional Canopy Water Balance'
         write(92,*)
         write(92,*)&
              '         New       Old                           Wet    Fraction'
         write(92,*)' Time   Canopy    Canopy   Total      Net      Canopy      Wet    Change in  Sum of     Water'
         write(92,*)' Step  Storage   Storage  Rainfall  Rainfall    Evap     Canopy    Storage   Fluxes    Balance'
         write(92,*)'         (mm)      (mm)    (mm/h)    (mm/h)    (mm/h)                (m)       (m)       (m)'

      endif
         
      if (iprn(93).eq.1) then

         write(93,*)' Regional Precipitation/Infiltration/Runoff'
         write(93,*)
         write(93,*)'                                                       Saturation  Infilt'
         write(93,*)' Time   Total      Net               Surface    Total    Excess    Excess'
         write(93,*)' Step Rainfall  Rainfall   Condens   Runoff    Infilt    Runoff    Runoff'
         write(93,*)'       (mm/h)    (mm/h)     (mm/h)   (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif

      if (iprn(94).eq.1) then

         write(94,*)' Regional Evapotranspiration Rates'
         write(94,*)
         write(94,*)'                  Bare       Dry       Wet     Frac   Frac'
         write(94,*)' Time   Total     Soil     Canopy    Canopy    Wet    Bare'
         write(94,*)' Step   Evap      Evap      Trans     Evap    Canopy  Soil'
         write(94,*)'       (mm/h)    (mm/h)    (mm/h)    (mm/h)'

      endif
     
      if (iprn(95).eq.1) then

         write(95,*)' Regional Root and Transmission Zone Water Balance'
         write(95,*)
         write(95,*)&
              '            Root Zone                      Transmision Zone'
         write(95,*)&
               ' Time  Soil  Change in  Sum of         Soil  Change in   Sum of'
         write(95,*)&
               ' Step Moist   Storage   Fluxes        Moist   Storage    Fluxes'
         write(95,*)&
               '               (mm)      (mm)                   (mm)       (mm)'

      endif
         
      if (iprn(96).eq.1) then

! NWC 06/13/11
!         write(96,*)' Regional Water Table Balance and Vertical Fluxes'
!         write(96,*)
!         write(96,*)'        New       Old'
!         write(96,*)' Time   Avg       Avg     Capillary       Drainage       Evap
!                       Drainage      Drainage      Diffusion  Diffusion  '
!         write(96,*)' Step  Depth     Depth    Rise fr WT       to WT        from WT
!              Baseflow  fr RZ        fr TZ          to RZ       to TZ   '
!         write(96,*)'        (mm)      (mm)     (mm/h)         (mm/h)         (mm/h)   
!             (mm/h)     (mm/h)      (mm/h)         (mm/h)      (mm/h)  '

      endif

      if (iprn(97).eq.1) then

         write(97,*)' Regional Fractional Saturation States'
         write(97,*)
         write(97,*)&
              ' Time  Region 3         Region 2                        Region 1 '
         write(97,*)' Step                  Saturated   Unsat         Saturated  TZ Sat  RZ Sat   Unsat'

      endif

      if (iprn(98).eq.1) then

         write(98,*)' Regional Evapotranspiration Controls and'
         write(98,*)'    Infiltration Mechanisms'
         write(98,*)
         write(98,*)&
               ' Time    Total  Atmos  Atmos   Land         Net    Sat   Infilt'
         write(98,*)&
                ' Step     Evap   Sat   Unsat  Unsat      Rainfall Exces  Exces'
         write(98,*)&
                              '         (mm/h)                           (mm/h)'

      endif

      if (iprn(110).eq.1) then

         write (110,*) 'Results for the moss layer, actual fluxes '
         write (110,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin      Ttopsoil        Tmid   theta_mois'

      endif

      if (iprn(111).eq.1) then

         write (111,*) 'Results for the understory layer, actual fluxes'
         write (111,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin        Tmid'

      endif

      if (iprn(112).eq.1) then

         write (112,*) 'Results for the overstory layer, actual fluxes'
          write (112,*)& 
       'Timestep Rnet      LE           H           G           ds        Tskin        Tmid'

      endif

      return

      end
