! ====================================================================
!
!			subroutine transv
!
! ====================================================================
!
! Subroutine to find actual transpiration from dry canopy.
!
! ====================================================================

      subroutine transv(epetd,epetd_us,i_und,iopveg,f1par,f3vpd,f4temp,&
       rescan,inc_frozen,ivgtyp,rzsm,rzsm_u,tc,tw,smcond,tzsm,tzsm_u,&
       tc_us,tw_us,smcond_us,f1par_us,f3vpd_us,f4temp_us,rescan_us,&
       vegcap,ravd,vegcap_us,ravd_us,zrz,srzrel,thetas,thetar,psisoi,&
       psic,bcbeta,ikopt,xksrz,xk0,ff,ressoi,rtact,rtdens,psicri,&
       respla,xkrz,ztz,stzrel,xkstz,xktz,Swq,evtact,ievcon,Swq_us,evtact_us,&
       ievcon_us,bsdew,i,ipix)

      implicit none
      include "help/transv.h"
      integer :: i,ipix

! ====================================================================
! Set potential transpiration to zero if negative for both under s.or.&
! and over story.
! ====================================================================

      if (epetd.le.zero) epetd=zero
      if ( (i_und.gt.0) .and. (epetd_us.le.zero) ) epetd_us=zero

! ====================================================================
! Calculate the maximum water vapor flux out of the plants possible
! by using a soil moisture conductance parameterization.
! ====================================================================

      if (iopveg.eq.0) then

! --------------------------------------------------------------------&
! Set the multiplier for the linear interpolation between 
! critical and wilting soil moisture condition based on 
! soil moisture in either the upper or lower layer.
! --------------------------------------------------------------------&

         resist = f1par*f3vpd*f4temp*rescan

         if (inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! Assuming no difference between ice particles and unfrozen water.
! --------------------------------------------------------------------&

            if (ivgtyp.eq.1) then

! ....................................................................
! Vegetation roots are in the upper layer.  Linear increase with
! soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(rzsm,tc,smcond,one,tw,zero)

            else

! ....................................................................
! Vegetation roots are in the lower layer.  Linear increase with
! soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(tzsm,tc,smcond,one,tw,zero)

            endif

         else

! --------------------------------------------------------------------&
! Treating the frozen soil water as soil particles.
! --------------------------------------------------------------------&

            if (ivgtyp.eq.1) then

! ....................................................................
! Vegetation roots are in the upper layer.  Linear increase with
! liquid soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(rzsm_u,tc,smcond,one,tw,zero)

            else

! ....................................................................
! Vegetation roots are in the lower layer.  Linear increase with
! liquid soil moisture between wilting and critical soil moisture, 1 if
! the soil moisture is higher than the wilting soil moisture and 0
! if the soil moisture is lower than the critical soil moisture.
! ....................................................................

               call calcsmcond(tzsm_u,tc,smcond,one,tw,zero)

            endif

         endif

! ====================================================================
! Calculate the soil moisture effect on stomatal resistance for
! the under story.  The roots of the under story are always assumed
! to be in the upper soil layer.
! ====================================================================

         if (inc_frozen.eq.0) then

! --------------------------------------------------------------------&
! Assuming no difference between ice particles and unfrozen water.
! --------------------------------------------------------------------&

            if (i_und.gt.0) then

! ....................................................................
! Linear increase with liquid soil moisture between wilting and critical
! soil moisture, 1 if the soil moisture is higher than the wilting soil
! moisture and 0 if the soil moisture is lower than the critical soil
! moisture.
! ....................................................................

               call calcsmcond(rzsm,tc_us,smcond_us,one,tw_us,zero)

               resist_us=f1par_us*f3vpd_us*f4temp_us*rescan_us

           endif

        else

! --------------------------------------------------------------------&
! Treating the frozen soil water as soil particles.
! --------------------------------------------------------------------&

           if (i_und.gt.0) then

! ....................................................................
! Linear increase with liquid soil moisture between wilting and critical
! soil moisture, 1 if the soil moisture is higher than the wilting soil
! moisture and 0 if the soil moisture is lower than the critical soil
! moisture.
! ....................................................................

               call calcsmcond(rzsm_u,tc_us,smcond_us,one,tw_us,zero)

               resist_us=f1par_us*f3vpd_us*f4temp_us*rescan_us

            endif

         endif

! ====================================================================
! Calculate vegetation capacity for the over story using linear
! interpolation for the canopy resistance.
! ====================================================================

         call calcvegcap(smcond,zero,vegcap,epetd,resist,ravd,smcond)

! ====================================================================
! Calculate vegetation capacity for the under story using linear
! interpolation for the canopy resistance.
! ====================================================================

         if (i_und.gt.0) then

            call calcvegcap(smcond_us,zero,vegcap_us,&
                            epetd_us,resist_us,ravd_us,smcond_us)

         endif

! ====================================================================
! Calculate the maximum water vapor flux out of the plants possible
! by using a root resistivity parameterization.
! ====================================================================

      else

! ====================================================================
! Calculate the maximum plant evaporation for the under story.
! ====================================================================

         if (i_und.gt.0) then

            call maxplevap(zrz,0.d0,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap_us,psicri,respla,xkrz)

         endif

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         if (ivgtyp.eq.1) then

            call maxplevap(zrz,0.d0,epetd,inc_frozen,srzrel,rzsm,thetas,&
       thetar,rzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xksrz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla,xkrz)

! ====================================================================
! Calculate the maximum plant evaporation for the over story for
! vegetation with its roots in the upper soil layer.
! ====================================================================

         else

            call maxplevap(ztz,zrz,epetd,inc_frozen,stzrel,tzsm,thetas,&
       thetar,tzsm_u,zero,one,psisoi,psic,bcbeta,ikopt,xkstz,xk0,ff,&
       two,three,ressoi,rtact,rtdens,vegcap,psicri,respla,xktz)

         endif

      endif

! ====================================================================
! Set actual transpiration to the minimum of potential 
! transpiration or vegetation capacity.
! ====================================================================

      call acttrans(Swq,vegcap,epetd,evtact,ievcon,zrz)

      if (i_und.gt.0) then

        call acttrans(Swq_us,vegcap_us,epetd_us,evtact_us,ievcon_us,zrz)

      endif 
! ====================================================================
! Since condensation occurs in canopy, set bare soil condensation
! to zero.
! ====================================================================

      bsdew = zero

      return

      end
