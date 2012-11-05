! ====================================================================
!
!			subroutine states
!
! ====================================================================
!
! Subroutine to update water table depths and root and transmission
! zone soil moisture to the conditions at the end of the previous
! time step.
!
! ====================================================================

      subroutine states(zw0,inc_frozen,i_moss,tkmid_moss,&
       r_mossm_u,r_mossm_f,r_mossm,zw,zbar,ff,atanb,xlamda,psic,&
       zrz,ztz,rzsm1,tzsm1,thetas,zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,&
       tzsm1_u,rzsm1_f,tzsm1_f,tsoilold,bulk_dens,a_ice,b_ice,&
       row,rzsmold,tzsmold,r_mossmold,rzsm,tzsm,r_mossm1,zmoss,r_moss_depth,&
       thetas_moss,rzsm_u,rzsm_f,tzsm_u,tzsm_f,r_mossm1_u,&
       r_mossm1_f,i,a_ice_moss,b_ice_moss,bulk_dens_moss)

      implicit none
      include "help/states.h"
    
! ====================================================================
! Calculate the frozen moss content.
! ====================================================================

      if (inc_frozen.eq.1) then

         if (i_moss.eq.1) then

            if (tkmid_moss.lt.(273.15d0)) then

               r_mossm_u=0.d0
               r_mossm_u=1000.d0*bulk_dens_moss*a_ice_moss*(273.15d0-tkmid_moss)**b_ice_moss/row

               if (r_mossm_u.ge.r_mossm) then

                  r_mossm_u=r_mossm

               endif

               r_mossm_f=r_mossm-r_mossm_u

            endif

            if (tkmid_moss.ge.(273.15d0)) then

               r_mossm_f=0.d0
               r_mossm_u=r_mossm

            endif

         endif

      endif

! ====================================================================
! Update local water table depth.
! ====================================================================

      zw0 = zw - psic
      zw = zbar-(one/ff)*(atanb-xlamda)
      
!cw! minimum size for root zone and transmission zone is 1 cm

      if (zw-psic.lt.(zrzmax+0.01).and.zw-psic.gt.zrzmax) then

         zw=zrzmax + psic - .00001

      endif
	 
      if (zw-psic.lt.(0.01).and.zw-psic.gt.zero) then

         zw=psic

      endif
!cw!
! ====================================================================
! First time step call initsm to get spatially variable initial
! conditions based on brooks-corey or user specified values.
! ====================================================================

      if (i.eq.1) then

         call initsm(zw,psic,zrz,ztz,rzsm1,tzsm1,thetas,&
       zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,&
       inc_frozen,tsoilold,bulk_dens,a_ice,b_ice,row)

         rzsm_u=rzsm1_u
         tzsm_u=tzsm1_u
         rzsm_f=rzsm1_f
         tzsm_f=tzsm1_f

      endif

      if (inc_frozen.eq.0) then

! ====================================================================
! Update the values for the old and new timestep for in the option
! that the frozen soil water is treated as liquid water.
! ====================================================================

         rzsmold = rzsm
         tzsmold = tzsm
         r_mossmold = r_mossm

! --------------------------------------------------------------------&
! Update root and transmission zone depths and soil moisture.
! Also update upper and lower soil layer soil moistures.
! --------------------------------------------------------------------&

         if ((zw-psic).le.zero) then

! ....................................................................
! For saturated areas.
! ....................................................................

            zrz = zero
            rzsm = thetas
            ztz = zero
            tzsm = thetas


         else if ((zw-psic).lt.zrzmax) then

! ....................................................................
! For not saturated areas where the ground water level reaches
! the root zone.
! ....................................................................

            zrz = zw-psic
            rzsm = rzsm1
            ztz = zero
            tzsm = thetas

         else

! ....................................................................
! For not saturated areas where the ground water level is in the
! transmission zone.
! ....................................................................

            zrz = zrzmax
            ztz = zw-psic-zrz
            rzsm = rzsm1
            tzsm = tzsm1

         endif

     
! ....................................................................
! Update the values for moss moisture content of the old and new
! timestep.
! ....................................................................

         r_mossm1=r_mossm

         if (thetas_moss.gt.0.d0) then

            zmoss=r_moss_depth*r_mossm/thetas_moss

         else

            zmoss=0.d0

         endif

      else

! ====================================================================
! Update the values for the old and new timestep for in the option
! that the frozen soil water is treated as solid soil particles.
! ====================================================================

         rzsmold = rzsm
         tzsmold = tzsm
         r_mossmold = r_mossm

! --------------------------------------------------------------------&
! Update root and transmission zone depths and soil moisture.
! Also update upper and lower soil layer soil moistures.
! --------------------------------------------------------------------&

         if ((zw-psic).le.zero) then

! ....................................................................
! For saturated areas.
! ....................................................................

            zrz = zero
            rzsm_u = thetas-rzsm_f
            ztz = zero
            tzsm_u = thetas-tzsm_f
            rzsm=rzsm_u+rzsm_f
            tzsm=tzsm_u+tzsm_f

         else if ((zw-psic).lt.zrzmax) then

! ....................................................................
! For not saturated areas where the ground water level reaches
! the root zone.
! ....................................................................

            zrz = zw-psic
            rzsm_u = rzsm1_u
            rzsm=rzsm_u+rzsm_f
            ztz = zero
            tzsm_u = thetas-tzsm_f
            tzsm=tzsm_u+tzsm_f

         else

! ....................................................................
! For not saturated areas where the ground water level is in the
! transmission zone.
! ....................................................................

            zrz = zrzmax
            ztz = zw-psic-zrz
            rzsm_u = rzsm1_u
            tzsm_u = tzsm1_u
            rzsm=rzsm_u+rzsm_f
            tzsm=tzsm_u+tzsm_f

         endif

      
  
! ....................................................................
! Update the values for moss moisture content of the old and new
! timestep.
! ....................................................................

         r_mossm1_u=r_mossm_u
         r_mossm=r_mossm_u+r_mossm_f
!CVAL         zmoss=r_moss_depth*r_mossm/(thetas_moss-r_mossm_f)

         if (thetas_moss.gt.0.d0) then

            zmoss=r_moss_depth*r_mossm/(thetas_moss)

         else

            zmoss=0.d0

         endif

         r_mossm1_f=r_mossm_f

      endif

!CVAL      write (129,*) real(r_mossm_f/r_mossm),real(rzsm_f/rzsm),real(tzsm_f/tzsm)
!CVAL      write (130,*) real(zw-psic)

      return

      end
