! ====================================================================
!
!			subroutine initsm
!
! ====================================================================
!
! Subroutine to assign spatially variable initial conditions in
! the root and transmission zones for the first time step in the
! storm.
!
! ====================================================================

      subroutine initsm(zw,psic,zrz,ztz,rzsm1,tzsm1,thetas,&
       zrzmax,iopsmini,thetar,bcbeta,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,&
       inc_frozen,tsoilold,bulk_dens,a_ice,b_ice,row)

      implicit none
      include "help/initsm.h"

! ====================================================================
! Update root and transmission zone depths and soil moisture.
! ====================================================================

      if ((zw-psic).le.zero) then

! --------------------------------------------------------------------
! If the root zone is saturated (Region 3).
! --------------------------------------------------------------------

         zrz = zero
         ztz = zero
         rzsm1 = thetas
         tzsm1 = thetas

      else if ((zw-psic).lt.zrzmax) then

! --------------------------------------------------------------------
! If the transmission zone is saturated and root zone is
! unsaturated (Region 2).
! --------------------------------------------------------------------

         zrz = zw-psic

         if (iopsmini.eq.0) then

            rzsm1 = thetar+(thetas-thetar)*((psic/zw)**bcbeta)

         endif

         ztz = zero
         tzsm1 = thetas

      else

! --------------------------------------------------------------------
! If the transmission and root zone are both unsaturated (Region 1).
! --------------------------------------------------------------------

         zrz = zrzmax
         ztz = zw-psic-zrz
 
         if (iopsmini.eq.0) then

            rzsm1 = thetar+(thetas-thetar)*((psic/(zw-0.5*zrz))**bcbeta)
            tzsm1 = thetar+(thetas-thetar)*((psic/(0.5*ztz+psic))** bcbeta)

         endif

      endif

! ====================================================================
! Calculate the frozen water and liquid water fractions in the
! soil if requested.
! ====================================================================

      rzsm1_u=0.d0
      tzsm1_u=0.d0
      rzsm1_f=0.d0
      tzsm1_f=0.d0

      if (inc_frozen.eq.1) then

         if (tsoilold.gt.(273.15d0)) then

! --------------------------------------------------------------------
! If the soil is not frozen all soil water is liquid.
! --------------------------------------------------------------------

            rzsm1_u=rzsm1
            tzsm1_u=tzsm1
            rzsm1_f=0.d0
            tzsm1_f=0.d0

         else

! --------------------------------------------------------------------
! In case of frozen soil calculate the unfrozen soil water.
! --------------------------------------------------------------------

            ttt=273.15d0-tsoilold
            rzsm1_u=1000.d0*bulk_dens*a_ice*ttt**b_ice/row
            tzsm1_u=1000.d0*bulk_dens*a_ice*ttt**b_ice/row

            if (rzsm1_u.gt.rzsm1) then

! ....................................................................
! Check whether the unfrozen soil water in the root zone
! is not higher than the root zone soil moisture.  Adapt
! if necessary.
! ....................................................................

               rzsm1_u=rzsm1
               rzsm1_f=0.d0

            else

! ....................................................................
! Frozen water content = total water content minus unfrozen water
! content.
! ....................................................................

               rzsm1_f=rzsm1-rzsm1_u

            endif

! ....................................................................
! Check whether the unfrozen soil water in the transmission zone
! is not higher than the transmission zone soil moisture.  Adapt
! if necessary.
! ....................................................................

            if (tzsm1_u.gt.tzsm1) then

               tzsm1_u=tzsm1
               tzsm1_f=0.d0

            else

! ....................................................................
! Frozen water content = total water content minus unfrozen water
! content.
! ....................................................................

               tzsm1_f=tzsm1-tzsm1_u

            endif

         endif

! --------------------------------------------------------------------
! Frozen water content = total water content minus unfrozen water
! content.  (This, in some cases, has already been calculated
! here above but doint this calculation again dos not affect the
! results.
! --------------------------------------------------------------------

         rzsm1_u=rzsm1-rzsm1_f
         tzsm1_u=tzsm1-tzsm1_f

         if ( (thetas-rzsm1_f).le.(thetar+0.d0)) then

! ....................................................................
! Check whether the frozen water content in the root zone is not
! higher than the saturated water content.
! ....................................................................

            rtdif=0.0d0+thetar-(thetas-rzsm1_f)
            rzsm1_f=rzsm1_f-rtdif
            rzsm1_u=rzsm1_u+rtdif

         endif

         if ( (thetas-tzsm1_f).le.(thetar)+0.d0) then

! ....................................................................
! Check whether the frozen water content in the transmission zone is not
! higher than the saturated water content.
! ....................................................................

            rtdif=0.0d0+thetar-(thetas-tzsm1_f)
            tzsm1_f=tzsm1_f-rtdif
            tzsm1_u=tzsm1_u+rtdif

         endif

         if (rzsm1_u.le.(thetar+0.d0)) then

! ....................................................................
! Check whether the unfrozen water content in the root zone is not
! lower than the residual water content.
! ....................................................................

            rtdif=thetar+0.01d0-rzsm1_u
            rzsm1_u=rzsm1_u+rtdif
            rzsm1_f=rzsm1_f-rtdif

         endif

         if (tzsm1_u.le.(thetar+0.d0)) then

! ....................................................................
! Check whether the unfrozen water content in the transmission zone is
! not lower than the residual water content.
! ....................................................................

            rtdif=thetar+0.01d0-tzsm1_u
            tzsm1_u=tzsm1_u+rtdif
            tzsm1_f=tzsm1_f-rtdif

         endif

      endif

      return

      end
