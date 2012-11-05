! ====================================================================
!
!            subroutine calcrsoil
!
! ====================================================================
!
! Calculate the bare soil resistance to evaporation.
!
! ====================================================================

      subroutine calcrsoil(irestype,rsoil,&
                           srespar1,srespar2,srespar3,&
                           thetas,rzsm,tkact)

      implicit none
      include "help/calcrsoil.h"

      if (irestype.eq.2) then

! ....................................................................
! Formulation of Sun (1982).
! ....................................................................

         rsoil = (srespar1*(((thetas)/(rzsm))**srespar2)) + srespar3

      endif

      if (irestype.eq.3) then

! ....................................................................
! Formulation of Kondo (1990).
! ....................................................................

         D0 = 0.229d-04
         Datm = D0 * (tkact/273.16)**srespar3
         rsoil = (srespar1* ((thetas - rzsm) **srespar2)) / Datm

      endif

      if (irestype.eq.4) then

! ....................................................................
! Formulation of Camillo and Gurney (1986).
! ....................................................................

         rsoil = srespar1*(thetas-rzsm) - srespar2

         if (rsoil.lt.(0.d0)) rsoil=0.d0

      endif

      if (irestype.eq.5) then

! ....................................................................
! Formulation of Passerat (1986).
! ....................................................................

         rsoil = srespar1* dexp(srespar2*rzsm/ srespar3)

      endif

      return

      end
