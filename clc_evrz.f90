! ====================================================================
!
!                  subroutine clc_evrz
!
! ====================================================================
!
! Calculate the evaporation/condensation into the root zone.
!
! ====================================================================

      subroutine clc_evrz(evrz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
       i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,dummy,f_und)

      implicit none
      include "help/clc_evrz.h"

      evrz=zero

      if ((Swq+Swq_us).le.(0.d0)) then

! ....................................................................
! If the vegetation is bare soil and there is no snow layer present
! all the evaporative demand comes from the bare soil and under story
! or moss is not represented.
! ....................................................................

         if (ivgtyp.eq.0) evrz=evtact*dc

         if (ivgtyp.eq.1) then

            if ( (i_und.eq.0).and.&
                 (i_moss.eq.0) ) then

! ....................................................................
! In case of vegetation with lower roots all the evaporative demand
! for the soil comes from the over story if there is no under story
! represented and if there is no snow.
! ....................................................................

               evrz=evtact*dc*(1.d0-fw)
               dummy=evrz

            else

               if (i_und.gt.0) then

! ....................................................................
! If there is under story and no snow part of the evporative demand
! comes from the over story and part from the under story if under
! story is represented.
! ....................................................................

                  evrz=(evtact*dc*(1.d0-fw))+ &
                       f_und*(evtact_us*dc_us*(1.d0-fw_us))

                  dummy=evrz

               endif

               if (i_moss.gt.0) then

! ....................................................................
! If there is a moss layer and no snow then all the evaporative demand
! for the soil comes from the over story layer.  This is justified
! by the fact that roots under a moss layer  re sufficiently
! vertically distributed to make this assumption hold.
! ....................................................................

                  evrz=(evtact*dc*(1.d0-fw))
                  dummy=evrz

               endif

            endif

         endif

         if (ivgtyp.eq.2) then

            if ( (i_und.eq.0).and.&
                 (i_moss.eq.0) ) then

! ....................................................................
! In case of lower layer roots and no under story
! than there is no evaporative demand for the root zone.
! ....................................................................

               evrz=zero
               dummy=evrz

            else

               if (i_und.gt.0) then

! ....................................................................
! In case of lower layer roots and an under story layer than
! the evaporative demand for the root zone comes entirely from
! the under story.
! ....................................................................

                  evrz=zero+&
                       f_und*(evtact_us*dc_us*(1.d0-fw_us))
                  dummy=evrz

               endif

               if (i_moss.gt.0) then

! ....................................................................
! In case of lower layer roots and a moss layer
! than there is no evaporative demand for the root zone.
! ....................................................................

                  evrz=zero
                  dummy=evrz

               endif

            endif

         endif

      else

! ....................................................................
! In case of snow on top of the under story all the evaporative demand
! comes from the over story.
! ....................................................................

         evrz=zero

         if (Swq.le.0.d0) then

            if (ivgtyp.eq.1) then

               evrz=evrz+evtact*dc*(1.d0-fw)*(1.-f_und)

            endif

            if (ivgtyp.eq.0) evrz=evtact*dc

         endif
! ....................................................................
! In case of snow on top of the over story all the evaporative demand
! comes from the under story.
! ....................................................................

         if (Swq_us.le.0.d0) then

            if (i_und.gt.0) then

               evrz=evrz+evtact_us*dc_us*(1.d0-fw_us)*f_und

            endif

            if (i_moss.gt.0) then

               evrz=evrz+zero

            endif

         endif

         dummy=evrz

      endif

! ....................................................................
! Water for evaporation of the over story is for 50 % extracted
! from soil water and 50 % from water in the moss layer.
! ....................................................................

      if (i_moss.eq.1) then

         evrz_moss=0.5d0*evrz
         evrz=0.5d0*evrz

      else

         evrz_moss=0.d0

      endif

      return

      end
