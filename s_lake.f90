! ====================================================================
!
!                        subroutine s_lake
!
! ====================================================================
!
! Change the lake variables of the statistical and distributed modes
! with each other.
!
! ====================================================================

      subroutine s_lake(ii,int,landc,ipix,s_preca_a,preca_a,&
       s_tp_in,tp_in,s_hice_in,hice_in,s_hsnw_in,hsnw_in,s_temp_a,temp_a,&
       s_tempi_a,tempi_a,s_hice_a,hice_a,s_hsnow_a,hsnow_a,s_mixmax,mixmax,&
       s_xlat_a,xlat_a,s_xlong_a,xlong_a,s_eta_a,eta_a,s_numnod,numnod,&
       s_surface,surface,s_fraci_a,fraci_a,s_precacc,precacc)

      implicit none
      include "SNOW.h"
      include "LAKE.h"
      include "wgtpar.h"
      include "help/s_lake.h"

      if (LAK_FLG.eq.1) then

         s_lc=landc
         s_px=ipix
         nodmax=MAX_NOD

      endif

      if (LAK_FLG.eq.0) then

         s_lc=1
         s_px=1
         nodmax=1

      endif

      if (ii.eq.0) then

         preca_a(s_px)=s_preca_a(int,s_lc,s_px)
         tp_in(s_px)=s_tp_in(int,s_lc,s_px)
         hice_in(s_px)=s_hice_in(int,s_lc,s_px)
         hsnw_in(s_px)=s_hsnw_in(int,s_lc,s_px)
         tempi_a(s_px)=s_tempi_a(int,s_lc,s_px)
         hice_a(s_px)=s_hice_a(int,s_lc,s_px)
         hsnow_a(s_px)=s_hsnow_a(int,s_lc,s_px)
         mixmax(s_px)=s_mixmax(int,s_lc,s_px)
         xlat_a(s_px)=s_xlat_a(int,s_lc,s_px)
         xlong_a(s_px)=s_xlong_a(int,s_lc,s_px)
         eta_a(s_px)=s_eta_a(int,s_lc,s_px)
         numnod(s_px)=s_numnod(int,s_lc,s_px)
         fraci_a(s_px)=s_fraci_a(int,s_lc,s_px)
         precacc(s_px)=s_precacc(int,s_lc,s_px)

         do j=1,nodmax

            surface(s_px,j)=s_surface(int,s_lc,s_px,j)
            temp_a(s_px,j)=s_temp_a(int,s_lc,s_px,j)

         enddo

      endif


      if (ii.eq.1) then

         s_preca_a(int,s_lc,s_px)=preca_a(s_px)
         s_tp_in(int,s_lc,s_px)=tp_in(s_px)
         s_hice_in(int,s_lc,s_px)=hice_in(s_px)
         s_hsnw_in(int,s_lc,s_px)=hsnw_in(s_px)
         s_tempi_a(int,s_lc,s_px)=tempi_a(s_px)
         s_hice_a(int,s_lc,s_px)=hice_a(s_px)
         s_hsnow_a(int,s_lc,s_px)=hsnow_a(s_px)
         s_mixmax(int,s_lc,s_px)=mixmax(s_px)
         s_xlat_a(int,s_lc,s_px)=xlat_a(s_px)
         s_xlong_a(int,s_lc,s_px)=xlong_a(s_px)
         s_eta_a(int,s_lc,s_px)=eta_a(s_px)
         s_numnod(int,s_lc,s_px)=numnod(s_px)
         s_fraci_a(int,s_lc,s_px)=fraci_a(s_px)
         s_precacc(int,s_lc,s_px)=precacc(s_px)

         do j=1,nodmax

            s_surface(int,s_lc,s_px,j)=surface(s_px,j)
            s_temp_a(int,s_lc,s_px,j)=temp_a(s_px,j)

         enddo

      endif

      return

      end
