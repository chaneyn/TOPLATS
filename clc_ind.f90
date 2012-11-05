! ====================================================================
!
!             subroutine clc_ind
!
! ====================================================================
!
! Calculate the matrix indexes.
!
! ====================================================================

      subroutine clc_ind (ilandc,ipix,SNOW_RUN,sw_lc,sw_px,SNW_FLG,s_lc,s_px)

      implicit none

      include 'help/clc_ind.h'

! --------------------------------------------------------------------
! Snow indexes, overstory.
! --------------------------------------------------------------------

      if (SNOW_RUN.eq.1) then

         sw_lc=ilandc
         sw_px=ipix

      endif

      if (SNOW_RUN.eq.0) then

         sw_lc=1
         sw_px=1

      endif

! --------------------------------------------------------------------
! Snow indexes, understory.
! --------------------------------------------------------------------

      if (SNW_FLG.eq.1) then

         s_lc=ilandc
         s_px=ipix

      endif

      if (SNW_FLG.eq.0) then

         s_lc=1
         s_px=1

      endif

      return

      end
