! ====================================================================
!
!			subroutine inittk
!
! ====================================================================
!
! Subroutine to assign initial soil temperatures
! and set heat storage variables for the first time step in the
! simulation.
!
! ====================================================================

      subroutine inittk(tdeep,tmid0,tmid0_moss,tkmid,&
       tkmid_us,tkmid_moss,tkel,tk0moss,tkact,tkact_us,&
       tkact_moss,tskinact_moss,dshact,dshact_us,dshact_moss,tkpet,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,dspet,dspet_us,dspet_moss,&
       Tsurf,Tpack,Tsurf_us,Tpack_us)

      implicit none
      include "help/inittk.h"

! ====================================================================
! Initialize average intermediate soil temperatures.
! ====================================================================

      ttemp=tdeep
      tkmidtmp=tmid0
      tkmidtmp_us=tmid0
      tkmidtmp_moss=tmid0_moss

! ====================================================================
! Initialize distributed intermediate soil temperatures.
! ====================================================================

      tkmid = tkmidtmp
      tkmid_us = tkmidtmp_us
      tkmid_moss = tkmidtmp_moss
      tkact = tkel
      tkact_us = tkel
      tkact_moss = tk0moss
      tskinact_moss = tkel

! ====================================================================
! Initialize storages.
! ====================================================================

      dshact= zero
      dshact_us= zero
      dshact_moss= zero

! ====================================================================
! Initialize skin temperatures.
! ====================================================================

      tkpet=ttemp
      tkmidpet = tkmidtmp
      tkmidpet_us = tkmidtmp_us
      tkmidpet_moss = tkmidtmp_moss

! ====================================================================
! Initialize storages.
! ====================================================================

      dspet=zero
      dspet_us=zero
      dspet_moss=zero

! ====================================================================
! Initialize snow temperatures.
! ====================================================================

      Tsurf=tkel-273.15d0
      Tpack=tkel-273.15d0
      Tsurf=-10.d0
      Tpack=-10.d0
      Tsurf_us=tkel-273.15d0
      Tpack_us=tkel-273.15d0

      return

      end
