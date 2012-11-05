! ====================================================================
!
!                       subroutine calc_rs
!
! ====================================================================
!
! Calculate the incoming solar radiation for the under and over story
! under the assumption of only one reflection.
!
! ====================================================================

      subroutine calc_rs(canclos,extinct,i_und,i_moss,Swq_us,&
                         albd_us,alb_moss,alb_snow,rsd,rs_over,rs_under)

      implicit none
      include "help/calc_rs.h"

! ====================================================================
! Calculate the incoming solar radiation for the under story and over
! story layers.
! ====================================================================

      refus=0.d0
      ccc=canclos
      thr=extinct

! --------------------------------------------------------------------
! Determine what albedo of the understory is : moss, snow or normal
! vegetation.
! --------------------------------------------------------------------

      if (i_und.gt.0) refus=albd_us
      if (i_moss.gt.0) refus=alb_moss
      if (Swq_us.gt.(0.d0)) refus=alb_snow

! --------------------------------------------------------------------
! Calculate the incoming radiation under the assumption of only
! one reflection.
! --------------------------------------------------------------------

      if ( (i_und.gt.0).or.(i_moss.gt.0) ) then

         rs_over=(1.d0+refus*thr)*rsd
         rs_under=rsd*(thr*ccc+1.d0-ccc)

      endif

      if ( (i_und.eq.0).and.(i_moss.eq.0) ) then

         rs_over=rsd
         rs_under=rsd

      endif

      return

      end
