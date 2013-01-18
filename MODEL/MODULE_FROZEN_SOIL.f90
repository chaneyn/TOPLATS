MODULE MODULE_FROZEN_SOIL

USE MODULE_VARIABLES

contains

! ====================================================================
!
!                       subroutine ice_change
!
! ====================================================================
!
! Calculate the change in ice content in the soil over time using
! Osterkamp (1987).
!
! ====================================================================

      subroutine ice_change(ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold)

      implicit none
      integer ipix
      real*8 rzdthetaidt,tzdthetaidt,f_moss,f_und,tkmidpet,tkmidpet_us
      real*8 tkmidpet_moss,rzsm1_f,tzsm1_f,bulk_dens,a_ice,b_ice
      real*8 row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,thetas,thetar
      real*8 rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold
      real*8 ttt,zsm,tsm,rtdif
      integer wrong1,wrong2

! --------------------------------------------------------------------
! Check the soil moistures fed into the program.
! --------------------------------------------------------------------

      wrong1=0
      wrong2=0

      call check_in_ice (rzsm1,1,wrong1,ipix)
      call check_in_ice (rzsm_f,2,wrong1,ipix)
      call check_in_ice (rzsm1_f,3,wrong1,ipix)
      call check_in_ice (rzsm1_u,4,wrong1,ipix)
      call check_in_ice (tzsm1,5,wrong2,ipix)
      call check_in_ice (tzsm_f,6,wrong2,ipix)
      call check_in_ice (tzsm1_f,7,wrong2,ipix)
      call check_in_ice (tzsm1_u,8,wrong2,ipix)

      if (wrong1.eq.1) then

!         write (*,*) 'ICE_CHANGE : input root soil moistures out of bound '
!         write (*,*) rzsm1,rzsm1_u,rzsm1_f,rzsm_f
         call correct (thetas,thetar,rzsm1,rzsm1_u,rzsm1_f,rzsm_f)

      endif

      if (wrong2.eq.1) then

!         write (*,*) 'ICE_CHANGE : input trans soil moistures out of bound '
!         write (*,*) tzsm1,tzsm1_u,tzsm1_f,tzsm_f
         call correct (thetas,thetar,tzsm1,tzsm1_u,tzsm1_f,tzsm_f)

      endif
      wrong1=0
      wrong2=0

! --------------------------------------------------------------------
! Temperature used for the representation of the soil temperature is
! a weighed averge of the temperature under vegetation, understory
! or moss.
! --------------------------------------------------------------------

      ttt = tkmidpet

! --------------------------------------------------------------------
! If the soil temperature is above 0C then there is no frozen soil
! water.
! --------------------------------------------------------------------

      if (ttt.gt.(273.15d0)) then

         rzsm1_f=0.d0
         tzsm1_f=0.d0

      else

! --------------------------------------------------------------------
! Calculate the soil unfrozen water content using the formula from
! Osterkamp (1987).
! --------------------------------------------------------------------

         ttt=273.15d0-ttt
         zsm=1000.d0*bulk_dens*a_ice*ttt**b_ice/row
         tsm=1000.d0*bulk_dens*a_ice*ttt**b_ice/row

         if (zsm.gt.rzsm1) then

! --------------------------------------------------------------------
! If the unfrozen water content is higher than the soil moisture there
! is no frozen water in the root zone.
! --------------------------------------------------------------------

            rzsm1_f=0.d0

         else

! --------------------------------------------------------------------
! If the unfrozen water content is lower than the soil moisture the
! frozen water content is the soil moisture - the unfrozen water
! content in the root zone.
! --------------------------------------------------------------------

            rzsm1_f=rzsm1-zsm

         endif


         if (tsm.gt.tzsm1) then

! --------------------------------------------------------------------
! If the unfrozen water content is higher than the soil moisture there
! is no frozen water in the transmission zone.
! --------------------------------------------------------------------

            tzsm1_f=0.d0

         else

! --------------------------------------------------------------------
! If the unfrozen water content is lower than the soil moisture the
! frozen water content is the soil moisture - the unfrozen water
! content in the transmission zone.
! --------------------------------------------------------------------

            tzsm1_f=tzsm1-tsm

         endif

      endif

! --------------------------------------------------------------------
! Calculate the unfrozen water content.
! --------------------------------------------------------------------

      rzsm1_u=rzsm1-rzsm1_f
      tzsm1_u=tzsm1-tzsm1_f

      if ((thetas-rzsm1_f).le.&
             (thetar+0.d0)) then

! --------------------------------------------------------------------
! If the frozen water content in the root zone is higher than
! the saturated - residual soil moisture update the frozen and
! unfrozen water content.
! --------------------------------------------------------------------

         rtdif=0.0d0+thetar-(thetas-rzsm1_f)
         rzsm1_f=rzsm1_f-rtdif
         rzsm1_u=rzsm1_u+rtdif

      endif

! --------------------------------------------------------------------
! If the frozen water content in the transmission zone is higher than
! the saturated - residual soil moisture update the frozen and
! unfrozen water content.
! --------------------------------------------------------------------

      if ((thetas-tzsm1_f).le.&
             (thetar+0.d0)) then

         rtdif=0.0d0+thetar-(thetas-tzsm_f)
         tzsm1_f=tzsm1_f-rtdif
         tzsm1_u=tzsm1_u+rtdif

      endif

! --------------------------------------------------------------------
! If the unfrozen soil water content in the root zone is lower
! than the residual soil moisture uptdate both frozen and
! unfrozen soil moisture.
! --------------------------------------------------------------------

      if (rzsm1_u.le.(thetar+0.d0)) then

         rtdif=thetar+0.0d0-rzsm1_u
         rzsm1_u=rzsm1_u+rtdif
         rzsm1_f=rzsm1_f-rtdif

      endif

! --------------------------------------------------------------------
! If the unfrozen soil water content in the transmission zone is lower
! than the residual soil moisture uptdate both frozen and
! unfrozen soil moisture.
! --------------------------------------------------------------------

      if (tzsm1_u.le.(thetar+0.d0)) then

         rtdif=thetar+0.0d0-tzsm1_u
         tzsm1_u=tzsm1_u+rtdif
         tzsm1_f=tzsm1_f-rtdif

      endif

! --------------------------------------------------------------------
! Calculate the change in frozen water content over time.
! --------------------------------------------------------------------

      rzdthetaidt=(rzsm1_f-rzsm_f)/dt
      tzdthetaidt=(tzsm1_f-tzsm_f)/dt
      rzdthetaidt=rzdthetaidt*roi/row
      tzdthetaidt=tzdthetaidt*roi/row

! --------------------------------------------------------------------
! Calculate the change in unfrozen water content in the soil versus
! soil temperature.
! --------------------------------------------------------------------

      ttt = tkmidpet

      if (ttt.eq.tsoilold) then

         rzdthetaudtemp=0.d0

      endif

      if (ttt.ne.tsoilold) then

         rzdthetaudtemp=-(rzsm1_f-rzsm_f)/(ttt-tsoilold)

         if (rzdthetaudtemp.lt.0.d0) then

            if (rzdthetaudtemp.gt.-0.01d0) rzdthetaudtemp=0.01d0

         endif

         if (rzdthetaudtemp.gt.0.d0) then

            if (rzdthetaudtemp.lt.0.01d0) rzdthetaudtemp=0.01d0

         endif

      endif

! --------------------------------------------------------------------
! Check the newly calculated soil moistures.
! --------------------------------------------------------------------

      wrong1=0
      wrong2=0

      call check_ice (rzsm1,1,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong1)

      call check_ice (tzsm1,2,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong2)

      call check_ice (rzsm1_u,3,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong1)

      call check_ice (rzsm1_f,4,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong1)

      call check_ice (tzsm1_u,5,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong2)

      call check_ice (tzsm1_f,6,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong2)

      if (wrong1.eq.1) then

!         write (*,*) 'ICE_CHANGE : results root soil moistures out of bound '
!         write (*,*) rzsm1,rzsm1_u,rzsm1_f,rzsm_f
         call correct (thetas,thetar,rzsm1,rzsm1_u,rzsm1_f,rzsm_f)

      endif

      if (wrong2.eq.1) then

!         write (*,*) 'ICE_CHANGE : results trans soil moistures out of bound '
!         write (*,*) tzsm1,tzsm1_u,tzsm1_f,tzsm_f
         call correct (thetas,thetar,tzsm1,tzsm1_u,tzsm1_f,tzsm_f)

      endif

      return

      end subroutine ice_change

! ====================================================================
!
!                       subroutine correct
!
! ====================================================================
!
! Correct the frozen and unfrozen water contents.
!
! ====================================================================

      subroutine correct (ts,tr,tz1,tzu,tzf,tzof)

      implicit none

      real *8 ts,tr,tz1,tzu,tzf,tzof

      if (tz1.gt.ts) then

         tz1=ts
         tzu=ts
         tzf=0.d0
         tzof=0.d0

      endif

      if (tz1.lt.tr) then

         tz1=tr
         tzu=tr
         tzf=0.d0
         tzof=0.d0

      endif

      return

      end subroutine correct

! ====================================================================
!
!                       subroutine check_in_ice
!
! ====================================================================
!
! Check the input frozen and unfrozen water contents.
!
! ====================================================================

      subroutine check_in_ice (soilm,val,wrong,ipix)

      integer val,wrong,ipix

      real*8 soilm

      if ( (soilm.ge.0.d0).and.(soilm.le.1.d0) ) then

         soilm=soilm

      else

!         write (*,*) 'pixel : ',ipix
!         if (val.eq.1) write (*,*) 'ICE_CHANGE : rzsm1 out of bounds ',soilm
!         if (val.eq.2) write (*,*) 'ICE_CHANGE : rzsm_f out of bounds ',soilm
!         if (val.eq.3) write (*,*) 'ICE_CHANGE : rzsm1_f out of bounds ',soilm
!         if (val.eq.4) write (*,*) 'ICE_CHANGE : rzsm1_u out of bounds ',soilm
!         if (val.eq.5) write (*,*) 'ICE_CHANGE : tzsm1 out of bounds ',soilm
!         if (val.eq.6) write (*,*) 'ICE_CHANGE : tzsm_f out of bounds ',soilm
!         if (val.eq.7) write (*,*) 'ICE_CHANGE : tzsm1_f out of bounds ',soilm
!         if (val.eq.8) write (*,*) 'ICE_CHANGE : tzsm1_u out of bounds ',soilm
         wrong=1

      endif

      return

      end subroutine check_in_ice

! ====================================================================
!
!                       subroutine check_ice
!
! ====================================================================
!
! Check the calculated frozen and unfrozen water contents.
!
! ====================================================================

      subroutine check_ice (soilm,nbad,ttt,zsm,tsm,rtdif,&
       ipix,rzdthetaidt,tzdthetaidt,f_moss,f_und,&
       tkmidpet,tkmidpet_us,tkmidpet_moss,rzsm1_f,tzsm1_f,&
       bulk_dens,a_ice,b_ice,row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,&
       thetas,thetar,rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold,wrong)

      implicit none
      integer ipix
      real*8 rzdthetaidt,tzdthetaidt,f_moss,f_und,tkmidpet,tkmidpet_us
      real*8 tkmidpet_moss,rzsm1_f,tzsm1_f,bulk_dens,a_ice,b_ice
      real*8 row,roi,rzsm1_u,rzsm1,tzsm1_u,tzsm1,thetas,thetar
      real*8 rzdthetaudtemp,dt,rzsm_f,tzsm_f,tsoilold
      real*8 ttt,zsm,tsm,rtdif

      integer wrong,nbad
      real*8 soilm

      if ( (soilm.ge.0.d0).and.(soilm.le.1.d0) ) then

         soilm=soilm

      else

!         if (nbad.eq.1) write (*,*) 'ICE_CHANGE : rzsm1 out of bounds ',soilm
!         if (nbad.eq.2) write (*,*) 'ICE_CHANGE : tzsm1 out of bounds ',soilm
!         if (nbad.eq.3) write (*,*) 'ICE_CHANGE : rzsm1_u out of bounds ',soilm
!         if (nbad.eq.4) write (*,*) 'ICE_CHANGE : rzsm1_f out of bounds ',soilm
!         if (nbad.eq.5) write (*,*) 'ICE_CHANGE : tzsm1_u out of bounds ',soilm
!         if (nbad.eq.6) write (*,*) 'ICE_CHANGE : tzsm1_f out of bounds ',soilm
!         write (*,*) 'A',ipix,rzdthetaidt,tzdthetaidt
!         write (*,*) 'B',f_moss,f_und
!         write (*,*) 'C',tkmidpet,tkmidpet_us,tkmidpet_moss
!         write (*,*) 'D',rzsm1_f,tzsm1_f
!         write (*,*) 'E',bulk_dens,a_ice,b_ice,row,roi
!         write (*,*) 'F',rzsm1_u,rzsm1,tzsm1_u,tzsm1
!         write (*,*) 'G',thetas,thetar
!         write (*,*) 'H',rzdthetaudtemp,dt
!         write (*,*) 'I',rzsm_f,tzsm_f
!         write (*,*) 'J',tsoilold,ttt,zsm,tsm,rtdif
         wrong=1

      endif

      return

      end subroutine check_ice
END MODULE MODULE_FROZEN_SOIL
