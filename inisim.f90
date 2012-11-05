! ====================================================================
!
!			subroutine inisim
!
! ====================================================================
!
! Subroutine to initialize simulation total water balance variables.
!
! ====================================================================

      subroutine inisim(iopsmini,nrow,ncol,ipixnum,ilandc,npix,inc_frozen,&
       istorm,intstm,istmst,intstp,istorm_moss,intstm_moss,istmst_moss,&
       intstp_moss,isoil,idifind,smpet0,r_mossmpet0,endstm,rzsm1,&
       tzsm1,r_mossm1,r_mossm,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f,r_mossm1_u,&
       r_mossm_u,r_mossm1_f,r_mossm_f,rzdthetaidt,tzdthetaidt,zmoss,&
       r_moss_depth,thetas_moss,xintst,xintst_moss,cuminf,xk0,psic,thetas,&
       thetar,bcgamm,bcbeta,sorp,cc,dt,sesq,corr,par,&
       PackWater_us,SurfWater_us,Swq_us,VaporMassFlux_us,r_MeltEnergy_us,&
       Outflow_us,PackWater,SurfWater,Swq,VaporMassFlux,r_MeltEnergy,Outflow)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/inisim.h"

      if (iopsmini.eq.1) then

! ====================================================================
! If initial root and transmission zone soil moistures are    
! input as images then read the images and assign to the
! initial condition.
! ====================================================================

! --------------------------------------------------------------------
! Alert the user.
! --------------------------------------------------------------------

	 print *, "Reading initial soil moisture images"

         call rdimgr(rzsm1,14,nrow,ncol,ipixnum)
         call rdimgr(tzsm1,15,nrow,ncol,ipixnum)

	 close(14)
	 close(15)

      else

! ====================================================================
! If initial root zone is not entered, then set the root
! zone soil moisture for calculation of thermal conductivity.
! This is changed later (initsm) to a new initial condition      
! based on Brooks-Corey and local water table depth.
! ====================================================================

         do 50 kk=1,npix

            if (MOS_FLG.eq.1) then

               m_kk=kk
               v_kk=ilandc(kk)

            endif

            if (MOS_FLG.eq.0) then

               m_kk=1
               v_kk=1

            endif

            rzsm1(kk) = smpet0
            tzsm1(kk) = smpet0
            r_mossm1(m_kk) = r_mossmpet0(v_kk)
            r_mossm(m_kk) = r_mossmpet0(v_kk)
            rzsm1_u(kk) = smpet0
            tzsm1_u(kk) = smpet0
            rzsm1_f(kk) = 0.d0
            tzsm1_f(kk) = 0.d0
            r_mossm1_u(m_kk) = r_mossmpet0(v_kk)
            r_mossm_u(m_kk) = r_mossmpet0(v_kk)
            r_mossm1_f(m_kk) = 0.d0
            r_mossm_f(m_kk) = 0.d0
            rzdthetaidt(kk)=0.d0
            tzdthetaidt(kk)=0.d0

50       continue

      endif

      do kk=1,npix

         if (MOS_FLG.eq.1) then

            m_kk=kk
            v_kk=ilandc(kk)

         endif

         if (MOS_FLG.eq.0) then

            m_kk=1
            v_kk=1

         endif

         if (inc_frozen.eq.0) then

          if(thetas_moss(v_kk).eq.0.)then
            zmoss(m_kk)=0.
          else
            zmoss(m_kk)=r_moss_depth(v_kk)*r_mossm(m_kk)/&
                        thetas_moss(v_kk)
          endif
         else
          if(thetas_moss(v_kk).eq.0.)then
            zmoss(m_kk)=0.
          else

            zmoss(m_kk)=r_moss_depth(v_kk)*r_mossm_u(m_kk)/&
                        thetas_moss(v_kk)
          endif

         endif

      enddo

! ====================================================================
! Read data to tell how program will initialize the storm
! and interstorm event flags and times.
! ====================================================================

      read(1000,*) iopflg

      if (iopflg.eq.0) then

! --------------------------------------------------------------------
! If one event flag value is used then set all flags and times
! accordingly.
! --------------------------------------------------------------------

         read(1000,*) istflg

         if (istflg.eq.1) then

! ....................................................................
! If the event is a storm event.
! ....................................................................

            do 100 kk=1,npix

               if (MOS_FLG.eq.1) m_kk=kk
               if (MOS_FLG.eq.0) m_kk=1

               istorm(kk) = 1
               intstm(kk) = 0
               istmst(kk) = 0
               intstp(kk) = 0
               xintst(kk) = 0.0
               istorm_moss(m_kk) = 1
               intstm_moss(m_kk) = 0
               istmst_moss(m_kk) = 0
               intstp_moss(m_kk) = 0
               xintst_moss(m_kk) = 0.0

100         continue

         else

! ....................................................................
! If the event is an interstorm event.
! ....................................................................

            do 200 kk=1,npix

               if (MOS_FLG.eq.1) m_kk=kk
               if (MOS_FLG.eq.0) m_kk=1

               istorm(kk) = 0
               intstm(kk) = 1
               istmst(kk) = 0
               intstp(kk) = 0
               xintst(kk) = endstm
               istorm_moss(m_kk) = 0
               intstm_moss(m_kk) = 1
               istmst_moss(m_kk) = 0
               intstp_moss(m_kk) = 0
               xintst_moss(m_kk) = endstm

200         continue

         endif

      else

! --------------------------------------------------------------------
! If spatially varying data is used read image data.
! --------------------------------------------------------------------

         call rdimgi(istorm,16,nrow,ncol,ipixnum)
         call rdimgi(istep,17,nrow,ncol,ipixnum)
         call rdimgr(cumdep,18,nrow,ncol,ipixnum)
         call rdimgr(smbeg,19,nrow,ncol,ipixnum)

! --------------------------------------------------------------------
! Set event flags, time of events, cumulative values.
! --------------------------------------------------------------------

         do 300 kk=1,npix

! ....................................................................
! For pixels under storm event.
! ....................................................................

            if (istorm(kk).eq.1) then

               intstm(kk) = 0
               istmst(kk) = istep(kk)
               intstp(kk) = 0
               xintst(kk) = 0.0
               cuminf(kk) = cumdep(kk)

! ....................................................................
! Find philip's equation parameters.
! ....................................................................

               sorp(kk) = (((two*xk0(isoil(kk))*&
                          ((thetas(isoil(kk))-smbeg(kk))**two)&
                          *psic(isoil(kk)))/&
                          (thetas(isoil(kk))-thetar(isoil(kk))))*&
                          ((one/(bcgamm(isoil(kk))+&
                               0.5d0*bcbeta(isoil(kk))-one)) +&
                          ((thetas(isoil(kk))-thetar(isoil(kk))) /&
                           (thetas(isoil(kk))-smbeg(kk))))) ** 0.5d0
               deltrz = smbeg(kk) - thetar(isoil(kk))
               if (deltrz.le.zero) deltrz=zero
               cc(kk) = 0.5d0 *&
                       (one+((deltrz/&
                              (thetas(isoil(kk))-thetar(isoil(kk)))) **&
                             (bcgamm(isoil(kk))/bcbeta(isoil(kk)))))

! ....................................................................
! For pixels under interstorm event.
! ....................................................................

            else

               intstm(kk) = 1
               istmst(kk) = 0
               intstp(kk) = istep(kk)
               xintst(kk) = intstp(kk)*dt + endstm
               
               relrze = (smbeg(kk) - thetar(isoil(kk))) /&
                        (thetas(isoil(kk)) - thetar(isoil(kk)))

               if (relrze.le.zero) relrze=zero
               if (relrze.ge.one) relrze=one


            endif

300      continue

      endif

! ====================================================================
! Initialize snow pack variables.
! ====================================================================

      do kk=1,npix

         if (SNW_FLG.eq.1) s_kk=kk
         if (SNW_FLG.eq.0) s_kk=1
         if (SNOW_RUN.eq.1) sw_kk=kk
         if (SNOW_RUN.eq.0) sw_kk=1

         PackWater_us(s_kk)=zero
         SurfWater_us(s_kk)=zero
         Swq_us(s_kk)=zero
         VaporMassFlux_us(s_kk)=zero
         r_MeltEnergy_us(s_kk)=zero
         Outflow_us(s_kk)=zero
         PackWater(sw_kk)=zero
         SurfWater(sw_kk)=zero
         Swq(sw_kk)=zero
         VaporMassFlux(sw_kk)=zero
         r_MeltEnergy(sw_kk)=zero
         Outflow(sw_kk)=zero

      enddo

      return

      end
