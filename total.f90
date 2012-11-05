! ====================================================================
!
!                              subroutine total
!
! ====================================================================
!
! Subroutine to calculate the average flux or state variable over a
! pixel given the distribution of values in it.
!
! ====================================================================

      subroutine total(tkact,gact,hact,xleact,rnact,s_tkact,s_gact,s_hact,&
       s_xleact,s_rnact,s_nr_tkact,s_nr_gact,s_nr_hact,s_nr_xleact,s_nr_rnact,&
       runtot,xinact,irntyp,zw,tzsm1,rzsm1,r_mossm,s_r_mossm,s_nr_r_mossm,&
       s_runtot,s_xinact,s_irntyp,s_zw,s_tzsm1,s_rzsm1,s_nr_runtot,&
       s_nr_xinact,s_nr_irntyp,s_nr_zw,s_nr_tzsm1,s_nr_rzsm1,tkpet,wcip1,&
       gpet,hpet,xlepet,rnpet,s_tkpet,s_wcip1,s_gpet,s_hpet,s_xlepet,s_rnpet,&
       s_nr_tkpet,s_nr_wcip1,s_nr_gpet,s_nr_hpet,s_nr_xlepet,s_nr_rnpet,&
       s_nr_rzsm,s_nr_ievcon,rzsm,s_rzsm,ievcon,s_ievcon,Swq,s_Swq,s_nr_Swq,&
       Swq_us,s_Swq_us,s_nr_Swq_us,atb_pdf,veg_pdf,natb,npix,nlcs,&
       i,FRC,ivgtyp,FRCOV,veg)

      implicit none
      include "SNOW.h"
      include "wgtpar.h"
      include "help/total.h"

      do ipix=1,npix

         if (MOS_FLG.eq.1) m_px=ipix
         if (MOS_FLG.eq.0) m_px=1
         s_px=1+SNOW_RUN*(ipix-1)
         sw_px=1+SNW_FLG*(ipix-1)

         rnact(ipix)=0.d0
         xleact(ipix)=0.d0
         hact(ipix)=0.d0
         gact(ipix)=0.d0
         tkact(ipix)=0.d0
         xinact(ipix)=0.d0
         runtot(ipix)=0.d0
         zw(ipix)=0.d0
         rzsm(ipix)=0.d0
         rzsm1(ipix)=0.d0
         tzsm1(ipix)=0.d0
         rnpet(ipix)=0.d0
         xlepet(ipix)=0.d0
         hpet(ipix)=0.d0
         gpet(ipix)=0.d0
         tkpet(ipix)=0.d0
         wcip1(ipix)=0.d0
         t1(ipix)=0.d0
         t2(ipix)=0.d0
         r_mossm(m_px)=0.d0
         Swq(s_px)=0.d0
         Swq_us(sw_px)=0.d0

      enddo

      do ipix=1,npix

         sum(ipix)=0.d0

      enddo

      do ipix=1,npix

         do landc=1,nlcs

            if (ivgtyp(veg(ipix,landc)).ge.0) then

               sum(ipix)=sum(ipix)+veg_pdf(ipix,landc)

            endif

         enddo

      enddo

      do int=1,natb

         do landc=1,nlcs

            do ipix=1,npix

               if (MOS_FLG.eq.1) m_lc=landc
               if (MOS_FLG.eq.0) m_lc=1
               if (MOS_FLG.eq.1) m_px=ipix
               if (MOS_FLG.eq.0) m_px=1
               s_px=1+SNOW_RUN*(ipix-1)
               sw_px=1+SNW_FLG*(ipix-1)
               s_lc=1+SNOW_RUN*(landc-1)
               sw_lc=1+SNW_FLG*(landc-1)

               if (FRCOV.eq.0) FRC_VALUE=1.d0
               if (FRCOV.eq.1) FRC_VALUE=FRC(ipix,i)

               if ( (atb_pdf(ipix,int)*veg_pdf(ipix,landc).gt.(0.d0)).and.&
                    (FRCOV.eq.1) ) then

                  Swq(s_px)=Swq(s_px)+s_Swq(int,s_lc,s_px)*&
                            atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                            s_nr_Swq(int,s_lc,s_px)*&
                            atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                            (1.d0-FRC_VALUE)
                  Swq_us(s_px)=Swq_us(s_px)+s_Swq_us(int,s_lc,s_px)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                               s_nr_Swq_us(int,s_lc,s_px)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                               (1.d0-FRC_VALUE)
                  rnact(ipix)=rnact(ipix)+s_rnact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_rnact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)
                  xleact(ipix)=xleact(ipix)+s_xleact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                               s_nr_xleact(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                               (1.d0-FRC_VALUE)
                  hact(ipix)=hact(ipix)+s_hact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                             s_nr_hact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                             (1.d0-FRC_VALUE)
                  gact(ipix)=gact(ipix)+s_gact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                             s_nr_gact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                             (1.d0-FRC_VALUE)
                  tkact(ipix)=tkact(ipix)+s_tkact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_tkact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)
                  xinact(ipix)=xinact(ipix)+s_xinact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                               s_nr_xinact(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                               (1.d0-FRC_VALUE)
                  runtot(ipix)=runtot(ipix)+s_runtot(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                               s_nr_runtot(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                               (1.d0-FRC_VALUE)

                  if (ivgtyp(veg(ipix,landc)).ge.0) then

                     zw(ipix)=zw(ipix)+s_zw(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_zw(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)

                     rzsm(ipix)=rzsm(ipix)+s_rzsm(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                                s_nr_rzsm(int,landc,ipix)*&
                                atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                                (1.d0-FRC_VALUE)
                     rzsm1(ipix)=rzsm1(ipix)+s_rzsm1(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                                 s_nr_rzsm1(int,landc,ipix)*&
                                 atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                                 (1.d0-FRC_VALUE)
                     tzsm1(ipix)=tzsm1(ipix)+s_tzsm1(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                                 s_nr_tzsm1(int,landc,ipix)*&
                                 atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                                 (1.d0-FRC_VALUE)

                     r_mossm(m_px)=r_mossm(m_px)+s_r_mossm(int,m_lc,m_px)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                                   s_nr_r_mossm(int,m_lc,m_px)*&
                                   atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                                   (1.d0-FRC_VALUE)

                  endif

                  rnpet(ipix)=rnpet(ipix)+s_rnpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_rnpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)
                  xlepet(ipix)=xlepet(ipix)+s_xlepet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                               s_nr_xlepet(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                               (1.d0-FRC_VALUE)
                  hpet(ipix)=hpet(ipix)+s_hpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                             s_nr_hpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                             (1.d0-FRC_VALUE)
                  gpet(ipix)=gpet(ipix)+s_gpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                             s_nr_gpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                             (1.d0-FRC_VALUE)
                  tkpet(ipix)=tkpet(ipix)+s_tkpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_tkpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)
                  wcip1(ipix)=wcip1(ipix)+s_wcip1(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                              s_nr_wcip1(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                              (1.d0-FRC_VALUE)

                  t1(ipix)=t1(ipix)+real(s_irntyp(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                           real(s_nr_irntyp(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                           (1.d0-FRC_VALUE)
                  t2(ipix)=t2(ipix)+real(s_ievcon(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)*FRC_VALUE+&
                           real(s_nr_ievcon(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)*&
                           (1.d0-FRC_VALUE)

               endif

               if ( (atb_pdf(ipix,int)*veg_pdf(ipix,landc).gt.(0.d0)).and.&
                    (FRCOV.eq.0) ) then

                  Swq(s_px)=Swq(s_px)+s_Swq(int,s_lc,s_px)*&
                            atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  Swq_us(s_px)=Swq_us(s_px)+s_Swq_us(int,s_lc,s_px)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  rnact(ipix)=rnact(ipix)+s_rnact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  xleact(ipix)=xleact(ipix)+s_xleact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  hact(ipix)=hact(ipix)+s_hact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  gact(ipix)=gact(ipix)+s_gact(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  tkact(ipix)=tkact(ipix)+s_tkact(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  xinact(ipix)=xinact(ipix)+s_xinact(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  runtot(ipix)=runtot(ipix)+s_runtot(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)

                  if (ivgtyp(veg(ipix,landc)).ge.0) then

                     zw(ipix)=zw(ipix)+s_zw(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)

                     rzsm(ipix)=rzsm(ipix)+s_rzsm(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                     rzsm1(ipix)=rzsm1(ipix)+s_rzsm1(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                     tzsm1(ipix)=tzsm1(ipix)+s_tzsm1(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)

                     r_mossm(m_px)=r_mossm(m_px)+s_r_mossm(int,m_lc,m_px)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)

                  endif

                  rnpet(ipix)=rnpet(ipix)+s_rnpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  xlepet(ipix)=xlepet(ipix)+s_xlepet(int,landc,ipix)*&
                               atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  hpet(ipix)=hpet(ipix)+s_hpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  gpet(ipix)=gpet(ipix)+s_gpet(int,landc,ipix)*&
                             atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  tkpet(ipix)=tkpet(ipix)+s_tkpet(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  wcip1(ipix)=wcip1(ipix)+s_wcip1(int,landc,ipix)*&
                              atb_pdf(ipix,int)*veg_pdf(ipix,landc)

                  t1(ipix)=t1(ipix)+real(s_irntyp(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)
                  t2(ipix)=t2(ipix)+real(s_ievcon(int,landc,ipix))*&
                           atb_pdf(ipix,int)*veg_pdf(ipix,landc)

               endif

            enddo

         enddo

      enddo

      do ipix=1,npix

         irntyp(ipix)=t1(ipix)
         ievcon(ipix)=t2(ipix)

         if (sum(ipix).gt.0.d0) then

            zw(ipix)=zw(ipix)/sum(ipix)

         endif

      enddo

      return

      end
