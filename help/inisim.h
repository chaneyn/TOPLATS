      integer iopsmini,nrow,ncol,ipixnum(MAX_ROW,MAX_COL)
      integer ilandc(MAX_PIX),npix,inc_frozen
      integer istorm(MAX_PIX),intstm(MAX_PIX),istmst(MAX_PIX)
      integer intstp(MAX_PIX),istorm_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstm_moss(1+MOS_FLG*(MAX_PIX-1))
      integer istmst_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstp_moss(1+MOS_FLG*(MAX_PIX-1)),isoil(MAX_PIX)
      integer idifind(MAX_SOI),s_kk,sw_kk

      real*8 smpet0,r_mossmpet0(1+MOS_FLG*(MAX_VEG-1)),endstm
      real*8 rzsm1(MAX_PIX),tzsm1(MAX_PIX)
      real*8 r_mossm1(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm(1+MOS_FLG*(MAX_PIX-1))
      real*8 rzsm1_u(MAX_PIX),tzsm1_u(MAX_PIX)
      real*8 rzsm1_f(MAX_PIX),tzsm1_f(MAX_PIX)
      real*8 r_mossm1_u(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm_u(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm1_f(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm_f(1+MOS_FLG*(MAX_PIX-1))
      real*8 rzdthetaidt(MAX_PIX),tzdthetaidt(MAX_PIX)
      real*8 zmoss(1+MOS_FLG*(MAX_VST-1)),r_moss_depth(1+MOS_FLG*(MAX_PIX-1))
      real*8 thetas_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 xintst(MAX_PIX),xintst_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 cuminf(MAX_PIX),xk0(MAX_SOI),psic(MAX_SOI)
      real*8 thetas(MAX_SOI),thetar(MAX_SOI),bcgamm(MAX_SOI)
      real*8 bcbeta(MAX_SOI),sorp(MAX_PIX),cc(MAX_PIX),dt
      real*8 sesq(MAX_PIX),corr(MAX_SOI),par(MAX_SOI)
      real*8 PackWater_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 SurfWater_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Swq_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 VaporMassFlux_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 r_MeltEnergy_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Outflow_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 PackWater(1+SNOW_RUN*(MAX_PIX-1))
      real*8 SurfWater(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Swq(1+SNOW_RUN*(MAX_PIX-1))
      real*8 VaporMassFlux(1+SNOW_RUN*(MAX_PIX-1))
      real*8 r_MeltEnergy(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Outflow(1+SNOW_RUN*(MAX_PIX-1))
      real*8 zero,one,two,three,four,five,six
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/


      integer istep(MAX_PIX),kk,iopflg,istflg,m_kk,v_kk

      real*8 cumdep(MAX_PIX),smbeg(MAX_PIX),deltrz,relrze
