      integer s_istmst(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_intstm(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_istmst_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_intstm_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_intstp(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_istorm(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_intstp_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_istorm_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer v_m,u_m,v_u,u_u

      real*8 s_xintst_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_xintst(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_wcip1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_wsc(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_wcip1_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_wsc_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))

      integer istmst(MAX_PIX)
      integer intstm(MAX_PIX)
      integer istmst_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstm_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstp(MAX_PIX)
      integer istorm(MAX_PIX)
      integer intstp_moss(1+MOS_FLG*(MAX_PIX-1))
      integer istorm_moss(1+MOS_FLG*(MAX_PIX-1))

      real*8 xintst_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 xintst(MAX_PIX)
      real*8 wcip1(MAX_PIX)
      real*8 wsc(MAX_VEG)
      real*8 wcip1_us(1+UST_FLG*(MAX_PIX-1))
      real*8 wsc_us(1+UST_FLG*(MAX_VEG-1))

      integer nlcs,npix,natb,t,u,v
