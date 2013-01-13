      integer iwppt(MAX_PIX,MAX_SPP),iwpet(MAX_PIX,MAX_SPP)
      integer iwta(MAX_PIX,MAX_SPP),iwhu(MAX_PIX,MAX_SPP)
      integer iwpa(MAX_PIX,MAX_SPP),iwws(MAX_PIX,MAX_SPP)
      integer iwsw(MAX_PIX,MAX_SPP),iwlw(MAX_PIX,MAX_SPP)
      integer iwrn(MAX_PIX,MAX_SPP),iwgb(MAX_PIX,MAX_SPP)
      integer nppt(MAX_PIX,MAX_SPP),npet(MAX_PIX,MAX_SPP)
      integer nta(MAX_PIX,MAX_SPP),nhu(MAX_PIX,MAX_SPP)
      integer npa(MAX_PIX,MAX_SPP),nws(MAX_PIX,MAX_SPP)
      integer nsw(MAX_PIX,MAX_SPP),nlw(MAX_PIX,MAX_SPP)
      integer nrn(MAX_PIX,MAX_SPP),ngb(MAX_PIX,MAX_SPP)
      integer nsta_ppt,nsta_pet,nsta_ta,nsta_hu,nsta_pa,nsta_ws
      integer nsta_sw,nsta_lw,nsta_rn,nsta_gb
      integer nsp_ppt,nsp_pet,nsp_ta,nsp_hu,nsp_pa,nsp_ws
      integer nsp_sw,nsp_lw,nsp_rn,nsp_gb
      integer npix,iprn(MAX_FIL)

      real*8 ppt_min,pet_min,ta_min,hu_min,pa_min,ws_min
      real*8 sw_min,rlw_min,rn_min,gb_min
      real*8 ppt_max,pet_max,ta_max,hu_max,pa_max,ws_max
      real*8 sw_max,rlw_max,rn_max,gb_max 

      integer kk,jj
