      integer npix,nrow,ncol,i,iopwv,ipixnum(MAX_ROW,MAX_COL)
      integer ioppet,img_opt
      integer nsta_ppt,nsta_pet,nsta_ta,nsta_hu,nsta_pa
      integer nsta_ws,nsta_sw,nsta_lw,nsta_rn,nsta_gb
      integer nsp_ppt,nsp_pet,nsp_ta,nsp_hu,nsp_pa
      integer nsp_ws,nsp_sw,nsp_lw,nsp_rn,nsp_gb
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

      real*8 dt
      real*8 ppt_max,pet_max,ta_max,hu_max,pa_max
      real*8 ws_max,sw_max,rlw_max,rn_max,gb_max
      real*8 ppt_min,pet_min,ta_min,hu_min,pa_min
      real*8 ws_min,sw_min,rlw_min,rn_min,gb_min
      real*8 tdry(MAX_PIX),twet(MAX_PIX),rh(MAX_PIX)
      real*8 rnetpn(MAX_PIX),rsd(MAX_PIX),rld(MAX_PIX)
      real*8 press(MAX_PIX),uzw(MAX_PIX)
      real*8 gbspen(MAX_PIX),ebspot(MAX_PIX),pptms(MAX_PIX)

      character*100 fnimg(MAX_FIL)

      integer ii,jj,kk,day,iyear,iday,ihour
      integer iwel(MAX_PIX),iihour,iiday,iimonth,iiyear,jday

      real*8 tmpval,tmprl,tmprs,tmprn,tmpg
      real*8 tmppet,tmpt,tmph,tmpp,tmpw,hour,r_minstep
      real*8 wslp(MAX_PIX,2),rlatitude,rlongitude,rlng_merid

