      integer ndata,npix,nrow,ncol,ipixnum(MAX_ROW,MAX_COL)
      integer nsta_ppt,nsp_ppt,NEWSTORM(MAX_PP2,MAX_TST)
      integer iwppt(MAX_PP2,MAX_SPP),nppt(MAX_PP2,MAX_SPP),img_opt

      real*8 ppt_max,ppt_min,FRC(MAX_PP2,MAX_TST)
      real*8 pptms(MAX_PP2),rainfall(MAX_PP2,MAX_TST),frcbeta

      character*100 fnimg(MAX_FIL)

      integer i,ii,jj,kk,istorm(MAX_TST),j,iyear,iday,ihour
      integer iwel(MAX_PIX),iihour,iiday,iimonth,iiyear,jday

      real*8 tmpval,storm(MAX_TST),max,dt
      real*8 wslp(MAX_PIX,2), rlatitude,rlongitude,rlng_merid
