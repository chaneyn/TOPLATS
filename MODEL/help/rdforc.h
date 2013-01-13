      integer iweight(MAX_PIX,MAX_SPP),istanum(MAX_PIX,MAX_SPP)
      integer ipixnum(MAX_ROW,MAX_COL),istaunit,imgunit,npix
      integer nsta,nsp,nrow,ncol,i,img_opt,iyear,iday,ihour
      integer iihour,iiday,iimonth,iiyear,jday
      integer iwel(MAX_PIX)

      real*8 pixdat(MAX_PIX),tsval,rminval,rmaxval,dt
      real*8 rlatitude,rlongitude,rlng_merid
      real*8 wslp(MAX_PIX,2)

      character*100 fnimg(MAX_FIL)

      integer kk
