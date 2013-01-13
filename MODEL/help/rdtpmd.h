      integer icount(MAX_CAT),iopbf,iopwt0,ncatch,nrow,ncol,npix
      integer ipixnum(MAX_ROW,MAX_COL),ixpix(MAX_PIX),iypix(MAX_PIX)
      integer icatch(MAX_PIX),iprn(MAX_FIL)
      integer ioptlr,iopslp

      real*8 pixsiz,q0(MAX_CAT),ff(MAX_CAT),dd(MAX_CAT)
      real*8 area(MAX_CAT),dtil(MAX_CAT),atanb(MAX_PIX)
      real*8 xlength(MAX_CAT),basink(MAX_CAT)
      real*8 zbar1(MAX_CAT),xlamda(MAX_CAT)

      integer kk,jj
      real*8 lat_deg,lat_min,lng_deg,lng_min,lng_mer
      integer iwel(MAX_PIX)

      real*8 atb(MAX_PIX),ti(MAX_PIX),hbar0
      real*8 zbar0(MAX_CAT),sumatb(MAX_CAT),sumlti(MAX_CAT)
      real*8 qb0(MAX_CAT),lte(MAX_CAT)
      real*8 rlatitude,rlongitude,rlng_merid
      real*8 wslp(MAX_PIX,2)
