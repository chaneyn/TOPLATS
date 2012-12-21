      integer nsoil,irestype,ikopt,iopsmini
      integer isoil(MAX_PIX),nrow,ncol,ifcoarse(MAX_SOI),inc_frozen
      integer idifind(MAX_SOI),ipixnum(MAX_ROW,MAX_COL),icatch(MAX_PIX)
      integer ncatch,npix,iprn(MAX_FIL)

      real*8 zrzmax,smpet0,bcbeta(MAX_SOI),psic(MAX_SOI)
      real*8 thetas(MAX_SOI),thetar(MAX_SOI),xk0(MAX_SOI),zdeep(MAX_SOI)
      real*8 tdeep(MAX_SOI),zmid(MAX_SOI),tmid0(MAX_SOI),rocpsoil(MAX_SOI)
      real*8 quartz(MAX_SOI),srespar1(MAX_SOI),srespar2(MAX_SOI)
      real*8 srespar3(MAX_SOI),a_ice(MAX_SOI),b_ice(MAX_SOI)
      real*8 bulk_dens(MAX_SOI),amp(MAX_SOI),phase(MAX_SOI)
      real*8 shift(MAX_SOI),bcgamm(MAX_SOI),par(MAX_SOI),corr(MAX_SOI)
      real*8 pixsiz,area(MAX_CAT),psicav(MAX_CAT)
      real*8 tc(MAX_SOI),tw(MAX_SOI)
      real*8 zero,one,two,three,four,five,six

      integer icount(MAX_SOI,MAX_CAT+1),jj,kk,nn

      real*8 frsoil(MAX_SOI,MAX_CAT+1),tempsum,dtaken
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
