      integer i,ic,iopbf,ipix,npix
      integer ilandc(MAX_PIX),ivgtyp(MAX_VEG)
      integer icatch(MAX_PIX),isoil(MAX_PIX)
      integer iprn(MAX_FIL),nlcs,natb,int,landc
      integer veg(MAX_PP1,MAX_VST)

      real*8 q0,ff,zbar,dtil
      real*8 basink,dd,xlength
      real*8 gwtsum,capsum,area
      real*8 r_lakearea,dt,etwtsum
      real*8 rzpsum,tzpsum,psicav
      real*8 zrzmax,zbar1,qbreg,zbar1rg
      real*8 psic(MAX_SOI)
      real*8 atb_pdf(MAX_PP1,MAX_ATB),veg_pdf(MAX_PP1,MAX_VST)
      real*8 s_tzsm1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_zw(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 thetas(MAX_SOI),FRT
      real*8 zero,one,two,three,four,five,six

      integer mm

      real*8 qb,hbar,zbrflx,zbrpor,qzbar,dzbar
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
