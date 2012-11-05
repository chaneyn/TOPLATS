      integer i,ic,iopbf,npix
      integer ilandc(MAX_PIX),ivgtyp(MAX_VEG)
      integer icatch(MAX_PIX),isoil(MAX_PIX)
      integer iprn(MAX_FIL)

      real*8 q0,ff,zbar,dtil
      real*8 basink,dd,xlength
      real*8 gwtsum,capsum,area
      real*8 r_lakearea,dt,etwtsum
      real*8 rzpsum,tzpsum,psicav
      real*8 zrzmax,zbar1,qbreg,zbar1rg
      real*8 psic(MAX_SOI),tzsm1(MAX_PIX),rzsm1(MAX_PIX)
      real*8 zw(MAX_PIX),thetas(MAX_SOI)
      real*8 zero,one,two,three,four,five,six
      real*8 pixsiz

      integer mm

      real*8 qb,hbar,zbrflx,zbrpor,qzbar,dzbar
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
