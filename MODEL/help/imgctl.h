      real*8 tkact(MAX_PIX),gact(MAX_PIX),hact(MAX_PIX)
      real*8 xleact(MAX_PIX),rnact(MAX_PIX),etpix(MAX_PIX)
      real*8 runtot(MAX_PIX),xinact(MAX_PIX)
      real*8 zw(MAX_PIX),tzsm1(MAX_PIX),rzsm1(MAX_PIX)
      real*8 tkpet(MAX_PIX),wcip1(MAX_PIX),gpet(MAX_PIX)
      real*8 hpet(MAX_PIX),xlepet(MAX_PIX),rnpet(MAX_PIX)
      real*8 r_mossm(1+MOS_FLG*(MAX_PIX-1)),pptms(MAX_PIX),pptsumrg
      real*8 Swq(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Swq_us(1+SNW_FLG*(MAX_PIX-1))

      integer ievcon(MAX_PIX),ipixnum(MAX_ROW,MAX_COL),ncol,nrow
      integer icurser(MAX_FIL),nseries(MAX_FIL),iouten(MAX_FIL,MAX_SER)
      integer ioutst(MAX_FIL,MAX_SER),iprn(MAX_FIL),i
      integer ioutsp(MAX_FIL,MAX_SER),irntyp(MAX_PIX)
      integer img_opt,iyear,iday,ihour

      character*100 fnimg(MAX_FIL)

      integer imgprn(MAX_FIL),iu
