      integer iwater,lakpix,numnod

      real*8 qbot,qw,sw,lnet,Qe,Qh,eta,dt,surface(MAX_PIX,MAX_NOD)
      real*8 T(MAX_NOD,2), de(MAX_NOD), dnsty(MAX_NOD)
      real*8 cpz(MAX_NOD), z(MAX_NOD), zhalf(MAX_NOD)
      real*8 a(MAX_NOD), b(MAX_NOD), c(MAX_NOD), d(MAX_NOD)
      real*8 told(MAX_NOD), tnew(MAX_NOD)
      real*8 work1(MAX_NOD), work2(MAX_NOD)

      integer k

      real*8 dist12,surface_1,surface_2,T1,cnextra,swtop,bot,top
