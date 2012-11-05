      integer MAX_NOD
      real*8 dz_lake,delta,rhowat,rhosnow,rhoice,rhosurf,Le,Lei
      real*8 fusion,surf,fracmin,fraclim,qwtau,cpw_ice,beta,dm,pi
      real*8 snocrit
      parameter (dz_lake=.1)
      parameter ( MAX_NOD = 5 )

      logical iceflag
      parameter (iceflag = .true.)

      logical ar_flag
      parameter (ar_flag = .true.)
 
      parameter (delta=5.67e-8)   ! s-b constant
      parameter (rhowat = 1000.)
      parameter (rhosnow = 330.)  ! densities of water and snow 
      parameter (rhoice = 917.)   ! ice density
      parameter (rhosurf=1.275)   ! surface air density
      parameter (Le = 2.25e6)
      parameter (Lei = 2.5e6)     ! latent heats 
      parameter (fusion=3.34e5)   ! heat of fusion
      parameter (surf = 0.1)     ! surface thickness for E-B computation
      parameter (fracmin= 0.01)   ! min ice thickness in meters
      parameter (fraclim = 0.02)  ! lower limit on fractional ice cover
      parameter (qwtau=86400./2.) ! D. Pollard sub-ice time constant
      parameter (cpw_ice = 4200.) ! specific heat of ice
      parameter (beta = 0.4)      ! portion of s-w absorbed at the surface 
      real*8 kv 
      parameter (kv= 0.4)         ! vonkarman constant
      parameter (dm=1.38889E-07)  ! molecular diffusivity of water
      parameter (pi= 3.141592654)
      parameter (snocrit = 0.05)  ! for albedo, in m
