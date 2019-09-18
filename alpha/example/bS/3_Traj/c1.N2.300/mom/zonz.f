      subroutine zonz(ttt,bmax,m1,m2,zz,zhs,zlj)
      implicit double precision (a-h,o-z)
      double precision mu1,mu2,mu
      character*5 m1,m2

      avagadro=6.0221415d23
      pi=dacos(-1.d0)

      include 'param.inc'

      mu = (mu1*mu2/(mu1+mu2))/(avagadro*1000.d0)       ! red mass is kg
      pref = pi*dsqrt(8.d0*1.380603d-23*ttt/(pi*mu))*1.d6*1.d-20
      zhs = bmax**2*pref     ! hard sphere collision rate for whatever bmax was used in the trajectories
      zlj = si**2*pref/(0.7d0+0.52d0*dlog(0.69502d0*ttt/ep)/dlog(10.d0))       !  LJ collision rate
      zz = zhs/zlj

      write(7,*)"!  ",m1," + ",m2
      write(7,*)"!  T       = ",ttt," K"
      write(7,*)"!  sigma   = ",si," A"
      write(7,*)"!  epsilon = ",ep," cm-1"
      write(7,*)"!  ZLJ     = ",zlj," cm3/s"

      return
      end
