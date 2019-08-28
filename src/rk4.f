      subroutine rk4(y,dydx,nv,hstep,nsurf)
c     Fourth order Runge-Kutta integrator. From numerical recipes.

      implicit none
      include 'param.f'

      integer nv
      double precision y(mnyarray),yt(mnyarray),dydx(mnyarray)

      integer i,j,nsurf
      double precision hh,h6,dum,hstep
      double precision dym(mnyarray),dyt(mnyarray)

      hh=hstep/2.d0
      h6=hstep/6.d0
      dum=0.d0

c     first step
      do 11 i=1,nv
      yt(i)=y(i)+hh*dydx(i)
 11   continue

c     second step
      call derivs(dum,yt,dyt,nv,nsurf)
      do 12 i=1,nv
      yt(i)=y(i)+hh*dyt(i)
 12   continue

c     third step
      call derivs(dum,yt,dym,nv,nsurf)
      do 13 i=1,nv
      yt(i)=y(i)+hstep*dym(i)
      dym(i)=dyt(i)+dym(i)
 13   continue

c     fourth step
      call derivs(dum,yt,dyt,nv,nsurf)
c     accumulate increments with proper weights
      do 14 i=1,nv
      y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
 14   continue

      return
      end
