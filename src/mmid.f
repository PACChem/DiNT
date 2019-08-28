      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,nsurf,nsurft)
c Modified midpoint method from Numerical Recipes
c Used by the Bulirsch-Stoer integrator.
      implicit none
      include 'param.f'
      INTEGER nstep,nvar,nsurf,nsurft
      DOUBLE PRECISION htot,xs,dydx(mnyarray),y(mnyarray),yout(mnyarray)
      INTEGER i,n
      DOUBLE PRECISION h,h2,swap,x,ym(mnyarray),yn(mnyarray)
      h=htot/dble(nstep)
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout,nvar,nsurf)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout,nvar,nsurf)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
