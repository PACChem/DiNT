      subroutine no2arbeam(zz)

      implicit none

      include 'param.f'
      double precision r1,r2,x1,x2,w,y1,y2,zz

c convert uniformly distributed numbers to a gaussian distributed set
c polar form of the Box-Muller transformation
c transformed random numbers are from a distribution that has zero mean
c and unit standard deviation
 10      continue
         call ranno(r1,0.d0,1.d0)
         call ranno(r2,0.d0,1.d0)
         x1=2.d0*r1-1.d0
         x2=2.d0*r2-1.d0
         w = x1*x1+x2*x2
         if (w.ge.1.d0.or.w.eq.0.d0) go to 10
         w = dsqrt(-2.d0*dlog(w)/w)
         y1 = x1*w
         y2 = x2*w

         zz=560.d0+y1*60.d0
         zz=max(zz,10.d0)
c         print *,'hihi',zz
         zz=zz/autocmi


      return

      end
