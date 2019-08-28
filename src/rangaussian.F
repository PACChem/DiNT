      subroutine rangaussian(xx,xe,alpha)

c Generate momenta randomly according to a Gaussian distribution.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c input/output
      double precision xx,alpha,xe

c local
      double precision r1,r2,x1,x2,w,y1,y2

c convert uniformly distributed numbers to a gaussian distributed set
c polar form of the Box-Muller transformation
c transformed random numbers are from a distribution that has zero mean
c and unit standard deviation
 10      continue
         r1 = sprng(rng_stream)
         r2 = sprng(rng_stream)
         x1=2.d0*r1-1.d0
         x2=2.d0*r2-1.d0
         w = x1*x1+x2*x2
         if (w.ge.1.d0.or.w.eq.0.d0) go to 10
         w = dsqrt(-2.d0*dlog(w)/w)
         y1 = x1*w
         y2 = x2*w
c you automatically get two random numbers (y1 and y2)
c i'm just going to use one for now, inefficient but simpler
       xx=y1*alpha+xe

      if (xx.gt.1.d0) go to 10
      if (xx.lt.0.d0) go to 10

      return

      end
