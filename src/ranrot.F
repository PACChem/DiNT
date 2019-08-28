      subroutine ranrot(pp,xx,mmm,nclu,temp0)

c Generate momenta randomly according to a Maxwell distribution.
c Subtract off momentum that is not rotational.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu
      double precision pp(3,mnat),mmm(mnat),temp0,pp2(3,mnat),
     &   pp1(3,mnat),xx(3,mnat)

c local
      integer i,ix,iy,iz,j
      real ran(2)
      double precision r1,r2,x1,x2,w,y1,y2

c convert uniformly distributed numbers to a gaussian distributed set
c polar form of the Box-Muller transformation
c transformed random numbers are from a distribution that has zero mean
c and unit standard deviation
      do i=1,nclu
      do j=1,3
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
c        transform to velocity kb*tempt, mmm, and pp are in a.u.
         pp2(j,i)=y1*dsqrt( kb * temp0 * mmm(i) )
         pp1(j,i)=pp2(j,i)
      enddo
      enddo

      call noang(xx,pp1,mmm,nclu)
      do i=1,nclu
      do j=1,3
        pp2(j,i)=pp2(j,i)-pp1(j,i)
      enddo
      enddo

      do i=1,nclu
      do j=1,3
        pp(j,i)=pp(j,i)+pp2(j,i)
      enddo
      enddo

      return

      end
