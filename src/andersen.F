      subroutine andersen(pp,mmm,nclu,step,lhit)

c Andersen thermostat, H. C. Andersen, JCP 72, 2348, 1980.
c Collisions occur with a probability of probhit*timestep.
c In general, you can choose a single atom as a "particle" or the whole system or
c anything in between.  We choose the entire system.
c If a random number is < the collision prob, then the momenta for each component
c of every atom is replaced with momenta selected from a Boltzmann distribution.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu
      double precision pp(3,mnat),mmm(mnat),temp,probhit,step
      logical lhit

c local
      integer i,ix,iy,iz,j
      real ran(2),rtmp
      double precision r1,r2,x1,x2,w,y1,y2,tempnow,kenow,pei,pe

c frequency of collisions
      probhit = andersen_freq
c temperature
      temp = andersen_temp

c get random number
      rtmp = sprng(rng_stream)
c check for collision
      IF (dble(rtmp).le.probhit*step) THEN
      write(6,*)"hit!"
      lhit = .true.

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
         pp(j,i)=y1*dsqrt( kb * temp * mmm(i) )
         enddo
      enddo

      ENDIF


      return

      end
