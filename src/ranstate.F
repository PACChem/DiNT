      subroutine ranstate(freq,nmqn,temp)

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

c local
      integer nmax,i
      double precision temp,freq,ptot,ei,pp,p
      double precision r1,nmqn

      nmax=100
      ptot=0.d0
      do i=0,nmax
        ei=freq*(dble(i)+0.5d0)
        ptot=ptot+dexp(-ei/(temp*kb))
      enddo

      r1 = sprng(rng_stream)

      pp=0.d0
      nmqn=-1.d0
      do i=0,nmax
        ei=freq*(dble(i)+0.5d0)
        pp=pp+dexp(-ei/(temp*kb))
        p = pp/ptot
        if (r1.lt.p.and.nmqn.lt.0.d0) nmqn = dble(i)
      enddo
      if (nmqn.lt.0.d0) then
c     something went wrong. this can happen when the temperature is small and all the boltzmann weights are ~zero.
      nmqn=0.d0
      endif

      return

      end
