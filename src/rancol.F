      subroutine rancol(temp0,ecol)

c Generate collision energies from a Maxwell distribution.
c The probability distribution is
c    P(E)dE = E/kT exp(-E/kT) * dE/kT

c This is taken from gamdev from Numerical Recipes
c We select from 
c    p_a(x)dx = x exp(-x) dx/gamma(2)  (i.e., ia = 2)
c and set x = E/kT to obtain E

      implicit none
      include 'param.f'
c      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c in/out
      double precision temp0,ecol

c local
      double precision x,r
      integer ia,i

      ia = 2
      x=1.d0
      do i=1,ia
         r=sprng(rng_stream)
         x=x*r
      enddo
      x=-dlog(x)

      ecol=x*kb*temp0   ! in au

      return

      end
