      subroutine ranno(x,a,b)

c picks a random number X evenly from A to B

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

      double precision x,a,b,r

      r = sprng(rng_stream)
      x=a+r*(b-a)

      return
      end

