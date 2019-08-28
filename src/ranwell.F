      subroutine ranwell(whwell,rel1)

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer whwell
      double precision rel1,r

      r=sprng(rng_stream)
      if (r.lt.rel1) then
        whwell=1
      else if (r.gt.rel1) then
        whwell=2
      else
        write (6,*) "ranwell: equal prob 1 or 2"
        whwell=2
      endif

      return
      end
  
