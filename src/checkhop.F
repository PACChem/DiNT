      subroutine checkhop(phop,nsurf,lhop,newsurf)

c Used for surface hopping methods.
c Draw a random number and compare it to PHOP to check for a hop.
c If a hop occurs, LHOP = .TRUE. and NEWSURF = whichever surface
c the system is trying to hop to.  The hop does not actually 
c occur in this subroutine.  Hops (frustrated and successful) 
c occur in HOP.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c input
      double precision phop(mnsurf)
      integer nsurf

c output
      logical lhop
      integer newsurf

c local
      integer i
      double precision dxi,tmp

c get random number
      dxi = sprng(rng_stream)

c check for hop
      newsurf = 0
      tmp = 0.d0
      do i=1,nsurft
      tmp = tmp + phop(i)
      if (dxi.le.tmp.and.newsurf.eq.0) newsurf = i
      enddo
      if (newsurf.eq.0) then
        write(6,*)"NEWSURF = 0 in HOPCHECK!"
        stop
      endif

      if (newsurf.ne.nsurf) then
         lhop = .true.
      endif

      return

      end
