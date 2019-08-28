      subroutine decocheck(phop,nsurf,lhop)

c Check for a switch of the decoherent state for the SCDM and CSDM 
c methods.  Draw and random number and compare it with PHOP.  If
c a switch occurs, set NSURF to the new surface and set LHOP = .TRUE.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c input
      integer nclu
      double precision phop(mnsurf)

c input/output
      integer nsurf
      logical lhop

c local
      integer newsurf,i
      double precision dxi,tmp

c get random number
      dxi = sprng(rng_stream)

c check for switch
      newsurf = 0
      tmp = 0.d0
      do i=1,nsurft
      tmp = tmp + phop(i)
      if (dxi.le.tmp.and.newsurf.eq.0) newsurf = i
c      print *,i,dxi,tmp,newsurf
      enddo
      if (newsurf.eq.0) then
        write(6,*)"NEWSURF = 0 in DECOCHECK!"
        stop
      endif

c switch!
      if (newsurf.ne.nsurf) then
        write(6,*)"DECO SWITCH from ",nsurf,"-->",newsurf
        nsurf = newsurf
        lhop = .true.
      endif

      return
      end
