      subroutine stodeco(ithistraj,istep,cre,cim,pem,nsurf,time,
     &    stodecotime,stodecotau,step)

c This subroutine checks for a stochastic decoherence (SD or STODECO) event, 
c and reinitializes the electronic wave function if such an event is called for.

c A random number X is compared with the probability for decoherence P.
c T is the time since the last hop or frustrated hop; the time at which the last hop or
c frustrated hop occured is stored as STODECOTIME.
c The formula for P is from Jasper and Truhlar, "Non-Born-Oppenheimer
c molecular dynamics of Na...FH photodissociation," J. Chem. Phys, in press (2007).
c TAU is the characteristic decoherence time and may be assigned to a constant
c or to some function. In the present version, TAU is set equal to STODECOTAU, which
c is the formula given in J. Chem. Phys. 123, 064103 (2005). This formula is evaluated
c at every step (for diagnostic purposes) via a call from DRIVER to ELECDECO.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nsurf,i,ithistraj,istep
      double precision cre(mnsurf),cim(mnsurf),pem(mnsurf)
      double precision stodecotime,stodecotau,time,step

      double precision t,tau,p,x

c      t = time - stodecotime
c     the above formula was in the original Na..FH paper, but it is wrong
c     the below formula is correct
      t = step
      tau = stodecotau
      p = dexp(-t/tau)
      x = sprng(rng_stream)

      if (x.gt.p.and.stodecotime.gt.0.d0) then  ! if stodecotime < 0, we have 
                                                ! already decohered once since the last hop and 
                                                ! we don't want to decohere again
!      reinitialize
          write(6,*)"REINITIALIZE to state ",nsurf
          write(33,1033)ithistraj,istep,time*autofs,"D",nsurf,
     &      nsurf,pem(nsurf)*autoev,pem(nsurf)*autoev,t*autofs
 1033     format(i5,i10,f12.3,1x,a,2i5,5f12.3)
       do i=1,nsurft
       cre(i) = 0.d0
       cim(i) = 0.d0
       enddo
       cre(nsurf) = 1.d0
!       stodecotime=-time   ! set to a negative value, so that we don't 
                           ! decohere more than once per hop
      endif

      return

      end
