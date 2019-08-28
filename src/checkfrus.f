      subroutine checkfrus(nclu,nsurf,pp,pem,dvec,frusflag,tu_maxt)

c Used for the FSTU method.
c At every step, check to see of a hop would be energetically allowed
c to all other states.  If a hop is allowed to state K, set FRUSFLAG(K) = 0.
c Otherwise, FRUSFLAG(K) = 1.

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer nclu,frusflag(mnsurf)
      integer nsurf

      integer is,i,j
      double precision pem(mnsurf),gpem(3,mnat,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),pp(3,mnat),
     & sectrm,a,b,a2,ab,radi,scr,deltae,tu_maxt(mnsurf),
     & talongd

c check for forbidden hops
      do is = 1,nsurft
       if (is.ne.nsurf) then
        deltae = pem(is) - pem(nsurf)
        a=0.0d0
        b=0.0d0
        do i=1,3
        do j=1,nclu
         scr=dvec(i,j,nsurf,is)/mm(j)
c        a = momentum.dot.dvec/mass
c        b = dvec.dot.dvec/mass
         a=a+scr*pp(i,j)
         b=b+scr*dvec(i,j,nsurf,is)
        enddo
        enddo
        a2=a*a
        ab=a/b
c       sectrm = deltae/kinetic energy along d
        if (b.gt.0.d0) then
c       DVEC is nonzero
        sectrm = 2.0d0*deltae*b/a2
        talongd = a2/b/2.d0 ! kinetic energy along d
        tu_maxt(is) = 0.5d0/(deltae - talongd)
        radi = 1.0d0-sectrm
        else
c       DVEC is zero, mark as frustrated
        radi = -1.d0
        tu_maxt(is) = 0.d0
        endif
c       if radi > 0 then delta > kinetic energy along d, so we can hop
c       note that for a hop down, sectrm is negative, so we can always hop
        if(radi.ge.0.0d0) then
c        successful hop
         frusflag(is) = 0
        else
c        frustrated hop
         frusflag(is) = 1
        endif
       else
        frusflag(is) = 0
       endif
      enddo

      return
      end
