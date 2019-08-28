      subroutine atomdiatom(arr,rho,ecol,bmin,bmax,jj,
     &   rin,rout,tau,ppreli,mmm,xx,pp)

c     monte carlo selection of initial conditions for 
c     the atom-diatom initial conditions (INITx = 3)
c     Some of these routines were borrowed (with modifications)
c     from NAT8.1.

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>
c input
      integer arr,jj
      double precision rho,ecol,rin,rout,tau,ppreli
      double precision mmm(mnat),mi,mj,xx(3,mnat),pp(3,mnat)

c local
      integer i,irel,ii,ij
      real ran(5)
      double precision tin,tout,dist,rr,eta,cmom,thetap,phip,cosn
      double precision thetar,t,arg,xvec(5),tmp(12),r(3),r1(3),r2(3),
     & r3(3),rd(3)
     
      double precision phi,snt,cst,snp,csp,sne,cse,bmax,bmin,bsq
      double precision theta,mdum,rmass,xji,ppo,xxo

      do i=1,5
      xvec(i) = sprng(rng_stream)
      enddo

      xji = dble(jj)+0.5d0

        xxo=2.0d0*xvec(1)-1.0d0
        phi=2.0d0*pi*xvec(2)
        eta=2.0d0*pi*xvec(3)
        bsq = bmin**2 + xvec(4)*(bmax**2-bmin**2)
        t = tau*xvec(5)

      mdum=(mmm(1)+mmm(2)+mmm(3))
      if (arr.eq.1) rmass = mmm(3)*(mmm(1)+mmm(2))/mdum
      if (arr.eq.2) rmass = mmm(1)*(mmm(2)+mmm(3))/mdum
      if (arr.eq.3) rmass = mmm(2)*(mmm(3)+mmm(1))/mdum

        theta=acos(xxo)
        tin = t
        tout = t-tau/2.0d0                                              2/4YL93
        if (tout.ge.0.0d0) go to 10                                     2/4YL93
        rr = rin
        dist = ppreli*tin/rmass
        go to 20
   10   continue
        rr = rout
        dist = ppreli*tout/rmass
   20 continue
        snt=sin(theta)
        cst=cos(theta)
        snp=sin(phi)
        csp=cos(phi)
        sne=sin(eta)
        cse=cos(eta)
c ...   initialize coordinates
        tmp(1)=rr*snt*csp
        tmp(2)=rr*snt*snp
        tmp(3)=rr*cst
        tmp(4)=0.0d0                                                      2/4YL93
        tmp(5) = sqrt(bsq)
        tmp(6) = -sqrt(rho**2-bsq)-dist
c ...   initialize momenta
        ppo=xji/rr
        tmp(7)=-ppo*(snp*cse+csp*cst*sne)
        tmp(8)=ppo*(csp*cse-snp*cst*sne)
        tmp(9)=ppo*snt*sne
        tmp(10)=0.0d0                                                     2/4YL93
        tmp(11)=0.0d0                                                     2/4YL93
        tmp(12) = sqrt(2.d0*rmass*ecol)

c transform Jacobis to x,y,z
c tmp(4-6) and tmp(10-12) correspond to the atom-diatom relative motion
c in the new coordinates we place the center of mass of the diatom
c at the origin, so tmp(4-6) and tmp(10-12) do not need to be changed
        if (arr.eq.1) then
           irel = 3
           ii = 1
           ij = 2
        endif
        if (arr.eq.2) then
           irel = 1
           ii = 2
           ij = 3
        endif
        if (arr.eq.3) then
           irel = 2
           ii = 3
           ij = 1
        endif
        xx(1,irel) = tmp(4)
        xx(2,irel) = tmp(5)
        xx(3,irel) = tmp(6)
        pp(1,irel) = tmp(10)
        pp(2,irel) = tmp(11)
        pp(3,irel) = tmp(12)

        mi = mmm(ij)/(mmm(ii)+mmm(ij))
        mj = mmm(ii)/(mmm(ii)+mmm(ij))
        xx(1,ii) = -tmp(1)*mi
        xx(2,ii) = -tmp(2)*mi
        xx(3,ii) = -tmp(3)*mi
        xx(1,ij) = tmp(1)*mj
        xx(2,ij) = tmp(2)*mj
        xx(3,ij) = tmp(3)*mj

        pp(1,ii) = -tmp(7)-mj*tmp(10)
        pp(2,ii) = -tmp(8)-mj*tmp(11)
        pp(3,ii) = -tmp(9)-mj*tmp(12)
        pp(1,ij) = tmp(7)-mi*tmp(10)
        pp(2,ij) = tmp(8)-mi*tmp(11)
        pp(3,ij) = tmp(9)-mi*tmp(12)

        call rtox(xx,r,1)

        return

        end
