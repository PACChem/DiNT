      subroutine preatomdiatom(escat,ji,vi,arri,ri,mm,rin,rout,
     &  tau,bmax,bmin,ppreli,ecoli,nsurf,nsurft,easym,rasym)

c Precompute some information for the atom-diatom initial conditions
c (INITx = 3).

c This subroutine assumes three atoms
c INPUT
c escat = total energy
c ji = initial rot state of diatom
c vi = initial vib state of the diatom
c arri = initial arrangement
c ri = initial atom-diatom separation
c mm = masses of the atoms
c OUTPUT
c xx = nuclear coords
c pp = nuclear momenta

      implicit none
      include 'param.f'

      integer ji,vi,arri,nsurf,nsurft,ia,ib
      double precision xx(3,mnat),pp(3,mnat),mm(mnat),escat,ri,
     & einti,eroti,evibi,eelec,rmin,rin,rout,ecoli,rmass,term,
     & mdum,ppreli,bmax,bmin,tau,easym(mnsurf,3),etmp,rasym(mnsurf,3)

      write(6,*)"Asymptotic analysis"
      write(6,*)" surf  arr      Emin (eV)      Rmin (A)"
      do ib=1,3
      do ia=1,nsurft
      call diamin(rmin,etmp,-1.d0,ib,mm,ia)
      easym(ia,ib) = etmp
      rasym(ia,ib) = rmin
      write(6,50)ia,ib,etmp*autoev,rmin*autoang
      if (ia.eq.nsurf.and.ib.eq.arri) rasym(ia,ib) = rmin
      enddo
      enddo
 50   format(2i5,2f15.5)

      write(6,*)"Preparing the diatom in arrangement ",arri,
     &  "with v = ",vi," j = ",ji," surface = ",nsurf
      write(6,*)
      eelec = easym(nsurf,arri)
      write(6,*)"Energetics (relative to zero of E)"
      write(6,*)"At the classical minimum in this arrangement:"
      write(6,100)"Re    = ",rasym(nsurf,arri)*autoang," A"
      write(6,100)"Eelec = ",eelec*autoev," eV"
      call ewkb(0,vi,arri,mm,evibi,rin,rout,nsurf)
      call ewkb(ji,vi,arri,mm,einti,rin,rout,nsurf)
      eroti = einti-evibi
      ecoli = escat - einti

      write(6,100)"Escat = ",escat*autoev," eV"
      write(6,102)"Eint  = E(v=",vi,",j=",ji,") =",einti*autoev," eV"
      write(6,102)"Evib  = E(v=",vi,",j=",0,") =",evibi*autoev," eV"
      write(6,103)"Erot  = Eint  - Evib   = ",eroti*autoev," eV"
      write(6,103)"Ecol  = Escat - Eint   = ",ecoli*autoev," eV"
      write(6,101)"Turning points = ",rin*autoang," and ",
     & rout*autoang," A"
 100  format(1x,a,20x,f13.5,a)
 101  format(1x,a,f13.5,a,f13.5,a)
 102  format(1x,a,i3,a,i3,a,4x,f13.5,a)
 103  format(1x,a,3x,f13.5,a)

      if (escat.le.einti) then
        write(6,*)"Total energy ",escat*autoev," <  energy of ",
     & "initial diatomic state ",einti*autoev
        stop
      endif

      call period(einti,ji,arri,mm,rin,rout,tau,nsurf)
      write(6,103)"Period for the vibration = ",tau*autofs," fs"

c     for J=0, l=j, so the range of impact parameters is
      mdum = mm(1)+mm(2)+mm(3)
      if (arri.eq.1) rmass = mm(3)*(mm(1)+mm(2))/mdum
      if (arri.eq.2) rmass = mm(1)*(mm(2)+mm(3))/mdum
      if (arri.eq.3) rmass = mm(2)*(mm(3)+mm(1))/mdum
      term = 1.d0/dsqrt(2.d0*rmass*ecoli)
      bmax=dble(ji+1)*term
      bmin=dble(ji)*term
      write(6,101)"Range of impact parameters = ",bmin*autoang,
     &  " to ",bmax*autoang," A"

c     initial relative momentum
      ppreli = dsqrt(2.d0*ecoli*rmass)
      write(6,103)"Initial relative momentum = ",ppreli," au"
      write(6,*)

      return

      end
