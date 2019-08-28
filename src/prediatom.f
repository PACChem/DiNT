      subroutine prediatom(vi,ji,rguess,mm,nsurf,symb,
     &    tau,rin,rout)

c Precompute some information for the atom-diatom initial conditions
c (INITx = 4).

      implicit none
      include 'param.f'

      integer nsurf,arri
      double precision mm(mnat),escat,rguess,
     & einti,eroti,evibi,eelec,rmin,rin,rout,ecoli,rmass,term,
     & mdum,ppreli,bmax,bmin,tau,etmp,
     & ji,vi
      character*2 symb(mnat)

      write(6,*)"Asymptotic analysis"
      write(6,*)"   Emin (eV)      Rmin (A)"
      call diamin2(rmin,etmp,-1.d0,mm,nsurf,symb,rguess)
      write(6,50)etmp*autoev,rmin*autoang
 50   format(2f15.5)

      write(6,*)"Preparing the diatom ",
     &  "with v = ",vi," j = ",ji," surface = ",nsurf
      write(6,*)
      eelec = etmp
      write(6,*)"Energetics (relative to zero of E)"
      write(6,*)"At the classical minimum:"
      write(6,100)"Re    = ",rmin*autoang," A"
      write(6,100)"Eelec = ",eelec*autoev," eV"
      call ewkb2(0.d0,vi,mm,evibi,rin,rout,nsurf,symb,rmin)
      call ewkb2(ji,vi,mm,einti,rin,rout,nsurf,symb,rmin)
      eroti = einti-evibi

      write(6,102)"Eint  = E(v=",vi,",j=",ji,") =",einti*autoev," eV"
      write(6,102)"Evib  = E(v=",vi,",j=",0.,") =",evibi*autoev," eV"
c     & ,2.d0*evibi*autocmi," cm-1"
      write(6,103)"Erot  = Eint  - Evib   = ",eroti*autoev," eV"
      write(6,101)"Turning points = ",rin*autoang," and ",
     & rout*autoang," A"
 100  format(1x,a,20x,f13.5,a)
 101  format(1x,a,f13.5,a,f13.5,a)
 102  format(1x,a,f6.2,a,f6.2,a,4x,f13.5,a,f13.5,a)
 103  format(1x,a,3x,f13.5,a)

      call period2(einti,ji,mm,rin,rout,tau,nsurf,symb)
      write(6,103)"Period for the vibration = ",tau*autofs," fs"
      write(6,*)

      return

      end
