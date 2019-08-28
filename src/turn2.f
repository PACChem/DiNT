      subroutine turn2(mm,edia,xj,rin,rout,nsurf,symb)

c For the atom-diatom initial conditions method (INITx = 3).
c Computes the classical turning points for the energy edia and the
c rotational state xj.  From NAT8.1.

      implicit none
      include 'param.f'
      integer nsurf
      double precision mm(2),edia,rin,rout,pl,pr,rtbis_turn2,emin,
     & rmin,xj,rguess
      character*2 symb(2)

c      print *,"aj in turn2"

c to find the minimum of the effective potential, we now call diamin
      rguess=1.
      call diamin2(rmin,emin,xj,mm,nsurf,symb,rguess)

c      write(6,*)"minimum energy for j=",xj," is ",emin*27.211d0,
c     & ' at rmin=',rmin

      if (edia .lt. emin) then
        write(*,*)'diatomic energy',edia,' is less than the',
     & ' minimum vib',
     &  ' energy ',emin,' for this value of j (',xj,')'
        rin = rmin
        rout = rmin
      elseif (edia .eq. emin) then
        rin = rmin
        rout = rmin
      else
c the turning points will now be found by calling rtbis.  For a given
c energy, we need only pick values of r which are large (and small)
c enough.  Should these values fail, we print them out, and the user
c can fine tune them.
        pl =  rin
        pr =  rout
        rin =  rtbis_turn2(pl,rmin,edia,xj,mm,nsurf,symb)
        rout = rtbis_turn2(rmin,pr,edia,xj,mm,nsurf,symb)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION rtbis_turn2(x1,x2,val,xj,mmm,nsurf,symb)

C     This routine is copied from
C     W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
C     Numerical Recipes in FORTRAN, The Art of Scientific Computing,
C     Second Edition, Chapter 9, $9.1, p.347, Cambridge University
C     Press, 1994.
C     Using bisection, find the root of a function func known
C     to lie between x1 and x2. Ths root, returned as rtbis,
C     will be refined until its accuracy is +-xacc.

      implicit none
      INTEGER JMAX,NSURF,IMOL
      DOUBLE PRECISION X1,X2,XACC,VAL
      DOUBLE PRECISION dv(3),mmm(2)
C     Maximum allowed number of bisections
      PARAMETER (JMAX = 500)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid
      double precision el,dl,xj
      character*2 symb(2)

      call diapot2(x1,el,xj,mmm,nsurf,symb)
      f = el - val

      call diapot2(x2,el,xj,mmm,nsurf,symb)
      fmid= el - val

      if (f*fmid.ge.0.d0) then
          write(6,1000)
          write(6,*)x1,(f + val)
          write(6,*)x2,(fmid + val)
          write(6,*)'nsurf = ',nsurf
          return                                                        02Y22V99
      end if
      if (f.lt.0.d0) then
          rtbis_turn2=x1
          dx=x2-x1
      else
          rtbis_turn2=x2
          dx=x1-x2
      end if
      do 11 j = 1, JMAX
          dx=dx*0.5d0
          xmid=rtbis_turn2+dx

          call diapot2(xmid,el,xj,mmm,nsurf,symb)
          fmid= el - val

          if (fmid.le.0.d0) rtbis_turn2=xmid
          if (abs(dx).lt.1.d-10.or.fmid.eq.0.d0) return
   11 continue
      write(6,1010)
      stop
      return                                                            02Y22V99
C-----------------------------------------------------------------------
 1000 format(2x,'root must be bracketed in rtbis_turn2')
 1010 format(2x,'Too many bisections in rtbis_turn2')
C-----------------------------------------------------------------------
      END

