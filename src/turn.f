      subroutine turn(arr,mm,edia,xj,rin,rout,nsurf)

c For the atom-diatom initial conditions method (INITx = 3).
c Computes the classical turning points for the energy edia and the
c rotational state xj.  From NAT8.1.

      implicit none
      include 'param.f'
      integer arr,nsurf
      double precision mm(mnat),edia,rin,rout,pl,pr,rtbis_turn,emin,
     & rmin,xj

c to find the minimum of the effective potential, we now call diamin
      call diamin(rmin,emin,xj,arr,mm,nsurf)

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
        pl =  0.5d0
        pr = 30.0d0
        pr = 4.0d0  ! IMLS
        rin =  rtbis_turn(pl,rmin,edia,arr,xj,mm,nsurf)
        rout = rtbis_turn(rmin,pr,edia,arr,xj,mm,nsurf)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION rtbis_turn(x1,x2,val,arr,xj,mmm,nsurf)

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
      DOUBLE PRECISION dv(3),mmm(3)
C     Maximum allowed number of bisections
      PARAMETER (JMAX = 500)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid
      integer arr
      double precision el,dl,xj

      call diapot(x1,arr,el,xj,mmm,nsurf)
      f = el - val

      call diapot(x2,arr,el,xj,mmm,nsurf)
      fmid= el - val

      if (f*fmid.ge.0.d0) then
          write(6,1000)
          write(6,*)x1,(f + val)
          write(6,*)x2,(fmid + val)
          write(6,*)'nsurf = ',nsurf,' arr = ',arr
          return                                                        02Y22V99
      end if
      if (f.lt.0.d0) then
          rtbis_turn=x1
          dx=x2-x1
      else
          rtbis_turn=x2
          dx=x1-x2
      end if
      do 11 j = 1, JMAX
          dx=dx*0.5d0
          xmid=rtbis_turn+dx

          call diapot(xmid,arr,el,xj,mmm,nsurf)
          fmid= el - val

          if (fmid.le.0.d0) rtbis_turn=xmid
          if (abs(dx).lt.1.d-10.or.fmid.eq.0.d0) return
   11 continue
      write(6,1010)
      stop
      return                                                            02Y22V99
C-----------------------------------------------------------------------
 1000 format(2x,'root must be bracketed in rtbis_turn')
 1010 format(2x,'Too many bisections in rtbis_turn')
C-----------------------------------------------------------------------
      END

