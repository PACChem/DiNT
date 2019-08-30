      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension x3(maxatom),y3(maxatom),z3(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      dimension dvdx3(maxatom),dvdy3(maxatom),dvdz3(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      character*2 symb(maxatom),symb2(maxatom),symb3(maxatom)
      integer at(maxatom),at2(maxatom),at3(maxatom)

      COMMON/USROCM/ PENGYGS,PENGYES(5),
     +               PENGYIJ(15),
     +               DGSCART(25,3),DESCART(25,3,5),
     +               DIJCART(25,3,15)
C
      COMMON/USRICM/ CART(25,3),ANUZERO,
     +               NULBL(25),NFLAG(20),
     +               NASURF(5+1,5+1),NDER
C
C READ INPUT POSITIONS,DERIVATIVES, EXCITED STATE, NFLAG1-5 INFO

       n=7
       ncart=21
       natoms=8
       nder=1
       nasurf=0
       nflag(1)=0
       nflag(2)=0
       nflag(3)=0
       nflag(4)=0
       nflag(5)=0
       do i = 1,natom
         cart(i,1)=x(i)
         cart(i,2)=y(i)
         cart(i,3)=z(i)
       enddo

       if (natom.eq.2) then
         cart(1,1)= x(1) -50.38170503 / 0.52917706
         cart(1,2)= y(1) -51.47183481 / 0.52917706
         cart(1,3)= z(1) -51.84669857 / 0.52917706
         cart(2,1)= x(1) -51.47034004 / 0.52917706
         cart(2,2)= y(1) -51.47034004 / 0.52917706
         cart(2,3)= z(1) -51.95453475 / 0.52917706
         cart(3,1)= x(1) -51.93489083 / 0.52917706
         cart(3,2)= y(1) -51.47263978 / 0.52917706
         cart(3,3)= z(1) -50.96410840 / 0.52917706
         cart(4,1)= x(1) -51.78140649 / 0.52917706
         cart(4,2)= y(1) -50.57522911 / 0.52917706
         cart(4,3)= z(1) -52.50113437 / 0.52917706
         cart(5,1)= x(1) -51.78335782 / 0.52917706
         cart(5,2)= y(1) -52.36165647 / 0.52917706
         cart(5,3)= z(1) -52.50619765 / 0.52917706
         cart(6,1)= x(1)
         cart(6,2)= y(1)
         cart(6,3)= z(1)
         cart(7,1)= x(2)
         cart(7,2)= y(2)
         cart(7,3)= z(2)
       else if (natom.eq.5) then
         cart(6,1) = x(2) + 48.52965994/ 0.52917706
         cart(6,2) = y(2) + 48.52965994/ 0.52917706
         cart(6,3) = z(2) + 49.04372221/ 0.52917706
         cart(7,1) = x(2) + 48.52965994/ 0.52917706
         cart(7,2) = y(2) + 48.52965994/ 0.52917706
         cart(7,3) = z(2) + 48.07312832/ 0.52917706
       endif

       call potlibpot
 
       v=pengygs           
       do i=1,natom
         dvdx(i)=dgscart(i,1)
         dvdy(i)=dgscart(i,2)
         dvdz(i)=dgscart(i,3)
       enddo

       if (natom.eq.2) then
         do i=3,7
           dvdx(i)=0.d0
           dvdy(i)=0.d0
           dvdz(i)=0.d0
         enddo
           dvdx(1)=dgscart(6,1)
           dvdy(1)=dgscart(6,2)
           dvdz(1)=dgscart(6,3)
           dvdx(2)=dgscart(7,1)
           dvdy(2)=dgscart(7,2)
           dvdz(2)=dgscart(7,3)
       else if (natom.eq.5) then
         do i=6,7
           dvdx(i)=0.d0
           dvdy(i)=0.d0
           dvdz(i)=0.d0
         enddo
       endif

      return
      end



c**************************************************************************
C   System:           CH4OH
C   Number of electronic surfaces: 1
C   Number of derivatives: 0
C   Number of bodies: 7
C   Interface:        potlib2001
C   Common name:      CH4OH
C   Notes:
C
C   References:  J. Espinosa-Garcia, J. C. Corchado, J. Chem. Phys.,
C                Vol. 112, p. 5731, 2000.
C
C   Calling conventions:  POLYRATE interface

      SUBROUTINE POTLIBPOT
C
C        This potential is written such that:
C
C                       X(1)  - X(3)  : X, Y, Z for H1
C                       X(4)  - X(6)  : X, Y, Z for C
C                       X(7)  - X(9)  : X, Y, Z for H3
C                       X(10) - X(12) : X, Y, Z for H4
C                       X(13) - X(15) : X, Y, Z for H2
C                       X(16) - X(18) : X, Y, Z for O
C                       X(19) - X(21) : X, Y, Z for  H(O)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
      CALL CARTOU
      CALL CARTTOR
C
C Changing units from angstrom to bohr and initialize
C
      DO 10 I = 1, 21
         q(I) = R(I)*0.52918d0
         pdot(I)=0.D0
 10   CONTINUE
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0

       call coorden
c
c  calculate switching functions
c
       call switchf
c
c  calculate reference angles and their derivatives
c
       call refangles
c
c  calculate stretching potential
c
       call stretch(vstr)
c
c  calculate out of plane bending potential
c
       call opbend(vop)
c
c  calculate in plane bending potential
c
       call ipbend(vip)
c
c  total potential energy is vstr+vop+vip
c
       en=vstr+vop+vip
c
c 0.03812 conversion factor from 10(5) j/mol to hartrees
c
       en = en*0.03812D0
       ENGYGS = en
c
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
c
c 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
c
         do i=1,21
            DEGSDR(i)=pdot(i)*0.0201723d0
         enddo
         CALL RTOCART
         CALL DEDCOU
      ENDIF
C
       return
       end
c
c******************************************************
c
       subroutine coorden
c
c  calculates relative coordinates and bond lengths
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
c  calculate relative coordinates
c
       do ind=1,3
         tcb(ind)=q(nc(ind))-q(nhb(ind))
         tno(ind)=q(no(ind))-q(nhb(ind))
         do i=1,4
           tch(i,ind)=q(nc(ind))-q(nh(i,ind))
           tbh(i,ind)=q(nhb(ind))-q(nh(i,ind))
         enddo
       enddo
c
c  calculate bond lengths
c
       rcb=sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
       rno=sqrt(tno(1)*tno(1)+tno(2)*tno(2)+tno(3)*tno(3))
       do i=1,4
         rch(i)=sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+
     *                tch(i,3)*tch(i,3))
         rbh(i)=sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+
     *                tbh(i,3)*tbh(i,3))
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine refangles
c
c  subroutine calculates reference angles for the "in-plane" potential
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       tau=acos(-1.0d0/3.0d0)
C       pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
c
c  set diagonal elements to zero
c
       do i=1,4
         theta0(i,i)=0.0d0
         do k=1,4
           dtheta0(i,i,k)=0.0d0
         enddo
       enddo
c
c  calculate reference angles
c
       theta0(1,2)=tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
       theta0(1,3)=tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
       theta0(1,4)=tau+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
       theta0(2,3)=tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
       theta0(2,4)=tau+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
       theta0(3,4)=tau+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c
c  calculate the derivatives of theta0(i,j) in terms of rch(k)
c  quantity calulated is dtheta0(i,j,k)
c
c  derivatives wrt rch(1)
c
       dtheta0(1,2,1)=(tau-halfpi)*dsphi(1)*sphi(2)
       dtheta0(1,3,1)=(tau-halfpi)*dsphi(1)*sphi(3)
       dtheta0(1,4,1)=(tau-halfpi)*dsphi(1)*sphi(4)
       dtheta0(2,3,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(4)
       dtheta0(2,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(3)
       dtheta0(3,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(2)
c
c  derivatives wrt rch(2)
c
       dtheta0(1,2,2)=(tau-halfpi)*sphi(1)*dsphi(2)
       dtheta0(1,3,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(4)
       dtheta0(1,4,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(3)
       dtheta0(2,3,2)=(tau-halfpi)*dsphi(2)*sphi(3)
       dtheta0(2,4,2)=(tau-halfpi)*dsphi(2)*sphi(4)
       dtheta0(3,4,2)=(tau-twopi/3.0d0)*stheta(1)*dstheta(2)
c
c  derivatives wrt rch(3)
c
       dtheta0(1,2,3)=(tau-twopi/3.0d0)*dstheta(3)*stheta(4)
       dtheta0(1,3,3)=(tau-halfpi)*sphi(1)*dsphi(3)
       dtheta0(1,4,3)=(tau-twopi/3.0d0)*stheta(2)*dstheta(3)
       dtheta0(2,3,3)=(tau-halfpi)*sphi(2)*dsphi(3)
       dtheta0(2,4,3)=(tau-twopi/3.0d0)*stheta(1)*dstheta(3)
       dtheta0(3,4,3)=(tau-halfpi)*dsphi(3)*sphi(4)
c
c  derivatives wrt rch(4)
c
       dtheta0(1,2,4)=(tau-twopi/3.0d0)*stheta(3)*dstheta(4)
       dtheta0(1,3,4)=(tau-twopi/3.0d0)*stheta(2)*dstheta(4)
       dtheta0(1,4,4)=(tau-halfpi)*sphi(1)*dsphi(4)
       dtheta0(2,3,4)=(tau-twopi/3.0d0)*stheta(1)*dstheta(4)
       dtheta0(2,4,4)=(tau-halfpi)*sphi(2)*dsphi(4)
       dtheta0(3,4,4)=(tau-halfpi)*sphi(3)*dsphi(4)
c
c  fill in the other half of the matrix
c
        do i=1,3
          do j=i+1,4
            theta0(j,i)=theta0(i,j)
            do k=1,4
              dtheta0(j,i,k)=dtheta0(i,j,k)
            enddo
          enddo
        enddo
       return
       end
c
c******************************************************
c
c
       subroutine stretch(vstr)
c
c  subroutine to calculate leps-type stretching potential and its
c  derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       dimension vqch(4),vjch(4),vqbh(4),vjbh(4),vq(4),vj(4),
     *           achdc(3),achdh(4,3)
c
c  calculate avergage bond length for the methane moiety
c
       rav=(rch(1)+rch(2)+rch(3)+rch(4))/4.0d0
c
c  initialise:
c
       vstr=0.0d0
c
c  calculate d3cb and dd3cb
c
       texp = exp(-(4.d0*(rav-rcbsp)/b3cb)**4.d0)
       d3cb = (d3cbi-a3cb) + a3cb * texp
       dd3cb = -4.d0*a3cb*texp*(rav-rcbsp)**3.d0*(4.d0/b3cb)**4.d0
c
c  ach:
c
c  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
c  in to avoid overflow/underflow errors
c
       arga=c1ch*(rav-r0ch)
       if(arga.lt.19.0d0)then
         ach=a1ch+b1ch*(tanh(arga)+1.0d0)*0.5d0
         dumach=b1ch*c1ch/(2.0d0*cosh(arga)**2)
       else
         ach=a1ch+b1ch
         dumach=0.0d0
       endif
c
c  calculate singlet: e1, triplet: e3 energies and vq and vj
c  terms for each bond
c
       e1=d1cb*(exp(-2.0d0*acb*(rcb-r0cb))-2.0d0*exp(-acb*(rcb-r0cb)))
       e3=d3cb*(exp(-2.0d0*acb*(rcb-r0cb))+2.0d0*exp(-acb*(rcb-r0cb)))
       vqcb=(e1+e3)*0.5d0
       vjcb=(e1-e3)*0.5d0
c
c O-H term (simple morse term)
c
       dt=(rno-r0hh)
       expterm=exp(-ahh*dt)
       vno=d1hh*(1.d0-expterm)**2.d0
c
       do i=1,4
         e1=d1ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              -2.0d0*exp(-ach*(rch(i)-r0ch)))
         e3=d3ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              +2.0d0*exp(-ach*(rch(i)-r0ch)))
         vqch(i)=(e1+e3)*0.5d0
         vjch(i)=(e1-e3)*0.5d0
         e1=d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         e3=d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         vqbh(i)=(e1+e3)*0.5d0
         vjbh(i)=(e1-e3)*0.5d0
c
c  calculate 3 body potential
c
         vq(i)=vqch(i)+vqcb+vqbh(i)
         vj(i)=-sqrt(((vjch(i)-vjcb)**2+(vjcb-vjbh(i))**2
     *                 +(vjbh(i)-vjch(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
c  add OH contribution
c
       vstr=vstr+vno
c
c  partial derivatives
c
c  derivative of OH morse term with respect to dt:
c
       deddt = 2.0d0 * ahh * d1hh * (1.0d0-expterm) *
     *         expterm
c
c  chain rule: derivatives with respect to cartesian coordinates.
c
       de = deddt / rno
       do i=1,3
          ded(i) = de * tno(i)
       enddo
c
c  now we need the derivative of ach:
c
       do ind=1,3
         achdc(ind)=dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2)
     *            +tch(3,ind)/rch(3)+tch(4,ind)/rch(4))/4.0d0
         do i=1,4
           achdh(i,ind)=-dumach*tch(i,ind)/rch(i)/4.0d0
         enddo
       enddo
       dumqcb=-acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))-
     *         (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
c
c  calculate cartesian derivatives:
c  looping over ch(i) and bh(i)
c
       do i=1,4
         dumqbh=-ahh*((d1hh+d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))-
     *           (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
c
c  adding the derv of D3cb wrt r
c
         addd3 = dd3cb*(exp(-2.0d0*acb*(rcb-r0cb))+
     *         2.0d0*exp(-acb*(rcb-r0cb)))
         addd3 = addd3/4.d0 * 0.5d0 / rcb
         dumqbh = dumqbh + addd3
c
         factj=0.5d0/vj(i)
         dumjcb=-acb*((d1cb-d3cb)*exp(-2.0d0*acb*(rcb-r0cb))
     *            -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
         dumjbh=-ahh*((d1hh-d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *            -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))*factj/rbh(i)
c
c  adding derv of D3cb wrt r
c
         addd3 = dd3cb*(exp(-2.0d0*acb*(rcb-r0cb))+
     *           2.0d0*exp(-acb*(rcb-r0cb)))
         addd3 = addd3/4.d0 * 0.5d0 * factj / rcb
         dumjbh = dumjbh - addd3
c
         do ind=1,3
c
c  deriv wrt hb:
c
                  pdot(nhb(ind))=pdot(nhb(ind))
     *             -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh
     *            +(vjch(i)-vjcb)*(dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind))
     *            +(vjbh(i)-vjch(i))*dumjbh*tbh(i,ind)
c
c  dvqch(i)/dc
c
           dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
               pdot(nc(ind))=pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
c
c  dvqch(i)/dh(i)
c
           dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
              pdot(nh(i,ind))=pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
c
c  dvjch(i)/dc
c
           dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnc(ind)
c
           pdot(nc(ind))=pdot(nc(ind))
     *            +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*dumjcb*tcb(ind)
     *            -(vjbh(i)-vjch(i))*dumjch
c
c  dvjch(i)/dh(i)
c
           dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnh(i,ind)
c
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *            +(vjch(i)-vjcb)*dumjhi
     *            +(vjcb-vjbh(i))*dumjbh*tbh(i,ind)
     *            +(vjbh(i)-vjch(i))*(-dumjbh*tbh(i,ind)-dumjhi)
c
c  dv(i)/dh(j)
c
           do k=1,3
             j=i+k
             if(j.gt.4)j=j-4
             dumqhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                    -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
             dumjhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                  -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
             pdot(nh(j,ind))=pdot(nh(j,ind))+dumqhj
     *            +(vjch(i)-vjcb)*dumjhj
     *            -(vjbh(i)-vjch(i))*dumjhj
           enddo
         enddo
       enddo
c
c add derivatives of OH morse term
c
       do ind=1,3
          pdot(nhb(ind))=pdot(nhb(ind)) - ded(ind)
          pdot(no(ind))=pdot(no(ind)) + ded(ind)
       enddo
c
       return
       end
c
c******************************************************
c
c
       subroutine opbend(vop)
c
c  subroutine calculates symmetrized vop potential and derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       double precision norma
       dimension sumd2(4),sumd4(4)
       dimension in(3),a(3),b(3),axb(3),c(4,3),argd(4)
c
c
       vop=0.0d0
c
c  calculate force constants and their derivatives
c
       call opforce
c
c  calculate out-of-plane angle and derivatives
c
       do i=1,4
         j=i+1
         if(j.gt.4)j=j-4
         k=j+1
         if(k.gt.4)k=k-4
         l=k+1
         if(l.gt.4)l=l-4
c
c  modification to ensure that the set of methane CH bond vector
c  (rj,rk,rl) is a right-handed set
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(k,ind))-q(nh(j,ind))
         b(ind)=q(nh(l,ind))-q(nh(j,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
c
c  if argd > 0 we need to switch vectors k and l around
c
         if (argd(in(ii)).gt.0.d0) then
             itemp=k
             k=l
             l=itemp
         endif
       enddo
c
c  subroutine performs sum over j, k, l
c  sum2 = sum delta**2
c  sum4 = sum delta**4
c
         call calcdelta(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
       enddo
       do i=1,4
         do j=1,4
c
c  overall derivatives of force constants i wrt the bond-length rch(j)
c
           ddr=dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
c
c  calculate derivatives in terms of cartesian coordinates:
c
           do ind=1,3
             pdot(nh(j,ind))=pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
             pdot(nc(ind))=pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
           enddo
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipbend(vip)
c
c  subroutine calculates symmetrised in plane bend term
c  and its derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       dimension costh(4,4),theta(4,4),dth(4,4)
c
c  initialise
c
       vip=0.0d0
c
c  calculate force constants: fk0(i,j), f1(i)
c  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
c
       call ipforce
c
c  calculate theta(i,j) and in plane bend potential
c
       do i=1,3
         do j=i+1,4
           costh(i,j)=tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2)
     *                       +tch(i,3)*tch(j,3)
           costh(i,j)=costh(i,j)/rch(i)/rch(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
c
c  calculate partial derivatives wrt cartesian coordinates
c
c  calculate pdots wrt theta:
c
           termth=-1.0d0/sqrt(1.0d0-costh(i,j)*costh(i,j))
           do ind=1,3
             dthi=-tch(j,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
             dthi=dthi*termth
             dthj=-tch(i,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
             dthj=dthj*termth
             dthc=-(dthi+dthj)
             pdot(nh(i,ind))=pdot(nh(i,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
             pdot(nh(j,ind))=pdot(nh(j,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
             pdot(nc(ind))=pdot(nc(ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
             do k=1,4
c
c  calculate pdots wrt force constants and wrt theta0
c
               dth0k=-dtheta0(i,j,k)*tch(k,ind)/rch(k)
               dth0c=-dth0k
               pdot(nh(k,ind))=pdot(nh(k,ind))
     *                -0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                -0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0k*dth(i,j)
               pdot(nc(ind))=pdot(nc(ind))
     *                +0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0c*dth(i,j)
               pdot(nhb(ind))=pdot(nhb(ind))
     *                +0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
             enddo
           enddo
         enddo
       enddo
c
c  Now, calculate h2o bending energies.
c
c  First, calculate the angles:
c
       do i = 1,4
         dot = 0.d0
         do j = 1,3
           dot = dot - tno(j)*tbh(i,j)
         enddo
         cosine = dot / (rno*rbh(i))
         cosine = min(1.d0, max(-1.d0, cosine))
         angh2o(i) = acos(cosine)
       enddo
c
c  Now, calculate each force constant fkh2o(i) as a function of the O-H(I)
c  distance, rbh(i)
c
       do i=1,4
          arga= alph2o * (rbh(i) - r0hh)
          if(arga.lt.19.0d0)then
            fkh2o(i) = fkh2oeq * (1 - tanh(arga))
          else
            fkh2o(i) = 0.0d0
          endif
       enddo
c
c  Now calculate the contribution to the energy for each H(I)-O-H(O)
c  harmonic bending term
c
      do i=1,4
           dang = (angh2o(i) - anh2oeq)
           vip=vip+0.5d0* fkh2o(i) * dang * dang
      enddo
c
c  Calculate the derivatives of the bending hamonic terms
c  (in order to make the code more efficient, some of these terms
c  can be calculated at the same time as the terms for the energy,
c  but in order to make it more clear we will keep them separated
c  even though it means recalculating some stuff)
c
c  d vip                                   d fkh2o(i)    d r
c ------ = 1/2 * (angh2o(i)-anh2oeq)**2 * ---------- * -----  +
c  d x                                       d r         d x
c
c                                            d angh2o(i)
c         + fkh2o(i) * (angh2o(i)-anh2oeq) * -----------
c                                               d x
c
c
c  First, calculate d fkh2o(i) / d r
c
      do i=1,4
         deno = cosh(alph2o * (rbh(i) - r0hh))
         dkdr(i) = - (fkh2oeq * alph2o) / (deno*deno)
      enddo
c
c  Now we calculate d fhk2o(i) / d x
c
      do i=1,4
         dkdr(i) = dkdr(i)/rbh(i)
         do j=1,3
            dkdx(i,j) = dkdr(i) * tbh(i,j)
         enddo
      enddo
c
c  We add the contribution of the first half of d vip / d x to pdot
c
      do i=1,4
         do j=1,3
            term1 = 0.5d0 * (angh2o(i)-anh2oeq)**2.d0
            pdot(nhb(j)) = pdot(nhb(j)) + dkdx(i,j) * term1
            pdot(nh(i,j)) = pdot(nh(i,j)) - dkdx(i,j) * term1
         enddo
      enddo
c
c  Now calculate the second half:
c
      do i=1,4
        dstda = fkh2o(i) * (angh2o(i)-anh2oeq)
c
c  contributions of d angle / d x adapted from the tinker code (eangle1.f)
c
        xp = tno(2)*tbh(i,3) - tno(3)*tbh(i,2)
        yp = tno(3)*tbh(i,1) - tno(1)*tbh(i,3)
        zp = tno(1)*tbh(i,2) - tno(2)*tbh(i,1)
        rp = sqrt (xp*xp + yp*yp + zp*zp)
c
c  check of lineality
c
        if (rp.lt.1.d-6) rp=1.d-6
c
        terma = dstda / (rbh(i)*rbh(i)*rp)
        termc = dstda / (rno*rno*rp)
        dedhi(1) = - terma * (tbh(i,2)*zp-tbh(i,3)*yp)
        dedhi(2) = - terma * (tbh(i,3)*xp-tbh(i,1)*zp)
        dedhi(3) = - terma * (tbh(i,1)*yp-tbh(i,2)*xp)
        dedho(1) = - termc * (tno(2)*zp-tno(3)*yp)
        dedho(2) = - termc * (tno(3)*xp-tno(1)*zp)
        dedho(3) = - termc * (tno(1)*yp-tno(2)*xp)
        dedo(1) = -dedhi(1) - dedho(1)
        dedo(2) = -dedhi(2) - dedho(2)
        dedo(3) = -dedhi(3) - dedho(3)
c
c  add these contributions to pdot
c  add these contributions to pdot
c
        do j=1,3
           pdot(nhb(j)) = pdot(nhb(j)) + dedo(j)
           pdot(nh(i,j)) = pdot(nh(i,j)) + dedhi(j)
           pdot(no(j)) = pdot(no(j)) + dedho(j)
        enddo
      enddo
c
       return
       end
c
c*************************************************************************
c
       subroutine calcdelta(i,j,k,l,sum2,sum4)
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       dimension  delta(4),in(3),a(3),b(3),axb(3),c(4,3),argd(4),
     *            daxb(4,3,3),cdot(4,3,3),atemp2(3)
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  set j,k,l indices
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(k,ind))-q(nh(j,ind))
         b(ind)=q(nh(l,ind))-q(nh(j,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
c        write(*,*) 'theta,delta ',theta0(i,in(ii)),delta(in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
c  derivatives of axb wrt hj:
c
       daxb(j,1,1)=0.0d0
       daxb(j,1,2)=b(3)-a(3)
       daxb(j,1,3)=-b(2)+a(2)
       daxb(j,2,1)=-b(3)+a(3)
       daxb(j,2,2)=0.0d0
       daxb(j,2,3)=b(1)-a(1)
       daxb(j,3,1)=b(2)-a(2)
       daxb(j,3,2)=-b(1)+a(1)
       daxb(j,3,3)=0.0d0
c
c  derivatives of axb wrt hk:
c
       daxb(k,1,1)=0.0d0
       daxb(k,1,2)=-b(3)
       daxb(k,1,3)=b(2)
       daxb(k,2,1)=b(3)
       daxb(k,2,2)=0.0d0
       daxb(k,2,3)=-b(1)
       daxb(k,3,1)=-b(2)
       daxb(k,3,2)=b(1)
       daxb(k,3,3)=0.0d0
c
c  derivatives of axb wrt hl:
c
       daxb(l,1,1)=0.0d0
       daxb(l,1,2)=a(3)
       daxb(l,1,3)=-a(2)
       daxb(l,2,1)=-a(3)
       daxb(l,2,2)=0.0d0
       daxb(l,2,3)=a(1)
       daxb(l,3,1)=a(2)
       daxb(l,3,2)=-a(1)
       daxb(l,3,3)=0.0d0
c
c   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
c   wrt h(in(ii),jind) with components jind
c
       do ii=1,3
c
c  deriv of cdot(in(ii),x) wrt x, y, z
c
         cdot(in(ii),1,1)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,2)=tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,3)=tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),y) wrt x, y, z
c
         cdot(in(ii),2,1)=tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,2)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,3)=tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),z) wrt x, y, z
c
         cdot(in(ii),3,1)=tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,2)=tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,3)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
       enddo
c
       do ii=1,3
         do ind=1,3
            deldot=-dtheta0(i,in(ii),i)
c
c  derivative wrt h(i,ind)
c  for  rch(i) only terms are from the derivatives of theta0
c
            deldot=-deldot*tch(i,ind)/rch(i)
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  derivative wrt c(ind)
c
            deldot=-deldot
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
           do jj=1,3
c
c  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
c
c   atemp1 is axb dot daxb wrt h(in(jj))
c
            atemp1=axb(1)*daxb(in(jj),ind,1)
     *            +axb(2)*daxb(in(jj),ind,2)
     *            +axb(3)*daxb(in(jj),ind,3)
            atemp1=atemp1/(norma**3)
c
c  atemp2 is deriv of normalised axb
c
            atemp2(1)=daxb(in(jj),ind,1)/norma-atemp1*axb(1)
            atemp2(2)=daxb(in(jj),ind,2)/norma-atemp1*axb(2)
            atemp2(3)=daxb(in(jj),ind,3)/norma-atemp1*axb(3)
c
c  atemp3 is daxb dot c(in(ii))

            atemp3=atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2)
     *                             +atemp2(3)*c(in(ii),3)
c
c  atemp4 is axb dot cdot
c
            atemp4=0.0d0
            if(ii.eq.jj)then
c
c  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
c
              atemp4=axb(1)*cdot(in(ii),1,ind)
     *                     +axb(2)*cdot(in(ii),2,ind)
     *                     +axb(3)*cdot(in(ii),3,ind)
              atemp4=atemp4/norma
            endif
c
c  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
c
            atemp5=-dtheta0(i,in(ii),in(jj))
c
c  deriv wrt h(in(jj)),ind):
c
            atemp5=-atemp5*tch(in(jj),ind)/rch(in(jj))
            deldot=atemp3+atemp4
            deldot=-1.0d0/sqrt(1.0d0-argd(in(ii))**2)*deldot
            deldot=deldot+atemp5
            pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  for carbon the only contributions are from axb dot cdot term and
c  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
c
            deldot=1.0d0/sqrt(1.0d0-argd(in(ii))**2)*atemp4
            deldot=deldot-atemp5
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
          enddo
        enddo
       enddo
       return
       end

c******************************************************
c
       subroutine opforce
c
c  calculates the out-of-plane bending force constants
c  and their derivatives
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       dimension switch(4),dswitch(4,4)
c
c  calculate switching functions:
c
       switch(1)=(1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
       switch(2)=(1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
       switch(3)=(1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
       switch(4)=(1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
c
c  calculate derivatives:
c  derivative of switch(1) wrt the 4 rch bond lengths
c
       dswitch(1,1)=-ds3(1)*s3(2)*s3(3)*s3(4)
       dswitch(1,2)=(1.0d0-s3(1))*ds3(2)*s3(3)*s3(4)
       dswitch(1,3)=(1.0d0-s3(1))*s3(2)*ds3(3)*s3(4)
       dswitch(1,4)=(1.0d0-s3(1))*s3(2)*s3(3)*ds3(4)
c
c  derivative of switch(2) wrt the 4 rch bond lengths
c
       dswitch(2,1)=(1.0d0-s3(2))*s3(3)*s3(4)*ds3(1)
       dswitch(2,2)=-ds3(2)*s3(3)*s3(4)*s3(1)
       dswitch(2,3)=(1.0d0-s3(2))*ds3(3)*s3(4)*s3(1)
       dswitch(2,4)=(1.0d0-s3(2))*s3(3)*ds3(4)*s3(1)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(3,1)=(1.0d0-s3(3))*s3(4)*ds3(1)*s3(2)
       dswitch(3,2)=(1.0d0-s3(3))*s3(4)*s3(1)*ds3(2)
       dswitch(3,3)=-ds3(3)*s3(4)*s3(1)*s3(2)
       dswitch(3,4)=(1.0d0-s3(3))*ds3(4)*s3(1)*s3(2)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(4,1)=(1.0d0-s3(4))*ds3(1)*s3(2)*s3(3)
       dswitch(4,2)=(1.0d0-s3(4))*s3(1)*ds3(2)*s3(3)
       dswitch(4,3)=(1.0d0-s3(4))*s3(1)*s3(2)*ds3(3)
       dswitch(4,4)=-ds3(4)*s3(1)*s3(2)*s3(3)
c
c  calculate the force constants and their derivatives
c
       do i=1,4
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
         do j=1,4
           dfdelta(i,j)=dswitch(i,j)*fch3
           dhdelta(i,j)=dswitch(i,j)*hch3
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipforce
c
c  calculates the symmetrised in plane bend force constants and
c  all partial derivatives involving them
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       dimension dfk0(4,4,4),df1dc(4),df1dh(4)
c
c  set force constant at asymptotes
c
       f0=fkinf+ak
       f2=fkinf
c
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
       fk0(1,4)=f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
       fk0(2,4)=f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
       fk0(3,4)=f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
c
c  derivative of fk0
c
       dfk0(1,2,1)=f0*ds1(1)*s1(2)
       dfk0(1,2,2)=f0*s1(1)*ds1(2)
       dfk0(1,2,3)=(f0-f2)*ds2(3)*s2(4)
       dfk0(1,2,4)=(f0-f2)*s2(3)*ds2(4)
c
       dfk0(1,3,1)=f0*ds1(1)*s1(3)
       dfk0(1,3,2)=(f0-f2)*ds2(2)*s2(4)
       dfk0(1,3,3)=f0*s1(1)*ds1(3)
       dfk0(1,3,4)=(f0-f2)*s2(2)*ds2(4)
c
       dfk0(1,4,1)=f0*ds1(1)*s1(4)
       dfk0(1,4,2)=(f0-f2)*ds2(2)*s2(3)
       dfk0(1,4,3)=(f0-f2)*s2(2)*ds2(3)
       dfk0(1,4,4)=f0*s1(1)*ds1(4)
c
       dfk0(2,3,1)=(f0-f2)*ds2(1)*s2(4)
       dfk0(2,3,2)=f0*ds1(2)*s1(3)
       dfk0(2,3,3)=f0*s1(2)*ds1(3)
       dfk0(2,3,4)=(f0-f2)*s2(1)*ds2(4)
c
       dfk0(2,4,1)=(f0-f2)*ds2(1)*s2(3)
       dfk0(2,4,2)=f0*ds1(2)*s1(4)
       dfk0(2,4,3)=(f0-f2)*s2(1)*ds2(3)
       dfk0(2,4,4)=f0*s1(2)*ds1(4)
c
       dfk0(3,4,1)=(f0-f2)*ds2(1)*s2(2)
       dfk0(3,4,2)=(f0-f2)*s2(1)*ds2(2)
       dfk0(3,4,3)=f0*ds1(3)*s1(4)
       dfk0(3,4,4)=f0*s1(3)*ds1(4)
c
       do i=1,4
c
c  calculate the terms f1(i)
c
         arga1=aa1*rbh(i)*rbh(i)
         arga2=aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rch(i)-r0ch)**2)
c
c  and calculate the derivatives wrt rch(i) and rbh(i)
c
         duma1=2.0d0*aa1*rbh(i)*exp(-arga1)
         duma2=-2.0d0*aa3*aa4*(rbh(i)-r0hh)*exp(-arga2)
         df1dc(i)=-2.0d0*(rch(i)-r0ch)*a1*a2*exp(-a2*(rch(i)-r0ch)**2)
         df1dh(i)=duma1*exp(-a2*(rch(i)-r0ch)**2)
     *             -duma2*(rch(i)-r0ch)**2*a1*exp(-a2*(rch(i)-r0ch)**2)
       enddo
c
c  derivative of total force constant f(i,j) wrt bond length rch(k)
c  is given by dfdc(i,j,k)
c
      dfdc(1,2,1)=dfk0(1,2,1)*f1(1)*f1(2)+fk0(1,2)*df1dc(1)*f1(2)
      dfdc(1,2,2)=dfk0(1,2,2)*f1(1)*f1(2)+fk0(1,2)*f1(1)*df1dc(2)
      dfdc(1,2,3)=dfk0(1,2,3)*f1(1)*f1(2)
      dfdc(1,2,4)=dfk0(1,2,4)*f1(1)*f1(2)
c
      dfdc(1,3,1)=dfk0(1,3,1)*f1(1)*f1(3)+fk0(1,3)*df1dc(1)*f1(3)
      dfdc(1,3,2)=dfk0(1,3,2)*f1(1)*f1(3)
      dfdc(1,3,3)=dfk0(1,3,3)*f1(1)*f1(3)+fk0(1,3)*f1(1)*df1dc(3)
      dfdc(1,3,4)=dfk0(1,3,4)*f1(1)*f1(3)
c
      dfdc(1,4,1)=dfk0(1,4,1)*f1(1)*f1(4)+fk0(1,4)*df1dc(1)*f1(4)
      dfdc(1,4,2)=dfk0(1,4,2)*f1(1)*f1(4)
      dfdc(1,4,3)=dfk0(1,4,3)*f1(1)*f1(4)
      dfdc(1,4,4)=dfk0(1,4,4)*f1(1)*f1(4)+fk0(1,4)*f1(1)*df1dc(4)
c
      dfdc(2,3,1)=dfk0(2,3,1)*f1(2)*f1(3)
      dfdc(2,3,2)=dfk0(2,3,2)*f1(2)*f1(3)+fk0(2,3)*df1dc(2)*f1(3)
      dfdc(2,3,3)=dfk0(2,3,3)*f1(2)*f1(3)+fk0(2,3)*f1(2)*df1dc(3)
      dfdc(2,3,4)=dfk0(2,3,4)*f1(2)*f1(3)
c
      dfdc(2,4,1)=dfk0(2,4,1)*f1(2)*f1(4)
      dfdc(2,4,2)=dfk0(2,4,2)*f1(2)*f1(4)+fk0(2,4)*df1dc(2)*f1(4)
      dfdc(2,4,3)=dfk0(2,4,3)*f1(2)*f1(4)
      dfdc(2,4,4)=dfk0(2,4,4)*f1(2)*f1(4)+fk0(2,4)*f1(2)*df1dc(4)
c
      dfdc(3,4,1)=dfk0(3,4,1)*f1(3)*f1(4)
      dfdc(3,4,2)=dfk0(3,4,2)*f1(3)*f1(4)
      dfdc(3,4,3)=dfk0(3,4,3)*f1(3)*f1(4)+fk0(3,4)*df1dc(3)*f1(4)
      dfdc(3,4,4)=dfk0(3,4,4)*f1(3)*f1(4)+fk0(3,4)*f1(3)*df1dc(4)
c
c  derivative of total force constant f(i,j) wrt bond length rbh(k)
c  is given by dfdh(i,j,k)
c
c  only non-zero derivatives are those from rbh(i) and rbh(j)
c
       dfdh(1,2,1)=fk0(1,2)*df1dh(1)*f1(2)
       dfdh(1,2,2)=fk0(1,2)*f1(1)*df1dh(2)
       dfdh(1,2,3)=0.0d0
       dfdh(1,2,4)=0.0d0
c
       dfdh(1,3,1)=fk0(1,3)*df1dh(1)*f1(3)
       dfdh(1,3,2)=0.0d0
       dfdh(1,3,3)=fk0(1,3)*f1(1)*df1dh(3)
       dfdh(1,3,4)=0.0d0
c
       dfdh(1,4,1)=fk0(1,4)*df1dh(1)*f1(4)
       dfdh(1,4,2)=0.0d0
       dfdh(1,4,3)=0.0d0
       dfdh(1,4,4)=fk0(1,4)*f1(1)*df1dh(4)
c
       dfdh(2,3,1)=0.0d0
       dfdh(2,3,2)=fk0(2,3)*df1dh(2)*f1(3)
       dfdh(2,3,3)=fk0(2,3)*f1(2)*df1dh(3)
       dfdh(2,3,4)=0.0d0
c
       dfdh(2,4,1)=0.0d0
       dfdh(2,4,2)=fk0(2,4)*df1dh(2)*f1(4)
       dfdh(2,4,3)=0.0d0
       dfdh(2,4,4)=fk0(2,4)*f1(2)*df1dh(4)
c
       dfdh(3,4,1)=0.0d0
       dfdh(3,4,2)=0.0d0
       dfdh(3,4,3)=fk0(3,4)*df1dh(3)*f1(4)
       dfdh(3,4,4)=fk0(3,4)*f1(3)*df1dh(4)
c
       return
       end
c
c******************************************************
c
c
       subroutine switchf
c
c  calculates switching functions: s3,sphi,stheta
c  and their derivatives ds3,dsphi,dstheta
c
       implicit double precision (a-h,o-z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
       a1s=1.5313681d-7
       b1s=-4.6696246d0
       a2s=1.0147402d-7
       b2s=-12.363798d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
       do i=1,4
         args1=a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
         if(args1.lt.argmax)then
           s1(i)=1.0d0-tanh(args1)
           ds1(i)=a1s*((rch(i)-b1s)**8
     *                 +8.0d0*(rch(i)-r0ch)*(rch(i)-b1s)**7)
           ds1(i)=-ds1(i)/cosh(args1)**2
         else
           s1(i)=0.0d0
           ds1(i)=0.0d0
         endif
c
         args2=a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
         if(args2.lt.argmax)then
           s2(i)=1.0d0-tanh(args2)
           ds2(i)=a2s*((rch(i)-b2s)**6
     *                 +6.0d0*(rch(i)-r0ch)*(rch(i)-b2s)**5)
           ds2(i)=-ds2(i)/cosh(args2)**2
         else
           s2(i)=0.0d0
           ds2(i)=0.0d0
         endif
c
c  calculate s3 and ds3
c
         args3=a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i)=1.0d0-tanh(args3)
           ds3(i)=a3s*(3.0d0*rch(i)**2-2.0d0*rch(i)*(r0ch+2.0d0*b3s)
     *          +b3s*(b3s+2.0d0*r0ch))
           ds3(i)=-ds3(i)/cosh(args3)**2
         else
           s3(i)=0.0d0
           ds3(i)=0.0d0
         endif
c
c  calculate sphi and dsphi
c
c  condition here is on the bondlength rch(i)
c  st argsphi is lt approx 19.0d0
c
         if(rch(i).lt.3.8d0)then
           argsphi=aphi*(rch(i)-r0ch)*exp(bphi*(rch(i)-cphi)**3)
           sphi(i)=1.0d0-tanh(argsphi)
           dsphi(i)=aphi*(1.0d0+3.0d0*bphi*(rch(i)-r0ch)
     *                      *(rch(i)-cphi)**2)
           dsphi(i)=dsphi(i)*exp(bphi*(rch(i)-cphi)**3)
           dsphi(i)=-dsphi(i)/cosh(argsphi)**2
         else
           sphi(i)=0.0d0
           dsphi(i)=0.0d0
         endif
c
c  calculate stheta and dstheta
c
         if(rch(i).lt.3.8d0)then
           argstheta=atheta*(rch(i)-r0ch)*exp(btheta*(rch(i)-ctheta)**3)
           stheta(i)=1.0d0-tanh(argstheta)
           dstheta(i)=atheta*(1.0d0+3.0d0*btheta*(rch(i)-r0ch)
     *           *(rch(i)-ctheta)**2)
           dstheta(i)=dstheta(i)*exp(btheta*(rch(i)-ctheta)**3)
           dstheta(i)=-dstheta(i)/cosh(argstheta)**2
         else
           stheta(i)=0.0d0
           dstheta(i)=0.0d0
         endif
       enddo
       return
       end
C
      SUBROUTINE PREPOT
C
C   N3TMMN = 3 * NATOMS
C   NATOMS = the number of atoms represented by this potential function
C
C   The variable N3TMMN is the minimum value of N3TM allowed to be
C   passed by the calling routine for the number of cartesian
C   coordinates needed to represent the full system represented by this
C   potential energy surface routine.
C   N3TM must be greater than or equal to N3TMMN.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
      PARAMETER (N3TMMN = 21)
C
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
C
C
C  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
C
      IF(NATOMS.GT.25) THEN
         WRITE(NFLAG(18),1111)
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')
         STOP
      END IF
C
       DO I=1,5
          REF(I) = ' '
       END DO
C
       REF(1)='J. Espinosa-Garcia and J. C. Corchado'
       REF(2)='J. Chem. Phys., Vol. 112, p. 5731, 2000'
C
      INDEXES(1) = 1
      INDEXES(2) = 6
      INDEXES(3) = 1
      INDEXES(4) = 1
      INDEXES(5) = 1
      INDEXES(6) = 8
      INDEXES(7) = 1
C
C
C
      IRCTNT=6
C
      CALL POTINFO
C
      CALL ANCVRT
c
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc(ind)=3*nnc+icount
         nhb(ind)=3*nnb+icount
         do i=1,4
           nh(i,ind)=3*nnh(i)+icount
           no(ind)=3*nno+icount
         enddo
       enddo
c
c  convert to appropriate units:
c
c  energy   in 1.0d+05 j/mol
c  time     in 1.0d-14 s
c  distance in 1.0d-10 m
c  angles   in radians
c  mass     in amu
c
       fact1=0.041840d0
       fact2=6.022045d0
       fact3=2.d0*PI/360.d0
c
      d1ch=d1ch*fact1
      d3ch=d3ch*fact1
      d1cb=d1cb*fact1
      d3cbi=d3cbi*fact1
      a3cb=a3cb*fact1
      d1hh=d1hh*fact1
      d3hh=d3hh*fact1
      fch3=fch3*fact2
      hch3=hch3*fact2
      fkinf=fkinf*fact2
      ak=ak*fact2
      fkh2oeq=fkh2oeq*fact2
      anh2oeq=anh2oeq*fact3
C
      RETURN
      END
C
      BLOCK DATA PTPACM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(5),nno,
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,
     +               fkh2oeq,alph2o,anh2oeq,
     +               a3cb,b3cb,rcbsp
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rno,rch(4),rbh(4)
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
       common /addh2o/ded(3),angh2o(4),fkh2o(4),
     +                dkdr(4),dkdx(4,3),dedo(3),dedhi(3),dedho(3)
C
      DATA NASURF /1,35*0/
      DATA NDER /1/
       DATA NFLAG /1,1,15*0,6,0,0/
C
      DATA ANUZERO /0.0D0/
      DATA ICARTR,MSURF,MDER/1,0,1/
      DATA NULBL /25*0/
      DATA NATOMS /7/
C
       DATA nnc     /2/
       DATA nnb     /6/
       DATA nnh     /3,4,5,1,0/
       DATA nno     /7/
       DATA r0ch    /   1.09397d0/
       DATA d1ch    / 112.17000d0/
       DATA d3ch    /  32.65328d0/
       DATA a1ch    /   1.78000d0/
       DATA b1ch    /   0.15000d0/
       DATA c1ch    /  15.00000d0/
       DATA r0hh    /   0.97060d0/
       DATA d1hh    / 125.44000d0/
       DATA d3hh    /  20.41017d0/
       DATA ahh     /   2.15000d0/
       DATA r0cb    /   1.49092d0/
       DATA d1cb    /  91.47526d0/
       DATA d3cbi   / 112.69509d0/
       DATA acb     /   2.98688d0/
       DATA a3s     /   0.1419100d0/
       DATA b3s     /  -0.3068400d0/
       DATA aphi    /   0.5287900d0/
       DATA bphi    /   0.4006600d0/
       DATA cphi    /   1.9209900d0/
       DATA atheta  /   0.9078700d0/
       DATA btheta  /   0.3548900d0/
       DATA ctheta  /   1.8915500d0/
       DATA fch3    /   0.0740000d0/
       DATA hch3    /   0.1915000d0/
       DATA fkinf   /   0.4400000d0/
       DATA ak      /   0.1260000d0/
       DATA bk      /  10.7132d0/
       DATA aa1     /   0.303746d0/
       DATA aa2     /   1.599960d0/
       DATA aa3     /   3.216595d0/
       DATA aa4     /  11.569980d0/
       DATA fkh2oeq /   0.7300000d0/
       DATA alph2o  /   1.1080000d0/
       DATA anh2oeq / 104.7132000d0/
       DATA a3cb    /   0.000000d0/
       DATA b3cb    /   0.900000d0/
       DATA rcbsp   /   2.606485d0/
C
       END


      SUBROUTINE POTINFO

C***********************************************************************
C   VERSION HISTORY
C
C   Version 1.0
C   Author:              R. J. Duchovic 
C   Date:                99-10-15   
C                        99-11-05   
C                        99-11-09   
C                        99-11-10   
C                        99-11-17   
C                        99-11-19   
C                        99-12-09   
C                        99-12-10   
C                        99-12-14   
C                        01-07-18   
C                                  
C   Version 1.0.1
C   Corrections by:      Jose Corchado
C   Date:                May 13, 2003
C   Summary of changes:  This new version of utility.f corrects two
C                        errors in the original version.  
C                        Corrected lines in are marked 'JC05132003' in 
C                        the margin.  
C                        The corrections are:
C                        1) In the part of the subroutine CARTTOR where
C                           the code handles ICARTR = 2, the counter 
C                           "I" was added.  Previously, not all of the 
C                           components of R(I) were getting the correct
C                           values.
C                        2) In the part of the subroutine RTOCART where
C                           the code handles ICARTR = 2, the DO loop
C                           was changed to run from 1 to NATOMS.  
C                           Previously, some of the derivatives with 
C                           resepect to the Cartesian coordinates were 
C                           incorrectly returned as zero.
C
C***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      PARAMETER (PI = 3.141592653589793D0)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
C LIST CONTENTS OF REF
C
      write(NFLAG(18),96)
 96   format(/)
      do i =1,5
         write(NFLAG(18),97) REF(i)
      end do
 97   format(2x,a75)
C
C LIST AND CHECK ON NUMBER OF SURFACES
C
      WRITE(NFLAG(18),96)
      KMAX = 0
      DO I = 1,ISURF+1
         DO J = 1,ISURF+1
            IF(NASURF(I,J).NE.0.AND.KMAX.LT.MAX(I,J)) KMAX = MAX(I,J)
         ENDDO
      ENDDO
      WRITE(NFLAG(18),101) MSURF,KMAX-1
101   FORMAT(2x,' MAX. AND ACTUAL NO. OF EXCITED SURFACES: ',I3,5x,I3)
      IF(KMAX-1.GT.MSURF) THEN
         WRITE(6,*) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
         STOP
      ENDIF
      KSDIAG = 0
      KEDIAG = 0
      DO I = 2,ISURF+1
         IF(NASURF(I,I).NE.0) THEN
            KEDIAG = I-1
            IF(KSDIAG.EQ.0) KSDIAG = I-1
         ENDIF
      ENDDO
      KSOFFD = 0
      KEOFFD = 0
      K = 0
      DO I = 1,ISURF
         DO J = I+1,ISURF+1
            K = K+1
            IF(NASURF(I,J)+NASURF(J,I).NE.0) THEN
               KEOFFD = K
               IF(KSOFFD.EQ.0) KSOFFD = K
            ENDIF
         ENDDO
      ENDDO
C
C LIST AND CHECK ON ORDER OF DERIVATIVES
C
      WRITE(NFLAG(18),103) MDER,NDER
103   FORMAT(2x,' MAX. AND ACTUAL ORDER OF DERIVATIVES:    ',I3,5x,I3)
      IF(NDER.GT.MDER) THEN
         WRITE(6,*) ' WRONG INPUT ON ORDER OF DERIVATIVES'
         STOP
      ENDIF
C
C**********************************************
C                                             *
C   NFLAG(19) is used to determine whether    *
C   or not a detailed description of the      *
C   options is printed.                       *
C                                             *
C**********************************************
C
      IF(NFLAG(19).EQ.1) THEN
         write(NFLAG(18),100)
 100     format(/)
         write(NFLAG(18),120)
 120     format(2x,'Cartesian coordinates are supplied by',/,
     +          2x,'the user in the array CART.',//)
         write(NFLAG(18),125)
 125     format(2x,'Provide cartesian coordinates in the',/,
     +          2x,'following order using the array CART',//,
     +          2x,' CART(1,1)...CART(1,3)   => ATOM 1',/,
     +          2x,' CART(2,1)...CART(2,3)   => ATOM 2',/,
     +          2x,' CART(3,1)...CART(3,3)   => ATOM 3',/,
     +          2x,' CART(N,1)...CART(N,3)   => ATOM N',/,
     +          2x,'CART(25,1)...CART(25,3)  => ATOM 25',/)
C
C**********************************************
C                                             *
C        Cartesian Coordinates                *
C                                             *
C   CART(1,1)  ... CART(1,3)  => ATOM 1       *
C   CART(2,1)  ... CART(2,3)  => ATOM 2       *
C   CART(3,1)  ... CART(3,3)  => ATOM 3       *
C       .               .                     *
C       .               .                     *
C       .               .                     *
C   CART(N,1)  ... CART(N,3)  => ATOM N       *
C       .               .                     *
C       .               .                     *
C       .               .                     *
C   CART(25,1) ... CART(25,3) => ATOM 25      *
C                                             *
C**********************************************
C
         write(NFLAG(18),130)
 130     format(2x,'If the user wishes to relabel the atoms,',/,
     +          2x,'set the variable IREORDER equal to 1',/,
     +          2x,'in the PARAMETER statement.  The user',/,
     +          2x,'must also specify the new labeling',/,
     +          2x,'scheme.  This is done using the array',/,
     +          2x,'NULBL in the following manner:',//,
     +          2x,'NULBL(i) = j',/,
     +          2x,'where:  i => old index',/,
     +          2x,'        j => new index',//)
         write(NFLAG(18),150)
 150     format(2x,'Cartesian coordinates can be provided to',/,
     +          2x,'the potential routine in a variety of units.',/,
     +          2x,'The input units will be converted to Bohr',/,
     +          2x,'based on the following values of the NFLAG',/,
     +          2x,'variable:',//,
     +          2x,'NFLAG(1)  =  1  =>  CARTESIANS IN BOHR (no',/,
     +          2x,'                    conversion required)',/,
     +          2x,'NFLAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS',//)
         write(NFLAG(18),160)
 160     format(2x,'The value of the energy and derivatives',/,
     +          2x,'(if computed) can be reported in a variety',/,
     +          2x,'units.  A units conversion will take place',/,
     +          2x,'as specified by the following values of the',/,
     +          2x,'NFLAG variable:',//,
     +          2x,'NFLAG(2) = 1 =>  ENERGIES REPORTED IN HARTEEE',/,
     +          2x,'NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTREE',/,
     +          2x,'NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV',/,
     +          2x,'NFLAG(2) = 4 =>  ENERGIES REPORTED IN kcal/mol',/,
     +          2x,'NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1',//)
         write(NFLAG(18),165)
 165     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/BOHR',//)
         write(NFLAG(18),170)
 170     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/ANGSTROM',//)
      ENDIF
      RETURN
      END
 
 
 
 
C
      SUBROUTINE ANCVRT
C
C***********************************
C                                  *
C   WRITTEN BY: R. J. DUCHOVIC     *
C               A. F. WAGNER       *
C                                  *            
C   DATE:               99-10-15   *
C                       99-11-09   *
C                       99-11-10   *
C                       99-11-15   *
C                       99-11-16   *
C                       99-11-17   *
C                       99-11-19   *
C                       99-12-09   *
C                       00-01-23   *
C                       01-07-18   *
C                                  *
C***********************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*75 REF(5)
      CHARACTER*3 PERIODIC_1(7,32)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      CHARACTER*2 NAME1(NATOM)
      CHARACTER*2 NAME2(NATOM)
      CHARACTER*1 IBLANK
      CHARACTER*20 DISTANCE
      CHARACTER*20 UNITS
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION IANUM(7,32)
      DIMENSION ISAVE(NATOM),JSAVE(NATOM)
C
C***********************************************************************
C                                                                      *
C  PI:                          3.141592653589793                      *
C                                                                      *
C  FUNDAMENTAL  CONSTANTS (1998 CODATA values)                         *
C                                                                      *
C  SPPED OF LIGHT:              2.99792458 x 10**8   ms**-1            *
C  PERMEABILITY  OF VACUUM:     4.0D0*PI*1.0 x 10**-07  NA**-2         *
C  PERMITIVITY  OF VACUUM:      1.0/(CMU_0*CLIGHT**2) Fm**-1           *
C  ELEMENTARY  CHARGE:          1.602176462 x 10**-19  C               *
C  PLANCK'S CONSTANT            6.62606876 x 10**-34  Js               *
C  ELECTRON MASS:               9.10938188 x 10**-31  kg               *        
C  ANGSTROM:                    1.0 X 10**-10 m                        *
C  AVOGADRO'S CONSTANT:         6.02214199 x 10**23 mol**-1            *
C  KILOCALORIE:                 4.184  x 10**10 ergs                   *
C                                                                      *
C  BOHR     = 2.0D0*CPLANCK**2/(2.0D0*PI*CM_E*CLIGHT**2*CMU_0*CE**2)   *
C           = 5.291772083 x 10**-11 m                                  *
C  HTOMILLH = 1000 milliHartree                                        *
C  HTOEV    = CE/(4.0D0*PI*CEPSILON_0*BOHR)                            *
C           = 27.2113834 eV                                            *
C  HTOKCAL  = (CE**2*CAVOGADRO)/(4.0D0*PI*CEPSILON_0*BOHR*CKCAL)       *
C           = 627.509470 kcal/mole                                     *
C  HTOWAVE  = (CE**2)/(4.0D0*PI*CEPSILON_0*BOHR*CPLANCK*CLIGHT)        *
C           = 219474.631 Hartree                                       *
C  HTOKJ    = HTOKCAL*4.184                                            *
C           = 2625.49962                                               *
C                                                                      *
C***********************************************************************
C
      PARAMETER(        PI = 3.141592653589793D0)
      PARAMETER(    CLIGHT = 2.99792458D08)
      PARAMETER(     CMU_0 = 4.0D0*PI*1.0D-07)
      PARAMETER(CEPSILON_0 = 1.0D0/(CMU_0*CLIGHT**2))
      PARAMETER(        CE = 1.602176462D-19)
      PARAMETER(   CPLANCK = 6.62606876D-34)
      PARAMETER(      CM_E = 9.10938188D-31)
      PARAMETER(      CANG = 1.0D-10)
      PARAMETER( CAVOGADRO = 6.02214199D23)
      PARAMETER(     CKCAL = 4.184D10)
      PARAMETER(  HTOMILLH = 1000.D0)
      PARAMETER(     HTOEV = 27.2113834D0)
      PARAMETER(   HTOKCAL = 627.509470D0)
      PARAMETER(   HTOWAVE = 219474.631D0)
      PARAMETER(     HTOKJ = 2625.49962D0)
      PARAMETER(    BOHR_A = .5291772083D0)
C
      DO I=1,7
         DO J=1,32
            IANUM(I,J)=0
            PERIODIC_1(I,J)=' '
         END DO
      END DO
C
      DISTANCE = 'BOHR                '
      UNITS    = 'HARTREE             '
C
      IANUM(1,1)  =  1
      IANUM(1,32) =  2
      IANUM(2,1)  =  3
      IANUM(2,2)  =  4
      IANUM(2,27) =  5
      IANUM(2,28) =  6
      IANUM(2,29) =  7
      IANUM(2,30) =  8
      IANUM(2,31) =  9
      IANUM(2,32) = 10
      IANUM(3,1)  = 11
      IANUM(3,2)  = 12
      IANUM(3,27) = 13
      IANUM(3,28) = 14
      IANUM(3,29) = 15
      IANUM(3,30) = 16
      IANUM(3,31) = 17
      IANUM(3,32) = 18
      IANUM(4,1)  = 19
      IANUM(4,2)  = 20
      IANUM(4,17) = 21
      IANUM(4,18) = 22
      IANUM(4,19) = 23
      IANUM(4,20) = 24
      IANUM(4,21) = 25
      IANUM(4,22) = 26
      IANUM(4,23) = 27
      IANUM(4,24) = 28
      IANUM(4,25) = 29
      IANUM(4,26) = 30
      IANUM(4,27) = 31
      IANUM(4,28) = 32
      IANUM(4,29) = 33
      IANUM(4,30) = 34
      IANUM(4,31) = 35
      IANUM(4,32) = 36
      IANUM(5,1)  = 37
      IANUM(5,2)  = 38
      IANUM(5,17) = 39
      IANUM(5,18) = 40
      IANUM(5,19) = 41
      IANUM(5,20) = 42
      IANUM(5,21) = 43
      IANUM(5,22) = 44
      IANUM(5,23) = 45
      IANUM(5,24) = 46
      IANUM(5,25) = 47
      IANUM(5,26) = 48
      IANUM(5,27) = 49
      IANUM(5,28) = 50
      IANUM(5,29) = 51
      IANUM(5,30) = 52
      IANUM(5,31) = 53
      IANUM(5,32) = 54
      IANUM(6,1)  = 55
      IANUM(6,2)  = 56
      IANUM(6,3)  = 57
      IANUM(6,4)  = 58
      IANUM(6,5)  = 59
      IANUM(6,6)  = 60
      IANUM(6,7)  = 61
      IANUM(6,8)  = 62
      IANUM(6,9)  = 63
      IANUM(6,10) = 64
      IANUM(6,11) = 65
      IANUM(6,12) = 66
      IANUM(6,13) = 67
      IANUM(6,14) = 68
      IANUM(6,15) = 69
      IANUM(6,16) = 70
      IANUM(6,17) = 71
      IANUM(6,18) = 72
      IANUM(6,19) = 73
      IANUM(6,20) = 74
      IANUM(6,21) = 75
      IANUM(6,22) = 76
      IANUM(6,23) = 77
      IANUM(6,24) = 78
      IANUM(6,25) = 79
      IANUM(6,26) = 80
      IANUM(6,27) = 81
      IANUM(6,28) = 82
      IANUM(6,29) = 83
      IANUM(6,30) = 84
      IANUM(6,31) = 85
      IANUM(6,32) = 86
      IANUM(7,1)  = 87
      IANUM(7,2)  = 88
      IANUM(7,3)  = 89
      IANUM(7,4)  = 90
      IANUM(7,5)  = 91
      IANUM(7,6)  = 92
      IANUM(7,7)  = 93
      IANUM(7,8)  = 94
      IANUM(7,9)  = 95
      IANUM(7,10) = 96
      IANUM(7,11) = 97
      IANUM(7,12) = 98
      IANUM(7,13) = 99
      IANUM(7,14) = 100
      IANUM(7,15) = 101
      IANUM(7,16) = 102
      IANUM(7,17) = 103
      IANUM(7,18) = 104
      IANUM(7,19) = 105
      IANUM(7,20) = 106
      IANUM(7,21) = 107
      IANUM(7,22) = 108
      IANUM(7,23) = 109
      IANUM(7,24) = 110
      IANUM(7,25) = 111
      IANUM(7,26) = 112
      IANUM(7,27) = 113
      IANUM(7,28) = 114
      IANUM(7,29) = 115
      IANUM(7,30) = 116
      IANUM(7,31) = 117
      IANUM(7,32) = 120
C
      PERIODIC_1(1,1)   = 'H  '
      PERIODIC_1(1,32)  = 'He '
      PERIODIC_1(2,1)   = 'Li '
      PERIODIC_1(2,2)   = 'Be '
      PERIODIC_1(2,27)  = 'B  '
      PERIODIC_1(2,28)  = 'C  '
      PERIODIC_1(2,29)  = 'N  '
      PERIODIC_1(2,30)  = 'O  '
      PERIODIC_1(2,31)  = 'F  '
      PERIODIC_1(2,32)  = 'Ne '
      PERIODIC_1(3,1)   = 'Na '
      PERIODIC_1(3,2)   = 'Mg '
      PERIODIC_1(3,27)  = 'Al '
      PERIODIC_1(3,28)  = 'Si '
      PERIODIC_1(3,29)  = 'P  '
      PERIODIC_1(3,30)  = 'S  '
      PERIODIC_1(3,31)  = 'Cl '
      PERIODIC_1(3,32)  = 'Ar '
      PERIODIC_1(4,1)   = 'K  '
      PERIODIC_1(4,2)   = 'Ca '
      PERIODIC_1(4,17)  = 'Sc '
      PERIODIC_1(4,18)  = 'Ti '
      PERIODIC_1(4,19)  = 'V  '
      PERIODIC_1(4,20)  = 'Cr '
      PERIODIC_1(4,21)  = 'Mn '
      PERIODIC_1(4,22)  = 'Fe '
      PERIODIC_1(4,23)  = 'Co '
      PERIODIC_1(4,24)  = 'Ni '
      PERIODIC_1(4,25)  = 'Cu '
      PERIODIC_1(4,26)  = 'Zn '
      PERIODIC_1(4,27)  = 'Ga '
      PERIODIC_1(4,28)  = 'Ge '
      PERIODIC_1(4,29)  = 'As '
      PERIODIC_1(4,30)  = 'Se '
      PERIODIC_1(4,31)  = 'Br '
      PERIODIC_1(4,32)  = 'Kr '
      PERIODIC_1(5,1)   = 'Rb '
      PERIODIC_1(5,2)   = 'Sr '
      PERIODIC_1(5,17)  = 'Y  '
      PERIODIC_1(5,18)  = 'Zr '
      PERIODIC_1(5,19)  = 'Nb '
      PERIODIC_1(5,20)  = 'Mo '
      PERIODIC_1(5,21)  = 'Tc '
      PERIODIC_1(5,22)  = 'Ru '
      PERIODIC_1(5,23)  = 'Rh '
      PERIODIC_1(5,24)  = 'Pd '
      PERIODIC_1(5,25)  = 'Ag '
      PERIODIC_1(5,26)  = 'Cd '
      PERIODIC_1(5,27)  = 'In '
      PERIODIC_1(5,28)  = 'Sn '
      PERIODIC_1(5,29)  = 'Sb '
      PERIODIC_1(5,30)  = 'Te '
      PERIODIC_1(5,31)  = 'I  '
      PERIODIC_1(5,32)  = 'Xe '
      PERIODIC_1(5,32)  = 'Xe '
C
C
      DO I=1,NATOMS
         ISAVE(I)=0
         JSAVE(I)=0
         NAME1(I)='  '
         NAME2(I)='  '
      END DO
      IBLANK=' '
C
      DO IND=1,NATOMS
         DO I=1,7
            DO J=1,32
               IF(INDEXES(IND).EQ.IANUM(I,J)) THEN
                  ISAVE(IND)=I
                  JSAVE(IND)=J
               END IF
            END DO
         END DO
      END DO
      write(NFLAG(18),100)
 100  format(/,2x,'The potential routine in POTLIB assumes a',/,
     +         2x,'default labeling of the atoms.  The user has',/,
     +         2x,'made the following selection of atomic labels:',/)
      write(NFLAG(18),125)
 125  format(15x,'DEFAULT',25x,'USER-SELECTED',/)
 
      DO IND=1,NATOMS
         IND2=NULBL(IND)
         IF(IND2.EQ.0) IND2=IND
         WRITE(NFLAG(18),1000) IND,
     +                 PERIODIC_1(ISAVE(IND),JSAVE(IND)),
     +                 PERIODIC_1(ISAVE(IND2),JSAVE(IND2))
 1000    FORMAT(2x,'Atom ',I2,' is ',6x,a2,34x,a2)
      END DO
      WRITE(NFLAG(18),1500)
 1500 FORMAT(/)
C
C*******************************
C                              *
C   Determine Reactant #1      *
C                              *
C*******************************
C
      INC1=0
      DO IND=1,IRCTNT-1
         INC1=INC1+1
         NAME1(INC1)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
         WRITE(NFLAG(18),2000) IND,
     +                    PERIODIC_1(ISAVE(IND),JSAVE(IND))
 2000    FORMAT(2x,'Atom ',I2,2x,a2,' is a member of Reactant #1')
      END DO
C
C*******************************
C                              *
C   Determine Reactant #2      *
C                              *
C*******************************
C
      INC2=0
      DO IND=IRCTNT,NATOMS
         INC2=INC2+1
         NAME2(INC2)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
         WRITE(NFLAG(18),2005) IND,
     +                    PERIODIC_1(ISAVE(IND),JSAVE(IND))
 2005    FORMAT(2x,'Atom ',I2,2x,a2,' is a member of Reactant #2')
      END DO
      write(NFLAG(18),140)
 140  format(2x,/,2x,'The Potential Energy is Zero for the following',/,
     +       2x,'geometric arrangement:  The reactants are at',/,
     +       2x,'their equilibrium geometries, widely separated.',/)
      write(NFLAG(18),145)
 145  format(2x,/,2x,'Reactant #1')
      write(NFLAG(18),*) IBLANK,(NAME1(I),I=1,INC1)
      write(NFLAG(18),147)
 147  format(2x,/,2x,'Reactant #2')
      write(NFLAG(18),*) IBLANK,(NAME2(I),I=1,INC2)
      write(NFLAG(18),148)
 148  format(2x,/,2x,'The default units are:',/,
     +       5x,'BOHR FOR DISTANCE',/,
     +       5x,'HARTREE FOR ENERGY',/)
C
      IF(NFLAG(1).EQ.2) DISTANCE = 'ANGSTROMS           '
      IF(NFLAG(2).EQ.2) THEN
         UNITS = 'MILLIHARTREE        '
      ELSEIF(NFLAG(2).EQ.3) THEN
         UNITS = 'EV                  '
      ELSEIF(NFLAG(2).EQ.4) THEN
         UNITS = 'KCAL PER MOLE       '
      ELSEIF(NFLAG(2).EQ.5) THEN
         UNITS = 'WAVENUMBERS         '
      ELSEIF(NFLAG(2).EQ.6) THEN
         UNITS = 'KILOJOULES PER MOLE '
      ENDIF
      write(NFLAG(18),149) DISTANCE,UNITS
 149  format(2x,/,2x,'The user has chosen units :',/,
     +       5x,a20,1x,'FOR DISTANCE',/,
     +       5x,a20,1x,'FOR ENERGY',/)
C
      write(NFLAG(18),150)
 
 150  format(2x,/,2x,'To reset the zero-of-energy, assign a',/,
     +            2x,'non-zero value to the variable ANUZERO.',/,
     +            2x,'The value of ANUZERO will be SUBTRACTED from',/,
     +            2x,'the energy calculated by this routine.',/)
C
C*********************************************
C                                            *
C   SET UP CONSTANTS FOR UNITS CONVERSIONS   *
C                                            *
C*********************************************
C
C*****************************************************************
C                                                                *
C   NFLAG(1) = 1    => UNITS IN BOHR, NO CONVERSION REQUIRED     *
C   NFLAG(1) = 2    => UNITS IN ANGSTROMS, CONVERSION REQUIRED   *
C                                                                *
C   NFLAG(2) = 1 => UNITS IN HARTREE, NO CONVERSION REQUIRED     *
C            = 2 => UNITS IN MILLIHARTREE, CONVERSION REQUIRED   *
C            = 3 => UNITS IN EV, CONVERSION REQUIRED             *
C            = 4 => UNITS IN KCLAL, CONVERSION REQUIRED          *
C            = 5 => UNITS IN WAVENUMBER, CONVERSION REQUIRED     *
C            = 6 => UNITS IN KILOJOULES/MOLE, CONVERSION REQUIRED*
C                                                                *
C*****************************************************************
C
      CNVRTD = 1.D0
      CNVRTE = 1.D0
      CNVRTDE = 1.D0
C
C**************************
C                         *
C   DISTANCE CONVERSION   *
C                         *
C**************************
C
      IF(NFLAG(1).EQ.2) CNVRTD = BOHR_A
C
C**************************
C                         *
C   ENERGY CONVERSION     *
C                         *
C**************************
C
      IF(NFLAG(2).EQ.2) THEN
         CNVRTE = CNVRTE*HTOMILLH
      ELSEIF(NFLAG(2).EQ.3) THEN
         CNVRTE = CNVRTE*HTOEV
      ELSEIF(NFLAG(2).EQ.4) THEN
         CNVRTE = CNVRTE*HTOKCAL
      ELSEIF(NFLAG(2).EQ.5) THEN
         CNVRTE = CNVRTE*HTOWAVE
      ELSEIF(NFLAG(2).EQ.6) THEN
         CNVRTE = CNVRTE*HTOKJ
      ENDIF
C
C******************************
C                             *
C   DERIVATIVE CONVERSION     *
C                             *
C******************************
C
      CNVRTDE = CNVRTE/CNVRTD
C
C*********************************
C                                *
C   ASSIGN VALUE TO "IREORDER"   *
C                                *
C*********************************
C
      ISUM = 0
      DO INU=1,25
         ISUM=ISUM + NULBL(INU)
      END DO
      IREORDER = 0
      IF(ISUM.NE.0) IREORDER = 1
C
      RETURN
      END
C
      SUBROUTINE CARTOU
C
C********************************************************
C                                                       *
C   WRITTEN BY: R. J. DUCHOVIC,A. F. WAGNER             *
C                                                       *
C   DATE:               00-01-22                        *
C                       01-07-18                        *
C                                                       *
C   THIS SUBROUTINE IS DESIGNED TO CONVERT              *
C   POSITIONS CART IN USER ORDER AND UNITS INTO         *
C   POSITIONS CARTNU IN DEVELOPER ORDER AND UNITS       *
C   (DEVELOPER UNITS ARE FIXED AT BOHR).                *
C   CONVERSION FACTOR DONE IN ANCVRT.F                  *
C                                                       *
C********************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
c      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
c     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
c     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C**********************************************
C                                             *
C   The following code will convert units     *
C   and relabel the                           *
C   atoms to match the the numbering          *
C   convention chosen by the user             *
C                                             *
C   If the user wishes to relabel the         *
C   atoms set the variable IREORDER equal     *
C   to 1 in the PARAMETER statement.  The     *
C   user must also specify the new labeling   *
C   scheme.  This is done using the array     *
C   NULBL in the following manner:            *
C                                             *
C   NULBL(i) = j                              *
C                                             *
C      where:  i => old index                 *
C              j => new index                 *
C                                             *
C**********************************************
C
      IF (IREORDER.EQ.1) THEN
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(NULBL(I),J)=CART(I,J)/CNVRTD
             END DO
          END DO
      ELSE
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(I,J)=CART(I,J)/CNVRTD
             END DO
          END DO
      END IF
C
      RETURN
      END
 
C
      SUBROUTINE CARTTOR
C
C*****************************************************************
C                                                                *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                     *
C                                                                *
C   DATE:               00-01-22                                 *
C                       00-08-03                                 *
C                       00-08-11                                 *
C                       01-07-18                                 * 
C                                                                *
C   THIS SUBROUTINE IS DESIGNED TO CONVERT FROM                  *
C   POSITIONS CARTNU IN DEVELOPER UNITS AND ORDER                *
C   TO DIFFERENT SETS OF INTERNAL COORDINATES R                  *
C   AS CONTROLLED BY ICARTR                                      *
C                                                                *
C   CURRENT ICARTR OPTIONS ARE:                                  *
C                                                                *
C   ICARTR = 0                                                   *
C     NOTHING HAPPENS (R VECTOR IS NOT FILLED)                   *
C                                                                *
C   ICARTR = 1                                                   *
C     ARRAY CARTNU WILL BE STORED IN THE VECTOR R AS             *
C      R(1) ....... R(3)  CARTESIANS (X,Y,Z) ATOM 1              *
C      R(4) ....... R(6)  CARTESIANS (X,Y,Z) ATOM 2              *
C      R(7) ....... R(9)  CARTESIANS (X,Y,Z) ATOM 3              *
C       .            .                                           *
C      R(3*N-2) ... R(3*N)  CARTESIANS (X,Y,Z) ATOM N            *
C       .            .                                           *
C      R(73) ...... R(75)  CARTESIANS (X,Y,Z) ATOM 25            *
C                                                                *
C   ICARTR = 2                                                   *
C     ARRAY CARTNU WILL BE CONVERTED INTO DISTANCE AS            *
C      R(1) .................... DISTANCE OF 1 TO 2              *
C      R(2) .................... DISTANCE OF 1 TO 3              *
C       .                        .                               *
C      R(NATOMS) ............... DISTANCE OF 1 TO NATOMS         *
C      R(NATOMS+1) ............. DISTANCE OF 2 TO 3              *
C      R(NATOMS+2) ............. DISTANCE OF 2 TO 4              *
C       .                        .                               *
C      R(2*NATOMS-1) ........... DISTANCE OF 2 TO NATOMS         *
C      R(2*NATOMS)   ........... DISTANCE OF 3 TO 4              *
C       .                        .                               *
C      R(NATOMS*(NATOMS-1)/2 ... DISTANCE OF NATOMS-1 TO NATOMS  *
C                                                                *
C   ICARTR = 3                                                   *
C     !!!FOR NATOMS = 3 ONLY !!!                                 *
C      R(1) .................... DISTANCE OF 1 TO 2              *
C      R(2) .................... DISTANCE OF 2 TO 3              *
C      R(3) .................... DISTANCE OF 1 TO 3              *
C                                                                *
C   ICARTR = 4                                                   *
C     Transformation for PES requiring Distances and Angles      *
C      HF Dimer PESs                                             *
C      R(1) .................... Center-of-Mass DISTANCE         *
C      R(2) .................... Bond DISTANCE First Diatom      *
C      R(3) .................... Bond DISTANCE Second Diatom     *
C      R(4) .................... Theta one angle                 *
C      R(5) .................... Theta two angle                 *
C      R(6) .................... Dihedral Tau angle              *
C                                                                *
C   NATOMS  =>  NUMBER OF ATOMS                                  *
C                                                                *
C*****************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C************************************************************
C                                                           *
C   SELECT CARTNU -> R OPTION WITH ICARTR                   *
C                                                           *
C************************************************************
C
      IF(ICARTR.EQ.1) THEN
C
C************************************************************
C                                                           *
C   R = VECTOR MAP OF CARTNU ARRAY                          *
C                                                           *
C************************************************************
C
         DO I=1,NATOMS
            IND=3*I-2
            R(IND)   = CARTNU(I,1)
            R(IND+1) = CARTNU(I,2)
            R(IND+2) = CARTNU(I,3)
         END DO
      ELSEIF(ICARTR.EQ.2) THEN
C
C************************************************************
C                                                           *
C   R = INTERPARTICLE DISTANCES IN CANONICAL ORDER          *
C                                                           *
C   R(1) 1-2 DISTANCE                                       *
C   R(2) 1-3 DISTANCE                                       *
C   ....                                                    *
C   R(2) 2-3 DISTANCE                                       *
C   ....                                                    *
C                                                           *
C************************************************************
C
         I = 1
cJC05132003
         DO K=1,NATOMS-1
            DO L = K+1,NATOMS
cJC05132003
               R(I) = SQRT( (CARTNU(K,1)-CARTNU(L,1))**2 +
     +                      (CARTNU(K,2)-CARTNU(L,2))**2 +
     +                      (CARTNU(K,3)-CARTNU(L,3))**2 )
               I = I + 1
cJC05132003
            END DO
         ENDDO
      ELSEIF(ICARTR.EQ.3) THEN
C
C************************************************************
C                                                           *
C   !!! FOR NATOMS = 3 ONLY!!!                              *
C   R = INTERPARTICLE DISTANCES IN THE ORDER                *
C       R(1)=1-2 DIST., R(2)=2-3 DIST., R(3)=1-3 DIST.      *
C                                                           *
C************************************************************
C
         R(1) = SQRT( (CARTNU(1,1)-CARTNU(2,1))**2 +
     +                (CARTNU(1,2)-CARTNU(2,2))**2 +
     +                (CARTNU(1,3)-CARTNU(2,3))**2 )
         R(2) = SQRT( (CARTNU(2,1)-CARTNU(3,1))**2 +
     +                (CARTNU(2,2)-CARTNU(3,2))**2 +
     +                (CARTNU(2,3)-CARTNU(3,3))**2 )
         R(3) = SQRT( (CARTNU(1,1)-CARTNU(3,1))**2 +
     +                (CARTNU(1,2)-CARTNU(3,2))**2 +
     +                (CARTNU(1,3)-CARTNU(3,3))**2 )
C
      ELSEIF(ICARTR.EQ.4) THEN
C
C*************************************************************
C                                                            *
C   ICARTR = 4                                               *
C     Transformation for PES requiring Distances and Angles  *
C      HF Dimer PESs                                         *
C                                                            *
C*************************************************************
C
C      COMMON//V,XB(30),DXB(30)
C      COMMON/SFCOM/IFLAG1,ICOR,ISFCT
C
C*********************************************
C                                            *
C   THE CARTESIAN COORDINATES WERE READ IN   *
C   Calculate INTERNAL COORDINATES           * 
C                                            *
C   FIND CENTER-OF-MASS OF EACH MONOMER      *
C                                            *
C   MASS OF HYDROGEN AND OF FLUORINE         *
C                                            *
C*********************************************
C
      FLM=18.99840D0
      HYM=1.007825D0
C
C************************************************************************
C                                                                       *
C   CENTER OF MASS OF FIRST MONOMER, DO X, Y, AND Z COORDINATE IN TURN
C   *
C                                                                       *
C************************************************************************
C
C      XCM1=(HYM*XB(1)+FLM*XB(4))/(FLM+HYM)
C      YCM1=(HYM*XB(2)+FLM*XB(5))/(FLM+HYM)
C      ZCM1=(HYM*XB(3)+FLM*XB(6))/(FLM+HYM)
C
      XCM1=(HYM*CARTNU(1,1)+FLM*CARTNU(2,1))/(FLM+HYM)
      YCM1=(HYM*CARTNU(1,2)+FLM*CARTNU(2,2))/(FLM+HYM)
      ZCM1=(HYM*CARTNU(1,3)+FLM*CARTNU(2,3))/(FLM+HYM)
C
C****************************************
C                                       *
C   CENTER OF MASS OF SECOND MONOMER    *
C                                       *
C****************************************
C
C      XCM2=(HYM*XB(7)+FLM*XB(10))/(FLM+HYM)
C      YCM2=(HYM*XB(8)+FLM*XB(11))/(FLM+HYM)
C      ZCM2=(HYM*XB(9)+FLM*XB(12))/(FLM+HYM)
C
      XCM2=(HYM*CARTNU(3,1)+FLM*CARTNU(4,1))/(FLM+HYM)
      YCM2=(HYM*CARTNU(3,2)+FLM*CARTNU(4,2))/(FLM+HYM)
      ZCM2=(HYM*CARTNU(3,3)+FLM*CARTNU(4,3))/(FLM+HYM)
C
C*******************************************
C                                          *
C   CENTER-OF-MASS SEPARATION OF MONOMERS  *
C                                          *
C*******************************************
C
      XCM3=XCM2-XCM1
      YCM3=YCM2-YCM1
      ZCM3=ZCM2-ZCM1
C
C******************************************************
C                                                     *
C   DEFINE THE VECTOR FROM CENTER OF MASS OF FIRST    *
C   MONOMER (XCM1,YCM1,ZCM1) TO ATOM H1               *
C                                                     *
C******************************************************
C                        
C      XRM1=XB(1)-XCM1
C      YRM1=XB(2)-YCM1
C      ZRM1=XB(3)-ZCM1
C
      XRM1=CARTNU(1,1)-XCM1
      YRM1=CARTNU(1,2)-YCM1
      ZRM1=CARTNU(1,3)-ZCM1
C
C*************************
C                        *
C   DEFINE COS(THETA1)   *
C                        *
C*************************
C
      THETA1=(XRM1*XCM3+YRM1*YCM3+ZRM1*ZCM3)
      THETA1=THETA1/(SQRT(XRM1**2+YRM1**2+ZRM1**2))
      THETA1=THETA1/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA1.GT.1.0D0)THETA1=1.0D0
      IF(THETA1.LT.-1.0D0)THETA1=-1.0D0
      THETA1=ACOS(THETA1)
C
C******************************************************
C                                                     *
C   DEFINE THE VECTOR FROM CENTER OF MASS OF SECOND   *
C   MONOMER (XCM2,YCM2,ZCM2) TO ATOM H2               *
C                                                     *
C******************************************************
C
C      XRM2=XB(7)-XCM2
C      YRM2=XB(8)-YCM2
C      ZRM2=XB(9)-ZCM2
C
      XRM2=CARTNU(3,1)-XCM2
      YRM2=CARTNU(3,2)-YCM2
      ZRM2=CARTNU(3,3)-ZCM2
C
C*************************
C                        *
C   DEFINE COS(THETA2)   *
C                        *
C*************************
C
      THETA2=(XRM2*(-XCM3)+YRM2*(-YCM3)+ZRM2*(-ZCM3))
      THETA2=THETA2/(SQRT(XRM2**2+YRM2**2+ZRM2**2))
      THETA2=THETA2/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA2.GT.1.0D0)THETA2=1.0D0
      IF(THETA2.LT.-1.0D0)THETA2=-1.0D0
      THETA2=ACOS(THETA2)
      PI=ACOS(-1.0D0)
      THETA2=PI-THETA2
C
C*****************
C                *
C   DEFINE PHI   *
C                *
C*****************
C
      Q1=SQRT(XRM1**2+YRM1**2+ZRM1**2)
      Q2=SQRT(XRM2**2+YRM2**2+ZRM2**2)
      CMM=(XCM3**2+YCM3**2+ZCM3**2)
      CMM=SQRT(CMM)
C      HHD=(XB(1)-XB(7))**2+(XB(2)-XB(8))**2+(XB(3)-XB(9))**2
      HHD=(CARTNU(1,1)-CARTNU(3,1))**2 +
     +    (CARTNU(1,2)-CARTNU(3,2))**2 +
     +    (CARTNU(1,3)-CARTNU(3,3))**2
      HHD=SQRT(HHD)
      Q=CMM-Q1*COS(THETA1)+Q2*COS(THETA2)
      Q3=SQRT(ABS(HHD**2-Q**2))
      Q1=Q1*SIN(THETA1)
      Q2=Q2*SIN(THETA2)
      CPHI=(Q1**2+Q2**2-Q3**2)/(2.*Q1*Q2)
      IF(CPHI.LT.-1.0D0)CPHI=-1.0D0
      IF(CPHI.GT.1.0D0)CPHI=1.0D0
      PHI=ACOS(CPHI)
C
C
C*************************************************
C                                                *
C   ASSIGN INTERNAL COORDINATES TO THE ARRAY R   *
C                                                *
C*************************************************
C
C      RCM=SQRT(XCM3**2+YCM3**2+ZCM3**2)
C      R1=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
C      R2=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
C     WRITE(6,2001)RCM,R1,R2,THETA1,THETA2,PHI
 2001 FORMAT(6F12.8)
C
C************************************************************
C                                                           *
C   Assign internal coordinates to elements of the R array  *
C                                                           *
C************************************************************
C
      R(1)=SQRT(XCM3**2+YCM3**2+ZCM3**2)
      R(2)=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
      R(3)=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
      R(4)=THETA1
      R(5)=THETA2
      R(6)=PHI
C
C      RETURN
C      END
C
      ELSEIF(ICARTR.NE.0) THEN
C
         WRITE(NFLAG(18),1000) ICARTR
 1000    FORMAT(2X,'WRONG ICARTR FOR CARTNU; ICARTR =',I5//)
         STOP
C
      ENDIF
      RETURN
      END
 
 
C
      SUBROUTINE EUNITZERO
C
C************************************************************
C                                                           *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                *
C                                                           *
C   DATE:               00-01-23                            *
C                       01-07-18                            *
C                                                           *
C   THIS PROGRAM CONVERTS OUTPUT ENERGIES FROM              *
C   DEVELOPER UNITS (NAMELY, ATOMIC UNITS) TO USER UNITS.   *
C   CONVERSION CALCULATED IN ANCVRT.F.                      *
C                                                           *
C************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
 
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      PENGYGS = ENGYGS * CNVRTE - ANUZERO
      IF(KSDIAG.NE.0) THEN
         DO I=KSDIAG,KEDIAG
            PENGYES(I) = ENGYES(I) * CNVRTE - ANUZERO
         END DO
      ENDIF
      IF(KSOFFD.NE.0) THEN
         DO J=KSOFFD,KEOFFD
            PENGYIJ(J) = ENGYIJ(J) * CNVRTE
         END DO
      ENDIF
C
      RETURN
      END
 
C
      SUBROUTINE RTOCART
C
C************************************************************
C                                                           *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                *
C                                                           *
C   DATE:               00-01-23                            *
C                       00-08-04                            *
C                       00-08-11                            *
C                       00-08-14                            *
C                       01-07-18                            *
C                                                           *
C   THIS PROGRAM IS DESIGNED TO CONVERT                     *
C   DEGSDR,DEESDR,DEIJDR IN INTERNAL COORDINATES R INTO     *
C   DGSCARTNU,DESCARTNU,DIJCARTNU IN DEVELOPER UNITS AND    *
C   ORDER.                                                  *
C   THIS CODE IS ONLY CALLED IF DERIVATIVES ARE REQUIRED.   *
C   THIS CODE IS THE INVERSE OF CARTTOR.                    *
C                                                           *
C************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
C
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION YGS(N3ATOM),YES(N3ATOM,ISURF),YIJ(N3ATOM,JSURF)
C
C************************************************************
C                                                           *
C   SELECT R -> CARTNU OPTION WITH ICARTR                   *
C                                                           *
C************************************************************
C
      IF(ICARTR.EQ.1) THEN
C
C********************************************************
C                                                       *
C   CARTNU = ARRAY MAP OF VECTOR R                      *
C   Do in order Ground State, Excited State, Coupling   *
C                                                       *
C********************************************************
C
         DO I = 1, NATOMS
            IND=3*I-2
            DGSCARTNU(I,1) = DEGSDR(IND)
            DGSCARTNU(I,2) = DEGSDR(IND+1)
            DGSCARTNU(I,3) = DEGSDR(IND+2)
            IF(KSDIAG.NE.0) THEN
               DO J = KSDIAG,KEDIAG
                  DESCARTNU(I,1,J) = DEESDR(IND,J)
                  DESCARTNU(I,2,J) = DEESDR(IND+1,J)
                  DESCARTNU(I,3,J) = DEESDR(IND+2,J)
               END DO
            ENDIF
            IF(KEOFFD.NE.0) THEN
               DO K = KSOFFD,KEOFFD
                  DIJCARTNU(I,1,K) = DEIJDR(IND,K)
                  DIJCARTNU(I,2,K) = DEIJDR(IND+1,K)
                  DIJCARTNU(I,3,K) = DEIJDR(IND+2,K)
               END DO
            ENDIF
         END DO
C
      ELSEIF(ICARTR.EQ.2) THEN
C
C**************************************************************
C                                                             *
C   DERIVATIVE WRT CARTNU DERIVED FROM DV/DR VIA CHAIN RULE   *
C     WHERE R= INTERPARTICLE DISTANCES IN CANONICAL ORDER:    *
C   R(1)=1-2 DIST.,R(2)=1-3 DIST.,...,R(NATOMS) 2-3 DIST.,... *
C   Do in order Ground State, Excited State, Coupling         *
C                                                             *
C**************************************************************

         DO I = 1, NATOMS
cJC05132003
            DGSCARTNU(I,1) = 0.D0
            DGSCARTNU(I,2) = 0.D0
            DGSCARTNU(I,3) = 0.D0
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                  DESCARTNU(I,1,J1) = 0.D0
                  DESCARTNU(I,2,J1) = 0.D0
                  DESCARTNU(I,3,J1) = 0.D0
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                  DIJCARTNU(I,1,J2) = 0.D0
                  DIJCARTNU(I,2,J2) = 0.D0
                  DIJCARTNU(I,3,J2) = 0.D0
               ENDDO
            ENDIF
            DO J = 1,NATOMS
               IF(J.LT.I) THEN
                  M1 = NATOMS*(J-1) - (J*(J-1))/2 + I-J
               ELSEIF(J.GT.I) THEN
                  M1 = NATOMS*(I-1) - (I*(I-1))/2 + J-I
               ELSE
                  GO TO 20
               ENDIF
               Y = DEGSDR(M1)
               TERMX = (CARTNU(I,1)-CARTNU(J,1))/R(M1)
               TERMY = (CARTNU(I,2)-CARTNU(J,2))/R(M1)
               TERMZ = (CARTNU(I,3)-CARTNU(J,3))/R(M1)
               DGSCARTNU(I,1) = DGSCARTNU(I,1) + TERMX*Y
               DGSCARTNU(I,2) = DGSCARTNU(I,2) + TERMY*Y
               DGSCARTNU(I,3) = DGSCARTNU(I,3) + TERMZ*Y
               IF(KSDIAG.GT.0) THEN
                  Y = DEESDR(M1,J1)
                  DO J1=KSDIAG,KEDIAG
                     DESCARTNU(I,1,J1)=DESCARTNU(I,1,J1) + TERMX*Y
                     DESCARTNU(I,2,J1)=DESCARTNU(I,2,J1) + TERMY*Y
                     DESCARTNU(I,3,J1)=DESCARTNU(I,3,J1) + TERMZ*Y
                  ENDDO
               ELSEIF(KSOFFD.GT.0) THEN
                  DO J2=KSOFFD,KEOFFD
                     Y = DEIJDR(M1,J2)
                     DIJCARTNU(I,1,J2)=DIJCARTNU(I,1,J2) + TERMX*Y
                     DIJCARTNU(I,2,J2)=DIJCARTNU(I,2,J2) + TERMY*Y
                     DIJCARTNU(I,3,J2)=DIJCARTNU(I,3,J2) + TERMZ*Y
                  ENDDO
               ENDIF
20             CONTINUE
            ENDDO
         ENDDO
C
      ELSEIF(ICARTR.EQ.3) THEN
C
C************************************************************
C                                                           *
C   !!!! FOR 3 ATOMS ONLY !!!                               *
C   DERIVATIVE WRT CARTNU DERIVED FROM DV/DR VIA CHAIN RULE *
C     WHERE R(1)=1-2 DIST.,R(2)=2-3 DIST.,R(3)=1-3 DIST.    *
C                                                           *
C   First do prepartory work                                *
C                                                           *
C************************************************************
C
         DO I = 1, NATOMS
            YGS(I) = DEGSDR(I)/R(I)
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  YES(I,J) = DEESDR(I,J)/R(I)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  YIJ(I,K) = DEIJDR(I,K)/R(I)
               ENDDO
            ENDIF
         ENDDO
C
C****************************************************************
C                                                               *
C   Second, do in order Ground State, Excited State, Coupling   *
C                                                               *
C****************************************************************
C
         DO K = 1,3
            TERM12 = CARTNU(1,K)-CARTNU(2,K)
            TERM23 = CARTNU(2,K)-CARTNU(3,K)
            TERM13 = CARTNU(1,K)-CARTNU(3,K)
            DGSCARTNU(1,K) = TERM12*YGS(1) + TERM13*YGS(3)
            DGSCARTNU(2,K) =-TERM12*YGS(1) + TERM23*YGS(2)
            DGSCARTNU(3,K) =-TERM13*YGS(3) - TERM23*YGS(2)
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                 DESCARTNU(1,K,J1) = TERM12*YES(1,J1) + TERM13*YES(3,J1)
                 DESCARTNU(2,K,J1) =-TERM12*YES(1,J1) + TERM23*YES(2,J1)
                 DESCARTNU(3,K,J1) =-TERM13*YES(3,J1) - TERM23*YES(2,J1)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                 DIJCARTNU(1,K,J2) = TERM12*YIJ(1,J2) + TERM13*YIJ(3,J2)
                 DIJCARTNU(2,K,J2) =-TERM12*YIJ(1,J2) + TERM23*YIJ(2,J2)
                 DIJCARTNU(3,K,J2) =-TERM13*YIJ(3,J2) - TERM23*YIJ(2,J2)
               ENDDO
            ENDIF
         ENDDO
C
      ELSEIF(ICARTR.EQ.4) THEN
C
C************************************************************
C                                                           *
C   ICARTR = 4                                              *
C                                                           *
C   Transformation for PES requiring Distances and Angles   *
C                                                           *
C   HF Dimer PESs                                           *
C                                                           *
C************************************************************
C                                                           *
C   CALCULATE THE ATOM-ATOM DISTANCES IN THE HF DIMER       *
C   ARBITRARILY SET ONE HF IN THE X-Z PLANE AND THE OTHER   *
C   ALONG THE X-AXIS, WITH THE FIRST HF'S CENTER OF MASS    *
C   AT THE ORIGIN. FIND THE CARTESIANS OF EACH INDIVIDUAL   * 
C   ATOM AND THEN FIND THE SEPARATIONS OF EACH PAIR OF      *
C   ATOMS AS DESIRED                                        *
C                                                           *     
C   FOR THE CODE BELOW, "A", "B", "C", AND "D"              *
C   REFER TO THE INDIVIDUAL ATOMS OF THE COMPLEX.           *
C                                                           *
C   A IS THE FIRST FLUORINE                                 *
C   B IS THE FIRST HYDRDOGEN                                *
C   C IS THE SECOND FLUORINE                                *
C   D IS THE SECOND HYDROGEN                                *
C                                                           *
C   RCM IS SEPARATION OF THE TWO CENTERS OF MASS            *
C   R1 IS FIRST HF BOND LENGTH                              *
C   R2 IS SECOND HF BOND LENGTH                             *
C   U1 IS COSINE OF THETA1                                  *
C   U2 IS COSINE OF THETA2                                  *
C   U3 IS COSINE OF PHI                                     *
C   SS1 IS SINE OF THETA1                                   *
C   SS2 IS SINE OF THETA2                                   *
C   SS3 IS SINE OF PHI                                      *
C   WH IS MASS OF HYDROGEN IN AMU'S                         *
C   WF IS MASS OF FLUORINE IN AMU'S                         *
C                                                           *
C************************************************************
C
      WH=1.007825D0
      WF=18.99840D0
C
      SUM=WH+WF
      EPS=WF/SUM
      EPSP=WH/SUM
C
      U1=COS(THETA1)
      U2=COS(THETA2)
      U3=COS(PHI)
C
      SS1=SIN(THETA1)
      SS2=SIN(THETA2)
      SS3=SIN(PHI)
C
      YA=0.0D0
      YB=0.0D0
C
      T0=R1*U1
      ZA=-EPSP*T0
      ZB=EPS*T0
C
      T0=R1*SS1
      XA=-EPSP*T0
      XBB=EPS*T0
C
      T0=R2*SS2
      T1=T0*U3
      XC=-EPSP*T1
      XD=EPS*T1
C
      T1=T0*SS3
      YC=-EPSP*T1
      YD=EPS*T1
C
      T0=R2*U2
      ZC=-EPSP*T0+RCM
      ZD=EPS*T0+RCM
C
      RFF=SQRT((XA-XC)**2+YC**2+(ZA-ZC)**2)
C
C      RETURN
C      END
C
      ELSE
C
          WRITE(NFLAG(18),1000) ICARTR
1000      FORMAT(2X,' WRONG ICARTR FOR DERIVATIVE; ICARTR =',I5//)
          STOP
C
      ENDIF
C
      RETURN
      END
 
C
      SUBROUTINE DEDCOU
C
C********************************************************************
C                                                                   *
C   WRITTEN BY: R. J. DUCHOVIC, A. F. WAGNER                        *
C                                                                   *
C   DATE:               00-01-23                                    *
C                       01-07-18                                    * 
C                                                                   *
C   THIS PROGRAM CONVERTS DERIVATIVES                               *
C   DGSCARTNU,DESCARTNU,DIJCARTNU IN DEVELOPER ORDER AND UNITS INTO *
C   DCSCART,DESCART,DIJCART IN USER ODER AND UNITS                  *
C   (DEVELOPER UNITS ARE ATOMIC UNITS)                              *
C   CONVERSION CALCULATED IN ANCVRT.F                               *
C                                                                   *
C********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C***********************************************************************
C                                                                      *
C   Do in order, derivatives of Ground State, Excited State, Couplings *
C                                                                      *
C***********************************************************************
C
      IF (IREORDER.EQ.1) THEN
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(NULBL(I),1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(NULBL(I),2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(NULBL(I),3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(NULBL(I),1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(NULBL(I),2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(NULBL(I),3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(NULBL(I),1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(NULBL(I),2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(NULBL(I),3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ELSE
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(I,1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(I,2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(I,3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(I,1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(I,2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(I,3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(I,1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(I,2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(I,3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ENDIF
C
      RETURN
      END
 
