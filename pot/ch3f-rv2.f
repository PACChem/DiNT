      subroutine pot(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)
c      program ch3f

      implicit none

      integer mnsurf,nsurf,mnat,nat,n3,i,ii,j,k,n,irow,icol,infoflag 
c      parameter(mnsurf=12)
c      parameter(mnat=5)
      parameter(n=12)
      double precision a(n,n),v(n,n),d(n)
      double precision upperPacked(n*(n+1)/2),work(3*n)
      double precision autoang
      parameter(autoang=0.52917706d0)
      character*2 symb(mnat)
      double precision sing3(3,3),gsing3(3,mnat,3,3),trip3(3,3),
     & gtrip3(3,mnat,3,3),sing3b(3,3)
      double precision sing3bis(3,3)
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf)
      double precision x(mnat),y(mnat),z(mnat),xx(20),yy(20),
     & zz(20)
      double precision autokcal,singsig,tripsig1,ajsigma,
     & ajpi,c1,c2,trip3a,trip3b,singpix,singpiy,trippix1,trippiy1
      double precision ajsigma0
      double precision c1a,c1b,c2a,c2b
      double precision singpi,trippi
      double precision singsig2,tripsig2,singpix2,singpiy2,
     &trippix2,trippiy2
      double precision tripsig,trippix,trippiy
      double precision singsig1,singpix1,singpiy1,singpi1,
     &singpi2
      double precision pi,a1,a2,a3,a4,a5,a6,a7,a8,th,phi
      complex*16 SOC(12,12),omega(12,12),SOCadj(12,12),interm(12,12)
      double precision lambda(12,12),omega2(12,12),delta,aa,int
      double precision z1,z2,z3,z4,z5
      double precision Fp(17),coefF(60),coefCl(60),coefBr(60)
      double precision xcmh1,ycmh1,zcmh1,xcmh2,ycmh2,zcmh2,xcmh3,ycmh3,
     & zcmh3,intx,inty,intz,nx,ny,nz,px,py,pz,rcmh1,rcmh2,rcmh3,phi1,
     & phi2,phi3,pmod,xh1h2,yh1h2,zh1h2,xh2h3,yh2h3,zh2h3,th1,xfcm,
     & yfcm,zfcm,rfcm,xcm,ycm,zcm,rh1h2,rh2h3
      double precision angh1cmh2,angh1cmh3,angh2cmh3
      double precision gajsigma(3,mnat),gajpi(3,mnat),gtrip3a(3,mnat),
     &gtrip3b(3,mnat),inter(3,mnat)
      double precision ech3,gech3(3,mnat)
      integer itimes,l,m
      double precision glambda(3,mnat,12,12),ginterm(3,mnat,12,12) 
      double precision plgndr

      common /pot_coe/ coefF,coefCl,coefBr
      common /SOC_matrix/ SOC

      pi = 3.14159274
      autokcal = 627.509470

      call prepot

      nat=4
      nsurf=3

c geometry for CH3 + F in A
c      x(1)=    0.000000000
c      y(1)=    0.000000000
c      z(1)=    0.000000000
c      x(2)=    1.0895667439
c      y(2)=    0.0000000000
c      z(2)=    0.0000000000
c      x(3)=   -0.5447833549
c      y(3)=   -0.9435924877
c      z(3)=    0.0000000000
c      x(4)=   -0.5447833910
c      y(4)=    0.9435924877
c      z(4)=    0.0000000000
c      x(5) =  -0.520945
c      y(5) =   0.513030
c      z(5) =  -2.909539

c      do i=1,nat+1
c        x(i)=x(i)/autoang
c        y(i)=y(i)/autoang
c        z(i)=z(i)/autoang
c      enddo

c compute and return the three singlet energies and gradients
c     sing3 contains the 3 singlet energies
c     gsing3 contians the gradients of the 3 singlet energies
c     the off diagonal elements of sing3 and gsing3 are not used

c First compute energy of CH3 at current CH3F geometry

      symb(1)="c"
      symb(2)="h"
      symb(3)="h"
      symb(4)="h"

      call system("cp qcch3.m06 qc.m06")

      n3=1
      call dd1(symb,x,y,z,sing3,gsing3,nat,mnat,n3,n3)
      
      ech3 = sing3(1,1) + 39.69314153

      do i=1,3
       do j=1,nat 
        gech3(i,j) = gsing3(i,j,1,1) 
       enddo
      enddo

C     The Cartesian components of the gradient for F are zero
 
      do i=1,3
      gech3(i,5) = 0.d0
      enddo

c Now do CH3F
c these will be replaced with passed variables once this program is transformed into a subroutine

      nat=5
      nsurf=3
      
      symb(1)="c"
      symb(2)="h"
      symb(3)="h"
      symb(4)="h"
      symb(5)="f"

      call system("cp qcch4f.m06 qc.m06")

      n3=3
      call dd(symb,x,y,z,sing3,gsing3,nat,mnat,n3,n3)

c print x,y,z components of the gradient for each atom for surface 1
c     do i=1,mnat
c      print *,"gE1(",i,") = ",gsing3(1,i,1,1),gsing3(2,i,1,1),
c    &                          gsing3(3,i,1,1)," hartree/bohr"
c     enddo
c     print *

c print x,y,z components of the gradient for each atom for surface 2
c     do i=1,mnat
c       print *,"gE2(",i,") = ",gsing3(1,i,2,2),gsing3(2,i,2,2),
c    &                          gsing3(3,i,2,2)," hartree/bohr"
c     enddo
c     print *

c print x,y,z components of the gradient for each atom for surface 3
c     do i=1,mnat
c       print *,"gE3(",i,") = ",gsing3(1,i,3,3),gsing3(2,i,3,3),
c    &                          gsing3(3,i,3,3)," hartree/bohr"
c     enddo

c Do your magic here....

c Eventually, package everything into the 12x12 diabatic matrix PEMD and its gradient GPEMD
C
C       Calculate general expression of theta, phi for arbitrary geometry
C       First determine normal to plane of three H

        xh1h2 = x(3) - x(2) 
        yh1h2 = y(3) - y(2)
        zh1h2 = z(3) - z(2)
        xh2h3 = x(4) - x(3)
        yh2h3 = y(4) - y(3)
        zh2h3 = z(4) - z(3)
        rh1h2 = dsqrt(xh1h2**2 + yh1h2**2 + zh1h2**2)
        rh2h3 = dsqrt(xh2h3**2 + yh2h3**2 + zh2h3**2)
        th1 = acos((xh1h2*xh2h3 + yh1h2*yh2h3 + zh1h2*zh2h3)/
     &  (rh1h2*rh2h3))
        nx = (yh1h2*zh2h3 - zh1h2*yh2h3)/(rh1h2*rh2h3*sin(th1))
        ny = -(xh1h2*zh2h3 - zh1h2*xh2h3)/(rh1h2*rh2h3*sin(th1))
        nz = (xh1h2*yh2h3 - yh1h2*xh2h3)/(rh1h2*rh2h3*sin(th1))

C       Then determine center of mass of three H and the vector from F (Cl,Br) to CM
C       Theta is the angle between the vector normal to the plane and the F-CM vector

        xcm = (x(2) + x(3) + x(4))/3.d0
        ycm = (y(2) + y(3) + y(4))/3.d0
        zcm = (z(2) + z(3) + z(4))/3.d0

        xfcm = x(5) - xcm
        yfcm = y(5) - ycm
        zfcm = z(5) - zcm
        rfcm = dsqrt(xfcm**2 + yfcm**2 + zfcm**2)

        th = 360.d0/(2.d0*pi)*acos((nx*xfcm + ny*yfcm + nz*zfcm)/rfcm)

C       Determine Phi1, Phi2, Phi3. They are the angles between the projection of the F-CM vector
C       on the plane of the three Hs and the lines CM-H1, CM-H2, CM-H3.
C       Final Phi is the average of Phi1, Phi2, Phi3 

        xcmh1 = x(2) - xcm 
        ycmh1 = y(2) - ycm
        zcmh1 = z(2) - zcm
  
        xcmh2 = x(3) - xcm
        ycmh2 = y(3) - ycm
        zcmh2 = z(3) - zcm
 
        xcmh3 = x(4) - xcm
        ycmh3 = y(4) - ycm
        zcmh3 = z(4) - zcm

        rcmh1 = dsqrt(xcmh1**2 + ycmh1**2 + zcmh1**2)
        rcmh2 = dsqrt(xcmh2**2 + ycmh2**2 + zcmh2**2)
        rcmh3 = dsqrt(xcmh3**2 + ycmh3**2 + zcmh3**2) 

C       Projection of F-CM vector on the plane is n^(rfcm^n)

        intx = yfcm*nz - zfcm*ny
        inty = -(xfcm*nz - zfcm*nx)
        intz = xfcm*ny - yfcm*nx

        px = ny*intz - nz*inty
        py = -(nx*intz - nz*intx)
        pz = nx*inty - ny*intx

        pmod = dsqrt(px**2 + py**2 + pz**2)

      phi1 = 360.d0/(2.d0*pi)*acos((px*xcmh1+py*ycmh1+pz*zcmh1)/
     &  (pmod*rcmh1))
      phi2 = 360.d0/(2.d0*pi)*acos((px*xcmh2+py*ycmh2+pz*zcmh2)/
     &  (pmod*rcmh2))
      phi3 = 360.d0/(2.d0*pi)*acos((px*xcmh3+py*ycmh3+pz*zcmh3)/
     &  (pmod*rcmh3))

C     With this formula the three phi are not measured always in the same
C     rotation sense (clockwise or anticlockwise). To correct for this
C     we determine which two phi add up to one of the angles on the plane
C     of the three Hs. For example the sum of phi2 and phi3 can be equal to H2-CM-H3.
C     Then we know that the projection of F-CM on the plane of the three Hs is
C     between the projection of CM-H2 and CM-H3 on the plane and that
C     phi2 and phi3 are measured in opposite rotation senses. To determine the
C     sense of phi1 we see if phi1-phi3 is equal to the angle H1-CM-H3 or if
C     phi1-phi2 is equal to the angle H1-CM-H2. In the first case phi1 and phi3
C     have the same sense, in the second phi1 and phi2 have the same sense.
C     In the first case phi2 is transformed into 360 - phi2, and in the second
C     case phi3 is transformed into 360 - phi3.
C     The same procedure is followed if the sum of phi1 and phi2 is equal to H1-CM-H2
C     and if the sum of phi1 and phi3 is equal to H1-CM-H3

      angh1cmh2 = 360.d0/(2.d0*pi)*acos((xcmh1*xcmh2 + ycmh1*ycmh2 
     & + zcmh1*zcmh2)/(rcmh1*rcmh2))

      angh2cmh3 = 360.d0/(2.d0*pi)*acos((xcmh2*xcmh3 + ycmh2*ycmh3 
     & + zcmh2*zcmh3)/(rcmh2*rcmh3))

      angh1cmh3 = 360.d0/(2.d0*pi)*acos((xcmh1*xcmh3 + ycmh1*ycmh3 
     & + zcmh1*zcmh3)/(rcmh1*rcmh3))

      if(abs((phi1+phi2)-angh1cmh2).lt.0.01) then
       if(abs((phi3-phi1)-angh1cmh3).lt.0.01) then
         phi2 = 360.d0 - phi2
       elseif(abs((phi3-phi2)-angh2cmh3).lt.0.01) then
         phi1 = 360.d0 - phi1
       endif
      endif

      if(abs((phi2+phi3)-angh2cmh3).lt.0.01) then
       if(abs((phi1-phi2)-angh1cmh2).lt.0.01) then
         phi3 = 360.d0 - phi3
       elseif(abs((phi1-phi3)-angh1cmh3).lt.0.01) then
         phi2 = 360.d0 - phi2
       endif
      endif

      if(abs((phi1+phi3)-angh1cmh3).lt.0.01) then
       if(abs((phi2-phi1)-angh1cmh2).lt.0.01) then
         phi3 = 360.d0 - phi3
       elseif(abs((phi2-phi3)-angh2cmh3).lt.0.01) then
         phi1 = 360.d0 - phi1
       endif
      endif

        phi = (phi1 + phi2 + phi3)/3.d0

C     Fitted singlet energies at C-F = 2.3 Angs

        z1 = cos(2.0*pi*th/360.0)
        z2 = cos(6.0*pi*phi/360.0)
        z3 = sin(6.0*pi*phi/360.0)
        z4 = cos(12.0*pi*phi/360.0)
        z5 = sin(12.0*pi*phi/360.0)
        Fp(1) = 1.0
        Fp(2) = plgndr(2,0,z1)
        Fp(3) = plgndr(3,3,z1) * z2
        Fp(4) = plgndr(3,3,z1) * z3
        Fp(5) = plgndr(4,0,z1)
        Fp(6) = plgndr(5,3,z1) * z2
        Fp(7) = plgndr(5,3,z1) * z3
        Fp(8) = plgndr(6,0,z1)
        Fp(9) = plgndr(6,6,z1) * z4
        Fp(10) = plgndr(6,6,z1) * z5

        singsig = 0.d0
        singpix = 0.d0
        singpiy = 0.d0
        tripsig = 0.d0
        trippix = 0.d0
        trippiy = 0.d0

        if(symb(5).eq."f") then

        do i=1,10
         singsig = singsig + Fp(i)*coefF(i)
        enddo

        do i=1,10
         singpix = singpix + Fp(i)*coefF(i+10)
        enddo

        do i=1,10
         singpiy = singpiy + Fp(i)*coefF(i+20)
        enddo

        do i=1,10
         tripsig = tripsig + Fp(i)*coefF(i+30)
        enddo

        do i=1,10
         trippix = trippix + Fp(i)*coefF(i+40)
        enddo

        do i=1,10
         trippiy = trippiy + Fp(i)*coefF(i+50)
        enddo

        elseif(symb(5).eq."cl") then

        do i=1,10
         singsig = singsig + Fp(i)*coefCl(i)
        enddo

        do i=1,10
         singpix = singpix + Fp(i)*coefCl(i+10)                         
        enddo

        do i=1,10
         singpiy = singpiy + Fp(i)*coefCl(i+20)
        enddo

        do i=1,10
         tripsig = tripsig + Fp(i)*coefCl(i+30)
        enddo

        do i=1,10
         trippix = trippix + Fp(i)*coefCl(i+40)
        enddo

        do i=1,10
         trippiy = trippiy + Fp(i)*coefCl(i+50)
        enddo

        elseif(symb(5).eq."br") then

        do i=1,10
         singsig = singsig + Fp(i)*coefBr(i)
        enddo

        do i=1,10
         singpix = singpix + Fp(i)*coefBr(i+10)                         
        enddo

        do i=1,10
         singpiy = singpiy + Fp(i)*coefBr(i+20)
        enddo

        do i=1,10
         tripsig = tripsig + Fp(i)*coefBr(i+30)
        enddo

        do i=1,10
         trippix = trippix + Fp(i)*coefBr(i+40)
        enddo

        do i=1,10
         trippiy = trippiy + Fp(i)*coefBr(i+50)
        enddo

        endif

C     Apply expression for triplets as functions of singlets

      singpi = (singpix + singpiy)/2.d0
      trippi = (trippix + trippiy)/2.d0

      c1 = (singsig + 2.d0*singpi - 7.d0*tripsig + 6.d0*trippi)/
     &  (5.d0*(singsig - tripsig))
      c2 = (singsig + 2.d0*singpi + 3.d0*tripsig - 4.d0*trippi)/
     &  (5.d0*(singpi - trippi))

        sing3bis(1,1) = sing3(1,1)
        sing3bis(2,2) = sing3(2,2)

        ajsigma = (3.d0*c2*(sing3bis(1,1)-ech3) - (4.d0*c2 - 2.d0)*
     &  (sing3bis(2,2) - ech3))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )
        ajpi = ( (1.d0+c1)*(sing3bis(2,2) - ech3) - (2.d0*c1 - 1.d0)*
     &  (sing3bis(1,1) - ech3))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )

        trip3(2,2)= sing3(2,2) - 2.d0*ajpi
        trip3a = sing3(1,1) - 2.d0*ajsigma

        sing3b(3,3) = sing3(2,2) + (sing3(3,3)-sing3(2,2))*
     &  (trippiy - trippix)/(singpiy - singpix)

        sing3bis(3,3) = sing3b(3,3)

        ajsigma = (3.d0*c2*(sing3bis(1,1) - ech3)-(4.d0*c2 - 2.d0)*
     &  (sing3bis(3,3) - ech3))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )
        ajpi = ( (1.d0+c1)*(sing3bis(3,3)-ech3) - (2.d0*c1 - 1.d0)*
     &  (sing3bis(1,1) - ech3))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )

        trip3(3,3)= sing3b(3,3) - 2.d0*ajpi
        trip3b = sing3(1,1) - 2.d0*ajsigma

        trip3(1,1) = (trip3a + trip3b)/2.d0

C     Calculate gradients of triplet states

      do i=1,3
       do j=1,mnat

        gajsigma(i,j) = (3.d0*c2*(gsing3(i,j,1,1) - gech3(i,j)) 
     &  - (4.d0*c2 - 2.d0)*
     &  (gsing3(i,j,2,2) - gech3(i,j)))/( 3.d0*c2*(c1+1.d0)-
     &  (2.d0*c1 - 1.d0)* (4.d0*c2 - 2.d0))

        gajpi(i,j) = ( (1.d0+c1)*(gsing3(i,j,2,2) - gech3(i,j)) - 
     &  (2.d0*c1 - 1.d0)*(gsing3(i,j,1,1) - gech3(i,j)))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     &  (4.d0*c2 - 2.d0) )

        gtrip3(i,j,2,2)= gsing3(i,j,2,2) - 2.d0*gajpi(i,j)
        gtrip3a(i,j) = gsing3(i,j,1,1) - 2.d0*gajsigma(i,j)

        gajsigma(i,j) = (3.d0*c2*(gsing3(i,j,1,1) - gech3(i,j)) - 
     &  (4.d0*c2 - 2.d0)*(gsing3(i,j,3,3) - gech3(i,j)))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     &  (4.d0*c2 - 2.d0) )

        gajpi(i,j) = ( (1.d0+c1)*(gsing3(i,j,3,3) - gech3(i,j)) - 
     &  (2.d0*c1 - 1.d0)*(gsing3(i,j,1,1) - gech3(i,j)))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     &  (4.d0*c2 - 2.d0) )

        gtrip3(i,j,3,3)= gsing3(i,j,3,3) - 2.d0*gajpi(i,j)
        gtrip3b(i,j) = gsing3(i,j,1,1) - 2.d0*gajsigma(i,j)

        gtrip3(i,j,1,1) = (gtrip3a(i,j) + gtrip3b(i,j))/2.d0

       enddo
      enddo

C     Transform singlets and triplets to diabatic spin-coupled representation

       do i=1,12
        do j=1,12
         SOCadj(i,j) = conjg(SOC(j,i))
        enddo
       enddo

       do i=1,12
        do j=1,12
        lambda(i,j) = (0.d0,0.d0)
        enddo
       enddo

C      Diabatic energies

        lambda(1,1) = sing3(1,1)
        lambda(2,2) = trip3(2,2)
        lambda(3,3) = trip3(3,3)
        lambda(4,4) = sing3(2,2)
        lambda(5,5) = trip3(1,1)
        lambda(6,6) = trip3(3,3)
        lambda(7,7) = sing3(3,3)
        lambda(8,8) = trip3(1,1)
        lambda(9,9) = trip3(2,2)
        lambda(10,10) = trip3(1,1)
        lambda(11,11) = trip3(2,2)
        lambda(12,12) = trip3(3,3)

C      Diabatic gradients

      do i=1,3
       do j=1,mnat

         glambda(i,j,1,1) = gsing3(i,j,1,1) 
         glambda(i,j,2,2) = gtrip3(i,j,2,2)
         glambda(i,j,3,3) = gtrip3(i,j,3,3)
         glambda(i,j,4,4) = gsing3(i,j,2,2)
         glambda(i,j,5,5) = gtrip3(i,j,1,1)
         glambda(i,j,6,6) = gtrip3(i,j,3,3)
         glambda(i,j,7,7) = gsing3(i,j,3,3)
         glambda(i,j,8,8) = gtrip3(i,j,1,1)
         glambda(i,j,9,9) = gtrip3(i,j,2,2)
         glambda(i,j,10,10) = gtrip3(i,j,1,1)
         glambda(i,j,11,11) = gtrip3(i,j,2,2)
         glambda(i,j,12,12) = gtrip3(i,j,3,3)

        enddo
       enddo

       do i=1,12
        do j=1,12
        interm(i,j) = (0.d0,0.d0)
        omega(i,j) = (0.d0,0.d0)
        do k=1,3
         do l=1,mnat
         ginterm(k,l,i,j) = (0.d0,0.d0)
         gpemd(k,l,i,j) = (0.d0,0.d0) 
         enddo
        enddo
       enddo
       enddo

C     Apply SOC transformation matrix to lambda matrix and add SOC

      do j=1,12
       do i=1,12
        do k=1,12
          interm(i,j) = interm(i,j) + SOCadj(i,k) * lambda(k,j)
        enddo
       enddo
      enddo

      do j=1,12
       do i=1,12
        do k=1,12
          omega(i,j) = omega(i,j) + interm(i,k) * SOC(k,j)
        enddo
       enddo
      enddo

C     Apply the same transformation to the spin-free gradients

      do l=1,3
      do m=1,mnat
      do j=1,12
       do i=1,12
        do k=1,12
          ginterm(l,m,i,j) = ginterm(l,m,i,j) + SOCadj(i,k) * 
     &    glambda(l,m,k,j)
        enddo
       enddo
      enddo
      enddo
      enddo

      do l=1,3
      do m=1,mnat
      do j=1,12
       do i=1,12
        do k=1,12
          gpemd(l,m,i,j) = gpemd(l,m,i,j) + ginterm(l,m,i,k) * SOC(k,j)
        enddo
       enddo
      enddo
      enddo
      enddo

c print x,y,z components of the gradient for each atom for surface 1
c     do i=1,mnat
c       print *,"gE1(",i,") = ",gpemd(1,i,1,1),gpemd(2,i,1,1),
c    &                          gpemd(3,i,1,1)," hartree/bohr"
c     enddo
c     print *

c print x,y,z components of the gradient for each atom for surface 2
c     do i=1,mnat
c       print *,"gE2(",i,") = ",gpemd(1,i,2,2),gpemd(2,i,1,1),
c    &                          gpemd(3,i,2,2)," hartree/bohr"
c     enddo
c     print *

c print x,y,z components of the gradient for each atom for surface 3
c     do i=1,mnat
c       print *,"gE3(",i,") = ",gpemd(1,i,3,3),gpemd(2,i,3,3),
c    &                          gpemd(3,i,3,3)," hartree/bohr"
c     enddo

C     Add experimental SOC constant to diagonal elements of diabatic spin-coupled matrix
C     Fluorine: SOC = 404.1 cm-1 = 1.16 kcal/mol
C     Chlorine: SOC = 882.4 cm-1 = 2.52 kcal/mol
C     Bromine:  SOC = 3685.2 cm-1 = 10.54 kcal/mol
C     Source is NIST Chemistry Webbook
C     delta = 1/3 of SOC constant

      if(symb(5).eq."f") then
      delta = 0.385/autokcal

      elseif(symb(5).eq."cl") then
      delta = 0.841/autokcal

      elseif(symb(5).eq."br") then
      delta = 3.512/autokcal

      endif

C     Add SOC constant to diagonal elements of diabatic spin-coupled matrix

      omega(1,1) = omega(1,1) - delta
      omega(2,2) = omega(2,2) - delta
      omega(3,3) = omega(3,3) + 2.d0*delta
      omega(4,4) = omega(4,4) - delta
      omega(5,5) = omega(5,5) + 2.d0*delta
      omega(6,6) = omega(6,6) - delta
      omega(7,7) = omega(7,7) - delta
      omega(8,8) = omega(8,8) + 2.d0*delta
      omega(9,9) = omega(9,9) - delta
      omega(10,10) = omega(10,10) + 2.d0*delta
      omega(11,11) = omega(11,11) - delta
      omega(12,12) = omega(12,12) - delta

       do i=1,12
        do j=1,12
         pemd(i,j) = dble(omega(i,j))
        enddo
       enddo

      nsurf=12
      nat=5

      end

       FUNCTION plgndr(l,m,x)
       integer l,m
       real*8 plgndr,x
       integer i,ll
       real*8 fact,pll,pmm,pmmp1,somx2
       if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0)
     &  pause 'bad arguments in plgndr'
       pmm=1.0
       if(m.gt.0) then
        somx2=sqrt((1.0-x)*(1.0+x))
        fact=1.0
        do i=1,m
         pmm = -pmm * fact* somx2
         fact = fact + 2.0
        enddo
       endif
       if(l.eq.m) then
        plgndr = pmm
       else
        pmmp1 = x*(2*m + 1) * pmm
        if(l.eq.m+1) then
         plgndr = pmmp1
        else
          do ll = m+2,l 
           pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm = pmmp1 
           pmmp1 = pll 
          enddo
           plgndr = pll
        endif 
       endif 
       return
       end    
c **********************************************************************
c **********************************************************************
c     DD: A common interface for direct dynamics calls to the Gaussian & 
c     Molpro packages. Jasper, Dec 2007
c
c     Some of this has been borrowed from the POLYRATE HOOKS programs
c **********************************************************************
c **********************************************************************

      subroutine dd(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)
c     this subroutine has been specialized for CH3 + F
c     three singlet energies and gradients are computed and returned

      implicit none
      double precision autoang,autocmi
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)

c in/out
      integer nat,mnat,nsurf,mnsurf
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     & gpemd(3,mnat,mnsurf,mnsurf),pemd(mnsurf,mnsurf)

c local
      integer mnc,ifail
      parameter(mnc=10) ! max number of QC calls per geometry
      integer ic,nc,ns(mnc),np(mnc),ntmp,ictot,i,j,k,l,ii,jj
      double precision gptmp(3,mnat,mnsurf,mnsurf),ptmp(mnsurf,mnsurf)
      integer it1,it2
      character*8 tname(mnc),tnameB(mnc),tnam

c SO coupling
      integer icc,iff
      double precision x1,x2,somax,so,del,rcent,rcom,arg
      double precision dx,dy,dz,darg,ddel,dso

c -------------------
c SET THESE
c     NC = number of separate QC calls per geom
c     NS(i) = number of states for call #i
c     NP(i) = QC package to be used for call #i
c           = 1 for G03
c           = 2 for Molpro 2006
      nc    = 1
      ns(1) = 3
      np(1) = 2
      tname(1) = "qc.m06"
c      tnameB(1) = "qcB.m06"
c -------------------

      if (nc.gt.mnc) then
          print *,"nc (",nc,") is bigger than mnc (",mnc,")" 
          stop
      endif

c zero
      do i=1,mnsurf
      do j=1,mnsurf
        pemd(i,j)=0.d0
        do k=1,3
        do l=1,mnat
          gpemd(k,l,i,j)=0.d0
        enddo
        enddo
      enddo
      enddo

c make calls
      ifail=0
      ictot=0
      do ic = 1,nc
        ntmp = ns(ic)
        tnam = tname(ic)
        if (np(ic).eq.1) then
c         call G03
          call dd_g03(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
        elseif (np(ic).eq.2) then
c         call Molpro 2006
          call system_clock(it1)
          call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
          call system_clock(it2)
          write(83,*)dble(it2-it1)/1.d2,ifail
c this part has been turned off
c          IF (.FALSE.) THEN
c this part has been turned on
          IF (.TRUE.) THEN
c         handle failures
          if (ifail.eq.0) then
c           everything OK
          elseif (ifail.eq.10000) then
c           retry with numerical gradient in backup template
            ifail=0
            do i=1,mnsurf
            do j=1,mnsurf
              pemd(i,j)=0.d0
              do k=1,3
              do l=1,mnat
                gpemd(k,l,i,j)=0.d0
              enddo
              enddo
            enddo
            enddo
            tnam=tnameB(ic)
            ifail=0
            call system_clock(it1)
            call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
            call system_clock(it2)
            write(83,*)dble(it2-it1)/1.d2,ifail
            if (ifail.ne.0) then
c             failed again
              write(6,*)"QC backup file failure"
              stop
            endif
          else
c           can't fix
            write(6,*)"QC failure"
            stop
          endif
          ENDIF
        else
          print *,"np(",ic,") = ",np(ic),", which is not allowed"
          stop
        endif
        do i=1,ntmp
        do j=1,ntmp
          ii=i+ictot
          jj=j+ictot
          pemd(ii,jj)=ptmp(i,j)
          do k=1,3
          do l=1,nat
            gpemd(k,l,ii,jj)=gptmp(k,l,i,j)
          enddo
          enddo
        enddo
        enddo
        ictot=ictot+ntmp
      enddo

      return

      end

      subroutine dd1(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)
c     this subroutine has been specialized for CH3 + F
c     three singlet energies and gradients are computed and returned

      implicit none
      double precision autoang,autocmi
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)

c in/out
      integer nat,mnat,nsurf,mnsurf
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     & gpemd(3,mnat,mnsurf,mnsurf),pemd(mnsurf,mnsurf)

c local
      integer mnc,ifail
      parameter(mnc=10) ! max number of QC calls per geometry
      integer ic,nc,ns(mnc),np(mnc),ntmp,ictot,i,j,k,l,ii,jj
      double precision gptmp(3,mnat,mnsurf,mnsurf),ptmp(mnsurf,mnsurf)
      integer it1,it2
      character*8 tname(mnc),tnameB(mnc),tnam

c SO coupling
      integer icc,iff
      double precision x1,x2,somax,so,del,rcent,rcom,arg
      double precision dx,dy,dz,darg,ddel,dso

c -------------------
c SET THESE
c     NC = number of separate QC calls per geom
c     NS(i) = number of states for call #i
c     NP(i) = QC package to be used for call #i
c           = 1 for G03
c           = 2 for Molpro 2006
      nc    = 1
      ns(1) = 1
      np(1) = 2
      tname(1) = "qc.m06"
c      tnameB(1) = "qcB.m06"
c -------------------

      if (nc.gt.mnc) then
          print *,"nc (",nc,") is bigger than mnc (",mnc,")" 
          stop
      endif

c zero
      do i=1,mnsurf
      do j=1,mnsurf
        pemd(i,j)=0.d0
        do k=1,3
        do l=1,mnat
          gpemd(k,l,i,j)=0.d0
        enddo
        enddo
      enddo
      enddo

c make calls
      ifail=0
      ictot=0
      do ic = 1,nc
        ntmp = ns(ic)
        tnam = tname(ic)
        if (np(ic).eq.1) then
c         call G03
          call dd_g03(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
        elseif (np(ic).eq.2) then
c         call Molpro 2006
          call system_clock(it1)
          call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
          call system_clock(it2)
          write(83,*)dble(it2-it1)/1.d2,ifail
c this part has been turned off
c          IF (.FALSE.) THEN
c this part has been turned on
          IF (.TRUE.) THEN
c         handle failures
          if (ifail.eq.0) then
c           everything OK
          elseif (ifail.eq.1000) then
c           retry with numerical gradient in backup template
            ifail=0
            do i=1,mnsurf
            do j=1,mnsurf
              pemd(i,j)=0.d0
              do k=1,3
              do l=1,mnat
                gpemd(k,l,i,j)=0.d0
              enddo
              enddo
            enddo
            enddo
            tnam=tnameB(ic)
            ifail=0
            call system_clock(it1)
            call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
            call system_clock(it2)
            write(83,*)dble(it2-it1)/1.d2,ifail
            if (ifail.ne.0) then
c             failed again
              write(6,*)"QC backup file failure"
              stop
            endif
          else
c           can't fix
            write(6,*)"QC failure"
            stop
          endif
          ENDIF
        else
          print *,"np(",ic,") = ",np(ic),", which is not allowed"
          stop
        endif
        do i=1,ntmp
        do j=1,ntmp
          ii=i+ictot
          jj=j+ictot
          pemd(ii,jj)=ptmp(i,j)
          do k=1,3
          do l=1,nat
            gpemd(k,l,ii,jj)=gptmp(k,l,i,j)
          enddo
          enddo
        enddo
        enddo
        ictot=ictot+ntmp
      enddo

      return

      end

c **********************************************************************
c **********************************************************************
      subroutine prepot
     
      double precision coefF(60),coefCl(60),coefBr(60)
      complex*16 SOC(12,12)
      common /pot_coe/ coefF,coefCl,coefBr
      common /SOC_matrix/ SOC

C       Coefficients for CH3F

        coefF( 1 ) =  -5.66085005574342848
        coefF( 2 ) =  -30.4526634110748695
        coefF( 3 ) =  -0.283339724865703912
        coefF( 4 ) =  0.909144408983299235E-01
        coefF( 5 ) =  11.2519578553738455
        coefF( 6 ) =  0.155429650424575936E-01
        coefF( 7 ) =  -0.167968903563195660E-01
        coefF( 8 ) =  -4.01197981608117171
        coefF( 9 ) =  0.496438759321592900E-04
        coefF( 10 ) =  0.112954592552487292E-03

        coefF( 11 ) =  18.3255793525378827
        coefF( 12 ) =  -9.94297659807739187
        coefF( 13 ) =  -0.890684189636730683
        coefF( 14 ) =  0.155384408888473971E-03
        coefF( 15 ) =  1.47339035807007068
        coefF( 16 ) =  0.207476442847030207E-01
        coefF( 17 ) =  0.626148924425637258E-03
        coefF( 18 ) =  0.647909103693702804
        coefF( 19 ) =  0.975957240767974656E-04
        coefF( 20 ) =  -0.417566258918837223E-04

        coefF( 21 ) =  21.4359404051685907
        coefF( 22 ) =  -15.6876073603005342
        coefF( 23 ) =  -0.988269471915655950
        coefF( 24 ) =  0.675225779691994432E-01
        coefF( 25 ) =  6.68269656748551810
        coefF( 26 ) =  0.213656675130531395E-01
        coefF( 27 ) =  -0.108799232729098593E-01
        coefF( 28 ) =  -2.01975402256017045
        coefF( 29 ) =  0.839832495730722089E-04
        coefF( 30 ) =  0.322413949556970045E-04

        coefF( 31 ) =  6.91306946781652965
        coefF( 32 ) =  1.77904847263440957
        coefF( 33 ) =  -0.166358992116615145
        coefF( 34 ) =  -0.297634972110518391E-01
        coefF( 35 ) =  -0.589078849492522871
        coefF( 36 ) =  0.686717510694524486E-02
        coefF( 37 ) =  0.517590194980125772E-02
        coefF( 38 ) =  0.329942580669068508
        coefF( 39 ) =  0.247301418375610871E-04
        coefF( 40 ) =  -0.535849718325920571E-05

        coefF( 41 ) =  20.4782061386153913
        coefF( 42 ) =  -16.8775446738371890
        coefF( 43 ) =  -0.993215701797845352
        coefF( 44 ) =  0.761466860447098243E-01
        coefF( 45 ) =  6.79924213377671638
        coefF( 46 ) =  0.199529327570264411E-01
        coefF( 47 ) =  -0.125532006490565937E-01
        coefF( 48 ) =  -2.01514826386988899
        coefF( 49 ) =  0.957755326252302132E-04
        coefF( 50 ) =  0.293335322338308182E-04

        coefF( 51 ) =  21.0203209011237107
        coefF( 52 ) =  -16.2402111428953226
        coefF( 53 ) =  -1.00009872772126296
        coefF( 54 ) =  0.593420895673896151E-01
        coefF( 55 ) =  6.16560991245408374
        coefF( 56 ) =  0.210470850892566404E-01
        coefF( 57 ) =  -0.995189146073606305E-02
        coefF( 58 ) =  -1.41649784779081700
        coefF( 59 ) =  0.105818199437863436E-03
        coefF( 60 ) =  0.339186475065961953E-04

C       Coefficients for CH3Cl

        coefCl( 1 ) =  -6.19487133690628511
        coefCl( 2 ) =  -29.6140919555237474
        coefCl( 3 ) =  -0.295124601673081410
        coefCl( 4 ) =  0.615849422909814023E-01
        coefCl( 5 ) =  9.36718118507110198
        coefCl( 6 ) =  0.106269289035111486E-01
        coefCl( 7 ) =  -0.117836912964743076E-01
        coefCl( 8 ) =  -3.15785872919283284
        coefCl( 9 ) =  0.462877168084204602E-04
        coefCl( 10 ) =  0.776363519749927443E-04

        coefCl( 11 ) =  16.3830870679290364
        coefCl( 12 ) =  -7.66677085235455724
        coefCl( 13 ) =  -0.633226485826246077
        coefCl( 14 ) =  -0.393070909291043987E-02
        coefCl( 15 ) =  -0.836949836657998603E-01
        coefCl( 16 ) =  0.110168315452027800E-01
        coefCl( 17 ) =  0.150824303049820036E-02
        coefCl( 18 ) =  0.881618536264790364
        coefCl( 19 ) =  0.481211785132343143E-04
        coefCl( 20 ) =  -0.409104694234156455E-04

        coefCl( 21 ) =  18.7135851598321103
        coefCl( 22 ) =  -12.0071725597013632
        coefCl( 23 ) =  -0.696126329768568808
        coefCl( 24 ) =  0.331584083333504587E-01
        coefCl( 25 ) =  3.73688592263619501
        coefCl( 26 ) =  0.115335760662151395E-01
        coefCl( 27 ) =  -0.532606388417097756E-02
        coefCl( 28 ) =  -1.00194907052305204
        coefCl( 29 ) =  0.350793673471011232E-04
        coefCl( 30 ) =  0.178718937907040788E-04

        coefCl( 31 ) =  6.97070212448414850
        coefCl( 32 ) =  0.239324002996673835
        coefCl( 33 ) =  -0.167765520410741303
        coefCl( 34 ) =  0.426765987996296364E-02
        coefCl( 35 ) =  0.498343900330858791
        coefCl( 36 ) =  0.556317845435584917E-02
        coefCl( 37 ) =  -0.406111681215035134E-03
        coefCl( 38 ) =  -0.355788342182287276
        coefCl( 39 ) =  0.308622579496519578E-04
        coefCl( 40 ) =  0.401584499620915860E-05

        coefCl( 41 ) =  17.7330900312878157
        coefCl( 42 ) =  -13.3640606590466913
        coefCl( 43 ) =  -0.703449875330017682
        coefCl( 44 ) =  0.404241813901437749E-01
        coefCl( 45 ) =  3.97351038038355098
        coefCl( 46 ) =  0.104015001646859517E-01
        coefCl( 47 ) =  -0.650282143136824298E-02
        coefCl( 48 ) =  -1.05699997391849165
        coefCl( 49 ) =  0.373049859457823133E-04
        coefCl( 50 ) =  0.185789141245987270E-04

        coefCl( 51 ) =  18.1254655402286886
        coefCl( 52 ) =  -12.8344037916175324
        coefCl( 53 ) =  -0.698917858752339471
        coefCl( 54 ) =  0.281408560341282581E-01
        coefCl( 55 ) =  3.61348917129132330
        coefCl( 56 ) =  0.119769059165607920E-01
        coefCl( 57 ) =  -0.498127357334230293E-02
        coefCl( 58 ) =  -0.596754605130297988
        coefCl( 59 ) =  0.551949690664402050E-04
        coefCl( 60 ) =  0.156089903962416113E-04

C       Coefficients for CH3Br

        coefBr( 1 ) =  -4.57580255476205888
        coefBr( 2 ) =  -19.6014780840912515
        coefBr( 3 ) =  -0.210290832408275452
        coefBr( 4 ) =  0.319945750381702837E-01
        coefBr( 5 ) =  5.34718652348615375
        coefBr( 6 ) =  0.773031875767682831E-02
        coefBr( 7 ) =  -0.607318405917784375E-02
        coefBr( 8 ) =  -1.63091656610339597
        coefBr( 9 ) =  0.257565867995306501E-04
        coefBr( 10 ) =  0.379391509804299443E-04

        coefBr( 11 ) =  9.37016365556973341
        coefBr( 12 ) =  -4.50672966536395236
        coefBr( 13 ) =  -0.381264434305581712
        coefBr( 14 ) =  -0.468823558580684901E-02
        coefBr( 15 ) =  0.250616148533752484E-01
        coefBr( 16 ) =  0.607737728209738638E-02
        coefBr( 17 ) =  0.142320478500531419E-02
        coefBr( 18 ) =  0.439705349641294096
        coefBr( 19 ) =  0.284409836247023851E-04
        coefBr( 20 ) =  -0.228635256124437782E-04

        coefBr( 21 ) =  10.3598754795014969
        coefBr( 22 ) =  -6.53267359444823992
        coefBr( 23 ) =  -0.380835921229914665
        coefBr( 24 ) =  0.153772537374428017E-01
        coefBr( 25 ) =  1.87678143833666078
        coefBr( 26 ) =  0.601068355319096456E-02
        coefBr( 27 ) =  -0.263696635309326006E-02
        coefBr( 28 ) =  -0.413554079860463442
        coefBr( 29 ) =  0.874279264861569298E-05
        coefBr( 30 ) =  0.101212829393572074E-04

        coefBr( 31 ) =  3.96558687007669475
        coefBr( 32 ) =  -0.281781176181172488
        coefBr( 33 ) =  -0.148931882734184257
        coefBr( 34 ) =  0.249466208631324250E-02
        coefBr( 35 ) =  0.664720284786069593
        coefBr( 36 ) =  0.448251007974247336E-02
        coefBr( 37 ) =  -0.457178938586405895E-03
        coefBr( 38 ) =  -0.950073852590082313E-01
        coefBr( 39 ) =  0.176323085810440722E-04
        coefBr( 40 ) =  0.282168785160503574E-05

        coefBr( 41 ) =  9.71453830469256552
        coefBr( 42 ) =  -7.30951537510202076
        coefBr( 43 ) =  -0.383375621356974616
        coefBr( 44 ) =  0.167015890140731760E-01
        coefBr( 45 ) =  1.90297873608476431
        coefBr( 46 ) =  0.579622924272952308E-02
        coefBr( 47 ) =  -0.270671866424681759E-02
        coefBr( 48 ) =  -0.425815259434083770
        coefBr( 49 ) =  0.238020342924025799E-04
        coefBr( 50 ) =  0.576506789425034131E-05

        coefBr( 51 ) =  9.81662404177448522
        coefBr( 52 ) =  -7.16535094799607997
        coefBr( 53 ) =  -0.392231939501513971
        coefBr( 54 ) =  0.953875130337036164E-02
        coefBr( 55 ) =  1.64534337675065911
        coefBr( 56 ) =  0.579389987835394072E-02
        coefBr( 57 ) =  -0.145799286022983937E-02
        coefBr( 58 ) =  -0.354716058601023232
        coefBr( 59 ) =  0.152541519427672754E-04
        coefBr( 60 ) =  0.825921411444731959E-05

C      Spin-orbit transformation matrix

       SOC(1,1) = (.81649662d0,0.d0)
       SOC(2,1) = (.40824826d0,0.d0)
       SOC(3,1) = (0.d0,-.40824826d0)
       SOC(4,1) = (0.d0,0.d0)
       SOC(5,1) = (0.d0,0.d0)
       SOC(6,1) = (0.d0,0.d0)
       SOC(7,1) = (0.d0,0.d0)
       SOC(8,1) = (0.d0,0.d0)
       SOC(9,1) = (0.d0,0.d0)
       SOC(10,1) = (0.d0,0.d0)
       SOC(11,1) = (0.d0,0.d0) 
       SOC(12,1) = (0.d0,0.d0)    

       SOC(1,2) = (0.d0,0.d0)
       SOC(2,2) = (.70710672d0,0.d0)
       SOC(3,2) = (0.d0,.70710672d0)
       SOC(4,2) = (0.d0,0.d0)
       SOC(5,2) = (0.d0,0.d0)
       SOC(6,2) = (0.d0,0.d0)
       SOC(7,2) = (0.d0,0.d0)
       SOC(8,2) = (0.d0,0.d0)
       SOC(9,2) = (0.d0,0.d0)
       SOC(10,2) = (0.d0,0.d0)
       SOC(11,2) = (0.d0,0.d0)
       SOC(12,2) = (0.d0,0.d0) 

       SOC(1,3) = (-.57735027d0,0.d0)
       SOC(2,3) = (.57735027d0,0.d0)
       SOC(3,3) = (0.d0,-.57735027d0)
       SOC(4,3) = (0.d0,0.d0)
       SOC(5,3) = (0.d0,0.d0)
       SOC(6,3) = (0.d0,0.d0)
       SOC(7,3) = (0.d0,0.d0)
       SOC(8,3) = (0.d0,0.d0)
       SOC(9,3) = (0.d0,0.d0)
       SOC(10,3) = (0.d0,0.d0)
       SOC(11,3) = (0.d0,0.d0)
       SOC(12,3) = (0.d0,0.d0)    

       SOC(1,4) = (0.d0,0.d0)
       SOC(2,4) = (0.d0,0.d0)
       SOC(3,4) = (0.d0,0.d0)
       SOC(4,4) = (.62952645d0,0.d0)
       SOC(5,4) = (.76506570d0,0.d0)
       SOC(6,4) = (0.d0,.13553938d0)
       SOC(7,4) = (0.d0,0.d0)
       SOC(8,4) = (0.d0,0.d0)
       SOC(9,4) = (0.d0,0.d0)
       SOC(10,4) = (0.d0,0.d0)
       SOC(11,4) = (0.d0,0.d0)
       SOC(12,4) = (0.d0,0.d0)    
 
       SOC(1,5) = (0.d0,0.d0)
       SOC(2,5) = (0.d0,0.d0)
       SOC(3,5) = (0.d0,0.d0)
       SOC(4,5) = (-.57735027d0,0.d0)
       SOC(5,5) = (.57735027d0,0.d0)
       SOC(6,5) = (0.d0,-.57735027d0)
       SOC(7,5) = (0.d0,0.d0)
       SOC(8,5) = (0.d0,0.d0)
       SOC(9,5) = (0.d0,0.d0)
       SOC(10,5) = (0.d0,0.d0)
       SOC(11,5) = (0.d0,0.d0)
       SOC(12,5) = (0.d0,0.d0)    

       SOC(1,6) = (0.d0,0.d0)
       SOC(2,6) = (0.d0,0.d0)
       SOC(3,6) = (0.d0,0.d0)
       SOC(4,6) = (-.51996450d0,0.d0)
       SOC(5,6) = (.28520349d0,0.d0)
       SOC(6,6) = (0.d0,.80516823d0)
       SOC(7,6) = (0.d0,0.d0)
       SOC(8,6) = (0.d0,0.d0)
       SOC(9,6) = (0.d0,0.d0)
       SOC(10,6) = (0.d0,0.d0)
       SOC(11,6) = (0.d0,0.d0)
       SOC(12,6) = (0.d0,0.d0)    

       SOC(1,7) = (0.d0,0.d0)
       SOC(2,7) = (0.d0,0.d0)
       SOC(3,7) = (0.d0,0.d0)
       SOC(4,7) = (0.d0,0.d0)
       SOC(5,7) = (0.d0,0.d0)
       SOC(6,7) = (0.d0,0.d0)
       SOC(7,7) = (.62952637d0,0.d0)
       SOC(8,7) = (0.d0,-.76506575d0)
       SOC(9,7) = (0.d0,-.13553951d0)
       SOC(10,7) = (0.d0,0.d0)
       SOC(11,7) = (0.d0,0.d0)
       SOC(12,7) = (0.d0,0.d0)    
 
       SOC(1,8) = (0.d0,0.d0)
       SOC(2,8) = (0.d0,0.d0)
       SOC(3,8) = (0.d0,0.d0)
       SOC(4,8) = (0.d0,0.d0)
       SOC(5,8) = (0.d0,0.d0)
       SOC(6,8) = (0.d0,0.d0)
       SOC(7,8) = (.57735027d0,0.d0)
       SOC(8,8) = (0.d0,.57735027d0)
       SOC(9,8) = (0.d0,-.57735027d0)
       SOC(10,8) = (0.d0,0.d0)
       SOC(11,8) = (0.d0,0.d0)
       SOC(12,8) = (0.d0,0.d0)    

       SOC(1,9) = (0.d0,0.d0)
       SOC(2,9) = (0.d0,0.d0)
       SOC(3,9) = (0.d0,0.d0)
       SOC(4,9) = (0.d0,0.d0)
       SOC(5,9) = (0.d0,0.d0)
       SOC(6,9) = (0.d0,0.d0)
       SOC(7,9) = (.51996461d0,0.d0)
       SOC(8,9) = (0.d0,.28520336d0)
       SOC(9,9) = (0.d0,.80516821d0)
       SOC(10,9) = (0.d0,0.d0)
       SOC(11,9) = (0.d0,0.d0)
       SOC(12,9) = (0.d0,0.d0)    

       SOC(1,10) = (0.d0,0.d0)
       SOC(2,10) = (0.d0,0.d0)
       SOC(3,10) = (0.d0,0.d0)
       SOC(4,10) = (0.d0,0.d0)
       SOC(5,10) = (0.d0,0.d0)
       SOC(6,10) = (0.d0,0.d0)
       SOC(7,10) = (0.d0,0.d0)
       SOC(8,10) = (0.d0,0.d0)
       SOC(9,10) = (0.d0,0.d0)
       SOC(10,10) = (.57735027d0,0.d0)
       SOC(11,10) = (.57735027d0,0.d0)
       SOC(12,10) = (0.d0,-.57735027d0)

       SOC(1,11) = (0.d0,0.d0)
       SOC(2,11) = (0.d0,0.d0)
       SOC(3,11) = (0.d0,0.d0)
       SOC(4,11) = (0.d0,0.d0)
       SOC(5,11) = (0.d0,0.d0)
       SOC(6,11) = (0.d0,0.d0)
       SOC(7,11) = (0.d0,0.d0)
       SOC(8,11) = (0.d0,0.d0)
       SOC(9,11) = (0.d0,0.d0)
       SOC(10,11) = (0.d0,0.d0)
       SOC(11,11) = (.70710672d0,0.d0)
       SOC(12,11) = (0.d0,.70710672d0)

       SOC(1,12) = (0.d0,0.d0)
       SOC(2,12) = (0.d0,0.d0)
       SOC(3,12) = (0.d0,0.d0)
       SOC(4,12) = (0.d0,0.d0)
       SOC(5,12) = (0.d0,0.d0)
       SOC(6,12) = (0.d0,0.d0)
       SOC(7,12) = (0.d0,0.d0)
       SOC(8,12) = (0.d0,0.d0)
       SOC(9,12) = (0.d0,0.d0)
       SOC(10,12) = (.81649669d0,0.d0)
       SOC(11,12) = (-.40824819d0,0.d0)
       SOC(12,12) = (0.d0,.40824819d0)

      return
      end
c **********************************************************************
c **********************************************************************




c **********************************************************************
c **********************************************************************
      subroutine dd_g03(symbol,tnam,x,y,z,pemd,gpemd,nclu,
     $                                      mnclu,nsurf,mnsurf,ifail)


c INPUT
c
c SYMBOL(MNCLU) : Array of atomic symbols (H, C, etc...)
c X,Y,Z(MNCLU) :  Arrays of cartesian coordinates in bohr
c NCLU :          Number of atoms
c MNCLU :         Max number of atoms for declaring arrays
c 
c OUTPUT:
c
c V:               Potential energy in hartree
c DX,DY,DZ(MNCLU): Gradients in hartree/bohr
c IFAIL :          Flag giving info about QC failures
c                  (Not yet implemented for Gaussian)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,nsurf,mnsurf,ifail
      character*2 symbol(mnclu)
      character*80 string,tmpstr
      character*12 Estr
      character*7 Gstr
      double precision x(mnclu),y(mnclu),z(mnclu),xtmp(mnclu*3),v
      double precision dx(mnclu),dy(mnclu),dz(mnclu)
      double precision cfloat
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnclu,mnsurf,mnsurf)
      character*8 tnam

c read template & write QC input file
      open(unit=7,file=tnam)       ! template file
      open(unit=10,file='qc.in')   ! temporary input file
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
          do j=1,nclu
            write(10,fmt='(a2,2x,3f20.10)')symbol(j),x(j),y(j),z(j)
          enddo
        else
          write(10,fmt='(a80)')string
        endif
      enddo
      write(6,*)"QC TEMPLATE IS TOO LONG (> 100 lines)"
      stop
 100  close(7)
      close(10)

c do gaussian calculation
      call system('./g.x ')

c read the formatted checkpoint file
      open (8,file='Test.FChk')

c     get energy
      Estr='Total Energy'
      lstr=len(Estr)
 200  read (8,fmt='(a80)',end=220) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Estr) then
          tmpstr=string(46:80)
          v=cfloat(tmpstr)
          goto 299
        endif
      enddo
      goto 200

 220  write(6,*)"Couldn't find string '",estr,"' in Test.FChk"
      stop

 299  continue

c     get gradients
      Gstr='Cartesian Gradient'
      lstr=len(Gstr)
 300  read (8,fmt='(a80)',end=320) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Gstr) then
          j=0
          do while (j.lt.nclu*3)
            if (nclu*3-j.gt.5) then
              read(8,*)(xtmp(j+k),k=1,5)
              j=j+5
            else
c I think it goes x1,y1,z1,x2,...,zN
              read(8,*)(xtmp(j+k),k=1,nclu*3-j)
              j=nclu*3
            endif
          enddo
          goto 399
        endif
      enddo
      goto 300

 320  write(6,*)"Couldn't find string '",gstr,"' in Test.FChk"
      stop

 399  continue
      close(8)

      do i=1,nclu
        j=(i-1)*3
        dx(i)=xtmp(j+1)
        dy(i)=xtmp(j+2)
        dz(i)=xtmp(j+3)
      enddo

c organize things
      do i=1,nsurf
        pemd(i,i)=v
        do j=1,nclu
          gpemd(1,j,i,i)=dx(j)
          gpemd(2,j,i,i)=dy(j)
          gpemd(3,j,i,i)=dz(j)
        enddo
      enddo

      end
c **********************************************************************
c **********************************************************************



c **********************************************************************
c **********************************************************************
      subroutine dd_m06(symbol,tnam,x,y,z,pemd,gpemd,nclu,
     &                                      mnclu,nsurf,mnsurf,ifail)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,lstr2,idum,ithis(mnclu),nsurf,
     &   mnsurf,isurf
      integer ifail
      character*2 symbol(mnclu),cdum
      character*80 string,tmpstr
      character*21 Estr
      character*18 Gstr,Gstr2
      character*18 Astr
      double precision x(mnclu),y(mnclu),z(mnclu),v(mnsurf)
      double precision xk,dx1,dy1,dz1
      double precision gpemd(3,mnclu,mnsurf,mnsurf),pemd(mnsurf,mnsurf)
      double precision xtmp,ytmp,ztmp,total,dum
      double precision dx(mnclu,mnsurf),dy(mnclu,mnsurf),
     &                                  dz(mnclu,mnsurf)
      double precision cfloat
      double precision autoang,so,autoev
      character*8 tnam
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.211d0)

c read template & write QC input file
      open(unit=7,file=tnam)
      open(unit=10,file='qc.in')
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
c          write(10,*)"geomtyp=xyz"
c          write(10,*)"geometry"
c          write(10,*)"nosym"
c          write(10,*)"noorient"
          write(10,*)nclu
          write(10,*)"ANT Direct Dynamics Calculation"
          do j=1,nclu
            write(10,fmt='(a2,2x,3f20.10)')symbol(j),
     &       x(j)*autoang,y(j)*autoang,z(j)*autoang   ! default is Angstroms for geomtyp=xyz
          enddo
c          write(10,*)"end"
        else
          write(10,fmt='(a80)')string
        endif
      enddo
      write(6,*)"QC TEMPLATE IS TOO LONG (> 100 lines)"
      stop
 100  close(7)
      close(10)

c do molpro calculation
      call system('./m.x ')

c read the formatted checkpoint file
      open (8,file='qc.out')

c     molpro reorders the atoms, and I can't figure out how to make it stop, so...
      Astr='ATOMIC COORDINATES'
      lstr=len(Astr)
 150  read (8,fmt='(a80)',end=170) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Astr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            ithis(j)=-1
          enddo
          do j=1,nclu
            read(8,*)idum,cdum,dum,xtmp,ytmp,ztmp  ! read where Molpro reprints the reordered atoms
            do k=1,nclu
             total=dabs(xtmp-x(k))+dabs(ytmp-y(k))+dabs(ztmp-z(k))  ! check each (x,y,z) against input (both are in bohr)
             if (total.le.1.d-4) ithis(j)=k ! if the difference is small, then this is it 
c                                             (obviously this will fail for weird geometries with very small distances)
c                                             if everything works, then line number J in the molpro output corresponds
c                                             to ANT atom K
            enddo
          enddo
          goto 199
        endif
      enddo
      goto 150

 170  write(6,*)"Couldn't find string '",astr,"' in qc.out"
      ifail=1
      go to 999

 199  continue
      do j=1,nclu
        if (ithis(j).eq.-1) then
          write(6,*)"Problem with atom ordering in Molpro output file"
          ifail=2
          go to 999
        endif
      enddo

      DO ISURF=1,NSURF

c     get energy
      Estr='SETTING MOLPRO_ENERGY'
      lstr=len(Estr)
 200  read (8,fmt='(a80)',end=220) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Estr) then
          tmpstr=string(28:46)
          v(isurf)=cfloat(tmpstr)
          goto 299
        endif
      enddo
      goto 200

 220  write(6,*)"Couldn't find string '",estr,"' in qc.out"
      ifail=3
      go to 999

 299  continue

c     get gradients
      Gstr='MOLGRAD'
      lstr=len(Gstr)
      Gstr2='Total Energy'
      lstr2=len(Gstr2)
 300  read (8,fmt='(a80)',end=320) string
      do i=1,10
        if (string(i:i+lstr-1).eq.Gstr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)xk,dx1,dy1,dz1
            k=nint(xk)
            dx(ithis(k),isurf)=dx1
            dy(ithis(k),isurf)=dy1
            dz(ithis(k),isurf)=dz1
          enddo
          goto 399
        elseif (string(i:i+lstr2-1).eq.Gstr2) then
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)xk,dx1,dy1,dz1
            k=nint(xk)
            dx(ithis(k),isurf)=dx1
            dy(ithis(k),isurf)=dy1
            dz(ithis(k),isurf)=dz1
          enddo
          goto 399
        endif
      enddo
      goto 300

 320  continue
      ifail=4
      go to 999

 399  continue

      ENDDO

c organize things
      do i=1,nsurf
        pemd(i,i)=v(i)
        do j=1,nclu
          gpemd(1,j,i,i)=dx(j,i)
          gpemd(2,j,i,i)=dy(j,i)
          gpemd(3,j,i,i)=dz(j,i)
        enddo
      enddo

c off-diagonal elements (set to zero here)
      do i=1,nsurf
      do j=i+1,nsurf
        pemd(i,j)=0.d0/autoev
        pemd(j,i)=pemd(i,j)
        do k=1,nclu
          gpemd(1,k,i,j)=0.d0
          gpemd(2,k,i,j)=0.d0
          gpemd(3,k,i,j)=0.d0
          gpemd(1,k,j,i)=0.d0
          gpemd(2,k,j,i)=0.d0
          gpemd(3,k,j,i)=0.d0
        enddo
      enddo
      enddo

 999  continue
      close(8)
      return

      end
C**********************************************************************
C**********************************************************************




C**********************************************************************
C**********************************************************************
C CFLOAT
C
      double precision function cfloat(string)
C
      implicit double precision(a-h,o-z)
      character*80 string,numbe
      character ch
      logical lexp,ldec

c AJ
      integer fu6
      fu6=6
C

      LEXP = .FALSE.
      LDEC = .FALSE.
      LENGTH = LEN(STRING)
      IF (LENGTH .EQ. 0) THEN
         CFLOAT = 0.0D0
         RETURN
      ENDIF
C
C     Find the first nonblank character
C
      I = 1
10    IF (STRING(I:I) .EQ. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 10
      ENDIF
C
C     If it is a blank string set function to zero
C
      IF (I .GT. LENGTH) THEN
         CFLOAT = 0.0D0
         RETURN
      ENDIF
      IBEG = I
C
C     Find the first blank character after the number
C
      I = IBEG+1
20    IF (STRING(I:I) .NE. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 20
      ENDIF
      IEND = I-1
C
C     Stripe the blanks before and after the number
C
      NUMBE = STRING(IBEG:IEND)
      LENGTH = IEND - IBEG + 1
C   
C     Make sure there is no blank left
C
      IF (INDEX(NUMBE,' ') .LE. LENGTH) THEN
         WRITE(FU6,1000) STRING
         STOP 'CFLOAT 1'
      ENDIF
C
C     Find the decimal point
C
      IDEC = INDEX(NUMBE,'.')
      IF (IDEC .NE. 0) LDEC = .TRUE.
C
C     Find the exponential symbol
C
      IUE = INDEX(NUMBE,'E')
      ILE = INDEX(NUMBE,'e')
      IUD = INDEX(NUMBE,'D')
      ILD = INDEX(NUMBE,'d')
      ISUM = IUE + ILE + IUD + ILD
      IEXP = MAX0(IUE,ILE,IUD,ILD)
      IF (ISUM .GT. IEXP) THEN
         WRITE(FU6,1000) STRING
         STOP 'CFLOAT 2'
      ENDIF
      IF (IEXP .NE. 0) THEN
         LEXP = .TRUE.
      ELSE
         IEXP = LENGTH + 1
      ENDIF
C
      IF (.NOT. LDEC) IDEC = IEXP
C
C     Get the number before decimal
C
      IBEG = 2
      IF (NUMBE(1:1) .EQ. '+') THEN
         SIGN = 1.0D0
      ELSEIF(NUMBE(1:1) .EQ. '-') THEN
         SIGN = -1.0D0
      ELSE
         SIGN = 1.0D0
         IBEG = 1
      ENDIF
      IF (IBEG .EQ. IEXP) THEN
         F1 = 1.0D0
      ELSE
         F1 = 0.0D0
      ENDIF
      DO 50 I = IBEG,IDEC-1
         CH = NUMBE(I:I)
         IF (CH .GE. '0' .AND. CH .LE. '9') THEN
            N = ICHAR(CH) - ICHAR('0')
            F1 = F1 * 10.0D0 + DBLE(N)
         ELSE
            WRITE(FU6,1000) STRING
            STOP 'CFLOAT 3'
         ENDIF
50    CONTINUE
C
C     Get the number after decimal 
C
      F2 = 0.0D0
      IF (LDEC) THEN
         J = 0
         DO 60 I = IDEC+1,IEXP-1
            CH = NUMBE(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F2 = F2 * 10.0D0 + DBLE(N)
               J = J + 1
            ELSE
               WRITE(FU6,1000) STRING
               STOP 'CFLOAT 4'
            ENDIF
60       CONTINUE
         F2 = F2 / 10.0D0 ** DBLE(J)
      ENDIF
C
C    Get the exponent
C
      ESIGN = 1.0D0
      F3 = 0.0D0
      IF (LEXP) THEN 
         IBEG = IEXP + 2
         IF (NUMBE(IEXP+1:IEXP+1) .EQ. '+') THEN
            ESIGN = 1.0D0
         ELSEIF(NUMBE(IEXP+1:IEXP+1) .EQ. '-') THEN
            ESIGN = -1.0D0
         ELSE
            ESIGN = 1.0D0
            IBEG = IEXP + 1
         ENDIF
         DO 70 I = IBEG,LENGTH
            CH = NUMBE(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F3 = F3 * 10.0D0 + DBLE(N)
            ELSE
               WRITE(FU6,1000) STRING
               STOP 'CFLOAT 5'
            ENDIF
70       CONTINUE
      ENDIF 
C
      CFLOAT = (SIGN * (F1 + F2)) * 10.0D0 ** (ESIGN*F3)
C
      RETURN
C
1000  FORMAT(/1X,'Illegal number: ',A80)
C
      END

C**********************************************************************
C**********************************************************************
