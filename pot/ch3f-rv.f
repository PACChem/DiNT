      subroutine pot(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)

      implicit none

      integer mnsurf,nsurf,mnat,nat,n3,i,ii,j,k
      double precision autoang
      parameter(autoang=0.52917706d0)
      character*2 symb(mnat)
      double precision sing3(3,3),gsing3(3,mnat,3,3),trip3(3,3),
     & gtrip3(3,mnat,3,3)
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf)
      double precision x(mnat),y(mnat),z(mnat)
      double precision autokcal,singsig,singpi,tripsig,trippi,ajsigma,
     & ajpi,c1,c2,trip3a,trip3b
      double precision pi,a1,a2,a3,a4,a5,a6,a7,a8,th,phi
      complex*16 SOC(12,12),omega(12,12),SOCadj(12,12),interm(12,12)
      double precision lambda(12,12),aa,inte
      double precision delta
      double precision z2,z3,zz2,zz3,zz4,zz5,zz6,zz7
      double precision Fp(11),Fp2(11),Fp3(11),coefF(44),coefCl(44),
     & coefBr(44)
      double precision xcmh1,ycmh1,zcmh1,xcmh2,ycmh2,zcmh2,xcmh3,ycmh3,
     & zcmh3,intx,inty,intz,nx,ny,nz,px,py,pz,rcmh1,rcmh2,rcmh3,phi1,
     & phi2,phi3,pmod,xh1h2,yh1h2,zh1h2,xh2h3,yh2h3,zh2h3,th1,xfcm,
     & yfcm,zfcm,rfcm,xcm,ycm,zcm,rh1h2,rh2h3
      double precision angh1cmh2,angh1cmh3,angh2cmh3
      double precision gajsigma(3,mnat),gajpi(3,mnat),gtrip3a(3,mnat),
     &gtrip3b(3,mnat),inter(3,mnat)
      integer itimes,l,m
      double precision glambda(3,mnat,12,12),ginterm(3,mnat,12,12) 

      common /pot_coe/ coefF,coefCl,coefBr
      common /SOC_matrix/ SOC

      pi = 3.14159274
      autokcal = 627.509470

c compute and return the three singlet energies and gradients
c     sing3 contains the 3 singlet energies
c     gsing3 contians the gradients of the 3 singlet energies
c     the off diagonal elements of sing3 and gsing3 are not used
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
c       print *,"gE2(",i,") = ",gsing3(1,i,2,2),gsing3(2,i,1,1),
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

C       Compute energies of triplet states from singlet energies
C       First decide if sing3(1,1) is "Sigma", and sing3(2,2),sing3(3,3) are "Pi" OR
C       sing3(1,1),sing3(2,2) are "Pi", and sing3(3,3) is "Sigma"

        itimes = 0

        if(abs(sing3(3,3)-sing3(2,2)).gt.
     &   abs(sing3(2,2)-sing3(1,1))) then
         inte = sing3(1,1)
         sing3(1,1) = sing3(3,3)
         sing3(3,3) = inte
        do l=1,3
         do m=1,mnat
          inter(l,m) = gsing3(l,m,1,1) 
          gsing3(l,m,1,1) = gsing3(l,m,3,3)
          gsing3(l,m,3,3) = inter(l,m)
         enddo
        enddo
        endif

22      if(itimes.eq.1) then
          inte = sing3(1,1)
          sing3(1,1) = sing3(3,3)
          sing3(3,3) = inte
        do l=1,3
         do m=1,mnat
          inter(l,m) = gsing3(l,m,1,1) 
          gsing3(l,m,1,1) = gsing3(l,m,3,3)
          gsing3(l,m,3,3) = inter(l,m)
         enddo
        enddo
        endif

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

C     With the last formula the three phi are not measured always in the same
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

C     Calculate energy of singlets and triplets from the fits

        z2 = (-cos(4.0*pi*th/360.0))
        z3 = (-cos(4.0*pi*th/360.0) + 1.0)
        zz2 = cos(6.0*pi*phi/360.0)
        zz3 = cos(12.0*pi*phi/360.0)
        zz4 = cos(6.0*pi*(phi+40.0)/360.0)
        zz5 = cos(12.0*pi*(phi+40.0)/360.0)
        zz6 = cos(6.0*pi*(phi-40.0)/360.0)
        zz7 = cos(12.0*pi*(phi-40.0)/360.0)

        Fp(1) = z2
        Fp(2) = 0.5*(5.0*z2**3 - 3.0*z2) * zz2
        Fp(3) = (1.0/8.0)*(63.0*z2**5 - 70.0*z2**3 + 15.0*z2)*zz2
        Fp(4) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz2
        Fp(5) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz3
        Fp(6) = z3
        Fp(7) = 0.5*(5.0*z3**3 - 3.0*z3) * zz2
        Fp(8) = (1.0/8.0)*(63.0*z3**5 - 70.0*z3**3 + 15.0*z3)*zz2
        Fp(9) = (1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz2
        Fp(10)=(1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz3
        Fp(11)=(-cos((4.0*pi*th/360.0))**2.0)

        Fp2(1) = z2
        Fp2(2) = 0.5*(5.0*z2**3 - 3.0*z2) * zz4
        Fp2(3) = (1.0/8.0)*(63.0*z2**5 - 70.0*z2**3 + 15.0*z2)*zz4
        Fp2(4) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz4
        Fp2(5) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz5
        Fp2(6) = z3
        Fp2(7) = 0.5*(5.0*z3**3 - 3.0*z3) * zz4
        Fp2(8) = (1.0/8.0)*(63.0*z3**5 - 70.0*z3**3 + 15.0*z3)*zz4
        Fp2(9) = (1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz4
        Fp2(10)=(1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz5
        Fp2(11)=(-cos((4.0*pi*th/360.0))**2.0)

        Fp3(1) = z2
        Fp3(2) = 0.5*(5.0*z2**3 - 3.0*z2) * zz6
        Fp3(3) = (1.0/8.0)*(63.0*z2**5 - 70.0*z2**3 + 15.0*z2)*zz6
        Fp3(4) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz6
        Fp3(5) = (1.0/16.0)*(429.0*z2**7 - 693.0*z2**5 + 315.0*z2**3
     &  -35.0*z2)*zz7
        Fp3(6) = z3
        Fp3(7) = 0.5*(5.0*z3**3 - 3.0*z3) * zz6
        Fp3(8) = (1.0/8.0)*(63.0*z3**5 - 70.0*z3**3 + 15.0*z3)*zz6
        Fp3(9) = (1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz6
        Fp3(10)=(1.0/16.0)*(429.0*z3**7 - 693.0*z3**5 + 315.0*z3**3
     &  -35.0*z3)*zz7
        Fp3(11)=(-cos((4.0*pi*th/360.0))**2.0)

        singsig = 0.d0
        singpi = 0.d0
        tripsig = 0.d0
        trippi = 0.d0

        if(symb(5).eq."f") then

        do i=1,11
       singsig = singsig + 1.1d0*Fp(i)*coefF(i)*
     & sin(2.0*pi*th/360.0)**2 +
     & (Fp2(i)*coefF(i) + Fp(i)*coefF(i) + Fp3(i)*coefF(i))
     & *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        singpi = singpi + Fp(i)*coefF(i+11)*sin(2.0*pi*th/360.0)**2 + 
     &  (Fp2(i)*coefF(i+11) + Fp(i)*coefF(i+11) + Fp3(i)*coefF(i+11))
     &  *cos(2.0*pi*th/360.0)**2/3.0  
        enddo

        do i=1,11
       tripsig=tripsig + Fp(i)*coefF(i+22)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefF(i+22) + Fp(i)*coefF(i+22) + Fp3(i)*coefF(i+22))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        trippi = trippi + Fp(i)*coefF(i+33)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefF(i+33) + Fp(i)*coefF(i+33) + Fp3(i)*coefF(i+33))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        elseif(symb(5).eq."cl") then

        do i=1,11
       singsig = singsig + Fp(i)*coefCl(i)*sin(2.0*pi*th/360.0)**2 +
     & (Fp2(i)*coefCl(i) + Fp(i)*coefCl(i) + Fp3(i)*coefCl(i))
     & *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        singpi = singpi + Fp(i)*coefCl(i+11)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefCl(i+11) + Fp(i)*coefCl(i+11) + Fp3(i)*coefCl(i+11))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
       tripsig=tripsig + Fp(i)*coefCl(i+22)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefCl(i+22) + Fp(i)*coefCl(i+22) + Fp3(i)*coefCl(i+22))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        trippi = trippi + Fp(i)*coefCl(i+33)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefCl(i+33) + Fp(i)*coefCl(i+33) + Fp3(i)*coefCl(i+33))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        elseif(symb(5).eq."br") then

        do i=1,11
       singsig = singsig + Fp(i)*coefBr(i)*
     & sin(2.0*pi*th/360.0)**2 +
     & (Fp2(i)*coefBr(i) + Fp(i)*coefBr(i) + Fp3(i)*coefBr(i))
     & *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        singpi = singpi + Fp(i)*coefBr(i+11)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefBr(i+11) + Fp(i)*coefBr(i+11) + Fp3(i)*coefBr(i+11))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
       tripsig=tripsig + Fp(i)*coefBr(i+22)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefBr(i+22) + Fp(i)*coefBr(i+22) + Fp3(i)*coefBr(i+22))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        do i=1,11
        trippi = trippi + Fp(i)*coefBr(i+33)*sin(2.0*pi*th/360.0)**2 +
     &  (Fp2(i)*coefBr(i+33) + Fp(i)*coefBr(i+33) + Fp3(i)*coefBr(i+33))
     &  *cos(2.0*pi*th/360.0)**2/3.0
        enddo

        endif

C     Apply VB expression for triplets as functions of singlets

      c1 = (singsig + 2.d0*singpi - 7.d0*tripsig + 6.d0*trippi)/
     &  (5.d0*(singsig - tripsig))
      c2 = (singsig + 2.d0*singpi + 3.d0*tripsig - 4.d0*trippi)/
     &  (5.d0*(singpi - trippi))

        ajsigma = (3.d0*c2*sing3(1,1) - (4.d0*c2 - 2.d0)*sing3(2,2))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )
        ajpi = ( (1.d0+c1)*sing3(2,2) - (2.d0*c1 - 1.d0)*sing3(1,1))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )

        trip3(2,2)= sing3(2,2) - 2.d0*ajpi
        trip3a = sing3(1,1) - 2.d0*ajsigma

        ajsigma = (3.d0*c2*sing3(1,1) - (4.d0*c2 - 2.d0)*sing3(3,3))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )
        ajpi = ( (1.d0+c1)*sing3(3,3) - (2.d0*c1 - 1.d0)*sing3(1,1))/
     &  ( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*(4.d0*c2 - 2.d0) )

        trip3(3,3)= sing3(3,3) - 2.d0*ajpi
        trip3b = sing3(1,1) - 2.d0*ajsigma

        trip3(1,1) = (trip3a + trip3b)/2.d0

C     If one or more of the triplets has a negative energy below -1 kcal/mol
C     revise criterion to classify 1Sigma, 1Pi

       if((trip3(1,1).lt.-0.0016. or .trip3(2,2).lt.-0.0016. or.
     & trip3(3,3).lt.-0.0016). and .itimes.eq.0) then
        itimes = itimes + 1
        goto 22
       else
        continue
       endif

C     Calculate gradients of triplet states

      do i=1,3
       do j=1,mnat

        gajsigma(i,j) = (3.d0*c2*gsing3(i,j,1,1) - (4.d0*c2 - 2.d0)*
     &  gsing3(i,j,2,2))/( 3.d0*c2*(c1+1.d0)-(2.d0*c1 - 1.d0)*
     &  (4.d0*c2 - 2.d0))
        gajpi(i,j) = ( (1.d0+c1)*gsing3(i,j,2,2) - (2.d0*c1 - 1.d0)*
     &  gsing3(i,j,1,1))/( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     &  (4.d0*c2 - 2.d0) )

        gtrip3(i,j,2,2)= gsing3(i,j,2,2) - 2.d0*gajpi(i,j)
        gtrip3a(i,j) = gsing3(i,j,1,1) - 2.d0*gajsigma(i,j)

        gajsigma(i,j) = (3.d0*c2*gsing3(i,j,1,1) - (4.d0*c2 - 2.d0)*
     & gsing3(i,j,3,3))/( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     & (4.d0*c2 - 2.d0) )
        gajpi(i,j) = ( (1.d0+c1)*gsing3(i,j,3,3) - (2.d0*c1 - 1.d0)*
     & gsing3(i,j,1,1))/( 3.d0*c2*(c1+1.d0) - (2.d0*c1 - 1.d0)*
     & (4.d0*c2 - 2.d0) )

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
        print *,(pemd(i,j)*27.211d0,j=1,12)
       enddo

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
          IF (.FALSE.) THEN
c         handle failures
          if (ifail.eq.0) then
c           everything OK
          elseif (ifail.eq.3) then
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





c **********************************************************************
c **********************************************************************
      subroutine prepot
     
      double precision coefF(44),coefCl(44),coefBr(44)
      complex*16 SOC(12,12)
      common /pot_coe/ coefF,coefCl,coefBr
      common /SOC_matrix/ SOC

C       Coefficients for CH3F

        coefF(1) = 34.2568984972709600
        coefF(2) = 0.815803987099408268
        coefF(3) = 1.97638309369391441
        coefF(4) = -0.467651031822151009
        coefF(5) = 0.442298752716083488
        coefF(6) = -16.6983932137480053
        coefF(7) = 0.0972974714932482443
        coefF(8) = 0.0295736634971844638
        coefF(9) = 0.00201476235105344827
        coefF(10) = -0.000394534344577274578
        coefF(11) = -8.17196524417184023

        coefF(12) = -8.01495034277106022
        coefF(13) = 4.65186302191481982
        coefF(14) = -3.03285247004536451
        coefF(15) = -0.761206819820244518
        coefF(16) = 0.389789912124984983
        coefF(17) = 14.7676201963799851
        coefF(18) = 4.51361229607966052
        coefF(19) = -0.689309678388396074
        coefF(20) = 0.0273193262919496534
        coefF(21) = 0.00128536970742070379
        coefF(22) = -3.86152414222647167

        coefF(23) = -8.35508189195861561
        coefF(24) = 0.618488042997546894
        coefF(25) = -0.0703274746469362971
        coefF(26) = -0.240785714840227150
        coefF(27) = 0.221654799589455609
        coefF(28) = 7.26005279951082372
        coefF(29) = 0.661228846021161387
        coefF(30) = -0.0800461550730983823
        coefF(31) = 0.00218531672986167292
        coefF(32) = 0.000174343026190504742
        coefF(33) = -0.583406765980881259

        coefF(34) = -4.90731606249476382
        coefF(35) = 5.69806550202457274
        coefF(36) = -2.52820616630996797
        coefF(37) = -0.552398026852378266
        coefF(38) = 0.245478368637441074
        coefF(39) = 13.9648150841097998
        coefF(40) = 5.00254308248387503
        coefF(41) = -0.801481681162784731
        coefF(42) = 0.0359223033624843438
        coefF(43) = 0.000812372534080395476
        coefF(44) = -6.51448459478853170

C       Coefficients for CH3Cl

       coefCl(1) =  34.4816438803044605
       coefCl(2) =  1.50543043375334840
       coefCl(3) =  1.14964895827507418
       coefCl(4) =  -0.769029936328295194
       coefCl(5) =  0.320727016371964935
       coefCl(6) =  -16.5745372528336929
       coefCl(7) =  0.773471854903561851
       coefCl(8) =  -0.100846009060196112
       coefCl(9) =  0.00691045795048445856
       coefCl(10) =  -0.000337363581357153636
       coefCl(11) =  -7.15847038720442885

       coefCl(12) =  -9.09748124995457630
       coefCl(13) =  3.51463448506935761
       coefCl(14) =  -4.14219636406842540
       coefCl(15) =  2.96113861194362427
       coefCl(16) =  -0.613130311023586838
       coefCl(17) =  14.3916561506331142
       coefCl(18) =  4.58415250374828531
       coefCl(19) =  -0.927505177355741450
       coefCl(20) =  0.0484160982123369027
       coefCl(21) =  0.000285793700553368956
       coefCl(22) =  -1.86658962673352891

       coefCl(23) =  -7.27211638118602988
       coefCl(24) =  0.790263160688133293
       coefCl(25) =  0.356780859375622983
       coefCl(26) =  -0.674123370530858912
       coefCl(27) =  0.317672312119389122
       coefCl(28) =  6.70995726938073922
       coefCl(29) =  0.477669292689389346
       coefCl(30) =  -0.0213914071848136239
       coefCl(31) =  -0.00146910939343289793
       coefCl(32) =  0.000362547871052952759
       coefCl(33) =  -0.904105579031015960

       coefCl(34) =  -5.22177617036604147
       coefCl(35) =  4.38357174869286670
       coefCl(36) =  -2.00366848950089516
       coefCl(37) =  -0.0685577036387190614
       coefCl(38) =  0.262192714482998135
       coefCl(39) =  12.8880349196509396
       coefCl(40) =  3.86114671102057860
       coefCl(41) =  -0.621503733103910383
       coefCl(42) =  0.0267199748498283055E-01
       coefCl(43) =  0.000524914567991625219
       coefCl(44) =  -4.81289522882179543

C       Coefficients for CH3Br

       coefBr(1) =  23.5531108157567566
       coefBr(2) =  0.990896905770226710
       coefBr(3) =  0.646188622445623917
       coefBr(4) =  -0.535585787506201116
       coefBr(5) =  0.208134161638168114
       coefBr(6) =  -11.2776932130959828
       coefBr(7) =  0.591402710549134292
       coefBr(8) =  -0.0755785688601180450
       coefBr(9) =  0.00449664431021431967
       coefBr(10) =  -0.0000675892928135768178
       coefBr(11) =  -4.44888048265373470

       coefBr(12) =  -4.83652900629097005
       coefBr(13) =  2.29469017157383304
       coefBr(14) =  -1.47528536671117516
       coefBr(15) =  0.392641468357126711
       coefBr(16) =  0.0120295615163227530
       coefBr(17) =  7.90794029057960213
       coefBr(18) =  2.46930086188693254
       coefBr(19) =  -0.444114956778328407
       coefBr(20) =  0.0208162945869615340
       coefBr(21) =  0.000552687284838069637
       coefBr(22) =  -1.71610585724866005

       coefBr(23) =  -3.71912319141572523
       coefBr(24) =  0.630234985570018225
       coefBr(25) =  -0.380189404614719950
       coefBr(26) =  -0.144218433925908690
       coefBr(27) =  0.110144379988511940
       coefBr(28) =  3.62212814669608418
       coefBr(29) =  0.608102009696999013
       coefBr(30) =  -0.0785145455047682045
       coefBr(31) =  0.00248300245884082536
       coefBr(32) =  0.000112038001449138553
       coefBr(33) =  -0.831411406860692814

       coefBr(34) =  -2.76620264583337283
       coefBr(35) =  2.36600993875504351
       coefBr(36) =  -1.14215806467275205
       coefBr(37) =  -0.0619528798245787965
       coefBr(38) =  0.182996607967677694
       coefBr(39) =  7.10941723108774948
       coefBr(40) =  2.12566700217732496
       coefBr(41) =  -0.340167646464717366
       coefBr(42) =  0.0143215934809430671
       coefBr(43) =  0.000322817111954396040
       coefBr(44) =  -2.55569206244497771

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

      print *,"PREPOT CALLED"

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
