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

      v=0.d0
      v2=0.d0
      v3=0.d0
      do i=1,natom
      symb2(i)="xx"
      symb3(i)="xx"
      x2(i)=0.d0
      y2(i)=0.d0
      z2(i)=0.d0
      x3(i)=0.d0
      y3(i)=0.d0
      z3(i)=0.d0
      dvdz2(i)=0.d0
      dvdy2(i)=0.d0
      dvdx2(i)=0.d0
      dvdz3(i)=0.d0
      dvdy3(i)=0.d0
      dvdx3(i)=0.d0
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do i=1,natom
      at(i)=0
      if ((symb(i).eq."H").or.
     &    (symb(i).eq."h"))  at(i)=1       ! hydrogen
      if ((symb(i).eq."C").or.
     &    (symb(i).eq."c"))  at(i)=2       ! carbon
      if ((symb(i).eq."O").or.
     &    (symb(i).eq."o"))  at(i)=4       ! oxygen
      if ((symb(i).eq."HO").or.
     &    (symb(i).eq."Ho"))  at(i)=5       ! OH H
      if ((symb(i).eq."He").or.
     &    (symb(i).eq."he").or.
     &    (symb(i).eq."HE")) at(i)=21      ! helium
      if ((symb(i).eq."Ar").or.
     &    (symb(i).eq."ar").or.
     &    (symb(i).eq."AR")) at(i)=23      ! argon
      if (at(i).eq.0) then ! atom not found
           write(6,*)"Atom # ",i," (",symb(i),") not found"
           stop
      endif
      enddo

      natom2=0
      natom3=0
      do i=1,natom
        if (at(i).le.20) then ! collect target atoms
          natom2=natom2+1
          x2(natom2)=x(i)
          y2(natom2)=y(i)
          z2(natom2)=z(i)
          symb2(natom2)=symb(i)
          at2(natom2)=at(i)
        else ! collect bath atoms
          natom3=natom3+1
          x3(natom3)=x(i)
          y3(natom3)=y(i)
          z3(natom3)=z(i)
          symb3(natom3)=symb(i)
          at3(natom3)=at(i)
        endif
      enddo
c      print *,natom,natom2,natom3
      if (natom3.ne.0) 
     &   call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
      if (natom2.ne.0) 
     &   call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
     &          natom2,maxatom)
      if (natom3.ne.0) 
     &   call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      v=v+v2+v3

      natom2=0
      natom3=0
      do i=1,natom
        if (at(i).le.20) then
        natom2=natom2+1
c        print *,i,dvdx2(i),dvdy2(i),dvdz2(i)
        dvdx(i)=dvdx(i)+dvdx2(natom2)
        dvdy(i)=dvdy(i)+dvdy2(natom2)
        dvdz(i)=dvdz(i)+dvdz2(natom2)
        else
        natom3=natom3+1
        dvdx(i)=dvdx(i)+dvdx3(natom3)
        dvdy(i)=dvdy(i)+dvdy3(natom3)
        dvdz(i)=dvdz(i)+dvdz3(natom3)
        endif
      enddo

      return

      end


      subroutine bath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      integer at(maxatom)

      v=0.
      do i=1,natom
        dvdx(i)=0.d0
        dvdy(i)=0.d0
        dvdz(i)=0.d0
      enddo

      if (natom.eq.1) return

      if (natom.gt.2) then
         print *,"Can't handle more than 2 bath atoms"
         stop
      endif

      if (natom.eq.2) then
      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)  ! au

        if (at(1).eq.25.and.at(2).eq.25) then
!       H2 bath
!       From Hack's fit (eq 8 in Hack, Truhlar, JCP 110, 4315 (1999))
!                        to Kolos and Wolniewicz JCP 43, 2429 (1965)
!       Rmin = 1.40121 au, Vmin = -4.74772265 eV relative to H+H
        c1=139.7160d0        ! eV
        c2=-123.8978d0       ! eV / bohr
        c3=3.4031d0          ! 1 / bohr
        c4=-6.8725d0         ! eV / bohr**2
        c5=-23.0440d0        ! eV / bohr**3
        c6=2.032d0           ! 1 / bohr

        v=(c1+c2*rr)*dexp(-c3*rr)
     &   +(c4+c5*rr)*dexp(-c6*rr)*rr**2
c       move zero from asymptote to minimum
        v=v+4.74772265
c        print *,rr,v
        v=v/autoev

        dvdr=((c1+c2*rr)*(-c3)+c2)*dexp(-c3*rr)
     &      +((c4+c5*rr)*(-c6)+c5)*dexp(-c6*rr)*rr**2
     &       +(c4+c5*rr)*dexp(-c6*rr)*rr*2.d0
        dvdr=dvdr/autoev

        elseif (at(1).eq.26.and.at(2).eq.26) then
!       N2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       agrees reasonably well with more complicated form of LeRoy (JCP 125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, PRL 53, 159 (1984), as quoted by LeRoy)
        re=1.097679d0 ! exp
        c1=2.68872341 ! my fit
        c2=0.240070803
        c3=0.472261727

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2-de    ! A and cm-1
        v=v/autocmi  ! convert to au

        dvdr = 2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)*beta
        dvdr=dvdr*autocmi/autocmi  ! convert to au

!       elseif (at(1).eq.??.and.at(2).eq.??) then
!       OTHER DIATOMIC BATHS HERE
        else
        print *,"Don't know this diatomic bath"
        stop
        endif

      dvdx(1) =  dvdr*dx/rr
      dvdx(2) = -dvdr*dx/rr
      dvdy(1) =  dvdr*dy/rr
      dvdy(2) = -dvdr*dy/rr
      dvdz(1) =  dvdr*dz/rr
      dvdz(2) = -dvdr*dz/rr

      endif

      return
      end



      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c Rare Gas exp6 potential subroutine
c loops over geometry and looks for Rg-X interactions
c returns the full Rg-target intermolecular potential and its derivatives

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)
      integer at(maxatom)
      logical troya,cutoff
      character dum
      dimension ccc(100)

      save iprepot2
      if (iprepot2.eq.0) then
c     set up the structure and mechanics calculation
      iprepot2=1
      open(22,file='coef.dat')
      read(22,*)ncoef
      do i=1,ncoef
      read(22,*)dum,ccc(i)
      enddo
      return
      endif

      v1=0.d0
      v=0.d0
      do i=1,natom
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do 1 i=1,natom
      do 2 j=i+1,natom

      m1=min(at(i),at(j))
      m2=max(at(i),at(j))
      troya=.false.   ! do or don't use Troya's form
      cutoff=.false.   ! do or don't use cutoff

      if (m1.ge.21) then ! two rare gases, skip this pair
         go to 2
      endif
      if (m2.le.20) then ! no rare gas, skip this pair
         go to 2
      endif

      if (m2.eq.21) then  ! He-
        if     (m1.eq.1) then !    H
        aa=   6.3939634388256481     
        bb=  0.22434766685995056     
        cc=   4.9925229651783809     
        rrc=   2.8015991698965426     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
        aa=   6.9955442976165045     
        bb=  0.22477797784356213     
        cc=   4.3569750053407397     
        rrc=   2.2914822840052489     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.4) then !    O
        aa=   7.7865840632343515     
        bb=  0.19954222235786007     
        cc=   6.6693929868465220     
        rrc=   3.5884578997161780     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.5) then !    OH H
        aa=   4.9993896298104801
        bb=  0.29424726096377452
        cc=   1.9019135105441449
        rrc=   19.194006164738916     
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a He-? interaction"
          stop
        endif
      elseif (m2.eq.23) then  ! Ar-
        if     (m1.eq.1) then !    H
        aa=   6.4374217963194678        ! original
        bb=  0.28550370799890135     
        cc=   7.3473311563463239     
        rrc=   2.8850978118228703     
        ix=1
        aa=ccc(ix)
        bb=ccc(ix+1)
        cc=ccc(ix+2)
        rrc=ccc(ix+3)
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
        aa=   7.5285195471053195     
        bb=  0.26272774437696461     
        cc=   6.7756889553514208     
        rrc=   3.1267433698538163     
c       aa=   8.4549059050348472
c       bb=  0.22670596532155032
c       cc=   8.3936305075403297
c       rrc=  0.87450288084102457
        ix=5
        aa=ccc(ix)
        bb=ccc(ix+1)
        cc=ccc(ix+2)
        rrc=ccc(ix+3)
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.4) then !    O
        aa=   8.9658192693868841     
        bb=  0.19998474074526201     
        cc=   7.5914487136448257     
        rrc=   2.7619861445966980     
c       aa=   6.9569004059739488
c       bb=  0.31582372778194251
c       cc=   8.9676265789396208
c       rrc=   2.4929106053279413
        ix=9
        aa=ccc(ix)
        bb=ccc(ix+1)
        cc=ccc(ix+2)
        rrc=ccc(ix+3)
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.5) then !    OH H
        aa=   5.9999084444715720     
        bb=  0.28738670003357036     
        cc=   7.5051118503372294     
        rrc=   2.9338053529465622     
c       aa=   6.0639100577323024
c       bb=  0.29138109446203270
c       cc=   5.8219878706767902
c       rrc=   2.1081727302431230
        ix=13
        aa=ccc(ix)
        bb=ccc(ix+1)
        cc=ccc(ix+2)
        rrc=ccc(ix+3)
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a Ar-? interaction"
          stop
        endif
      else
          write(6,*)"Cant find a ?-? interaction"
          stop
      endif

      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rra=rr*autoang

! NOTE CANNOT HAVE BOTH TROYA FORM AND CUTOFF FORM

      if (troya) then   ! Troya uses different form & units
        v=aa*dexp(-rra*bb)+cc/rra**6
        v=v/autokcal
        dvdr = -aa*bb*dexp(-rra*bb)-6.d0*cc/rra**7
        dvdr=dvdr/autokcal*autoang
      elseif (cutoff) then  ! cutoff 1/R**-6 at short distances
        v=aa*dexp(-rra/bb)-(cc**6/(rra**6+rrc**6))
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)
     &      +6.d0*(cc**6)*(rra**5)/(rra**6+rrc**6)**2
        dvdr=dvdr/autocmi*autoang
      else
        v=aa*dexp(-rra/bb)-(cc/rra)**6
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)+(6.d0/rra)*(cc/rra)**6
        dvdr=dvdr/autocmi*autoang
      endif
        v1=v1+v

c      print *,m1,m2,rra,v*autocmi,v1*autocmi

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + dvdr*dx/rr
      dvdx(j) = dvdx(j) - dvdr*dx/rr
      dvdy(i) = dvdy(i) + dvdr*dy/rr
      dvdy(j) = dvdy(j) - dvdr*dy/rr
      dvdz(i) = dvdz(i) + dvdr*dz/rr
      dvdz(j) = dvdz(j) - dvdr*dz/rr

    2 continue
    1 continue
      v=v1
c      print *,'rgexp',v1
      return
      end

c **********************************************************************
c **********************************************************************
      subroutine tinkerpot(symb,xx,yy,zz,pema,dxx,dyy,dzz,nat,mnat)

      use sizes
      use atoms
      use files
      use inform
      use iounit

      implicit none
      double precision autokcal,autoang
      parameter(autokcal=627.509d0)
      parameter(autoang=0.52917706d0)

      integer iprepot,mnat,nat,i
      double precision xx(mnat), yy(mnat), zz(mnat),
     &                dxx(mnat),dyy(mnat),dzz(mnat),
     &                pema,energy,derivs(3,maxatm)
      character*2 symb(mnat)

      save iprepot
      entry prepot
      if (iprepot.eq.0) then
c     set up the structure and mechanics calculation
      call initial
      call getxyz
      call mechanic
      iprepot=1
      return
      endif

      do i=1,n
      x(i)=xx(i)*autoang
      y(i)=yy(i)*autoang
      z(i)=zz(i)*autoang
      enddo

      call gradient (energy,derivs)
      pema=energy/autokcal

      do i=1,nat
        dxx(i)=derivs(1,i)/autokcal*autoang
        dyy(i)=derivs(2,i)/autokcal*autoang
        dzz(i)=derivs(3,i)/autokcal*autoang
      enddo

      return

      end
c **********************************************************************
c **********************************************************************




