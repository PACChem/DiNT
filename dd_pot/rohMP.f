      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c NEW "AUTOFIT" PESs

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
      if ((symb(i).eq."N").or.
     &    (symb(i).eq."N"))  at(i)=2       ! carbon
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
      if ((symb(i).eq."N2").or.
     &    (symb(i).eq."n2")) at(i)=26      ! N2  ! label your bath atoms N2, not N
      if ((symb(i).eq."Ow").or.
     &    (symb(i).eq."ow")) at(i)=31      ! O in H2O 
      if ((symb(i).eq."Hw").or.
     &    (symb(i).eq."Hw")) at(i)=32      ! H in H2O
      if (at(i).eq.0) then ! atom not found
           write(6,*)"Atom # ",i," (",symb(i),") not found"
           stop
      endif
      enddo

      natom2=0
      natom3=0
      nh=1
      do i=1,natom
        if (at(i).le.20) then ! collect target atoms
          natom2=natom2+1
          x2(natom2)=x(i)
          y2(natom2)=y(i)
          z2(natom2)=z(i)
          symb2(natom2)=symb(i)
          at2(natom2)=at(i)
        elseif (at(i).le.30) then ! collect atomic and diatomic bath atoms
          natom3=natom3+1
          x3(natom3)=x(i)
          y3(natom3)=y(i)
          z3(natom3)=z(i)
          symb3(natom3)=symb(i)
          at3(natom3)=at(i)
        else ! collect larger bath atoms
          natom3=natom3+1
          if (at(i).eq.31) then
           x3(1)=x(i)
           y3(1)=y(i)
           z3(1)=z(i)
           symb3(1)=symb(i)
          elseif (at(i).eq.32) then
           nh=nh+1
           x3(nh)=x(i)
           y3(nh)=y(i)
           z3(nh)=z(i)
           symb3(nh)=symb(i)
          endif
        endif
      enddo
c CATCH SPECIAL CASES
      if (natom2.eq.3.and.at(1).eq.5.and.at(2).eq.1.and.at(3).eq.1.and.
     & natom3.eq.3.and.at(4).eq.31.and.at(5).eq.32.and.at(6).eq.33) then
c      h2o+h2o
       call water(x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
       call water(x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
       call lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
      else
c      print *,natom,natom2,natom3
      if (natom3.eq.2) 
     &   call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
      if (natom3.eq.3) 
     &   call water(x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
      if (natom2.ne.0) 
     &   call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
     &          natom2,maxatom)
      if (natom3.ne.0.and.natom2.ne.0) 
     &   call lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
      endif

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


! ONE DIMENSIONAL DIATOMIC BATHS
      subroutine bath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      integer at(maxatom)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      parameter(autoev=27.2113961d0)

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
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

        if (at(1).eq.25.and.at(2).eq.25) then
! H2 bath
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
        v=v/autoev

        dvdr=((c1+c2*rr)*(-c3)+c2)*dexp(-c3*rr)
     &      +((c4+c5*rr)*(-c6)+c5)*dexp(-c6*rr)*rr**2
     &       +(c4+c5*rr)*dexp(-c6*rr)*rr*2.d0
        dvdr=dvdr/autoev

        elseif (at(1).eq.29.and.at(2).eq.29) then
! O2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       Jasper April 3, 2012
        de=42046.5d0 ! exp De in cm-1 
        re=1.2075d0 ! exp in A
        c1= 2.6938139d0 ! my fit
        c2= 0.384763939d0
        c3= 0.812506485d0

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif (at(1).eq.26.and.at(2).eq.26) then
! N2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       agrees reasonably well with more complicated form of LeRoy (JCP
!       125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, 
!                                     PRL 53, 159 (1984), as quoted by
!                                     LeRoy)
        re=1.097679d0 ! exp in A
        c1=2.68872341 ! my fit
        c2=0.240070803
        c3=0.472261727

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif ((at(1).eq.27.and.at(2).eq.28).or.
     &          (at(1).eq.28.and.at(2).eq.27)) then
! CO bath
!       Morse. Fit to RKR data of PAUL H. KRUPENIE and STANLEY WEISSMAN,
!       J. CHEM. PHYS. 43, 1529 (1965)
!       with De = 11.06 eV
        de=11.06d0    ! exp De in eV
        de=de/autoev
        re=1.128322d0 ! exp in A
        re=re/autoang
        beta=1.d0/0.428d0  ! my fit in 1/A
        beta=beta*autoang

        yy=rr-re
        v = de*(1.d0-dexp(-beta*yy))**2
        dvdr=2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)*beta

!       elseif (at(1).eq.??.and.at(2).eq.??) then
! OTHER DIATOMIC BATHS HERE
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

      integer nfitparams,mfitparams,ix
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams

      save/amfit/

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

      if (m2.eq.m2) then  ! ANYTHING
        ix=m1
        if (ix.gt.2) ix=ix-1
        aa = ccc((ix-1)*4+1)
        bb = ccc((ix-1)*4+2)
        cc = ccc((ix-1)*4+3)
        rrc = ccc((ix-1)*4+4)
        aa=(10.d0**aa)
        cutoff=.true.
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
      character*2 symb(mnat),dum

      integer nfitparams,mfitparams
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams

      save iprepot
      entry prepot
      if (iprepot.eq.0) then
      call prepot2
      call prepot3
c     set up the structure and mechanics calculation
      call initial
      call getxyz
      call getkey
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




********************************************

      subroutine prepot2

      implicit double precision(a-h,o-z)
      include 'paramls.inc'
      dimension iagroup(maxatom),
     &  ind(maxterm,maxpair),nfound(maxatom),
     &  iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair)
      character*2 symb(maxatom),dum
      logical lreadbasis
 
      common/foox/rrr,nncoef,natom1

      save npairs,nterms,ind,ibasis

ccc GENERATE BASIS ccc
      if (.false.) then
ccc GENERATE BASIS ccc


c     generate atom permutation lists
      do i=1,natom
      nngroup(i)=0
      enddo
      do i=1,natom
      if (iagroup(i).gt.ngroup) ngroup=iagroup(i)
      nngroup(iagroup(i))=nngroup(iagroup(i))+1
      enddo

      nn=0

      do i=1,ngroup

      n=0
      do k=1,natom
      if (iagroup(k).eq.i) then
      n=n+1
      idum(n)=k
      endif
      enddo
      
      npermute=0
      call heapp(idum,n,n,idum2,npermute)
      nperm0(i)=nn+1
      nperm1(i)=nn+npermute
      do k=1,npermute
      nn=nn+1
      m=0
      do j=1,natom
      idum3(nn,j)=0
      if (iagroup(j).eq.i) then
          m=m+1
          idum3(nn,j)=idum2(k,m)
      endif
      enddo
      enddo

      enddo

      ntmp=1
      do i=1,ngroup
      idum(i)=nperm0(i)
      print *,"Group ",i," has ",(nperm1(i)-nperm0(i)+1)," permutations"
      ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo
      print *,"For a total of ",ntmp," permutations"

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
        print *,"npermute (",npermute,") > maxperm (",maxperm,")"
        print *,"NOTE: maxperm needs to be at least npermute + 1"
        stop
        endif

        do i=1,natom
        iatom(npermute,i)=0
        do j=1,ngroup
        iatom(npermute,i)=iatom(npermute,i)+idum3(idum(j),i)
        enddo
        enddo

        idum(ngroup)=idum(ngroup)+1
 777    continue

        do i=1,ngroup
        if (idum(i).gt.nperm1(i)) then
        if (i.eq.1) go to 778
        idum(i)=nperm0(i)
        idum(i-1)=idum(i-1)+1
        go to 777
        endif
        enddo 

      enddo
 778  continue

      print *
      print *,'Atom permutations',npermute
      do i=1,min(npermute,100)
      print *,i,":",(iatom(i,j),j=1,natom)
      enddo

      ii=0
      do i=1,natom1
      do j=natom1+1,natom
      ii=ii+1
      index(i,j)=ii
      enddo
      enddo

      write(6,*)
      write(6,*)"Pair permutations"
      write(6,'(22x,100(a3,"- ",a3,4x))')
     &   ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1) 
      write(6,'(21x,100(i3," -",i3,4x))')((i,j,j=1+natom1,natom),
     &   i=1,natom1) 
      do ii=1,npermute
      iix=0
      do i=1,natom1
      do j=natom1+1,natom
      iix=iix+1
      ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
      enddo
      enddo
      if (ii.le.100) print *,ii,":",(ix(ii,iix),iix=1,npairs)
      enddo

c generate terms using individual power constraints
      ii=1
      do i=1,npairs
      ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
      print *,"number of terms (",ii,") > maxterm (",maxterm,")"
        stop
        endif

        do i=1,npairs
        ind(ii,i)=ind(ii-1,i)
        enddo
        ind(ii,npairs)=ind(ii,npairs)+1
 300    continue
        indtot=0
        do i=1,npairs
        indtot=indtot+ind(ii,i)
        if (ind(ii,i).gt.ipow.or.indtot.gt.ipowt) then ! ipow(i) would allow atom-atom-type-dependent limits
        if (i.eq.1) go to 400
        ind(ii,i)=0
        ind(ii,i-1)=ind(ii,i-1)+1
        go to 300
        endif
        enddo
      enddo
 400  continue
      nterms=ii-1

      print *
      print *,"Basis # (Group):  Powers"

c symmetrize
      nbasis=0
      DO ii=1,nterms
      ifail=0
      do i=1,ii-1
        do j=1,npermute
          ifail=1
          do k=1,npairs
            if (ind(i,k).ne.ind(ii,ix(j,k))) ifail=0
          enddo
          if (ifail.eq.1) go to 1010
        enddo
      enddo
 1010 continue

      if (ifail.eq.0) then
      nbasis=nbasis+1
      ibasis(ii)=nbasis
      else
      ibasis(ii)=ibasis(i)
      endif
      write(6,'(i5,"  (",i5,"):",100i8)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      ENDDO

      nncoef=nbasis
      print *,'nncoef = ',nncoef
 
      open(55,file="basis.dat")
      write(55,*)natom1,npairs,nncoef,nterms,
     & " ! atom pairs, coefficients, terms"
      write(55,*)" TERM GROUP :     EXPONENTS"
      do ii=1,nterms
      write(55,'(2i6," : ",1000i5)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)

ccc READ BASIS ccc
      else
      open(55,file="basis.dat")
      read(55,*)natom1,npairs,nncoef,nterms
      read(55,*)
      do i=1,nterms
      read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
      enddo
      close(55)

      endif

      return

      entry funcs1(iii,basis,ncoef)

c print *,ncoef,npairs,nterms

      do j=1,ncoef
      basis(j)=0.d0
      enddo

      do j=1,npairs
      r(j)=dexp(-rrr(iii,j)*autoang)
      enddo

      do i=1,nterms
      arg=1.d0
      do j=1,npairs
      arg=arg*(r(j)**ind(i,j))
      enddo
      basis(ibasis(i))=basis(ibasis(i))+arg
      enddo

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'paramls.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
         ii=ii+1
         do i=1,n
           iia(ii,i)=ia(i)
         enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif

      enddo

      end subroutine

***************************************************
      subroutine lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatomx)
c
      implicit double precision (a-h,o-z)
c
      include 'paramls.inc'
      dimension coef(maxterm),sig(maxdata)
      dimension basis(maxterm)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension vv2(maxdata),rrr2(maxdata,maxpair)
      dimension rcom2(maxdata),xprint(50,20)
      dimension x(maxatom),y(maxatom),z(maxatom)

      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      integer at(maxatom)


      character*2 dum

      common/foox/rrr,nncoef,natom1

      save iprepot3
      entry prepot3
      if (iprepot3.eq.0) then
      open(77,file="coef.dat")
      do k=1,nncoef
      read(77,*)i,coef(k)
      enddo
      iprepot3=1
      return
      endif

c      print *
c      do j=1,natom
c      print *,at(j),x(j)*autoang,y(j)*autoang,z(j)*autoang
c      enddo

      ii=0
      do j=1,natom1  ! molecule
      do k=natom1+1,natom  ! bath
      ii=ii+1
      rrr(1,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
      enddo  
      enddo  


      ncoef=nncoef
      call funcs1(1,basis,ncoef) 
      v=0.d0
      do j=1,ncoef
         v=v+coef(j)*basis(j)
      enddo
c      write(6,*)v,(rrr(1,j)*autoang,j=1,15)
      v=v/autocmi

      resp=0.0001
      ii=0
      do j=1,natom1  ! molecule
      do k=natom1+1,natom  ! bath
      dx=x(j)-x(k)
      dy=y(j)-y(k)
      dz=z(j)-z(k)
      ii=ii+1

      rrr(1,ii)=rrr(1,ii)+resp
      call funcs1(1,basis,ncoef) 
      vp=0.d0
      do l=1,ncoef
         vp=vp+coef(l)*basis(l)
      enddo
      vp=vp/autocmi

      rrr(1,ii)=rrr(1,ii)-2.d0*resp
      call funcs1(1,basis,ncoef) 
      vm=0.d0
      do l=1,ncoef
         vm=vm+coef(l)*basis(l)
      enddo
      vm=vm/autocmi
      rrr(1,ii)=rrr(1,ii)+resp

      dtmpdrr=(vp-vm)/(2.d0*resp)
      dvdx(j) = dvdx(j) + dtmpdrr*dx/rrr(1,ii)
      dvdx(k) = dvdx(k) - dtmpdrr*dx/rrr(1,ii)
      dvdy(j) = dvdy(j) + dtmpdrr*dy/rrr(1,ii)
      dvdy(k) = dvdy(k) - dtmpdrr*dy/rrr(1,ii)
      dvdz(j) = dvdz(j) + dtmpdrr*dz/rrr(1,ii)
      dvdz(k) = dvdz(k) - dtmpdrr*dz/rrr(1,ii)
      enddo
      enddo

      return
 
      end

****************************************************************
C****************************************************************
      subroutine vibpot(AJa,AJb,AJc,vvv,n)
      implicit real*8 (a-h,o-z)
c
c     pes for h2o,
c     Harry Partridge and David W. Schwenke, J. Chem. Phys.,
c     submitted Nov. 8, 1996.
c     rij(i,1)& rij(i,2) are oh distances in au
c     rij(i,3) is hoh angle in rad
c     v(i) is pes in au
c     n is number of geometries
c     mass dependent factors are included. the nuclear masses
c     should be passed to this program using the array xm in
c     common potmcm. xm(1) is the
c     mass of the hydrogen associated with rij(i,1), and xm(2)
c     is the mass of the hydrogen associated with rij(i,2).
c     all masses are in au.
c
      dimension rij(n,3),v(n),c5z(245),cbasis(245),ccore(245),
     $          crest(245),idx(245,3),fmat(15,3),cmass(9),idxm(9,3)
c      common/potrot/fact1,fact2,c1,s1,icoord,xm(2),xmx,iperm 7d27s90
c      common/potmcm/xm(2)
c
c     expansion indicies
c
       data (idx(i,1),i=1,245)/
     $ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
     $ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $ 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
     $ 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
     $ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     $ 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
     $ 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
     $ 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     $ 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5,
     $ 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,
     $ 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6,
     $ 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     $ 9, 9, 9, 9, 9/
       data (idx(i,2),i=1,245)/
     $ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $ 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $ 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
     $ 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,
     $ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $ 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
     $ 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
     $ 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,
     $ 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,
     $ 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,
     $ 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,
     $ 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,
     $ 1, 1, 1, 1, 1/
       data (idx(i,3),i=1,245)/
     $ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,
     $ 6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
     $12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,
     $ 6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,
     $ 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     $11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,
     $ 9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,
     $ 9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,
     $ 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,
     $ 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,
     $ 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,
     $ 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,
     $ 3, 4, 5, 6, 7/
c
c     expansion coefficients for 5z ab initio data
c
       data (c5z(i),i=1,245)/
     $ 4.2278462684916D+04, 4.5859382909906D-02, 9.4804986183058D+03,
     $ 7.5485566680955D+02, 1.9865052511496D+03, 4.3768071560862D+02,
     $ 1.4466054104131D+03, 1.3591924557890D+02,-1.4299027252645D+03,
     $ 6.6966329416373D+02, 3.8065088734195D+03,-5.0582552618154D+02,
     $-3.2067534385604D+03, 6.9673382568135D+02, 1.6789085874578D+03,
     $-3.5387509130093D+03,-1.2902326455736D+04,-6.4271125232353D+03,
     $-6.9346876863641D+03,-4.9765266152649D+02,-3.4380943579627D+03,
     $ 3.9925274973255D+03,-1.2703668547457D+04,-1.5831591056092D+04,
     $ 2.9431777405339D+04, 2.5071411925779D+04,-4.8518811956397D+04,
     $-1.4430705306580D+04, 2.5844109323395D+04,-2.3371683301770D+03,
     $ 1.2333872678202D+04, 6.6525207018832D+03,-2.0884209672231D+03,
     $-6.3008463062877D+03, 4.2548148298119D+04, 2.1561445953347D+04,
     $-1.5517277060400D+05, 2.9277086555691D+04, 2.6154026873478D+05,
     $-1.3093666159230D+05,-1.6260425387088D+05, 1.2311652217133D+05,
     $-5.1764697159603D+04, 2.5287599662992D+03, 3.0114701659513D+04,
     $-2.0580084492150D+03, 3.3617940269402D+04, 1.3503379582016D+04,
     $-1.0401149481887D+05,-6.3248258344140D+04, 2.4576697811922D+05,
     $ 8.9685253338525D+04,-2.3910076031416D+05,-6.5265145723160D+04,
     $ 8.9184290973880D+04,-8.0850272976101D+03,-3.1054961140464D+04,
     $-1.3684354599285D+04, 9.3754012976495D+03,-7.4676475789329D+04,
     $-1.8122270942076D+05, 2.6987309391410D+05, 4.0582251904706D+05,
     $-4.7103517814752D+05,-3.6115503974010D+05, 3.2284775325099D+05,
     $ 1.3264691929787D+04, 1.8025253924335D+05,-1.2235925565102D+04,
     $-9.1363898120735D+03,-4.1294242946858D+04,-3.4995730900098D+04,
     $ 3.1769893347165D+05, 2.8395605362570D+05,-1.0784536354219D+06,
     $-5.9451106980882D+05, 1.5215430060937D+06, 4.5943167339298D+05,
     $-7.9957883936866D+05,-9.2432840622294D+04, 5.5825423140341D+03,
     $ 3.0673594098716D+03, 8.7439532014842D+04, 1.9113438435651D+05,
     $-3.4306742659939D+05,-3.0711488132651D+05, 6.2118702580693D+05,
     $-1.5805976377422D+04,-4.2038045404190D+05, 3.4847108834282D+05,
     $-1.3486811106770D+04, 3.1256632170871D+04, 5.3344700235019D+03,
     $ 2.6384242145376D+04, 1.2917121516510D+05,-1.3160848301195D+05,
     $-4.5853998051192D+05, 3.5760105069089D+05, 6.4570143281747D+05,
     $-3.6980075904167D+05,-3.2941029518332D+05,-3.5042507366553D+05,
     $ 2.1513919629391D+03, 6.3403845616538D+04, 6.2152822008047D+04,
     $-4.8805335375295D+05,-6.3261951398766D+05, 1.8433340786742D+06,
     $ 1.4650263449690D+06,-2.9204939728308D+06,-1.1011338105757D+06,
     $ 1.7270664922758D+06, 3.4925947462024D+05,-1.9526251371308D+04,
     $-3.2271030511683D+04,-3.7601575719875D+05, 1.8295007005531D+05,
     $ 1.5005699079799D+06,-1.2350076538617D+06,-1.8221938812193D+06,
     $ 1.5438780841786D+06,-3.2729150692367D+03, 1.0546285883943D+04,
     $-4.7118461673723D+04,-1.1458551385925D+05, 2.7704588008958D+05,
     $ 7.4145816862032D+05,-6.6864945408289D+05,-1.6992324545166D+06,
     $ 6.7487333473248D+05, 1.4361670430046D+06,-2.0837555267331D+05,
     $ 4.7678355561019D+05,-1.5194821786066D+04,-1.1987249931134D+05,
     $ 1.3007675671713D+05, 9.6641544907323D+05,-5.3379849922258D+05,
     $-2.4303858824867D+06, 1.5261649025605D+06, 2.0186755858342D+06,
     $-1.6429544469130D+06,-1.7921520714752D+04, 1.4125624734639D+04,
     $-2.5345006031695D+04, 1.7853375909076D+05,-5.4318156343922D+04,
     $-3.6889685715963D+05, 4.2449670705837D+05, 3.5020329799394D+05,
     $ 9.3825886484788D+03,-8.0012127425648D+05, 9.8554789856472D+04,
     $ 4.9210554266522D+05,-6.4038493953446D+05,-2.8398085766046D+06,
     $ 2.1390360019254D+06, 6.3452935017176D+06,-2.3677386290925D+06,
     $-3.9697874352050D+06,-1.9490691547041D+04, 4.4213579019433D+04,
     $ 1.6113884156437D+05,-7.1247665213713D+05,-1.1808376404616D+06,
     $ 3.0815171952564D+06, 1.3519809705593D+06,-3.4457898745450D+06,
     $ 2.0705775494050D+05,-4.3778169926622D+05, 8.7041260169714D+03,
     $ 1.8982512628535D+05,-2.9708215504578D+05,-8.8213012222074D+05,
     $ 8.6031109049755D+05, 1.0968800857081D+06,-1.0114716732602D+06,
     $ 1.9367263614108D+05, 2.8678295007137D+05,-9.4347729862989D+04,
     $ 4.4154039394108D+04, 5.3686756196439D+05, 1.7254041770855D+05,
     $-2.5310674462399D+06,-2.0381171865455D+06, 3.3780796258176D+06,
     $ 7.8836220768478D+05,-1.5307728782887D+05,-3.7573362053757D+05,
     $ 1.0124501604626D+06, 2.0929686545723D+06,-5.7305706586465D+06,
     $-2.6200352535413D+06, 7.1543745536691D+06,-1.9733601879064D+04,
     $ 8.5273008477607D+04, 6.1062454495045D+04,-2.2642508675984D+05,
     $ 2.4581653864150D+05,-9.0376851105383D+05,-4.4367930945690D+05,
     $ 1.5740351463593D+06, 2.4563041445249D+05,-3.4697646046367D+03,
     $-2.1391370322552D+05, 4.2358948404842D+05, 5.6270081955003D+05,
     $-8.5007851251980D+05,-6.1182429537130D+05, 5.6690751824341D+05,
     $-3.5617502919487D+05,-8.1875263381402D+02,-2.4506258140060D+05,
     $ 2.5830513731509D+05, 6.0646114465433D+05,-6.9676584616955D+05,
     $ 5.1937406389690D+05, 1.7261913546007D+05,-1.7405787307472D+04,
     $-3.8301842660567D+05, 5.4227693205154D+05, 2.5442083515211D+06,
     $-1.1837755702370D+06,-1.9381959088092D+06,-4.0642141553575D+05,
     $ 1.1840693827934D+04,-1.5334500255967D+05, 4.9098619510989D+05,
     $ 6.1688992640977D+05, 2.2351144690009D+05,-1.8550462739570D+06,
     $ 9.6815110649918D+03,-8.1526584681055D+04,-8.0810433155289D+04,
     $ 3.4520506615177D+05, 2.5509863381419D+05,-1.3331224992157D+05,
     $-4.3119301071653D+05,-5.9818343115856D+04, 1.7863692414573D+03,
     $ 8.9440694919836D+04,-2.5558967650731D+05,-2.2130423988459D+04,
     $ 4.4973674518316D+05,-2.2094939343618D+05/
c
c     expansion coefficients for basis correction
c
       data (cbasis(i),i=1,245)/
     $ 6.9770019624764D-04,-2.4209870001642D+01, 1.8113927151562D+01,
     $ 3.5107416275981D+01,-5.4600021126735D+00,-4.8731149608386D+01,
     $ 3.6007189184766D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-7.7178474355102D+01,-3.8460795013977D+01,-4.6622480912340D+01,
     $ 5.5684951167513D+01, 1.2274939911242D+02,-1.4325154752086D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00,-6.0800589055949D+00,
     $ 8.6171499453475D+01,-8.4066835441327D+01,-5.8228085624620D+01,
     $ 2.0237393793875D+02, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 3.3525582670313D+02, 7.0056962392208D+01,-4.5312502936708D+01,
     $-3.0441141194247D+02, 2.8111438108965D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-1.2983583774779D+02, 3.9781671212935D+01,
     $-6.6793945229609D+01,-1.9259805675433D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-8.2855757669957D+02,-5.7003072730941D+01,
     $-3.5604806670066D+01, 9.6277766002709D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 8.8645622149112D+02,-7.6908409772041D+01,
     $ 6.8111763314154D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 2.5090493428062D+02,-2.3622141780572D+02, 5.8155647658455D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 2.8919570295095D+03,
     $-1.7871014635921D+02,-1.3515667622500D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-3.6965613754734D+03, 2.1148158286617D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00,-1.4795670139431D+03,
     $ 3.6210798138768D+02, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-5.3552886800881D+03, 3.1006384016202D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 1.6241824368764D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 4.3764909606382D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 1.0940849243716D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 3.0743267832931D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00/
c
c     expansion coefficients for core correction
c
       data (ccore(i),i=1,245)/
     $ 2.4332191647159D-02,-2.9749090113656D+01, 1.8638980892831D+01,
     $-6.1272361746520D+00, 2.1567487597605D+00,-1.5552044084945D+01,
     $ 8.9752150543954D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-3.5693557878741D+02,-3.0398393196894D+00,-6.5936553294576D+00,
     $ 1.6056619388911D+01, 7.8061422868204D+01,-8.6270891686359D+01,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00,-3.1688002530217D+01,
     $ 3.7586725583944D+01,-3.2725765966657D+01,-5.6458213299259D+00,
     $ 2.1502613314595D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 5.2789943583277D+02,-4.2461079404962D+00,-2.4937638543122D+01,
     $-1.1963809321312D+02, 2.0240663228078D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-6.2574211352272D+02,-6.9617539465382D+00,
     $-5.9440243471241D+01, 1.4944220180218D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-1.2851139918332D+03,-6.5043516710835D+00,
     $ 4.0410829440249D+01,-6.7162452402027D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 1.0031942127832D+03, 7.6137226541944D+01,
     $-2.7279242226902D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-3.3059000871075D+01, 2.4384498749480D+01,-1.4597931874215D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 1.6559579606045D+03,
     $ 1.5038996611400D+02,-7.3865347730818D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-1.9738401290808D+03,-1.4149993809415D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00,-1.2756627454888D+02,
     $ 4.1487702227579D+01, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-1.7406770966429D+03,-9.3812204399266D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-1.1890301282216D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 2.3723447727360D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-1.0279968223292D+03, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 5.7153838472603D+02, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00/
c
c     expansion coefficients for v rest
c
       data (crest(i),i=1,245)/
     $ 0.0000000000000D+00,-4.7430930170000D+00,-1.4422132560000D+01,
     $-1.8061146510000D+01, 7.5186735000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $-2.7962099800000D+02, 1.7616414260000D+01,-9.9741392630000D+01,
     $ 7.1402447000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00,-7.8571336480000D+01,
     $ 5.2434353250000D+01, 7.7696745000000D+01, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 1.7799123760000D+02, 1.4564532380000D+02, 2.2347226000000D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-4.3823284100000D+02,-7.2846553000000D+02,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00,-2.6752313750000D+02, 3.6170310000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,
     $ 0.0000000000000D+00, 0.0000000000000D+00/
c
c     expansion indicies for mass correction
c
       data idxm/1,2,1,1,3,2,1,2,1,
     $           2,1,1,3,1,2,2,1,1,
     $           1,1,2,1,1,1,2,2,3/
c
c     expansion coefficients for mass correction
c
       data cmass/ -8.3554183D+00,3.7036552D+01,-5.2722136D+00,
     $      1.6843857D+01,-7.0929741D+01,5.5380337D+00,-2.9962997D+01,
     $      1.3637682D+02,-3.0530195d+00/
c
c     two body parameters
c
       data reoh,thetae,b1,roh,alphaoh,deoh,phh1,phh2/0.958649d0,
     $      104.3475d0,2.0d0,0.9519607159623009d0,2.587949757553683d0,
     $      42290.92019288289d0,16.94879431193463d0,12.66426998162947d0/
c
c     scaling factors for contributions to emperical potential
c
       data f5z,fbasis,fcore,frest/0.99967788500000d0,
     $      0.15860145369897d0,-1.6351695982132d0,1d0/
      save
      data ifirst/0/

c AJ HACK
      rij(1,1)=AJa
      rij(1,2)=AJb
      rij(1,3)=AJc

      if(ifirst.eq.0)then
       ifirst=1
c      write(6,1)
    1  format(/1x,'pes for h2o',
     $        /1x,'by Harry Partridge and David W. Schwenke',
     $        /1x,'submitted to J. Chem. Phys. Nov. 8, 1996')
c      write(6,56)
   56  format(/1x,'parameters before adjustment')
c      write(6,55)phh1,phh2,deoh,alphaoh,roh
   55  format(/1x,'two body potential parameters:',
     $        /1x,'hh: phh1 = ',f10.1,' phh2 = ',f5.2,
     $        /1x,'oh: deoh = ',f10.1,' alpha = ',f7.4,
     $        ' re = ',f7.4)
c      write(6,4)reoh,thetae,b1
    4  format(/1x,'three body parameters:',
     $        /1x,'reoh = ',f10.4,' thetae = ',f10.4,
     $        /1x,'betaoh = ',f10.4,
     $        /1x,'    i    j    k',7x,'c5z',9x,'cbasis',10x,'ccore',
     $        10x,'crest')
       do 2 i=1,245
c       write(6,5)(idx(i,j)-1,j=1,3),c5z(i),cbasis(i),ccore(i),crest(i)
    5   format(1x,3i5,1p4e15.7)
    2  continue
c
c     remove mass correction from vrest
c
c       xmh=1836.152697d0
       xmh=1.0078250d0*1822.88853d0
       xmhi=1d0/xmh
       xmd=3670.483031d0
       fact=1d0/((1d0/xmd)-(1d0/xmh))
c      write(6,65)
   65  format(/1x,'parameters for delta v hdo ',
     $       /1x,'    i    j    k')
       do 60 i=1,9
c       write(6,5)(idxm(i,j)-1,j=1,3),cmass(i)
        cmass(i)=cmass(i)*fact
        corr=cmass(i)*xmhi
        if(idxm(i,1).eq.idxm(i,2))corr=corr*0.5d0
        do 61 j=1,245
         if(idx(j,1).eq.idxm(i,1).and.idx(j,2).eq.idxm(i,2).and.
     $      idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 62
         end if
   61   continue
   62   continue
        do 63 j=1,245
         if(idx(j,2).eq.idxm(i,1).and.idx(j,1).eq.idxm(i,2).and.
     $      idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 64
         end if
   63   continue
   64   continue
   60  continue
c      write(6,70)xm
   70  format(/1x,'masses used for mass correction: ',1p2e15.7)
       xm1=1d0/xmh
       xm2=1d0/xmh
c
c     adjust parameters using scale factors
c
c      write(6,57)f5z,fbasis,fcore,frest
   57  format(/1x,'adjusting parameters using scale factors ',
     $        /1x,'f5z =    ',f11.8,
     $        /1x,'fbasis = ',f11.8,
     $        /1x,'fcore =  ',f11.8,
     $        /1x,'frest =  ',f11.8)
       phh1=phh1*f5z
       deoh=deoh*f5z
       do 59 i=1,245
        c5z(i)=f5z*c5z(i)+fbasis*cbasis(i)+fcore*ccore(i)
     $       +frest*crest(i)
   59  continue
c      write(6,55)phh1,phh2,deoh,alphaoh,roh
c      write(6,58)reoh,thetae,b1,((idx(i,j)-1,j=1,3),c5z(i),i=1,245)
   58  format(/1x,'three body parameters:',
     $        /1x,'reoh = ',f10.4,' thetae = ',f10.4,
     $        /1x,'betaoh = ',f10.4,
     $        /1x,'    i    j    k   cijk',
     $        /(1x,3i5,1pe15.7))
       do 66 i=1,9
        cmass(i)=cmass(i)*frest
   66  continue
c      write(6,76)((idxm(i,j),j=1,3),cmass(i),i=1,9)
   76  format(/1x,'mass correction factors ',
     $        /1x,'    i    j    k   cijk',
     $        /(1x,3i5,1pe15.7))
c
c     convert parameters from 1/cm, angstrom to a.u.
c
       reoh=reoh/0.529177249d0
       b1=b1*0.529177249d0*0.529177249d0
       do 3 i=1,245
        c5z(i)=c5z(i)*4.556335d-6 
    3  continue
       do 67 i=1,9
        cmass(i)=cmass(i)*4.556335d-6
   67  continue
       rad=acos(-1d0)/1.8d2
       ce=cos(thetae*rad)
       phh1=phh1*exp(phh2)
       phh1=phh1*4.556335d-6
       phh2=phh2*0.529177249d0
       deoh=deoh*4.556335d-6
       roh=roh/0.529177249d0
       alphaoh=alphaoh*0.529177249d0
       c5z(1)=c5z(1)*2d0
      end if
      do 6 i=1,n
       x1=(rij(i,1)-reoh)/reoh
       x2=(rij(i,2)-reoh)/reoh
       x3=cos(rij(i,3))-ce
       rhh=sqrt(rij(i,1)**2+rij(i,2)**2
     $      -2d0*rij(i,1)*rij(i,2)*cos(rij(i,3)))
       vhh=phh1*exp(-phh2*rhh)
       ex=exp(-alphaoh*(rij(i,1)-roh))
       voh1=deoh*ex*(ex-2d0)
       ex=exp(-alphaoh*(rij(i,2)-roh))
       voh2=deoh*ex*(ex-2d0)
       fmat(1,1)=1d0
       fmat(1,2)=1d0
       fmat(1,3)=1d0
       do 10 j=2,15
        fmat(j,1)=fmat(j-1,1)*x1
        fmat(j,2)=fmat(j-1,2)*x2
        fmat(j,3)=fmat(j-1,3)*x3
   10  continue
       v(i)=0d0
       do 12 j=2,245
        term=c5z(j)*(fmat(idx(j,1),1)*fmat(idx(j,2),2)
     $                    +fmat(idx(j,2),1)*fmat(idx(j,1),2))
     $                    *fmat(idx(j,3),3)
        v(i)=v(i)+term
   12  continue
       v1=0d0
       v2=0d0
       do 13 j=1,9
        v1=v1+cmass(j)*fmat(idxm(j,1),1)*fmat(idxm(j,2),2)
     $       *fmat(idxm(j,3),3)
        v2=v2+cmass(j)*fmat(idxm(j,2),1)*fmat(idxm(j,1),2)
     $       *fmat(idxm(j,3),3)
   13  continue
       v(i)=v(i)+xm1*v1+xm2*v2
       v(i)=v(i)*exp(-b1*((rij(i,1)-reoh)**2+(rij(i,2)-reoh)**2))
     $      +c5z(1)
     $      +voh1+voh2+vhh
    6 continue
      vvv=v(1)
      return
      end



      subroutine water(x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)

      dimension r(3)

      r(1)=dsqrt((x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2)
      r(2)=dsqrt((x(1)-x(3))**2+(y(1)-y(3))**2+(z(1)-z(3))**2)
      r(3)=dsqrt((x(2)-x(3))**2+(y(2)-y(3))**2+(z(2)-z(3))**2)
      angle=-(r(3)**2-r(1)**2-r(2)**2)/r(1)/r(2)/2.d0
      angle=dacos(min(1.d0,max(angle,-1.d0)))

      call vibpot(r(1),r(2),angle,v,1)
      v=v+0.000052/27.211d0

c      print *,'xyz0',(r(k)*autoang,k=1,3),angle/3.14158*180.,v

      resp=0.000001d0

      do ij=1,3
      if (ij.eq.1) i=1 
      if (ij.eq.1) j=2 
      if (ij.eq.2) i=1 
      if (ij.eq.2) j=3 
      if (ij.eq.3) i=2 
      if (ij.eq.3) j=3 
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(ij)=r(ij)+resp
      angle=-(r(3)**2-r(1)**2-r(2)**2)/r(1)/r(2)/2.d0
      angle=dacos(min(1.d0,max(angle,-1.d0)))
c      print *,'xyz1',(r(k)*autoang,k=1,3),angle/3.14158*180.
      call vibpot(r(1),r(2),angle,vp,1)
      r(ij)=r(ij)-2.d0*resp
      angle=-(r(3)**2-r(1)**2-r(2)**2)/r(1)/r(2)/2.d0
      angle=dacos(min(1.d0,max(angle,-1.d0)))
c      print *,'xyz2',(r(k)*autoang,k=1,3),angle/3.14158*180.
      call vibpot(r(1),r(2),angle,vm,1)
      r(ij)=r(ij)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)
      dvdx(i) = dvdx(i) + dtmpdrr*dx/r(ij)
      dvdx(j) = dvdx(j) - dtmpdrr*dx/r(ij)
      dvdy(i) = dvdy(i) + dtmpdrr*dy/r(ij)
      dvdy(j) = dvdy(j) - dtmpdrr*dy/r(ij)
      dvdz(i) = dvdz(i) + dtmpdrr*dz/r(ij)
      dvdz(j) = dvdz(j) - dtmpdrr*dz/r(ij)
      enddo

      return
      end

