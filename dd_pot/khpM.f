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
     &   call aipot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
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
        aa=   6.4374217963194678     
        bb=  0.28550370799890135     
        cc=   7.3473311563463239     
        rrc=   2.8850978118228703     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
        aa=   7.5285195471053195     
        bb=  0.26272774437696461     
        cc=   6.7756889553514208     
        rrc=   3.1267433698538163     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.4) then !    O
        aa=   8.9658192693868841     
        bb=  0.19998474074526201     
        cc=   7.5914487136448257     
        rrc=   2.7619861445966980     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.5) then !    OH H
        aa=   5.9999084444715720     
        bb=  0.28738670003357036     
        cc=   7.5051118503372294     
        rrc=   2.9338053529465622     
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
c     DD: A common interface for direct dynamics calls to the Gaussian & 
c     Molpro packages. Jasper, Dec 2007
c **********************************************************************
c **********************************************************************

c      subroutine pot(symb,x,y,z,pema,gpema,dvec,nat,mnat,nsurf,mnsurf)
      subroutine aipot(symb,x,y,z,pema,dx,dy,dz,nat,mnat)

      implicit none

c one state hack
      integer nsurf,mnsurf
      parameter(nsurf=1)
      parameter(mnsurf=1)

c in/out
      integer nat,mnat
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     & gpema(3,mnat,mnsurf),pema(mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf)
      double precision dz(mnat),dy(mnat),dx(mnat)

c local
      integer mnc,ifail
      parameter(mnc=10) ! max number of QC calls per geometry
      integer ic,nc,ns(mnc),np(mnc),ntmp,ictot,i,j,k,l,ii,jj
      double precision gptmp(3,mnat,mnsurf),ptmp(mnsurf)
      double precision dptmp(3,mnat,mnsurf,mnsurf)
      character*8 tname(mnc),tnam

c -------------------
c SET THESE
c     NC = number of separate QC calls per geom
c     NS(i) = number of states for call #i
c     NP(i) = QC package to be used for call #i
c           = 1 for G03
c           = 2 for Molpro 2006
      nc    = 1
      ns(1) = 1
      np(1) = 1
      tname(1) = "qc.g09"
c -------------------

      if (nc.gt.mnc) then
          print *,"nc (",nc,") is bigger than mnc (",mnc,")" 
          stop
      endif

c zero
      do i=1,mnsurf
        pema(i)=0.d0
        do k=1,3
        do l=1,mnat
          gpema(k,l,i)=0.d0
          do j=1,mnsurf
          dvec(k,l,i,j)=0.d0
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
c       call G03
        call dd_g03(symb,tnam,x,y,z,ptmp,gptmp,dptmp,nat,
     &               mnat,ntmp,mnsurf,ifail)
        elseif (np(ic).eq.2) then
c       call Molpro 2006
        call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,dptmp,nat,
     &               mnat,ntmp,mnsurf,ifail)
c         handle failures
          if (ifail.eq.0) then
c           everything OK
          else
c           can't fix
            write(6,*)"QC failure"
            stop
          endif
        else
          print *,"np(",ic,") = ",np(ic),", which is not allowed"
          stop
        endif
        do i=1,ntmp
          ii=i+ictot
          pema(ii)=ptmp(i)
          do k=1,3
          do l=1,nat
            gpema(k,l,ii)=gptmp(k,l,i)
c                print *,"x",k,l,ii,gpema(k,l,ii)
          enddo
          enddo
        enddo
        do i=1,ntmp
        do j=1,ntmp
          ii=i+ictot
          jj=j+ictot
          do k=1,3
          do l=1,nat
            dvec(k,l,ii,jj)=dptmp(k,l,i,j)
          enddo
          enddo
        enddo
        enddo
        ictot=ictot+ntmp
      enddo

      do i=1,nat
        dx(i)=gpema(1,i,1)
        dy(i)=gpema(2,i,1)
        dz(i)=gpema(3,i,1)
c       print *,i,dx(i),dy(i),dz(i)
      enddo

      return

      end
c **********************************************************************
c **********************************************************************





c **********************************************************************
c **********************************************************************
      subroutine prepot
      return
      end
c **********************************************************************
c **********************************************************************




c **********************************************************************
c **********************************************************************
      subroutine dd_g03(symbol,tnam,x,y,z,pema,gpema,dvec,nclu,
     &                              mnclu,nsurf,mnsurf,ifail)

c NOTE: SINGLE SURFACE FOR NOW, NO NA COUPLING

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
      double precision pema(mnsurf),gpema(3,mnclu,mnsurf)
      double precision dvec(3,mnclu,mnsurf,mnsurf)
      character*8 tnam

c read template & write QC input file
      open(unit=7,file=tnam)       ! template file
      open(unit=10,file='qc.in')   ! temporary input file
 122  format(a8,2x,3f20.10)
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
          do j=1,nclu
c          if (symbol(j).eq."Ho") write(10,122)"H-HO",x(j),y(j),z(j)
c          if (symbol(j).eq."C") write(10,122)"C-CT",x(j),y(j),z(j)
c          if (symbol(j).eq."O") write(10,122)"O-OH",x(j),y(j),z(j)
c          if (symbol(j).eq."H") write(10,122)"H-H1",x(j),y(j),z(j)
          write(10,122)symbol(j),x(j),y(j),z(j)
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
c        print *,"in",i,dx(i),dy(i),dz(i)
      enddo

c organize things
      do i=1,nsurf
        pema(i)=v
        do j=1,nclu
          gpema(1,j,i)=dx(j)
          gpema(2,j,i)=dy(j)
          gpema(3,j,i)=dz(j)
c          print *,"in2",i,gpema(1,j,i),gpema(2,j,i),gpema(3,j,i)
        enddo
      enddo

      end
c **********************************************************************
c **********************************************************************



c **********************************************************************
c **********************************************************************
      subroutine dd_m06(symbol,tnam,x,y,z,pema,gpema,dvec,nclu,
     &                                 mnclu,nsurf,mnsurf,ifail)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,idum,ithis(mnclu),nsurf,mnsurf,
     &          isurf,jsurf,ifail
      character*2 symbol(mnclu),cdum
      character*80 string,tmpstr
      character*21 Estr
      character*18 Gstr
      character*18 Astr
      character*4 Nstr
      double precision xk,dx1,dy1,dz1
      double precision x(mnclu),y(mnclu),z(mnclu),v(mnsurf)
      double precision gpema(3,mnclu,mnsurf),pema(mnsurf)
      double precision dvec(3,mnclu,mnsurf,mnsurf)
      double precision xtmp,ytmp,ztmp,total,dum
      double precision dx(mnclu,mnsurf),dy(mnclu,mnsurf),
     &                                  dz(mnclu,mnsurf)
      double precision cx(mnclu,mnsurf,mnsurf),cy(mnclu,mnsurf,mnsurf),
     &                                         cz(mnclu,mnsurf,mnsurf)
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


c READ ENERGIES & GRADIENTS

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
        endif
      enddo
      goto 300

 320  write(6,*)"Couldn't find string '",gstr,"' in qc.out"
      ifail=4
      go to 999

 399  continue

      ENDDO

c READ NA COUPLINGS

      DO ISURF=1,NSURF
      DO JSURF=ISURF+1,NSURF

c     get coupling
      Nstr='MOLD'
      lstr=len(Nstr)
 400  read (8,fmt='(a80)',end=420) string
      do i=1,10
        if (string(i:i+lstr-1).eq.Nstr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)xk,dx1,dy1,dz1
            k=nint(xk)
            cx(ithis(k),isurf,jsurf)=dx1
            cy(ithis(k),isurf,jsurf)=dy1
            cz(ithis(k),isurf,jsurf)=dz1
          enddo
          goto 499
        endif
      enddo
      goto 400

 420  write(6,*)"Couldn't find string '",nstr,"' in qc.out"
      ifail=5
      go to 999


 499  continue

      ENDDO
      ENDDO

      close(8)

c organize things
      do i=1,nsurf
        pema(i)=v(i)
        do j=1,nclu
          gpema(1,j,i)=dx(j,i)
          gpema(2,j,i)=dy(j,i)
          gpema(3,j,i)=dz(j,i)
          dvec(1,j,i,i)=0.d0
          dvec(2,j,i,i)=0.d0
          dvec(3,j,i,i)=0.d0
          do k=i+1,nsurf
            dvec(1,j,i,k)=cx(j,i,k)
            dvec(2,j,i,k)=cy(j,i,k)
            dvec(3,j,i,k)=cz(j,i,k)
            dvec(1,j,k,i)=-cx(j,i,k)
            dvec(2,j,k,i)=-cy(j,i,k)
            dvec(3,j,k,i)=-cz(j,i,k)
          enddo
        enddo
      enddo

999   continue
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
