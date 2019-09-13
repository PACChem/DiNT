
! This subroutine returns energies and analytic gradients for a 
! "universal" TB+exp/6 potential for CxHy + He, Ne, Kr, Ar, H2, N2, O2

! The "separable" approximation is used, where V = V_CxHy + V_M + V_int
! The V_int potential is further approximated as a sum of pairwise exp/6
! interactions, with parameters based on the CH4 + M system.

! V_CxHy: The tight binding hydrocarbon CxHy PES is from
! (1) Wang, Y.; Mak, C. H. Chem. Phys. Lett. 1995, 235, 37.
! Liu, T.; Truhlar, D. G.TB, version 1.0.1; University of Minnesota:
! Minnesota, 2004; see comp.chem.umn.edu/tbpac.
! Liu, T.Ph.D. Thesis, University of Minnesota, 2000.

! V_int: The "separable pairwise" atom-atom exp/6 parametrizations are from
! (2) Theoretical study of the Ar-, Kr-, and Xe-CH4, -CF4 intermolecular
! Potential-Energy Surfaces
! W. A. Alexander and D. Troya, J. Phys. Chem. A 110, 10834 (2006).
! (3) Theoretical unimolecular kinetics for CH4 + M → CH3 + H + M in 
! eight baths, M = He, Ne, Ar, Kr, H2, CO, N2, and CH4
! A. W. Jasper and J. A. Miller, J. Phys. Chem. A 115, 6438 (2011).
! (4) The collision efficiency of water in the unimolecular reaction CH4
! (+H2O) --> CH3 + H (+H2O): One-dimensional and two-dimensional
! solutions of the low-pressure-limit master equation
! A. W. Jasper, J. A. Miller, and S. J. Klippenstein, J. Phys. Chem. A
! 117, 12243 (2013).
! (5) Temperature- and pressure-dependent rate coefficients for the HACA
! pathways from benzene to naphthalene
! A. M. Mebel, Y. Georgievskii, A. W. Jasper, and S. J. Klippenstein,
! Proc. Combust. Inst., in press (2017:w
! (6) Kinetics of propargyl radical dissociation
! S. J. Klippenstein, J. A. Miller, and A. W. Jasper, J. Phys. Chem. A
! 119, 7780-7791 (2015).
!
! V_M: One-dimensional potentials for the diatomic baths were obtained
! by fitting modified Morse curves to ab initio data.
!
! In addition to the validations given in the above references, the 
! present exp/6 interaction potentials have been validated for and used to 
! predict energy transfer in larger systems:
! (6) A. W. Jasper, C. M. Oana, and J. A. Miller, Proc. Combust. Inst.
! 34, 197-204 (2015).
! and used to predict diffusion coefficients:
! (7) A. W. Jasper and J. M. Miller, Combust. Flame, 161, 101-110 (2014).

      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

! INPUT: 
! This subroutine follows POTLIB's "HE-MM-1" convention (comp.chem.umn.edu/potlib)
! SYMB: Array of 1 or 2 character atom labels. The labels follow special
! conventions. See below.
! MAXATOM: dimension of the X,Y,Z,DVDX,DVDY,DVDZ arrays
! NATOM: number of atoms
! X(1,NATOM): Array of x coordinates for each atom in bohr
! Y(1,NATOM): Array of y coordinates for each atom in bohr
! Z(1,NATOM): Array of z coordinates for each atom in bohr

! OUTPUT:
! V: Energy in hartree
! DVDX(1,NATOM): Array of x components of the gradient for each atom in hartree/bohr
! DVDY(1,NATOM): Array of y components of the gradient for each atom in hartree/bohr
! DVDZ(1,NATOM): Array of z components of the gradient for each atom in hartree/bohr

! ATOM LABEL CONVENTIONS:
! Atom labels are used to distinguish "bath" atoms from "target" atoms.
! The CxHy TB PES is evaluated for "target" atoms only. The interaction 
! PES is evaluated for every intermolecular pair of atoms. The diatomic
! bath gas PES is evaluated for the two atoms in the diatomic bath.
!     "Target" atoms should be labeled
!        C or c for Carbon
!        H or h for Hydrogen
!        Ca or ca for the radical carbon in naphthyl when used with the Ha bath
!     "Bath" atoms should be labeled
!        He, HE, or he for Helium
!        Ne, NE, or ne for Neon
!        Ar, AR, or ar for Argon
!        Kr, KR, or kr for Krypton
!        H2 or h2 for each H atom in the H2 bath
!        C2 or c2 for the C atom in the CO bath
!        O2 or o2 for the O atom in the CO bath
!        N2 or o2 for each N atom in the N2 bath
!        O3 or o3 for each O atom in the O2 bath
!        Note: By default the RADIAL fits from Ref 3 are returned for
!        the diatomic baths.
!      Special cases
!        PAH + He and Ar:
!        Label the bath Ha or ha for naphthyl radical and similar systems + He
!          ** The radical C in naphthyl should be labeled Ca **
!        Label the bath Aa or aa for naphthyl radical and similar systems + Ar
!        C3H3 + He and Ar:
!        Label the bath Hp or hp for propargyl radical + He
!        Label the bath Ap or ap for propargyl radical + Ar

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

      dimension tmpprint(50)
      common/tmp/tmpprint

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
     &    (symb(i).eq."n"))  at(i)=3       ! nitrogen
      if ((symb(i).eq."O").or.
     &    (symb(i).eq."o"))  at(i)=4       ! oxygen
      if ((symb(i).eq."He").or.
     &    (symb(i).eq."he").or.
     &    (symb(i).eq."HE")) at(i)=21      ! helium
      if ((symb(i).eq."Ne").or.
     &    (symb(i).eq."ne").or.
     &    (symb(i).eq."NE")) at(i)=22      ! neon
      if ((symb(i).eq."Ar").or.
     &    (symb(i).eq."ar").or.
     &    (symb(i).eq."AR")) at(i)=23      ! argon
      if ((symb(i).eq."Kr").or.
     &    (symb(i).eq."kr").or.
     &    (symb(i).eq."KR")) at(i)=24      ! krypton
      if ((symb(i).eq."H2").or.
     &    (symb(i).eq."h2")) at(i)=25      ! H2  ! label your bath atoms H2, not H
      if ((symb(i).eq."N2").or.
     &    (symb(i).eq."n2")) at(i)=26      ! N2  ! label your bath atoms N2, not N
      if ((symb(i).eq."C2").or.
     &    (symb(i).eq."c2")) at(i)=27      ! CO bath ! label atoms C2 and O2
      if ((symb(i).eq."O2").or.
     &    (symb(i).eq."c2")) at(i)=28      ! CO bath ! label atoms C2 and O2
      if ((symb(i).eq."O3").or.
     &    (symb(i).eq."o3")) at(i)=29      ! O2 bath ! label atoms O3, not O (the numbers are just labels; yes, it is confusing)

! Special cases
!       C2H4 + Ar
      if ((symb(i).eq."Ae").or.
     &    (symb(i).eq."ae")) at(i)=37

!       A3 + Ar
      if ((symb(i).eq."Aa").or.
     &    (symb(i).eq."aa")) at(i)=35
!       A3 + He
      if ((symb(i).eq."Ca").or.
     &    (symb(i).eq."ca"))  at(i)=5
      if ((symb(i).eq."Ha").or.
     &    (symb(i).eq."ha"))  at(i)=36

!       propargyl + He
      if ((symb(i).eq."Hp").or.
     &    (symb(i).eq."hp")) at(i)=32
!       propargyl + Ar
      if ((symb(i).eq."Ap").or.
     &    (symb(i).eq."ap")) at(i)=33

      if ((symb(i).eq."Hx").or.
     &    (symb(i).eq."hx").or.
     &    (symb(i).eq."HX")) at(i)=31      ! helium for C2H3
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
          if (at2(natom2).eq.5) symb2(natom2)="C"
        else ! collect bath atoms
          natom3=natom3+1
          x3(natom3)=x(i)
          y3(natom3)=y(i)
          z3(natom3)=z(i)
          symb3(natom3)=symb(i)
          at3(natom3)=at(i)
        endif
      enddo

      if (natom3.ne.0) then
        call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom) ! 1D diatomic baths
        call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom) ! Interaction PES
      endif
      if (natom2.ne.0)
     &   call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
     &          natom2,maxatom)

      tmpprint(1)=v    ! interaction
      tmpprint(3)=v2   ! target internal
      tmpprint(4)=v3   ! bath internal

      v=v+v2+v3

      natom2=0
      natom3=0
      do i=1,natom
        if (at(i).le.20) then
        natom2=natom2+1
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
!       agrees reasonably well with more complicated form of LeRoy (JCP 125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, 
!                                     PRL 53, 159 (1984), as quoted by LeRoy)
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
!       Morse. Fit to RKR data of PAUL H. KRUPENIE and STANLEY WEISSMAN, J. CHEM. PHYS. 43, 1529 (1965)
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

      if (m1.ge.21) then ! two bath gases, skip this pair
         go to 2
      endif
      if (m2.le.20) then ! no bath gas, skip this pair
         go to 2
      endif

      if (m2.eq.21) then      ! He-
        if     (m1.eq.1) then !    H
c       tinkered with to fit qcisd(t)/av5z CH4+He (Jasper, 2008)
          aa=6.040254d0
          bb=0.286639d0
          cc=5.124943d0
      ascale=1.25d0
      bscale=0.88d0
      cscale=0.88d0
          aa=ascale*(10.d0**aa)
          bb=bscale*bb
          cc=cscale*cc
        elseif (m1.eq.2) then !    C
c       tinkered with to fit qcisd(t)/av5z CH4+He (Jasper, 2008)
          aa=6.393873d0
          bb=0.277317d0
          cc=5.987628d0
      ascale=1.05d0
      bscale=1.040
      cscale=1.030
          aa=ascale*(10.d0**aa)
          bb=bscale*bb
          cc=cscale*cc
        else
          write(6,*)"Cant find a ?-He interaction"
          stop
        endif

      elseif (m2.eq.22) then  ! Ne-
        if     (m1.eq.1) then !    H
c       fit to cc qcisd(t)/cbs(adz,atz) CH4+Ne (Jasper, 2009)
          aa=6.39071017000
          bb=0.26216342700
          cc=5.74940641000
          aa=(10.d0**aa)
        elseif (m1.eq.2) then !    C
c       fit to cc qcisd(t)/cbs(adz,atz) CH4+Ne (Jasper, 2009)
          aa=7.49964293000
          bb=0.23913852400
          cc=4.80220038000
          aa=(10.d0**aa)
        else
          write(6,*)"Cant find a Ne-? interaction"
          stop
        endif
      elseif (m2.eq.23) then  ! Ar-
        if     (m1.eq.1) then !    H
c       Troya's fit to CH4+Ar  (J. Phys. Chem. 110, 10834 (2006))
          aa=11426.51d0
          bb=3.385d0
          cc=-374.119d0
          troya=.true.
        elseif (m1.eq.2) then !    C
c       Troya's fit to CH4+Ar  (J. Phys. Chem. 110, 10834 (2006))
          aa=96594.54d0
          bb=3.608d0
          cc=-356.575d0
          troya=.true.
        else
          write(6,*)"Cant find a Ar-? interaction"
          stop
        endif
      elseif (m2.eq.24) then  ! Kr-
        if     (m1.eq.1) then !    H
c       Troya's fit to CH4+Kr  (J. Phys. Chem. 110, 10834 (2006))
          aa=13754.02d0
          bb=3.238d0
          cc=-621.784d0
          troya=.true.
c         reoptimized for i-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**4.19044465
c          bb=3.41540574
c          cc=-350.001526
c          troya=.true.
c         reoptimized for n-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**4.27774285
c          bb=3.47558214
c          cc=-457.029939
c          troya=.true.
        elseif (m1.eq.2) then !    C
c       Troya's fit to CH4+Kr  (J. Phys. Chem. 110, 10834 (2006))
          aa=112927.4d0
          bb=3.520d0
          cc=-268.460d0
          troya=.true.
c         reoptimized for i-proply + KR, QCISD(T)/CBS CC
c           aa=10.**5.6093936
c           bb=3.91874142
c           cc=-1102.49641
c          troya=.true.
c         reoptimized for n-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**5.10611591
c          bb=3.51407819
c          cc=-683.828242
c          troya=.true.
        else
          write(6,*)"Cant find a Kr-? interaction"
          stop
        endif
      elseif (m2.eq.25) then  ! H2-
        if     (m1.eq.1) then !    H
c       aa=4.91209143  ! compromise radial, perpendicular fit
c       bb=0.386484268
c       cc=6.27030549
c       rrc=2.50001526
       aa=5.23924375 ! radial only fit
       bb=0.31757561
       cc=5.53859066
       rrc=2.62514115
c       aa=5.68126469 ! perp only
c       bb=0.277195654
c       cc=5.45313883
c       rrc=2.55272073
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
c       aa=5.1396527 ! compromise radial, perp fit
c       bb=0.428481094
c       cc=3.24997711
c       rrc=2.94143498
       aa=6.91601306 ! radial only fit
       bb=0.224359874
       cc=5.47520371
       rrc=2.46021912
c       aa=6.34962615 ! perp only
c       bb=0.289641102
c       cc=1.51567125
c       rrc=3.
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a H2-? interaction"
          stop
        endif
      elseif (m2.eq.29) then  ! O2-
        if     (m1.eq.1) then !    H
         aa=6.04197821
         bb=0.288725394
         cc=5.99391156
         rrc=2.63286233
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
         aa=7.57335124
         bb=0.245461135
         cc=6.92338328
         rcc=2.03799554
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a N2-? interaction"
          stop
        endif
      elseif (m2.eq.26) then  ! N2-
        if     (m1.eq.1) then !    H
         aa=6.52842799 ! radial
         bb=0.275934629
         cc=6.45893735
         rrc=0.942960906
c         aa=6.35029145 ! perp only fit
c         bb=0.273452864
c         cc=0.0175786615
c         rrc=8.70967742
c         aa=6.50016785  ! compromise fit to radial and perp
c         bb=0.259999695
c         cc=4.49601733
c         rrc=3.09872738
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
         aa=7.14549699 ! radial
         bb=0.286131779
         cc=0.0288094729
         rrc=6.05728324
c         aa=6.74993133  ! perp only fit
c         bb=0.300780358
c         cc=8.43803217
c         rrc=2.84856716
c         aa=7.43794061  ! compromise fit to radial and per
c         bb=0.265233619
c         cc=8.18771935
c         rrc=3.0496231
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a N2-? interaction"
          stop
        endif
      elseif (m2.eq.27) then  ! C of a CO bath-
        if     (m1.eq.1) then !    H
         aa=  6.67976928  ! radial O-in
         bb=  0.292242805
         cc=  0.738883633
         rrc=1.45118564
c         aa=7.12680441  ! compromise
c         bb=0.22606708
c         cc=6.56242561
c         rrc=3.28461562
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
         aa=  7.1884518 ! radial O-in
         bb=  0.238867763
         cc=  8.74794763
         rrc=2.50135807
c         aa=7.500412   ! compromise 
c         bb=0.276250801
c         cc=1.68080691
c         rrc=7.13440352
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a CO-? interaction"
          stop
        endif
      elseif (m2.eq.28) then  ! O of a CO bath-
        if     (m1.eq.1) then !    H
         aa=  4.99996948 ! radial O-in
         bb=  0.413102512
         cc=  0.90664388
         rrc=8.74816126
c         aa=4.99996948  ! compromise
c         bb=0.398960234
c         cc=1.99975585
c         rrc=2.4582049
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
         aa=  6.9365215 ! radial O-in
         bb=  0.30321543
         cc=  8.98471023
         rrc=2.29291665
c         aa=5.99981689  ! compromise
c         bb=0.339911191
c         cc=8.74993133
c         rrc=3.00070193
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a CO-? interaction"
          stop
        endif

      elseif (m2.eq.31) then      ! He-C2H3
        if     (m1.eq.1) then !    H
        aa =  6.92322459
        bb =  0.193029115
        cc =  4.99956359
        rrc =  1.27176122
        cutoff=.true.
        aa=10.d0**aa
        elseif (m1.eq.2) then !    C
        aa =  6.35630665
        bb =  0.309332865
        cc =  5.90316477
        rrc =  2.76535539
        cutoff=.true.
        aa=10.d0**aa
        else
          write(6,*)"Cant find a ?-Hx interaction"
          stop
        endif

      elseif (m2.eq.32) then      ! He-C3H3
        if     (m1.eq.1) then !    propagyl H
          aa=6.7498397778252510
          bb=0.26839808343760491
          cc=5.2494888149662771
          rrc=2.2638019959105198
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    propagyl C
          aa=6.1233558153019807
          bb=0.21866359447004607
          cc=6.0000000000000000
          rrc=3.5030671102023376
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Hp interaction"
          stop
        endif

      elseif (m2.eq.33) then      ! Ar-C3H3
        if     (m1.eq.1) then !    propagyl H
          aa=7.4028138065736870
          bb=0.28616351817377239
          cc=8.2187261574144728
          rrc=2.8802758873256629
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    propagyl C
          aa=7.6640827661976987
          bb=0.19703238013855404
          cc=6.7231360820337533
          rrc=3.0937223426007874
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Ap interaction"
          stop
        endif

      elseif (m2.eq.37) then      ! Ar-C2H4
        if     (m1.eq.1) then !    H
          aa=  6.9052705465865047     
          bb=  0.25637623218482009     
          cc=  6.6610919522690510     
          rrc= 1.7641834772789697     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !    C
          aa=   8.0212714011047694     
          bb=  0.25264870143742180     
          cc=   7.9999694814905240     
          rrc=   3.4531693472090823     
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Ap interaction"
          stop
        endif

      elseif (m2.eq.36) then      ! He+A3
        if     (m1.eq.1) then !   H
           aa=6.2186193426313059     
           bb=0.22906155583361310     
           cc=5.1328775902584916     
           rrc=3.7504196295052949     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !   C
          aa=6.7790765099032564     
          bb=0.26019440290536211     
          cc=5.8884548478652299     
          rrc=0.46217230750450150     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.5) then !   C
          aa=6.5000152592547380     
          bb=0.29316537980285040     
          cc=2.9682607501449629     
          rrc=0.16724143192846461     
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Ha interaction"
          stop
        endif

      elseif (m2.eq.34) then      ! N2+A3
        if     (m1.eq.1) then !   H
c          aa=5.9743339335306862   ! first fit
c          bb=0.30863979003265479
c          cc=5.8125553147984252
c          rrc=1.7656483657338176

          aa=5.9833979308450580      ! 2nd fit to avg of cc & non-cc QM
          bb=0.30161137730033266     
          cc=5.3140354625080111     
          rrc=1.8628803369243445     
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !   C
c          aa=7.4881435590685754     ! first fit
c          bb=0.26138676107058934
c          cc=7.7470320749534594
c          rrc=1.9371929074983978

          aa=7.0429700613422037      ! 2nd fit to avg of cc & non-cc QM
          bb=0.28108615375225071     
          cc=7.7815485091708121     
          rrc=2.6391796624652852     
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Ap interaction"
          stop
        endif

      elseif (m2.eq.35) then      ! Ar+A3
        if     (m1.eq.1) then !   H
          aa=5.887490463
          bb=0.339483627
          cc=7.124729148
          rrc=2.500961333
          aa=(10.d0**aa)
          cutoff=.true.
        elseif (m1.eq.2) then !   C
          aa=8.011932737
          bb=0.254977569
          cc=8.37653737
          rrc=0.384533219
          aa=(10.d0**aa)
          cutoff=.true.
        else
          write(6,*)"Cant find a ?-Ap interaction"
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

