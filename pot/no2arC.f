      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)

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

c TMP
c      dimension ccc(15)
c TMP

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

c      print *,1,symb(1)
c      print *,2,symb(2)
c      print *,3,symb(3)
c      print *,4,symb(4)

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
      if (natom3.ne.0) then
c        call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
        call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
      endif
      if (natom2.ne.0) then
        call no2(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
      endif

      tmpprint(1)=v

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


c      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
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

c TMP
      dimension ccc(15)
c TMP

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

       ccc(1) =7.5047456282235174
       ccc(2) =  0.27335978270821254
       ccc(3) =   7.8613238929410691
       ccc(4) =   1.9245582445753349
       ccc(5) =   7.9206213568529318
       ccc(6) =  0.23934049501022370
       ccc(7) =   7.7472151860103153
       ccc(8) =   2.5243232520523700

       ccc(1) =7.5002594073305460     
       ccc(2) =0.27177220984527117     
       ccc(3) =7.6946867275002289     
       ccc(4) =1.7557359538560138     
       ccc(5) =8.0000305185094760     
       ccc(6) =0.23250129703665273     
       ccc(7) =7.8164311655018768     
       ccc(8) =2.8247993408001952     


      if (m2.eq.23) then  ! Ar-
        if (m1.eq.3) then !    N
        aa=ccc(1)
        bb=ccc(2)
        cc=ccc(3)
        rrc=ccc(4)
c        aa=7.75022431
c        bb=0.259142125
c        cc=7.87899411
c        rrc=2.21254921
        cutoff=.true.
        aa=10.d0**aa
        elseif (m1.eq.4) then !    O
        aa=ccc(5)
        bb=ccc(6)
        cc=ccc(7)
        rrc=ccc(8)
c        aa=7.72470473
c        bb=0.271240883
c        cc=8.65068819
c        rrc=2.25741142
        cutoff=.true.
        aa=10.d0**aa
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


**** NO2 ****

c      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
      subroutine no2(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom),r(3)
      dimension dvdx(maxatom),xx(3,maxatom)
      dimension dvdy(maxatom)
      dimension dvdz(maxatom)
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      character*2 symb(maxatom)

      v1=0.d0
      do i=1,natom
      dvdx(i)=0.d0
      dvdy(i)=0.d0
      dvdz(i)=0.d0
      enddo

      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      r(1)=dsqrt(dx*dx+dy*dy+dz*dz)
      r(1)=r(1)*autoang

      dx=x(2)-x(3)
      dy=y(2)-y(3)
      dz=z(2)-z(3)
      r(2)=dsqrt(dx*dx+dy*dy+dz*dz)
      r(2)=r(2)*autoang

      dx=x(3)-x(1)
      dy=y(3)-y(1)
      dz=z(3)-z(1)
      r(3)=dsqrt(dx*dx+dy*dy+dz*dz)
      r(3)=r(3)*autoang

c test
c       rno2=1.19369055d0
c       ano2=133.2153283d0
c       roo2=2.d0*rno2**2
c     & -2.d0*rno2**2*dcos(ano2/180.d0*dacos(-1.d0))
c       roo2=dsqrt(roo2)
c       r(1)=rno2
c       r(2)=rno2
c       r(3)=roo2

c      print *,r,v
      call no2fit(r,v)
      v=v/autocmi

c      print *,r,v,v*autocmi,v*27.211d0

      resp=0.0001d0

      do i1=1,3
      if (i1.eq.1) i=1
      if (i1.eq.1) j=2
      if (i1.eq.2) i=2
      if (i1.eq.2) j=3
      if (i1.eq.3) i=3
      if (i1.eq.3) j=1
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      dx=dx*autoang
      dy=dy*autoang
      dz=dz*autoang
      r(i1)=r(i1)+resp
      call no2fit(r,vp)
      r(i1)=r(i1)-2.d0*resp
      call no2fit(r,vm)
      r(i1)=r(i1)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)/autocmi*autoang
      dvdx(i) = dvdx(i) + dtmpdrr*dx/r(i1)
      dvdx(j) = dvdx(j) - dtmpdrr*dx/r(i1)
      dvdy(i) = dvdy(i) + dtmpdrr*dy/r(i1)
      dvdy(j) = dvdy(j) - dtmpdrr*dy/r(i1)
      dvdz(i) = dvdz(i) + dtmpdrr*dz/r(i1)
      dvdz(j) = dvdz(j) - dtmpdrr*dz/r(i1)
      enddo

      return

      end

c **************************************************
c zero at NO2
c r1 = no1
c r2 = no2
c r3 = oo
c A and cm-1

      subroutine no2fit(r,v)

      implicit double precision(a-h,o-z)
      dimension r(3),basis(1000),coef(1000)

       r1=r(1)
       r2=r(2)
       r3=r(3)

c ZERO
       re=1.1508d0 ! A
       de=53372. ! cm-1
       a1=2.78984954
       a2=0.36732078
       a3=0.24997711
       rno2=1.19369055d0
       ano2=133.2153283d0
       roo2=2.d0*rno2**2
     & -2.d0*rno2**2*dcos(ano2/180.d0*dacos(-1.d0))
       roo2=dsqrt(roo2)
       rr=r1-re
       aa=a1+a2*rr+a3*rr**2
       s1=de*(1.d0-dexp(-aa*rr))**2-de
       rr=r2-re
       aa=a1+a2*rr+a3*rr**2
       s2=de*(1.d0-dexp(-aa*rr))**2-de
       ss=s1+s2
       rr=rno2-re
       aa=a1+a2*rr+a3*rr**2
       sx=de*(1.d0-dexp(-aa*rr))**2-de
       eno2=-25617.2d0
       eno2x=eno2-sx*2.d0

      coef(    1) =    -0.2575717111E+07
       coef(    2) =     0.7617900894E+08
       coef(    3) =    -0.1781055816E+08
       coef(    4) =    -0.4423316030E+10
       coef(    5) =     0.1192299674E+11
       coef(    6) =     0.1514992055E+12
       coef(    7) =    -0.5347138232E+12
       coef(    8) =     0.5473551056E+08
       coef(    9) =    -0.1621443853E+10
       coef(   10) =     0.1978304192E+10
       coef(   11) =     0.8884192065E+11
       coef(   12) =    -0.5719972457E+12
       coef(   13) =    -0.3046483552E+12
       coef(   14) =     0.3953007177E+13
       coef(   15) =    -0.4619076759E+09
       coef(   16) =     0.1325148366E+11
       coef(   17) =    -0.2987306402E+11
       coef(   18) =    -0.5594130878E+12
       coef(   19) =     0.5474511653E+13
       coef(   20) =    -0.9324539425E+13
       coef(   21) =    -0.5530705142E+13
       coef(   22) =     0.1929095293E+10
       coef(   23) =    -0.5104187328E+11
       coef(   24) =     0.1779596512E+12
       coef(   25) =     0.8248928204E+12
       coef(   26) =    -0.1827843266E+14
       coef(   27) =     0.5143902881E+14
       coef(   28) =    -0.1172369299E+14
       coef(   29) =    -0.4050270797E+10
       coef(   30) =     0.8596439797E+11
       coef(   31) =    -0.4496589914E+12
       coef(   32) =     0.3692737860E+13
       coef(   33) =     0.1451435779E+14
       coef(   34) =    -0.8320872745E+14
       coef(   35) =     0.3312350838E+14
       coef(   36) =     0.3654443417E+10
       coef(   37) =    -0.2339327957E+11
       coef(   38) =     0.3718479647E+12
       coef(   39) =    -0.1458548775E+14
       coef(   40) =     0.3761108069E+14
       coef(   41) =    -0.1306767330E+14
       coef(   42) =     0.1404887031E+14
       coef(   43) =    -0.6299478585E+09
       coef(   44) =    -0.5990944612E+11
       coef(   45) =     0.8515930425E+11
       coef(   46) =     0.1493268304E+14
       coef(   47) =    -0.6208497430E+14
       coef(   48) =     0.1153908747E+15
       coef(   49) =    -0.1118185371E+15
       coef(   50) =     0.5473218152E+08
       coef(   51) =    -0.1621260049E+10
       coef(   52) =     0.1976973099E+10
       coef(   53) =     0.8884587191E+11
       coef(   54) =    -0.5722212359E+12
       coef(   55) =    -0.3030352347E+12
       coef(   56) =     0.3951303975E+13
       coef(   57) =    -0.1130466630E+10
       coef(   58) =     0.3270641572E+11
       coef(   59) =    -0.1249862019E+12
       coef(   60) =    -0.6625890061E+12
       coef(   61) =     0.6890134225E+13
       coef(   62) =     0.2201539682E+13
       coef(   63) =    -0.2712185984E+14
       coef(   64) =     0.9058669917E+10
       coef(   65) =    -0.2321306264E+12
       coef(   66) =     0.1234460947E+13
       coef(   67) =     0.5433517849E+12
       coef(   68) =    -0.3695704701E+14
       coef(   69) =     0.2149641415E+14
       coef(   70) =     0.4464868253E+14
       coef(   71) =    -0.3491538997E+11
       coef(   72) =     0.6233540078E+12
       coef(   73) =    -0.4069108345E+13
       coef(   74) =     0.8972009589E+13
       coef(   75) =     0.9241079861E+14
       coef(   76) =    -0.1483872522E+15
       coef(   77) =    -0.5239469882E+14
       coef(   78) =     0.6284670406E+11
       coef(   79) =     0.1323957004E+12
       coef(   80) =    -0.1833054322E+12
       coef(   81) =    -0.3061113383E+14
       coef(   82) =    -0.1163319799E+15
       coef(   83) =     0.3120753861E+15
       coef(   84) =     0.8161334516E+14
       coef(   85) =    -0.3459887042E+11
       coef(   86) =    -0.3489083638E+13
       coef(   87) =     0.2675848238E+14
       coef(   88) =     0.3733423787E+14
       coef(   89) =     0.4684407560E+14
       coef(   90) =    -0.2208480306E+15
       coef(   91) =     0.5018453994E+14
       coef(   92) =    -0.1776975037E+11
       coef(   93) =     0.4436922474E+13
       coef(   94) =    -0.3933152991E+14
       coef(   95) =    -0.1778608315E+14
       coef(   96) =     0.7765926221E+14
       coef(   97) =    -0.1205540664E+15
       coef(   98) =     0.2350235929E+15
       coef(   99) =    -0.4618379175E+09
       coef(  100) =     0.1324719377E+11
       coef(  101) =    -0.2983537118E+11
       coef(  102) =    -0.5593492317E+12
       coef(  103) =     0.5475833764E+13
       coef(  104) =    -0.9341829150E+13
       coef(  105) =    -0.5515477285E+13
       coef(  106) =     0.9058350984E+10
       coef(  107) =    -0.2320927214E+12
       coef(  108) =     0.1234014373E+13
       coef(  109) =     0.5420238855E+12
       coef(  110) =    -0.3695150381E+14
       coef(  111) =     0.2159342412E+14
       coef(  112) =     0.4458671187E+14
       coef(  113) =    -0.7176204299E+11
       coef(  114) =     0.1359313268E+13
       coef(  115) =    -0.1101601692E+14
       coef(  116) =     0.2099258340E+14
       coef(  117) =     0.9433258861E+14
       coef(  118) =     0.1712315857E+14
       coef(  119) =     0.4278167727E+14
       coef(  120) =     0.2849699209E+12
       coef(  121) =    -0.1977952765E+13
       coef(  122) =     0.3674929642E+14
       coef(  123) =    -0.9098765369E+14
       coef(  124) =    -0.3796725772E+14
       coef(  125) =    -0.1179773999E+15
       coef(  126) =    -0.8124984361E+14
       coef(  127) =    -0.5732095287E+12
       coef(  128) =    -0.8652807260E+13
       coef(  129) =    -0.2349475060E+14
       coef(  130) =     0.1780537242E+15
       coef(  131) =    -0.1058006785E+15
       coef(  132) =     0.1335671468E+15
       coef(  133) =    -0.2542614846E+15
       coef(  134) =     0.5166384456E+12
       coef(  135) =     0.3386522514E+14
       coef(  136) =    -0.1698471012E+15
       coef(  137) =    -0.9714000331E+14
       coef(  138) =    -0.3623460639E+13
       coef(  139) =    -0.1083066891E+15
       coef(  140) =    -0.2673106751E+15
       coef(  141) =    -0.1322054099E+12
       coef(  142) =    -0.3258413957E+14
       coef(  143) =     0.3145562053E+15
       coef(  144) =    -0.2106910365E+15
       coef(  145) =     0.8924692858E+14
       coef(  146) =     0.1644777045E+13
       coef(  147) =    -0.4080120066E+14
       coef(  148) =     0.1928517063E+10
       coef(  149) =    -0.5100202471E+11
       coef(  150) =     0.1774933072E+12
       coef(  151) =     0.8249938166E+12
       coef(  152) =    -0.1827819896E+14
       coef(  153) =     0.5151095465E+14
       coef(  154) =    -0.1179617908E+14
       coef(  155) =    -0.3491495477E+11
       coef(  156) =     0.6230266489E+12
       coef(  157) =    -0.4063532677E+13
       coef(  158) =     0.8973506576E+13
       coef(  159) =     0.9227716969E+14
       coef(  160) =    -0.1486603584E+15
       coef(  161) =    -0.5220161313E+14
       coef(  162) =     0.2849951462E+12
       coef(  163) =    -0.1977609897E+13
       coef(  164) =     0.3672900052E+14
       coef(  165) =    -0.9091012194E+14
       coef(  166) =    -0.3775184011E+14
       coef(  167) =    -0.1175318813E+15
       coef(  168) =    -0.8152244172E+14
       coef(  169) =    -0.1169603046E+13
       coef(  170) =    -0.8313955636E+13
       coef(  171) =    -0.1361709843E+15
       coef(  172) =     0.6963005938E+13
       coef(  173) =    -0.6933371284E+14
       coef(  174) =     0.2770325496E+15
       coef(  175) =     0.2799523730E+15
       coef(  176) =     0.2003748210E+13
       coef(  177) =     0.6225470954E+14
       coef(  178) =     0.2229193992E+15
       coef(  179) =     0.1137342230E+15
       coef(  180) =    -0.4532461410E+14
       coef(  181) =     0.2938661796E+15
       coef(  182) =     0.6696567865E+14
       coef(  183) =    -0.3258770125E+12
       coef(  184) =    -0.1018310054E+15
       coef(  185) =     0.1373448418E+15
       coef(  186) =     0.2528714285E+15
       coef(  187) =     0.8393427625E+14
       coef(  188) =     0.3965375568E+14
       coef(  189) =    -0.9778880569E+14
       coef(  190) =    -0.1798385406E+13
       coef(  191) =     0.1252550726E+14
       coef(  192) =    -0.5623637223E+15
       coef(  193) =     0.2257877240E+15
       coef(  194) =     0.1763621420E+15
       coef(  195) =     0.2877705452E+14
       coef(  196) =    -0.5089562139E+14
       coef(  197) =    -0.4047909401E+10
       coef(  198) =     0.8578401213E+11
       coef(  199) =    -0.4469927363E+12
       coef(  200) =     0.3684967546E+13
       coef(  201) =     0.1451624480E+14
       coef(  202) =    -0.8337737761E+14
       coef(  203) =     0.3333134630E+14
       coef(  204) =     0.6285749940E+11
       coef(  205) =     0.1334032063E+12
       coef(  206) =    -0.2086072097E+12
       coef(  207) =    -0.3055020520E+14
       coef(  208) =    -0.1158801906E+15
       coef(  209) =     0.3119893247E+15
       coef(  210) =     0.8166376151E+14
       coef(  211) =    -0.5734164327E+12
       coef(  212) =    -0.8650454111E+13
       coef(  213) =    -0.2343547666E+14
       coef(  214) =     0.1776415376E+15
       coef(  215) =    -0.1061767186E+15
       coef(  216) =     0.1335331280E+15
       coef(  217) =    -0.2543645454E+15
       coef(  218) =     0.2004091509E+13
       coef(  219) =     0.6224603648E+14
       coef(  220) =     0.2229375294E+15
       coef(  221) =     0.1140122786E+15
       coef(  222) =    -0.4583585583E+14
       coef(  223) =     0.2936298033E+15
       coef(  224) =     0.6707177285E+14
       coef(  225) =     0.5659151381E+13
       coef(  226) =    -0.2080481541E+15
       coef(  227) =    -0.4816123537E+15
       coef(  228) =    -0.4335931761E+15
       coef(  229) =    -0.2265486311E+15
       coef(  230) =     0.1340712920E+15
       coef(  231) =    -0.1027020553E+14
       coef(  232) =    -0.4319445280E+14
       coef(  233) =     0.1671942914E+15
       coef(  234) =     0.5481383785E+14
       coef(  235) =    -0.1669422374E+15
       coef(  236) =    -0.1602194518E+15
       coef(  237) =    -0.5252609774E+14
       coef(  238) =    -0.9518719293E+14
       coef(  239) =     0.6204259255E+14
       coef(  240) =     0.2280312110E+15
       coef(  241) =     0.3322014341E+14
       coef(  242) =     0.1022799448E+15
       coef(  243) =    -0.2299273793E+14
       coef(  244) =    -0.5547915317E+14
       coef(  245) =    -0.7000000514E+14
       coef(  246) =     0.3649700281E+10
       coef(  247) =    -0.2299826229E+11
       coef(  248) =     0.3649535770E+12
       coef(  249) =    -0.1455080850E+14
       coef(  250) =     0.3754170923E+14
       coef(  251) =    -0.1274224149E+14
       coef(  252) =     0.1365617924E+14
       coef(  253) =    -0.3464492284E+11
       coef(  254) =    -0.3490227148E+13
       coef(  255) =     0.2680695230E+14
       coef(  256) =     0.3710027828E+14
       coef(  257) =     0.4654208266E+14
       coef(  258) =    -0.2205483103E+15
       coef(  259) =     0.4981802539E+14
       coef(  260) =     0.5171775751E+12
       coef(  261) =     0.3385403681E+14
       coef(  262) =    -0.1698650464E+15
       coef(  263) =    -0.9660036073E+14
       coef(  264) =    -0.3666028399E+13
       coef(  265) =    -0.1079371993E+15
       coef(  266) =    -0.2672557850E+15
       coef(  267) =    -0.3266311662E+12
       coef(  268) =    -0.1018187163E+15
       coef(  269) =     0.1371938317E+15
       coef(  270) =     0.2533882676E+15
       coef(  271) =     0.8365456326E+14
       coef(  272) =     0.3976423924E+14
       coef(  273) =    -0.9753991829E+14
       coef(  274) =    -0.4319385719E+14
       coef(  275) =     0.1672455070E+15
       coef(  276) =     0.5452107717E+14
       coef(  277) =    -0.1669157732E+15
       coef(  278) =    -0.1602900813E+15
       coef(  279) =    -0.5240932137E+14
       coef(  280) =    -0.9504880210E+14
       coef(  281) =     0.2196855555E+15
       coef(  282) =    -0.3758282197E+14
       coef(  283) =     0.3136101926E+15
       coef(  284) =    -0.1302296764E+15
       coef(  285) =    -0.1906336287E+15
       coef(  286) =    -0.1229882194E+15
       coef(  287) =    -0.1021633021E+15
       coef(  288) =    -0.3069262517E+15
       coef(  289) =    -0.2270371591E+15
       coef(  290) =     0.2820854869E+15
       coef(  291) =     0.2368461380E+14
       coef(  292) =    -0.1011018654E+15
       coef(  293) =    -0.9310937021E+14
       coef(  294) =    -0.6842852048E+14
       coef(  295) =    -0.6262058322E+09
       coef(  296) =    -0.6024364409E+11
       coef(  297) =     0.9169811862E+11
       coef(  298) =     0.1489005855E+14
       coef(  299) =    -0.6195728126E+14
       coef(  300) =     0.1150175317E+15
       coef(  301) =    -0.1113877070E+15
       coef(  302) =    -0.1771664600E+11
       coef(  303) =     0.4437155335E+13
       coef(  304) =    -0.3936474494E+14
       coef(  305) =    -0.1755905872E+14
       coef(  306) =     0.7753295170E+14
       coef(  307) =    -0.1203555621E+15
       coef(  308) =     0.2347574465E+15
       coef(  309) =    -0.1326870152E+12
       coef(  310) =    -0.3257299752E+14
       coef(  311) =     0.3145217264E+15
       coef(  312) =    -0.2109236557E+15
       coef(  313) =     0.8913952780E+14
       coef(  314) =     0.1918125906E+13
       coef(  315) =    -0.4077510355E+14
       coef(  316) =    -0.1797461524E+13
       coef(  317) =     0.1252045204E+14
       coef(  318) =    -0.5622707645E+15
       coef(  319) =     0.2256764899E+15
       coef(  320) =     0.1761558746E+15
       coef(  321) =     0.2893317808E+14
       coef(  322) =    -0.5070904136E+14
       coef(  323) =     0.6204022956E+14
       coef(  324) =     0.2280046528E+15
       coef(  325) =     0.3324969564E+14
       coef(  326) =     0.1021545667E+15
       coef(  327) =    -0.2303427601E+14
       coef(  328) =    -0.5533883332E+14
       coef(  329) =    -0.6986678795E+14
       coef(  330) =    -0.3069281039E+15
       coef(  331) =    -0.2269879205E+15
       coef(  332) =     0.2822499717E+15
       coef(  333) =     0.2367860086E+14
       coef(  334) =    -0.1010794805E+15
       coef(  335) =    -0.9305287244E+14
       coef(  336) =    -0.6838564735E+14
       coef(  337) =     0.4321187037E+15
       coef(  338) =    -0.3620209768E+15
       coef(  339) =     0.2454971307E+15
       coef(  340) =     0.3662846248E+14
       coef(  341) =    -0.7399633249E+14
       coef(  342) =    -0.7002953234E+14
       coef(  343) =    -0.4534095563E+14

       call getbasis(r,basis,ncoef)

c       print *
c       print *
c       print *
c       print *
c       print *
c       print *
c       print *
c       print *
       v=0.d0
       do j=1,ncoef
       v=v+coef(j)*basis(j)
c       print *,j,coef(j),basis(j)
       enddo
       v=v+ss+eno2x-eno2
c tmp
c       ang=(r(1)**2+r(2)**2-r(3)**2)/(2.d0*r(1)*r(2))
c       ang=dacos(ang)/dacos(-1.d0)*180.d0
c       print *,ang,r(1),r(2),v,v+eno2

       end

********************************************

      subroutine getbasis(r,basis,ncoef)

      implicit double precision(a-h,o-z)
      dimension r(3),basis(1000)

      norder=6
      norder2=6

      r1=r(1)
      r2=r(2)
      r3=r(3)
c      xq=-(r3**2-r1**2-r2**2)/(2.d0*r1*r2)
c      xq=dacos(xq)/dacos(-1.d0)*180.d0

c NO
       re=1.1508d0 ! A
       de=53372. ! cm-1
       a1=2.78984954
       a2=0.36732078
       a3=0.24997711

c NO2 min (calc)
       rno2=1.19369055d0
       ano2=133.2153283d0
       roo2=2.d0*rno2**2
     & -2.d0*rno2**2*dcos(ano2/180.d0*dacos(-1.d0))
       roo2=dsqrt(roo2)

c asymptotes
       rr=r1-re
       aa=a1+a2*rr+a3*rr**2
       s1=de*(1.d0-dexp(-aa*rr))**2-de
       rr=r2-re
       aa=a1+a2*rr+a3*rr**2
       s2=de*(1.d0-dexp(-aa*rr))**2-de
       ss=s1+s2

c zero
       rr=rno2-re
       aa=a1+a2*rr+a3*rr**2
       sx=de*(1.d0-dexp(-aa*rr))**2-de
       eno2=-25617.2d0-sx*2.d0

c fit
      b1=1.d0
      b3=1.d0

c      e1=1.d0-dexp(-(r1-rno2)/b1)
c      e2=1.d0-dexp(-(r2-rno2)/b1)
c      e3=1.d0-dexp(-(r3-roo2)/b3)
c      e3=1.d0-dexp(-(xq-ano2)/b3)
      e1=dexp(-r1/b1)
      e2=dexp(-r2/b1)
      e3=dexp(-r3/b3)
c      e1=(r1-rno2)
c      e2=(r2-rno2)
c      e3=dcos((xq-ano2)/180.d0*dacos(-1.d0))
c      e3=(xq-ano2)

c      print *,"e1",rno2,r1,e1
c      print *,"e2",rno2,r2,e2
c      print *,"e3",roo2,r3,e3

      ii=0
      e=0.d0
      do i1=0,norder2
      x1=e1**i1
      do i2=0,norder2
      x2=x1*(e2**i2)
      do i3=0,norder
      x3=x2*(e3**i3)
c      print *,i1,i2,i3,x1,x2,x3
      ii=ii+1
      basis(ii)=1.d0
      if (i1+i2+i3.ne.0) basis(ii)=x3
      enddo
      enddo
      enddo

      ncoef=ii

      return 
      end

***************************************************
      subroutine prepot
      return
      end
