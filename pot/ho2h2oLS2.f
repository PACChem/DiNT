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
      dimension ccc(15)
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
      if ((symb(i).eq."Ow").or.
     &    (symb(i).eq."ow")) at(i)=30      ! water bath
      if ((symb(i).eq."Hw").or.
     &    (symb(i).eq."hw")) at(i)=31      ! water bath
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
      if (natom3.ne.0.and.natom2.ne.0) then
c        call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
c        call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
        call lsbath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
      endif
      if (natom2.ne.0) then
        call ho2(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
      endif
      if (natom3.ne.0) then
        call water(x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
      endif

      tmpprint(1)=v

c      print *,v*autocmi,v2*autocmi,v3*autocmi,natom,natom2,natom3

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


      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
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

       if (m2.eq.21) then
        ccc(1) =   6.2813196203497421     
        ccc(2) =  0.18875026703695791     
        ccc(3) =   4.8749046296578875     
        ccc(4) =   2.5312967314676351     
        ccc(5) =   6.7499313333536790     
        ccc(6) =  0.25311593981749930     
        ccc(7) =   5.4999542222357860     
        ccc(8) =   2.0362865077669605     
       elseif(m2.eq.23) then
        ccc(1)=   7.1964476454969937
        ccc(2)=  0.19953154087954345
        ccc(3)=   6.4964751121555224
        ccc(4)=   2.3697012237922301
        ccc(5)=   7.4057130649739067
        ccc(6)=  0.27445326090273753
        ccc(7)=   7.9421979430524612
        ccc(8)=   2.5391094698934902
       else
        print *,'bad bath in ho2he.f'
        stop
       endif

      if (m2.eq.21.or.m2.eq.23) then  ! He-
        if (m1.eq.1) then !    H
        aa=ccc(1)
        bb=ccc(2)
        cc=ccc(3)
        rrc=ccc(4)
        cutoff=.true.
        aa=10.d0**aa
        elseif (m1.eq.4) then !    O
        aa=ccc(5)
        bb=ccc(6)
        cc=ccc(7)
        rrc=ccc(8)
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


!-----------------------------------------------------------------------
      subroutine ho2(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

! O O H

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

      dx=x(3)-x(2)
      dy=y(3)-y(2)
      dz=z(3)-z(2)
      r(1)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oo
      r(1)=r(1)*autoang

      dx=x(2)-x(1)
      dy=y(2)-y(1)
      dz=z(2)-z(1)
      r(2)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oh1
      r(2)=r(2)*autoang

      dx=x(3)-x(1)
      dy=y(3)-y(1)
      dz=z(3)-z(1)
      r(3)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oh2
      r(3)=r(3)*autoang

c      print *,"ho2",(r(k),k=1,3)

      call ho2fit(r,v)
      v=v/autocmi

      resp=0.0001d0

      do i1=1,3
      if (i1.eq.1) i=3
      if (i1.eq.1) j=2
      if (i1.eq.2) i=2
      if (i1.eq.2) j=1
      if (i1.eq.3) i=3
      if (i1.eq.3) j=1
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      dx=dx*autoang
      dy=dy*autoang
      dz=dz*autoang
      r(i1)=r(i1)+resp
      call ho2fit(r,vp)
      r(i1)=r(i1)-2.d0*resp
      call ho2fit(r,vm)
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
c r1 = oo
c r2 = ho1
c r3 = ho2
c A and cm-1

      subroutine ho2fit(r,v)

      implicit double precision(a-h,o-z)
      dimension r(3),basis(1000),coef(1000)

       r1=r(1)
       r2=r(2)
       r3=r(3)

       coef(    1) =     0.2885755880E+07
       coef(    2) =    -0.5187870513E+08
       coef(    3) =     0.3913454355E+09
       coef(    4) =    -0.1565117811E+10
       coef(    5) =     0.3485101857E+10
       coef(    6) =    -0.4096488749E+10
       coef(    7) =     0.1990397144E+10
       coef(    8) =    -0.5187501892E+08
       coef(    9) =     0.1038949244E+10
       coef(   10) =    -0.8485105642E+10
       coef(   11) =     0.3608041911E+11
       coef(   12) =    -0.8424083243E+11
       coef(   13) =     0.1024915854E+12
       coef(   14) =    -0.5083830889E+11
       coef(   15) =     0.3912897735E+09
       coef(   16) =    -0.8484890012E+10
       coef(   17) =     0.7567156007E+11
       coef(   18) =    -0.3537875717E+12
       coef(   19) =     0.9124059122E+12
       coef(   20) =    -0.1233830930E+13
       coef(   21) =     0.6867787126E+12
       coef(   22) =    -0.1564785372E+10
       coef(   23) =     0.3608030860E+11
       coef(   24) =    -0.3538188152E+12
       coef(   25) =     0.1858845657E+13
       coef(   26) =    -0.5436257680E+13
       coef(   27) =     0.8425218741E+13
       coef(   28) =    -0.5440821210E+13
       coef(   29) =     0.3484120106E+10
       coef(   30) =    -0.8424611841E+11
       coef(   31) =     0.9126420581E+12
       coef(   32) =    -0.5437150914E+13
       coef(   33) =     0.1808696302E+14
       coef(   34) =    -0.3198733745E+14
       coef(   35) =     0.2365337013E+14
       coef(   36) =    -0.4095055159E+10
       coef(   37) =     0.1025087128E+12
       coef(   38) =    -0.1234410751E+13
       coef(   39) =     0.8427842898E+13
       coef(   40) =    -0.3199071923E+14
       coef(   41) =     0.6404172242E+14
       coef(   42) =    -0.5360346252E+14
       coef(   43) =     0.1989569008E+10
       coef(   44) =    -0.5085387177E+11
       coef(   45) =     0.6872486139E+12
       coef(   46) =    -0.5442966401E+13
       coef(   47) =     0.2365643909E+14
       coef(   48) =    -0.5360446417E+14
       coef(   49) =     0.5131882408E+14
       coef(   50) =    -0.6299627645E+08
       coef(   51) =     0.1080106688E+10
       coef(   52) =    -0.7745801965E+10
       coef(   53) =     0.2970615520E+11
       coef(   54) =    -0.6408633613E+11
       coef(   55) =     0.7355745552E+11
       coef(   56) =    -0.3503807167E+11
       coef(   57) =     0.1079988181E+10
       coef(   58) =    -0.1775566904E+11
       coef(   59) =     0.1230595352E+12
       coef(   60) =    -0.4404332763E+12
       coef(   61) =     0.8270964012E+12
       coef(   62) =    -0.7263202999E+12
       coef(   63) =     0.1933461150E+12
       coef(   64) =    -0.7743630753E+10
       coef(   65) =     0.1230391319E+12
       coef(   66) =    -0.8729146678E+12
       coef(   67) =     0.3218662189E+13
       coef(   68) =    -0.6039384311E+13
       coef(   69) =     0.4755381922E+13
       coef(   70) =    -0.4943945698E+12
       coef(   71) =     0.2969065771E+11
       coef(   72) =    -0.4402545713E+12
       coef(   73) =     0.3218514175E+13
       coef(   74) =    -0.1280249256E+14
       coef(   75) =     0.2621026888E+14
       coef(   76) =    -0.2422298382E+14
       coef(   77) =     0.6366899469E+13
       coef(   78) =    -0.6403275868E+11
       coef(   79) =     0.8264377134E+12
       coef(   80) =    -0.6039890284E+13
       coef(   81) =     0.2621811452E+14
       coef(   82) =    -0.6003961523E+14
       coef(   83) =     0.7387266913E+14
       coef(   84) =    -0.4378871279E+14
       coef(   85) =     0.7346776534E+11
       coef(   86) =    -0.7251902913E+12
       coef(   87) =     0.4758258363E+13
       coef(   88) =    -0.2424826077E+14
       coef(   89) =     0.7390622207E+14
       coef(   90) =    -0.1466172707E+15
       coef(   91) =     0.1473719204E+15
       coef(   92) =    -0.3497979472E+11
       coef(   93) =     0.1926099445E+12
       coef(   94) =    -0.4974867034E+12
       coef(   95) =     0.6388237713E+13
       coef(   96) =    -0.4381736105E+14
       coef(   97) =     0.1473776345E+15
       coef(   98) =    -0.2074958409E+15
       coef(   99) =     0.5800278032E+09
       coef(  100) =    -0.9329047338E+10
       coef(  101) =     0.6575012449E+11
       coef(  102) =    -0.2624346026E+12
       coef(  103) =     0.6156815904E+12
       coef(  104) =    -0.7906648391E+12
       coef(  105) =     0.4299032162E+12
       coef(  106) =    -0.9328528213E+10
       coef(  107) =     0.1210615235E+12
       coef(  108) =    -0.7064634557E+12
       coef(  109) =     0.2256831494E+13
       coef(  110) =    -0.3948542170E+13
       coef(  111) =     0.3357944004E+13
       coef(  112) =    -0.9530402922E+12
       coef(  113) =     0.6573814051E+11
       coef(  114) =    -0.7063122994E+12
       coef(  115) =     0.3769605575E+13
       coef(  116) =    -0.1097989210E+14
       coef(  117) =     0.1506848602E+14
       coef(  118) =    -0.1755316979E+13
       coef(  119) =    -0.1099782788E+14
       coef(  120) =    -0.2623313858E+12
       coef(  121) =     0.2255099854E+13
       coef(  122) =    -0.1097697048E+14
       coef(  123) =     0.3305619386E+14
       coef(  124) =    -0.4621926909E+14
       coef(  125) =    -0.5179611767E+13
       coef(  126) =     0.4980137306E+14
       coef(  127) =     0.6152766137E+12
       coef(  128) =    -0.3940708057E+13
       coef(  129) =     0.1505788294E+14
       coef(  130) =    -0.4623742528E+14
       coef(  131) =     0.4393558985E+14
       coef(  132) =     0.7675493552E+14
       coef(  133) =    -0.1066172476E+15
       coef(  134) =    -0.7899252952E+12
       coef(  135) =     0.3342394771E+13
       coef(  136) =    -0.1741184038E+13
       coef(  137) =    -0.5115695035E+13
       coef(  138) =     0.7667170223E+14
       coef(  139) =    -0.2039492167E+15
       coef(  140) =     0.8976746701E+14
       coef(  141) =     0.4293911666E+12
       coef(  142) =    -0.9417862195E+12
       coef(  143) =    -0.1100292004E+14
       coef(  144) =     0.4974879181E+14
       coef(  145) =    -0.1065681306E+15
       coef(  146) =     0.8978815222E+14
       coef(  147) =     0.1836124186E+15
       coef(  148) =    -0.2840444893E+10
       coef(  149) =     0.4265799645E+11
       coef(  150) =    -0.3043826020E+12
       coef(  151) =     0.1329893125E+13
       coef(  152) =    -0.3440430819E+13
       coef(  153) =     0.4723703176E+13
       coef(  154) =    -0.2676487825E+13
       coef(  155) =     0.4266271618E+11
       coef(  156) =    -0.3989172259E+12
       coef(  157) =     0.1855415280E+13
       coef(  158) =    -0.6159583199E+13
       coef(  159) =     0.1368777062E+14
       coef(  160) =    -0.1707469643E+14
       coef(  161) =     0.9304597067E+13
       coef(  162) =    -0.3044410397E+12
       coef(  163) =     0.1855526083E+13
       coef(  164) =    -0.5309612110E+13
       coef(  165) =     0.1045160931E+14
       coef(  166) =    -0.1424280220E+14
       coef(  167) =    -0.2895227641E+13
       coef(  168) =     0.2438823608E+14
       coef(  169) =     0.1330116137E+13
       coef(  170) =    -0.6156135065E+13
       coef(  171) =     0.1044792201E+14
       coef(  172) =    -0.9790798343E+13
       coef(  173) =     0.2988117173E+14
       coef(  174) =    -0.1240816315E+13
       coef(  175) =    -0.7785469466E+14
       coef(  176) =    -0.3440730689E+13
       coef(  177) =     0.1366054595E+14
       coef(  178) =    -0.1420714642E+14
       coef(  179) =     0.2989879467E+14
       coef(  180) =    -0.1268350522E+15
       coef(  181) =     0.1077000011E+15
       coef(  182) =    -0.2855188574E+14
       coef(  183) =     0.4723650074E+13
       coef(  184) =    -0.1700807976E+14
       coef(  185) =    -0.2976087112E+13
       coef(  186) =    -0.1269567100E+13
       coef(  187) =     0.1077236890E+15
       coef(  188) =    -0.2127941123E+14
       coef(  189) =    -0.3100240455E+14
       coef(  190) =    -0.2676175787E+13
       coef(  191) =     0.9251164799E+13
       coef(  192) =     0.2443666963E+14
       coef(  193) =    -0.7784872242E+14
       coef(  194) =    -0.2849143148E+14
       coef(  195) =    -0.3107928142E+14
       coef(  196) =    -0.5358522920E+14
       coef(  197) =     0.7758652059E+10
       coef(  198) =    -0.1085001335E+12
       coef(  199) =     0.7964588772E+12
       coef(  200) =    -0.3786954287E+13
       coef(  201) =     0.1003498132E+14
       coef(  202) =    -0.1282628934E+14
       coef(  203) =     0.6441317325E+13
       coef(  204) =    -0.1085444794E+12
       coef(  205) =     0.5816102562E+12
       coef(  206) =    -0.1547406925E+13
       coef(  207) =     0.8677643066E+13
       coef(  208) =    -0.2958212176E+14
       coef(  209) =     0.3954116675E+14
       coef(  210) =    -0.2300780301E+14
       coef(  211) =     0.7971279379E+12
       coef(  212) =    -0.1551049041E+13
       coef(  213) =    -0.6375456772E+13
       coef(  214) =     0.1541953334E+14
       coef(  215) =     0.1167591359E+14
       coef(  216) =    -0.1202794960E+14
       coef(  217) =    -0.1623867393E+14
       coef(  218) =    -0.3790840769E+13
       coef(  219) =     0.8688916725E+13
       coef(  220) =     0.1539372814E+14
       coef(  221) =    -0.1011717662E+15
       coef(  222) =     0.9301569060E+14
       coef(  223) =    -0.5225037392E+14
       coef(  224) =     0.1192797668E+15
       coef(  225) =     0.1004663008E+14
       coef(  226) =    -0.2956521114E+14
       coef(  227) =     0.1168904457E+14
       coef(  228) =     0.9297446324E+14
       coef(  229) =    -0.4310817190E+14
       coef(  230) =    -0.4891520413E+13
       coef(  231) =     0.7577271690E+14
       coef(  232) =    -0.1284388760E+14
       coef(  233) =     0.3943760247E+14
       coef(  234) =    -0.1197684873E+14
       coef(  235) =    -0.5225682763E+14
       coef(  236) =    -0.4894435872E+13
       coef(  237) =    -0.9187642408E+14
       coef(  238) =    -0.4603851749E+14
       coef(  239) =     0.6451748734E+13
       coef(  240) =    -0.2290438263E+14
       coef(  241) =    -0.1624829948E+14
       coef(  242) =     0.1193009698E+15
       coef(  243) =     0.7575914764E+14
       coef(  244) =    -0.4612207290E+14
       coef(  245) =    -0.6796529482E+14
       coef(  246) =    -0.1121292927E+11
       coef(  247) =     0.1448601620E+12
       coef(  248) =    -0.1093171788E+13
       coef(  249) =     0.5409505681E+13
       coef(  250) =    -0.1268539338E+14
       coef(  251) =     0.9538003778E+13
       coef(  252) =     0.3886786968E+12
       coef(  253) =     0.1449853812E+12
       coef(  254) =    -0.1049358110E+12
       coef(  255) =    -0.2236229495E+13
       coef(  256) =    -0.3779134094E+13
       coef(  257) =     0.2439675907E+14
       coef(  258) =     0.1170680909E+14
       coef(  259) =    -0.3199933755E+14
       coef(  260) =    -0.1095109649E+13
       coef(  261) =    -0.2225878334E+13
       coef(  262) =     0.2585569250E+14
       coef(  263) =    -0.1896817929E+14
       coef(  264) =    -0.1055770422E+15
       coef(  265) =     0.5462131890E+14
       coef(  266) =     0.7360365528E+13
       coef(  267) =     0.5421429778E+13
       coef(  268) =    -0.3824196121E+13
       coef(  269) =    -0.1892366523E+14
       coef(  270) =     0.9084602852E+14
       coef(  271) =     0.3694532578E+14
       coef(  272) =    -0.7339420729E+14
       coef(  273) =     0.2637316409E+14
       coef(  274) =    -0.1272363422E+14
       coef(  275) =     0.2445666055E+14
       coef(  276) =    -0.1056328884E+15
       coef(  277) =     0.3694447014E+14
       coef(  278) =    -0.5354465236E+14
       coef(  279) =    -0.1252590700E+14
       coef(  280) =     0.5528976920E+14
       coef(  281) =     0.9598936259E+13
       coef(  282) =     0.1173352538E+14
       coef(  283) =     0.5467462754E+14
       coef(  284) =    -0.7334352714E+14
       coef(  285) =    -0.1252716552E+14
       coef(  286) =     0.3307742222E+13
       coef(  287) =     0.2755978025E+14
       coef(  288) =     0.3512359718E+12
       coef(  289) =    -0.3208200634E+14
       coef(  290) =     0.7241072695E+13
       coef(  291) =     0.2638026379E+14
       coef(  292) =     0.5525696945E+14
       coef(  293) =     0.2750482317E+14
       coef(  294) =     0.1133597241E+14
       coef(  295) =     0.6728522777E+10
       coef(  296) =    -0.7884579248E+11
       coef(  297) =     0.6009183039E+12
       coef(  298) =    -0.2872782016E+13
       coef(  299) =     0.3918911423E+13
       coef(  300) =     0.7991960432E+13
       coef(  301) =    -0.1395512644E+14
       coef(  302) =    -0.7896567352E+11
       coef(  303) =    -0.3980782704E+12
       coef(  304) =     0.4072253370E+13
       coef(  305) =    -0.3555893516E+13
       coef(  306) =     0.7779564972E+13
       coef(  307) =    -0.1071672439E+15
       coef(  308) =     0.1180192909E+15
       coef(  309) =     0.6027717067E+12
       coef(  310) =     0.4063511218E+13
       coef(  311) =    -0.2504524927E+14
       coef(  312) =     0.9382231514E+13
       coef(  313) =     0.1242182663E+15
       coef(  314) =     0.4380269384E+14
       coef(  315) =    -0.8280241315E+14
       coef(  316) =    -0.2884434781E+13
       coef(  317) =    -0.3517810428E+13
       coef(  318) =     0.9396945706E+13
       coef(  319) =    -0.1033234024E+15
       coef(  320) =     0.4180202567E+13
       coef(  321) =    -0.1175448305E+15
       coef(  322) =    -0.1254519005E+15
       coef(  323) =     0.3957361841E+13
       coef(  324) =     0.7709797372E+13
       coef(  325) =     0.1241486997E+15
       coef(  326) =     0.4172806469E+13
       coef(  327) =     0.5700487477E+14
       coef(  328) =     0.4655506031E+14
       coef(  329) =     0.2529933821E+14
       coef(  330) =     0.7930028835E+13
       coef(  331) =    -0.1071392895E+15
       coef(  332) =     0.4388428195E+14
       coef(  333) =    -0.1175078399E+15
       coef(  334) =     0.4654701086E+14
       coef(  335) =     0.8483517390E+14
       coef(  336) =     0.6305474245E+14
       coef(  337) =    -0.1391745862E+14
       coef(  338) =     0.1180581082E+15
       coef(  339) =    -0.8280683684E+14
       coef(  340) =    -0.1254445313E+15
       coef(  341) =     0.2526975674E+14
       coef(  342) =     0.6302423885E+14
       coef(  343) =     0.4606668298E+14

       call getbasis(r,basis,ncoef)

       v=0.d0
       do j=1,ncoef
       v=v+coef(j)*basis(j)
       enddo

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

      b1=1.d0
      b3=1.d0

      e1=dexp(-r1/b1)
      e2=dexp(-r2/b3)
      e3=dexp(-r3/b3)

      ii=0
      e=0.d0
      do i1=0,norder
      x1=e1**i1
      do i2=0,norder2
      x2=x1*(e2**i2)
      do i3=0,norder2
      x3=x2*(e3**i3)
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

***************************************************

      subroutine lsbath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)
      integer at(maxatom)

      dimension r(9)

      ii=0
      do i=1,3
      do j=4,6
      ii=ii+1
      ij=ii
      if (ii.eq.5) ij=6
      if (ii.eq.6) ij=7
      if (ii.eq.7) ij=5
      r(ij)=dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
      enddo
      enddo

      call lsbathfit(r,v)

      resp=0.0001d0

      ii=0
      do i=1,3
      do j=4,6
      ii=ii+1
      ij=ii
      if (ii.eq.5) ij=6
      if (ii.eq.6) ij=7
      if (ii.eq.7) ij=5
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(ij)=r(ij)+resp
      call lsbathfit(r,vp)
      r(ij)=r(ij)-2.d0*resp
      call lsbathfit(r,vm)
      r(ij)=r(ij)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)
      dvdx(i) = dvdx(i) + dtmpdrr*dx/r(ij)
      dvdx(j) = dvdx(j) - dtmpdrr*dx/r(ij)
      dvdy(i) = dvdy(i) + dtmpdrr*dy/r(ij)
      dvdy(j) = dvdy(j) - dtmpdrr*dy/r(ij)
      dvdz(i) = dvdz(i) + dtmpdrr*dz/r(ij)
      dvdz(j) = dvdz(j) - dtmpdrr*dz/r(ij)
      enddo
      enddo

      return
      end

      subroutine lsbathfit(r,v)

      implicit real*8(a-h,o-z)

      parameter(autoang=0.529177249d0)
      parameter(autocmi=219474.63067d0)

      dimension r(9)
      dimension coef(1100),basis(1100)

        ncoef = 1057 
       coef(    1) =     0.3983694314E+00
       coef(    2) =    -0.3002771123E+03
       coef(    3) =     0.1664598547E+05
       coef(    4) =    -0.2023627874E+06
       coef(    5) =     0.1232232177E+07
       coef(    6) =     0.2118287061E+05
       coef(    7) =    -0.1129691306E+07
       coef(    8) =     0.6510976485E+07
       coef(    9) =    -0.3328913208E+08
       coef(   10) =     0.9149446197E+07
       coef(   11) =    -0.1023777403E+09
       coef(   12) =     0.1896565787E+09
       coef(   13) =     0.2220225970E+09
       coef(   14) =    -0.3136541985E+09
       coef(   15) =    -0.2247014991E+06
       coef(   16) =     0.1484610349E+08
       coef(   17) =    -0.1142789104E+09
       coef(   18) =     0.4531564683E+09
       coef(   19) =     0.1105980361E+08
       coef(   20) =     0.7414112849E+09
       coef(   21) =    -0.2432041386E+10
       coef(   22) =    -0.1076023289E+10
       coef(   23) =    -0.2483610217E+09
       coef(   24) =    -0.1311541569E+09
       coef(   25) =    -0.1035242199E+08
       coef(   26) =     0.2328122786E+09
       coef(   27) =    -0.4542134074E+09
       coef(   28) =     0.6209059935E+09
       coef(   29) =    -0.1257669745E+10
       coef(   30) =    -0.1921123136E+10
       coef(   31) =     0.5255715951E+09
       coef(   32) =     0.2591907529E+03
       coef(   33) =    -0.1209341182E+06
       coef(   34) =     0.9403360633E+06
       coef(   35) =    -0.3045590522E+07
       coef(   36) =    -0.1735632419E+08
       coef(   37) =     0.2687444398E+07
       coef(   38) =    -0.4569945164E+08
       coef(   39) =     0.3361230245E+09
       coef(   40) =    -0.3445106908E+09
       coef(   41) =     0.5529372991E+08
       coef(   42) =    -0.1165325705E+10
       coef(   43) =     0.1752869224E+10
       coef(   44) =     0.1247247621E+10
       coef(   45) =    -0.6989571607E+07
       coef(   46) =    -0.2875547838E+09
       coef(   47) =    -0.3925483942E+09
       coef(   48) =    -0.9291373698E+09
       coef(   49) =     0.3437613989E+09
       coef(   50) =     0.1144514369E+11
       coef(   51) =    -0.1830695242E+10
       coef(   52) =     0.4712375051E+08
       coef(   53) =     0.1335291337E+09
       coef(   54) =     0.2113759273E+10
       coef(   55) =     0.7076749770E+10
       coef(   56) =     0.1688982607E+05
       coef(   57) =     0.2526600363E+07
       coef(   58) =    -0.1053041772E+08
       coef(   59) =     0.2887877194E+08
       coef(   60) =     0.1262384314E+09
       coef(   61) =    -0.4567101635E+08
       coef(   62) =     0.5463835097E+09
       coef(   63) =    -0.2016517970E+10
       coef(   64) =     0.2666981784E+10
       coef(   65) =     0.9583006994E+09
       coef(   66) =    -0.8058530363E+10
       coef(   67) =     0.3392389221E+09
       coef(   68) =     0.5961361098E+09
       coef(   69) =    -0.5749117120E+10
       coef(   70) =    -0.1039513403E+11
       coef(   71) =    -0.7075025640E+09
       coef(   72) =    -0.5233759592E+10
       coef(   73) =    -0.8018945831E+06
       coef(   74) =    -0.4930775140E+07
       coef(   75) =    -0.6642869872E+07
       coef(   76) =     0.2245842885E+09
       coef(   77) =    -0.1136921556E+10
       coef(   78) =     0.2492375779E+08
       coef(   79) =    -0.2515650704E+10
       coef(   80) =     0.6229088514E+10
       coef(   81) =     0.8267803694E+10
       coef(   82) =    -0.1684187327E+10
       coef(   83) =     0.1490142237E+11
       coef(   84) =     0.3083852779E+10
       coef(   85) =     0.5050317941E+07
       coef(   86) =    -0.3330245070E+08
       coef(   87) =    -0.2702576503E+09
       coef(   88) =     0.1059126070E+10
       coef(   89) =     0.1534595895E+10
       coef(   90) =    -0.7335331150E+10
       coef(   91) =    -0.2830796962E+10
       coef(   92) =     0.6669175710E+05
       coef(   93) =    -0.1883915893E+07
       coef(   94) =     0.4818526887E+07
       coef(   95) =     0.9687677503E+07
       coef(   96) =     0.7637976402E+08
       coef(   97) =     0.4843757008E+08
       coef(   98) =    -0.2348816194E+09
       coef(   99) =     0.4253749777E+09
       coef(  100) =    -0.2004643975E+10
       coef(  101) =     0.4311736345E+09
       coef(  102) =    -0.1941432359E+10
       coef(  103) =    -0.2622077836E+09
       coef(  104) =     0.8014042529E+09
       coef(  105) =     0.1170040931E+10
       coef(  106) =     0.3234418295E+10
       coef(  107) =    -0.2671162011E+09
       coef(  108) =     0.6460134693E+10
       coef(  109) =    -0.2021167060E+07
       coef(  110) =     0.4496422487E+07
       coef(  111) =    -0.1965075438E+09
       coef(  112) =    -0.1340115624E+09
       coef(  113) =     0.5924049436E+09
       coef(  114) =    -0.2364631264E+09
       coef(  115) =     0.2193822615E+10
       coef(  116) =     0.3728362459E+10
       coef(  117) =    -0.1103348805E+09
       coef(  118) =     0.2983595151E+10
       coef(  119) =    -0.2175271914E+11
       coef(  120) =    -0.7722251322E+10
       coef(  121) =     0.3443839732E+08
       coef(  122) =     0.4918008692E+08
       coef(  123) =     0.1091597678E+10
       coef(  124) =    -0.3199416158E+10
       coef(  125) =    -0.2703510101E+10
       coef(  126) =    -0.3560683672E+10
       coef(  127) =     0.8655677094E+10
       coef(  128) =    -0.1862223148E+09
       coef(  129) =     0.1113796896E+10
       coef(  130) =     0.1013706526E+10
       coef(  131) =     0.2633182080E+09
       coef(  132) =     0.2407401540E+07
       coef(  133) =     0.1951501406E+09
       coef(  134) =    -0.3782921727E+09
       coef(  135) =     0.6011584030E+08
       coef(  136) =    -0.1179337757E+10
       coef(  137) =    -0.2025338980E+10
       coef(  138) =     0.2281570287E+10
       coef(  139) =    -0.1760149492E+09
       coef(  140) =    -0.3518342791E+09
       coef(  141) =     0.5289460104E+10
       coef(  142) =     0.1179492368E+11
       coef(  143) =     0.1076665053E+10
       coef(  144) =    -0.1061239466E+11
       coef(  145) =     0.1673541659E+09
       coef(  146) =    -0.3230607406E+10
       coef(  147) =     0.9973997956E+09
       coef(  148) =     0.2662256476E+03
       coef(  149) =    -0.9491619718E+04
       coef(  150) =     0.2485409987E+06
       coef(  151) =    -0.3800907365E+07
       coef(  152) =     0.1141832486E+08
       coef(  153) =     0.3692735500E+06
       coef(  154) =    -0.3199339831E+08
       coef(  155) =     0.1869425794E+09
       coef(  156) =    -0.5601205622E+09
       coef(  157) =     0.3308204317E+09
       coef(  158) =    -0.1267176867E+10
       coef(  159) =     0.2186743991E+10
       coef(  160) =    -0.1116971826E+10
       coef(  161) =     0.2456342303E+08
       coef(  162) =    -0.5435578848E+09
       coef(  163) =    -0.4448663731E+09
       coef(  164) =     0.2669639150E+10
       coef(  165) =     0.1004584381E+09
       coef(  166) =     0.5548932494E+10
       coef(  167) =     0.3283162564E+10
       coef(  168) =    -0.4373042408E+09
       coef(  169) =     0.7281233748E+10
       coef(  170) =    -0.5663513468E+10
       coef(  171) =     0.4003681393E+10
       coef(  172) =     0.3502702901E+05
       coef(  173) =    -0.2794703366E+07
       coef(  174) =     0.5484286284E+08
       coef(  175) =    -0.2013197916E+09
       coef(  176) =     0.5120617266E+09
       coef(  177) =     0.5897740946E+08
       coef(  178) =    -0.6940587032E+08
       coef(  179) =     0.4223700031E+09
       coef(  180) =    -0.5140378871E+10
       coef(  181) =    -0.1763437302E+10
       coef(  182) =     0.5835348783E+09
       coef(  183) =     0.9554793758E+09
       coef(  184) =    -0.5726073158E+10
       coef(  185) =     0.1630461347E+11
       coef(  186) =    -0.1284856591E+11
       coef(  187) =    -0.4996612457E+10
       coef(  188) =    -0.1479860673E+11
       coef(  189) =    -0.8241760494E+06
       coef(  190) =    -0.1971576065E+08
       coef(  191) =    -0.6928058555E+09
       coef(  192) =     0.1719992725E+10
       coef(  193) =    -0.1311562328E+10
       coef(  194) =    -0.1945247913E+10
       coef(  195) =     0.1165283068E+11
       coef(  196) =     0.6803607568E+10
       coef(  197) =     0.1222214177E+10
       coef(  198) =     0.8704389685E+10
       coef(  199) =    -0.2747018370E+11
       coef(  200) =     0.1408274548E+11
       coef(  201) =     0.1587772603E+08
       coef(  202) =     0.1108487407E+10
       coef(  203) =    -0.1260753956E+10
       coef(  204) =    -0.5530024589E+10
       coef(  205) =    -0.4690410974E+10
       coef(  206) =    -0.1188106525E+11
       coef(  207) =    -0.3510312426E+10
       coef(  208) =    -0.1688821108E+09
       coef(  209) =     0.1074763209E+09
       coef(  210) =     0.5508069984E+10
       coef(  211) =     0.3382262448E+10
       coef(  212) =     0.2693295624E+07
       coef(  213) =    -0.1063114063E+09
       coef(  214) =     0.2892522903E+09
       coef(  215) =    -0.1724612258E+09
       coef(  216) =     0.7029482942E+09
       coef(  217) =     0.6339148773E+09
       coef(  218) =    -0.2033340962E+10
       coef(  219) =    -0.8483141458E+10
       coef(  220) =     0.1828462835E+11
       coef(  221) =    -0.1803230745E+10
       coef(  222) =    -0.2821769921E+11
       coef(  223) =    -0.9561958791E+10
       coef(  224) =    -0.2947253577E+07
       coef(  225) =     0.1325382740E+10
       coef(  226) =    -0.3419248620E+10
       coef(  227) =     0.1267584100E+10
       coef(  228) =    -0.8779331419E+10
       coef(  229) =     0.4629585121E+10
       coef(  230) =     0.7638511235E+11
       coef(  231) =    -0.4052039892E+09
       coef(  232) =     0.2506980967E+10
       coef(  233) =     0.8046178907E+10
       coef(  234) =    -0.2137506819E+11
       coef(  235) =    -0.1902781795E+09
       coef(  236) =     0.4806733587E+09
       coef(  237) =    -0.1584036359E+09
       coef(  238) =     0.1041316496E+10
       coef(  239) =     0.7727928938E+10
       coef(  240) =    -0.2077456519E+10
       coef(  241) =     0.1008928694E+10
       coef(  242) =    -0.2210118427E+11
       coef(  243) =     0.4563917095E+10
       coef(  244) =     0.2836901720E+09
       coef(  245) =     0.3471701525E+04
       coef(  246) =     0.7904092403E+05
       coef(  247) =    -0.4918653010E+07
       coef(  248) =     0.4874443340E+08
       coef(  249) =    -0.8545768692E+08
       coef(  250) =    -0.9605246871E+07
       coef(  251) =    -0.1746432326E+09
       coef(  252) =     0.1199721565E+10
       coef(  253) =     0.1156616283E+10
       coef(  254) =     0.4123581266E+09
       coef(  255) =    -0.1995558795E+10
       coef(  256) =     0.4726052631E+09
       coef(  257) =    -0.2438146622E+10
       coef(  258) =    -0.4300573828E+10
       coef(  259) =    -0.3764680513E+10
       coef(  260) =    -0.2390887527E+10
       coef(  261) =     0.1719998496E+08
       coef(  262) =    -0.3610494917E+05
       coef(  263) =    -0.8880953322E+07
       coef(  264) =     0.1112277881E+09
       coef(  265) =    -0.1833514823E+10
       coef(  266) =     0.8952054949E+09
       coef(  267) =    -0.8775975193E+09
       coef(  268) =     0.4731904833E+10
       coef(  269) =    -0.1507962019E+11
       coef(  270) =     0.6652239892E+10
       coef(  271) =     0.7695187669E+10
       coef(  272) =     0.1549364642E+10
       coef(  273) =    -0.7702765942E+10
       coef(  274) =     0.6034890111E+07
       coef(  275) =     0.1483456017E+10
       coef(  276) =    -0.1337097042E+10
       coef(  277) =     0.4255224311E+10
       coef(  278) =    -0.1444039988E+11
       coef(  279) =     0.1456148246E+11
       coef(  280) =     0.1305999953E+11
       coef(  281) =    -0.5234858098E+09
       coef(  282) =     0.4050034255E+10
       coef(  283) =    -0.6001689758E+10
       coef(  284) =     0.2918311771E+09
       coef(  285) =     0.1119469272E+09
       coef(  286) =     0.2815755120E+09
       coef(  287) =     0.1465007397E+08
       coef(  288) =    -0.5278444433E+09
       coef(  289) =     0.4176179823E+09
       coef(  290) =     0.5480803321E+10
       coef(  291) =     0.6091493398E+10
       coef(  292) =    -0.9688513403E+10
       coef(  293) =     0.9590096720E+10
       coef(  294) =    -0.1732438409E+09
       coef(  295) =    -0.1683175279E+10
       coef(  296) =     0.3409280805E+10
       coef(  297) =    -0.2596421357E+11
       coef(  298) =     0.2392505013E+10
       coef(  299) =     0.2151290233E+10
       coef(  300) =    -0.4478150132E+10
       coef(  301) =     0.1880560132E+09
       coef(  302) =    -0.2859293187E+10
       coef(  303) =     0.2687373393E+10
       coef(  304) =    -0.1200169797E+06
       coef(  305) =     0.1195980840E+08
       coef(  306) =    -0.9966014672E+06
       coef(  307) =    -0.1393823364E+09
       coef(  308) =     0.8679271914E+07
       coef(  309) =    -0.1691566055E+09
       coef(  310) =     0.2977744112E+09
       coef(  311) =     0.1470006131E+09
       coef(  312) =    -0.9672354724E+09
       coef(  313) =     0.1126685381E+09
       coef(  314) =     0.1781346869E+10
       coef(  315) =     0.3464160929E+10
       coef(  316) =    -0.8972307696E+05
       coef(  317) =    -0.6585653773E+08
       coef(  318) =     0.8708134706E+09
       coef(  319) =     0.1494779138E+10
       coef(  320) =     0.9140009203E+09
       coef(  321) =     0.5554828191E+10
       coef(  322) =    -0.4814562770E+10
       coef(  323) =    -0.6327533271E+08
       coef(  324) =    -0.4694559446E+09
       coef(  325) =    -0.1714940266E+10
       coef(  326) =    -0.4291953006E+10
       coef(  327) =     0.6228758212E+09
       coef(  328) =     0.1701939113E+10
       coef(  329) =    -0.9120065025E+09
       coef(  330) =     0.3004837111E+08
       coef(  331) =     0.7520053392E+09
       coef(  332) =    -0.7552794644E+10
       coef(  333) =    -0.9528892234E+10
       coef(  334) =    -0.6095929037E+09
       coef(  335) =     0.1542980456E+11
       coef(  336) =    -0.9421557696E+09
       coef(  337) =    -0.2447082743E+10
       coef(  338) =     0.5601730008E+05
       coef(  339) =    -0.2045200542E+08
       coef(  340) =    -0.1328944958E+08
       coef(  341) =     0.6004204610E+08
       coef(  342) =     0.5904427186E+09
       coef(  343) =    -0.1102206847E+10
       coef(  344) =    -0.4292156359E+09
       coef(  345) =     0.2488817943E+08
       coef(  346) =    -0.6720962765E+09
       coef(  347) =    -0.1267827769E+10
       coef(  348) =    -0.2884110649E+09
       coef(  349) =     0.1117857947E+09
       coef(  350) =     0.2575767603E+10
       coef(  351) =    -0.1244318505E+10
       coef(  352) =    -0.3674274125E+08
       coef(  353) =     0.2104232949E+10
       coef(  354) =    -0.1263179069E+10
       coef(  355) =     0.2984699742E+04
       coef(  356) =     0.9744629783E+06
       coef(  357) =    -0.1132363289E+08
       coef(  358) =     0.8383796798E+08
       coef(  359) =    -0.1099329280E+09
       coef(  360) =    -0.1454077308E+08
       coef(  361) =    -0.5290262641E+08
       coef(  362) =     0.9264438449E+09
       coef(  363) =     0.3150174112E+09
       coef(  364) =     0.4206102521E+09
       coef(  365) =    -0.1501304628E+10
       coef(  366) =     0.3468708208E+09
       coef(  367) =    -0.7354993574E+09
       coef(  368) =    -0.1155455953E+11
       coef(  369) =     0.6511259536E+09
       coef(  370) =    -0.1699918158E+10
       coef(  371) =    -0.5295768091E+10
       coef(  372) =    -0.3322615990E+06
       coef(  373) =    -0.1652134067E+08
       coef(  374) =     0.6423150030E+07
       coef(  375) =    -0.1831153386E+10
       coef(  376) =     0.3253159081E+10
       coef(  377) =    -0.7151325043E+09
       coef(  378) =     0.2123546036E+10
       coef(  379) =    -0.3255635531E+10
       coef(  380) =    -0.4746179832E+10
       coef(  381) =     0.3695198917E+10
       coef(  382) =     0.1642291974E+11
       coef(  383) =    -0.5342477302E+10
       coef(  384) =     0.9843768006E+07
       coef(  385) =     0.1182200451E+10
       coef(  386) =    -0.1946854495E+09
       coef(  387) =    -0.5224028885E+10
       coef(  388) =    -0.7319540689E+10
       coef(  389) =     0.7668294800E+10
       coef(  390) =     0.3364417095E+10
       coef(  391) =    -0.3463877457E+09
       coef(  392) =     0.1927256641E+10
       coef(  393) =     0.1871277838E+10
       coef(  394) =    -0.3941881049E+10
       coef(  395) =    -0.1702405940E+09
       coef(  396) =     0.1640029147E+10
       coef(  397) =     0.1821017882E+08
       coef(  398) =    -0.2937260310E+09
       coef(  399) =     0.1067492525E+09
       coef(  400) =     0.1014804875E+10
       coef(  401) =     0.4923374950E+10
       coef(  402) =     0.6384865990E+10
       coef(  403) =     0.3579884467E+10
       coef(  404) =    -0.2247965722E+09
       coef(  405) =    -0.6094627634E+09
       coef(  406) =     0.9239494519E+10
       coef(  407) =    -0.2927024482E+11
       coef(  408) =     0.1012129006E+10
       coef(  409) =     0.6161036051E+10
       coef(  410) =    -0.4161778074E+10
       coef(  411) =    -0.2279712479E+09
       coef(  412) =    -0.9723210319E+10
       coef(  413) =     0.6863118733E+10
       coef(  414) =    -0.3226098499E+06
       coef(  415) =     0.7484727868E+07
       coef(  416) =     0.1852735032E+09
       coef(  417) =    -0.6648363827E+09
       coef(  418) =    -0.1797601288E+09
       coef(  419) =     0.1645576266E+09
       coef(  420) =    -0.4574497914E+10
       coef(  421) =     0.2413171233E+11
       coef(  422) =    -0.5417930473E+10
       coef(  423) =    -0.1871689869E+10
       coef(  424) =    -0.1544896293E+11
       coef(  425) =     0.9131058321E+10
       coef(  426) =     0.1597046484E+08
       coef(  427) =    -0.7834664299E+09
       coef(  428) =     0.5634005915E+10
       coef(  429) =    -0.1107378249E+11
       coef(  430) =     0.9865045336E+10
       coef(  431) =     0.9878240286E+10
       coef(  432) =     0.3762242374E+11
       coef(  433) =    -0.1152070060E+09
       coef(  434) =    -0.6313702157E+10
       coef(  435) =     0.1600157353E+11
       coef(  436) =    -0.2900874160E+11
       coef(  437) =     0.1728377471E+10
       coef(  438) =     0.6670520813E+10
       coef(  439) =    -0.3012437213E+10
       coef(  440) =     0.2201078991E+09
       coef(  441) =    -0.4552975164E+10
       coef(  442) =     0.5585282290E+09
       coef(  443) =    -0.1696672841E+11
       coef(  444) =     0.6352234089E+09
       coef(  445) =    -0.7937230498E+10
       coef(  446) =     0.6827595392E+10
       coef(  447) =    -0.3123580438E+09
       coef(  448) =    -0.2361642085E+07
       coef(  449) =    -0.2686690806E+08
       coef(  450) =    -0.3151447482E+09
       coef(  451) =    -0.4392677355E+09
       coef(  452) =     0.4487628320E+09
       coef(  453) =     0.1305432405E+10
       coef(  454) =     0.4128330065E+10
       coef(  455) =     0.5859181567E+08
       coef(  456) =     0.2071588505E+09
       coef(  457) =    -0.7046988245E+10
       coef(  458) =    -0.4686303650E+11
       coef(  459) =     0.2713091270E+09
       coef(  460) =     0.2243632808E+11
       coef(  461) =    -0.5457996623E+10
       coef(  462) =    -0.1476788832E+09
       coef(  463) =     0.1448692476E+11
       coef(  464) =    -0.3819103468E+09
       coef(  465) =     0.1079567906E+08
       coef(  466) =    -0.2331064432E+09
       coef(  467) =     0.4761899622E+09
       coef(  468) =     0.1428437989E+10
       coef(  469) =    -0.1214783327E+09
       coef(  470) =     0.5946670344E+10
       coef(  471) =    -0.7625644839E+09
       coef(  472) =    -0.2459911722E+10
       coef(  473) =    -0.1022572003E+07
       coef(  474) =     0.4182705767E+08
       coef(  475) =    -0.2940740877E+09
       coef(  476) =    -0.1011127626E+10
       coef(  477) =     0.2349492493E+09
       coef(  478) =     0.4112977553E+10
       coef(  479) =     0.3857394981E+10
       coef(  480) =     0.3767104337E+08
       coef(  481) =    -0.2937040562E+08
       coef(  482) =    -0.3273944533E+10
       coef(  483) =    -0.5220349283E+11
       coef(  484) =     0.1280770594E+09
       coef(  485) =     0.2261019657E+11
       coef(  486) =    -0.5142831867E+10
       coef(  487) =    -0.4107073758E+09
       coef(  488) =     0.1194074752E+11
       coef(  489) =     0.4415027911E+10
       coef(  490) =    -0.1380633292E+08
       coef(  491) =    -0.6919528455E+09
       coef(  492) =     0.1846247245E+10
       coef(  493) =     0.1168663236E+11
       coef(  494) =    -0.8068992037E+08
       coef(  495) =     0.9253752931E+10
       coef(  496) =    -0.1100214197E+10
       coef(  497) =    -0.2731397346E+10
       coef(  498) =     0.1502279826E+09
       coef(  499) =    -0.1474855926E+10
       coef(  500) =    -0.1228934893E+10
       coef(  501) =     0.8224385587E+08
       coef(  502) =    -0.1026140869E+10
       coef(  503) =    -0.4125445031E+09
       coef(  504) =    -0.4167359630E+08
       coef(  505) =    -0.3278633336E+03
       coef(  506) =     0.1106135700E+05
       coef(  507) =    -0.2641167358E+06
       coef(  508) =     0.3407751651E+07
       coef(  509) =    -0.8536125849E+07
       coef(  510) =     0.1805747120E+06
       coef(  511) =     0.3235212111E+08
       coef(  512) =    -0.1246897077E+09
       coef(  513) =     0.4267378852E+09
       coef(  514) =    -0.4178284282E+09
       coef(  515) =     0.1935613691E+10
       coef(  516) =    -0.3406322594E+10
       coef(  517) =     0.5463219826E+09
       coef(  518) =    -0.4641559701E+08
       coef(  519) =     0.5796631957E+09
       coef(  520) =    -0.3469952393E+09
       coef(  521) =    -0.1047659296E+10
       coef(  522) =    -0.2645252919E+09
       coef(  523) =    -0.1035527775E+10
       coef(  524) =     0.4496633222E+09
       coef(  525) =     0.4720375106E+09
       coef(  526) =    -0.8919239242E+10
       coef(  527) =     0.7994709429E+10
       coef(  528) =     0.4787348668E+10
       coef(  529) =    -0.3515565584E+05
       coef(  530) =    -0.2459393266E+06
       coef(  531) =    -0.3909283819E+08
       coef(  532) =     0.6040821712E+08
       coef(  533) =    -0.3220851221E+09
       coef(  534) =     0.1821374280E+08
       coef(  535) =    -0.5112729687E+09
       coef(  536) =     0.1388271233E+10
       coef(  537) =     0.3551641990E+10
       coef(  538) =     0.1444719190E+10
       coef(  539) =     0.2583649939E+09
       coef(  540) =    -0.1075371678E+10
       coef(  541) =     0.1237022182E+11
       coef(  542) =    -0.3278006756E+11
       coef(  543) =    -0.1318040104E+11
       coef(  544) =     0.8385510387E+10
       coef(  545) =    -0.7210984031E+10
       coef(  546) =     0.1923706332E+07
       coef(  547) =     0.9860339339E+07
       coef(  548) =     0.9355314856E+09
       coef(  549) =    -0.2174052197E+10
       coef(  550) =     0.1366759154E+10
       coef(  551) =     0.1968800337E+10
       coef(  552) =    -0.1646177665E+11
       coef(  553) =     0.7392509680E+10
       coef(  554) =     0.6632012274E+10
       coef(  555) =    -0.1598584608E+11
       coef(  556) =     0.8522309282E+11
       coef(  557) =    -0.8370766986E+10
       coef(  558) =    -0.2374292350E+08
       coef(  559) =    -0.1354223590E+10
       coef(  560) =     0.3323086706E+10
       coef(  561) =     0.1296195646E+10
       coef(  562) =     0.1200067408E+11
       coef(  563) =    -0.2463406451E+11
       coef(  564) =    -0.1561846871E+11
       coef(  565) =     0.2593206401E+09
       coef(  566) =    -0.2585575568E+10
       coef(  567) =    -0.7840901887E+09
       coef(  568) =     0.6837816481E+10
       coef(  569) =    -0.1209688404E+07
       coef(  570) =     0.8270973162E+08
       coef(  571) =     0.3036530674E+08
       coef(  572) =    -0.6891046776E+09
       coef(  573) =     0.1186382013E+09
       coef(  574) =    -0.8937867130E+09
       coef(  575) =     0.2772490335E+10
       coef(  576) =     0.1211925765E+11
       coef(  577) =    -0.2132063093E+11
       coef(  578) =     0.3497682182E+10
       coef(  579) =     0.1339702987E+11
       coef(  580) =     0.5198430049E+10
       coef(  581) =    -0.2798533262E+08
       coef(  582) =    -0.1091503286E+10
       coef(  583) =     0.1495081765E+10
       coef(  584) =    -0.2727788762E+09
       coef(  585) =     0.5495869301E+10
       coef(  586) =     0.9074279586E+10
       coef(  587) =    -0.6265607372E+11
       coef(  588) =     0.5577739194E+09
       coef(  589) =     0.8966880271E+09
       coef(  590) =    -0.1251374687E+11
       coef(  591) =     0.4135726905E+10
       coef(  592) =    -0.8664881784E+09
       coef(  593) =     0.4955863263E+10
       coef(  594) =     0.2174286054E+09
       coef(  595) =    -0.3294767036E+09
       coef(  596) =    -0.9557479281E+10
       coef(  597) =     0.3838163665E+10
       coef(  598) =    -0.2549729134E+10
       coef(  599) =     0.2616942195E+11
       coef(  600) =    -0.2516840359E+10
       coef(  601) =    -0.1408244406E+09
       coef(  602) =    -0.1109546406E+05
       coef(  603) =    -0.1893589681E+07
       coef(  604) =     0.3607387363E+08
       coef(  605) =    -0.2793761165E+09
       coef(  606) =     0.3019060406E+09
       coef(  607) =     0.3816370251E+08
       coef(  608) =     0.5049294853E+09
       coef(  609) =    -0.3370223748E+10
       coef(  610) =    -0.1943733982E+10
       coef(  611) =    -0.3763262491E+10
       coef(  612) =     0.1123399795E+11
       coef(  613) =    -0.1198454596E+10
       coef(  614) =     0.7986616937E+10
       coef(  615) =     0.1862669332E+11
       coef(  616) =    -0.1590706031E+11
       coef(  617) =     0.9601868139E+10
       coef(  618) =    -0.5302332543E+10
       coef(  619) =     0.7502464827E+06
       coef(  620) =     0.4502763474E+08
       coef(  621) =    -0.3648319975E+09
       coef(  622) =     0.7987842656E+10
       coef(  623) =    -0.9111677770E+10
       coef(  624) =     0.2256271509E+10
       coef(  625) =    -0.1217691132E+11
       coef(  626) =     0.3570737730E+11
       coef(  627) =    -0.3454234309E+10
       coef(  628) =    -0.3196544427E+11
       coef(  629) =     0.4141329111E+11
       coef(  630) =     0.2877214800E+11
       coef(  631) =    -0.2733407453E+08
       coef(  632) =    -0.4170948407E+10
       coef(  633) =     0.2878078483E+10
       coef(  634) =    -0.2980872398E+10
       coef(  635) =     0.5162990057E+11
       coef(  636) =    -0.8657121076E+11
       coef(  637) =    -0.6660602908E+11
       coef(  638) =     0.1419778374E+10
       coef(  639) =    -0.1922925970E+11
       coef(  640) =     0.2272785672E+11
       coef(  641) =     0.3942992458E+11
       coef(  642) =     0.1731270797E+10
       coef(  643) =    -0.1165228645E+11
       coef(  644) =    -0.7010375458E+08
       coef(  645) =     0.1975262465E+10
       coef(  646) =    -0.1948374237E+10
       coef(  647) =    -0.1477723716E+11
       coef(  648) =    -0.2099302974E+11
       coef(  649) =     0.1872558302E+11
       coef(  650) =     0.9305056423E+10
       coef(  651) =     0.7399693789E+09
       coef(  652) =     0.1865588993E+10
       coef(  653) =    -0.1947865913E+11
       coef(  654) =     0.6433591432E+11
       coef(  655) =    -0.5343000201E+10
       coef(  656) =     0.1216787077E+10
       coef(  657) =     0.1248183950E+11
       coef(  658) =     0.5576467334E+08
       coef(  659) =     0.1710505410E+11
       coef(  660) =    -0.1333840959E+11
       coef(  661) =     0.6415655218E+06
       coef(  662) =    -0.4421061048E+08
       coef(  663) =    -0.7447342817E+08
       coef(  664) =     0.8453109516E+09
       coef(  665) =    -0.4129020574E+09
       coef(  666) =     0.8731392170E+09
       coef(  667) =    -0.6775164871E+09
       coef(  668) =    -0.1442160914E+11
       coef(  669) =     0.9522280202E+10
       coef(  670) =    -0.3636669014E+10
       coef(  671) =     0.3646789276E+11
       coef(  672) =    -0.6040782971E+10
       coef(  673) =    -0.1007741764E+08
       coef(  674) =     0.2665879356E+09
       coef(  675) =    -0.5550867665E+10
       coef(  676) =     0.8549275785E+10
       coef(  677) =     0.4455635959E+09
       coef(  678) =    -0.5070568941E+11
       coef(  679) =    -0.3986803932E+11
       coef(  680) =     0.4155149937E+09
       coef(  681) =    -0.3012369268E+10
       coef(  682) =    -0.3502284793E+10
       coef(  683) =     0.7240407241E+11
       coef(  684) =    -0.5912764147E+09
       coef(  685) =    -0.1928208382E+11
       coef(  686) =     0.4521983322E+10
       coef(  687) =    -0.2131576779E+09
       coef(  688) =     0.4398090970E+09
       coef(  689) =     0.1791566807E+11
       coef(  690) =     0.2004640327E+11
       coef(  691) =     0.1742034626E+10
       coef(  692) =    -0.1859215878E+11
       coef(  693) =    -0.9835316834E+10
       coef(  694) =     0.7147385773E+10
       coef(  695) =     0.5078016017E+05
       coef(  696) =     0.2214955579E+08
       coef(  697) =     0.2138944536E+09
       coef(  698) =     0.6909774519E+09
       coef(  699) =    -0.2158812305E+10
       coef(  700) =     0.2028317728E+10
       coef(  701) =    -0.5797053608E+10
       coef(  702) =    -0.1044179930E+09
       coef(  703) =     0.2437342432E+10
       coef(  704) =     0.5501428531E+10
       coef(  705) =     0.2763584706E+11
       coef(  706) =    -0.6778567888E+09
       coef(  707) =    -0.2406934560E+11
       coef(  708) =     0.7614142628E+10
       coef(  709) =     0.2549426200E+09
       coef(  710) =    -0.1237261408E+11
       coef(  711) =     0.4146577574E+10
       coef(  712) =    -0.5598958234E+07
       coef(  713) =     0.2476858053E+09
       coef(  714) =     0.1402516916E+09
       coef(  715) =    -0.2351656581E+10
       coef(  716) =     0.1998912959E+08
       coef(  717) =    -0.9242162228E+09
       coef(  718) =     0.9701941234E+09
       coef(  719) =     0.1103164761E+10
       coef(  720) =     0.2907753176E+06
       coef(  721) =    -0.2272700591E+08
       coef(  722) =    -0.2325521316E+09
       coef(  723) =     0.9130954497E+09
       coef(  724) =    -0.3041356760E+09
       coef(  725) =    -0.6990688645E+08
       coef(  726) =     0.5783325598E+10
       coef(  727) =    -0.4249321892E+11
       coef(  728) =     0.1702696999E+11
       coef(  729) =    -0.6140113320E+08
       coef(  730) =     0.5574028293E+11
       coef(  731) =    -0.4168399148E+10
       coef(  732) =    -0.2665361098E+08
       coef(  733) =     0.1334679870E+10
       coef(  734) =    -0.8602211279E+10
       coef(  735) =     0.2697678745E+11
       coef(  736) =    -0.1011780397E+11
       coef(  737) =    -0.4104515103E+11
       coef(  738) =    -0.8152566041E+11
       coef(  739) =     0.9879306847E+08
       coef(  740) =     0.3384778954E+10
       coef(  741) =    -0.2805972908E+11
       coef(  742) =     0.7450003955E+11
       coef(  743) =    -0.4148621195E+09
       coef(  744) =    -0.1588461987E+11
       coef(  745) =     0.4501618316E+10
       coef(  746) =    -0.4294023775E+09
       coef(  747) =     0.8357522421E+10
       coef(  748) =    -0.6810781926E+10
       coef(  749) =     0.5207705448E+10
       coef(  750) =     0.2562825891E+09
       coef(  751) =     0.2474955859E+11
       coef(  752) =    -0.1718118888E+11
       coef(  753) =     0.3586678152E+10
       coef(  754) =     0.1479168440E+08
       coef(  755) =    -0.2832495368E+09
       coef(  756) =     0.1169850617E+10
       coef(  757) =     0.4042018524E+10
       coef(  758) =    -0.3131552717E+08
       coef(  759) =    -0.2584938849E+11
       coef(  760) =    -0.3230242031E+11
       coef(  761) =    -0.2667905243E+09
       coef(  762) =     0.2078224040E+10
       coef(  763) =     0.1373458810E+10
       coef(  764) =     0.2409688805E+12
       coef(  765) =    -0.6545106678E+09
       coef(  766) =    -0.1028400922E+12
       coef(  767) =     0.2480948539E+11
       coef(  768) =     0.1065880587E+10
       coef(  769) =    -0.4118206330E+11
       coef(  770) =    -0.1727501422E+11
       coef(  771) =     0.6443709448E+07
       coef(  772) =     0.2283111493E+10
       coef(  773) =     0.1179138861E+09
       coef(  774) =    -0.1201230303E+11
       coef(  775) =    -0.1137382292E+09
       coef(  776) =    -0.2176040723E+11
       coef(  777) =     0.1479107667E+10
       coef(  778) =     0.7713970891E+10
       coef(  779) =    -0.2295299205E+09
       coef(  780) =    -0.3320178323E+09
       coef(  781) =     0.1190766112E+10
       coef(  782) =     0.2997612156E+08
       coef(  783) =     0.1504164630E+10
       coef(  784) =    -0.7734339658E+09
       coef(  785) =    -0.2093760632E+11
       coef(  786) =     0.2453892561E+08
       coef(  787) =    -0.8578769305E+10
       coef(  788) =    -0.1067615865E+10
       coef(  789) =     0.1981346154E+10
       coef(  790) =    -0.6995476288E+09
       coef(  791) =    -0.6252652260E+09
       coef(  792) =     0.5496713181E+10
       coef(  793) =     0.6884491748E+09
       coef(  794) =     0.6168519181E+09
       coef(  795) =     0.5602544339E+04
       coef(  796) =     0.1049932064E+07
       coef(  797) =    -0.2232867279E+08
       coef(  798) =     0.1482501392E+09
       coef(  799) =    -0.9272735175E+08
       coef(  800) =    -0.3103198234E+08
       coef(  801) =     0.4260338778E+08
       coef(  802) =    -0.1432001508E+09
       coef(  803) =     0.3176521223E+09
       coef(  804) =     0.3063227231E+10
       coef(  805) =    -0.7672570110E+10
       coef(  806) =     0.7555480591E+09
       coef(  807) =    -0.5124995170E+10
       coef(  808) =    -0.2449959201E+10
       coef(  809) =     0.1269712111E+11
       coef(  810) =    -0.6148628469E+10
       coef(  811) =     0.1047376735E+11
       coef(  812) =    -0.6148687977E+06
       coef(  813) =     0.1085317971E+08
       coef(  814) =    -0.9130444260E+07
       coef(  815) =    -0.2654938397E+10
       coef(  816) =     0.4685391050E+10
       coef(  817) =    -0.1242184434E+10
       coef(  818) =     0.8670499712E+10
       coef(  819) =    -0.1982868140E+11
       coef(  820) =     0.3065928786E+09
       coef(  821) =     0.1715491580E+11
       coef(  822) =    -0.4799663456E+11
       coef(  823) =    -0.1215149051E+11
       coef(  824) =     0.6025081929E+07
       coef(  825) =     0.1580187562E+10
       coef(  826) =    -0.3791633530E+10
       coef(  827) =     0.2472798837E+10
       coef(  828) =    -0.2518226558E+11
       coef(  829) =     0.5095270524E+11
       coef(  830) =     0.3991832995E+11
       coef(  831) =    -0.4923560376E+09
       coef(  832) =     0.1064338952E+11
       coef(  833) =    -0.9182452891E+10
       coef(  834) =    -0.2230509989E+11
       coef(  835) =    -0.1025088716E+10
       coef(  836) =     0.3966130203E+10
       coef(  837) =     0.3197063379E+08
       coef(  838) =    -0.1041797288E+10
       coef(  839) =     0.4212932643E+09
       coef(  840) =     0.7268165927E+10
       coef(  841) =     0.9839740637E+10
       coef(  842) =    -0.1125803120E+11
       coef(  843) =    -0.2145596135E+11
       coef(  844) =    -0.1699577896E+09
       coef(  845) =     0.8311455099E+09
       coef(  846) =     0.1356629759E+11
       coef(  847) =    -0.4769947300E+10
       coef(  848) =     0.1095332311E+10
       coef(  849) =    -0.2020938303E+11
       coef(  850) =    -0.1722047177E+10
       coef(  851) =     0.3932572397E+08
       coef(  852) =    -0.1070175342E+11
       coef(  853) =     0.6728230649E+10
       coef(  854) =    -0.6911111431E+06
       coef(  855) =     0.7341840385E+08
       coef(  856) =     0.3072882459E+08
       coef(  857) =    -0.1160446628E+10
       coef(  858) =     0.2442601350E+10
       coef(  859) =    -0.1831464709E+10
       coef(  860) =     0.2850029550E+10
       coef(  861) =     0.4304763794E+11
       coef(  862) =    -0.2550686307E+11
       coef(  863) =     0.9862842082E+10
       coef(  864) =    -0.1238678181E+12
       coef(  865) =    -0.1873438391E+11
       coef(  866) =     0.2239381875E+08
       coef(  867) =    -0.6740899573E+09
       coef(  868) =     0.8055415541E+10
       coef(  869) =    -0.4214641988E+11
       coef(  870) =    -0.5648250984E+10
       coef(  871) =     0.1206216182E+12
       coef(  872) =     0.1529310693E+12
       coef(  873) =    -0.5276788027E+09
       coef(  874) =     0.1544049182E+11
       coef(  875) =     0.2730495209E+11
       coef(  876) =    -0.2026389164E+12
       coef(  877) =    -0.3270575397E+10
       coef(  878) =     0.5124227781E+11
       coef(  879) =    -0.9646639267E+10
       coef(  880) =     0.5654910545E+09
       coef(  881) =    -0.7256904467E+10
       coef(  882) =     0.2960026562E+10
       coef(  883) =     0.1113920574E+11
       coef(  884) =    -0.4191811663E+10
       coef(  885) =    -0.3193838016E+11
       coef(  886) =     0.3360309824E+11
       coef(  887) =    -0.9402480958E+10
       coef(  888) =    -0.5516211905E+07
       coef(  889) =     0.2964701373E+09
       coef(  890) =    -0.1135462920E+09
       coef(  891) =    -0.5143927398E+10
       coef(  892) =     0.4965924132E+09
       coef(  893) =     0.1053212251E+11
       coef(  894) =     0.2837905904E+11
       coef(  895) =     0.1628579664E+09
       coef(  896) =    -0.2306557834E+10
       coef(  897) =     0.5762979691E+10
       coef(  898) =    -0.1169897273E+12
       coef(  899) =    -0.1849580480E+08
       coef(  900) =     0.6983359293E+11
       coef(  901) =    -0.1742765937E+11
       coef(  902) =    -0.6694632169E+09
       coef(  903) =     0.1828646669E+11
       coef(  904) =     0.2949690532E+10
       coef(  905) =     0.1183828378E+08
       coef(  906) =    -0.1406960215E+10
       coef(  907) =    -0.2348350374E+10
       coef(  908) =     0.1217191779E+11
       coef(  909) =     0.3311419404E+09
       coef(  910) =     0.9469317304E+09
       coef(  911) =    -0.1884584583E+10
       coef(  912) =    -0.3134426833E+10
       coef(  913) =     0.6455623686E+08
       coef(  914) =     0.2506456484E+09
       coef(  915) =    -0.9005011268E+09
       coef(  916) =    -0.1154170447E+08
       coef(  917) =     0.5322702248E+09
       coef(  918) =    -0.5963986589E+09
       coef(  919) =    -0.3031108964E+10
       coef(  920) =    -0.1309779577E+10
       coef(  921) =     0.2007254497E+11
       coef(  922) =     0.2387204594E+11
       coef(  923) =     0.2575935599E+09
       coef(  924) =    -0.4705820984E+10
       coef(  925) =     0.2489845007E+11
       coef(  926) =    -0.1469779526E+12
       coef(  927) =     0.1114785445E+09
       coef(  928) =     0.6338158364E+11
       coef(  929) =    -0.1764725019E+11
       coef(  930) =    -0.7264492192E+09
       coef(  931) =     0.6332554235E+10
       coef(  932) =     0.1573482032E+11
       coef(  933) =    -0.1316367000E+09
       coef(  934) =    -0.5462817476E+10
       coef(  935) =    -0.3317220657E+10
       coef(  936) =     0.3549350193E+11
       coef(  937) =     0.1024028948E+10
       coef(  938) =     0.1456600799E+10
       coef(  939) =     0.8510137724E+10
       coef(  940) =    -0.1560908289E+10
       coef(  941) =     0.7319136604E+09
       coef(  942) =     0.2967952610E+10
       coef(  943) =    -0.2953820323E+10
       coef(  944) =    -0.9103383396E+08
       coef(  945) =     0.6979053886E+09
       coef(  946) =     0.7568354006E+10
       coef(  947) =    -0.5565934818E+10
       coef(  948) =    -0.2892280717E+10
       coef(  949) =     0.1494379628E+06
       coef(  950) =    -0.2581527227E+08
       coef(  951) =     0.8617308410E+08
       coef(  952) =     0.1769411741E+09
       coef(  953) =    -0.1432468938E+10
       coef(  954) =     0.8493938742E+09
       coef(  955) =    -0.4843452540E+10
       coef(  956) =    -0.6160525112E+10
       coef(  957) =     0.7005767387E+10
       coef(  958) =    -0.4987043315E+10
       coef(  959) =     0.5086824168E+11
       coef(  960) =     0.1693612668E+11
       coef(  961) =    -0.3236508320E+07
       coef(  962) =    -0.2686676056E+08
       coef(  963) =     0.1645672653E+09
       coef(  964) =     0.1285084790E+11
       coef(  965) =     0.5213563476E+10
       coef(  966) =    -0.5225261529E+11
       coef(  967) =    -0.5879384093E+11
       coef(  968) =     0.1658753918E+09
       coef(  969) =    -0.7513306334E+10
       coef(  970) =    -0.6301336897E+10
       coef(  971) =     0.7928823425E+11
       coef(  972) =     0.1542815584E+10
       coef(  973) =    -0.2088908374E+11
       coef(  974) =     0.3538606561E+10
       coef(  975) =    -0.1943050238E+09
       coef(  976) =     0.2511662280E+10
       coef(  977) =    -0.5163989451E+10
       coef(  978) =    -0.9054132653E+10
       coef(  979) =     0.1476957324E+10
       coef(  980) =     0.1407381705E+11
       coef(  981) =    -0.9152946072E+10
       coef(  982) =     0.2383145578E+10
       coef(  983) =     0.7553536726E+07
       coef(  984) =    -0.6476220590E+09
       coef(  985) =     0.1936762792E+09
       coef(  986) =     0.6192840761E+10
       coef(  987) =     0.5487540057E+10
       coef(  988) =    -0.2580969440E+11
       coef(  989) =    -0.2629538979E+11
       coef(  990) =    -0.1718891676E+09
       coef(  991) =     0.1443057208E+10
       coef(  992) =    -0.3145122989E+11
       coef(  993) =     0.1254653834E+12
       coef(  994) =     0.1465219069E+10
       coef(  995) =    -0.7686070493E+11
       coef(  996) =     0.1961415026E+11
       coef(  997) =     0.8970352227E+09
       coef(  998) =    -0.1540079395E+10
       coef(  999) =    -0.7489138355E+10
       coef( 1000) =     0.2096955401E+08
       coef( 1001) =     0.1741044788E+10
       coef( 1002) =     0.5010090166E+10
       coef( 1003) =    -0.2196033372E+11
       coef( 1004) =    -0.7073469175E+09
       coef( 1005) =     0.8276723430E+10
       coef( 1006) =    -0.1807509492E+10
       coef( 1007) =     0.5037026095E+09
       coef( 1008) =    -0.1910611383E+09
       coef( 1009) =    -0.9904888259E+08
       coef( 1010) =     0.2035862409E+10
       coef( 1011) =    -0.2061877144E+08
       coef( 1012) =     0.7526474668E+08
       coef( 1013) =     0.2225611971E+10
       coef( 1014) =     0.9116073986E+09
       coef( 1015) =    -0.1553990250E+11
       coef( 1016) =    -0.8068448800E+09
       coef( 1017) =     0.1547249348E+11
       coef( 1018) =    -0.6851358262E+10
       coef( 1019) =    -0.1748436767E+10
       coef( 1020) =    -0.6351555941E+09
       coef( 1021) =    -0.7835735383E+10
       coef( 1022) =     0.3237852358E+10
       coef( 1023) =     0.4051621054E+09
       coef( 1024) =     0.1432191850E+10
       coef( 1025) =    -0.8266795831E+06
       coef( 1026) =     0.1496940144E+09
       coef( 1027) =    -0.1094936939E+09
       coef( 1028) =    -0.1533653569E+10
       coef( 1029) =    -0.2189621873E+10
       coef( 1030) =     0.1109060122E+11
       coef( 1031) =     0.1475846421E+10
       coef( 1032) =     0.1861505803E+08
       coef( 1033) =     0.2279849825E+09
       coef( 1034) =     0.7189082511E+10
       coef( 1035) =    -0.2395224518E+11
       coef( 1036) =    -0.6253276857E+09
       coef( 1037) =     0.1803234116E+11
       coef( 1038) =    -0.4394718890E+10
       coef( 1039) =    -0.1428495025E+09
       coef( 1040) =    -0.4429106387E+09
       coef( 1041) =    -0.1390877307E+09
       coef( 1042) =    -0.2283096679E+08
       coef( 1043) =    -0.7278316205E+09
       coef( 1044) =    -0.1927994826E+10
       coef( 1045) =     0.1082373505E+11
       coef( 1046) =     0.5204863439E+09
       coef( 1047) =    -0.8969004027E+10
       coef( 1048) =     0.2803865821E+10
       coef( 1049) =     0.7118086827E+09
       coef( 1050) =     0.1387349712E+09
       coef( 1051) =     0.5203800559E+09
       coef( 1052) =    -0.1521838633E+10
       coef( 1053) =     0.4310195737E+07
       coef( 1054) =    -0.1085137265E+08
       coef( 1055) =     0.9595174481E+09
       coef( 1056) =     0.2384338918E+09
       coef( 1057) =    -0.2640816535E+09

      r1=r(1)*autoang
      r2=r(2)*autoang
      r3=r(3)*autoang
      r4=r(4)*autoang
      r5=r(5)*autoang
      r6=r(6)*autoang
      r7=r(7)*autoang
      r8=r(8)*autoang
      r9=r(9)*autoang

      b1=1.d0

      e1=dexp(-r1/b1)
      e2=dexp(-r2/b1)
      e3=dexp(-r3/b1)
      e4=dexp(-r4/b1)
      e5=dexp(-r5/b1)
      e6=dexp(-r6/b1)
      e7=dexp(-r7/b1)
      e8=dexp(-r8/b1)
      e9=dexp(-r9/b1)

      norder=4
      norder2=4
      norder3=4
      norder4=4
      maxorder=7

      ii=0
      e=0.d0
      do i1=0,norder
      x1=e1**i1
      do i2=0,norder2
      do i3=i2,norder2
      x2=x1*(e2**i2*e3**i3+e3**i2*e2**i3)
      do i4=0,norder3
      do i5=i4,norder3
      x3=x2*(e4**i4*e5**i5+e5**i4*e6**i5)
      do i6=0,norder4
      do i7=i6,norder4
      do i8=i7,norder4
      do i9=i8,norder4
      if ((i1+i2+i3+i4+i5+i6+i7+i8+i9).le.maxorder) then
      x4=x3*(
     & e6**i6*e7**i7*e8**i8*e9**i9 +
     & e6**i6*e7**i7*e9**i8*e8**i9 +
     & e6**i6*e8**i7*e9**i8*e7**i9 +
     & e6**i6*e8**i7*e7**i8*e9**i9 +
     & e6**i6*e9**i7*e8**i8*e7**i9 +
     & e6**i6*e9**i7*e7**i8*e8**i9 +
     & e7**i6*e8**i7*e9**i8*e6**i9 +
     & e7**i6*e8**i7*e6**i8*e9**i9 +
     & e7**i6*e9**i7*e6**i8*e8**i9 +
     & e7**i6*e9**i7*e8**i8*e6**i9 +
     & e7**i6*e6**i7*e9**i8*e8**i9 +
     & e7**i6*e6**i7*e8**i8*e9**i9 +
     & e8**i6*e9**i7*e6**i8*e7**i9 +
     & e8**i6*e9**i7*e7**i8*e6**i9 +
     & e8**i6*e6**i7*e7**i8*e9**i9 +
     & e8**i6*e6**i7*e9**i8*e7**i9 +
     & e8**i6*e7**i7*e6**i8*e9**i9 +
     & e8**i6*e7**i7*e9**i8*e6**i9 +
     & e9**i6*e6**i7*e7**i8*e8**i9 +
     & e9**i6*e6**i7*e8**i8*e7**i9 +
     & e9**i6*e7**i7*e8**i8*e6**i9 +
     & e9**i6*e7**i7*e6**i8*e8**i9 +
     & e9**i6*e8**i7*e7**i8*e6**i9 +
     & e9**i6*e8**i7*e6**i8*e7**i9 )
      ii=ii+1
      basis(ii)=1.d0
      if (i1+i2+i3+i4+i5+i6+i7+i8+i9.ne.0) basis(ii)=x4
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      v=0.d0
      do j=1,ncoef
        v=v+coef(j)*basis(j)
      enddo
      v=v/autocmi

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

