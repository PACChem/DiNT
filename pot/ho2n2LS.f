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
        call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
c        call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
        call lsbath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
      endif
      if (natom2.ne.0) then
        call ho2(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
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

      dimension r(6)

      ii=0
      do i=1,3
      do j=4,5
      ii=ii+1
      ij=ii
      r(ij)=dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
      enddo
      enddo

      call lsbathfit(r,v)

      resp=0.0001d0

      ii=0
      do i=1,3
      do j=4,5
      ii=ii+1
      ij=ii
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

      dimension r(6)
      dimension coef(3000),basis(3000)

       ncoef=  982 
       coef(    1) =    -0.20146179700128095114E+01
       coef(    2) =    -0.13061509849402082750E+03
       coef(    3) =     0.60724971336275622889E+05
       coef(    4) =    -0.77022876770858420059E+06
       coef(    5) =    -0.11142799125094839931E+09
       coef(    6) =     0.97358296776296305656E+09
       coef(    7) =    -0.19755632086051759720E+10
       coef(    8) =    -0.33673021067852561828E+05
       coef(    9) =    -0.18332042128086980432E+08
       coef(   10) =     0.11859079107210268974E+10
       coef(   11) =    -0.23953383400837860107E+10
       coef(   12) =    -0.44846277972343223572E+11
       coef(   13) =     0.16385182060432989502E+12
       coef(   14) =     0.42550698881224989891E+09
       coef(   15) =    -0.37563752899434127808E+11
       coef(   16) =     0.33471425504882470703E+12
       coef(   17) =    -0.10054040327441966553E+13
       coef(   18) =     0.34268326179568634033E+11
       coef(   19) =     0.32988927300896997070E+12
       coef(   20) =    -0.54559904712178457031E+13
       coef(   21) =     0.20391908750036589844E+14
       coef(   22) =    -0.17355631735940582031E+14
       coef(   23) =     0.19906269580103097656E+14
       coef(   24) =    -0.13382110882803243750E+15
       coef(   25) =     0.10649832553953085938E+15
       coef(   26) =     0.17173260482448796875E+15
       coef(   27) =    -0.91986599938611093750E+14
       coef(   28) =    -0.25565222414047156250E+15
       coef(   29) =     0.12773118758646558970E+08
       coef(   30) =    -0.13846268642974917889E+10
       coef(   31) =    -0.25559383173413452148E+11
       coef(   32) =     0.12146250532389486694E+12
       coef(   33) =     0.58703574146735876465E+12
       coef(   34) =    -0.25310232937173281250E+13
       coef(   35) =     0.62702418644818656921E+11
       coef(   36) =    -0.13510800682336973572E+12
       coef(   37) =    -0.72092857364349658203E+12
       coef(   38) =    -0.34369735868971845703E+13
       coef(   39) =     0.23015310218490242188E+14
       coef(   40) =    -0.11916638998874863281E+13
       coef(   41) =     0.52531623242009667969E+13
       coef(   42) =    -0.22909261398958761719E+14
       coef(   43) =    -0.17465129539404058594E+14
       coef(   44) =     0.29536698218469312500E+14
       coef(   45) =    -0.19089967285079808594E+14
       coef(   46) =    -0.75403088576383515625E+14
       coef(   47) =    -0.94391936555578125000E+13
       coef(   48) =    -0.39205574704226765625E+14
       coef(   49) =    -0.63564198727427026367E+12
       coef(   50) =     0.82384442689200878906E+13
       coef(   51) =    -0.18606438121877265625E+14
       coef(   52) =     0.26166860127591023438E+14
       coef(   53) =    -0.20806268440002566406E+14
       coef(   54) =    -0.22709429509690058594E+14
       coef(   55) =     0.59572298855821320312E+14
       coef(   56) =     0.92024306829175742188E+13
       coef(   57) =    -0.53792975242393187500E+14
       coef(   58) =    -0.58758664076709746094E+13
       coef(   59) =    -0.31862692545826390625E+14
       coef(   60) =    -0.39324053923625984375E+14
       coef(   61) =    -0.20877703916830304688E+14
       coef(   62) =     0.95759147322813300781E+13
       coef(   63) =    -0.26409469656063226562E+14
       coef(   64) =    -0.17501721660017707031E+14
       coef(   65) =    -0.25498552171798515625E+14
       coef(   66) =    -0.40060808741686367188E+14
       coef(   67) =    -0.17930581235077984375E+14
       coef(   68) =    -0.25615914553895476562E+14
       coef(   69) =     0.27350350648546826839E+09
       coef(   70) =    -0.10586692811103097916E+11
       coef(   71) =     0.42758446578491290283E+12
       coef(   72) =    -0.37853379356960878906E+13
       coef(   73) =     0.91323217475211230469E+13
       coef(   74) =    -0.42767270707607158203E+13
       coef(   75) =    -0.13721347170923989868E+12
       coef(   76) =    -0.29784732036332221680E+13
       coef(   77) =     0.42826805143450171875E+14
       coef(   78) =    -0.10758186456708329688E+15
       coef(   79) =     0.44860714760481218750E+14
       coef(   80) =     0.18993262398668554688E+14
       coef(   81) =    -0.14230883767609306250E+15
       coef(   82) =     0.19973152780698181250E+15
       coef(   83) =     0.17348090783466153125E+15
       coef(   84) =    -0.33252297560537460938E+14
       coef(   85) =     0.76638014289246562500E+14
       coef(   86) =     0.52772744895886335938E+14
       coef(   87) =     0.36053542553857898438E+14
       coef(   88) =     0.25919967838967250977E+13
       coef(   89) =    -0.43794865868901820312E+14
       coef(   90) =     0.50935123109404203125E+14
       coef(   91) =    -0.41376122729321968750E+14
       coef(   92) =    -0.39188091326890523438E+14
       coef(   93) =     0.87882771697059468750E+14
       coef(   94) =     0.12680004327779601562E+14
       coef(   95) =     0.17666350173764339844E+14
       coef(   96) =     0.63683397141833242188E+13
       coef(   97) =    -0.32364352784650761719E+13
       coef(   98) =     0.55075162968557128906E+13
       coef(   99) =     0.70963169135382656250E+13
       coef(  100) =    -0.18426343194536390625E+14
       coef(  101) =    -0.45644089115333417969E+13
       coef(  102) =    -0.12488733538253632812E+14
       coef(  103) =     0.48657172391970458984E+13
       coef(  104) =    -0.64692699791427988281E+13
       coef(  105) =     0.24828242076866695312E+14
       coef(  106) =    -0.17016986845003158203E+14
       coef(  107) =    -0.22677890992076812500E+14
       coef(  108) =    -0.56609313949507460938E+14
       coef(  109) =    -0.17467567823592583984E+14
       coef(  110) =    -0.49676047336500683594E+13
       coef(  111) =    -0.64038621902850244141E+13
       coef(  112) =    -0.48966545842702421875E+14
       coef(  113) =    -0.14651874877433568359E+14
       coef(  114) =    -0.23143150214045117188E+14
       coef(  115) =     0.40700095763392505432E+03
       coef(  116) =    -0.51750957297988439677E+05
       coef(  117) =     0.40524161068432601169E+07
       coef(  118) =    -0.43992659283050276339E+08
       coef(  119) =     0.69648456236928844452E+09
       coef(  120) =     0.84238576052194280624E+10
       coef(  121) =    -0.49991255525249290466E+11
       coef(  122) =     0.61226407352060936391E+07
       coef(  123) =    -0.93998292908547604084E+09
       coef(  124) =     0.66935938751908283234E+10
       coef(  125) =    -0.13183055041971403503E+12
       coef(  126) =     0.43904376255594702148E+12
       coef(  127) =     0.91172518387683502197E+11
       coef(  128) =     0.16270993915296646118E+11
       coef(  129) =    -0.32801571682120013428E+12
       coef(  130) =     0.19678910577193081055E+13
       coef(  131) =    -0.55318995050203964844E+13
       coef(  132) =     0.10257089120728695312E+14
       coef(  133) =     0.23980382098060673828E+13
       coef(  134) =    -0.22756918237874167969E+14
       coef(  135) =     0.31431831408764328125E+14
       coef(  136) =    -0.45498824318133828125E+14
       coef(  137) =     0.47321153381138781250E+14
       coef(  138) =    -0.75855642313741437500E+14
       coef(  139) =     0.32616293822372351562E+14
       coef(  140) =    -0.14151542651972201172E+14
       coef(  141) =     0.38906937865201140625E+14
       coef(  142) =     0.36047071890339642763E+09
       coef(  143) =     0.98802181656957321167E+10
       coef(  144) =     0.21096372122282687378E+12
       coef(  145) =     0.76886181842581811523E+12
       coef(  146) =    -0.87526833346772177734E+13
       coef(  147) =     0.99412936067054726562E+13
       coef(  148) =    -0.42088174854852069092E+12
       coef(  149) =    -0.24500039842403847656E+13
       coef(  150) =     0.75668506308292949219E+13
       coef(  151) =     0.47533444967199640625E+14
       coef(  152) =    -0.12217076345933210938E+15
       coef(  153) =     0.18510155671798355469E+14
       coef(  154) =    -0.87898350045484765625E+14
       coef(  155) =     0.13332870462791240625E+15
       coef(  156) =     0.97559838663801718750E+14
       coef(  157) =    -0.66139180724298781250E+14
       coef(  158) =     0.55713059582812570312E+14
       coef(  159) =     0.55213609298694828125E+14
       coef(  160) =     0.44705355303685562500E+14
       coef(  161) =     0.49221165101825986328E+13
       coef(  162) =    -0.43573639622179062500E+14
       coef(  163) =     0.22695460175368554688E+14
       coef(  164) =     0.50077878830562304688E+12
       coef(  165) =    -0.46161454058242781250E+14
       coef(  166) =     0.10970482986062881250E+15
       coef(  167) =    -0.72394160250042031250E+14
       coef(  168) =    -0.19735750308972031250E+14
       coef(  169) =    -0.45212299070103486328E+13
       coef(  170) =    -0.64868923514491734375E+14
       coef(  171) =    -0.12354997423289833984E+14
       coef(  172) =    -0.23240096657679562500E+14
       coef(  173) =    -0.70802732350954093750E+14
       coef(  174) =    -0.24717003058147496094E+14
       coef(  175) =    -0.37605697767416562500E+14
       coef(  176) =    -0.12248514808438083649E+11
       coef(  177) =     0.17989989124773901367E+12
       coef(  178) =    -0.24996312560818515625E+13
       coef(  179) =     0.35392136343924628906E+13
       coef(  180) =     0.11623712616962634766E+14
       coef(  181) =    -0.56632431759461738281E+13
       coef(  182) =     0.12835608953992653809E+13
       coef(  183) =     0.20312645740139808594E+14
       coef(  184) =    -0.53212812705977359375E+14
       coef(  185) =    -0.32849598004250507812E+13
       coef(  186) =     0.37479806361965332031E+13
       coef(  187) =    -0.23970608774886785156E+14
       coef(  188) =     0.47992629430429429688E+14
       coef(  189) =     0.82530107936460093750E+14
       coef(  190) =     0.48498072966785929688E+14
       coef(  191) =     0.50594300176709164062E+14
       coef(  192) =     0.41411746067714093750E+14
       coef(  193) =    -0.23691590685199375000E+14
       coef(  194) =     0.46073888108510265625E+14
       coef(  195) =     0.12063317049401675781E+14
       coef(  196) =     0.54017774857461386719E+13
       coef(  197) =     0.97209903865274414062E+11
       coef(  198) =    -0.56364825248964546875E+14
       coef(  199) =    -0.17199969968645470703E+14
       coef(  200) =     0.36363073940007065430E+13
       coef(  201) =    -0.34895681178768935547E+13
       coef(  202) =    -0.59488623744678835938E+14
       coef(  203) =    -0.19540698839905695312E+14
       coef(  204) =     0.50788837568162726562E+14
       coef(  205) =    -0.83877548237595312500E+14
       coef(  206) =    -0.25039174267042984375E+14
       coef(  207) =    -0.47550992354510898438E+13
       coef(  208) =    -0.70140836792093773438E+14
       coef(  209) =    -0.18465165654849648438E+14
       coef(  210) =    -0.32177536200240210938E+14
       coef(  211) =    -0.20870568785086645221E+05
       coef(  212) =     0.12908136290376207326E+07
       coef(  213) =    -0.15438805415493947268E+09
       coef(  214) =     0.88841116809139289856E+10
       coef(  215) =    -0.12516408000118840027E+12
       coef(  216) =     0.43879058222667028809E+12
       coef(  217) =    -0.16980038583336218262E+12
       coef(  218) =     0.13075427700077532232E+09
       coef(  219) =    -0.14104216057565176010E+11
       coef(  220) =     0.24747140147003939819E+12
       coef(  221) =    -0.23487599585835197449E+11
       coef(  222) =    -0.90727787987765795898E+12
       coef(  223) =    -0.68618940769885087891E+13
       coef(  224) =     0.13897575340297192383E+12
       coef(  225) =    -0.19435168823378759766E+13
       coef(  226) =    -0.90401917212784257812E+13
       coef(  227) =     0.30096794605981992188E+14
       coef(  228) =    -0.21861748448534964844E+14
       coef(  229) =     0.72154620360184550781E+13
       coef(  230) =     0.22321529009757054688E+14
       coef(  231) =    -0.47844400589571929688E+14
       coef(  232) =     0.51455840479573453125E+14
       coef(  233) =    -0.53472302035765859375E+14
       coef(  234) =    -0.30542479933022148438E+13
       coef(  235) =     0.53978556993888484375E+14
       coef(  236) =     0.30908444256603191406E+14
       coef(  237) =     0.22195890942028074265E+10
       coef(  238) =    -0.38662014729607623291E+12
       coef(  239) =    -0.28474826340564660645E+12
       coef(  240) =     0.10629070091800365234E+14
       coef(  241) =    -0.25622181278644820312E+14
       coef(  242) =     0.61683101852293531250E+14
       coef(  243) =     0.23294498219244306641E+13
       coef(  244) =    -0.25076806965485322266E+13
       coef(  245) =     0.11856403646159888672E+14
       coef(  246) =    -0.10604245124946101562E+15
       coef(  247) =    -0.51922970024364468750E+14
       coef(  248) =    -0.41614763756075500000E+14
       coef(  249) =     0.69391255676196437500E+14
       coef(  250) =     0.52432664013064257812E+14
       coef(  251) =     0.35725685210595984375E+14
       coef(  252) =     0.70999771423366718750E+14
       coef(  253) =     0.45434002939491500000E+14
       coef(  254) =    -0.88104851002199394531E+13
       coef(  255) =     0.49738541235079773438E+14
       coef(  256) =     0.75713374584968171875E+14
       coef(  257) =     0.58902722717657519531E+13
       coef(  258) =    -0.44654513322871230469E+13
       coef(  259) =    -0.35089522357948113281E+14
       coef(  260) =     0.46492706929049042969E+13
       coef(  261) =     0.77760192661296640625E+13
       coef(  262) =     0.57863802024990605469E+13
       coef(  263) =    -0.54159026559958171875E+14
       coef(  264) =    -0.16137121664939082031E+14
       coef(  265) =     0.26794246101135437012E+12
       coef(  266) =    -0.28915532459473999023E+13
       coef(  267) =     0.61050385459457031250E+13
       coef(  268) =    -0.41628661730108750000E+13
       coef(  269) =    -0.35877423855804328125E+14
       coef(  270) =     0.18482581213894839844E+14
       coef(  271) =     0.15425324973382365234E+14
       coef(  272) =    -0.98786601765243156250E+14
       coef(  273) =     0.67219425710086242188E+14
       coef(  274) =     0.38161769667064429688E+14
       coef(  275) =     0.20692987820589054688E+14
       coef(  276) =    -0.40684118926453835938E+14
       coef(  277) =     0.41914538798456828125E+14
       coef(  278) =     0.31765661130968101562E+14
       coef(  279) =     0.29556058659710539062E+14
       coef(  280) =     0.85784709399619589844E+13
       coef(  281) =     0.10605969679483160938E+15
       coef(  282) =     0.76329243775203078125E+14
       coef(  283) =     0.29585305734628515625E+14
       coef(  284) =     0.39294551013179281250E+14
       coef(  285) =     0.19940680372486300781E+14
       coef(  286) =     0.72611070598007246094E+13
       coef(  287) =    -0.46060088160800921875E+14
       coef(  288) =     0.23024982665638742188E+14
       coef(  289) =     0.20504065718049570312E+14
       coef(  290) =     0.61175509562450605469E+13
       coef(  291) =     0.49011078270960698137E+06
       coef(  292) =    -0.66673007227690786123E+08
       coef(  293) =     0.36731190804481754303E+10
       coef(  294) =    -0.21364226211604513550E+12
       coef(  295) =     0.26490490993773974609E+13
       coef(  296) =    -0.94084184611269316406E+13
       coef(  297) =     0.70001187066273349609E+13
       coef(  298) =     0.13867132637579495907E+10
       coef(  299) =     0.11513504295470982361E+12
       coef(  300) =    -0.21232870594639602051E+13
       coef(  301) =    -0.40303834875790258789E+13
       coef(  302) =     0.15089367443735251953E+14
       coef(  303) =     0.70806803014189406250E+14
       coef(  304) =    -0.22188577727031235352E+13
       coef(  305) =     0.30364596806938734375E+14
       coef(  306) =     0.22278557268924242188E+14
       coef(  307) =    -0.10318999622074041016E+14
       coef(  308) =     0.22019195896871435547E+12
       coef(  309) =    -0.10016232235688189062E+15
       coef(  310) =     0.31178728009056300781E+14
       coef(  311) =    -0.16166577134629433594E+14
       coef(  312) =     0.62066329808103359375E+13
       coef(  313) =     0.66437410262322015625E+14
       coef(  314) =     0.29371463611935691406E+14
       coef(  315) =    -0.30227203859663562775E+11
       coef(  316) =     0.58756947081082285156E+13
       coef(  317) =    -0.22108635655376367188E+14
       coef(  318) =    -0.25886838411892523438E+14
       coef(  319) =    -0.55878357466230136719E+13
       coef(  320) =     0.94891074136672617188E+13
       coef(  321) =    -0.13883969838507807617E+13
       coef(  322) =     0.38401654804563812500E+14
       coef(  323) =    -0.46271847201993718750E+14
       coef(  324) =    -0.66409200544035718750E+14
       coef(  325) =    -0.29726903306849542969E+14
       coef(  326) =     0.67083742526889687500E+14
       coef(  327) =     0.47553089590213039062E+14
       coef(  328) =     0.12653196299474253906E+14
       coef(  329) =     0.38174517972699078125E+14
       coef(  330) =    -0.20698022662500011719E+14
       coef(  331) =    -0.12835083616085492188E+15
       coef(  332) =    -0.44189152534316710938E+14
       coef(  333) =    -0.21117902824644765625E+14
       coef(  334) =    -0.70002611461686984375E+14
       coef(  335) =    -0.12800317433034597656E+14
       coef(  336) =    -0.35764862648603164062E+14
       coef(  337) =    -0.28170096013419189453E+13
       coef(  338) =    -0.10133934430179307861E+13
       coef(  339) =     0.10921087492441395312E+15
       coef(  340) =    -0.45204367432022968750E+14
       coef(  341) =    -0.55014912144515343750E+14
       coef(  342) =    -0.18973123365654085938E+14
       coef(  343) =    -0.12673683197692648438E+14
       coef(  344) =    -0.21644725133636039062E+14
       coef(  345) =     0.34176605337142324219E+13
       coef(  346) =    -0.19291268952517883301E+13
       coef(  347) =     0.63707677884558027344E+13
       coef(  348) =     0.15373984783799277344E+14
       coef(  349) =    -0.30915184135177445312E+14
       coef(  350) =     0.18603701833469804688E+14
       coef(  351) =     0.15471743557478146484E+14
       coef(  352) =     0.11132885693556826172E+14
       coef(  353) =    -0.22991768521535730469E+14
       coef(  354) =     0.74309091770234648438E+13
       coef(  355) =     0.32107334375450462103E+07
       coef(  356) =    -0.34840498539721906185E+08
       coef(  357) =     0.60229485517400093079E+10
       coef(  358) =     0.89130703338951049805E+12
       coef(  359) =    -0.12446451686108927734E+14
       coef(  360) =     0.46808668949869406250E+14
       coef(  361) =    -0.46652060055737226562E+14
       coef(  362) =    -0.31796504761644557953E+11
       coef(  363) =     0.27168630169781893921E+12
       coef(  364) =     0.30760675871999023438E+13
       coef(  365) =     0.33627796280624156250E+14
       coef(  366) =    -0.73581235736154625000E+14
       coef(  367) =    -0.19639892998368587500E+15
       coef(  368) =     0.16002826902664350586E+13
       coef(  369) =    -0.35584793127444515625E+14
       coef(  370) =    -0.94206720516240203125E+14
       coef(  371) =    -0.10916158847962356250E+15
       coef(  372) =    -0.87554442406863593750E+14
       coef(  373) =     0.13656941657148978125E+15
       coef(  374) =     0.10927266949029234375E+15
       coef(  375) =     0.15423721628162324219E+14
       coef(  376) =     0.85558014739183671875E+14
       coef(  377) =    -0.30088725983605603027E+12
       coef(  378) =    -0.15535873382374472656E+14
       coef(  379) =     0.36239388603259007812E+14
       coef(  380) =     0.67908737830307804688E+14
       coef(  381) =     0.81548227864377406250E+14
       coef(  382) =     0.13055368195874746094E+14
       coef(  383) =    -0.19284724853868441406E+14
       coef(  384) =     0.44029296620210328125E+14
       coef(  385) =     0.40627328226952304688E+13
       coef(  386) =    -0.93203140574777578125E+13
       coef(  387) =     0.80000472186222843750E+14
       coef(  388) =     0.37109027805404257812E+14
       coef(  389) =     0.90581395767513375000E+14
       coef(  390) =    -0.19515989451962082031E+14
       coef(  391) =    -0.11927188321234660156E+14
       coef(  392) =    -0.12399162712372914062E+14
       coef(  393) =     0.13786315449879537109E+14
       coef(  394) =    -0.25882888093827531250E+14
       coef(  395) =    -0.11419277408445657812E+15
       coef(  396) =    -0.55229636283846179688E+14
       coef(  397) =    -0.17598615885589941406E+14
       coef(  398) =     0.15346372257639085938E+14
       coef(  399) =    -0.22100131818320148438E+14
       coef(  400) =    -0.69919903093149160156E+13
       coef(  401) =     0.34413719299361010742E+13
       coef(  402) =     0.43019296144481679688E+13
       coef(  403) =     0.17417136775479423828E+13
       coef(  404) =    -0.29485110578170083008E+13
       coef(  405) =    -0.10507614256479880214E+09
       coef(  406) =     0.85263170359583206177E+10
       coef(  407) =    -0.26079878576371908569E+12
       coef(  408) =     0.29920989356598059082E+12
       coef(  409) =     0.15532060711511527344E+14
       coef(  410) =    -0.70636350666644671875E+14
       coef(  411) =     0.79704128557722703125E+14
       coef(  412) =     0.86534498699532089233E+11
       coef(  413) =    -0.11971519025825866699E+12
       coef(  414) =    -0.36721819583310097656E+13
       coef(  415) =    -0.77826398308117703125E+14
       coef(  416) =     0.19592926683112543750E+15
       coef(  417) =     0.54184346137543828125E+14
       coef(  418) =    -0.28806148996690776367E+13
       coef(  419) =     0.11688913945432078125E+14
       coef(  420) =     0.10586968677053287500E+15
       coef(  421) =     0.61353268322202125000E+14
       coef(  422) =    -0.58093713180601898438E+14
       coef(  423) =     0.39854691349657164062E+14
       coef(  424) =     0.35088002050745825195E+12
       coef(  425) =     0.39599563419950000000E+14
       coef(  426) =    -0.94414805660104218750E+14
       coef(  427) =     0.34244069582140562500E+14
       coef(  428) =     0.94876221327834828125E+14
       coef(  429) =    -0.22685000137188914062E+14
       coef(  430) =     0.58202119571295351562E+13
       coef(  431) =     0.17721432580530078125E+14
       coef(  432) =     0.24387837548928539062E+14
       coef(  433) =     0.96111835189792953125E+14
       coef(  434) =     0.13893458231152826172E+14
       coef(  435) =    -0.31267699102682710938E+14
       coef(  436) =     0.78642346736087046875E+14
       coef(  437) =    -0.81784540372366140625E+14
       coef(  438) =    -0.22246828590004339844E+14
       coef(  439) =     0.11504498080423330078E+14
       coef(  440) =    -0.16097786237530406250E+14
       coef(  441) =     0.59483577218670214844E+13
       coef(  442) =     0.41350146161719119549E+09
       coef(  443) =    -0.37742278065852851868E+11
       coef(  444) =     0.88770108213490295410E+12
       coef(  445) =    -0.67685634900463183594E+13
       coef(  446) =     0.24246467301736691406E+14
       coef(  447) =    -0.58638183175340515625E+14
       coef(  448) =     0.78740981880496656250E+14
       coef(  449) =     0.29309394078786462402E+12
       coef(  450) =    -0.93013300711318984375E+13
       coef(  451) =     0.52577613083181937500E+14
       coef(  452) =    -0.11033635701742046875E+15
       coef(  453) =     0.13662510255202378125E+15
       coef(  454) =     0.29284076611331218750E+14
       coef(  455) =    -0.73520021887513406250E+14
       coef(  456) =     0.63536465080334937500E+14
       coef(  457) =    -0.14574727623922890625E+15
       coef(  458) =    -0.96802073806764526367E+12
       coef(  459) =    -0.12609386705853591797E+14
       coef(  460) =     0.29754833143500480469E+14
       coef(  461) =     0.34924433801611511719E+14
       coef(  462) =    -0.97537306075957531250E+14
       coef(  463) =    -0.16767410802277994141E+14
       coef(  464) =     0.16712267824176820312E+14
       coef(  465) =     0.24099216130121878906E+14
       coef(  466) =    -0.10390744981838066406E+14
       coef(  467) =    -0.30691401195740738281E+14
       coef(  468) =    -0.26326025762403886719E+14
       coef(  469) =    -0.25293692401780424461E+03
       coef(  470) =    -0.45442475999839487486E+05
       coef(  471) =     0.40540454755092780106E+07
       coef(  472) =    -0.43983519733571842313E+08
       coef(  473) =     0.69647839654000020027E+09
       coef(  474) =     0.84238578718974037170E+10
       coef(  475) =    -0.49991259451593544006E+11
       coef(  476) =     0.61252374147418448702E+07
       coef(  477) =    -0.93998442397919344902E+09
       coef(  478) =     0.66935958521154336929E+10
       coef(  479) =    -0.13183054955192230225E+12
       coef(  480) =     0.43904376098758471680E+12
       coef(  481) =     0.91172515416058868408E+11
       coef(  482) =     0.16270991972833166122E+11
       coef(  483) =    -0.32801571473340447998E+12
       coef(  484) =     0.19678910559461909180E+13
       coef(  485) =    -0.55318995031461005859E+13
       coef(  486) =     0.10257089116826298828E+14
       coef(  487) =     0.23980382106157622070E+13
       coef(  488) =    -0.22756918235280089844E+14
       coef(  489) =     0.31431831408255390625E+14
       coef(  490) =    -0.45498824318129562500E+14
       coef(  491) =     0.47321153381392335938E+14
       coef(  492) =    -0.75855642312329500000E+14
       coef(  493) =     0.32616293821937718750E+14
       coef(  494) =    -0.14151542651868216797E+14
       coef(  495) =     0.38906937865134515625E+14
       coef(  496) =     0.36047163832995128632E+09
       coef(  497) =     0.98802164709960174561E+10
       coef(  498) =     0.21096372100150009155E+12
       coef(  499) =     0.76886181902509985352E+12
       coef(  500) =    -0.87526833347444531250E+13
       coef(  501) =     0.99412936065591796875E+13
       coef(  502) =    -0.42088174796636657715E+12
       coef(  503) =    -0.24500039836734682617E+13
       coef(  504) =     0.75668506313939365234E+13
       coef(  505) =     0.47533444966612304688E+14
       coef(  506) =    -0.12217076345893743750E+15
       coef(  507) =     0.18510155671562648438E+14
       coef(  508) =    -0.87898350046472421875E+14
       coef(  509) =     0.13332870462842668750E+15
       coef(  510) =     0.97559838662910765625E+14
       coef(  511) =    -0.66139180724138179688E+14
       coef(  512) =     0.55713059583675601562E+14
       coef(  513) =     0.55213609298185695312E+14
       coef(  514) =     0.44705355303999578125E+14
       coef(  515) =     0.49221165113513066406E+13
       coef(  516) =    -0.43573639622402242188E+14
       coef(  517) =     0.22695460175292214844E+14
       coef(  518) =     0.50077878837871093750E+12
       coef(  519) =    -0.46161454057623429688E+14
       coef(  520) =     0.10970482986086942188E+15
       coef(  521) =    -0.72394160249761328125E+14
       coef(  522) =    -0.19735750308914570312E+14
       coef(  523) =    -0.45212299074125253906E+13
       coef(  524) =    -0.64868923514452398438E+14
       coef(  525) =    -0.12354997423049746094E+14
       coef(  526) =    -0.23240096657329523438E+14
       coef(  527) =    -0.70802732351303437500E+14
       coef(  528) =    -0.24717003058147953125E+14
       coef(  529) =    -0.37605697767387312500E+14
       coef(  530) =    -0.12248514823749221802E+11
       coef(  531) =     0.17989989087797671509E+12
       coef(  532) =    -0.24996312560742758789E+13
       coef(  533) =     0.35392136344505849609E+13
       coef(  534) =     0.11623712617125363281E+14
       coef(  535) =    -0.56632431759369433594E+13
       coef(  536) =     0.12835608952538146973E+13
       coef(  537) =     0.20312645740340968750E+14
       coef(  538) =    -0.53212812705904921875E+14
       coef(  539) =    -0.32849598002985488281E+13
       coef(  540) =     0.37479806363611679688E+13
       coef(  541) =    -0.23970608774791429688E+14
       coef(  542) =     0.47992629430505929688E+14
       coef(  543) =     0.82530107936379968750E+14
       coef(  544) =     0.48498072966627000000E+14
       coef(  545) =     0.50594300176932500000E+14
       coef(  546) =     0.41411746067816562500E+14
       coef(  547) =    -0.23691590685230460938E+14
       coef(  548) =     0.46073888108548625000E+14
       coef(  549) =     0.12063317049379964844E+14
       coef(  550) =     0.54017774857621904297E+13
       coef(  551) =     0.97209903860148437500E+11
       coef(  552) =    -0.56364825248912007812E+14
       coef(  553) =    -0.17199969968579628906E+14
       coef(  554) =     0.36363073939946562500E+13
       coef(  555) =    -0.34895681178530781250E+13
       coef(  556) =    -0.59488623744764968750E+14
       coef(  557) =    -0.19540698839911347656E+14
       coef(  558) =     0.50788837568195328125E+14
       coef(  559) =    -0.83877548237574281250E+14
       coef(  560) =    -0.25039174267052019531E+14
       coef(  561) =    -0.47550992354669042969E+13
       coef(  562) =    -0.70140836792123664062E+14
       coef(  563) =    -0.18465165654797972656E+14
       coef(  564) =    -0.32177536200238871094E+14
       coef(  565) =     0.47454996281615778571E+05
       coef(  566) =    -0.63522128025911543518E+07
       coef(  567) =     0.48775456502316278219E+09
       coef(  568) =    -0.22710441911161247253E+11
       coef(  569) =     0.30465557618001306152E+12
       coef(  570) =    -0.15233657634424062500E+13
       coef(  571) =     0.24976509828449511719E+13
       coef(  572) =    -0.43889841304299646616E+09
       coef(  573) =     0.56118772028504653931E+11
       coef(  574) =    -0.61221555106768054199E+12
       coef(  575) =     0.25763525049014709473E+12
       coef(  576) =     0.12587475763106097656E+14
       coef(  577) =    -0.33907266555638464844E+14
       coef(  578) =    -0.43457855537490295410E+12
       coef(  579) =     0.76961606763051220703E+13
       coef(  580) =     0.67906600496389970703E+13
       coef(  581) =    -0.72428356783121890625E+14
       coef(  582) =     0.11847471645006735938E+15
       coef(  583) =    -0.46188049004192960938E+14
       coef(  584) =     0.15347634723067965625E+15
       coef(  585) =    -0.95624741974856984375E+14
       coef(  586) =     0.18340934842460820312E+14
       coef(  587) =    -0.61173288060369179688E+14
       coef(  588) =    -0.59471118385667320312E+14
       coef(  589) =     0.16967777669796902344E+14
       coef(  590) =    -0.53789618279649707031E+13
       coef(  591) =    -0.29097272588852615356E+11
       coef(  592) =     0.74682600990187145996E+12
       coef(  593) =    -0.14034767803247488403E+12
       coef(  594) =    -0.30456886864822691406E+14
       coef(  595) =     0.75117949846614109375E+14
       coef(  596) =     0.24979911693387363281E+14
       coef(  597) =    -0.70781277048952958984E+13
       coef(  598) =     0.55059715221816101562E+14
       coef(  599) =    -0.82094069566138687500E+14
       coef(  600) =    -0.61408284176633726562E+14
       coef(  601) =    -0.19310972873468789062E+14
       coef(  602) =    -0.25181141077522855469E+14
       coef(  603) =     0.26916392991199359375E+14
       coef(  604) =     0.37370861658094343750E+14
       coef(  605) =     0.28752291805131785156E+14
       coef(  606) =     0.30992425333497371094E+14
       coef(  607) =     0.26634981844925734375E+14
       coef(  608) =     0.87510385194352304688E+13
       coef(  609) =    -0.50076260011132351562E+14
       coef(  610) =     0.38010844391378171875E+14
       coef(  611) =     0.18011973693837531250E+14
       coef(  612) =     0.40356281384237773438E+13
       coef(  613) =    -0.66318736429715054688E+14
       coef(  614) =    -0.52023062493396132812E+13
       coef(  615) =     0.69121935569566210938E+13
       coef(  616) =     0.33063441162678027344E+12
       coef(  617) =    -0.55305896857649734375E+14
       coef(  618) =    -0.15434844463826898438E+14
       coef(  619) =    -0.21268261418659286499E+11
       coef(  620) =    -0.45378917910545336914E+12
       coef(  621) =     0.74654425805981572266E+13
       coef(  622) =     0.14027735617720210938E+14
       coef(  623) =    -0.90698354117092609375E+14
       coef(  624) =    -0.71844079661043656250E+14
       coef(  625) =    -0.53991725833209414062E+13
       coef(  626) =     0.19862452070581003906E+14
       coef(  627) =     0.51827228708954062500E+14
       coef(  628) =     0.10331519683969382812E+14
       coef(  629) =    -0.13243765045513012695E+13
       coef(  630) =    -0.11734650388938398438E+13
       coef(  631) =     0.29719229166607378906E+14
       coef(  632) =     0.20324457303546550781E+14
       coef(  633) =     0.19128178558970312500E+14
       coef(  634) =     0.20030129127458867188E+13
       coef(  635) =     0.82677470280295078125E+14
       coef(  636) =     0.49185301056400695312E+14
       coef(  637) =     0.15586491730546388672E+14
       coef(  638) =     0.23237728339664746094E+14
       coef(  639) =     0.11134183686021744141E+14
       coef(  640) =     0.11424355135870820312E+13
       coef(  641) =    -0.10029106697358082812E+15
       coef(  642) =    -0.40510103009266958008E+13
       coef(  643) =     0.82338257277192216797E+13
       coef(  644) =    -0.29476909286107207031E+13
       coef(  645) =     0.14052968499580078060E+06
       coef(  646) =     0.66805197010946579278E+08
       coef(  647) =    -0.48398688480602874756E+10
       coef(  648) =     0.14260898189309567261E+12
       coef(  649) =    -0.12766283165346022949E+13
       coef(  650) =     0.47459047510648955078E+13
       coef(  651) =    -0.56250613711980507812E+13
       coef(  652) =     0.36833255701896400452E+10
       coef(  653) =    -0.29922302610319860840E+12
       coef(  654) =     0.20639779188337602539E+13
       coef(  655) =    -0.41936692961402441406E+13
       coef(  656) =    -0.76315464204759443359E+13
       coef(  657) =     0.14754249527845121094E+14
       coef(  658) =     0.12859815639692460938E+13
       coef(  659) =    -0.15371452881871402344E+14
       coef(  660) =    -0.90158194996639140625E+13
       coef(  661) =     0.12740156789739816406E+14
       coef(  662) =     0.51666031620757210938E+14
       coef(  663) =     0.67516958906981437500E+14
       coef(  664) =    -0.60919941849658890625E+14
       coef(  665) =    -0.40026692160814976562E+14
       coef(  666) =     0.10287962595200044922E+14
       coef(  667) =    -0.67765755275636531250E+14
       coef(  668) =    -0.19050336929366578125E+14
       coef(  669) =     0.90840285410390991211E+11
       coef(  670) =     0.30667844577349822998E+12
       coef(  671) =     0.56921211892662890625E+12
       coef(  672) =     0.31062194355505156250E+14
       coef(  673) =    -0.22877042144120679688E+14
       coef(  674) =    -0.58649421137732968750E+13
       coef(  675) =    -0.78017043324959658203E+13
       coef(  676) =    -0.15625991056044107422E+14
       coef(  677) =     0.34351968114192179688E+14
       coef(  678) =     0.78340077210627949219E+13
       coef(  679) =     0.50105139999874931641E+13
       coef(  680) =    -0.92141117476475039062E+13
       coef(  681) =     0.30703091415568132812E+14
       coef(  682) =     0.19064641216734968750E+14
       coef(  683) =     0.21478883368602503906E+14
       coef(  684) =     0.13846574222801996094E+14
       coef(  685) =     0.23002372929274394531E+14
       coef(  686) =     0.48643032523256039062E+14
       coef(  687) =     0.20255680562424261719E+14
       coef(  688) =     0.19012869784822617188E+13
       coef(  689) =     0.12458768364001910156E+14
       coef(  690) =    -0.57512937333863535156E+13
       coef(  691) =    -0.86532759852196179199E+12
       coef(  692) =     0.14371324327568730469E+14
       coef(  693) =    -0.67386594324895828125E+14
       coef(  694) =     0.74352895898350296875E+14
       coef(  695) =     0.49213272127243710938E+13
       coef(  696) =    -0.71362577542233105469E+13
       coef(  697) =    -0.20032677934696445312E+14
       coef(  698) =    -0.73623386351088427734E+13
       coef(  699) =     0.30983475162681031250E+14
       coef(  700) =     0.11186472054099242188E+14
       coef(  701) =    -0.44638648214255458984E+13
       coef(  702) =     0.10625993168809757812E+14
       coef(  703) =     0.14076048739538125000E+14
       coef(  704) =     0.26745134691042847656E+14
       coef(  705) =     0.18143482905697078125E+14
       coef(  706) =     0.94836590564172890625E+13
       coef(  707) =    -0.16923399325387269531E+14
       coef(  708) =     0.37037473134873857422E+13
       coef(  709) =    -0.32117501775620233268E+08
       coef(  710) =    -0.10941015065252690315E+10
       coef(  711) =     0.10518161274848971558E+12
       coef(  712) =    -0.21259042668483857422E+13
       coef(  713) =     0.13017878584257628906E+14
       coef(  714) =    -0.34883760416380382812E+14
       coef(  715) =     0.34919931653789402344E+14
       coef(  716) =    -0.11109598441299419403E+11
       coef(  717) =     0.24389906812812827148E+13
       coef(  718) =    -0.11558433504022994141E+14
       coef(  719) =     0.36242617169803976562E+14
       coef(  720) =     0.14461288531457685547E+13
       coef(  721) =    -0.75449663232077468750E+14
       coef(  722) =     0.15509830240534912109E+12
       coef(  723) =    -0.34009113623723750000E+13
       coef(  724) =     0.41391162398946656250E+14
       coef(  725) =     0.27367432888591152344E+14
       coef(  726) =     0.54794701584819921875E+12
       coef(  727) =    -0.11271088965760021875E+15
       coef(  728) =    -0.45956771322631601562E+14
       coef(  729) =    -0.95101972195507832031E+13
       coef(  730) =    -0.20569881144124382812E+14
       coef(  731) =    -0.14277883378443803711E+13
       coef(  732) =    -0.13649146278102875000E+14
       coef(  733) =     0.41845323395963843750E+14
       coef(  734) =    -0.12932459259959501562E+15
       coef(  735) =    -0.75403218415601437500E+14
       coef(  736) =    -0.37131169886796468750E+14
       coef(  737) =     0.86680466170484718750E+14
       coef(  738) =    -0.83107696341023312500E+14
       coef(  739) =    -0.36796499173855828125E+14
       coef(  740) =    -0.13451296043376123047E+14
       coef(  741) =    -0.46128786745500343750E+14
       coef(  742) =    -0.36164304762577465820E+13
       coef(  743) =    -0.68731046908904101562E+14
       coef(  744) =    -0.48597104858888929688E+14
       coef(  745) =    -0.59660035308337832031E+13
       coef(  746) =    -0.17591446373196867188E+14
       coef(  747) =     0.97363612900527597656E+13
       coef(  748) =    -0.59551481173180226562E+14
       coef(  749) =     0.10164649337706401562E+15
       coef(  750) =     0.33876508247845765625E+14
       coef(  751) =    -0.36403223031340415039E+13
       coef(  752) =     0.12974711816471893750E+15
       coef(  753) =     0.32768690680214117188E+14
       coef(  754) =     0.10411275110854005859E+14
       coef(  755) =     0.36277043715769492188E+13
       coef(  756) =     0.62970707827033093750E+14
       coef(  757) =     0.16346688249383990234E+14
       coef(  758) =     0.15899502209416042969E+14
       coef(  759) =     0.59121551063910531998E+09
       coef(  760) =    -0.19363509101604141235E+11
       coef(  761) =    -0.38170522516834826660E+12
       coef(  762) =     0.83010550319545898438E+13
       coef(  763) =    -0.34669406428452468750E+14
       coef(  764) =     0.40611372800213359375E+14
       coef(  765) =    -0.19636845062510289062E+14
       coef(  766) =     0.33285532797548773193E+12
       coef(  767) =    -0.11360622234340951172E+14
       coef(  768) =     0.15222738863856302734E+14
       coef(  769) =    -0.49769194097584648438E+14
       coef(  770) =     0.95829088004929515625E+14
       coef(  771) =     0.17142826633301218750E+14
       coef(  772) =     0.16770066037159902344E+13
       coef(  773) =     0.87136456786414656250E+14
       coef(  774) =     0.11271648131665773438E+15
       coef(  775) =     0.66832945836824429688E+14
       coef(  776) =    -0.72941700642709000000E+14
       coef(  777) =    -0.86299112101869179688E+13
       coef(  778) =     0.82619034227631455078E+13
       coef(  779) =     0.17351145938554445312E+14
       coef(  780) =    -0.98442909633922265625E+12
       coef(  781) =    -0.50408023063640960938E+14
       coef(  782) =    -0.53745955363496611328E+13
       coef(  783) =    -0.56180960582016179688E+14
       coef(  784) =    -0.53117309914628960938E+14
       coef(  785) =    -0.14065883509842564453E+14
       coef(  786) =    -0.25938753002023523438E+14
       coef(  787) =    -0.40493069397852382812E+14
       coef(  788) =    -0.22545647459297285156E+14
       coef(  789) =    -0.35251748743414281250E+14
       coef(  790) =    -0.24477285067416320312E+14
       coef(  791) =     0.69074132249913109375E+14
       coef(  792) =     0.22151961988264406250E+14
       coef(  793) =     0.42680155735393226562E+14
       coef(  794) =     0.13370994957974347656E+14
       coef(  795) =     0.23796295930892699219E+14
       coef(  796) =    -0.30675010961178936958E+10
       coef(  797) =     0.15304880115427731323E+12
       coef(  798) =     0.31471746257849578857E+12
       coef(  799) =    -0.16237144807002169922E+14
       coef(  800) =     0.52283692667516062500E+14
       coef(  801) =     0.11871980220930687500E+14
       coef(  802) =    -0.25338581432786527344E+14
       coef(  803) =    -0.26151929807547412109E+13
       coef(  804) =     0.44240825968554015625E+14
       coef(  805) =    -0.43260556616216921875E+14
       coef(  806) =    -0.68879355592759984375E+14
       coef(  807) =     0.64631512668782945312E+14
       coef(  808) =    -0.55001442945007445312E+14
       coef(  809) =     0.21492543918687242188E+14
       coef(  810) =     0.65570030566253742188E+14
       coef(  811) =    -0.61415874660886109375E+14
       coef(  812) =    -0.17785066354759585938E+14
       coef(  813) =    -0.38369716760240921875E+14
       coef(  814) =     0.29473460962050500000E+14
       coef(  815) =    -0.14073961752861269531E+13
       coef(  816) =    -0.64807379583635578125E+14
       coef(  817) =    -0.21482340973003605469E+14
       coef(  818) =    -0.17308802218390990234E+14
       coef(  819) =     0.10030044882897934375E+15
       coef(  820) =     0.40496402162661376953E+13
       coef(  821) =     0.40323985467292125000E+14
       coef(  822) =     0.46602278845609052734E+13
       coef(  823) =     0.47686315141695146561E+10
       coef(  824) =    -0.26584054904835202026E+12
       coef(  825) =     0.29351797410640637207E+12
       coef(  826) =     0.11558404167029542969E+14
       coef(  827) =    -0.22341125378882054688E+14
       coef(  828) =    -0.60428713705562281250E+14
       coef(  829) =     0.46566493729573525391E+13
       coef(  830) =    -0.55648410341690507812E+14
       coef(  831) =     0.53747032175657164062E+14
       coef(  832) =    -0.36147656377849437500E+14
       coef(  833) =     0.53617388665484906250E+14
       coef(  834) =     0.17120930046407929688E+14
       coef(  835) =     0.72238764041575488281E+13
       coef(  836) =     0.10862005508259725000E+15
       coef(  837) =     0.80379234969656703125E+14
       coef(  838) =    -0.27870248951971826172E+13
       coef(  839) =    -0.10527317054685178125E+15
       coef(  840) =    -0.10485616198960808594E+14
       coef(  841) =    -0.20876785791607209831E+05
       coef(  842) =     0.12907917591761196963E+07
       coef(  843) =    -0.15438806482243725657E+09
       coef(  844) =     0.88841116833034248352E+10
       coef(  845) =    -0.12516408000628300476E+12
       coef(  846) =     0.43879058226468847656E+12
       coef(  847) =    -0.16980038582159735107E+12
       coef(  848) =     0.13075431370516654849E+09
       coef(  849) =    -0.14104216077549655914E+11
       coef(  850) =     0.24747140147191549683E+12
       coef(  851) =    -0.23487599584456783295E+11
       coef(  852) =    -0.90727787987550061035E+12
       coef(  853) =    -0.68618940769886728516E+13
       coef(  854) =     0.13897575343705322266E+12
       coef(  855) =    -0.19435168823405593262E+13
       coef(  856) =    -0.90401917212569316406E+13
       coef(  857) =     0.30096794605987632812E+14
       coef(  858) =    -0.21861748448541457031E+14
       coef(  859) =     0.72154620360022597656E+13
       coef(  860) =     0.22321529009735445312E+14
       coef(  861) =    -0.47844400589564882812E+14
       coef(  862) =     0.51455840479563007812E+14
       coef(  863) =    -0.53472302035774132812E+14
       coef(  864) =    -0.30542479932793911133E+13
       coef(  865) =     0.53978556993894859375E+14
       coef(  866) =     0.30908444256613343750E+14
       coef(  867) =     0.22195890852487788200E+10
       coef(  868) =    -0.38662014728151953125E+12
       coef(  869) =    -0.28474826341943103027E+12
       coef(  870) =     0.10629070091830757812E+14
       coef(  871) =    -0.25622181278656101562E+14
       coef(  872) =     0.61683101852285570312E+14
       coef(  873) =     0.23294498219330561523E+13
       coef(  874) =    -0.25076806965805449219E+13
       coef(  875) =     0.11856403646170269531E+14
       coef(  876) =    -0.10604245124947593750E+15
       coef(  877) =    -0.51922970024363304688E+14
       coef(  878) =    -0.41614763756074703125E+14
       coef(  879) =     0.69391255676215687500E+14
       coef(  880) =     0.52432664013058609375E+14
       coef(  881) =     0.35725685210602671875E+14
       coef(  882) =     0.70999771423366406250E+14
       coef(  883) =     0.45434002939477429688E+14
       coef(  884) =    -0.88104851002326953125E+13
       coef(  885) =     0.49738541235082726562E+14
       coef(  886) =     0.75713374584961671875E+14
       coef(  887) =     0.58902722717688281250E+13
       coef(  888) =    -0.44654513322901386719E+13
       coef(  889) =    -0.35089522357939875000E+14
       coef(  890) =     0.46492706928977275391E+13
       coef(  891) =     0.77760192661384033203E+13
       coef(  892) =     0.57863802025050537109E+13
       coef(  893) =    -0.54159026559949062500E+14
       coef(  894) =    -0.16137121664940841797E+14
       coef(  895) =     0.26794246101781536865E+12
       coef(  896) =    -0.28915532459537656250E+13
       coef(  897) =     0.61050385459496738281E+13
       coef(  898) =    -0.41628661730083115234E+13
       coef(  899) =    -0.35877423855803218750E+14
       coef(  900) =     0.18482581213898562500E+14
       coef(  901) =     0.15425324973385087891E+14
       coef(  902) =    -0.98786601765239718750E+14
       coef(  903) =     0.67219425710088625000E+14
       coef(  904) =     0.38161769667060390625E+14
       coef(  905) =     0.20692987820586953125E+14
       coef(  906) =    -0.40684118926455718750E+14
       coef(  907) =     0.41914538798453796875E+14
       coef(  908) =     0.31765661130968984375E+14
       coef(  909) =     0.29556058659714179688E+14
       coef(  910) =     0.85784709399597558594E+13
       coef(  911) =     0.10605969679483239062E+15
       coef(  912) =     0.76329243775200562500E+14
       coef(  913) =     0.29585305734626039062E+14
       coef(  914) =     0.39294551013180593750E+14
       coef(  915) =     0.19940680372486304688E+14
       coef(  916) =     0.72611070598080546875E+13
       coef(  917) =    -0.46060088160801921875E+14
       coef(  918) =     0.23024982665636593750E+14
       coef(  919) =     0.20504065718048628906E+14
       coef(  920) =     0.61175509562475449219E+13
       coef(  921) =     0.14052968499580078060E+06
       coef(  922) =     0.66805197010946579278E+08
       coef(  923) =    -0.48398688480602874756E+10
       coef(  924) =     0.14260898189309567261E+12
       coef(  925) =    -0.12766283165346022949E+13
       coef(  926) =     0.47459047510648955078E+13
       coef(  927) =    -0.56250613711980507812E+13
       coef(  928) =     0.36833255701896400452E+10
       coef(  929) =    -0.29922302610319860840E+12
       coef(  930) =     0.20639779188337602539E+13
       coef(  931) =    -0.41936692961402441406E+13
       coef(  932) =    -0.76315464204759443359E+13
       coef(  933) =     0.14754249527845121094E+14
       coef(  934) =     0.12859815639692460938E+13
       coef(  935) =    -0.15371452881871402344E+14
       coef(  936) =    -0.90158194996639140625E+13
       coef(  937) =     0.12740156789739816406E+14
       coef(  938) =     0.51666031620757210938E+14
       coef(  939) =     0.67516958906981437500E+14
       coef(  940) =    -0.60919941849658890625E+14
       coef(  941) =    -0.40026692160814976562E+14
       coef(  942) =     0.10287962595200044922E+14
       coef(  943) =    -0.67765755275636531250E+14
       coef(  944) =    -0.19050336929366578125E+14
       coef(  945) =     0.90840285410390991211E+11
       coef(  946) =     0.30667844577349822998E+12
       coef(  947) =     0.56921211892662890625E+12
       coef(  948) =     0.31062194355505156250E+14
       coef(  949) =    -0.22877042144120679688E+14
       coef(  950) =    -0.58649421137732968750E+13
       coef(  951) =    -0.78017043324959658203E+13
       coef(  952) =    -0.15625991056044107422E+14
       coef(  953) =     0.34351968114192179688E+14
       coef(  954) =     0.78340077210627949219E+13
       coef(  955) =     0.50105139999874931641E+13
       coef(  956) =    -0.92141117476475039062E+13
       coef(  957) =     0.30703091415568132812E+14
       coef(  958) =     0.19064641216734968750E+14
       coef(  959) =     0.21478883368602503906E+14
       coef(  960) =     0.13846574222801996094E+14
       coef(  961) =     0.23002372929274394531E+14
       coef(  962) =     0.48643032523256039062E+14
       coef(  963) =     0.20255680562424261719E+14
       coef(  964) =     0.19012869784822617188E+13
       coef(  965) =     0.12458768364001910156E+14
       coef(  966) =    -0.57512937333863535156E+13
       coef(  967) =    -0.86532759852196179199E+12
       coef(  968) =     0.14371324327568730469E+14
       coef(  969) =    -0.67386594324895828125E+14
       coef(  970) =     0.74352895898350296875E+14
       coef(  971) =     0.49213272127243710938E+13
       coef(  972) =    -0.71362577542233105469E+13
       coef(  973) =    -0.20032677934696445312E+14
       coef(  974) =    -0.73623386351088427734E+13
       coef(  975) =     0.30983475162681031250E+14
       coef(  976) =     0.11186472054099242188E+14
       coef(  977) =    -0.44638648214255458984E+13
       coef(  978) =     0.10625993168809757812E+14
       coef(  979) =     0.14076048739538125000E+14
       coef(  980) =     0.26745134691042847656E+14
       coef(  981) =     0.18143482905697078125E+14
       coef(  982) =     0.94836590564172890625E+13

      r1=r(1)*autoang
      r2=r(2)*autoang
      r3=r(3)*autoang
      r4=r(4)*autoang
      r5=r(5)*autoang
      r6=r(6)*autoang

      b1=1.d0

      e1=dexp(-r1/b1)
      e2=dexp(-r2/b1)
      e3=dexp(-r3/b1)
      e4=dexp(-r4/b1)
      e5=dexp(-r5/b1)
      e6=dexp(-r6/b1)

      norder=6
      norder2=6
      maxorder=12

      ii=0
      e=0.d0
      do i1=0,norder
      do i2=0,norder
      x1=e1**i1*e2**i2+e1**i2*e2**i1
      do i3=0,norder2
      do i4=i3,norder2
      do i5=i4,norder2
      do i6=i5,norder2
      if ((i1+i2+i3+i4+i5+i6).le.maxorder) then
      x2=x1*(
     & e3**i3*e4**i4*e5**i5*e6**i6 +
     & e3**i3*e4**i4*e6**i5*e5**i6 +
     & e3**i3*e5**i4*e6**i5*e4**i6 +
     & e3**i3*e5**i4*e4**i5*e6**i6 +
     & e3**i3*e6**i4*e5**i5*e4**i6 +
     & e3**i3*e6**i4*e4**i5*e5**i6 +
     & e4**i3*e5**i4*e6**i5*e3**i6 +
     & e4**i3*e5**i4*e3**i5*e6**i6 +
     & e4**i3*e6**i4*e3**i5*e5**i6 +
     & e4**i3*e6**i4*e5**i5*e3**i6 +
     & e4**i3*e3**i4*e6**i5*e5**i6 +
     & e4**i3*e3**i4*e5**i5*e6**i6 +
     & e5**i3*e6**i4*e3**i5*e4**i6 +
     & e5**i3*e6**i4*e4**i5*e3**i6 +
     & e5**i3*e3**i4*e4**i5*e6**i6 +
     & e5**i3*e3**i4*e6**i5*e4**i6 +
     & e5**i3*e4**i4*e3**i5*e6**i6 +
     & e5**i3*e4**i4*e6**i5*e3**i6 +
     & e6**i3*e3**i4*e4**i5*e5**i6 +
     & e6**i3*e3**i4*e5**i5*e4**i6 +
     & e6**i3*e4**i4*e5**i5*e3**i6 +
     & e6**i3*e4**i4*e3**i5*e5**i6 +
     & e6**i3*e5**i4*e4**i5*e3**i6 +
     & e6**i3*e5**i4*e3**i5*e4**i6 )
      ii=ii+1
      basis(ii)=1.d0
      if (i1+i2+i3+i4+i5+i6.ne.0) basis(ii)=x2
      endif
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
c        print *,rr,v
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

