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

       ncoef=1535
       coef(    1) =    -0.82753270477642049130E+00
       coef(    2) =    -0.19188880193000096597E+03
       coef(    3) =    -0.39976008634628538857E+05
       coef(    4) =     0.21997549039769386873E+07
       coef(    5) =    -0.59869767625177867711E+08
       coef(    6) =    -0.26380396564509165287E+09
       coef(    7) =     0.54924124664944410324E+10
       coef(    8) =    -0.32755037416067775726E+11
       coef(    9) =     0.60957393750258087158E+11
       coef(   10) =     0.10151417156504071318E+06
       coef(   11) =    -0.54686019272987935692E+07
       coef(   12) =     0.30947495755970585346E+09
       coef(   13) =     0.10634256256786628962E+10
       coef(   14) =     0.11797735692908117294E+11
       coef(   15) =    -0.37621422542478637695E+11
       coef(   16) =     0.27110755729008856201E+12
       coef(   17) =    -0.83700180054133239746E+12
       coef(   18) =    -0.17722122648091754317E+09
       coef(   19) =    -0.91022845621515512466E+09
       coef(   20) =     0.20351409416812862396E+11
       coef(   21) =    -0.11663114414444318848E+13
       coef(   22) =     0.58739533652519150391E+13
       coef(   23) =    -0.78975770766694345703E+13
       coef(   24) =    -0.21559840687816220703E+13
       coef(   25) =     0.18821675940491455078E+11
       coef(   26) =    -0.14444255598203735352E+11
       coef(   27) =     0.10682068534215443359E+14
       coef(   28) =    -0.67171338716366281250E+14
       coef(   29) =     0.11766003528863089062E+15
       coef(   30) =    -0.20904274779757011719E+14
       coef(   31) =    -0.49738193963141445312E+13
       coef(   32) =     0.19626210100534164062E+14
       coef(   33) =     0.51466662072815843750E+14
       coef(   34) =    -0.55110834085907171875E+14
       coef(   35) =    -0.14140547309548678125E+15
       coef(   36) =    -0.45314296728528179688E+14
       coef(   37) =     0.19062781340436992188E+14
       coef(   38) =    -0.11404856577642771484E+14
       coef(   39) =    -0.41422962384827218750E+14
       coef(   40) =     0.43859119265053320312E+14
       coef(   41) =     0.24543325643024761719E+14
       coef(   42) =    -0.75680064302416536957E+07
       coef(   43) =     0.61168630989513814449E+09
       coef(   44) =    -0.21149529057388183594E+11
       coef(   45) =     0.96579666692120513916E+11
       coef(   46) =    -0.88374291944176379395E+12
       coef(   47) =     0.21680408274659306641E+13
       coef(   48) =    -0.62478397474002656250E+13
       coef(   49) =     0.16460561661717650391E+14
       coef(   50) =     0.64776942969758868217E+09
       coef(   51) =     0.12739225862300213623E+12
       coef(   52) =     0.48321195903701542969E+13
       coef(   53) =    -0.25075765590305750000E+14
       coef(   54) =     0.92996575559203406250E+14
       coef(   55) =    -0.18275523048517746875E+15
       coef(   56) =     0.10701044957145329688E+15
       coef(   57) =    -0.37612443400816313477E+13
       coef(   58) =     0.11709304864710320312E+14
       coef(   59) =     0.46971323346890742188E+13
       coef(   60) =    -0.27720985774001769531E+14
       coef(   61) =    -0.57488308294830914062E+14
       coef(   62) =     0.20672567451995128906E+14
       coef(   63) =    -0.24175843730208570312E+14
       coef(   64) =     0.44607944538587796875E+14
       coef(   65) =     0.40361951519378531250E+14
       coef(   66) =     0.10213677448819093750E+14
       coef(   67) =     0.58509952456306904297E+13
       coef(   68) =    -0.11220715397619941406E+14
       coef(   69) =    -0.27620631536434853516E+13
       coef(   70) =     0.99637908643927221680E+12
       coef(   71) =    -0.29955592206179433594E+13
       coef(   72) =     0.19200952342320129395E+11
       coef(   73) =    -0.14884241759966020508E+13
       coef(   74) =    -0.17440887431535156250E+14
       coef(   75) =     0.44068824228117500000E+14
       coef(   76) =    -0.46762754511037484375E+14
       coef(   77) =    -0.81662522825377500000E+13
       coef(   78) =     0.11622712804181071875E+15
       coef(   79) =     0.26874208072091636719E+14
       coef(   80) =    -0.38216813794353413086E+13
       coef(   81) =    -0.42623150597477742188E+14
       coef(   82) =     0.15831613198028566406E+14
       coef(   83) =     0.48482975798567750000E+14
       coef(   84) =     0.56372467882725796875E+14
       coef(   85) =    -0.17155779870489351562E+14
       coef(   86) =    -0.25108840904720742188E+14
       coef(   87) =     0.91389776574814023438E+13
       coef(   88) =     0.19234300698894457031E+14
       coef(   89) =    -0.24784464942000070312E+14
       coef(   90) =    -0.57571984969355605469E+13
       coef(   91) =    -0.58804449048206195312E+14
       coef(   92) =     0.78214470053029125000E+14
       coef(   93) =     0.18261788415272640625E+14
       coef(   94) =     0.24226638514868562500E+14
       coef(   95) =     0.27193518215087046875E+14
       coef(   96) =     0.62338687043995062500E+14
       coef(   97) =     0.15793869339002433594E+14
       coef(   98) =     0.93870779480112734375E+13
       coef(   99) =     0.21230100061881616211E+13
       coef(  100) =     0.37203321663810953125E+14
       coef(  101) =     0.10568110943081953125E+14
       coef(  102) =    -0.30224638648040622473E+08
       coef(  103) =     0.71457340103734636307E+09
       coef(  104) =     0.12781380448441915894E+12
       coef(  105) =    -0.31211327106174086914E+13
       coef(  106) =     0.32320471909483820312E+14
       coef(  107) =    -0.13117287850953640625E+15
       coef(  108) =     0.25537707233435581250E+15
       coef(  109) =    -0.23524283167393709375E+15
       coef(  110) =    -0.51141240295067947388E+11
       coef(  111) =     0.34344285161275517578E+13
       coef(  112) =    -0.40539045024003585938E+14
       coef(  113) =     0.76461921590510750000E+14
       coef(  114) =    -0.26592356544857359375E+14
       coef(  115) =     0.87808955573323687500E+14
       coef(  116) =     0.71387368156343437500E+14
       coef(  117) =     0.11177540075649292969E+14
       coef(  118) =    -0.10409250670330531250E+14
       coef(  119) =     0.16202678459770509766E+14
       coef(  120) =    -0.74725925257850976562E+12
       coef(  121) =     0.24947583649755179688E+14
       coef(  122) =     0.27435323904665703125E+14
       coef(  123) =    -0.29165094410723210938E+14
       coef(  124) =    -0.10629951003111664062E+14
       coef(  125) =     0.24809743261805400391E+13
       coef(  126) =     0.93017604004907011719E+13
       coef(  127) =    -0.15687924998888408203E+14
       coef(  128) =    -0.50075781299693818359E+13
       coef(  129) =    -0.23761847199998037109E+13
       coef(  130) =     0.16928006276485976562E+14
       coef(  131) =     0.73590114702393921875E+14
       coef(  132) =    -0.17859562794009181250E+15
       coef(  133) =    -0.13735624342474195312E+15
       coef(  134) =    -0.17571874168239281250E+14
       coef(  135) =     0.21617363509571832031E+14
       coef(  136) =    -0.10331679233537912500E+15
       coef(  137) =     0.56714609819829335938E+14
       coef(  138) =    -0.12472048999379429688E+14
       coef(  139) =    -0.13999100187050658203E+14
       coef(  140) =     0.40173689075016308594E+13
       coef(  141) =     0.55127614171793507812E+14
       coef(  142) =     0.14945488359374078125E+14
       coef(  143) =     0.47944634975779609375E+13
       coef(  144) =     0.45745522916883251953E+13
       coef(  145) =     0.88947745497092984375E+14
       coef(  146) =     0.85498952278025359375E+14
       coef(  147) =     0.23029749624207906250E+14
       coef(  148) =     0.58345822419438388672E+13
       coef(  149) =     0.44036294899448234375E+14
       coef(  150) =     0.13043148684870037109E+14
       coef(  151) =     0.17707865617804296875E+14
       coef(  152) =     0.50889136257497539062E+13
       coef(  153) =    -0.44543034366975515625E+14
       coef(  154) =     0.15666432335200900000E+15
       coef(  155) =    -0.20615807998396718750E+14
       coef(  156) =    -0.36827136936695828125E+14
       coef(  157) =    -0.69973120874931328125E+13
       coef(  158) =    -0.28040659004748203125E+14
       coef(  159) =     0.48834517894033328125E+14
       coef(  160) =     0.74216843295468681641E+13
       coef(  161) =    -0.19118813028501459961E+13
       coef(  162) =     0.27853379651719070312E+14
       coef(  163) =     0.79136717765676074219E+13
       coef(  164) =     0.89126918850871640625E+13
       coef(  165) =     0.20329526006989660156E+14
       coef(  166) =     0.58188800711375156250E+13
       coef(  167) =     0.95135063306172148438E+13
       coef(  168) =     0.52728724647659218750E+13
       coef(  169) =     0.54509900044279736328E+13
       coef(  170) =     0.34606103252326056463E+03
       coef(  171) =    -0.50369777421015540313E+05
       coef(  172) =     0.63861037185736894608E+07
       coef(  173) =    -0.33366911325972288847E+09
       coef(  174) =     0.86039316634857063293E+10
       coef(  175) =    -0.68928483926863540649E+11
       coef(  176) =     0.29875991426586621094E+12
       coef(  177) =    -0.66638544127537524414E+12
       coef(  178) =     0.53835454345100909424E+12
       coef(  179) =     0.24036390570724927820E+07
       coef(  180) =    -0.30937514307694387436E+09
       coef(  181) =     0.28345153121747970581E+07
       coef(  182) =    -0.26699843570031042480E+12
       coef(  183) =     0.24237188578794760742E+13
       coef(  184) =    -0.10696524458528146484E+14
       coef(  185) =     0.24133197489907347656E+14
       coef(  186) =    -0.20637526752884128906E+14
       coef(  187) =     0.95837382813488941193E+10
       coef(  188) =    -0.57635030538178596497E+11
       coef(  189) =     0.19982679331585126953E+13
       coef(  190) =    -0.15928053045859501953E+14
       coef(  191) =     0.66342786189053593750E+14
       coef(  192) =    -0.13841152251246976562E+15
       coef(  193) =     0.11855155979043557812E+15
       coef(  194) =    -0.52191891682675549316E+12
       coef(  195) =     0.14338017596272907715E+13
       coef(  196) =     0.49474475472084365234E+13
       coef(  197) =    -0.18734719837799726562E+14
       coef(  198) =    -0.49705046201453601562E+14
       coef(  199) =     0.22549415978763898438E+14
       coef(  200) =    -0.26849118067314824219E+13
       coef(  201) =     0.12631857335724642578E+14
       coef(  202) =     0.62234984490664203125E+14
       coef(  203) =     0.42790939394521351562E+14
       coef(  204) =     0.29084133341962843750E+14
       coef(  205) =    -0.54928295787007343750E+14
       coef(  206) =    -0.76502376876825341797E+13
       coef(  207) =     0.94334640411798925781E+13
       coef(  208) =    -0.50601941567893564453E+13
       coef(  209) =     0.17604968712697350979E+09
       coef(  210) =     0.50228090027932131290E+09
       coef(  211) =     0.49783105448380932617E+12
       coef(  212) =    -0.14909220370426335449E+13
       coef(  213) =    -0.12786658931384785156E+13
       coef(  214) =     0.16440796776263859375E+14
       coef(  215) =    -0.34091630417418210938E+14
       coef(  216) =     0.13020957952335367188E+14
       coef(  217) =    -0.47239180982933099365E+12
       coef(  218) =    -0.19588860872790595703E+13
       coef(  219) =     0.16976744472333652344E+14
       coef(  220) =    -0.55655311060630953125E+14
       coef(  221) =     0.62385535064531148438E+14
       coef(  222) =     0.14769419482407316406E+14
       coef(  223) =     0.22321328497896957031E+14
       coef(  224) =     0.10020182604654496094E+14
       coef(  225) =    -0.39226831913542656250E+14
       coef(  226) =     0.18558762226039976562E+14
       coef(  227) =     0.46607651821358062500E+14
       coef(  228) =     0.24643939663311445312E+14
       coef(  229) =     0.15966973247338005859E+14
       coef(  230) =    -0.11265718127148119141E+14
       coef(  231) =     0.11606861250640859375E+14
       coef(  232) =     0.23165446489078218750E+14
       coef(  233) =     0.15394056926538300781E+14
       coef(  234) =    -0.59359957210302514648E+12
       coef(  235) =     0.35648520125434667969E+13
       coef(  236) =     0.43424950135269091797E+13
       coef(  237) =    -0.37700380979296437500E+14
       coef(  238) =     0.57641899202718742188E+14
       coef(  239) =    -0.16197631547869574219E+14
       coef(  240) =    -0.25639806330298570312E+14
       coef(  241) =    -0.22984352472133906250E+14
       coef(  242) =    -0.67905008908051093750E+13
       coef(  243) =     0.68920254281126031250E+14
       coef(  244) =    -0.57622832839478726562E+14
       coef(  245) =    -0.52416822319136023438E+14
       coef(  246) =    -0.21511010345618242188E+14
       coef(  247) =    -0.69221614323222480469E+13
       coef(  248) =    -0.63196405069364617188E+14
       coef(  249) =    -0.29266131218549710938E+14
       coef(  250) =    -0.79313938324833798828E+13
       coef(  251) =    -0.11505593594258820312E+14
       coef(  252) =    -0.60337391040580253906E+13
       coef(  253) =    -0.47686976122860507812E+14
       coef(  254) =    -0.29606492137280519531E+14
       coef(  255) =    -0.10872349627543511719E+14
       coef(  256) =    -0.26911733534239246094E+14
       coef(  257) =    -0.11017704349830296875E+14
       coef(  258) =    -0.94479820844239980469E+13
       coef(  259) =    -0.60714524229112472534E+10
       coef(  260) =     0.88880342040296234131E+11
       coef(  261) =    -0.36724381530371879883E+13
       coef(  262) =     0.14414920979074572266E+14
       coef(  263) =    -0.24952810358940031250E+14
       coef(  264) =     0.21985236607192468750E+14
       coef(  265) =     0.10053369926770132812E+14
       coef(  266) =    -0.18617409245204640625E+14
       coef(  267) =     0.30934626984900166016E+13
       coef(  268) =     0.61180891527595253906E+13
       coef(  269) =    -0.21789873066353496094E+14
       coef(  270) =     0.15000095570922386719E+14
       coef(  271) =    -0.78945453665746425781E+13
       coef(  272) =    -0.15808360394845005859E+14
       coef(  273) =    -0.11307941558120400391E+14
       coef(  274) =    -0.31183389940745710938E+14
       coef(  275) =     0.84411571209770296875E+14
       coef(  276) =     0.63790438395961109375E+14
       coef(  277) =     0.24880518607867246094E+14
       coef(  278) =     0.72208095255287382812E+13
       coef(  279) =     0.68107887727116382812E+14
       coef(  280) =     0.32328853112085875000E+14
       coef(  281) =     0.12889600625494199219E+14
       coef(  282) =     0.11082159270236582031E+14
       coef(  283) =    -0.25434182724971410156E+14
       coef(  284) =     0.93894652280598765625E+14
       coef(  285) =    -0.75032834990037156250E+14
       coef(  286) =    -0.68354080054360531250E+14
       coef(  287) =    -0.40117003326744367188E+14
       coef(  288) =    -0.18645851048013460938E+14
       coef(  289) =    -0.80769153711131015625E+14
       coef(  290) =    -0.49994827288559656250E+14
       coef(  291) =    -0.20857908597097699219E+14
       coef(  292) =    -0.83625788279082705078E+13
       coef(  293) =    -0.14729182209412373047E+14
       coef(  294) =    -0.35657701853048115234E+13
       coef(  295) =    -0.76087503154026234375E+14
       coef(  296) =    -0.28383925350705507812E+14
       coef(  297) =    -0.88560371983973105469E+13
       coef(  298) =    -0.81929347498476992188E+13
       coef(  299) =     0.42325324674798296875E+14
       coef(  300) =    -0.75828917333559640625E+14
       coef(  301) =    -0.71035774676882328125E+14
       coef(  302) =    -0.38903861288897875000E+14
       coef(  303) =    -0.17676444118636429688E+14
       coef(  304) =    -0.76337487467105218750E+14
       coef(  305) =    -0.30159962264604195312E+14
       coef(  306) =    -0.10520964299065785156E+14
       coef(  307) =    -0.80950951242620693359E+13
       coef(  308) =    -0.34541262049753976562E+14
       coef(  309) =    -0.10148459189973273438E+14
       coef(  310) =    -0.11043284803066744141E+14
       coef(  311) =    -0.14028054237088826994E+05
       coef(  312) =     0.34021215241761151701E+07
       coef(  313) =    -0.29615267429545253515E+09
       coef(  314) =     0.14726038113801223755E+11
       coef(  315) =    -0.30587156572881274414E+12
       coef(  316) =     0.22294060135319941406E+13
       coef(  317) =    -0.78917444083963408203E+13
       coef(  318) =     0.13031217054673433594E+14
       coef(  319) =    -0.67492143006973359375E+13
       coef(  320) =     0.63841430434388481081E+08
       coef(  321) =    -0.11356877745293428421E+11
       coef(  322) =     0.58544912570439416504E+12
       coef(  323) =    -0.20958543098300981445E+13
       coef(  324) =    -0.40030405063537890625E+12
       coef(  325) =     0.24267555398410750000E+14
       coef(  326) =    -0.63524837721014718750E+14
       coef(  327) =     0.35124744160597785156E+14
       coef(  328) =    -0.26316505764455032349E+12
       coef(  329) =     0.13592156187238439941E+13
       coef(  330) =    -0.15528015557717556641E+14
       coef(  331) =     0.26688636460396683594E+14
       coef(  332) =     0.30871980420893609375E+14
       coef(  333) =    -0.24298703592652582031E+14
       coef(  334) =     0.83811529537173652344E+13
       coef(  335) =     0.79239881679641250000E+13
       coef(  336) =    -0.17972093355626828125E+14
       coef(  337) =    -0.69855675601095625000E+13
       coef(  338) =    -0.23665931905170578125E+14
       coef(  339) =    -0.27105124873071335938E+14
       coef(  340) =    -0.94177505117498828125E+13
       coef(  341) =     0.24824707863426863281E+14
       coef(  342) =    -0.15136459443580390625E+14
       coef(  343) =    -0.80041437143085546875E+13
       coef(  344) =    -0.24727610242401015625E+13
       coef(  345) =    -0.29868109487718183594E+14
       coef(  346) =    -0.11854256590864792969E+14
       coef(  347) =     0.38881744387422380447E+10
       coef(  348) =    -0.25775048573638946533E+12
       coef(  349) =    -0.44308086607930263672E+13
       coef(  350) =     0.27190474418009476562E+14
       coef(  351) =    -0.56798203429101359375E+14
       coef(  352) =     0.36405825328252062500E+14
       coef(  353) =     0.26357810121045273438E+13
       coef(  354) =    -0.11728200637970320312E+14
       coef(  355) =     0.39355901461688984375E+13
       coef(  356) =    -0.20073722272049375000E+13
       coef(  357) =     0.18186983387122023438E+14
       coef(  358) =    -0.59280877985707250000E+14
       coef(  359) =    -0.68736741030430097656E+13
       coef(  360) =    -0.87198550644583593750E+12
       coef(  361) =    -0.81987146328468359375E+12
       coef(  362) =    -0.38775902547736210938E+14
       coef(  363) =     0.78493325830041421875E+14
       coef(  364) =     0.88783717877549589844E+13
       coef(  365) =    -0.15726921049993203125E+13
       coef(  366) =    -0.16943094659917917480E+13
       coef(  367) =     0.50548235561022187500E+14
       coef(  368) =     0.75420447193561523438E+13
       coef(  369) =     0.52806902184368603516E+12
       coef(  370) =    -0.10011374845951665039E+13
       coef(  371) =    -0.89408979358014746094E+13
       coef(  372) =     0.28129632467699382812E+14
       coef(  373) =     0.18117520184676574219E+14
       coef(  374) =    -0.23567782665192644531E+14
       coef(  375) =    -0.13081984584424839844E+14
       coef(  376) =    -0.59687048616690468750E+13
       coef(  377) =    -0.67022323885808300781E+12
       coef(  378) =    -0.63058165093885957031E+13
       coef(  379) =    -0.12823381162618216797E+14
       coef(  380) =    -0.67429204573333251953E+13
       coef(  381) =    -0.74459181457501806641E+13
       coef(  382) =    -0.53785525031789062500E+13
       coef(  383) =    -0.22562844231915382812E+14
       coef(  384) =    -0.13717815813906398438E+14
       coef(  385) =    -0.72431363468112265625E+13
       coef(  386) =    -0.64166042066104824219E+13
       coef(  387) =     0.19450884357898718262E+12
       coef(  388) =    -0.21779294503238430176E+13
       coef(  389) =     0.19239956019950585938E+14
       coef(  390) =    -0.63419959479118304688E+14
       coef(  391) =     0.50801981050817406250E+14
       coef(  392) =     0.38141140864581328125E+14
       coef(  393) =    -0.46116409157757812500E+13
       coef(  394) =    -0.17978648290780585938E+14
       coef(  395) =     0.33356023103213271484E+13
       coef(  396) =    -0.78261125264975968750E+14
       coef(  397) =     0.58727226060571484375E+14
       coef(  398) =     0.51773472923423367188E+14
       coef(  399) =     0.20261309310550464844E+14
       coef(  400) =     0.22774336701297519531E+13
       coef(  401) =    -0.36675735525987843750E+14
       coef(  402) =     0.28000127053841421875E+14
       coef(  403) =     0.17943061312948125000E+14
       coef(  404) =     0.62590127943417968750E+13
       coef(  405) =     0.16525428859810281250E+14
       coef(  406) =     0.61537514040247910156E+13
       coef(  407) =     0.59706691442821621094E+13
       coef(  408) =     0.16688698279933737500E+15
       coef(  409) =     0.65431695680093546875E+14
       coef(  410) =     0.17982755535367136719E+14
       coef(  411) =     0.32432525094871113281E+13
       coef(  412) =     0.58271004271006421875E+14
       coef(  413) =     0.14340098925646791016E+14
       coef(  414) =     0.28657154311447148438E+13
       coef(  415) =     0.27315300643625034180E+13
       coef(  416) =     0.11304051757053531250E+14
       coef(  417) =     0.10571333618651945801E+13
       coef(  418) =    -0.19950035167026671875E+14
       coef(  419) =     0.42893123600732265625E+14
       coef(  420) =     0.13654977147097597656E+14
       coef(  421) =     0.23545721305523989258E+13
       coef(  422) =     0.98855787566464570312E+13
       coef(  423) =     0.14661162177213537598E+13
       coef(  424) =     0.43525752024360131836E+12
       coef(  425) =    -0.10670567199179180898E+07
       coef(  426) =    -0.14944397935084164143E+08
       coef(  427) =     0.16271463438631825447E+10
       coef(  428) =    -0.22987550818431457520E+12
       coef(  429) =     0.46141058630239326172E+13
       coef(  430) =    -0.28959262802176066406E+14
       coef(  431) =     0.90885927126742656250E+14
       coef(  432) =    -0.13850186786811478125E+15
       coef(  433) =     0.62196028438374015625E+14
       coef(  434) =     0.85141677767883157730E+09
       coef(  435) =     0.18460932493517562866E+12
       coef(  436) =    -0.80015873650899853516E+13
       coef(  437) =     0.30396245139663382812E+14
       coef(  438) =    -0.38396048023735796875E+14
       coef(  439) =    -0.40602335072038062500E+14
       coef(  440) =     0.97786727698597750000E+14
       coef(  441) =     0.12107061251237567188E+15
       coef(  442) =     0.27508992167249130859E+13
       coef(  443) =    -0.10940583377402501953E+14
       coef(  444) =     0.78548648901869015625E+14
       coef(  445) =    -0.64525937147133843750E+14
       coef(  446) =    -0.14127798177276209375E+15
       coef(  447) =    -0.85468823898213890625E+14
       coef(  448) =    -0.30818636457232613281E+14
       coef(  449) =    -0.36669847361688851562E+14
       coef(  450) =     0.43514031517506703125E+14
       coef(  451) =    -0.68922294526825595703E+13
       coef(  452) =    -0.40810185066056265625E+14
       coef(  453) =    -0.29023844346677421875E+14
       coef(  454) =     0.19196375709998644531E+14
       coef(  455) =    -0.10192014089705925781E+14
       coef(  456) =    -0.12672085909126912109E+14
       coef(  457) =    -0.12364972904854425781E+14
       coef(  458) =    -0.13861183448496722412E+12
       coef(  459) =     0.70355826775237509766E+13
       coef(  460) =     0.33307260247424228516E+13
       coef(  461) =    -0.61880066774022328125E+14
       coef(  462) =     0.45400825799728570312E+14
       coef(  463) =     0.10734378729792589062E+15
       coef(  464) =     0.86026361764678171875E+14
       coef(  465) =     0.44186360090936937500E+14
       coef(  466) =    -0.36445140746742320312E+14
       coef(  467) =     0.49650481581232492188E+14
       coef(  468) =     0.50432602823700390625E+12
       coef(  469) =    -0.56495249456341960938E+14
       coef(  470) =    -0.28310725182309011719E+14
       coef(  471) =    -0.78111431182453281250E+13
       coef(  472) =    -0.62043635087444375000E+13
       coef(  473) =     0.26560677646368179688E+14
       coef(  474) =    -0.59331873450107929688E+13
       coef(  475) =    -0.81750077637378496094E+13
       coef(  476) =     0.18130679068169171875E+14
       coef(  477) =     0.13269246736842851562E+13
       coef(  478) =     0.29199267882316050781E+14
       coef(  479) =    -0.76519380858789578125E+14
       coef(  480) =    -0.40702748456591171875E+14
       coef(  481) =    -0.31998331743998710938E+14
       coef(  482) =    -0.14620796445176367188E+14
       coef(  483) =    -0.49507039857902640625E+14
       coef(  484) =    -0.13760319522868345703E+14
       coef(  485) =    -0.76225627125369658203E+13
       coef(  486) =    -0.30123938515253266602E+13
       coef(  487) =    -0.24192402483023648438E+14
       coef(  488) =    -0.67458108244832382812E+13
       coef(  489) =    -0.21782865528462568359E+13
       coef(  490) =    -0.42753758592346689453E+13
       coef(  491) =     0.41380142173316093750E+14
       coef(  492) =     0.30749789978597699219E+14
       coef(  493) =     0.61115792771461562500E+14
       coef(  494) =     0.50595585135455664062E+14
       coef(  495) =     0.25169710194699597656E+14
       coef(  496) =     0.80571868496147546875E+14
       coef(  497) =    -0.11154087567623615625E+15
       coef(  498) =    -0.32151153619256718750E+14
       coef(  499) =    -0.42338414558391479492E+13
       coef(  500) =     0.17756943746637749023E+13
       coef(  501) =    -0.52641103513098609375E+14
       coef(  502) =    -0.68109445726264785156E+13
       coef(  503) =    -0.59331083530734326172E+12
       coef(  504) =     0.98996706018458813477E+12
       coef(  505) =    -0.18712037254233631250E+15
       coef(  506) =    -0.37791713485122437500E+14
       coef(  507) =    -0.67474329327994335938E+13
       coef(  508) =    -0.18135877851048044434E+13
       coef(  509) =    -0.51701488081943886719E+13
       coef(  510) =    -0.52334261697887792969E+12
       coef(  511) =    -0.53923086361358740234E+12
       coef(  512) =    -0.10889573223619476562E+15
       coef(  513) =    -0.13991289743052925781E+14
       coef(  514) =    -0.18259966356190268555E+13
       coef(  515) =    -0.14365043096666674805E+13
       coef(  516) =     0.47755730408008731902E+08
       coef(  517) =    -0.20455694440961506367E+10
       coef(  518) =     0.77282874796423233032E+11
       coef(  519) =     0.23337360953217529297E+11
       coef(  520) =    -0.10860141802297824219E+14
       coef(  521) =     0.41138736094614218750E+14
       coef(  522) =    -0.44212408111175289062E+14
       coef(  523) =    -0.18689048185468437500E+13
       coef(  524) =     0.11187461857592670312E+15
       coef(  525) =    -0.13420595621863845825E+11
       coef(  526) =    -0.10934642557188989258E+13
       coef(  527) =     0.30999279750876921875E+14
       coef(  528) =    -0.89582305091993500000E+14
       coef(  529) =     0.84778922383620796875E+14
       coef(  530) =    -0.64374358115922109375E+14
       coef(  531) =    -0.57992637789196757812E+13
       coef(  532) =     0.28441833871876531250E+14
       coef(  533) =    -0.57587733002667548828E+13
       coef(  534) =     0.25811657957213816406E+14
       coef(  535) =    -0.55887378863376117188E+14
       coef(  536) =     0.85898936732914609375E+14
       coef(  537) =     0.10043547571020212891E+14
       coef(  538) =    -0.84792693283213789062E+13
       coef(  539) =    -0.27126123263017734375E+14
       coef(  540) =     0.33924013848234867188E+14
       coef(  541) =     0.47366069874190468750E+14
       coef(  542) =     0.11325766151314962891E+14
       coef(  543) =     0.25976828163309933594E+14
       coef(  544) =     0.12507215054109779297E+14
       coef(  545) =     0.26972184067014111328E+12
       coef(  546) =    -0.20422153576103367188E+14
       coef(  547) =    -0.18734266518147093750E+14
       coef(  548) =    -0.10804872518292605469E+14
       coef(  549) =    -0.58581636431400312500E+13
       coef(  550) =    -0.19652405307044921875E+12
       coef(  551) =     0.11240645037738835938E+14
       coef(  552) =     0.79103107095527000000E+14
       coef(  553) =     0.62230952984058671875E+14
       coef(  554) =     0.55048290024466968750E+14
       coef(  555) =     0.82753597565474580078E+13
       coef(  556) =    -0.35245434881090205078E+13
       coef(  557) =     0.31507140633588765625E+14
       coef(  558) =     0.30745887966332750000E+14
       coef(  559) =     0.85652605689225166016E+13
       coef(  560) =     0.14238055491097644531E+14
       coef(  561) =    -0.73805569250166015625E+12
       coef(  562) =    -0.27707033252988902344E+14
       coef(  563) =    -0.21425588199846218262E+13
       coef(  564) =    -0.42319492774972060547E+13
       coef(  565) =    -0.16684892716028554688E+14
       coef(  566) =    -0.16748756817575239258E+13
       coef(  567) =    -0.82235082939149804688E+13
       coef(  568) =     0.11440808094053759766E+14
       coef(  569) =    -0.36688939857024851562E+14
       coef(  570) =    -0.81867521905216968750E+14
       coef(  571) =    -0.44361673779863312500E+14
       coef(  572) =    -0.14204105794580429688E+14
       coef(  573) =     0.24946514370884765625E+12
       coef(  574) =     0.11473350710098154688E+15
       coef(  575) =    -0.25597053270095839844E+14
       coef(  576) =    -0.12915468914771792969E+14
       coef(  577) =    -0.48829195850642519531E+13
       coef(  578) =    -0.12483092933781878906E+14
       coef(  579) =    -0.11393970581716484375E+13
       coef(  580) =    -0.32585315126155546875E+14
       coef(  581) =    -0.13979189920373886719E+14
       coef(  582) =    -0.40452895922070854492E+13
       coef(  583) =    -0.29252716396333642578E+13
       coef(  584) =    -0.27188169306003492188E+14
       coef(  585) =    -0.47185664749659140625E+13
       coef(  586) =    -0.54700474590531980991E+09
       coef(  587) =     0.20925741907972488403E+11
       coef(  588) =    -0.63947510067948217773E+12
       coef(  589) =     0.34229218323343369141E+13
       coef(  590) =     0.17506280764426802734E+14
       coef(  591) =    -0.37873001259533687500E+14
       coef(  592) =    -0.50815534460415023438E+14
       coef(  593) =    -0.55830216814496140625E+14
       coef(  594) =    -0.19278729169402785156E+14
       coef(  595) =     0.29705748576538262939E+12
       coef(  596) =    -0.36978168444904394531E+12
       coef(  597) =    -0.36610016001490367188E+14
       coef(  598) =     0.75255110417005781250E+14
       coef(  599) =     0.80550172490578093750E+14
       coef(  600) =    -0.43775867676102625000E+14
       coef(  601) =    -0.29247852873429171875E+14
       coef(  602) =    -0.12473062852556265625E+14
       coef(  603) =    -0.16004185031395546875E+14
       coef(  604) =    -0.59232165202511242188E+14
       coef(  605) =     0.48698510001278914062E+14
       coef(  606) =     0.17545275904309296875E+14
       coef(  607) =     0.24895123481936804688E+14
       coef(  608) =     0.25313318031336789062E+14
       coef(  609) =     0.30495879240047011719E+14
       coef(  610) =     0.17811135259210636719E+14
       coef(  611) =    -0.18927555785439526367E+13
       coef(  612) =     0.53089081824400726562E+14
       coef(  613) =     0.14624755508985765625E+14
       coef(  614) =     0.10983458859500989062E+15
       coef(  615) =     0.40060226216778703125E+14
       coef(  616) =     0.14218053471358828125E+13
       coef(  617) =    -0.17012206041468096875E+15
       coef(  618) =    -0.96920284735088281250E+13
       coef(  619) =     0.43532304242844726562E+14
       coef(  620) =     0.18676948401210898438E+14
       coef(  621) =     0.12873815345885697266E+14
       coef(  622) =     0.17606184372495273438E+14
       coef(  623) =     0.41099828189721304688E+14
       coef(  624) =     0.14443539296680156250E+14
       coef(  625) =     0.11756204461423146484E+14
       coef(  626) =     0.32407568101936420898E+13
       coef(  627) =    -0.15824574637667824219E+14
       coef(  628) =    -0.47034453025390125000E+14
       coef(  629) =     0.10075945417724728125E+15
       coef(  630) =     0.47365404856736125000E+14
       coef(  631) =     0.92803405950756953125E+13
       coef(  632) =     0.49065945186998726562E+14
       coef(  633) =     0.22258369539283078125E+14
       coef(  634) =     0.92437939011782812500E+13
       coef(  635) =     0.54530268098291992188E+13
       coef(  636) =     0.13463989210253515625E+14
       coef(  637) =     0.50059105843344013672E+13
       coef(  638) =     0.10248209801942011719E+13
       coef(  639) =     0.31288836617353243828E+10
       coef(  640) =    -0.14636356190816149902E+12
       coef(  641) =     0.34028132407943256836E+13
       coef(  642) =    -0.23564596728988496094E+14
       coef(  643) =     0.10857064458494246094E+14
       coef(  644) =    -0.12076355808441699219E+14
       coef(  645) =     0.14216934302535475000E+15
       coef(  646) =     0.90833099938381343750E+14
       coef(  647) =    -0.82135246230663037109E+12
       coef(  648) =     0.42587276120228662109E+13
       coef(  649) =    -0.30829175628796406250E+13
       coef(  650) =    -0.32465755867634710938E+14
       coef(  651) =    -0.19151692763190617188E+14
       coef(  652) =    -0.31658126328036593750E+14
       coef(  653) =     0.74509207816174187500E+14
       coef(  654) =    -0.71348512402916757812E+13
       coef(  655) =    -0.51923264725418664062E+14
       coef(  656) =     0.53025516463809375000E+13
       coef(  657) =     0.26439359438265082031E+14
       coef(  658) =     0.10542779399307669922E+14
       coef(  659) =     0.10673178529437646484E+13
       coef(  660) =     0.10251148134170306641E+14
       coef(  661) =    -0.65731524425867820312E+14
       coef(  662) =     0.44050371285506625000E+14
       coef(  663) =     0.19639957122002039062E+14
       coef(  664) =    -0.90417664783136531250E+14
       coef(  665) =    -0.18078052548428980469E+14
       coef(  666) =     0.18033649089888328125E+14
       coef(  667) =     0.31136110019840751953E+13
       coef(  668) =     0.75236016558631781250E+14
       coef(  669) =     0.22754226521333660156E+14
       coef(  670) =     0.21941397369466406250E+14
       coef(  671) =     0.61404156080090859375E+13
       coef(  672) =     0.10898107307758732812E+15
       coef(  673) =     0.49678286764990976562E+14
       coef(  674) =     0.29982383500396335938E+14
       coef(  675) =     0.23939873391485632812E+14
       coef(  676) =     0.16870888516442781250E+14
       coef(  677) =    -0.65661283826072463989E+10
       coef(  678) =     0.32917309085040399170E+12
       coef(  679) =    -0.81494686214636796875E+13
       coef(  680) =     0.58298030126752406250E+14
       coef(  681) =    -0.50684165276970812500E+14
       coef(  682) =    -0.68490291010509523438E+14
       coef(  683) =     0.10359308042261637500E+15
       coef(  684) =     0.31053855770190156250E+13
       coef(  685) =    -0.25947247774721277344E+14
       coef(  686) =     0.42300040301502296875E+14
       coef(  687) =    -0.57406144136148296875E+14
       coef(  688) =    -0.45458530654316945312E+14
       coef(  689) =    -0.25160758164886132812E+14
       coef(  690) =    -0.23976567076517539062E+14
       coef(  691) =    -0.36924214221307945312E+14
       coef(  692) =     0.89613911267183632812E+13
       coef(  693) =    -0.14162277735751214844E+14
       coef(  694) =     0.52770999736575097656E+13
       coef(  695) =    -0.77799180516294250000E+14
       coef(  696) =    -0.61950893987525573730E+12
       coef(  697) =    -0.69771474556780400391E+13
       coef(  698) =    -0.86433829685770273438E+13
       coef(  699) =     0.65647319784080320312E+14
       coef(  700) =    -0.47546311455148994141E+13
       coef(  701) =     0.25988052096113218750E+14
       coef(  702) =     0.63792455934946531250E+14
       coef(  703) =     0.22527179894990625000E+14
       coef(  704) =     0.48409915728468856812E+10
       coef(  705) =    -0.29712481425290686035E+12
       coef(  706) =     0.82176376433453417969E+13
       coef(  707) =    -0.58774439906860218750E+14
       coef(  708) =     0.67522234666562945312E+14
       coef(  709) =     0.19205088714314808594E+14
       coef(  710) =    -0.32123134672632207031E+13
       coef(  711) =     0.16030019083753623047E+14
       coef(  712) =     0.73169806003407531250E+14
       coef(  713) =    -0.13801922741115146484E+14
       coef(  714) =    -0.69627811385917531250E+14
       coef(  715) =    -0.22410795426777515625E+14
       coef(  716) =     0.19861437431167589844E+14
       coef(  717) =     0.22189000751928644531E+14
       coef(  718) =    -0.46105303271808585938E+14
       coef(  719) =     0.23139153730648312500E+14
       coef(  720) =    -0.27927048520300949219E+14
       coef(  721) =     0.23878097874097855469E+14
       coef(  722) =    -0.22420763800461480741E+03
       coef(  723) =    -0.51950266372330523154E+05
       coef(  724) =     0.63857241582316523418E+07
       coef(  725) =    -0.33366955692816811800E+09
       coef(  726) =     0.86039316698797340393E+10
       coef(  727) =    -0.68928484106209152222E+11
       coef(  728) =     0.29875991371062835693E+12
       coef(  729) =    -0.66638544142870031738E+12
       coef(  730) =     0.53835454303419628906E+12
       coef(  731) =     0.24036504928493141197E+07
       coef(  732) =    -0.30937577697204315662E+09
       coef(  733) =     0.28353742150430679321E+07
       coef(  734) =    -0.26699843559204309082E+12
       coef(  735) =     0.24237188582592905273E+13
       coef(  736) =    -0.10696524458865519531E+14
       coef(  737) =     0.24133197490199343750E+14
       coef(  738) =    -0.20637526753872394531E+14
       coef(  739) =     0.95837393434877605438E+10
       coef(  740) =    -0.57635030296458534241E+11
       coef(  741) =     0.19982679337229692383E+13
       coef(  742) =    -0.15928053045556953125E+14
       coef(  743) =     0.66342786188882937500E+14
       coef(  744) =    -0.13841152251251635938E+15
       coef(  745) =     0.11855155979039343750E+15
       coef(  746) =    -0.52191891694110302734E+12
       coef(  747) =     0.14338017600707092285E+13
       coef(  748) =     0.49474475471678154297E+13
       coef(  749) =    -0.18734719837833195312E+14
       coef(  750) =    -0.49705046201531335938E+14
       coef(  751) =     0.22549415979091488281E+14
       coef(  752) =    -0.26849118068165551758E+13
       coef(  753) =     0.12631857335933224609E+14
       coef(  754) =     0.62234984490752531250E+14
       coef(  755) =     0.42790939394368523438E+14
       coef(  756) =     0.29084133341702687500E+14
       coef(  757) =    -0.54928295787027601562E+14
       coef(  758) =    -0.76502376876866933594E+13
       coef(  759) =     0.94334640412035468750E+13
       coef(  760) =    -0.50601941565945595703E+13
       coef(  761) =     0.17604991146590074897E+09
       coef(  762) =     0.50228072404702204466E+09
       coef(  763) =     0.49783105438509765625E+12
       coef(  764) =    -0.14909220367935798340E+13
       coef(  765) =    -0.12786658932945859375E+13
       coef(  766) =     0.16440796776470585938E+14
       coef(  767) =    -0.34091630417193039062E+14
       coef(  768) =     0.13020957952464265625E+14
       coef(  769) =    -0.47239180978991326904E+12
       coef(  770) =    -0.19588860873980947266E+13
       coef(  771) =     0.16976744472439976562E+14
       coef(  772) =    -0.55655311060659312500E+14
       coef(  773) =     0.62385535064576609375E+14
       coef(  774) =     0.14769419482524878906E+14
       coef(  775) =     0.22321328497968632812E+14
       coef(  776) =     0.10020182604846691406E+14
       coef(  777) =    -0.39226831913480093750E+14
       coef(  778) =     0.18558762225959566406E+14
       coef(  779) =     0.46607651821334781250E+14
       coef(  780) =     0.24643939663386988281E+14
       coef(  781) =     0.15966973247295802734E+14
       coef(  782) =    -0.11265718127125117188E+14
       coef(  783) =     0.11606861250521949219E+14
       coef(  784) =     0.23165446489133023438E+14
       coef(  785) =     0.15394056926595917969E+14
       coef(  786) =    -0.59359957205498754883E+12
       coef(  787) =     0.35648520124821870117E+13
       coef(  788) =     0.43424950134728481445E+13
       coef(  789) =    -0.37700380979338968750E+14
       coef(  790) =     0.57641899202721429688E+14
       coef(  791) =    -0.16197631547899441406E+14
       coef(  792) =    -0.25639806330293953125E+14
       coef(  793) =    -0.22984352472198730469E+14
       coef(  794) =    -0.67905008908268281250E+13
       coef(  795) =     0.68920254281120687500E+14
       coef(  796) =    -0.57622832839441679688E+14
       coef(  797) =    -0.52416822319164109375E+14
       coef(  798) =    -0.21511010345650496094E+14
       coef(  799) =    -0.69221614323036611328E+13
       coef(  800) =    -0.63196405069387054688E+14
       coef(  801) =    -0.29266131218528828125E+14
       coef(  802) =    -0.79313938325482089844E+13
       coef(  803) =    -0.11505593594257554688E+14
       coef(  804) =    -0.60337391040600517578E+13
       coef(  805) =    -0.47686976122837117188E+14
       coef(  806) =    -0.29606492137273191406E+14
       coef(  807) =    -0.10872349627564220703E+14
       coef(  808) =    -0.26911733534275195312E+14
       coef(  809) =    -0.11017704349792501953E+14
       coef(  810) =    -0.94479820843836660156E+13
       coef(  811) =    -0.60714523892500944138E+10
       coef(  812) =     0.88880342052271331787E+11
       coef(  813) =    -0.36724381530065434570E+13
       coef(  814) =     0.14414920979063939453E+14
       coef(  815) =    -0.24952810358962722656E+14
       coef(  816) =     0.21985236607223441406E+14
       coef(  817) =     0.10053369926771414062E+14
       coef(  818) =    -0.18617409245191675781E+14
       coef(  819) =     0.30934626984747753906E+13
       coef(  820) =     0.61180891527610195312E+13
       coef(  821) =    -0.21789873066357550781E+14
       coef(  822) =     0.15000095570931796875E+14
       coef(  823) =    -0.78945453665746054688E+13
       coef(  824) =    -0.15808360394826724609E+14
       coef(  825) =    -0.11307941558123414062E+14
       coef(  826) =    -0.31183389940749773438E+14
       coef(  827) =     0.84411571209787468750E+14
       coef(  828) =     0.63790438395965640625E+14
       coef(  829) =     0.24880518607861433594E+14
       coef(  830) =     0.72208095255231445312E+13
       coef(  831) =     0.68107887727113148438E+14
       coef(  832) =     0.32328853112081109375E+14
       coef(  833) =     0.12889600625477560547E+14
       coef(  834) =     0.11082159270227519531E+14
       coef(  835) =    -0.25434182724981175781E+14
       coef(  836) =     0.93894652280609875000E+14
       coef(  837) =    -0.75032834990029625000E+14
       coef(  838) =    -0.68354080054376117188E+14
       coef(  839) =    -0.40117003326747750000E+14
       coef(  840) =    -0.18645851048025519531E+14
       coef(  841) =    -0.80769153711137156250E+14
       coef(  842) =    -0.49994827288564515625E+14
       coef(  843) =    -0.20857908597081882812E+14
       coef(  844) =    -0.83625788279166015625E+13
       coef(  845) =    -0.14729182209410396484E+14
       coef(  846) =    -0.35657701852945415039E+13
       coef(  847) =    -0.76087503154026640625E+14
       coef(  848) =    -0.28383925350705476562E+14
       coef(  849) =    -0.88560371984075273438E+13
       coef(  850) =    -0.81929347498434667969E+13
       coef(  851) =     0.42325324674808828125E+14
       coef(  852) =    -0.75828917333564031250E+14
       coef(  853) =    -0.71035774676889421875E+14
       coef(  854) =    -0.38903861288897875000E+14
       coef(  855) =    -0.17676444118641117188E+14
       coef(  856) =    -0.76337487467109906250E+14
       coef(  857) =    -0.30159962264596058594E+14
       coef(  858) =    -0.10520964299064509766E+14
       coef(  859) =    -0.80950951242693808594E+13
       coef(  860) =    -0.34541262049761781250E+14
       coef(  861) =    -0.10148459189961609375E+14
       coef(  862) =    -0.11043284803060230469E+14
       coef(  863) =     0.47119303046054701554E+05
       coef(  864) =    -0.10541551275568570942E+08
       coef(  865) =     0.84878117173955416679E+09
       coef(  866) =    -0.31236444103822383881E+11
       coef(  867) =     0.50576844826024414062E+12
       coef(  868) =    -0.34481294122457729492E+13
       coef(  869) =     0.10455386325802570312E+14
       coef(  870) =    -0.12088312733114480469E+14
       coef(  871) =     0.12239941353455058594E+13
       coef(  872) =    -0.37508582174566709995E+09
       coef(  873) =     0.30959267336045745850E+11
       coef(  874) =    -0.93178279405905151367E+12
       coef(  875) =     0.37049105816155214844E+13
       coef(  876) =     0.40049348512314135742E+13
       coef(  877) =    -0.41242567249494656250E+14
       coef(  878) =     0.12924004611370019531E+14
       coef(  879) =     0.99240890209049500000E+14
       coef(  880) =     0.30309816693972827148E+12
       coef(  881) =    -0.12512917209554230957E+13
       coef(  882) =     0.26437956816929042969E+14
       coef(  883) =    -0.74489104114218265625E+14
       coef(  884) =     0.47144971314505312500E+14
       coef(  885) =     0.79615882705792519531E+13
       coef(  886) =     0.14674026885182132812E+14
       coef(  887) =    -0.18752094570353976562E+14
       coef(  888) =     0.30371919240069222656E+14
       coef(  889) =     0.66676955272494234375E+14
       coef(  890) =     0.23020048150399433594E+14
       coef(  891) =    -0.14869470987881273438E+14
       coef(  892) =    -0.12260145264157937500E+14
       coef(  893) =    -0.37489922365718531250E+14
       coef(  894) =    -0.66532583230470937500E+13
       coef(  895) =     0.49156206918197636719E+13
       coef(  896) =     0.23884705632097045898E+12
       coef(  897) =    -0.11645396212697605469E+14
       coef(  898) =    -0.14653096625758188477E+13
       coef(  899) =    -0.12305647382661199570E+11
       coef(  900) =     0.68848192866590527344E+12
       coef(  901) =     0.35390275044895615234E+13
       coef(  902) =    -0.49629862529652515625E+14
       coef(  903) =     0.13132724712766260938E+15
       coef(  904) =    -0.75039455671948496094E+13
       coef(  905) =    -0.14855006782326912500E+15
       coef(  906) =    -0.13908038774592048438E+15
       coef(  907) =    -0.97392319047426562500E+13
       coef(  908) =     0.65845399529995898438E+14
       coef(  909) =    -0.13644678955213021875E+15
       coef(  910) =     0.67785278196545867188E+14
       coef(  911) =     0.67449144058038406250E+14
       coef(  912) =     0.71874000305856796875E+13
       coef(  913) =    -0.12296904793096572266E+14
       coef(  914) =    -0.10150798980477222656E+14
       coef(  915) =     0.12689673682750648438E+14
       coef(  916) =     0.35825557495601437500E+14
       coef(  917) =     0.21906134276428687500E+14
       coef(  918) =     0.51240460688965683594E+13
       coef(  919) =     0.16882221958276562500E+13
       coef(  920) =     0.39456642516759584961E+13
       coef(  921) =     0.40373749130404506836E+13
       coef(  922) =    -0.12192957784233007812E+12
       coef(  923) =     0.11569274082195347656E+14
       coef(  924) =    -0.80714931797335875000E+14
       coef(  925) =     0.42600933536134218750E+14
       coef(  926) =     0.59136977986113531250E+14
       coef(  927) =     0.29394731489178363281E+14
       coef(  928) =     0.63698459886821943359E+13
       coef(  929) =     0.93019410629180351562E+13
       coef(  930) =     0.14963617136723082031E+14
       coef(  931) =     0.98621061268859003906E+13
       coef(  932) =     0.43025130343243896484E+13
       coef(  933) =    -0.69245296188451904297E+12
       coef(  934) =    -0.74013696764062182617E+12
       coef(  935) =     0.10242566003699558594E+14
       coef(  936) =     0.17689425415628095703E+13
       coef(  937) =     0.25882803484875585938E+12
       coef(  938) =    -0.12380321561982314453E+13
       coef(  939) =    -0.15944536980247351074E+12
       coef(  940) =    -0.51945511395217480469E+12
       coef(  941) =     0.14037946064130322266E+13
       coef(  942) =     0.97305749365336328125E+13
       coef(  943) =     0.91548302315334707031E+13
       coef(  944) =    -0.11175202422305595312E+15
       coef(  945) =    -0.13254629081055325000E+15
       coef(  946) =    -0.91262218380661171875E+14
       coef(  947) =     0.18854092976609335938E+14
       coef(  948) =    -0.46103544515648679688E+14
       coef(  949) =     0.73697171310182125000E+14
       coef(  950) =     0.49762459289997812500E+14
       coef(  951) =     0.14035130679191210938E+12
       coef(  952) =    -0.16061295873626222656E+14
       coef(  953) =     0.14621181826078601562E+14
       coef(  954) =     0.36924690641713265625E+14
       coef(  955) =     0.20570827367769523438E+14
       coef(  956) =     0.52967710835897392578E+13
       coef(  957) =     0.15385965769568511719E+14
       coef(  958) =     0.60858490249165458984E+13
       coef(  959) =    -0.56616754763374453125E+14
       coef(  960) =     0.15359982890621975000E+15
       coef(  961) =     0.70317678996754531250E+14
       coef(  962) =     0.20073365335475246094E+14
       coef(  963) =     0.11505476235443398438E+13
       coef(  964) =     0.64011617007056304688E+14
       coef(  965) =     0.18064002173567816406E+14
       coef(  966) =     0.43705115058656264648E+13
       coef(  967) =     0.38667471737670429688E+13
       coef(  968) =     0.15821605618010908203E+14
       coef(  969) =     0.28840385204450136719E+13
       coef(  970) =    -0.81856485925582187500E+14
       coef(  971) =     0.29961721663864554688E+14
       coef(  972) =     0.13427414565394453125E+14
       coef(  973) =     0.25854253742607436523E+13
       coef(  974) =     0.86373146379669052734E+13
       coef(  975) =     0.18887191956555034180E+13
       coef(  976) =     0.56800200971033105469E+12
       coef(  977) =     0.53701637360473582521E+06
       coef(  978) =     0.64294278391538411379E+08
       coef(  979) =    -0.44307303758825645447E+10
       coef(  980) =     0.18646712887953179932E+12
       coef(  981) =    -0.21431627903627512207E+13
       coef(  982) =     0.11532075835627615234E+14
       coef(  983) =    -0.32687515166100500000E+14
       coef(  984) =     0.54004123038855335938E+14
       coef(  985) =    -0.42740860907219601562E+14
       coef(  986) =     0.22852734499881172180E+10
       coef(  987) =    -0.25273607017151742554E+12
       coef(  988) =     0.34698026639051992188E+13
       coef(  989) =    -0.14190183868348281250E+14
       coef(  990) =     0.12975752255913691406E+14
       coef(  991) =     0.51275433648910781250E+13
       coef(  992) =     0.30453316562393710938E+14
       coef(  993) =     0.11433233381293773438E+14
       coef(  994) =     0.14303587942838323975E+12
       coef(  995) =    -0.49059438274171484375E+13
       coef(  996) =    -0.14337048254234306641E+14
       coef(  997) =     0.20543408354311933594E+14
       coef(  998) =    -0.16042555426586972656E+14
       coef(  999) =    -0.26148443860693015625E+14
       coef( 1000) =    -0.17963338275828472656E+14
       coef( 1001) =     0.46942379945665421875E+14
       coef( 1002) =    -0.43275646947811976562E+14
       coef( 1003) =     0.29560037196563266602E+13
       coef( 1004) =    -0.34650529179754179688E+13
       coef( 1005) =    -0.84329007176829628906E+13
       coef( 1006) =    -0.58503406180654570312E+14
       coef( 1007) =    -0.17142048685197019531E+14
       coef( 1008) =    -0.38079270403582910156E+13
       coef( 1009) =    -0.73740300789225830078E+13
       coef( 1010) =    -0.55225432409759082794E+10
       coef( 1011) =    -0.19674996863116149902E+12
       coef( 1012) =    -0.67155257695122617188E+13
       coef( 1013) =     0.46034981361242273438E+14
       coef( 1014) =    -0.48459451213260296875E+14
       coef( 1015) =    -0.72659130522292441406E+13
       coef( 1016) =    -0.19749578497683750000E+13
       coef( 1017) =    -0.97023850401833046875E+13
       coef( 1018) =     0.89184449054707666016E+12
       coef( 1019) =    -0.22550752314422246094E+14
       coef( 1020) =    -0.11821965000273037109E+14
       coef( 1021) =    -0.24710312211518796875E+14
       coef( 1022) =    -0.98585487809070234375E+13
       coef( 1023) =    -0.49039929584986513672E+13
       coef( 1024) =    -0.10735978961002675781E+14
       coef( 1025) =     0.80283837447616015625E+13
       coef( 1026) =     0.38441151271382275391E+12
       coef( 1027) =    -0.89165562051921411133E+12
       coef( 1028) =     0.39691402925852900391E+13
       coef( 1029) =     0.51319327507855859375E+12
       coef( 1030) =     0.21627130849872523438E+14
       coef( 1031) =    -0.71998350161179199219E+13
       coef( 1032) =     0.88128742488609160156E+13
       coef( 1033) =    -0.66380380776393041992E+12
       coef( 1034) =    -0.12732764831942104492E+13
       coef( 1035) =    -0.53478140418219199219E+13
       coef( 1036) =     0.18895852011094472656E+13
       coef( 1037) =     0.12715788119999340820E+12
       coef( 1038) =     0.47577971745077807617E+12
       coef( 1039) =    -0.40433241050699306641E+13
       coef( 1040) =    -0.64473225070283544922E+12
       coef( 1041) =     0.76131321040805834961E+12
       coef( 1042) =    -0.15789731593968564453E+13
       coef( 1043) =    -0.46492953889469707031E+13
       coef( 1044) =    -0.19438409400750976562E+14
       coef( 1045) =    -0.16462809434617603516E+14
       coef( 1046) =    -0.11446134648780958984E+14
       coef( 1047) =    -0.11345801749402957031E+14
       coef( 1048) =    -0.25518280079346453125E+14
       coef( 1049) =     0.91857407720733046875E+13
       coef( 1050) =     0.61577609689761396484E+13
       coef( 1051) =     0.44057340222892578125E+11
       coef( 1052) =    -0.17996442410991669922E+13
       coef( 1053) =     0.13320271774975421875E+14
       coef( 1054) =     0.73426169296255898438E+13
       coef( 1055) =     0.22428978304674130859E+13
       coef( 1056) =     0.30872164938728466797E+13
       coef( 1057) =     0.49915035693668789062E+14
       coef( 1058) =     0.57632313538126062500E+14
       coef( 1059) =     0.17469957589151623047E+14
       coef( 1060) =     0.33319194421761679688E+13
       coef( 1061) =     0.22891462010676957031E+14
       coef( 1062) =     0.57083061623259570312E+13
       coef( 1063) =     0.66073050267828652344E+13
       coef( 1064) =     0.26363726065523312500E+14
       coef( 1065) =     0.20683556464372707031E+14
       coef( 1066) =     0.56630867690622753906E+13
       coef( 1067) =     0.63844823001807675781E+13
       coef( 1068) =    -0.10836448530801504850E+09
       coef( 1069) =     0.51949221871807804108E+10
       coef( 1070) =    -0.85060360800480590820E+11
       coef( 1071) =    -0.14349950884808930664E+13
       coef( 1072) =     0.11564160027109660156E+14
       coef( 1073) =    -0.14509943259160671875E+14
       coef( 1074) =    -0.63912006112720382812E+14
       coef( 1075) =     0.54700220230511937500E+14
       coef( 1076) =     0.10650100230698082031E+14
       coef( 1077) =     0.34030622197459831238E+11
       coef( 1078) =     0.35920057422233242188E+13
       coef( 1079) =    -0.26328816799539898438E+14
       coef( 1080) =     0.70103274736783906250E+14
       coef( 1081) =     0.69488780622105968750E+14
       coef( 1082) =    -0.20747602891073210938E+14
       coef( 1083) =     0.29494428272803164062E+13
       coef( 1084) =     0.56028312031585283203E+13
       coef( 1085) =    -0.78466978591002324219E+13
       coef( 1086) =     0.92436695767507773438E+13
       coef( 1087) =    -0.99462710900734671875E+14
       coef( 1088) =     0.85474249128694326172E+13
       coef( 1089) =    -0.30565866596955605469E+13
       coef( 1090) =    -0.63942331697216113281E+13
       coef( 1091) =    -0.17157792988626867188E+14
       coef( 1092) =    -0.45553474316600289062E+14
       coef( 1093) =     0.56629418625714941406E+12
       coef( 1094) =     0.11089069540516259766E+13
       coef( 1095) =    -0.27306222590224976562E+14
       coef( 1096) =    -0.52531268128315615234E+13
       coef( 1097) =    -0.82591881350707165527E+12
       coef( 1098) =     0.90929331500833032227E+12
       coef( 1099) =     0.89729675692250921875E+14
       coef( 1100) =    -0.33816768953478886719E+14
       coef( 1101) =    -0.10946130999228364062E+15
       coef( 1102) =    -0.55056427438224796875E+14
       coef( 1103) =    -0.14340351276168433594E+14
       coef( 1104) =     0.67131002115360791016E+13
       coef( 1105) =    -0.29096953315260992188E+14
       coef( 1106) =    -0.22728168589461734375E+14
       coef( 1107) =    -0.26398059455318007812E+14
       coef( 1108) =    -0.13319871068350609375E+14
       coef( 1109) =    -0.18613119964998246094E+14
       coef( 1110) =    -0.17114830894659912109E+13
       coef( 1111) =    -0.21561191399561801758E+13
       coef( 1112) =     0.13497308112099814453E+13
       coef( 1113) =    -0.89799400139715000000E+14
       coef( 1114) =    -0.45228581993548218750E+14
       coef( 1115) =    -0.11803211902753662109E+14
       coef( 1116) =    -0.65890686107601240234E+13
       coef( 1117) =    -0.17633923930914839844E+14
       coef( 1118) =    -0.32736460224683925781E+13
       coef( 1119) =    -0.63261634423015917969E+13
       coef( 1120) =    -0.67335729431858154297E+13
       coef( 1121) =     0.45959931831488710938E+14
       coef( 1122) =    -0.13374659943808476562E+13
       coef( 1123) =    -0.65450212512601132812E+14
       coef( 1124) =    -0.47659799203762695312E+14
       coef( 1125) =    -0.18142413667967648438E+14
       coef( 1126) =     0.94575232025752832031E+13
       coef( 1127) =    -0.16169652224849158203E+14
       coef( 1128) =    -0.16527646371450009766E+14
       coef( 1129) =    -0.93314246287105390625E+13
       coef( 1130) =    -0.57077391848787421875E+13
       coef( 1131) =    -0.18480126835225410156E+13
       coef( 1132) =    -0.25721613829721416016E+13
       coef( 1133) =     0.23760870017854223633E+13
       coef( 1134) =    -0.10695401318593688965E+13
       coef( 1135) =     0.19058475502563647461E+13
       coef( 1136) =     0.11035830914746518555E+13
       coef( 1137) =     0.24882636843509033203E+13
       coef( 1138) =     0.10390180032105814219E+10
       coef( 1139) =    -0.58357077551278182983E+11
       coef( 1140) =     0.97756639504330615234E+12
       coef( 1141) =     0.72128840484013085938E+13
       coef( 1142) =    -0.52481573035519539062E+14
       coef( 1143) =     0.85399731273359578125E+14
       coef( 1144) =     0.10940305787063160156E+14
       coef( 1145) =     0.62242533799079734375E+14
       coef( 1146) =     0.44425337932781867188E+14
       coef( 1147) =    -0.79291684411215722656E+12
       coef( 1148) =    -0.57489772524505292969E+13
       coef( 1149) =    -0.19374210224102273438E+14
       coef( 1150) =    -0.63807468422767093750E+14
       coef( 1151) =    -0.48874940311027609375E+14
       coef( 1152) =    -0.60333966623153062500E+14
       coef( 1153) =    -0.20364505632492242188E+14
       coef( 1154) =     0.51634650078883789062E+14
       coef( 1155) =     0.70778619888589812500E+14
       coef( 1156) =    -0.53627706914519921875E+14
       coef( 1157) =     0.46742590507646679688E+12
       coef( 1158) =    -0.80339934685609741211E+12
       coef( 1159) =     0.63451867132702539062E+14
       coef( 1160) =     0.47061159844825527344E+13
       coef( 1161) =     0.82500149130695839844E+13
       coef( 1162) =     0.24896155932631152344E+12
       coef( 1163) =     0.82362504372485156250E+13
       coef( 1164) =    -0.66057045623644187500E+14
       coef( 1165) =    -0.56312694551184140625E+14
       coef( 1166) =    -0.38173056707279359375E+14
       coef( 1167) =    -0.54957323626716078125E+14
       coef( 1168) =    -0.31500731403093191406E+14
       coef( 1169) =    -0.40941120883670595703E+13
       coef( 1170) =    -0.35196897575740210938E+14
       coef( 1171) =    -0.12532187374293640625E+14
       coef( 1172) =    -0.91246546542991914062E+13
       coef( 1173) =    -0.13625333757247324219E+14
       coef( 1174) =    -0.97834450240658154297E+12
       coef( 1175) =     0.30584776205886945312E+14
       coef( 1176) =    -0.64149962078108808594E+13
       coef( 1177) =    -0.24616103062588437500E+13
       coef( 1178) =    -0.43050854544920561523E+13
       coef( 1179) =     0.30483834015973310547E+13
       coef( 1180) =    -0.39658025350166914062E+14
       coef( 1181) =     0.38174625968356968750E+14
       coef( 1182) =    -0.93313658907868398438E+13
       coef( 1183) =    -0.16710960926576501953E+14
       coef( 1184) =    -0.12512948459473054688E+14
       coef( 1185) =    -0.39244099848858769531E+13
       coef( 1186) =    -0.45891956034088339844E+13
       coef( 1187) =    -0.21877015089203447266E+13
       coef( 1188) =    -0.45103514611310488281E+13
       coef( 1189) =    -0.50642366289316845703E+12
       coef( 1190) =    -0.11217223770530124512E+13
       coef( 1191) =    -0.44455802103393917084E+10
       coef( 1192) =     0.43081409267232342529E+12
       coef( 1193) =    -0.86644707459532880859E+13
       coef( 1194) =     0.49117764672749287109E+13
       coef( 1195) =     0.59347930722204335938E+14
       coef( 1196) =    -0.40015227819391265625E+14
       coef( 1197) =    -0.50152403798047203125E+14
       coef( 1198) =     0.19864493826854580078E+13
       coef( 1199) =     0.20100572903167399902E+13
       coef( 1200) =     0.30027768169507636719E+14
       coef( 1201) =     0.99065862206378359375E+14
       coef( 1202) =     0.50524352915927625000E+14
       coef( 1203) =    -0.88736332121145937500E+13
       coef( 1204) =    -0.29559596476593363281E+14
       coef( 1205) =    -0.15444085748622725000E+15
       coef( 1206) =    -0.24324922156121289062E+14
       coef( 1207) =    -0.33250193354784148438E+14
       coef( 1208) =    -0.29257686658005043945E+13
       coef( 1209) =     0.21148950293945062500E+14
       coef( 1210) =     0.37537168589091845703E+13
       coef( 1211) =    -0.12401670741614988281E+14
       coef( 1212) =     0.13826213027771287500E+15
       coef( 1213) =    -0.63648278763645136719E+13
       coef( 1214) =     0.83238812435091699219E+13
       coef( 1215) =    -0.90683653641717871094E+13
       coef( 1216) =     0.70049570318013296875E+14
       coef( 1217) =    -0.85230222402177392578E+13
       coef( 1218) =    -0.14706123757854062500E+13
       coef( 1219) =    -0.54559125873362949219E+13
       coef( 1220) =     0.76711752374154125000E+14
       coef( 1221) =     0.11292283288594375000E+14
       coef( 1222) =     0.80010995147486093750E+13
       coef( 1223) =    -0.28972180074107359375E+14
       coef( 1224) =     0.46205590988540875000E+14
       coef( 1225) =     0.14475851017462548828E+14
       coef( 1226) =    -0.87224598496226894531E+13
       coef( 1227) =     0.33522680602469765625E+13
       coef( 1228) =    -0.43864091698219104004E+12
       coef( 1229) =    -0.32817248533276729584E+10
       coef( 1230) =    -0.84462524532937023926E+12
       coef( 1231) =     0.27837974920524742188E+14
       coef( 1232) =    -0.60561061781641312500E+14
       coef( 1233) =     0.57479705119389531250E+13
       coef( 1234) =    -0.68546963678442812500E+14
       coef( 1235) =    -0.36050527022278875000E+14
       coef( 1236) =    -0.77868352756604609375E+13
       coef( 1237) =    -0.66954595018029578125E+14
       coef( 1238) =    -0.43862122705519671875E+14
       coef( 1239) =     0.73339651582112109375E+13
       coef( 1240) =    -0.75101801743839804688E+13
       coef( 1241) =     0.45866926675241953125E+14
       coef( 1242) =    -0.11956269581988033203E+14
       coef( 1243) =    -0.17560292444682853516E+14
       coef( 1244) =     0.63678018698222363281E+13
       coef( 1245) =    -0.17758617035910394531E+14
       coef( 1246) =    -0.15193778007887167969E+14
       coef( 1247) =    -0.58099837156778304688E+14
       coef( 1248) =    -0.45194756190423613281E+13
       coef( 1249) =     0.34940642889154796875E+14
       coef( 1250) =    -0.72289149768106093750E+13
       coef( 1251) =     0.50303992653750695312E+14
       coef( 1252) =     0.82591490208804658203E+13
       coef( 1253) =    -0.40488762838213164062E+14
       coef( 1254) =     0.18317589543321433594E+14
       coef( 1255) =    -0.92607280492681386719E+13
       coef( 1256) =     0.46489773685702674866E+11
       coef( 1257) =    -0.38196291380582141113E+12
       coef( 1258) =    -0.35381476967243734375E+14
       coef( 1259) =     0.75251022168518546875E+14
       coef( 1260) =     0.18913446656654566406E+14
       coef( 1261) =    -0.35827954068278304688E+14
       coef( 1262) =     0.25779637741305890625E+14
       coef( 1263) =     0.34987861511586789062E+14
       coef( 1264) =    -0.53169652614772492188E+14
       coef( 1265) =    -0.10625949821277691406E+14
       coef( 1266) =     0.14912926816287240625E+15
       coef( 1267) =     0.11257359126014130859E+14
       coef( 1268) =     0.36836883923678179688E+14
       coef( 1269) =    -0.78509498092325625000E+14
       coef( 1270) =    -0.63899705373681734375E+14
       coef( 1271) =     0.67989953511523720703E+13
       coef( 1272) =     0.51021167320464824219E+13
       coef( 1273) =    -0.34468057213143566406E+14
       coef( 1274) =    -0.55339320865870971680E+11
       coef( 1275) =     0.14132324903468344727E+13
       coef( 1276) =     0.14443915055228130859E+14
       coef( 1277) =    -0.15017487824307781250E+14
       coef( 1278) =    -0.12913624899930138672E+14
       coef( 1279) =    -0.24445861992389671875E+14
       coef( 1280) =    -0.29830569090730371094E+13
       coef( 1281) =    -0.48182693599492320312E+14
       coef( 1282) =     0.10601267469390696875E+15
       coef( 1283) =    -0.10620925346671050781E+14
       coef( 1284) =    -0.79740225445961093750E+14
       coef( 1285) =    -0.12615306402718800781E+14
       coef( 1286) =    -0.14032736863918555173E+05
       coef( 1287) =     0.34021204345322018489E+07
       coef( 1288) =    -0.29615267079704391956E+09
       coef( 1289) =     0.14726038115450288773E+11
       coef( 1290) =    -0.30587156573044415283E+12
       coef( 1291) =     0.22294060135328281250E+13
       coef( 1292) =    -0.78917444084017998047E+13
       coef( 1293) =     0.13031217054674226562E+14
       coef( 1294) =    -0.67492143006977714844E+13
       coef( 1295) =     0.63841419717761859298E+08
       coef( 1296) =    -0.11356877748861391068E+11
       coef( 1297) =     0.58544912569905993652E+12
       coef( 1298) =    -0.20958543098335195312E+13
       coef( 1299) =    -0.40030405063754931641E+12
       coef( 1300) =     0.24267555398413003906E+14
       coef( 1301) =    -0.63524837721016687500E+14
       coef( 1302) =     0.35124744160598484375E+14
       coef( 1303) =    -0.26316505763963830566E+12
       coef( 1304) =     0.13592156187230107422E+13
       coef( 1305) =    -0.15528015557711324219E+14
       coef( 1306) =     0.26688636460396179688E+14
       coef( 1307) =     0.30871980420892695312E+14
       coef( 1308) =    -0.24298703592650914062E+14
       coef( 1309) =     0.83811529537159453125E+13
       coef( 1310) =     0.79239881679618632812E+13
       coef( 1311) =    -0.17972093355630667969E+14
       coef( 1312) =    -0.69855675601104570312E+13
       coef( 1313) =    -0.23665931905169066406E+14
       coef( 1314) =    -0.27105124873067328125E+14
       coef( 1315) =    -0.94177505117452382812E+13
       coef( 1316) =     0.24824707863426562500E+14
       coef( 1317) =    -0.15136459443580539062E+14
       coef( 1318) =    -0.80041437143042392578E+13
       coef( 1319) =    -0.24727610242394448242E+13
       coef( 1320) =    -0.29868109487719886719E+14
       coef( 1321) =    -0.11854256590863089844E+14
       coef( 1322) =     0.38881744433080043793E+10
       coef( 1323) =    -0.25775048573559289551E+12
       coef( 1324) =    -0.44308086607925605469E+13
       coef( 1325) =     0.27190474418007710938E+14
       coef( 1326) =    -0.56798203429105796875E+14
       coef( 1327) =     0.36405825328252687500E+14
       coef( 1328) =     0.26357810121084912109E+13
       coef( 1329) =    -0.11728200637970111328E+14
       coef( 1330) =     0.39355901461715087891E+13
       coef( 1331) =    -0.20073722272032246094E+13
       coef( 1332) =     0.18186983387124355469E+14
       coef( 1333) =    -0.59280877985709218750E+14
       coef( 1334) =    -0.68736741030423808594E+13
       coef( 1335) =    -0.87198550644737597656E+12
       coef( 1336) =    -0.81987146328458618164E+12
       coef( 1337) =    -0.38775902547737054688E+14
       coef( 1338) =     0.78493325830041031250E+14
       coef( 1339) =     0.88783717877544648438E+13
       coef( 1340) =    -0.15726921050003247070E+13
       coef( 1341) =    -0.16943094659911401367E+13
       coef( 1342) =     0.50548235561019140625E+14
       coef( 1343) =     0.75420447193571972656E+13
       coef( 1344) =     0.52806902184441601562E+12
       coef( 1345) =    -0.10011374845933437500E+13
       coef( 1346) =    -0.89408979357975019531E+13
       coef( 1347) =     0.28129632467697101562E+14
       coef( 1348) =     0.18117520184677273438E+14
       coef( 1349) =    -0.23567782665193964844E+14
       coef( 1350) =    -0.13081984584424726562E+14
       coef( 1351) =    -0.59687048616701660156E+13
       coef( 1352) =    -0.67022323885610009766E+12
       coef( 1353) =    -0.63058165093874121094E+13
       coef( 1354) =    -0.12823381162617927734E+14
       coef( 1355) =    -0.67429204573335888672E+13
       coef( 1356) =    -0.74459181457505615234E+13
       coef( 1357) =    -0.53785525031779433594E+13
       coef( 1358) =    -0.22562844231914660156E+14
       coef( 1359) =    -0.13717815813907898438E+14
       coef( 1360) =    -0.72431363468108740234E+13
       coef( 1361) =    -0.64166042066130195312E+13
       coef( 1362) =     0.19450884357914166260E+12
       coef( 1363) =    -0.21779294503236103516E+13
       coef( 1364) =     0.19239956019950371094E+14
       coef( 1365) =    -0.63419959479118539062E+14
       coef( 1366) =     0.50801981050817921875E+14
       coef( 1367) =     0.38141140864579843750E+14
       coef( 1368) =    -0.46116409157753613281E+13
       coef( 1369) =    -0.17978648290780402344E+14
       coef( 1370) =     0.33356023103205356445E+13
       coef( 1371) =    -0.78261125264975781250E+14
       coef( 1372) =     0.58727226060571703125E+14
       coef( 1373) =     0.51773472923422812500E+14
       coef( 1374) =     0.20261309310549613281E+14
       coef( 1375) =     0.22774336701310410156E+13
       coef( 1376) =    -0.36675735525987031250E+14
       coef( 1377) =     0.28000127053841558594E+14
       coef( 1378) =     0.17943061312947414062E+14
       coef( 1379) =     0.62590127943421982422E+13
       coef( 1380) =     0.16525428859811781250E+14
       coef( 1381) =     0.61537514040254697266E+13
       coef( 1382) =     0.59706691442815527344E+13
       coef( 1383) =     0.16688698279933740625E+15
       coef( 1384) =     0.65431695680093484375E+14
       coef( 1385) =     0.17982755535366875000E+14
       coef( 1386) =     0.32432525094883657227E+13
       coef( 1387) =     0.58271004271005320312E+14
       coef( 1388) =     0.14340098925646804688E+14
       coef( 1389) =     0.28657154311449384766E+13
       coef( 1390) =     0.27315300643619360352E+13
       coef( 1391) =     0.11304051757053152344E+14
       coef( 1392) =     0.10571333618658228760E+13
       coef( 1393) =    -0.19950035167026742188E+14
       coef( 1394) =     0.42893123600731671875E+14
       coef( 1395) =     0.13654977147097816406E+14
       coef( 1396) =     0.23545721305520214844E+13
       coef( 1397) =     0.98855787566466464844E+13
       coef( 1398) =     0.14661162177206933594E+13
       coef( 1399) =     0.43525752024328710938E+12
       coef( 1400) =     0.53701637360473582521E+06
       coef( 1401) =     0.64294278391538411379E+08
       coef( 1402) =    -0.44307303758825645447E+10
       coef( 1403) =     0.18646712887953179932E+12
       coef( 1404) =    -0.21431627903627512207E+13
       coef( 1405) =     0.11532075835627615234E+14
       coef( 1406) =    -0.32687515166100500000E+14
       coef( 1407) =     0.54004123038855335938E+14
       coef( 1408) =    -0.42740860907219601562E+14
       coef( 1409) =     0.22852734499881172180E+10
       coef( 1410) =    -0.25273607017151742554E+12
       coef( 1411) =     0.34698026639051992188E+13
       coef( 1412) =    -0.14190183868348281250E+14
       coef( 1413) =     0.12975752255913691406E+14
       coef( 1414) =     0.51275433648910781250E+13
       coef( 1415) =     0.30453316562393710938E+14
       coef( 1416) =     0.11433233381293773438E+14
       coef( 1417) =     0.14303587942838323975E+12
       coef( 1418) =    -0.49059438274171484375E+13
       coef( 1419) =    -0.14337048254234306641E+14
       coef( 1420) =     0.20543408354311933594E+14
       coef( 1421) =    -0.16042555426586972656E+14
       coef( 1422) =    -0.26148443860693015625E+14
       coef( 1423) =    -0.17963338275828472656E+14
       coef( 1424) =     0.46942379945665421875E+14
       coef( 1425) =    -0.43275646947811976562E+14
       coef( 1426) =     0.29560037196563266602E+13
       coef( 1427) =    -0.34650529179754179688E+13
       coef( 1428) =    -0.84329007176829628906E+13
       coef( 1429) =    -0.58503406180654570312E+14
       coef( 1430) =    -0.17142048685197019531E+14
       coef( 1431) =    -0.38079270403582910156E+13
       coef( 1432) =    -0.73740300789225830078E+13
       coef( 1433) =    -0.55225432409759082794E+10
       coef( 1434) =    -0.19674996863116149902E+12
       coef( 1435) =    -0.67155257695122617188E+13
       coef( 1436) =     0.46034981361242273438E+14
       coef( 1437) =    -0.48459451213260296875E+14
       coef( 1438) =    -0.72659130522292441406E+13
       coef( 1439) =    -0.19749578497683750000E+13
       coef( 1440) =    -0.97023850401833046875E+13
       coef( 1441) =     0.89184449054707666016E+12
       coef( 1442) =    -0.22550752314422246094E+14
       coef( 1443) =    -0.11821965000273037109E+14
       coef( 1444) =    -0.24710312211518796875E+14
       coef( 1445) =    -0.98585487809070234375E+13
       coef( 1446) =    -0.49039929584986513672E+13
       coef( 1447) =    -0.10735978961002675781E+14
       coef( 1448) =     0.80283837447616015625E+13
       coef( 1449) =     0.38441151271382275391E+12
       coef( 1450) =    -0.89165562051921411133E+12
       coef( 1451) =     0.39691402925852900391E+13
       coef( 1452) =     0.51319327507855859375E+12
       coef( 1453) =     0.21627130849872523438E+14
       coef( 1454) =    -0.71998350161179199219E+13
       coef( 1455) =     0.88128742488609160156E+13
       coef( 1456) =    -0.66380380776393041992E+12
       coef( 1457) =    -0.12732764831942104492E+13
       coef( 1458) =    -0.53478140418219199219E+13
       coef( 1459) =     0.18895852011094472656E+13
       coef( 1460) =     0.12715788119999340820E+12
       coef( 1461) =     0.47577971745077807617E+12
       coef( 1462) =    -0.40433241050699306641E+13
       coef( 1463) =    -0.64473225070283544922E+12
       coef( 1464) =     0.76131321040805834961E+12
       coef( 1465) =    -0.15789731593968564453E+13
       coef( 1466) =    -0.46492953889469707031E+13
       coef( 1467) =    -0.19438409400750976562E+14
       coef( 1468) =    -0.16462809434617603516E+14
       coef( 1469) =    -0.11446134648780958984E+14
       coef( 1470) =    -0.11345801749402957031E+14
       coef( 1471) =    -0.25518280079346453125E+14
       coef( 1472) =     0.91857407720733046875E+13
       coef( 1473) =     0.61577609689761396484E+13
       coef( 1474) =     0.44057340222892578125E+11
       coef( 1475) =    -0.17996442410991669922E+13
       coef( 1476) =     0.13320271774975421875E+14
       coef( 1477) =     0.73426169296255898438E+13
       coef( 1478) =     0.22428978304674130859E+13
       coef( 1479) =     0.30872164938728466797E+13
       coef( 1480) =     0.49915035693668789062E+14
       coef( 1481) =     0.57632313538126062500E+14
       coef( 1482) =     0.17469957589151623047E+14
       coef( 1483) =     0.33319194421761679688E+13
       coef( 1484) =     0.22891462010676957031E+14
       coef( 1485) =     0.57083061623259570312E+13
       coef( 1486) =     0.66073050267828652344E+13
       coef( 1487) =     0.26363726065523312500E+14
       coef( 1488) =     0.20683556464372707031E+14
       coef( 1489) =     0.56630867690622753906E+13
       coef( 1490) =     0.63844823001807675781E+13
       coef( 1491) =     0.39856930365428663790E+08
       coef( 1492) =    -0.43232452537219276428E+10
       coef( 1493) =     0.85388743931469558716E+11
       coef( 1494) =    -0.80100110444092187500E+12
       coef( 1495) =     0.46180364882827744141E+13
       coef( 1496) =    -0.21781207343360250000E+14
       coef( 1497) =     0.53603455194872875000E+14
       coef( 1498) =     0.21513939308665609375E+14
       coef( 1499) =    -0.12089541614711039062E+15
       coef( 1500) =     0.85434644662753250122E+11
       coef( 1501) =    -0.22335991184961054688E+13
       coef( 1502) =     0.14201507202991531250E+14
       coef( 1503) =    -0.20694071999685847656E+14
       coef( 1504) =    -0.42779617089832085938E+14
       coef( 1505) =    -0.54932508943116101562E+14
       coef( 1506) =    -0.19797884254793113281E+14
       coef( 1507) =    -0.15838686302654175781E+14
       coef( 1508) =     0.52755442660255976562E+13
       coef( 1509) =    -0.38836820770802695312E+14
       coef( 1510) =     0.13759218884355685938E+15
       coef( 1511) =     0.11464403872989357812E+15
       coef( 1512) =     0.34994984297825605469E+14
       coef( 1513) =     0.46771736196772558594E+13
       coef( 1514) =    -0.58793186500246820312E+14
       coef( 1515) =    -0.17240456810090097656E+14
       coef( 1516) =     0.16858144127423498047E+14
       coef( 1517) =     0.82913119861334462891E+13
       coef( 1518) =    -0.24048833184696726562E+14
       coef( 1519) =    -0.38209134856368886719E+13
       coef( 1520) =    -0.24174274976271377563E+12
       coef( 1521) =     0.85854288878134179688E+13
       coef( 1522) =    -0.66448924480989375000E+14
       coef( 1523) =     0.10673149430532646875E+15
       coef( 1524) =    -0.51860708659688875000E+14
       coef( 1525) =    -0.41939554920848359375E+14
       coef( 1526) =    -0.14711950587084136719E+14
       coef( 1527) =    -0.14152496419037679688E+14
       coef( 1528) =     0.10371281472043923438E+15
       coef( 1529) =     0.70786512507022156250E+14
       coef( 1530) =     0.74858534804145195312E+13
       coef( 1531) =    -0.33841673928317182617E+13
       coef( 1532) =     0.45914851095492632812E+14
       coef( 1533) =     0.20255952820454117188E+14
       coef( 1534) =     0.45431183980998164062E+13
       coef( 1535) =     0.56573821262348808594E+13

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

      norder=8
      norder2=8
      maxorder=13

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

