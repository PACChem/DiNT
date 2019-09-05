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

       ncoef=1084
       coef(    1) =    -0.24434545821422699596E+01
       coef(    2) =     0.57120541750657473301E+02
       coef(    3) =    -0.14724778803770761442E+05
       coef(    4) =    -0.20976019935127627105E+06
       coef(    5) =    -0.84660731800853520632E+08
       coef(    6) =     0.61341187732794916630E+09
       coef(    7) =    -0.33892245508591442108E+10
       coef(    8) =     0.11873244231931465149E+11
       coef(    9) =    -0.29342982804734332603E+05
       coef(   10) =     0.54027775880456408486E+07
       coef(   11) =     0.68238618226302933693E+09
       coef(   12) =    -0.12123594088461918831E+10
       coef(   13) =     0.11246184006822006226E+11
       coef(   14) =    -0.43369309548264068604E+11
       coef(   15) =    -0.17850827260525103760E+12
       coef(   16) =    -0.57232796599280166626E+09
       coef(   17) =     0.45036184803717365265E+10
       coef(   18) =    -0.66159799134840728760E+11
       coef(   19) =    -0.59821827291764526367E+11
       coef(   20) =     0.11906203442738867188E+13
       coef(   21) =    -0.11863979444896123047E+13
       coef(   22) =    -0.55177952991845970154E+11
       coef(   23) =     0.16809146098836228027E+13
       coef(   24) =    -0.56091334034951142578E+13
       coef(   25) =     0.11926381795024943359E+14
       coef(   26) =    -0.13854564975192031250E+14
       coef(   27) =    -0.74184502301071943359E+13
       coef(   28) =     0.38469835576492179688E+14
       coef(   29) =    -0.84263299487482500000E+14
       coef(   30) =     0.99707099920992218750E+14
       coef(   31) =    -0.40931448158180937500E+13
       coef(   32) =    -0.22702812034459859375E+14
       coef(   33) =     0.66919900130373390625E+14
       coef(   34) =    -0.15500488430089527344E+14
       coef(   35) =    -0.15597096759631303139E+07
       coef(   36) =    -0.88799166598568212986E+09
       coef(   37) =    -0.21145142164970626831E+11
       coef(   38) =     0.36967720994899551392E+11
       coef(   39) =     0.29161057735781103516E+12
       coef(   40) =    -0.77158444840203515625E+12
       coef(   41) =     0.32041765794239179688E+13
       coef(   42) =     0.52164339378167671204E+11
       coef(   43) =    -0.12591883901698854065E+12
       coef(   44) =     0.19020193635681589355E+13
       coef(   45) =    -0.10319556403088867188E+14
       coef(   46) =     0.70972766768957607422E+13
       coef(   47) =    -0.59369238312538291016E+13
       coef(   48) =    -0.27542296780399130859E+13
       coef(   49) =     0.12794662095494125000E+14
       coef(   50) =     0.16156077878591599609E+14
       coef(   51) =    -0.49794029449666750000E+14
       coef(   52) =    -0.65837586734844531250E+12
       coef(   53) =    -0.12328723304658347656E+14
       coef(   54) =    -0.52513234992619843750E+13
       coef(   55) =    -0.69983759976984453125E+13
       coef(   56) =     0.19859248420000843750E+14
       coef(   57) =    -0.45003054840383492188E+14
       coef(   58) =    -0.27859449138908707031E+14
       coef(   59) =    -0.51053832181463360596E+12
       coef(   60) =     0.65899754247642861328E+13
       coef(   61) =    -0.21332997901679136719E+14
       coef(   62) =     0.44664077705040820312E+14
       coef(   63) =    -0.84490900110395830078E+13
       coef(   64) =     0.11926852381939447266E+14
       coef(   65) =    -0.14058653745859439453E+14
       coef(   66) =     0.20809617638903078125E+14
       coef(   67) =    -0.14799329240713117188E+14
       coef(   68) =    -0.30599796515312105469E+14
       coef(   69) =    -0.99879312901233417969E+13
       coef(   70) =     0.53943684104254773438E+14
       coef(   71) =     0.21426021475291160156E+14
       coef(   72) =     0.23400675328445761719E+13
       coef(   73) =     0.63277592822931044922E+13
       coef(   74) =     0.17938340440290882812E+14
       coef(   75) =    -0.60592304832063031250E+14
       coef(   76) =    -0.41154993909107515625E+14
       coef(   77) =    -0.22048941989304812500E+14
       coef(   78) =    -0.28717833074216234375E+14
       coef(   79) =    -0.87771401253409062500E+13
       coef(   80) =    -0.12422365754945703125E+14
       coef(   81) =     0.62873442512897109985E+09
       coef(   82) =    -0.26903314884366886139E+11
       coef(   83) =     0.46188066260250579834E+12
       coef(   84) =    -0.33464265631012910156E+13
       coef(   85) =     0.10152315576695726562E+14
       coef(   86) =    -0.16588792516390576172E+14
       coef(   87) =     0.24637120292115966797E+12
       coef(   88) =     0.13097617590588394165E+12
       coef(   89) =    -0.29696576930349135742E+13
       coef(   90) =     0.20051543184804863281E+14
       coef(   91) =    -0.27204952898978207031E+14
       coef(   92) =     0.30947941183982242188E+14
       coef(   93) =     0.62163569879231734375E+14
       coef(   94) =     0.13410638029702365234E+14
       coef(   95) =    -0.88034127680337968750E+14
       coef(   96) =     0.59275623520853054688E+14
       coef(   97) =     0.71661678953713750000E+14
       coef(   98) =     0.48192560608167187500E+14
       coef(   99) =     0.39241497329891796875E+14
       coef(  100) =     0.55716614479711781250E+14
       coef(  101) =     0.35846480061686046875E+14
       coef(  102) =     0.17212340865310470703E+14
       coef(  103) =    -0.62171216604072998047E+12
       coef(  104) =    -0.42903035420491035156E+13
       coef(  105) =    -0.78146559837801474609E+13
       coef(  106) =    -0.90156538772585437500E+14
       coef(  107) =    -0.41660893416127632812E+14
       coef(  108) =    -0.15846260930048789062E+13
       coef(  109) =     0.22052790348128152344E+14
       coef(  110) =    -0.17536628856904369141E+14
       coef(  111) =    -0.26174620000016156250E+14
       coef(  112) =    -0.10661754312512757812E+14
       coef(  113) =     0.37529824546136621094E+12
       coef(  114) =    -0.89825491051228906250E+12
       coef(  115) =     0.95107018743568218750E+14
       coef(  116) =     0.16587220950466041016E+14
       coef(  117) =    -0.56363039029939062500E+13
       coef(  118) =     0.45181757809425966797E+13
       coef(  119) =     0.68451041685026708984E+13
       coef(  120) =    -0.43123123019995593750E+14
       coef(  121) =     0.14324065451813375000E+15
       coef(  122) =     0.43200833111877890625E+14
       coef(  123) =     0.81525265322242460938E+13
       coef(  124) =    -0.62912597606135867188E+14
       coef(  125) =     0.72390300448716796875E+13
       coef(  126) =     0.23786464601925893555E+13
       coef(  127) =     0.13931303177022470703E+13
       coef(  128) =    -0.31963544044290226562E+14
       coef(  129) =    -0.37288967060587070312E+13
       coef(  130) =    -0.16351187402856238281E+14
       coef(  131) =    -0.23063165356605651368E+03
       coef(  132) =    -0.15015053042460056531E+05
       coef(  133) =     0.31153039889453519136E+07
       coef(  134) =    -0.36321756997665174305E+08
       coef(  135) =     0.14017564664959228039E+10
       coef(  136) =     0.45305431137886209488E+10
       coef(  137) =    -0.68346203455056411743E+11
       coef(  138) =     0.75169701718110305786E+11
       coef(  139) =     0.60713986148382478859E+05
       coef(  140) =    -0.47660008275270235538E+09
       coef(  141) =    -0.43441733681350374222E+10
       coef(  142) =    -0.10335902471887097168E+12
       coef(  143) =     0.55777707586499291992E+12
       coef(  144) =     0.84197374017827270508E+11
       coef(  145) =    -0.95599198108888232422E+12
       coef(  146) =     0.11836152982653999329E+11
       coef(  147) =    -0.10191376007626126099E+12
       coef(  148) =     0.18070570260072133789E+13
       coef(  149) =    -0.10135764767065703125E+14
       coef(  150) =     0.23631457155352511719E+14
       coef(  151) =    -0.13143950241193791016E+14
       coef(  152) =    -0.16822189400608535767E+12
       coef(  153) =    -0.16723565627438051758E+13
       coef(  154) =     0.55071459907391611328E+13
       coef(  155) =    -0.37919279989613984375E+14
       coef(  156) =     0.10164339670346222656E+14
       coef(  157) =     0.59153752950215673828E+13
       coef(  158) =     0.67925428628596103516E+13
       coef(  159) =    -0.25855907212699085938E+14
       coef(  160) =    -0.27660940167436113281E+13
       coef(  161) =     0.19752188717250750000E+14
       coef(  162) =    -0.15055322030809785156E+13
       coef(  163) =     0.39745271395296072960E+09
       coef(  164) =     0.72172927731051187515E+10
       coef(  165) =     0.40519508390780346680E+12
       coef(  166) =    -0.61549080103531555176E+12
       coef(  167) =    -0.61592723164436171875E+13
       coef(  168) =     0.13135550543721242188E+14
       coef(  169) =    -0.63834158637315244141E+13
       coef(  170) =    -0.58143438973041772461E+12
       coef(  171) =    -0.13818876961926853027E+13
       coef(  172) =     0.85508201710823378906E+13
       coef(  173) =     0.71395030408702285156E+13
       coef(  174) =    -0.30871278479246503906E+14
       coef(  175) =    -0.40218689607460757812E+14
       coef(  176) =     0.18987441205827703125E+14
       coef(  177) =    -0.89732759071208531250E+14
       coef(  178) =     0.16354004433787540625E+15
       coef(  179) =     0.11815853198175798438E+15
       coef(  180) =     0.49825710317830343750E+14
       coef(  181) =    -0.67766378374639406250E+14
       coef(  182) =     0.50687266017084468750E+14
       coef(  183) =     0.42000363881964734375E+14
       coef(  184) =     0.26606748092843542969E+14
       coef(  185) =     0.52706201334068750000E+13
       coef(  186) =    -0.48620162201415742188E+14
       coef(  187) =     0.73751045523685828125E+14
       coef(  188) =    -0.80977324316876578125E+14
       coef(  189) =    -0.83309688176622703125E+14
       coef(  190) =    -0.45564013500233109375E+14
       coef(  191) =     0.91081928934928312500E+14
       coef(  192) =    -0.43471380564332875000E+14
       coef(  193) =    -0.28746687789964234375E+14
       coef(  194) =    -0.12866465388785859375E+14
       coef(  195) =    -0.41262852187038421875E+14
       coef(  196) =    -0.93457201592028574219E+13
       coef(  197) =    -0.31363257075244031250E+14
       coef(  198) =    -0.58391217494748390625E+14
       coef(  199) =    -0.24995392271181621094E+14
       coef(  200) =    -0.28599813722157828125E+14
       coef(  201) =    -0.12043692973183826447E+11
       coef(  202) =     0.17696455080589862061E+12
       coef(  203) =    -0.36530148181980571289E+13
       coef(  204) =     0.11743454653746855469E+14
       coef(  205) =    -0.37280162897736411133E+13
       coef(  206) =    -0.16566185805425085938E+14
       coef(  207) =     0.26464982373558710938E+14
       coef(  208) =     0.26407401650700825195E+13
       coef(  209) =     0.89034128105414316406E+13
       coef(  210) =    -0.26251828605217460938E+14
       coef(  211) =    -0.32961388410785839844E+13
       coef(  212) =    -0.20612955919099878906E+14
       coef(  213) =    -0.11866157110247921875E+14
       coef(  214) =    -0.19138636067636761719E+14
       coef(  215) =     0.10208318481453886719E+14
       coef(  216) =     0.45898752551279164062E+14
       coef(  217) =     0.23656643369930652344E+14
       coef(  218) =     0.14088360452325912109E+14
       coef(  219) =     0.18877006081382781250E+14
       coef(  220) =    -0.25251549897590808594E+14
       coef(  221) =     0.61948404701083296875E+14
       coef(  222) =     0.14730525744722171875E+14
       coef(  223) =    -0.91069184667713242188E+13
       coef(  224) =    -0.14076535006663531250E+14
       coef(  225) =    -0.55892244025243140625E+14
       coef(  226) =    -0.27361101513144140625E+14
       coef(  227) =    -0.69117557911294052734E+13
       coef(  228) =    -0.11178323140131324219E+14
       coef(  229) =    -0.56214003218872320312E+14
       coef(  230) =    -0.20884463870885332031E+14
       coef(  231) =     0.50832998618173046875E+14
       coef(  232) =    -0.90768071372000500000E+14
       coef(  233) =    -0.25823823600709371094E+14
       coef(  234) =    -0.78204421247116005859E+13
       coef(  235) =    -0.75007302302132640625E+14
       coef(  236) =    -0.21442309966779984375E+14
       coef(  237) =    -0.33537929377331285156E+14
       coef(  238) =    -0.31776585830421699939E+05
       coef(  239) =     0.42054419139068117365E+07
       coef(  240) =    -0.47125379714033091068E+09
       coef(  241) =     0.15159512399813467026E+11
       coef(  242) =    -0.22018088182351855469E+12
       coef(  243) =     0.12300503305561784668E+13
       coef(  244) =    -0.26019126119262031250E+13
       coef(  245) =     0.21952798202525749512E+13
       coef(  246) =     0.42683136463122856617E+09
       coef(  247) =    -0.20848756096053710938E+11
       coef(  248) =     0.58133792324219580078E+12
       coef(  249) =    -0.36339653553253500977E+13
       coef(  250) =     0.12521281986024048828E+14
       coef(  251) =    -0.36286508029539437500E+14
       coef(  252) =     0.37041946098173734375E+14
       coef(  253) =    -0.64978690023755073547E+11
       coef(  254) =    -0.13955518988561666870E+12
       coef(  255) =    -0.13735842340352091797E+14
       coef(  256) =     0.47456062396277812500E+14
       coef(  257) =    -0.38794191384630843750E+14
       coef(  258) =    -0.12734524455038480469E+14
       coef(  259) =     0.79891819176434052734E+13
       coef(  260) =    -0.16461373909704248047E+12
       coef(  261) =    -0.65299447516924218750E+13
       coef(  262) =     0.12733235631947609375E+14
       coef(  263) =     0.24880030018742937500E+14
       coef(  264) =    -0.19597251500619777344E+14
       coef(  265) =    -0.14995226584386566406E+14
       coef(  266) =     0.27748438553006250000E+12
       coef(  267) =    -0.84481714942221386719E+13
       coef(  268) =    -0.12770863128256678581E+10
       coef(  269) =    -0.22218202349797933960E+12
       coef(  270) =    -0.88932429363774243164E+12
       coef(  271) =     0.16623807289266687500E+14
       coef(  272) =    -0.25552493258777945312E+14
       coef(  273) =     0.26158798392986265625E+14
       coef(  274) =     0.13514002872837666016E+13
       coef(  275) =     0.11972479098038500977E+13
       coef(  276) =     0.49807323582848261719E+13
       coef(  277) =    -0.35183024356298851562E+14
       coef(  278) =    -0.14950574577385800781E+14
       coef(  279) =    -0.14291378789391578125E+14
       coef(  280) =    -0.13140518510865832031E+14
       coef(  281) =    -0.44667815641379421875E+14
       coef(  282) =     0.43504489339298476562E+14
       coef(  283) =     0.53463403180332093750E+14
       coef(  284) =     0.26848579720331691406E+14
       coef(  285) =     0.37943397994276945312E+14
       coef(  286) =     0.24515077273789140625E+14
       coef(  287) =    -0.92569070606546152344E+13
       coef(  288) =     0.57492258096934445312E+14
       coef(  289) =     0.60514144968449429688E+14
       coef(  290) =     0.11567428033076904297E+13
       coef(  291) =    -0.12637962338052681641E+14
       coef(  292) =    -0.17952828840393187500E+14
       coef(  293) =     0.16672961486546301270E+13
       coef(  294) =     0.13558918664846865234E+13
       coef(  295) =     0.80564472554769055176E+12
       coef(  296) =    -0.38655531176047976562E+14
       coef(  297) =    -0.13323032885712541016E+14
       coef(  298) =     0.23699314078090615845E+12
       coef(  299) =    -0.31322236979435439453E+13
       coef(  300) =     0.93652329907958144531E+13
       coef(  301) =    -0.32283596446632406250E+14
       coef(  302) =     0.85131780734266904297E+13
       coef(  303) =    -0.32664055810769257812E+13
       coef(  304) =    -0.54794343996263046875E+13
       coef(  305) =     0.18225537915535656250E+14
       coef(  306) =    -0.97780815019481171875E+14
       coef(  307) =     0.81821500211233546875E+14
       coef(  308) =     0.57093551929166101562E+14
       coef(  309) =     0.15664680753598464844E+14
       coef(  310) =     0.27079212436751820312E+14
       coef(  311) =     0.54071766108184226562E+14
       coef(  312) =     0.29910668816570406250E+14
       coef(  313) =     0.24459364955181812500E+14
       coef(  314) =    -0.41199569297910878906E+13
       coef(  315) =     0.60919239875843218750E+14
       coef(  316) =     0.58238261524922484375E+14
       coef(  317) =     0.23548926914673007812E+14
       coef(  318) =     0.24913838363338156250E+14
       coef(  319) =     0.13265104457752810547E+14
       coef(  320) =     0.18999406390518503418E+13
       coef(  321) =    -0.19533633441447753906E+14
       coef(  322) =     0.40539221013886494141E+13
       coef(  323) =     0.11101423984637375000E+14
       coef(  324) =    -0.26305864594910063477E+13
       coef(  325) =     0.10167270235168052604E+07
       coef(  326) =    -0.25544717695321702957E+09
       coef(  327) =     0.16716361364124353409E+11
       coef(  328) =    -0.44152907526370782471E+12
       coef(  329) =     0.51334294593118144531E+13
       coef(  330) =    -0.25745162578227230469E+14
       coef(  331) =     0.55161936865306140625E+14
       coef(  332) =    -0.47663592627919171875E+14
       coef(  333) =    -0.57075456206532783508E+10
       coef(  334) =     0.15546125310142114258E+12
       coef(  335) =    -0.56003395564837470703E+13
       coef(  336) =     0.28655767470416078125E+14
       coef(  337) =    -0.49458512432237960938E+14
       coef(  338) =     0.67132309203607187500E+14
       coef(  339) =     0.45994590347535750000E+14
       coef(  340) =     0.11243943952713459473E+13
       coef(  341) =     0.55969536141671289062E+13
       coef(  342) =     0.44723819404269234375E+14
       coef(  343) =    -0.37169847852443015625E+14
       coef(  344) =    -0.96099161445175875000E+14
       coef(  345) =    -0.59149273351982140625E+14
       coef(  346) =    -0.78640407151534312500E+14
       coef(  347) =     0.43348378956572375000E+14
       coef(  348) =    -0.21341956292767234375E+14
       coef(  349) =    -0.31486553098092453125E+14
       coef(  350) =     0.46000543724101039062E+14
       coef(  351) =     0.28434463636770292969E+13
       coef(  352) =     0.13801789528735248566E+11
       coef(  353) =     0.48676063750975917969E+13
       coef(  354) =    -0.25222747725729359375E+14
       coef(  355) =    -0.53184691283452648438E+14
       coef(  356) =     0.40577826458003312500E+14
       coef(  357) =     0.41025285583214109375E+14
       coef(  358) =     0.77519846661446367188E+13
       coef(  359) =     0.16270210611514472656E+13
       coef(  360) =     0.48400249407763195312E+14
       coef(  361) =    -0.36766341482269781250E+14
       coef(  362) =    -0.35950116513357265625E+14
       coef(  363) =    -0.23498543213161238281E+14
       coef(  364) =     0.12234845561019148438E+15
       coef(  365) =     0.56224875800765804688E+14
       coef(  366) =     0.13260743911600052734E+14
       coef(  367) =     0.29370893068059246094E+14
       coef(  368) =    -0.20862596158036082031E+14
       coef(  369) =    -0.17511484717590743750E+15
       coef(  370) =    -0.61755495593437085938E+14
       coef(  371) =    -0.24165345546107187500E+14
       coef(  372) =    -0.84254904458298718750E+14
       coef(  373) =    -0.18720407772047238281E+14
       coef(  374) =    -0.38455833336167093750E+14
       coef(  375) =    -0.18880289356236672363E+13
       coef(  376) =    -0.68877035550895087891E+13
       coef(  377) =     0.12290383258803256250E+15
       coef(  378) =    -0.26508509806711289062E+14
       coef(  379) =    -0.28321718418778695312E+14
       coef(  380) =    -0.17812326362384312500E+14
       coef(  381) =    -0.11235810677918835938E+14
       coef(  382) =    -0.10128524823558534375E+15
       coef(  383) =    -0.23525890125013968750E+14
       coef(  384) =    -0.58708388981115781250E+13
       coef(  385) =    -0.11179427926027205078E+14
       coef(  386) =     0.82724269408661474609E+13
       coef(  387) =     0.84645736527465341797E+13
       coef(  388) =    -0.82617687145764960938E+13
       coef(  389) =     0.21676122936379892578E+13
       coef(  390) =     0.72317556813783496094E+12
       coef(  391) =     0.89133906083915800781E+13
       coef(  392) =     0.23753070570908032227E+13
       coef(  393) =     0.41577237363544330001E+07
       coef(  394) =     0.16729975245070166588E+10
       coef(  395) =    -0.90310033176977523804E+11
       coef(  396) =     0.21668559609573286133E+13
       coef(  397) =    -0.20245280456670476562E+14
       coef(  398) =     0.73787350707227500000E+14
       coef(  399) =    -0.98512681221430125000E+14
       coef(  400) =     0.29403483399099226562E+14
       coef(  401) =     0.43834192850252746582E+11
       coef(  402) =    -0.33350425489295428467E+12
       coef(  403) =     0.15054589148690160156E+14
       coef(  404) =    -0.44720882137751671875E+14
       coef(  405) =    -0.67761007136343515625E+13
       coef(  406) =    -0.27041193442275722656E+14
       coef(  407) =    -0.21068274106950515625E+14
       coef(  408) =    -0.97303546652139042969E+13
       coef(  409) =     0.19483478916159269531E+14
       coef(  410) =    -0.39425182931673648438E+14
       coef(  411) =    -0.47038245369189656250E+14
       coef(  412) =    -0.62166794308530421875E+14
       coef(  413) =     0.87538402197222968750E+14
       coef(  414) =     0.71772516707840984375E+14
       coef(  415) =     0.65746359119424941406E+13
       coef(  416) =     0.50129410320819906250E+14
       coef(  417) =    -0.95259758970360119629E+12
       coef(  418) =    -0.84290907802180859375E+13
       coef(  419) =     0.44232153213026257812E+14
       coef(  420) =     0.63453481099571914062E+14
       coef(  421) =     0.83916731965539312500E+14
       coef(  422) =     0.36059357995184062500E+14
       coef(  423) =    -0.16870349387461826172E+14
       coef(  424) =     0.23662563228137625000E+14
       coef(  425) =     0.68714090763314511719E+13
       coef(  426) =    -0.81543763120177441406E+12
       coef(  427) =     0.76977606332322484375E+14
       coef(  428) =     0.32731612877138171875E+14
       coef(  429) =     0.10962132108294615625E+15
       coef(  430) =    -0.47355065156919468750E+14
       coef(  431) =    -0.20458217270233476562E+14
       coef(  432) =    -0.25891453439056890625E+14
       coef(  433) =     0.10115535433196876953E+14
       coef(  434) =    -0.63677330533787050781E+13
       coef(  435) =    -0.11916688528118082812E+15
       coef(  436) =    -0.51849175626595437500E+14
       coef(  437) =    -0.15435618057787908203E+14
       coef(  438) =    -0.10833540105615078125E+13
       coef(  439) =    -0.60513641163245031250E+14
       coef(  440) =    -0.20497201682619777344E+14
       coef(  441) =    -0.11162516040668095703E+14
       coef(  442) =     0.30463541507598457031E+14
       coef(  443) =    -0.51739454522211738281E+13
       coef(  444) =     0.15389137950586234375E+14
       coef(  445) =    -0.16670925009462952614E+09
       coef(  446) =    -0.34116207151314654350E+10
       coef(  447) =     0.22112186751683883667E+12
       coef(  448) =    -0.53094260178457988281E+13
       coef(  449) =     0.36812963554369015625E+14
       coef(  450) =    -0.65984890000126484375E+14
       coef(  451) =    -0.96546421447342753906E+13
       coef(  452) =     0.86707312814213562500E+14
       coef(  453) =    -0.17001579084072882080E+12
       coef(  454) =     0.20179433331733349609E+13
       coef(  455) =    -0.23448752832880714844E+14
       coef(  456) =     0.15017100804040390625E+13
       coef(  457) =     0.90407979862627890625E+14
       coef(  458) =     0.36009342351173171875E+14
       coef(  459) =     0.21360433870192804688E+14
       coef(  460) =    -0.67279822991402078125E+14
       coef(  461) =    -0.24413789019098335938E+14
       coef(  462) =     0.46283499982424365234E+12
       coef(  463) =    -0.44907134569010757812E+14
       coef(  464) =     0.48080466389551093750E+13
       coef(  465) =     0.31266838262006376953E+13
       coef(  466) =     0.15593716907788089844E+14
       coef(  467) =    -0.75060783178203937500E+14
       coef(  468) =     0.28893867362808847656E+14
       coef(  469) =     0.56516551149752062500E+14
       coef(  470) =    -0.81810650478385109375E+14
       coef(  471) =     0.60504781843753183594E+13
       coef(  472) =     0.14467760606547267578E+14
       coef(  473) =     0.31668889733757343750E+14
       coef(  474) =     0.12410949317675903125E+15
       coef(  475) =     0.14764118277117628906E+14
       coef(  476) =    -0.23069377821698144531E+14
       coef(  477) =     0.31404130834846695312E+14
       coef(  478) =    -0.53987414403405789062E+14
       coef(  479) =    -0.12689127854740152344E+14
       coef(  480) =     0.24284969979956082031E+14
       coef(  481) =    -0.13723196309611433594E+14
       coef(  482) =     0.30122408856930398438E+14
       coef(  483) =     0.10513872006589440107E+10
       coef(  484) =    -0.16232592984107505798E+11
       coef(  485) =    -0.12670554495338572693E+12
       coef(  486) =     0.57039514108195507812E+13
       coef(  487) =    -0.26070843524019734375E+14
       coef(  488) =    -0.37041640971158648438E+14
       coef(  489) =     0.76260973882865843750E+14
       coef(  490) =     0.68578416148569775391E+12
       coef(  491) =    -0.13328096444319242188E+14
       coef(  492) =     0.58555545510862031250E+14
       coef(  493) =    -0.48085659377476562500E+11
       coef(  494) =     0.47188682412694828125E+14
       coef(  495) =     0.36582679061425703125E+13
       coef(  496) =     0.71568397635658789062E+13
       coef(  497) =     0.43880769251569628906E+13
       coef(  498) =    -0.32754337929884273438E+14
       coef(  499) =    -0.63811177577317197266E+13
       coef(  500) =     0.34302995883598359375E+14
       coef(  501) =     0.37730293775423266602E+12
       coef(  502) =     0.19886211548197347656E+14
       coef(  503) =    -0.54112505013097171875E+14
       coef(  504) =     0.12855512814515634766E+14
       coef(  505) =     0.68232342250086406250E+14
       coef(  506) =     0.18575846474884492188E+14
       coef(  507) =     0.18656518544503234375E+14
       coef(  508) =    -0.28827210491918935547E+13
       coef(  509) =     0.15847223051711710938E+14
       coef(  510) =    -0.14911759105976617336E+10
       coef(  511) =     0.13026026247946798325E+11
       coef(  512) =     0.22871335203487298584E+12
       coef(  513) =    -0.36475639877754956055E+13
       coef(  514) =     0.72847078668231933594E+13
       coef(  515) =     0.58129216624875867188E+14
       coef(  516) =    -0.26315717007843688965E+12
       coef(  517) =     0.70288676934032656250E+13
       coef(  518) =    -0.39037907448151203125E+14
       coef(  519) =    -0.34235742468032550781E+14
       coef(  520) =     0.27855941998808359375E+13
       coef(  521) =     0.24424363824546144531E+14
       coef(  522) =     0.24300888192228120117E+13
       coef(  523) =    -0.17540489417634867188E+14
       coef(  524) =     0.95293121196814687500E+13
       coef(  525) =    -0.32280627663110851562E+14
       coef(  526) =    -0.23840520753541015625E+13
       coef(  527) =    -0.37650827156371347656E+13
       coef(  528) =     0.26520369597119292848E+03
       coef(  529) =    -0.17365571825484548754E+05
       coef(  530) =     0.31167594100612355396E+07
       coef(  531) =    -0.36321087897476971149E+08
       coef(  532) =     0.14017566728073906898E+10
       coef(  533) =     0.45305418378086824417E+10
       coef(  534) =    -0.68346204718690353394E+11
       coef(  535) =     0.75169701432383224487E+11
       coef(  536) =     0.60527623346493171994E+05
       coef(  537) =    -0.47659984065916788578E+09
       coef(  538) =    -0.43441732935493831635E+10
       coef(  539) =    -0.10335902462553985596E+12
       coef(  540) =     0.55777707675352770996E+12
       coef(  541) =     0.84197375354710708618E+11
       coef(  542) =    -0.95599198106683007812E+12
       coef(  543) =     0.11836153557389921188E+11
       coef(  544) =    -0.10191376073401982117E+12
       coef(  545) =     0.18070570278994907227E+13
       coef(  546) =    -0.10135764767963539062E+14
       coef(  547) =     0.23631457155493722656E+14
       coef(  548) =    -0.13143950241070632812E+14
       coef(  549) =    -0.16822189433831036377E+12
       coef(  550) =    -0.16723565619161882324E+13
       coef(  551) =     0.55071459905543125000E+13
       coef(  552) =    -0.37919279991343796875E+14
       coef(  553) =     0.10164339670796224609E+14
       coef(  554) =     0.59153752961546367188E+13
       coef(  555) =     0.67925428617828886719E+13
       coef(  556) =    -0.25855907211536359375E+14
       coef(  557) =    -0.27660940164209365234E+13
       coef(  558) =     0.19752188717969261719E+14
       coef(  559) =    -0.15055322024792329102E+13
       coef(  560) =     0.39745230705706900358E+09
       coef(  561) =     0.72172932886561508179E+10
       coef(  562) =     0.40519508468009210205E+12
       coef(  563) =    -0.61549080144824755859E+12
       coef(  564) =    -0.61592723159667812500E+13
       coef(  565) =     0.13135550543967839844E+14
       coef(  566) =    -0.63834158633240322266E+13
       coef(  567) =    -0.58143438988368127441E+12
       coef(  568) =    -0.13818876958253449707E+13
       coef(  569) =     0.85508201710686416016E+13
       coef(  570) =     0.71395030410711572266E+13
       coef(  571) =    -0.30871278479357234375E+14
       coef(  572) =    -0.40218689607235515625E+14
       coef(  573) =     0.18987441205573164062E+14
       coef(  574) =    -0.89732759070983203125E+14
       coef(  575) =     0.16354004433820187500E+15
       coef(  576) =     0.11815853198190775000E+15
       coef(  577) =     0.49825710317769812500E+14
       coef(  578) =    -0.67766378375051578125E+14
       coef(  579) =     0.50687266016942390625E+14
       coef(  580) =     0.42000363882023609375E+14
       coef(  581) =     0.26606748092702460938E+14
       coef(  582) =     0.52706201334899716797E+13
       coef(  583) =    -0.48620162201485773438E+14
       coef(  584) =     0.73751045523628781250E+14
       coef(  585) =    -0.80977324316926906250E+14
       coef(  586) =    -0.83309688176547968750E+14
       coef(  587) =    -0.45564013500182757812E+14
       coef(  588) =     0.91081928934753500000E+14
       coef(  589) =    -0.43471380564379812500E+14
       coef(  590) =    -0.28746687790041609375E+14
       coef(  591) =    -0.12866465388592562500E+14
       coef(  592) =    -0.41262852187068148438E+14
       coef(  593) =    -0.93457201591295566406E+13
       coef(  594) =    -0.31363257075275109375E+14
       coef(  595) =    -0.58391217494781273438E+14
       coef(  596) =    -0.24995392271277980469E+14
       coef(  597) =    -0.28599813722140898438E+14
       coef(  598) =    -0.12043693016847681046E+11
       coef(  599) =     0.17696455077180792236E+12
       coef(  600) =    -0.36530148183072822266E+13
       coef(  601) =     0.11743454653805699219E+14
       coef(  602) =    -0.37280162897293505859E+13
       coef(  603) =    -0.16566185805395359375E+14
       coef(  604) =     0.26464982373545007812E+14
       coef(  605) =     0.26407401650512255859E+13
       coef(  606) =     0.89034128105353457031E+13
       coef(  607) =    -0.26251828605214734375E+14
       coef(  608) =    -0.32961388411164667969E+13
       coef(  609) =    -0.20612955919065785156E+14
       coef(  610) =    -0.11866157110234109375E+14
       coef(  611) =    -0.19138636067586253906E+14
       coef(  612) =     0.10208318481433591797E+14
       coef(  613) =     0.45898752551259835938E+14
       coef(  614) =     0.23656643369939503906E+14
       coef(  615) =     0.14088360452358419922E+14
       coef(  616) =     0.18877006081375312500E+14
       coef(  617) =    -0.25251549897553351562E+14
       coef(  618) =     0.61948404701082203125E+14
       coef(  619) =     0.14730525744681847656E+14
       coef(  620) =    -0.91069184667645332031E+13
       coef(  621) =    -0.14076535006691392578E+14
       coef(  622) =    -0.55892244025252312500E+14
       coef(  623) =    -0.27361101513119664062E+14
       coef(  624) =    -0.69117557911334238281E+13
       coef(  625) =    -0.11178323140132574219E+14
       coef(  626) =    -0.56214003218847601562E+14
       coef(  627) =    -0.20884463870898128906E+14
       coef(  628) =     0.50832998618191375000E+14
       coef(  629) =    -0.90768071371984437500E+14
       coef(  630) =    -0.25823823600692207031E+14
       coef(  631) =    -0.78204421247024179688E+13
       coef(  632) =    -0.75007302302169828125E+14
       coef(  633) =    -0.21442309966773085938E+14
       coef(  634) =    -0.33537929377332437500E+14
       coef(  635) =     0.71267080682723011705E+05
       coef(  636) =    -0.10353901621053382754E+08
       coef(  637) =     0.99811181343006336689E+09
       coef(  638) =    -0.32748082408294677734E+11
       coef(  639) =     0.45317689435100244141E+12
       coef(  640) =    -0.28422581329349208984E+13
       coef(  641) =     0.72595628050049785156E+13
       coef(  642) =    -0.54069089150845957031E+13
       coef(  643) =    -0.99614606517123806477E+09
       coef(  644) =     0.63972877295210876465E+11
       coef(  645) =    -0.10611568884239313965E+13
       coef(  646) =     0.59733540961928984375E+13
       coef(  647) =    -0.87705213002713720703E+13
       coef(  648) =    -0.64015967673410429688E+13
       coef(  649) =    -0.21954529150729046875E+14
       coef(  650) =    -0.12966445684162513733E+12
       coef(  651) =     0.38542038528615522461E+13
       coef(  652) =     0.88762469649636953125E+13
       coef(  653) =    -0.41921921815451617188E+14
       coef(  654) =     0.80346276106915625000E+14
       coef(  655) =     0.61180254400743460938E+14
       coef(  656) =    -0.32046817840215609375E+14
       coef(  657) =     0.86141197300667312500E+14
       coef(  658) =    -0.92787458065632281250E+14
       coef(  659) =    -0.30044872443435070312E+14
       coef(  660) =     0.55652974776773300781E+13
       coef(  661) =     0.26758852159807195312E+14
       coef(  662) =    -0.55322301069014695312E+14
       coef(  663) =    -0.29279389685792109375E+14
       coef(  664) =    -0.38828398238644835938E+14
       coef(  665) =    -0.20954941963083549500E+11
       coef(  666) =     0.63971409457692431641E+12
       coef(  667) =    -0.41423749949316142578E+13
       coef(  668) =    -0.14915385425667757812E+14
       coef(  669) =     0.61336282145876382812E+14
       coef(  670) =     0.45779602574926828125E+14
       coef(  671) =    -0.40873429215720796875E+14
       coef(  672) =    -0.33656969309942998047E+13
       coef(  673) =     0.43452799960449593750E+14
       coef(  674) =    -0.11349076148295629688E+15
       coef(  675) =    -0.20869321091968164062E+12
       coef(  676) =     0.16021015487871557617E+13
       coef(  677) =    -0.14366787685644462891E+14
       coef(  678) =    -0.43166961517377490234E+12
       coef(  679) =     0.44591434457440203125E+14
       coef(  680) =     0.47478445014741203125E+14
       coef(  681) =     0.21570373485435097656E+14
       coef(  682) =     0.29740511016197113281E+14
       coef(  683) =     0.15352674970527355469E+14
       coef(  684) =     0.74349133646081699219E+13
       coef(  685) =    -0.37833942415991781250E+14
       coef(  686) =    -0.91062319570283144531E+13
       coef(  687) =    -0.14047692773011906250E+14
       coef(  688) =    -0.15609667787291835938E+14
       coef(  689) =    -0.60036258021761484375E+14
       coef(  690) =    -0.18610291349199281250E+14
       coef(  691) =    -0.56706585502348632812E+13
       coef(  692) =    -0.72604789190272724609E+13
       coef(  693) =    -0.52842718888350320312E+14
       coef(  694) =    -0.18873352399090128906E+14
       coef(  695) =     0.67459584315578483582E+11
       coef(  696) =    -0.22124408106923505859E+13
       coef(  697) =     0.28717279067799964844E+14
       coef(  698) =    -0.77081514178236406250E+14
       coef(  699) =     0.26747018552145968750E+14
       coef(  700) =    -0.61933683632103046875E+13
       coef(  701) =    -0.21953062721118386719E+14
       coef(  702) =    -0.62779367599345000000E+13
       coef(  703) =     0.13363015808568933594E+14
       coef(  704) =     0.73714573813381468750E+14
       coef(  705) =     0.49377720634348359375E+14
       coef(  706) =     0.10528212804367052734E+14
       coef(  707) =     0.63105658327235429688E+14
       coef(  708) =     0.47172939041745625000E+14
       coef(  709) =     0.23628204433379203125E+14
       coef(  710) =     0.17425934714151634766E+14
       coef(  711) =    -0.17412793174381609375E+14
       coef(  712) =     0.34807908440421593750E+14
       coef(  713) =     0.33628004914286046875E+14
       coef(  714) =     0.13956747256854373047E+14
       coef(  715) =     0.67303681122946083984E+13
       coef(  716) =     0.44400676955783603516E+13
       coef(  717) =    -0.55277365882265175781E+13
       coef(  718) =    -0.25946332355540562500E+14
       coef(  719) =    -0.97702328102343046875E+13
       coef(  720) =     0.24914923821948720703E+13
       coef(  721) =    -0.82792376884500371094E+13
       coef(  722) =    -0.45313015048374421895E+06
       coef(  723) =     0.18026378482919526100E+09
       coef(  724) =    -0.10632251417157266617E+11
       coef(  725) =     0.24313253119098367310E+12
       coef(  726) =    -0.22871282238422861328E+13
       coef(  727) =     0.10491281109301746094E+14
       coef(  728) =    -0.21803846088048421875E+14
       coef(  729) =     0.17166679109512156250E+14
       coef(  730) =     0.67902737992296781540E+10
       coef(  731) =    -0.29881149072151068115E+12
       coef(  732) =     0.34620124238688681641E+13
       coef(  733) =    -0.15633758910267568359E+14
       coef(  734) =     0.69878199133358730469E+13
       coef(  735) =     0.46615086691062648438E+14
       coef(  736) =    -0.77862856188213906250E+13
       coef(  737) =     0.93681110626847686768E+11
       coef(  738) =    -0.43700740588072465820E+13
       coef(  739) =    -0.26668442204637164062E+14
       coef(  740) =     0.25485384829594937500E+14
       coef(  741) =    -0.13493946155086953125E+13
       coef(  742) =    -0.14531109835957484375E+14
       coef(  743) =     0.55198094483173039062E+14
       coef(  744) =    -0.38054917399788257812E+14
       coef(  745) =    -0.47873466608640179688E+14
       coef(  746) =    -0.25385516347504406250E+14
       coef(  747) =    -0.48619149709857851562E+14
       coef(  748) =    -0.33528954486165042969E+14
       coef(  749) =    -0.24260619648080497742E+11
       coef(  750) =     0.14125821648716254883E+13
       coef(  751) =     0.37180136742920214844E+13
       coef(  752) =     0.46898072585585382812E+14
       coef(  753) =    -0.81664307611515062500E+14
       coef(  754) =    -0.77863468662108359375E+14
       coef(  755) =    -0.55342139987929234375E+14
       coef(  756) =    -0.16907063363893300781E+14
       coef(  757) =    -0.10011781019930447266E+14
       coef(  758) =     0.73457075606287375000E+14
       coef(  759) =     0.23776401831587199219E+14
       coef(  760) =    -0.54327321180067167969E+13
       coef(  761) =    -0.57973317345873093750E+14
       coef(  762) =     0.15577454760355285156E+14
       coef(  763) =     0.11236286312381140625E+14
       coef(  764) =     0.90377013365590996094E+13
       coef(  765) =     0.16262650213658712891E+14
       coef(  766) =     0.45926360025127546875E+14
       coef(  767) =     0.52765208887692625000E+14
       coef(  768) =     0.16241555956359935547E+14
       coef(  769) =     0.17554043128318598633E+13
       coef(  770) =     0.76722497321268076172E+13
       coef(  771) =    -0.62236789245903437500E+13
       coef(  772) =    -0.92436621648287573242E+12
       coef(  773) =     0.15281648844716855469E+14
       coef(  774) =    -0.80771499880096765625E+14
       coef(  775) =     0.94741059205045218750E+14
       coef(  776) =     0.16450331919308443359E+14
       coef(  777) =    -0.17682035957910234375E+14
       coef(  778) =    -0.44414333741613320312E+13
       coef(  779) =    -0.17462071374838785156E+14
       coef(  780) =     0.39450131783130218750E+14
       coef(  781) =     0.16378031954150498047E+14
       coef(  782) =    -0.59971536085725859375E+13
       coef(  783) =     0.94513221718251503906E+13
       coef(  784) =    -0.10399671353263218750E+14
       coef(  785) =     0.73638991271557832031E+13
       coef(  786) =     0.13150925624907601562E+14
       coef(  787) =     0.18874565949415830078E+13
       coef(  788) =    -0.12281393357134191406E+14
       coef(  789) =    -0.13288664973672236328E+13
       coef(  790) =    -0.22978859190490759909E+08
       coef(  791) =    -0.42574714674619569778E+10
       coef(  792) =     0.18731153628698971558E+12
       coef(  793) =    -0.34980381867500932617E+13
       coef(  794) =     0.22347062189240214844E+14
       coef(  795) =    -0.54559402823019984375E+14
       coef(  796) =     0.56770892103006125000E+14
       coef(  797) =    -0.55397478351484578125E+14
       coef(  798) =     0.15108933536827711105E+11
       coef(  799) =     0.14980206674130068359E+13
       coef(  800) =    -0.19874955037363750000E+14
       coef(  801) =     0.67084471747104945312E+14
       coef(  802) =     0.59963712512539101562E+13
       coef(  803) =    -0.31147139016306507812E+14
       coef(  804) =    -0.41362068040238367188E+14
       coef(  805) =     0.59736981341227021484E+13
       coef(  806) =    -0.36634768421689906250E+14
       coef(  807) =     0.94235095120879781250E+14
       coef(  808) =     0.59306762714551117188E+14
       coef(  809) =     0.28061599695821845703E+13
       coef(  810) =    -0.10265463420531935938E+15
       coef(  811) =    -0.26077644925529597656E+14
       coef(  812) =    -0.76029013131785439453E+13
       coef(  813) =    -0.17647312177303164062E+14
       coef(  814) =     0.46855397743218457031E+12
       coef(  815) =    -0.20232365173469898438E+14
       coef(  816) =     0.13824393826516818359E+14
       coef(  817) =    -0.88057014408553437500E+14
       coef(  818) =    -0.90386882270599515625E+14
       coef(  819) =    -0.51918555592551148438E+14
       coef(  820) =     0.13159746796427639062E+15
       coef(  821) =    -0.66679261550335250000E+14
       coef(  822) =    -0.12009275125877769531E+14
       coef(  823) =    -0.66335686280349804688E+13
       coef(  824) =    -0.62025308784277218750E+14
       coef(  825) =    -0.89765934011246523438E+13
       coef(  826) =    -0.93886526199765625000E+14
       coef(  827) =    -0.49623142798164351562E+14
       coef(  828) =    -0.37760682108233056641E+13
       coef(  829) =    -0.21112548495864460938E+14
       coef(  830) =     0.72513222645646923828E+13
       coef(  831) =    -0.64433709553508460938E+14
       coef(  832) =     0.10943294490727696875E+15
       coef(  833) =     0.65961286225018656250E+14
       coef(  834) =     0.50179690108639853516E+13
       coef(  835) =     0.12154496112259400000E+15
       coef(  836) =     0.29449586931428214844E+14
       coef(  837) =     0.15045208167615060547E+14
       coef(  838) =    -0.15828611176285839844E+12
       coef(  839) =     0.53243204288059414062E+14
       coef(  840) =     0.11253202793599537109E+14
       coef(  841) =     0.15934376350846097656E+14
       coef(  842) =     0.69452393350711917877E+09
       coef(  843) =    -0.67832202439179821014E+10
       coef(  844) =    -0.46279388027904827881E+12
       coef(  845) =     0.12522137987377496094E+14
       coef(  846) =    -0.65348750717632226562E+14
       coef(  847) =     0.54830855754472726562E+14
       coef(  848) =     0.50561150496063820312E+14
       coef(  849) =    -0.64094969280667968750E+12
       coef(  850) =    -0.41162117703335003662E+12
       coef(  851) =    -0.15720007681752392578E+12
       coef(  852) =     0.31551365698416367188E+14
       coef(  853) =    -0.59648865761048335938E+14
       coef(  854) =    -0.50293223587485867188E+14
       coef(  855) =    -0.36509138002015796875E+14
       coef(  856) =    -0.66329033684405898438E+13
       coef(  857) =     0.51904970440012539062E+14
       coef(  858) =     0.52884547203810921875E+14
       coef(  859) =     0.23000928295800015625E+14
       coef(  860) =    -0.39874731267661859375E+14
       coef(  861) =    -0.10225948151543421875E+14
       coef(  862) =     0.15219724216272802734E+12
       coef(  863) =     0.16839065523812957031E+14
       coef(  864) =    -0.13394617507535939453E+14
       coef(  865) =    -0.62972233861574179688E+14
       coef(  866) =    -0.42884235848579984375E+14
       coef(  867) =    -0.88633438318772468750E+14
       coef(  868) =    -0.51920549503023414062E+14
       coef(  869) =    -0.13856425258423367188E+14
       coef(  870) =    -0.27662614504577089844E+14
       coef(  871) =    -0.66863343775162335938E+14
       coef(  872) =    -0.25972413146176578125E+14
       coef(  873) =    -0.16367608992396664062E+14
       coef(  874) =    -0.32066280453177976562E+14
       coef(  875) =     0.79391656267864562500E+14
       coef(  876) =     0.31585191155445468750E+14
       coef(  877) =     0.54659359676932710938E+14
       coef(  878) =     0.18778162998472500000E+14
       coef(  879) =     0.29786537389303496094E+14
       coef(  880) =    -0.44985190842215833664E+10
       coef(  881) =     0.20409239167772659302E+12
       coef(  882) =    -0.16680159017354306641E+13
       coef(  883) =    -0.14002011138972662109E+14
       coef(  884) =     0.10082532643005684375E+15
       coef(  885) =    -0.15739578346219457031E+14
       coef(  886) =    -0.10214023877490357422E+14
       coef(  887) =     0.10389763913495168457E+13
       coef(  888) =    -0.81230776922868291016E+13
       coef(  889) =    -0.29611150556234460938E+14
       coef(  890) =    -0.72083656706317880859E+13
       coef(  891) =    -0.14471814920160150391E+14
       coef(  892) =    -0.19690111894452945312E+14
       coef(  893) =     0.83929422561480343750E+14
       coef(  894) =     0.46070972167062117188E+14
       coef(  895) =     0.60898347355464140625E+13
       coef(  896) =     0.24266266077147128906E+13
       coef(  897) =     0.32475151475650210938E+14
       coef(  898) =     0.40401678589961617188E+14
       coef(  899) =    -0.18244238737408176270E+13
       coef(  900) =    -0.56882982178183750000E+14
       coef(  901) =    -0.76328540888089257812E+13
       coef(  902) =    -0.18393071725777875000E+14
       coef(  903) =     0.48332503829825148438E+14
       coef(  904) =    -0.25963179994402324219E+13
       coef(  905) =     0.47945978496611812500E+14
       coef(  906) =     0.18936204043764253906E+14
       coef(  907) =     0.86808937299732208252E+10
       coef(  908) =    -0.63303237382690649414E+12
       coef(  909) =     0.93652568995979238281E+13
       coef(  910) =    -0.26644998241188574219E+14
       coef(  911) =    -0.19376753847234472656E+14
       coef(  912) =    -0.91047144842303875000E+14
       coef(  913) =    -0.10998839137846115723E+13
       coef(  914) =     0.29122581691012515625E+14
       coef(  915) =    -0.10356428150267398438E+14
       coef(  916) =     0.52702691261554238281E+13
       coef(  917) =    -0.12635492626800445312E+14
       coef(  918) =     0.60087212887437195312E+14
       coef(  919) =    -0.29799109948210031250E+14
       coef(  920) =     0.96511491041920781250E+13
       coef(  921) =     0.40739823936467703125E+14
       coef(  922) =    -0.25443254053217125000E+14
       coef(  923) =    -0.30996307573899437500E+14
       coef(  924) =    -0.13712311457146464844E+14
       coef(  925) =    -0.98146945231986212730E+09
       coef(  926) =     0.44787792824014117432E+12
       coef(  927) =    -0.11011391866019251953E+14
       coef(  928) =     0.50666490275137968750E+14
       coef(  929) =    -0.36470344193770328125E+14
       coef(  930) =     0.23071917563572827148E+13
       coef(  931) =    -0.27711585492694902344E+14
       coef(  932) =    -0.61788303255837451172E+13
       coef(  933) =    -0.12042614731688921875E+14
       coef(  934) =     0.32194230866961968750E+14
       coef(  935) =     0.97456321367475156250E+13
       coef(  936) =    -0.37890892554588210938E+14
       coef(  937) =    -0.31738547158569843305E+05
       coef(  938) =     0.42054412580029210076E+07
       coef(  939) =    -0.47125379927866917849E+09
       coef(  940) =     0.15159512397988882065E+11
       coef(  941) =    -0.22018088179040344238E+12
       coef(  942) =     0.12300503305558466797E+13
       coef(  943) =    -0.26019126119306152344E+13
       coef(  944) =     0.21952798202495068359E+13
       coef(  945) =     0.42683139818230915070E+09
       coef(  946) =    -0.20848756104534301758E+11
       coef(  947) =     0.58133792321973828125E+12
       coef(  948) =    -0.36339653553132529297E+13
       coef(  949) =     0.12521281986020289062E+14
       coef(  950) =    -0.36286508029533273438E+14
       coef(  951) =     0.37041946098150789062E+14
       coef(  952) =    -0.64978690020646972656E+11
       coef(  953) =    -0.13955518990247879028E+12
       coef(  954) =    -0.13735842340358138672E+14
       coef(  955) =     0.47456062396270562500E+14
       coef(  956) =    -0.38794191384636906250E+14
       coef(  957) =    -0.12734524455051626953E+14
       coef(  958) =     0.79891819176273339844E+13
       coef(  959) =    -0.16461373910427392578E+12
       coef(  960) =    -0.65299447516990693359E+13
       coef(  961) =     0.12733235631954398438E+14
       coef(  962) =     0.24880030018755144531E+14
       coef(  963) =    -0.19597251500624449219E+14
       coef(  964) =    -0.14995226584389115234E+14
       coef(  965) =     0.27748438553759033203E+12
       coef(  966) =    -0.84481714942232402344E+13
       coef(  967) =    -0.12770863108166017532E+10
       coef(  968) =    -0.22218202348978399658E+12
       coef(  969) =    -0.88932429363683020020E+12
       coef(  970) =     0.16623807289270767578E+14
       coef(  971) =    -0.25552493258775785156E+14
       coef(  972) =     0.26158798392991437500E+14
       coef(  973) =     0.13514002872883027344E+13
       coef(  974) =     0.11972479097998637695E+13
       coef(  975) =     0.49807323582915703125E+13
       coef(  976) =    -0.35183024356304375000E+14
       coef(  977) =    -0.14950574577384828125E+14
       coef(  978) =    -0.14291378789394056641E+14
       coef(  979) =    -0.13140518510868119141E+14
       coef(  980) =    -0.44667815641372625000E+14
       coef(  981) =     0.43504489339297773438E+14
       coef(  982) =     0.53463403180329296875E+14
       coef(  983) =     0.26848579720334414062E+14
       coef(  984) =     0.37943397994279171875E+14
       coef(  985) =     0.24515077273785312500E+14
       coef(  986) =    -0.92569070606494492188E+13
       coef(  987) =     0.57492258096935578125E+14
       coef(  988) =     0.60514144968449140625E+14
       coef(  989) =     0.11567428033101386719E+13
       coef(  990) =    -0.12637962338054023438E+14
       coef(  991) =    -0.17952828840392558594E+14
       coef(  992) =     0.16672961486532866211E+13
       coef(  993) =     0.13558918664894707031E+13
       coef(  994) =     0.80564472554525476074E+12
       coef(  995) =    -0.38655531176048101562E+14
       coef(  996) =    -0.13323032885712220703E+14
       coef(  997) =     0.23699314078312768555E+12
       coef(  998) =    -0.31322236979455205078E+13
       coef(  999) =     0.93652329908017597656E+13
       coef( 1000) =    -0.32283596446630214844E+14
       coef( 1001) =     0.85131780734271738281E+13
       coef( 1002) =    -0.32664055810761689453E+13
       coef( 1003) =    -0.54794343996273496094E+13
       coef( 1004) =     0.18225537915535906250E+14
       coef( 1005) =    -0.97780815019481750000E+14
       coef( 1006) =     0.81821500211234781250E+14
       coef( 1007) =     0.57093551929168921875E+14
       coef( 1008) =     0.15664680753598951172E+14
       coef( 1009) =     0.27079212436754796875E+14
       coef( 1010) =     0.54071766108181132812E+14
       coef( 1011) =     0.29910668816570121094E+14
       coef( 1012) =     0.24459364955183535156E+14
       coef( 1013) =    -0.41199569297914746094E+13
       coef( 1014) =     0.60919239875842000000E+14
       coef( 1015) =     0.58238261524922804688E+14
       coef( 1016) =     0.23548926914673750000E+14
       coef( 1017) =     0.24913838363336941406E+14
       coef( 1018) =     0.13265104457753228516E+14
       coef( 1019) =     0.18999406390504450684E+13
       coef( 1020) =    -0.19533633441447687500E+14
       coef( 1021) =     0.40539221013887666016E+13
       coef( 1022) =     0.11101423984636871094E+14
       coef( 1023) =    -0.26305864594908198242E+13
       coef( 1024) =    -0.45313015048374421895E+06
       coef( 1025) =     0.18026378482919526100E+09
       coef( 1026) =    -0.10632251417157266617E+11
       coef( 1027) =     0.24313253119098367310E+12
       coef( 1028) =    -0.22871282238422861328E+13
       coef( 1029) =     0.10491281109301746094E+14
       coef( 1030) =    -0.21803846088048421875E+14
       coef( 1031) =     0.17166679109512156250E+14
       coef( 1032) =     0.67902737992296781540E+10
       coef( 1033) =    -0.29881149072151068115E+12
       coef( 1034) =     0.34620124238688681641E+13
       coef( 1035) =    -0.15633758910267568359E+14
       coef( 1036) =     0.69878199133358730469E+13
       coef( 1037) =     0.46615086691062648438E+14
       coef( 1038) =    -0.77862856188213906250E+13
       coef( 1039) =     0.93681110626847686768E+11
       coef( 1040) =    -0.43700740588072465820E+13
       coef( 1041) =    -0.26668442204637164062E+14
       coef( 1042) =     0.25485384829594937500E+14
       coef( 1043) =    -0.13493946155086953125E+13
       coef( 1044) =    -0.14531109835957484375E+14
       coef( 1045) =     0.55198094483173039062E+14
       coef( 1046) =    -0.38054917399788257812E+14
       coef( 1047) =    -0.47873466608640179688E+14
       coef( 1048) =    -0.25385516347504406250E+14
       coef( 1049) =    -0.48619149709857851562E+14
       coef( 1050) =    -0.33528954486165042969E+14
       coef( 1051) =    -0.24260619648080497742E+11
       coef( 1052) =     0.14125821648716254883E+13
       coef( 1053) =     0.37180136742920214844E+13
       coef( 1054) =     0.46898072585585382812E+14
       coef( 1055) =    -0.81664307611515062500E+14
       coef( 1056) =    -0.77863468662108359375E+14
       coef( 1057) =    -0.55342139987929234375E+14
       coef( 1058) =    -0.16907063363893300781E+14
       coef( 1059) =    -0.10011781019930447266E+14
       coef( 1060) =     0.73457075606287375000E+14
       coef( 1061) =     0.23776401831587199219E+14
       coef( 1062) =    -0.54327321180067167969E+13
       coef( 1063) =    -0.57973317345873093750E+14
       coef( 1064) =     0.15577454760355285156E+14
       coef( 1065) =     0.11236286312381140625E+14
       coef( 1066) =     0.90377013365590996094E+13
       coef( 1067) =     0.16262650213658712891E+14
       coef( 1068) =     0.45926360025127546875E+14
       coef( 1069) =     0.52765208887692625000E+14
       coef( 1070) =     0.16241555956359935547E+14
       coef( 1071) =     0.17554043128318598633E+13
       coef( 1072) =     0.76722497321268076172E+13
       coef( 1073) =    -0.62236789245903437500E+13
       coef( 1074) =    -0.92436621648287573242E+12
       coef( 1075) =     0.15281648844716855469E+14
       coef( 1076) =    -0.80771499880096765625E+14
       coef( 1077) =     0.94741059205045218750E+14
       coef( 1078) =     0.16450331919308443359E+14
       coef( 1079) =    -0.17682035957910234375E+14
       coef( 1080) =    -0.44414333741613320312E+13
       coef( 1081) =    -0.17462071374838785156E+14
       coef( 1082) =     0.39450131783130218750E+14
       coef( 1083) =     0.16378031954150498047E+14
       coef( 1084) =    -0.59971536085725859375E+13

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

      norder=7
      norder2=7
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

