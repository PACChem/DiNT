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
       coef(    1) =    -0.22232291764969556169E+01
       coef(    2) =     0.25001323570273292773E+02
       coef(    3) =    -0.51450292729420514661E+05
       coef(    4) =     0.13090531469563287683E+07
       coef(    5) =    -0.59382658081002220511E+08
       coef(    6) =    -0.67118962296864187717E+09
       coef(    7) =     0.15575111497311908722E+11
       coef(    8) =    -0.92134865494977218628E+11
       coef(    9) =     0.18914575195986212158E+12
       coef(   10) =     0.33975069289802588173E+05
       coef(   11) =     0.62129924456147952005E+07
       coef(   12) =     0.12557240907003444433E+09
       coef(   13) =     0.11615212340012817383E+11
       coef(   14) =    -0.14970929255776043701E+12
       coef(   15) =     0.59622073068294018555E+12
       coef(   16) =    -0.33114366695379394531E+11
       coef(   17) =    -0.34090761241814423828E+13
       coef(   18) =    -0.47225144840597969294E+09
       coef(   19) =     0.85119136788860321045E+10
       coef(   20) =    -0.12492906090550352478E+12
       coef(   21) =     0.15319860146298811035E+13
       coef(   22) =    -0.10540880330327369141E+14
       coef(   23) =     0.29113729778963523438E+14
       coef(   24) =    -0.20053419520833796875E+14
       coef(   25) =    -0.12971740501539068604E+12
       coef(   26) =     0.21398135258446262207E+13
       coef(   27) =    -0.73764657396531396484E+13
       coef(   28) =     0.30201998732886953125E+14
       coef(   29) =    -0.89149931863095406250E+14
       coef(   30) =     0.75729164034552156250E+14
       coef(   31) =    -0.79373950019877011719E+13
       coef(   32) =     0.32888550214418101562E+14
       coef(   33) =    -0.16670745254470625000E+14
       coef(   34) =    -0.24529362101213437500E+13
       coef(   35) =     0.53279714124618250000E+14
       coef(   36) =    -0.24776563159177609375E+14
       coef(   37) =    -0.70488454668613320312E+13
       coef(   38) =     0.37520474109140515625E+14
       coef(   39) =     0.44360307485777867188E+14
       coef(   40) =    -0.70604158376065976562E+13
       coef(   41) =     0.12061603448633908203E+14
       coef(   42) =    -0.77401102413144344464E+07
       coef(   43) =     0.29560444096442341805E+09
       coef(   44) =    -0.38452767341301139832E+11
       coef(   45) =     0.27383372341365252686E+12
       coef(   46) =    -0.17845735282848632812E+12
       coef(   47) =    -0.30379979637085556641E+13
       coef(   48) =    -0.98899839730895605469E+12
       coef(   49) =     0.35525114454123648438E+14
       coef(   50) =     0.20387369384136478424E+11
       coef(   51) =     0.48618839658765222168E+12
       coef(   52) =    -0.59666225220060214844E+13
       coef(   53) =     0.24898294550342492188E+14
       coef(   54) =    -0.19892422882165433594E+14
       coef(   55) =    -0.56606447018584976562E+14
       coef(   56) =    -0.57608166533597921875E+14
       coef(   57) =    -0.25155222084377304688E+13
       coef(   58) =     0.17957730059564085938E+14
       coef(   59) =    -0.49983328343898859375E+14
       coef(   60) =     0.53260684990110695312E+14
       coef(   61) =     0.86313092660731025391E+13
       coef(   62) =    -0.46921873297968906250E+13
       coef(   63) =     0.37499004784290117188E+13
       coef(   64) =    -0.24808004702035582031E+14
       coef(   65) =     0.23480910923872644531E+14
       coef(   66) =     0.18105162793897500000E+14
       coef(   67) =     0.91544444179497304688E+13
       coef(   68) =    -0.69391560549945859375E+14
       coef(   69) =    -0.24714356319673132812E+14
       coef(   70) =    -0.14079577282426967773E+13
       coef(   71) =    -0.14809929535041175781E+14
       coef(   72) =    -0.40757870154388562012E+12
       coef(   73) =     0.48150214357997548828E+13
       coef(   74) =    -0.46585111702972031250E+13
       coef(   75) =    -0.30563793710218921875E+14
       coef(   76) =     0.90286911001436625000E+14
       coef(   77) =     0.41714007391253226562E+14
       coef(   78) =    -0.17635782032466250000E+13
       coef(   79) =    -0.12813290006640023438E+14
       coef(   80) =     0.24100131388158914062E+14
       coef(   81) =    -0.26090956727248058594E+14
       coef(   82) =     0.91124323201773085938E+13
       coef(   83) =     0.42282853850795175781E+13
       coef(   84) =    -0.20735421309373652344E+13
       coef(   85) =     0.59733894278241343750E+14
       coef(   86) =     0.30584653076184570312E+13
       coef(   87) =     0.34734582512022587891E+13
       coef(   88) =     0.26202558818671547852E+13
       coef(   89) =    -0.10511908950406312500E+14
       coef(   90) =    -0.42471155270092812500E+13
       coef(   91) =     0.81202219495477695312E+13
       coef(   92) =    -0.35797959407740023438E+14
       coef(   93) =    -0.44157648381066992188E+14
       coef(   94) =    -0.15981548663658769531E+14
       coef(   95) =    -0.53328679935803925781E+13
       coef(   96) =    -0.25841083976486882812E+14
       coef(   97) =    -0.17081580318599029297E+14
       coef(   98) =    -0.58512024225275068359E+13
       coef(   99) =    -0.72069142338992294922E+13
       coef(  100) =    -0.16055532150485941406E+14
       coef(  101) =    -0.68923256063157324219E+13
       coef(  102) =     0.14085011751662456989E+09
       coef(  103) =    -0.26620013242844672203E+10
       coef(  104) =     0.80066513877163177490E+11
       coef(  105) =    -0.20732746444953598633E+13
       coef(  106) =     0.10777864048372740234E+14
       coef(  107) =    -0.15095354407905232422E+14
       coef(  108) =     0.64296027262695234375E+13
       coef(  109) =    -0.75188633677155781250E+14
       coef(  110) =     0.21066309766146087646E+11
       coef(  111) =     0.31009263357761533203E+13
       coef(  112) =    -0.56561113915845751953E+13
       coef(  113) =    -0.54464565233006687500E+14
       coef(  114) =     0.12725504999710345312E+15
       coef(  115) =     0.12869671417815339062E+15
       coef(  116) =     0.66075221309470453125E+14
       coef(  117) =     0.52669602494679335938E+13
       coef(  118) =    -0.46884841900021640625E+14
       coef(  119) =     0.43761359985687078125E+14
       coef(  120) =     0.94937574866106828125E+14
       coef(  121) =     0.65122844532831742188E+14
       coef(  122) =     0.30971954990927800781E+14
       coef(  123) =     0.47792524051406484375E+14
       coef(  124) =     0.30495297960904949219E+14
       coef(  125) =     0.32265282813395222656E+14
       coef(  126) =     0.20971126403903414062E+14
       coef(  127) =    -0.41594682498765742188E+13
       coef(  128) =    -0.36116777575553271484E+12
       coef(  129) =    -0.29269717164735361328E+13
       coef(  130) =    -0.41849101036197578125E+13
       coef(  131) =     0.12176603801118493750E+15
       coef(  132) =    -0.21700022707804828125E+15
       coef(  133) =    -0.12659035389846115625E+15
       coef(  134) =    -0.34065851761506816406E+14
       coef(  135) =    -0.32765643913073857422E+13
       coef(  136) =    -0.31092815016148234375E+14
       coef(  137) =    -0.31354772019878914062E+14
       coef(  138) =    -0.92839404633154015625E+14
       coef(  139) =    -0.46830495589114335938E+14
       coef(  140) =    -0.13915348968764484375E+14
       coef(  141) =    -0.20969204325381214844E+14
       coef(  142) =    -0.24615292417770843750E+14
       coef(  143) =    -0.10418261528792421875E+14
       coef(  144) =    -0.10661664415286097656E+14
       coef(  145) =     0.13835405820423360938E+15
       coef(  146) =     0.24980925765504234375E+14
       coef(  147) =    -0.23890008300548113281E+14
       coef(  148) =    -0.15349190419059404297E+14
       coef(  149) =     0.40435016569588657227E+13
       coef(  150) =    -0.53383580822352011719E+13
       coef(  151) =     0.10767333712036154785E+13
       coef(  152) =     0.14949497278642906250E+14
       coef(  153) =    -0.93288921106234859375E+14
       coef(  154) =     0.17306287027682050000E+15
       coef(  155) =    -0.36999595225851945312E+14
       coef(  156) =    -0.47780834831599437500E+14
       coef(  157) =    -0.19098020129707109375E+14
       coef(  158) =    -0.20397444889414812500E+14
       coef(  159) =     0.19657783422746578125E+14
       coef(  160) =    -0.19800736206310976562E+14
       coef(  161) =    -0.15120519610734632812E+14
       coef(  162) =     0.14613360561370439453E+13
       coef(  163) =    -0.57393244866626621094E+13
       coef(  164) =     0.19324627341890320312E+14
       coef(  165) =     0.91069092153326445312E+13
       coef(  166) =    -0.40832443798102519531E+13
       coef(  167) =     0.22054526621384062500E+13
       coef(  168) =     0.10822976543787191406E+14
       coef(  169) =     0.34118117170501918945E+13
       coef(  170) =     0.18802491711596800883E+03
       coef(  171) =    -0.21557593047847134585E+05
       coef(  172) =     0.36922546765964906663E+07
       coef(  173) =    -0.75708329405135363340E+08
       coef(  174) =     0.20908119374217214584E+10
       coef(  175) =    -0.11862161101039104462E+10
       coef(  176) =    -0.13917980900587203979E+11
       coef(  177) =    -0.17520175167806051636E+12
       coef(  178) =     0.29374042909338964844E+12
       coef(  179) =    -0.45817341064970416483E+06
       coef(  180) =    -0.39593856336661159992E+09
       coef(  181) =    -0.39606286576935682297E+10
       coef(  182) =    -0.14111184103405902100E+12
       coef(  183) =     0.70926927273435205078E+12
       coef(  184) =     0.13060937153189978027E+12
       coef(  185) =     0.17097882145942050781E+13
       coef(  186) =    -0.60518585288168281250E+13
       coef(  187) =     0.10364152778042295456E+11
       coef(  188) =    -0.85943357110941925049E+11
       coef(  189) =     0.20008598806003427734E+13
       coef(  190) =    -0.15014299724212718750E+14
       coef(  191) =     0.47013836562495093750E+14
       coef(  192) =    -0.89678576086983593750E+14
       coef(  193) =     0.98716562902856546875E+14
       coef(  194) =    -0.32282920925675354004E+10
       coef(  195) =    -0.10940818292134243164E+13
       coef(  196) =     0.17956821781936171875E+13
       coef(  197) =    -0.54171458325521562500E+13
       coef(  198) =    -0.46371279575556140625E+14
       coef(  199) =     0.13853839550335960938E+14
       coef(  200) =     0.29633489096456982422E+13
       coef(  201) =     0.63840554395166953125E+13
       coef(  202) =     0.87673959199780312500E+13
       coef(  203) =    -0.94003011950116796875E+12
       coef(  204) =     0.68156268573799804688E+13
       coef(  205) =    -0.46830686635990781250E+13
       coef(  206) =    -0.38863667024977380371E+12
       coef(  207) =     0.11412308377795908203E+13
       coef(  208) =    -0.14332313745977785645E+13
       coef(  209) =     0.33500599817549681664E+09
       coef(  210) =     0.78286991099532394409E+10
       coef(  211) =     0.42710435410386999512E+12
       coef(  212) =     0.13132975190657614136E+12
       coef(  213) =    -0.10878173889271841797E+14
       coef(  214) =     0.66891481206554589844E+13
       coef(  215) =     0.34537030250868214844E+14
       coef(  216) =    -0.35311823309791609375E+14
       coef(  217) =    -0.53241368353433630371E+12
       coef(  218) =    -0.45583768724576191406E+13
       coef(  219) =     0.25256224452207554688E+14
       coef(  220) =     0.29184912503419414062E+13
       coef(  221) =    -0.31243571178804234375E+14
       coef(  222) =    -0.52711889992671046875E+14
       coef(  223) =    -0.31937604266776937500E+14
       coef(  224) =     0.21053954408398687500E+14
       coef(  225) =    -0.10157063260106950000E+15
       coef(  226) =     0.12645361979872101562E+15
       coef(  227) =     0.88759983708423296875E+14
       coef(  228) =     0.26711984209672171875E+14
       coef(  229) =     0.61654542394411318359E+13
       coef(  230) =    -0.32496167087814652344E+14
       coef(  231) =     0.52949768458660296875E+14
       coef(  232) =     0.38685814863400882812E+14
       coef(  233) =     0.15198312988942867188E+14
       coef(  234) =     0.26303287555177609375E+14
       coef(  235) =     0.12019336658956906250E+14
       coef(  236) =     0.63351611760779218750E+13
       coef(  237) =    -0.52277825165058859375E+14
       coef(  238) =     0.24832810934653281250E+14
       coef(  239) =    -0.25385405079289843750E+14
       coef(  240) =    -0.45295666025850281250E+14
       coef(  241) =    -0.34978998427363992188E+14
       coef(  242) =    -0.17503849577984671875E+14
       coef(  243) =     0.12538044660750881250E+15
       coef(  244) =    -0.33569801440599351562E+14
       coef(  245) =    -0.17882854909708621094E+14
       coef(  246) =    -0.82020013688669218750E+13
       coef(  247) =    -0.51043474502789990234E+13
       coef(  248) =    -0.41213538822876570312E+14
       coef(  249) =    -0.11317296846544796875E+14
       coef(  250) =    -0.16247725551511323242E+13
       coef(  251) =    -0.28983244201937158203E+13
       coef(  252) =    -0.64127910038529875000E+14
       coef(  253) =    -0.72554646669517796875E+14
       coef(  254) =    -0.27473352634474781250E+14
       coef(  255) =    -0.83822149627682324219E+13
       coef(  256) =    -0.36505492154402945312E+14
       coef(  257) =    -0.11090529055975251953E+14
       coef(  258) =    -0.14057028055598253906E+14
       coef(  259) =    -0.11490185403009384155E+11
       coef(  260) =     0.12989783457056585693E+12
       coef(  261) =    -0.35562211837026552734E+13
       coef(  262) =     0.12104188929075646484E+14
       coef(  263) =     0.24009379851835742188E+12
       coef(  264) =    -0.22787061102518519531E+14
       coef(  265) =     0.53474245312642636719E+13
       coef(  266) =     0.10810624260862296875E+14
       coef(  267) =     0.32592092383789106445E+13
       coef(  268) =     0.74727018371513427734E+13
       coef(  269) =    -0.37754511890597648438E+14
       coef(  270) =     0.24131058156366558594E+14
       coef(  271) =    -0.75277966954720625000E+13
       coef(  272) =    -0.11268409694102830078E+14
       coef(  273) =    -0.44147550519760439453E+13
       coef(  274) =    -0.26461029764760523438E+14
       coef(  275) =     0.46330791204053535156E+13
       coef(  276) =     0.40067635684805085938E+14
       coef(  277) =     0.19568992009772406250E+14
       coef(  278) =     0.57990008716352509766E+13
       coef(  279) =     0.83776288933271621094E+13
       coef(  280) =     0.13925018726195500000E+14
       coef(  281) =     0.76705565432788535156E+13
       coef(  282) =     0.52247893989886289062E+13
       coef(  283) =    -0.29237422494630921875E+14
       coef(  284) =     0.78122497977807453125E+14
       coef(  285) =     0.27887045898179484375E+14
       coef(  286) =    -0.41163219651525957031E+13
       coef(  287) =    -0.14644217869066046875E+14
       coef(  288) =    -0.93281355650760546875E+13
       coef(  289) =    -0.64358921023884187500E+14
       coef(  290) =    -0.34433756383383570312E+14
       coef(  291) =    -0.11784457573043064453E+14
       coef(  292) =    -0.49000490797772763672E+13
       coef(  293) =    -0.16820610042668609375E+14
       coef(  294) =    -0.48196321344279121094E+13
       coef(  295) =    -0.69347034967794656250E+14
       coef(  296) =    -0.26898010766290394531E+14
       coef(  297) =    -0.81831799179216054688E+13
       coef(  298) =    -0.94143655953861835938E+13
       coef(  299) =     0.62269254466469359375E+14
       coef(  300) =    -0.10781593674694376562E+15
       coef(  301) =    -0.31198152458698820312E+14
       coef(  302) =    -0.12230953340683439453E+14
       coef(  303) =    -0.70816442304709580078E+13
       coef(  304) =    -0.88177486520849625000E+14
       coef(  305) =    -0.26469626754222167969E+14
       coef(  306) =    -0.73653633352977050781E+13
       coef(  307) =    -0.84424310424961015625E+13
       coef(  308) =    -0.38976271978332695312E+14
       coef(  309) =    -0.10492718526048080078E+14
       coef(  310) =    -0.12646204418129435547E+14
       coef(  311) =    -0.26651187431818325422E+04
       coef(  312) =     0.31589293388367057778E+07
       coef(  313) =    -0.31694702680603575706E+09
       coef(  314) =     0.14083430290664489746E+11
       coef(  315) =    -0.22242099452061932373E+12
       coef(  316) =     0.14306579816473903809E+13
       coef(  317) =    -0.50298538300749472656E+13
       coef(  318) =     0.10905171737821939453E+14
       coef(  319) =    -0.92602432682114160156E+13
       coef(  320) =     0.14454197019401288033E+09
       coef(  321) =    -0.19393162045808391571E+11
       coef(  322) =     0.57459083121687866211E+12
       coef(  323) =    -0.30164594295612592773E+13
       coef(  324) =     0.10235902074367550781E+14
       coef(  325) =    -0.26986042097899593750E+14
       coef(  326) =    -0.37645642593212773438E+13
       coef(  327) =     0.44636568860999898438E+14
       coef(  328) =    -0.72099673221685287476E+11
       coef(  329) =    -0.61890344312850540161E+11
       coef(  330) =    -0.14214166583511748047E+14
       coef(  331) =     0.47168232409690593750E+14
       coef(  332) =     0.16126878313942773438E+12
       coef(  333) =    -0.77861308031274437500E+14
       coef(  334) =    -0.41156199673596882812E+14
       coef(  335) =     0.50477270317233437500E+13
       coef(  336) =     0.98231980472294667969E+13
       coef(  337) =    -0.92080173550719375000E+13
       coef(  338) =     0.38713838242628671875E+14
       coef(  339) =     0.17017376596740017578E+14
       coef(  340) =     0.75961227941287490234E+13
       coef(  341) =    -0.34269444937019957031E+14
       coef(  342) =    -0.27797869532715808594E+14
       coef(  343) =     0.93193222929655898438E+13
       coef(  344) =     0.10111320935936197266E+14
       coef(  345) =    -0.18527324535008214844E+14
       coef(  346) =    -0.16627165694653774414E+13
       coef(  347) =     0.72647586696214218140E+10
       coef(  348) =    -0.43710553215229534912E+12
       coef(  349) =    -0.12566185074696416016E+13
       coef(  350) =     0.81646328422611601562E+13
       coef(  351) =     0.99420487443119687500E+13
       coef(  352) =     0.42179190613807656250E+13
       coef(  353) =     0.14272669917283978516E+14
       coef(  354) =     0.19418002135653578125E+14
       coef(  355) =     0.33070314019257617188E+13
       coef(  356) =     0.15730625678320046875E+14
       coef(  357) =    -0.60770451669025765625E+14
       coef(  358) =    -0.52384182875495546875E+14
       coef(  359) =    -0.50407131080814195312E+14
       coef(  360) =    -0.32564437634626085938E+14
       coef(  361) =    -0.13050987557855359375E+14
       coef(  362) =    -0.62684160900142390625E+14
       coef(  363) =     0.89956473675261437500E+14
       coef(  364) =     0.64179639550413593750E+14
       coef(  365) =     0.23374053346115503906E+14
       coef(  366) =     0.53542114153880253906E+13
       coef(  367) =     0.77844882832000062500E+14
       coef(  368) =     0.35733629850394031250E+14
       coef(  369) =     0.13150802850546574219E+14
       coef(  370) =     0.12712332484904910156E+14
       coef(  371) =    -0.19081909117849031250E+14
       coef(  372) =     0.82475009385671406250E+14
       coef(  373) =     0.66428000750693437500E+14
       coef(  374) =     0.56313244721089921875E+13
       coef(  375) =    -0.15670072943675773438E+14
       coef(  376) =    -0.11956035721514699219E+14
       coef(  377) =    -0.22221151164387867188E+14
       coef(  378) =     0.11162271537802169922E+14
       coef(  379) =     0.67919637671832578125E+13
       coef(  380) =     0.41490149148295330811E+12
       coef(  381) =     0.71162234729992656250E+13
       coef(  382) =     0.34753936771921611328E+13
       coef(  383) =    -0.61770456471607750000E+14
       coef(  384) =    -0.17462248594091304688E+14
       coef(  385) =    -0.33384445142718354492E+13
       coef(  386) =    -0.59107874625654531250E+13
       coef(  387) =     0.26287152538230975342E+12
       coef(  388) =    -0.35156753168960390625E+13
       coef(  389) =     0.17636781668439878906E+14
       coef(  390) =    -0.45579102826322828125E+14
       coef(  391) =     0.14983137273954078125E+14
       coef(  392) =    -0.98233972238646523438E+13
       coef(  393) =    -0.68183532295283447266E+13
       coef(  394) =     0.89363347705272790527E+12
       coef(  395) =     0.14318155410643363281E+14
       coef(  396) =    -0.12679640042631815625E+15
       coef(  397) =     0.69828209563522218750E+14
       coef(  398) =     0.58872497714900765625E+14
       coef(  399) =     0.14648324058405023438E+14
       coef(  400) =     0.75897676699861572266E+12
       coef(  401) =     0.84188515954918578125E+14
       coef(  402) =     0.84940830329860953125E+14
       coef(  403) =     0.39327060469539312500E+14
       coef(  404) =     0.12071796443359783203E+14
       coef(  405) =     0.38714135294051781250E+14
       coef(  406) =     0.13676052915150910156E+14
       coef(  407) =     0.16101404509643667969E+14
       coef(  408) =     0.46242343493916968750E+14
       coef(  409) =     0.53554828708518671875E+14
       coef(  410) =     0.24171660360892066406E+14
       coef(  411) =     0.56247043107669101562E+13
       coef(  412) =     0.22559917714594210938E+14
       coef(  413) =     0.14209814975479171875E+14
       coef(  414) =     0.59687163757741933594E+13
       coef(  415) =     0.43647058685561601562E+13
       coef(  416) =    -0.23169596102813486328E+13
       coef(  417) =    -0.40232872508633874512E+12
       coef(  418) =    -0.27196610682888148438E+14
       coef(  419) =    -0.10099205834982669922E+14
       coef(  420) =     0.68276207335456279297E+13
       coef(  421) =     0.42884815742622744141E+13
       coef(  422) =    -0.82614382554654804688E+13
       coef(  423) =    -0.81387427768913708496E+12
       coef(  424) =    -0.48473222840516269531E+13
       coef(  425) =    -0.24075751779994186945E+07
       coef(  426) =    -0.21907071674724135548E+08
       coef(  427) =     0.64768168758226766586E+10
       coef(  428) =    -0.34492532882179040527E+12
       coef(  429) =     0.36867497092551318359E+13
       coef(  430) =    -0.16279773822425634766E+14
       coef(  431) =     0.42330364645294773438E+14
       coef(  432) =    -0.71448064833377437500E+14
       coef(  433) =     0.40866116294390859375E+14
       coef(  434) =     0.10510148057384524345E+10
       coef(  435) =     0.22657382990971893311E+12
       coef(  436) =    -0.51817797428183876953E+13
       coef(  437) =     0.25329206960278562500E+14
       coef(  438) =    -0.64303931446624609375E+14
       coef(  439) =     0.84318888114076593750E+14
       coef(  440) =     0.10668326088684196875E+15
       coef(  441) =     0.86135410185577156250E+14
       coef(  442) =     0.12008847512508276367E+12
       coef(  443) =     0.10061390755258398438E+14
       coef(  444) =     0.47474280555878203125E+14
       coef(  445) =    -0.96951167135532781250E+14
       coef(  446) =    -0.92814341304299687500E+14
       coef(  447) =    -0.55607507947773750000E+14
       coef(  448) =    -0.21070614347765183594E+14
       coef(  449) =    -0.65930473488766609375E+14
       coef(  450) =     0.39606177996297804688E+14
       coef(  451) =    -0.33897417290858085938E+14
       coef(  452) =    -0.24075698861798429688E+14
       coef(  453) =    -0.10879362703173910156E+14
       coef(  454) =     0.37774693855753320312E+14
       coef(  455) =    -0.31293677876725166016E+13
       coef(  456) =    -0.33418397660719521484E+13
       coef(  457) =    -0.44797854508260380859E+13
       coef(  458) =    -0.21076011235671185303E+12
       coef(  459) =     0.83774577746402060547E+13
       coef(  460) =    -0.23640593733024339844E+14
       coef(  461) =    -0.18794858611528914062E+14
       coef(  462) =    -0.11979631895063085938E+14
       coef(  463) =    -0.25203694469315867188E+14
       coef(  464) =    -0.10123026731869386719E+14
       coef(  465) =     0.22870104065342578125E+13
       coef(  466) =    -0.25295972547938820312E+14
       coef(  467) =     0.29460431859163863281E+14
       coef(  468) =     0.60844084265436750000E+14
       coef(  469) =    -0.65091893055950722656E+13
       coef(  470) =    -0.25886523061207140625E+14
       coef(  471) =    -0.16045169101230519531E+14
       coef(  472) =     0.88100805041912296875E+14
       coef(  473) =     0.92040437704356015625E+14
       coef(  474) =     0.28500598798219058594E+14
       coef(  475) =     0.37627269223912788086E+13
       coef(  476) =     0.50476353475937617188E+14
       coef(  477) =     0.15192940356449279297E+14
       coef(  478) =     0.21230887513244042969E+14
       coef(  479) =    -0.19651887237774171875E+15
       coef(  480) =    -0.33592252597779796875E+14
       coef(  481) =    -0.10084123041788123047E+14
       coef(  482) =    -0.82084491989246289062E+13
       coef(  483) =    -0.10960907242875243750E+15
       coef(  484) =    -0.15536997827220675781E+14
       coef(  485) =    -0.13336127890786186523E+13
       coef(  486) =    -0.85344669664483886719E+11
       coef(  487) =    -0.55131677053943601562E+14
       coef(  488) =    -0.11352405094172023438E+14
       coef(  489) =    -0.17215923619849697266E+13
       coef(  490) =    -0.10094494964767445312E+14
       coef(  491) =     0.10462211297923582812E+15
       coef(  492) =     0.24459474099799074219E+14
       coef(  493) =    -0.47193452319446210938E+13
       coef(  494) =    -0.18322852358816226562E+14
       coef(  495) =    -0.10764352093909203125E+14
       coef(  496) =     0.45018702209467921875E+14
       coef(  497) =    -0.18758844411821268750E+15
       coef(  498) =    -0.30120253289243058594E+14
       coef(  499) =    -0.10844211296514160156E+13
       coef(  500) =    -0.25847187441464824219E+13
       coef(  501) =    -0.28732943245350328125E+14
       coef(  502) =     0.14384439081934062500E+14
       coef(  503) =     0.78318515746502851562E+13
       coef(  504) =     0.10514354679094890625E+14
       coef(  505) =    -0.21330932294076914062E+14
       coef(  506) =    -0.46903074711441039062E+14
       coef(  507) =    -0.44185894541316220703E+13
       coef(  508) =     0.13895502446991813965E+13
       coef(  509) =    -0.10312071376481296875E+14
       coef(  510) =     0.76601671551390820312E+12
       coef(  511) =    -0.35352248504672568359E+13
       coef(  512) =    -0.21066666451651710938E+14
       coef(  513) =    -0.11856408281432173828E+14
       coef(  514) =    -0.10274181203754060059E+13
       coef(  515) =    -0.31886561157000878906E+13
       coef(  516) =     0.77197622410656765103E+08
       coef(  517) =    -0.29006368106894860268E+10
       coef(  518) =     0.41956126186586471558E+11
       coef(  519) =     0.17345156100104741211E+13
       coef(  520) =    -0.14317337671003189453E+14
       coef(  521) =     0.21396067550834453125E+14
       coef(  522) =     0.24142871852881500000E+14
       coef(  523) =    -0.50767329880519156250E+14
       coef(  524) =     0.19361530808772609375E+14
       coef(  525) =    -0.17324070331266475677E+11
       coef(  526) =    -0.10234471198046511230E+13
       coef(  527) =     0.31861003976523188477E+13
       coef(  528) =    -0.25619741577991289062E+13
       coef(  529) =     0.47868734444322675781E+13
       coef(  530) =    -0.37602309099844515625E+14
       coef(  531) =    -0.20396056763817125000E+14
       coef(  532) =     0.31199064632590380859E+13
       coef(  533) =    -0.32492265361668984375E+13
       coef(  534) =     0.23579552439079042969E+14
       coef(  535) =    -0.29575870756544339844E+14
       coef(  536) =    -0.35436635046269257812E+14
       coef(  537) =    -0.47036124297469171875E+14
       coef(  538) =    -0.29195136221502914062E+14
       coef(  539) =     0.32367767561745796875E+14
       coef(  540) =     0.56107592247006390625E+14
       coef(  541) =     0.72501250467966845703E+13
       coef(  542) =    -0.75847729170483583984E+13
       coef(  543) =     0.42521201876560703125E+14
       coef(  544) =     0.99237333657183300781E+13
       coef(  545) =     0.47381512568172943115E+12
       coef(  546) =    -0.14678138804673333984E+14
       coef(  547) =     0.38930128372014664062E+14
       coef(  548) =    -0.16045896268045880859E+14
       coef(  549) =     0.21449641544115812500E+14
       coef(  550) =    -0.62347647194351171875E+13
       coef(  551) =    -0.81407534658061328125E+13
       coef(  552) =     0.56697841805165859375E+14
       coef(  553) =    -0.20479409060131132812E+13
       coef(  554) =     0.41712719508200046875E+14
       coef(  555) =     0.13727345577004656250E+14
       coef(  556) =    -0.38290929021599228516E+13
       coef(  557) =     0.53624273430761984375E+14
       coef(  558) =     0.43769965989503031250E+14
       coef(  559) =     0.13499109899403009766E+14
       coef(  560) =     0.21020585008157671875E+14
       coef(  561) =     0.79116945811925968750E+14
       coef(  562) =    -0.69744070347862812500E+14
       coef(  563) =    -0.10201305106102968750E+14
       coef(  564) =    -0.11210658562849345703E+13
       coef(  565) =    -0.38828595754304640625E+14
       coef(  566) =    -0.53571571968992851562E+13
       coef(  567) =    -0.18812526638646097656E+14
       coef(  568) =     0.74190412157281103516E+13
       coef(  569) =    -0.22445035960909476562E+14
       coef(  570) =    -0.11400485922384142188E+15
       coef(  571) =    -0.50974887540151347656E+13
       coef(  572) =     0.58383586439364619141E+13
       coef(  573) =    -0.33801455395832763672E+13
       coef(  574) =     0.80561560721778554688E+13
       coef(  575) =    -0.10391961323625450000E+15
       coef(  576) =    -0.20990226387287304688E+14
       coef(  577) =    -0.24199051791873139648E+13
       coef(  578) =    -0.25631364290984382812E+14
       coef(  579) =     0.25060327233735253906E+12
       coef(  580) =     0.71303670175553554688E+13
       coef(  581) =    -0.24084048331757375000E+14
       coef(  582) =    -0.53363344196307910156E+13
       coef(  583) =    -0.67039133511195175781E+13
       coef(  584) =     0.26076202166597656250E+11
       coef(  585) =    -0.45244463790586005859E+13
       coef(  586) =    -0.94590442196306967735E+09
       coef(  587) =     0.43994502949880836487E+11
       coef(  588) =    -0.93684361095715966797E+12
       coef(  589) =    -0.30657072062868710938E+13
       coef(  590) =     0.43356216192607421875E+14
       coef(  591) =    -0.56378606124913617188E+14
       coef(  592) =    -0.85797777598192968750E+12
       coef(  593) =    -0.33372982298735199219E+14
       coef(  594) =    -0.11639501979783957031E+14
       coef(  595) =     0.10974301983216909790E+12
       coef(  596) =     0.10943709182899953125E+14
       coef(  597) =     0.82497244663451113281E+13
       coef(  598) =    -0.58953229317341531250E+14
       coef(  599) =     0.53892375735549031250E+14
       coef(  600) =    -0.55210061191935673828E+13
       coef(  601) =    -0.16154981503356722656E+14
       coef(  602) =    -0.27502619927426148438E+14
       coef(  603) =    -0.33297358971581101562E+14
       coef(  604) =    -0.55783131430725078125E+14
       coef(  605) =    -0.86699107638901406250E+13
       coef(  606) =    -0.12566175893448232422E+14
       coef(  607) =    -0.16030963557877566406E+14
       coef(  608) =     0.34243467331887343750E+13
       coef(  609) =     0.24921093377453134766E+13
       coef(  610) =     0.97478635755955839844E+13
       coef(  611) =    -0.46329555643206542969E+13
       coef(  612) =     0.40830416881277988281E+13
       coef(  613) =    -0.25589922900829750000E+14
       coef(  614) =     0.25387895848460109375E+14
       coef(  615) =     0.39554278733017484375E+14
       coef(  616) =     0.92353220444431074219E+13
       coef(  617) =    -0.14009207769399543750E+15
       coef(  618) =    -0.68550409294083964844E+13
       coef(  619) =     0.26617950353913757812E+14
       coef(  620) =     0.13491289472619527344E+14
       coef(  621) =     0.25205193405320910156E+14
       coef(  622) =     0.17214355625294189453E+14
       coef(  623) =     0.42324677342199039062E+14
       coef(  624) =    -0.38655872560532148438E+13
       coef(  625) =     0.34303787529159990234E+13
       coef(  626) =    -0.54977489432117216797E+13
       coef(  627) =     0.39298440975029804688E+13
       coef(  628) =     0.30713492227234601562E+14
       coef(  629) =    -0.30449333442639675781E+14
       coef(  630) =     0.19846732531044558594E+14
       coef(  631) =     0.13594236750500277344E+14
       coef(  632) =     0.28047199547258289062E+14
       coef(  633) =    -0.26415045334803296875E+14
       coef(  634) =    -0.31273579205373085938E+13
       coef(  635) =    -0.75663574140066738281E+13
       coef(  636) =     0.21543063901359957031E+14
       coef(  637) =    -0.40356984286324633789E+13
       coef(  638) =     0.65725347629988398438E+13
       coef(  639) =     0.50855359752138500214E+10
       coef(  640) =    -0.24150958283525100708E+12
       coef(  641) =     0.54182098546515322266E+13
       coef(  642) =    -0.14480585752913791016E+14
       coef(  643) =    -0.30408751742001257812E+14
       coef(  644) =    -0.84760390713590634766E+13
       coef(  645) =     0.38569867068280765625E+14
       coef(  646) =     0.10270730975200601562E+14
       coef(  647) =    -0.82787774587954296875E+12
       coef(  648) =    -0.34294142704565078125E+14
       coef(  649) =    -0.11962338393369304688E+14
       coef(  650) =     0.35111351719038589844E+14
       coef(  651) =     0.64596131755142726562E+14
       coef(  652) =     0.12227809681399878906E+14
       coef(  653) =     0.10400225207358242188E+15
       coef(  654) =     0.23474532563161437500E+14
       coef(  655) =    -0.33202464313071746094E+14
       coef(  656) =    -0.49628048747782314453E+13
       coef(  657) =    -0.14812878044527386719E+14
       coef(  658) =    -0.11668136565601134766E+14
       coef(  659) =     0.16489656448354056641E+14
       coef(  660) =     0.75622475195368859375E+14
       coef(  661) =     0.12510864955068984375E+13
       coef(  662) =     0.30124003750510402344E+14
       coef(  663) =     0.25732389553681562500E+14
       coef(  664) =    -0.43441191940259734375E+14
       coef(  665) =     0.11550295702689937500E+14
       coef(  666) =     0.14677700984530898438E+14
       coef(  667) =     0.13187679389429769531E+14
       coef(  668) =     0.34079508890035921875E+14
       coef(  669) =     0.95680622896863359375E+13
       coef(  670) =    -0.43504053769739273438E+14
       coef(  671) =     0.28750112852797449219E+14
       coef(  672) =     0.30885494576877382812E+13
       coef(  673) =     0.16528217978936960938E+14
       coef(  674) =     0.23003077288973777344E+14
       coef(  675) =    -0.24308812770348232422E+13
       coef(  676) =     0.14581788819015378906E+14
       coef(  677) =    -0.11395031248568449020E+11
       coef(  678) =     0.49423739060109167480E+12
       coef(  679) =    -0.12176725107760353516E+14
       coef(  680) =     0.54564995280402781250E+14
       coef(  681) =    -0.27204818347053414062E+14
       coef(  682) =     0.34611052618214710938E+14
       coef(  683) =     0.54342079895068429688E+14
       coef(  684) =     0.47368753405057978516E+13
       coef(  685) =     0.83853762219000878906E+13
       coef(  686) =     0.78715682459996875000E+13
       coef(  687) =     0.41988268973554867188E+14
       coef(  688) =     0.40235614793455195312E+14
       coef(  689) =    -0.10197987144087341797E+14
       coef(  690) =    -0.21647898219991550781E+14
       coef(  691) =    -0.30832866090245437500E+14
       coef(  692) =    -0.23980012422279351562E+14
       coef(  693) =    -0.33100210527747085938E+14
       coef(  694) =    -0.88480603862958320312E+13
       coef(  695) =    -0.33954759988644601562E+14
       coef(  696) =     0.53716749366826064453E+13
       coef(  697) =    -0.27497471033054070312E+14
       coef(  698) =     0.13774025615496977539E+13
       coef(  699) =     0.13702562745251464844E+14
       coef(  700) =     0.86038675172731933594E+13
       coef(  701) =    -0.29087943945503710938E+13
       coef(  702) =    -0.34263229273873320312E+13
       coef(  703) =     0.64580554514737666016E+13
       coef(  704) =     0.96769432734124832153E+10
       coef(  705) =    -0.38216890660750268555E+12
       coef(  706) =     0.99756009154737148438E+13
       coef(  707) =    -0.48782588011825554688E+14
       coef(  708) =     0.78539820736025781250E+12
       coef(  709) =     0.69985239195413164062E+14
       coef(  710) =    -0.60406542740613447266E+13
       coef(  711) =     0.40568473116059460938E+14
       coef(  712) =    -0.45020949255149296875E+13
       coef(  713) =     0.24677709551794320312E+14
       coef(  714) =    -0.14242844511321668750E+15
       coef(  715) =    -0.57148730878089000000E+14
       coef(  716) =     0.19408715602746914062E+14
       coef(  717) =    -0.72663311685540968750E+14
       coef(  718) =    -0.48741208275602812500E+14
       coef(  719) =    -0.33261801270246125000E+14
       coef(  720) =     0.88431124291308484375E+14
       coef(  721) =    -0.75660207507392265625E+13
       coef(  722) =    -0.14287513463115857348E+03
       coef(  723) =    -0.21557858923307350778E+05
       coef(  724) =     0.36922485745783876628E+07
       coef(  725) =    -0.75708674652054026723E+08
       coef(  726) =     0.20908120484373142719E+10
       coef(  727) =    -0.11862163570234518051E+10
       coef(  728) =    -0.13917982027503349304E+11
       coef(  729) =    -0.17520175116250335693E+12
       coef(  730) =     0.29374042817088757324E+12
       coef(  731) =    -0.45919368192603980424E+06
       coef(  732) =    -0.39593850713104450703E+09
       coef(  733) =    -0.39606279912058629990E+10
       coef(  734) =    -0.14111184168828549194E+12
       coef(  735) =     0.70926927301657739258E+12
       coef(  736) =     0.13060937147124322510E+12
       coef(  737) =     0.17097882150776865234E+13
       coef(  738) =    -0.60518585292481689453E+13
       coef(  739) =     0.10364152281726320267E+11
       coef(  740) =    -0.85943356837164596558E+11
       coef(  741) =     0.20008598803300837402E+13
       coef(  742) =    -0.15014299724105074219E+14
       coef(  743) =     0.47013836562802210938E+14
       coef(  744) =    -0.89678576086767765625E+14
       coef(  745) =     0.98716562903084531250E+14
       coef(  746) =    -0.32282918083258666992E+10
       coef(  747) =    -0.10940818287853645020E+13
       coef(  748) =     0.17956821783159375000E+13
       coef(  749) =    -0.54171458326785781250E+13
       coef(  750) =    -0.46371279575107742188E+14
       coef(  751) =     0.13853839550434101562E+14
       coef(  752) =     0.29633489100118740234E+13
       coef(  753) =     0.63840554394055068359E+13
       coef(  754) =     0.87673959201346445312E+13
       coef(  755) =    -0.94003011956025781250E+12
       coef(  756) =     0.68156268572726328125E+13
       coef(  757) =    -0.46830686637181191406E+13
       coef(  758) =    -0.38863667043084185791E+12
       coef(  759) =     0.11412308377593310547E+13
       coef(  760) =    -0.14332313745542172852E+13
       coef(  761) =     0.33500599866767436266E+09
       coef(  762) =     0.78286991238266525269E+10
       coef(  763) =     0.42710435392759155273E+12
       coef(  764) =     0.13132975169702209473E+12
       coef(  765) =    -0.10878173889138523438E+14
       coef(  766) =     0.66891481205823935547E+13
       coef(  767) =     0.34537030250872062500E+14
       coef(  768) =    -0.35311823309639242188E+14
       coef(  769) =    -0.53241368356578417969E+12
       coef(  770) =    -0.45583768723449423828E+13
       coef(  771) =     0.25256224452234945312E+14
       coef(  772) =     0.29184912504036601562E+13
       coef(  773) =    -0.31243571178809464844E+14
       coef(  774) =    -0.52711889992694218750E+14
       coef(  775) =    -0.31937604266811925781E+14
       coef(  776) =     0.21053954408320250000E+14
       coef(  777) =    -0.10157063260096621875E+15
       coef(  778) =     0.12645361979858575000E+15
       coef(  779) =     0.88759983708532125000E+14
       coef(  780) =     0.26711984209734851562E+14
       coef(  781) =     0.61654542394675087891E+13
       coef(  782) =    -0.32496167087727367188E+14
       coef(  783) =     0.52949768458794179688E+14
       coef(  784) =     0.38685814863408671875E+14
       coef(  785) =     0.15198312988998275391E+14
       coef(  786) =     0.26303287555140054688E+14
       coef(  787) =     0.12019336658976968750E+14
       coef(  788) =     0.63351611759934570312E+13
       coef(  789) =    -0.52277825165083937500E+14
       coef(  790) =     0.24832810934686046875E+14
       coef(  791) =    -0.25385405079239601562E+14
       coef(  792) =    -0.45295666025883992188E+14
       coef(  793) =    -0.34978998427360621094E+14
       coef(  794) =    -0.17503849578009978516E+14
       coef(  795) =     0.12538044660752840625E+15
       coef(  796) =    -0.33569801440611562500E+14
       coef(  797) =    -0.17882854909689390625E+14
       coef(  798) =    -0.82020013688795585938E+13
       coef(  799) =    -0.51043474502787109375E+13
       coef(  800) =    -0.41213538822897062500E+14
       coef(  801) =    -0.11317296846554541016E+14
       coef(  802) =    -0.16247725551734064941E+13
       coef(  803) =    -0.28983244201870131836E+13
       coef(  804) =    -0.64127910038554101562E+14
       coef(  805) =    -0.72554646669522109375E+14
       coef(  806) =    -0.27473352634450058594E+14
       coef(  807) =    -0.83822149627654873047E+13
       coef(  808) =    -0.36505492154434382812E+14
       coef(  809) =    -0.11090529056008703125E+14
       coef(  810) =    -0.14057028055610181641E+14
       coef(  811) =    -0.11490185361564630508E+11
       coef(  812) =     0.12989783457767950439E+12
       coef(  813) =    -0.35562211837347944336E+13
       coef(  814) =     0.12104188929104851562E+14
       coef(  815) =     0.24009379852072656250E+12
       coef(  816) =    -0.22787061102514789062E+14
       coef(  817) =     0.53474245312686523438E+13
       coef(  818) =     0.10810624260872691406E+14
       coef(  819) =     0.32592092383494531250E+13
       coef(  820) =     0.74727018371354257812E+13
       coef(  821) =    -0.37754511890622421875E+14
       coef(  822) =     0.24131058156336812500E+14
       coef(  823) =    -0.75277966955016015625E+13
       coef(  824) =    -0.11268409694109279297E+14
       coef(  825) =    -0.44147550519778955078E+13
       coef(  826) =    -0.26461029764746242188E+14
       coef(  827) =     0.46330791204015898438E+13
       coef(  828) =     0.40067635684781093750E+14
       coef(  829) =     0.19568992009757605469E+14
       coef(  830) =     0.57990008716526171875E+13
       coef(  831) =     0.83776288933017578125E+13
       coef(  832) =     0.13925018726201429688E+14
       coef(  833) =     0.76705565432995468750E+13
       coef(  834) =     0.52247893989738378906E+13
       coef(  835) =    -0.29237422494628148438E+14
       coef(  836) =     0.78122497977792484375E+14
       coef(  837) =     0.27887045898183777344E+14
       coef(  838) =    -0.41163219651600415039E+13
       coef(  839) =    -0.14644217869070976562E+14
       coef(  840) =    -0.93281355650938203125E+13
       coef(  841) =    -0.64358921023881156250E+14
       coef(  842) =    -0.34433756383372914062E+14
       coef(  843) =    -0.11784457573049728516E+14
       coef(  844) =    -0.49000490797859355469E+13
       coef(  845) =    -0.16820610042659929688E+14
       coef(  846) =    -0.48196321344328544922E+13
       coef(  847) =    -0.69347034967800546875E+14
       coef(  848) =    -0.26898010766300531250E+14
       coef(  849) =    -0.81831799179263671875E+13
       coef(  850) =    -0.94143655953947421875E+13
       coef(  851) =     0.62269254466468328125E+14
       coef(  852) =    -0.10781593674695164062E+15
       coef(  853) =    -0.31198152458696117188E+14
       coef(  854) =    -0.12230953340681058594E+14
       coef(  855) =    -0.70816442304642841797E+13
       coef(  856) =    -0.88177486520853640625E+14
       coef(  857) =    -0.26469626754218730469E+14
       coef(  858) =    -0.73653633352926376953E+13
       coef(  859) =    -0.84424310425009814453E+13
       coef(  860) =    -0.38976271978333390625E+14
       coef(  861) =    -0.10492718526053197266E+14
       coef(  862) =    -0.12646204418132171875E+14
       coef(  863) =     0.49670116825495088051E+04
       coef(  864) =    -0.67329113896426539868E+07
       coef(  865) =     0.64801298335973381996E+09
       coef(  866) =    -0.30713089667266914368E+11
       coef(  867) =     0.48348485922768640137E+12
       coef(  868) =    -0.33544610573637729492E+13
       coef(  869) =     0.10350872381841480469E+14
       coef(  870) =    -0.14631868695806746094E+14
       coef(  871) =     0.11822107224336769531E+14
       coef(  872) =    -0.32574027174298804998E+09
       coef(  873) =     0.58686036525235420227E+11
       coef(  874) =    -0.10150387766450244141E+13
       coef(  875) =     0.37093069329619130859E+13
       coef(  876) =     0.10643344603698347656E+14
       coef(  877) =    -0.73772535555640906250E+14
       coef(  878) =     0.66896871679406546875E+14
       coef(  879) =    -0.60820955178888328125E+14
       coef(  880) =    -0.71118334089906768799E+11
       coef(  881) =     0.33498684681562871094E+13
       coef(  882) =     0.14931717061111253906E+14
       coef(  883) =    -0.74441394406241250000E+14
       coef(  884) =     0.17011160911036906250E+15
       coef(  885) =     0.85450214660421031250E+14
       coef(  886) =     0.12618868766782000000E+14
       coef(  887) =    -0.32452989188461750000E+14
       coef(  888) =     0.86457613798431500000E+14
       coef(  889) =    -0.16798018156173640625E+15
       coef(  890) =    -0.44934063857185445312E+14
       coef(  891) =    -0.84261590339967333984E+13
       coef(  892) =    -0.24858990089028041992E+13
       coef(  893) =     0.94962871856928750000E+14
       coef(  894) =    -0.46764745943405828125E+14
       coef(  895) =    -0.22392903454447699219E+14
       coef(  896) =    -0.59922157012964619141E+13
       coef(  897) =    -0.34163409897066941406E+14
       coef(  898) =    -0.12817721567528095703E+14
       coef(  899) =    -0.39062070171951721191E+11
       coef(  900) =     0.10149213287636715088E+13
       coef(  901) =    -0.14949155532944882812E+13
       coef(  902) =    -0.29743165904423527344E+14
       coef(  903) =     0.99680215191563031250E+14
       coef(  904) =     0.45892616407141156250E+14
       coef(  905) =    -0.41334528395393750000E+14
       coef(  906) =    -0.70476514471552984375E+14
       coef(  907) =    -0.89595275505965371094E+13
       coef(  908) =     0.70144668575934328125E+14
       coef(  909) =    -0.15023878206813553125E+15
       coef(  910) =    -0.17199396757769876953E+14
       coef(  911) =    -0.15857997621783994141E+14
       coef(  912) =    -0.29502531454445007812E+14
       coef(  913) =    -0.23188695431468160156E+14
       coef(  914) =    -0.48592640841146351562E+14
       coef(  915) =     0.11248785037658595703E+14
       coef(  916) =     0.23895808818789800781E+14
       coef(  917) =     0.85220272919480839844E+13
       coef(  918) =    -0.13718992467395456543E+13
       coef(  919) =     0.19165047247477156250E+14
       coef(  920) =     0.90610639340012851562E+13
       coef(  921) =     0.32896131595524780273E+13
       coef(  922) =     0.19863834145278972168E+13
       coef(  923) =     0.73865003268427539062E+13
       coef(  924) =    -0.22973337695914968750E+14
       coef(  925) =     0.41066549595554117188E+14
       coef(  926) =     0.22278339885330164062E+14
       coef(  927) =    -0.46449151762373369141E+13
       coef(  928) =    -0.10073033488187609375E+14
       coef(  929) =    -0.53877790416380671875E+14
       coef(  930) =    -0.43062363929605322266E+13
       coef(  931) =     0.28301902427021508789E+13
       coef(  932) =    -0.50757626730851641846E+12
       coef(  933) =    -0.31853022779937436523E+13
       coef(  934) =    -0.57678775487573339844E+12
       coef(  935) =    -0.59422874662846390625E+14
       coef(  936) =    -0.17758867024942480469E+14
       coef(  937) =    -0.37850627842411005859E+13
       coef(  938) =    -0.64682354304850410156E+13
       coef(  939) =     0.24691143352532379150E+11
       coef(  940) =    -0.14695926641910090332E+13
       coef(  941) =     0.13636598684921798828E+14
       coef(  942) =    -0.77141696929325718750E+14
       coef(  943) =     0.28522588638146675781E+14
       coef(  944) =    -0.20772001046720921875E+14
       coef(  945) =    -0.35695963336291632812E+14
       coef(  946) =    -0.24784674272532878906E+14
       coef(  947) =     0.52377812125198789062E+13
       coef(  948) =     0.96414460529683046875E+14
       coef(  949) =     0.62948308675998328125E+14
       coef(  950) =     0.34960317683221199219E+14
       coef(  951) =     0.19309559806060125732E+12
       coef(  952) =    -0.81235145126318964844E+13
       coef(  953) =     0.14059651684129823438E+15
       coef(  954) =     0.61221555825721257812E+14
       coef(  955) =     0.22386348580951398438E+14
       coef(  956) =     0.50151793631770546875E+13
       coef(  957) =     0.20563158593683222656E+14
       coef(  958) =     0.61387211468519306641E+13
       coef(  959) =    -0.84279559403959078125E+14
       coef(  960) =     0.32489002554039742188E+14
       coef(  961) =     0.25852768281862878906E+14
       coef(  962) =     0.96097806867591152344E+13
       coef(  963) =    -0.24989861406548950195E+12
       coef(  964) =     0.88684771225875039062E+13
       coef(  965) =     0.36900758147899067383E+13
       coef(  966) =     0.14435206162287717285E+13
       coef(  967) =     0.21223354450170745850E+12
       coef(  968) =    -0.76740286658225996094E+13
       coef(  969) =    -0.28272735068458828125E+13
       coef(  970) =    -0.89315773983969406250E+14
       coef(  971) =    -0.28109839673757742188E+14
       coef(  972) =    -0.41596848681481372070E+13
       coef(  973) =    -0.42577415629192724609E+12
       coef(  974) =    -0.14126573807906605469E+14
       coef(  975) =    -0.35262181464740883789E+13
       coef(  976) =    -0.62619630436775312500E+13
       coef(  977) =     0.16457306050914905500E+07
       coef(  978) =    -0.12003693658087346703E+08
       coef(  979) =    -0.41884498239852671623E+10
       coef(  980) =     0.18203767082257040405E+12
       coef(  981) =    -0.17029988748056091309E+13
       coef(  982) =     0.80678935523036083984E+13
       coef(  983) =    -0.21176751633663050781E+14
       coef(  984) =     0.44914285486881500000E+14
       coef(  985) =    -0.51970964286422468750E+14
       coef(  986) =     0.32942216872651052475E+10
       coef(  987) =    -0.37605500079439227295E+12
       coef(  988) =     0.37259383845181469727E+13
       coef(  989) =    -0.16998612507288214844E+14
       coef(  990) =     0.82024611262562500000E+13
       coef(  991) =     0.21696205004859375000E+14
       coef(  992) =     0.44399341724556312500E+14
       coef(  993) =     0.12327368489921109375E+14
       coef(  994) =     0.72452802756998315430E+12
       coef(  995) =    -0.93808508041049003906E+13
       coef(  996) =    -0.17440988272368517578E+14
       coef(  997) =     0.22231932081742863281E+14
       coef(  998) =     0.10600633038189109375E+14
       coef(  999) =    -0.10743860636907400391E+14
       coef( 1000) =    -0.13132454438848679688E+14
       coef( 1001) =     0.66423319001862015625E+14
       coef( 1002) =    -0.50001288072596851562E+14
       coef( 1003) =    -0.71284820015896250000E+14
       coef( 1004) =    -0.31894325201322074219E+14
       coef( 1005) =    -0.13292714902502349609E+14
       coef( 1006) =    -0.44323338237859867188E+14
       coef( 1007) =    -0.35264142839499851562E+14
       coef( 1008) =    -0.14074221610256755859E+14
       coef( 1009) =    -0.16075275883438597656E+14
       coef( 1010) =     0.14204638413768252563E+12
       coef( 1011) =    -0.53104465301811547852E+12
       coef( 1012) =    -0.29240974915894628906E+12
       coef( 1013) =     0.64439507415401062500E+14
       coef( 1014) =    -0.98009226775112437500E+14
       coef( 1015) =    -0.99682638439516750000E+14
       coef( 1016) =    -0.54537619659010523438E+14
       coef( 1017) =    -0.25303531722825382812E+14
       coef( 1018) =    -0.61018907855594921875E+12
       coef( 1019) =    -0.46748097063315406250E+14
       coef( 1020) =     0.11648654056683025000E+15
       coef( 1021) =     0.37398674879677382812E+14
       coef( 1022) =    -0.68951429647078886719E+13
       coef( 1023) =    -0.11604732944044054688E+14
       coef( 1024) =    -0.25971944691942914062E+14
       coef( 1025) =     0.42179720950745460938E+14
       coef( 1026) =     0.16443053678413472656E+14
       coef( 1027) =     0.16779340050174858398E+13
       coef( 1028) =     0.20246041502273703125E+14
       coef( 1029) =     0.54589861502835214844E+13
       coef( 1030) =     0.24245936928070023438E+14
       coef( 1031) =    -0.57309442282479789062E+14
       coef( 1032) =     0.44926795033278625000E+14
       coef( 1033) =     0.21179675345057492188E+14
       coef( 1034) =     0.23355142142905888672E+13
       coef( 1035) =    -0.51774036695847359375E+14
       coef( 1036) =    -0.21161497976152062988E+12
       coef( 1037) =     0.29851152328660444336E+13
       coef( 1038) =     0.95181107973740722656E+12
       coef( 1039) =    -0.31637390487198535156E+14
       coef( 1040) =    -0.59249489814623476562E+13
       coef( 1041) =    -0.10420893373334082031E+13
       coef( 1042) =     0.15452977129798376953E+14
       coef( 1043) =    -0.71371119550734328125E+14
       coef( 1044) =     0.10552409022168360938E+15
       coef( 1045) =     0.20788289408318484375E+14
       coef( 1046) =    -0.19641927485820304688E+14
       coef( 1047) =    -0.17015991425286238281E+14
       coef( 1048) =    -0.46754371912050492188E+14
       coef( 1049) =    -0.50802492311341078125E+14
       coef( 1050) =     0.27168282604401859375E+14
       coef( 1051) =     0.12590405491584548828E+14
       coef( 1052) =    -0.39229412275205419922E+12
       coef( 1053) =     0.16394201591365664062E+13
       coef( 1054) =     0.13354567902683296875E+14
       coef( 1055) =     0.55481332581869218750E+13
       coef( 1056) =     0.56469466730527871094E+13
       coef( 1057) =     0.11916499911703084375E+15
       coef( 1058) =     0.27324511409608398438E+14
       coef( 1059) =     0.13545680522376128906E+14
       coef( 1060) =     0.45370238818226142578E+13
       coef( 1061) =     0.64183363989175683594E+13
       coef( 1062) =     0.27201102467723144531E+13
       coef( 1063) =    -0.82854545397324218750E+10
       coef( 1064) =     0.66522970345149531250E+14
       coef( 1065) =     0.12604065698566126953E+14
       coef( 1066) =     0.31948777015791450195E+13
       coef( 1067) =     0.16893543840889951172E+13
       coef( 1068) =    -0.10529739028779822588E+09
       coef( 1069) =     0.33022205725102000237E+10
       coef( 1070) =     0.23978790889556648254E+11
       coef( 1071) =    -0.24420118778795454102E+13
       coef( 1072) =     0.11140327862216070312E+14
       coef( 1073) =     0.20778907053942226562E+13
       coef( 1074) =    -0.68782756298745523438E+14
       coef( 1075) =    -0.33146531400096832031E+14
       coef( 1076) =    -0.19810613601345278320E+13
       coef( 1077) =     0.88012361693432296753E+11
       coef( 1078) =     0.11798158807043283691E+13
       coef( 1079) =     0.52212544321915722656E+13
       coef( 1080) =     0.19333275469055859375E+14
       coef( 1081) =     0.72914674008626625000E+14
       coef( 1082) =     0.10307566831292687500E+14
       coef( 1083) =     0.24155090942416699219E+13
       coef( 1084) =     0.30242522797073100586E+13
       coef( 1085) =     0.10359883388191833984E+14
       coef( 1086) =    -0.15207075564985200000E+15
       coef( 1087) =     0.14418275759198025000E+15
       coef( 1088) =     0.12159070676859368750E+15
       coef( 1089) =     0.36895613371166093750E+14
       coef( 1090) =     0.43296208603838017578E+13
       coef( 1091) =     0.13770965989380691406E+14
       coef( 1092) =     0.39635811672388093750E+14
       coef( 1093) =     0.15958481450756265625E+14
       coef( 1094) =     0.25346643716712832031E+13
       coef( 1095) =     0.73790727856803740234E+13
       coef( 1096) =    -0.97541372140769702148E+12
       coef( 1097) =    -0.20946705256003999023E+13
       coef( 1098) =    -0.10353225338138753906E+14
       coef( 1099) =     0.67753392614249667969E+13
       coef( 1100) =    -0.23342697696863695312E+14
       coef( 1101) =    -0.51331766244959953125E+14
       coef( 1102) =    -0.41983586086475710938E+14
       coef( 1103) =    -0.20001132683152816406E+14
       coef( 1104) =     0.43413511107943421875E+14
       coef( 1105) =     0.28190524890334582031E+14
       coef( 1106) =     0.66596756120113367188E+14
       coef( 1107) =     0.21385429883811695312E+14
       coef( 1108) =     0.19363008414457275391E+12
       coef( 1109) =     0.31927821435305339844E+14
       coef( 1110) =     0.29478236870645324219E+14
       coef( 1111) =     0.96652632074654843750E+13
       coef( 1112) =     0.11668563194208353516E+14
       coef( 1113) =    -0.27973247147605203125E+14
       coef( 1114) =    -0.29470662796310445312E+14
       coef( 1115) =     0.12408502654223798828E+14
       coef( 1116) =     0.63698886873389072266E+13
       coef( 1117) =    -0.18476220624582441406E+14
       coef( 1118) =    -0.11281783503281884766E+12
       coef( 1119) =    -0.10363333413972566406E+14
       coef( 1120) =     0.11361171921223906250E+14
       coef( 1121) =    -0.70950590380329765625E+13
       coef( 1122) =     0.57725231467435179688E+14
       coef( 1123) =     0.74543701099865562500E+14
       coef( 1124) =     0.17156453327727580078E+14
       coef( 1125) =    -0.34311324306602412109E+13
       coef( 1126) =     0.34860909075131015625E+14
       coef( 1127) =    -0.83796354044848906250E+13
       coef( 1128) =     0.10613659036839152344E+14
       coef( 1129) =     0.41529425817991586914E+13
       coef( 1130) =    -0.17824636310162353516E+13
       coef( 1131) =     0.33849816349837138672E+13
       coef( 1132) =     0.60067632045957515625E+14
       coef( 1133) =     0.86249002223629960938E+13
       coef( 1134) =     0.33439307707673491211E+13
       coef( 1135) =     0.13065996160945156250E+13
       coef( 1136) =     0.29180401644492585938E+14
       coef( 1137) =     0.49166444966900351562E+13
       coef( 1138) =     0.14361908290938534737E+10
       coef( 1139) =    -0.96280513634093383789E+11
       coef( 1140) =     0.81698761918767285156E+12
       coef( 1141) =     0.99114817256189375000E+13
       coef( 1142) =    -0.61282817414476093750E+14
       coef( 1143) =     0.60151741482730890625E+14
       coef( 1144) =     0.67806034637969804688E+14
       coef( 1145) =     0.33943998963627968750E+14
       coef( 1146) =     0.18460667623503585938E+14
       coef( 1147) =     0.25664106507594638062E+12
       coef( 1148) =    -0.21478462867051000000E+14
       coef( 1149) =     0.81306818287321093750E+12
       coef( 1150) =    -0.14649193454805518750E+15
       coef( 1151) =    -0.13656531486137160156E+14
       coef( 1152) =    -0.92667723552870839844E+13
       coef( 1153) =    -0.54766016103644238281E+13
       coef( 1154) =     0.72413270885413656250E+14
       coef( 1155) =     0.36754091202183234375E+14
       coef( 1156) =     0.10217552777953820312E+15
       coef( 1157) =     0.66328529539445148438E+14
       coef( 1158) =     0.19803094002989125000E+14
       coef( 1159) =     0.35807707301730031250E+14
       coef( 1160) =     0.30166560633056148438E+14
       coef( 1161) =     0.14385355460598300781E+14
       coef( 1162) =     0.81633542168307490234E+13
       coef( 1163) =     0.14411803168904925781E+14
       coef( 1164) =     0.68343063891309765625E+12
       coef( 1165) =    -0.39739922343480320312E+14
       coef( 1166) =    -0.33865210914969972656E+14
       coef( 1167) =    -0.20085225210359632812E+14
       coef( 1168) =    -0.13274835407032369141E+14
       coef( 1169) =    -0.44311426472918914062E+14
       coef( 1170) =     0.13290639946824138672E+14
       coef( 1171) =     0.24706237767426414062E+14
       coef( 1172) =     0.87766358676214589844E+13
       coef( 1173) =     0.17509813335470324219E+14
       coef( 1174) =     0.11420312734891447266E+14
       coef( 1175) =    -0.78053798493652539062E+13
       coef( 1176) =    -0.42826209509900888672E+13
       coef( 1177) =     0.45500302046508583984E+13
       coef( 1178) =    -0.35449389372504243164E+13
       coef( 1179) =    -0.62840313976622593750E+14
       coef( 1180) =    -0.26863463613573792969E+14
       coef( 1181) =     0.28367444813322757812E+14
       coef( 1182) =     0.29872778879846074219E+14
       coef( 1183) =     0.82284471853649550781E+13
       coef( 1184) =     0.16988320397571759766E+14
       coef( 1185) =    -0.40064239133631054688E+12
       coef( 1186) =     0.37910372236029648438E+13
       coef( 1187) =    -0.72752211352425781250E+12
       coef( 1188) =     0.22408399903568531250E+14
       coef( 1189) =     0.27770395217406743164E+13
       coef( 1190) =     0.97471063551196992188E+13
       coef( 1191) =    -0.55136806252221794128E+10
       coef( 1192) =     0.60933979153543371582E+12
       coef( 1193) =    -0.72805308930953779297E+13
       coef( 1194) =     0.36971596241872363281E+12
       coef( 1195) =     0.89232249437228343750E+14
       coef( 1196) =    -0.33053442961284949219E+14
       coef( 1197) =    -0.14039808835328730469E+14
       coef( 1198) =    -0.14830663524297656250E+12
       coef( 1199) =    -0.34622725493718666992E+13
       coef( 1200) =     0.68755916838725718750E+14
       coef( 1201) =     0.76634829498975953125E+14
       coef( 1202) =     0.57503019398700742188E+14
       coef( 1203) =     0.37325816985509812500E+14
       coef( 1204) =     0.63549937898040146484E+13
       coef( 1205) =    -0.20883095865832368750E+15
       coef( 1206) =     0.18137735390009988281E+14
       coef( 1207) =     0.51071914598036742188E+14
       coef( 1208) =     0.29096438909124972656E+14
       coef( 1209) =     0.94326523255976250000E+13
       coef( 1210) =     0.91735605696155136719E+13
       coef( 1211) =    -0.24962111379710320312E+14
       coef( 1212) =     0.14829802320486806641E+14
       coef( 1213) =     0.49370899002506109375E+14
       coef( 1214) =     0.26142651303475007812E+14
       coef( 1215) =     0.72928587808439023438E+13
       coef( 1216) =    -0.29169270287728671875E+14
       coef( 1217) =     0.14689617046168568359E+14
       coef( 1218) =     0.13155673640491867188E+14
       coef( 1219) =     0.83923622060095527344E+13
       coef( 1220) =     0.35982763356803242188E+13
       coef( 1221) =     0.31499659151866640625E+13
       coef( 1222) =     0.13853691763505495312E+15
       coef( 1223) =     0.30300632676072023438E+14
       coef( 1224) =     0.31012977842441609375E+14
       coef( 1225) =     0.18650539366779363281E+14
       coef( 1226) =     0.14832603199264847656E+14
       coef( 1227) =     0.35225310664990605469E+13
       coef( 1228) =     0.96541795369431875000E+13
       coef( 1229) =    -0.78022034519545755386E+10
       coef( 1230) =    -0.11735529625480861816E+13
       coef( 1231) =     0.20685524317656609375E+14
       coef( 1232) =    -0.64067273483394882812E+14
       coef( 1233) =    -0.29720850865520839844E+13
       coef( 1234) =    -0.10587078112820007812E+15
       coef( 1235) =    -0.53779400546348320312E+14
       coef( 1236) =     0.88891199806694238281E+13
       coef( 1237) =    -0.96330122580885468750E+14
       coef( 1238) =    -0.76827678234654218750E+14
       coef( 1239) =     0.24848331325531562500E+14
       coef( 1240) =     0.17611922675017781250E+14
       coef( 1241) =     0.70054595070754140625E+14
       coef( 1242) =     0.38961895846811421875E+14
       coef( 1243) =     0.21428685016624765625E+14
       coef( 1244) =     0.70776582654276953125E+13
       coef( 1245) =    -0.18626601585120613281E+14
       coef( 1246) =     0.17898962195631523438E+14
       coef( 1247) =     0.23422321204589210938E+14
       coef( 1248) =     0.18028664867636214844E+14
       coef( 1249) =     0.22028950585359296875E+13
       coef( 1250) =     0.98954389284986445312E+13
       coef( 1251) =     0.56521969970888525391E+13
       coef( 1252) =    -0.48148284497179703125E+14
       coef( 1253) =    -0.32675352482918232422E+13
       coef( 1254) =     0.12265539219586460938E+14
       coef( 1255) =     0.39136085993356352539E+13
       coef( 1256) =     0.77892081778948699951E+11
       coef( 1257) =    -0.55116834832233935547E+12
       coef( 1258) =    -0.15696892663005039062E+14
       coef( 1259) =     0.75481670438735625000E+14
       coef( 1260) =     0.37845049910759902344E+13
       coef( 1261) =    -0.68022233236009382812E+14
       coef( 1262) =    -0.33448815632267353516E+13
       coef( 1263) =     0.50609002341990171875E+14
       coef( 1264) =    -0.92199566859951296875E+14
       coef( 1265) =    -0.77352015396766308594E+13
       coef( 1266) =     0.18179301784495896875E+15
       coef( 1267) =     0.36546413846249890625E+14
       coef( 1268) =     0.68272840070101640625E+14
       coef( 1269) =     0.43180811612459003906E+13
       coef( 1270) =    -0.51247773382661875000E+13
       coef( 1271) =     0.59062369168214921875E+13
       coef( 1272) =    -0.12098859208396645312E+15
       coef( 1273) =    -0.22847132776731039062E+14
       coef( 1274) =    -0.10893464013261856079E+12
       coef( 1275) =     0.27230308201325244141E+13
       coef( 1276) =    -0.10221342915421585938E+14
       coef( 1277) =     0.20322836017519468750E+14
       coef( 1278) =     0.16636537885194536133E+13
       coef( 1279) =    -0.81846434557844404297E+13
       coef( 1280) =    -0.56048187736911367188E+13
       coef( 1281) =    -0.83465372039255593750E+14
       coef( 1282) =     0.10895840490035871875E+15
       coef( 1283) =    -0.47190937398970000000E+13
       coef( 1284) =    -0.26431490745812230469E+14
       coef( 1285) =    -0.10394979522478264062E+15
       coef( 1286) =    -0.26618474044452605085E+04
       coef( 1287) =     0.31589278851358033717E+07
       coef( 1288) =    -0.31694703237590312958E+09
       coef( 1289) =     0.14083430289359613419E+11
       coef( 1290) =    -0.22242099452810610962E+12
       coef( 1291) =     0.14306579816538017578E+13
       coef( 1292) =    -0.50298538300702509766E+13
       coef( 1293) =     0.10905171737817169922E+14
       coef( 1294) =    -0.92602432682131171875E+13
       coef( 1295) =     0.14454196767005902529E+09
       coef( 1296) =    -0.19393162041956367493E+11
       coef( 1297) =     0.57459083121857006836E+12
       coef( 1298) =    -0.30164594295609902344E+13
       coef( 1299) =     0.10235902074370255859E+14
       coef( 1300) =    -0.26986042097909832031E+14
       coef( 1301) =    -0.37645642593216035156E+13
       coef( 1302) =     0.44636568861000015625E+14
       coef( 1303) =    -0.72099673217291366577E+11
       coef( 1304) =    -0.61890344310796951294E+11
       coef( 1305) =    -0.14214166583514269531E+14
       coef( 1306) =     0.47168232409690703125E+14
       coef( 1307) =     0.16126878314027929688E+12
       coef( 1308) =    -0.77861308031275078125E+14
       coef( 1309) =    -0.41156199673597773438E+14
       coef( 1310) =     0.50477270317201191406E+13
       coef( 1311) =     0.98231980472294550781E+13
       coef( 1312) =    -0.92080173550727812500E+13
       coef( 1313) =     0.38713838242629984375E+14
       coef( 1314) =     0.17017376596740066406E+14
       coef( 1315) =     0.75961227941312119141E+13
       coef( 1316) =    -0.34269444937020351562E+14
       coef( 1317) =    -0.27797869532715570312E+14
       coef( 1318) =     0.93193222929661054688E+13
       coef( 1319) =     0.10111320935938240234E+14
       coef( 1320) =    -0.18527324535009226562E+14
       coef( 1321) =    -0.16627165694661591797E+13
       coef( 1322) =     0.72647586683646316528E+10
       coef( 1323) =    -0.43710553214949633789E+12
       coef( 1324) =    -0.12566185074706040039E+13
       coef( 1325) =     0.81646328422607939453E+13
       coef( 1326) =     0.99420487443081640625E+13
       coef( 1327) =     0.42179190613796679688E+13
       coef( 1328) =     0.14272669917284621094E+14
       coef( 1329) =     0.19418002135654433594E+14
       coef( 1330) =     0.33070314019268193359E+13
       coef( 1331) =     0.15730625678323146484E+14
       coef( 1332) =    -0.60770451669025406250E+14
       coef( 1333) =    -0.52384182875496796875E+14
       coef( 1334) =    -0.50407131080813789062E+14
       coef( 1335) =    -0.32564437634626929688E+14
       coef( 1336) =    -0.13050987557853585938E+14
       coef( 1337) =    -0.62684160900143062500E+14
       coef( 1338) =     0.89956473675262406250E+14
       coef( 1339) =     0.64179639550412281250E+14
       coef( 1340) =     0.23374053346112941406E+14
       coef( 1341) =     0.53542114153884199219E+13
       coef( 1342) =     0.77844882832000109375E+14
       coef( 1343) =     0.35733629850393812500E+14
       coef( 1344) =     0.13150802850545970703E+14
       coef( 1345) =     0.12712332484905382812E+14
       coef( 1346) =    -0.19081909117847585938E+14
       coef( 1347) =     0.82475009385672062500E+14
       coef( 1348) =     0.66428000750694437500E+14
       coef( 1349) =     0.56313244721098564453E+13
       coef( 1350) =    -0.15670072943674148438E+14
       coef( 1351) =    -0.11956035721512302734E+14
       coef( 1352) =    -0.22221151164388671875E+14
       coef( 1353) =     0.11162271537802791016E+14
       coef( 1354) =     0.67919637671817226562E+13
       coef( 1355) =     0.41490149148325823975E+12
       coef( 1356) =     0.71162234729992871094E+13
       coef( 1357) =     0.34753936771908959961E+13
       coef( 1358) =    -0.61770456471608125000E+14
       coef( 1359) =    -0.17462248594091109375E+14
       coef( 1360) =    -0.33384445142727773438E+13
       coef( 1361) =    -0.59107874625656484375E+13
       coef( 1362) =     0.26287152538143893433E+12
       coef( 1363) =    -0.35156753168970117188E+13
       coef( 1364) =     0.17636781668439851562E+14
       coef( 1365) =    -0.45579102826322679688E+14
       coef( 1366) =     0.14983137273953986328E+14
       coef( 1367) =    -0.98233972238654863281E+13
       coef( 1368) =    -0.68183532295280136719E+13
       coef( 1369) =     0.89363347705287158203E+12
       coef( 1370) =     0.14318155410642923828E+14
       coef( 1371) =    -0.12679640042631739062E+15
       coef( 1372) =     0.69828209563521468750E+14
       coef( 1373) =     0.58872497714899781250E+14
       coef( 1374) =     0.14648324058404839844E+14
       coef( 1375) =     0.75897676699845654297E+12
       coef( 1376) =     0.84188515954917578125E+14
       coef( 1377) =     0.84940830329860984375E+14
       coef( 1378) =     0.39327060469538468750E+14
       coef( 1379) =     0.12071796443360205078E+14
       coef( 1380) =     0.38714135294050890625E+14
       coef( 1381) =     0.13676052915150572266E+14
       coef( 1382) =     0.16101404509643337891E+14
       coef( 1383) =     0.46242343493917500000E+14
       coef( 1384) =     0.53554828708519453125E+14
       coef( 1385) =     0.24171660360892394531E+14
       coef( 1386) =     0.56247043107673574219E+13
       coef( 1387) =     0.22559917714594554688E+14
       coef( 1388) =     0.14209814975478861328E+14
       coef( 1389) =     0.59687163757745751953E+13
       coef( 1390) =     0.43647058685554384766E+13
       coef( 1391) =    -0.23169596102814790039E+13
       coef( 1392) =    -0.40232872508684722900E+12
       coef( 1393) =    -0.27196610682888210938E+14
       coef( 1394) =    -0.10099205834982535156E+14
       coef( 1395) =     0.68276207335457343750E+13
       coef( 1396) =     0.42884815742628842773E+13
       coef( 1397) =    -0.82614382554650371094E+13
       coef( 1398) =    -0.81387427768863037109E+12
       coef( 1399) =    -0.48473222840519443359E+13
       coef( 1400) =     0.16457306050914905500E+07
       coef( 1401) =    -0.12003693658087346703E+08
       coef( 1402) =    -0.41884498239852671623E+10
       coef( 1403) =     0.18203767082257040405E+12
       coef( 1404) =    -0.17029988748056091309E+13
       coef( 1405) =     0.80678935523036083984E+13
       coef( 1406) =    -0.21176751633663050781E+14
       coef( 1407) =     0.44914285486881500000E+14
       coef( 1408) =    -0.51970964286422468750E+14
       coef( 1409) =     0.32942216872651052475E+10
       coef( 1410) =    -0.37605500079439227295E+12
       coef( 1411) =     0.37259383845181469727E+13
       coef( 1412) =    -0.16998612507288214844E+14
       coef( 1413) =     0.82024611262562500000E+13
       coef( 1414) =     0.21696205004859375000E+14
       coef( 1415) =     0.44399341724556312500E+14
       coef( 1416) =     0.12327368489921109375E+14
       coef( 1417) =     0.72452802756998315430E+12
       coef( 1418) =    -0.93808508041049003906E+13
       coef( 1419) =    -0.17440988272368517578E+14
       coef( 1420) =     0.22231932081742863281E+14
       coef( 1421) =     0.10600633038189109375E+14
       coef( 1422) =    -0.10743860636907400391E+14
       coef( 1423) =    -0.13132454438848679688E+14
       coef( 1424) =     0.66423319001862015625E+14
       coef( 1425) =    -0.50001288072596851562E+14
       coef( 1426) =    -0.71284820015896250000E+14
       coef( 1427) =    -0.31894325201322074219E+14
       coef( 1428) =    -0.13292714902502349609E+14
       coef( 1429) =    -0.44323338237859867188E+14
       coef( 1430) =    -0.35264142839499851562E+14
       coef( 1431) =    -0.14074221610256755859E+14
       coef( 1432) =    -0.16075275883438597656E+14
       coef( 1433) =     0.14204638413768252563E+12
       coef( 1434) =    -0.53104465301811547852E+12
       coef( 1435) =    -0.29240974915894628906E+12
       coef( 1436) =     0.64439507415401062500E+14
       coef( 1437) =    -0.98009226775112437500E+14
       coef( 1438) =    -0.99682638439516750000E+14
       coef( 1439) =    -0.54537619659010523438E+14
       coef( 1440) =    -0.25303531722825382812E+14
       coef( 1441) =    -0.61018907855594921875E+12
       coef( 1442) =    -0.46748097063315406250E+14
       coef( 1443) =     0.11648654056683025000E+15
       coef( 1444) =     0.37398674879677382812E+14
       coef( 1445) =    -0.68951429647078886719E+13
       coef( 1446) =    -0.11604732944044054688E+14
       coef( 1447) =    -0.25971944691942914062E+14
       coef( 1448) =     0.42179720950745460938E+14
       coef( 1449) =     0.16443053678413472656E+14
       coef( 1450) =     0.16779340050174858398E+13
       coef( 1451) =     0.20246041502273703125E+14
       coef( 1452) =     0.54589861502835214844E+13
       coef( 1453) =     0.24245936928070023438E+14
       coef( 1454) =    -0.57309442282479789062E+14
       coef( 1455) =     0.44926795033278625000E+14
       coef( 1456) =     0.21179675345057492188E+14
       coef( 1457) =     0.23355142142905888672E+13
       coef( 1458) =    -0.51774036695847359375E+14
       coef( 1459) =    -0.21161497976152062988E+12
       coef( 1460) =     0.29851152328660444336E+13
       coef( 1461) =     0.95181107973740722656E+12
       coef( 1462) =    -0.31637390487198535156E+14
       coef( 1463) =    -0.59249489814623476562E+13
       coef( 1464) =    -0.10420893373334082031E+13
       coef( 1465) =     0.15452977129798376953E+14
       coef( 1466) =    -0.71371119550734328125E+14
       coef( 1467) =     0.10552409022168360938E+15
       coef( 1468) =     0.20788289408318484375E+14
       coef( 1469) =    -0.19641927485820304688E+14
       coef( 1470) =    -0.17015991425286238281E+14
       coef( 1471) =    -0.46754371912050492188E+14
       coef( 1472) =    -0.50802492311341078125E+14
       coef( 1473) =     0.27168282604401859375E+14
       coef( 1474) =     0.12590405491584548828E+14
       coef( 1475) =    -0.39229412275205419922E+12
       coef( 1476) =     0.16394201591365664062E+13
       coef( 1477) =     0.13354567902683296875E+14
       coef( 1478) =     0.55481332581869218750E+13
       coef( 1479) =     0.56469466730527871094E+13
       coef( 1480) =     0.11916499911703084375E+15
       coef( 1481) =     0.27324511409608398438E+14
       coef( 1482) =     0.13545680522376128906E+14
       coef( 1483) =     0.45370238818226142578E+13
       coef( 1484) =     0.64183363989175683594E+13
       coef( 1485) =     0.27201102467723144531E+13
       coef( 1486) =    -0.82854545397324218750E+10
       coef( 1487) =     0.66522970345149531250E+14
       coef( 1488) =     0.12604065698566126953E+14
       coef( 1489) =     0.31948777015791450195E+13
       coef( 1490) =     0.16893543840889951172E+13
       coef( 1491) =    -0.47343028259245734662E+07
       coef( 1492) =     0.80400291434611713886E+09
       coef( 1493) =     0.51902491872886772156E+10
       coef( 1494) =    -0.35421216315746331787E+12
       coef( 1495) =     0.44997697477473388672E+13
       coef( 1496) =    -0.21610324301968402344E+14
       coef( 1497) =     0.29940831374130792969E+14
       coef( 1498) =     0.31111361853032515625E+14
       coef( 1499) =     0.12546655315548691406E+13
       coef( 1500) =    -0.63139511376782974243E+11
       coef( 1501) =     0.22142935576795522461E+13
       coef( 1502) =    -0.16366357858452271484E+14
       coef( 1503) =     0.47397040281853960938E+14
       coef( 1504) =    -0.31971118034072921875E+14
       coef( 1505) =    -0.42621578077600765625E+14
       coef( 1506) =    -0.14483198167412832031E+14
       coef( 1507) =    -0.40489584091993945312E+13
       coef( 1508) =    -0.13084912755273130859E+14
       coef( 1509) =     0.10607711481518559375E+15
       coef( 1510) =    -0.88610814998245687500E+14
       coef( 1511) =    -0.10882666717855708984E+14
       coef( 1512) =    -0.77607085999230908203E+12
       coef( 1513) =    -0.32922631852569345703E+13
       coef( 1514) =    -0.11559068252604943750E+15
       coef( 1515) =    -0.11837692529166610938E+15
       coef( 1516) =    -0.44073777139338421875E+14
       coef( 1517) =    -0.12991473703788181641E+14
       coef( 1518) =    -0.62384006787907039062E+14
       coef( 1519) =    -0.21779552706587203125E+14
       coef( 1520) =    -0.18216643879389758301E+12
       coef( 1521) =     0.23944857356059858398E+13
       coef( 1522) =     0.23415630261202953125E+14
       coef( 1523) =    -0.87275836828867109375E+14
       coef( 1524) =    -0.11558650011137017188E+15
       coef( 1525) =    -0.68142607179832046875E+14
       coef( 1526) =    -0.27844561200827957031E+14
       coef( 1527) =    -0.16908183754213246094E+14
       coef( 1528) =     0.17169436431441246094E+14
       coef( 1529) =     0.33455638047147101562E+14
       coef( 1530) =     0.44625883777408427734E+13
       coef( 1531) =    -0.47040884591222558594E+13
       coef( 1532) =    -0.96024619372341113281E+12
       coef( 1533) =     0.92135346688822226562E+13
       coef( 1534) =     0.28192116944845131836E+13
       coef( 1535) =     0.35486501680246230469E+13

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

