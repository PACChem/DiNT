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

      resp=0.001d0

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

       ncoef=  540
       coef(    1) =    -0.13932968537648162677E+01
       coef(    2) =    -0.39888708523829293995E+02
       coef(    3) =     0.33455236187006055843E+05
       coef(    4) =     0.46196942171875335043E+06
       coef(    5) =    -0.15356968256780984998E+09
       coef(    6) =     0.13423572722665503025E+10
       coef(    7) =    -0.31988089096266560555E+10
       coef(    8) =    -0.61459677142324071610E+05
       coef(    9) =    -0.86015588652576711029E+07
       coef(   10) =     0.10931366293683180809E+10
       coef(   11) =    -0.42630928006393694878E+09
       coef(   12) =    -0.51475581648276321411E+11
       coef(   13) =     0.18217479085460620117E+12
       coef(   14) =    -0.14022031378255245090E+09
       coef(   15) =    -0.15717195745554462433E+11
       coef(   16) =     0.81686664766765350342E+11
       coef(   17) =    -0.44439711188324729919E+11
       coef(   18) =    -0.10534802033515794678E+13
       coef(   19) =     0.12497024875719430542E+12
       coef(   20) =    -0.10910356114568178711E+13
       coef(   21) =     0.39113604928235351562E+13
       coef(   22) =     0.27469716314730682373E+12
       coef(   23) =     0.77947035342364575195E+12
       coef(   24) =    -0.83018904711029746094E+13
       coef(   25) =     0.54948328444052734375E+12
       coef(   26) =     0.85454415686162500000E+13
       coef(   27) =     0.96052038602794706821E+07
       coef(   28) =    -0.14088001234705724716E+10
       coef(   29) =    -0.33052501070635101318E+11
       coef(   30) =     0.19943532082611053467E+12
       coef(   31) =     0.74191716553178863525E+11
       coef(   32) =    -0.17390347709856623535E+13
       coef(   33) =     0.71934931869884292603E+11
       coef(   34) =    -0.13966488598493814087E+12
       coef(   35) =     0.29853704890312457275E+12
       coef(   36) =    -0.40614928950260776367E+13
       coef(   37) =     0.19875720223320894531E+14
       coef(   38) =    -0.35279751397025488281E+13
       coef(   39) =     0.22641729016957460938E+14
       coef(   40) =    -0.61199152979299914062E+14
       coef(   41) =     0.12052764720716191406E+14
       coef(   42) =     0.84425583880457578125E+13
       coef(   43) =    -0.12814114783972500000E+13
       coef(   44) =    -0.70803911143165625000E+12
       coef(   45) =     0.95512556887867714844E+13
       coef(   46) =    -0.40061939557464773438E+14
       coef(   47) =     0.90382175312261500000E+14
       coef(   48) =    -0.74123768211963187500E+14
       coef(   49) =    -0.14502856204503699219E+14
       coef(   50) =     0.88003034295626062500E+14
       coef(   51) =    -0.10993546441219453125E+15
       coef(   52) =    -0.45276899610070093750E+14
       coef(   53) =    -0.18630177504767140625E+14
       coef(   54) =     0.25810787797972156250E+14
       coef(   55) =     0.49976653266676622629E+09
       coef(   56) =    -0.17775247980310379028E+11
       coef(   57) =     0.48577749287781469727E+12
       coef(   58) =    -0.42907348722957070312E+13
       coef(   59) =     0.11622593168501462891E+14
       coef(   60) =    -0.76609873926409189453E+13
       coef(   61) =    -0.16305020052604919434E+12
       coef(   62) =    -0.19222786531923647461E+13
       coef(   63) =     0.40084171893817429688E+14
       coef(   64) =    -0.12152386500079331250E+15
       coef(   65) =     0.77261127607929593750E+14
       coef(   66) =     0.19988247026473578125E+14
       coef(   67) =    -0.17482770673015378125E+15
       coef(   68) =     0.36311159623128800000E+15
       coef(   69) =    -0.96312795214452148438E+12
       coef(   70) =     0.32764147420984619141E+13
       coef(   71) =    -0.60516287534566312500E+14
       coef(   72) =     0.11132866177124378125E+15
       coef(   73) =    -0.13748957713603131250E+15
       coef(   74) =     0.13049725238887603125E+15
       coef(   75) =    -0.91943041475870968750E+14
       coef(   76) =     0.49720234971010546875E+12
       coef(   77) =     0.26482699346679418945E+13
       coef(   78) =     0.31888996165068445312E+14
       coef(   79) =    -0.90209697754822656250E+13
       coef(   80) =    -0.17928738173016950000E+15
       coef(   81) =    -0.11204165193344992585E+07
       coef(   82) =     0.10562023419014157844E+06
       coef(   83) =     0.46346001167789020110E+06
       coef(   84) =     0.38106423706642840989E+07
       coef(   85) =     0.15869877873209846020E+10
       coef(   86) =    -0.60508771910363025665E+10
       coef(   87) =    -0.21448938450421314240E+10
       coef(   88) =     0.38888307176824212074E+07
       coef(   89) =    -0.38900966316529494524E+09
       coef(   90) =    -0.77530258354755916595E+10
       coef(   91) =    -0.35529293426255577087E+11
       coef(   92) =     0.38136134474820227051E+12
       coef(   93) =    -0.62253629736938708496E+12
       coef(   94) =     0.10073911302508897781E+11
       coef(   95) =    -0.62667361608866027832E+11
       coef(   96) =     0.14075130138786288452E+12
       coef(   97) =    -0.14509325562037451172E+12
       coef(   98) =     0.47366645382230957031E+13
       coef(   99) =     0.74306934886365014648E+12
       coef(  100) =    -0.70803666457785859375E+13
       coef(  101) =    -0.22530995550808769531E+13
       coef(  102) =    -0.12829852562525941406E+14
       coef(  103) =     0.22789017010409941406E+14
       coef(  104) =    -0.20208162685568468750E+14
       coef(  105) =     0.18932980262961021066E+09
       coef(  106) =     0.17025544450572547913E+11
       coef(  107) =     0.31801429197689593506E+12
       coef(  108) =    -0.10894710116432200928E+13
       coef(  109) =    -0.26355440747446069336E+13
       coef(  110) =     0.82439614679452558594E+13
       coef(  111) =    -0.71497971028929162598E+12
       coef(  112) =     0.77554134537841638184E+12
       coef(  113) =     0.89110212491022089844E+13
       coef(  114) =     0.11047967303407835938E+14
       coef(  115) =    -0.79300073040309046875E+14
       coef(  116) =     0.14687493610147734375E+14
       coef(  117) =    -0.11773742660041296875E+15
       coef(  118) =     0.22169740296028240625E+15
       coef(  119) =     0.10655379637866250000E+14
       coef(  120) =     0.61906087603869394531E+13
       coef(  121) =    -0.66165197519916859375E+14
       coef(  122) =     0.90147874446247406250E+14
       coef(  123) =    -0.82721547844149843750E+14
       coef(  124) =     0.17253288354650606250E+15
       coef(  125) =    -0.21856554569246100000E+15
       coef(  126) =    -0.78489253889762515625E+14
       coef(  127) =    -0.17856508455593997955E+11
       coef(  128) =     0.48587325018586279297E+12
       coef(  129) =    -0.48792945597833652344E+13
       coef(  130) =     0.16926325743326175781E+14
       coef(  131) =    -0.16853431195010001953E+14
       coef(  132) =     0.39853533824659082031E+12
       coef(  133) =    -0.22483565634590043945E+13
       coef(  134) =     0.39456256827561578125E+14
       coef(  135) =    -0.11637004006658353125E+15
       coef(  136) =     0.10712596921176159375E+15
       coef(  137) =    -0.25096705697074832031E+14
       coef(  138) =     0.83730707575176031250E+14
       coef(  139) =    -0.23625789213133212891E+13
       coef(  140) =    -0.63138446653040367188E+14
       coef(  141) =     0.11111065545011215625E+15
       coef(  142) =    -0.48581935100839781250E+14
       coef(  143) =     0.15787642693899593750E+14
       coef(  144) =     0.61786849223399539062E+14
       coef(  145) =    -0.12988617506765331200E+04
       coef(  146) =    -0.34305620374234248884E+07
       coef(  147) =     0.20010742631976366043E+09
       coef(  148) =    -0.54576147560163676739E+08
       coef(  149) =    -0.51396650323546104431E+11
       coef(  150) =     0.28036294552503564453E+12
       coef(  151) =    -0.31500878366055419922E+12
       coef(  152) =     0.86289032140694230795E+08
       coef(  153) =    -0.99197199133578643799E+10
       coef(  154) =     0.25508218203364831543E+12
       coef(  155) =     0.45249046103725891113E+11
       coef(  156) =    -0.44148910724153027344E+13
       coef(  157) =     0.59608464106937890625E+13
       coef(  158) =     0.17048595610792289734E+11
       coef(  159) =    -0.93610945202134899902E+12
       coef(  160) =    -0.93725872175437734375E+13
       coef(  161) =     0.27324350320965050781E+14
       coef(  162) =    -0.21292614647119800781E+14
       coef(  163) =     0.66463472430775166016E+13
       coef(  164) =     0.10657065877975207031E+14
       coef(  165) =     0.97918137600439257812E+13
       coef(  166) =    -0.54959785220328125000E+14
       coef(  167) =     0.52788256655893731117E+09
       coef(  168) =    -0.28467830421541406250E+12
       coef(  169) =    -0.12084169082583591309E+13
       coef(  170) =     0.16915674848682851562E+14
       coef(  171) =    -0.32841459061992964844E+14
       coef(  172) =     0.25143535219723496094E+14
       coef(  173) =     0.22980464761764853516E+13
       coef(  174) =    -0.80655191581782480469E+13
       coef(  175) =     0.15405047014640679688E+14
       coef(  176) =    -0.86773713170533781250E+14
       coef(  177) =    -0.46227468111417070312E+14
       coef(  178) =     0.15083994887215815625E+15
       coef(  179) =    -0.75839451295603847656E+13
       coef(  180) =     0.74152867067728406250E+14
       coef(  181) =     0.16421027227353552734E+14
       coef(  182) =    -0.18390827974735740625E+15
       coef(  183) =     0.50826780285126123047E+12
       coef(  184) =    -0.59166269715071835938E+13
       coef(  185) =     0.29890458329844210938E+14
       coef(  186) =    -0.13228178425227254688E+15
       coef(  187) =     0.17861632112347387500E+15
       coef(  188) =     0.18279883273741222656E+14
       coef(  189) =    -0.13217460939773251562E+15
       coef(  190) =     0.87300223932101781250E+14
       coef(  191) =    -0.73013182515507470703E+13
       coef(  192) =     0.38132802640286414062E+14
       coef(  193) =     0.12347802079637825000E+15
       coef(  194) =    -0.11905752324834196875E+15
       coef(  195) =     0.37084868558805587236E+06
       coef(  196) =     0.33174701107169337571E+08
       coef(  197) =    -0.59685204650357007980E+10
       coef(  198) =     0.52753063960365570068E+11
       coef(  199) =     0.82876178243924209595E+11
       coef(  200) =    -0.64768781860471191406E+12
       coef(  201) =    -0.43862549969817871094E+12
       coef(  202) =     0.11917605337153856754E+10
       coef(  203) =     0.14992492931211959839E+12
       coef(  204) =    -0.38114746687100419922E+13
       coef(  205) =     0.57567951311057050781E+13
       coef(  206) =     0.21009837947722218750E+14
       coef(  207) =    -0.17090640656955800781E+14
       coef(  208) =    -0.47170653541641522217E+12
       coef(  209) =     0.20565823845154437500E+14
       coef(  210) =     0.29110632692906484375E+14
       coef(  211) =    -0.87475577047331484375E+14
       coef(  212) =    -0.97162061769833187500E+14
       coef(  213) =     0.10424640277648420312E+15
       coef(  214) =    -0.49419475754876914978E+11
       coef(  215) =     0.44349678577647226562E+13
       coef(  216) =    -0.17125185882038806641E+14
       coef(  217) =    -0.39623995184728367188E+14
       coef(  218) =     0.54487692933113734375E+14
       coef(  219) =     0.74745693513075957031E+13
       coef(  220) =     0.19878514761772460938E+14
       coef(  221) =    -0.16224463958781726562E+14
       coef(  222) =     0.41479182271911453125E+14
       coef(  223) =    -0.60016945507437257812E+14
       coef(  224) =    -0.10598488936748703125E+15
       coef(  225) =    -0.25333106365297934570E+13
       coef(  226) =     0.43444047191705268555E+13
       coef(  227) =     0.88489512202720359375E+14
       coef(  228) =    -0.12810430695516066406E+14
       coef(  229) =    -0.60377049074900367188E+14
       coef(  230) =     0.98235410559355312500E+14
       coef(  231) =     0.72948968671742031250E+14
       coef(  232) =    -0.24065656930828625336E+07
       coef(  233) =     0.50090012145424795151E+09
       coef(  234) =     0.18663000482223979950E+11
       coef(  235) =    -0.95067857443411529541E+11
       coef(  236) =    -0.94831101772894226074E+11
       coef(  237) =    -0.14711796782507534180E+13
       coef(  238) =     0.41410676515653393555E+13
       coef(  239) =    -0.33016201214958545685E+11
       coef(  240) =    -0.14330955436199157715E+11
       coef(  241) =     0.12894544601572869141E+14
       coef(  242) =    -0.32387026655464714844E+14
       coef(  243) =     0.77134436356347875977E+12
       coef(  244) =    -0.67161609030302060547E+13
       coef(  245) =    -0.34837548141488125000E+13
       coef(  246) =    -0.41830866640788500000E+14
       coef(  247) =     0.11946434227212075000E+15
       coef(  248) =     0.25021126048073406982E+12
       coef(  249) =    -0.15703700847656482422E+14
       coef(  250) =     0.31804831452453453125E+14
       coef(  251) =     0.35452313801581765625E+14
       coef(  252) =    -0.50119469176010253906E+13
       coef(  253) =     0.54717431041882515625E+14
       coef(  254) =     0.16642629919888512500E+15
       coef(  255) =     0.92287562589568710938E+13
       coef(  256) =     0.92014603844668164062E+12
       coef(  257) =    -0.18880598464089375000E+15
       coef(  258) =    -0.42975782608652796875E+14
       coef(  259) =    -0.61496621206146523356E+08
       coef(  260) =     0.17628917884488921165E+10
       coef(  261) =    -0.62194158298392921448E+11
       coef(  262) =    -0.45493015795961596680E+12
       coef(  263) =     0.34547406126998129883E+13
       coef(  264) =    -0.10750633747102568359E+13
       coef(  265) =     0.68112469215657447815E+11
       coef(  266) =     0.58500408235530541992E+12
       coef(  267) =    -0.16730387949801710938E+14
       coef(  268) =     0.35658053570142109375E+14
       coef(  269) =     0.16494557965258761719E+14
       coef(  270) =    -0.70899197348436796875E+14
       coef(  271) =    -0.96374633296929602051E+12
       coef(  272) =     0.29503164245141277344E+14
       coef(  273) =    -0.19605471511008265625E+14
       coef(  274) =    -0.77631893458462921875E+14
       coef(  275) =    -0.19020195388097257812E+14
       coef(  276) =     0.58090285179431007812E+14
       coef(  277) =     0.28330870504679548740E+09
       coef(  278) =    -0.16613201857817089081E+11
       coef(  279) =     0.22564280073314166260E+12
       coef(  280) =     0.37122214297963024902E+12
       coef(  281) =    -0.44871383475669755859E+13
       coef(  282) =     0.15565676136891992188E+12
       coef(  283) =    -0.59331554932551611328E+13
       coef(  284) =     0.16439868275384185547E+14
       coef(  285) =     0.10802132557820351562E+14
       coef(  286) =     0.14531783977711003418E+13
       coef(  287) =    -0.13992129453817490234E+14
       coef(  288) =     0.84094521494718515625E+13
       coef(  289) =     0.11204169782797354273E+07
       coef(  290) =    -0.10953031455489824293E+06
       coef(  291) =    -0.10720042439382482553E+06
       coef(  292) =     0.41293594707737723365E+07
       coef(  293) =     0.15843710717460374832E+10
       coef(  294) =    -0.60516604435616798401E+10
       coef(  295) =    -0.21440888081163425446E+10
       coef(  296) =     0.30556418580374312587E+07
       coef(  297) =    -0.38851534741210210323E+09
       coef(  298) =    -0.77534441370117893219E+10
       coef(  299) =    -0.35529872868088340759E+11
       coef(  300) =     0.38136083953443896484E+12
       coef(  301) =    -0.62253437598758154297E+12
       coef(  302) =     0.10074069635095077515E+11
       coef(  303) =    -0.62665506897353515625E+11
       coef(  304) =     0.14075143447879531860E+12
       coef(  305) =    -0.14509395830595214844E+12
       coef(  306) =     0.47366649630781972656E+13
       coef(  307) =     0.74306863215176232910E+12
       coef(  308) =    -0.70803664929517255859E+13
       coef(  309) =    -0.22530991787788261719E+13
       coef(  310) =    -0.12829852605810843750E+14
       coef(  311) =     0.22789016149662652344E+14
       coef(  312) =    -0.20208162455189914062E+14
       coef(  313) =     0.18956246630215251446E+09
       coef(  314) =     0.17025384128120721817E+11
       coef(  315) =     0.31801445156924279785E+12
       coef(  316) =    -0.10894710411574235840E+13
       coef(  317) =    -0.26355445188206098633E+13
       coef(  318) =     0.82439614483984511719E+13
       coef(  319) =    -0.71497979542597460938E+12
       coef(  320) =     0.77554159487776928711E+12
       coef(  321) =     0.89110212246418710938E+13
       coef(  322) =     0.11047967306605490234E+14
       coef(  323) =    -0.79300072927864687500E+14
       coef(  324) =     0.14687493555961119141E+14
       coef(  325) =    -0.11773742657450215625E+15
       coef(  326) =     0.22169740298490146875E+15
       coef(  327) =     0.10655379621090648438E+14
       coef(  328) =     0.61906087299637373047E+13
       coef(  329) =    -0.66165197588305164062E+14
       coef(  330) =     0.90147874477387656250E+14
       coef(  331) =    -0.82721547813389359375E+14
       coef(  332) =     0.17253288346070009375E+15
       coef(  333) =    -0.21856554582812675000E+15
       coef(  334) =    -0.78489253848085562500E+14
       coef(  335) =    -0.17856442837207656860E+11
       coef(  336) =     0.48587338560702075195E+12
       coef(  337) =    -0.48792946758718564453E+13
       coef(  338) =     0.16926325652902183594E+14
       coef(  339) =    -0.16853431231477890625E+14
       coef(  340) =     0.39853527511634472656E+12
       coef(  341) =    -0.22483565306859389648E+13
       coef(  342) =     0.39456256952503640625E+14
       coef(  343) =    -0.11637004000752104688E+15
       coef(  344) =     0.10712596929138554688E+15
       coef(  345) =    -0.25096705713790281250E+14
       coef(  346) =     0.83730707571063953125E+14
       coef(  347) =    -0.23625789493775019531E+13
       coef(  348) =    -0.63138446697050460938E+14
       coef(  349) =     0.11111065540699539062E+15
       coef(  350) =    -0.48581935090771265625E+14
       coef(  351) =     0.15787642705876371094E+14
       coef(  352) =     0.61786849195680312500E+14
       coef(  353) =     0.83705187646592912643E+04
       coef(  354) =     0.34680714619563547894E+07
       coef(  355) =    -0.22032478329608395696E+09
       coef(  356) =    -0.44694061588135318756E+10
       coef(  357) =     0.12221671591084042358E+12
       coef(  358) =    -0.78828264063684655762E+12
       coef(  359) =     0.14400782184993876953E+13
       coef(  360) =    -0.41065163756528133154E+09
       coef(  361) =     0.40003515507580360413E+11
       coef(  362) =    -0.47647265667733081055E+12
       coef(  363) =     0.70747623086957568359E+12
       coef(  364) =     0.84936609390840458984E+13
       coef(  365) =    -0.25801690329201367188E+14
       coef(  366) =    -0.20152039411667538452E+12
       coef(  367) =     0.33575021815104829102E+13
       coef(  368) =     0.16605640259269044922E+14
       coef(  369) =    -0.84525069342214343750E+14
       coef(  370) =     0.12748335563814196875E+15
       coef(  371) =    -0.30538123394476046875E+14
       coef(  372) =     0.10050891626129893750E+15
       coef(  373) =    -0.67463975519988695312E+14
       coef(  374) =    -0.53718913252368937500E+14
       coef(  375) =    -0.66712202636473369598E+10
       coef(  376) =     0.98158442955910186768E+11
       coef(  377) =     0.25202055430414433594E+13
       coef(  378) =    -0.36091060543112781250E+14
       coef(  379) =     0.73439958624446312500E+14
       coef(  380) =     0.10782035495487056641E+14
       coef(  381) =    -0.18403045161940637207E+13
       coef(  382) =     0.30956275051424261719E+14
       coef(  383) =    -0.56782020238278515625E+14
       coef(  384) =     0.39751163776675906250E+14
       coef(  385) =    -0.20843583257305285156E+14
       coef(  386) =     0.22280701259534375000E+14
       coef(  387) =    -0.46362832295527109375E+13
       coef(  388) =     0.27288424126373671875E+14
       coef(  389) =    -0.44975707314696406250E+14
       coef(  390) =    -0.14311371207865846875E+15
       coef(  391) =    -0.10950332880913900757E+12
       coef(  392) =     0.19434387262860644531E+13
       coef(  393) =    -0.66046918609274287109E+13
       coef(  394) =     0.74788617606212593750E+14
       coef(  395) =    -0.20619455234654025000E+15
       coef(  396) =    -0.18708097246813261719E+14
       coef(  397) =     0.29249405636456519531E+14
       coef(  398) =     0.72233145243632781250E+14
       coef(  399) =    -0.18167754250963851562E+14
       coef(  400) =     0.55049719526102210938E+14
       coef(  401) =     0.99587139024682609375E+14
       coef(  402) =    -0.17873753674446375000E+15
       coef(  403) =     0.48702491848787770141E+06
       coef(  404) =     0.78275081135236203671E+08
       coef(  405) =     0.10961610186613721848E+10
       coef(  406) =     0.30150643681563568115E+11
       coef(  407) =    -0.41290566626158532715E+12
       coef(  408) =     0.20739939482118647461E+13
       coef(  409) =    -0.19709007913051909180E+13
       coef(  410) =    -0.11165789503818228245E+10
       coef(  411) =    -0.32643968848985394287E+12
       coef(  412) =     0.25849107259914287109E+13
       coef(  413) =    -0.11679393364703082031E+14
       coef(  414) =     0.18488592826693007812E+14
       coef(  415) =    -0.16882990433969931641E+14
       coef(  416) =     0.18150174233815725098E+13
       coef(  417) =    -0.17473449614932697266E+14
       coef(  418) =    -0.43003887134656039062E+14
       coef(  419) =     0.70109093065492281250E+14
       coef(  420) =     0.10815402940464154688E+15
       coef(  421) =    -0.13359191575697895312E+15
       coef(  422) =     0.12390442003270101929E+12
       coef(  423) =     0.15292706973512858887E+13
       coef(  424) =    -0.16729411072785683594E+13
       coef(  425) =     0.73584189157337437500E+14
       coef(  426) =    -0.56174032876677609375E+14
       coef(  427) =    -0.29016448168483648438E+14
       coef(  428) =     0.13615081644176359375E+14
       coef(  429) =    -0.83045849404177281250E+14
       coef(  430) =     0.84401559930615828125E+14
       coef(  431) =     0.44177691486977953125E+14
       coef(  432) =     0.25212369155736910156E+14
       coef(  433) =    -0.15081575119220446777E+13
       coef(  434) =     0.28277340651603187500E+14
       coef(  435) =    -0.10916481954234246875E+15
       coef(  436) =     0.10609313415988653125E+15
       coef(  437) =    -0.45229092877258265625E+14
       coef(  438) =    -0.36044065956929031250E+14
       coef(  439) =     0.40452631491758421875E+14
       coef(  440) =    -0.34037838847370423377E+08
       coef(  441) =    -0.64117297717748439312E+09
       coef(  442) =     0.25311757775731708527E+11
       coef(  443) =    -0.42881842044413488770E+12
       coef(  444) =     0.24941305450746083984E+13
       coef(  445) =    -0.90557101967124296875E+13
       coef(  446) =     0.71744668857062470703E+13
       coef(  447) =     0.55770196630918006897E+11
       coef(  448) =     0.25340831570202493286E+12
       coef(  449) =     0.90676924062666320801E+12
       coef(  450) =     0.16094537124106070312E+14
       coef(  451) =    -0.35520418153220750000E+14
       coef(  452) =    -0.13096749581043110352E+13
       coef(  453) =    -0.86360399827939726562E+13
       coef(  454) =     0.11580402006847482812E+15
       coef(  455) =    -0.15477500290277456250E+15
       coef(  456) =    -0.98366486831471044922E+12
       coef(  457) =    -0.10309406520029046875E+14
       coef(  458) =     0.11801933490780273438E+14
       coef(  459) =    -0.94306762497378031250E+14
       coef(  460) =     0.14334667257944665625E+15
       coef(  461) =    -0.17969660569806212500E+15
       coef(  462) =    -0.12795840000566809375E+15
       coef(  463) =     0.69604040631678046875E+13
       coef(  464) =    -0.79622460367284328125E+14
       coef(  465) =     0.13121117532362478125E+15
       coef(  466) =     0.17349200882077234375E+15
       coef(  467) =     0.54042416954686725140E+09
       coef(  468) =    -0.17515478569809825897E+11
       coef(  469) =    -0.56360514548903648376E+11
       coef(  470) =     0.91426214786691320801E+12
       coef(  471) =    -0.30871194222752324219E+13
       coef(  472) =     0.11871197989938947266E+14
       coef(  473) =    -0.21511142080503211975E+11
       coef(  474) =     0.42763099832450000000E+12
       coef(  475) =    -0.17266572750732332031E+14
       coef(  476) =     0.37033171261579853516E+13
       coef(  477) =    -0.16873681203010947266E+13
       coef(  478) =     0.93220996846419593750E+14
       coef(  479) =     0.50635188729404052734E+13
       coef(  480) =    -0.15241585476245318359E+14
       coef(  481) =     0.45496740409019640625E+14
       coef(  482) =    -0.14729714914505825000E+15
       coef(  483) =    -0.11012206151309142578E+14
       coef(  484) =     0.31600497723313742188E+14
       coef(  485) =    -0.25228397826981945038E+10
       coef(  486) =     0.12122626112866676331E+12
       coef(  487) =    -0.91615106805686676025E+11
       coef(  488) =    -0.42187700128774810791E+11
       coef(  489) =    -0.39368446461820600586E+13
       coef(  490) =    -0.15878267920855693359E+13
       coef(  491) =     0.94176097113623574219E+13
       coef(  492) =    -0.77505887895667919922E+13
       coef(  493) =    -0.18338955508197480469E+14
       coef(  494) =    -0.54623524624235439453E+13
       coef(  495) =     0.55114622236705445312E+14
       coef(  496) =     0.85980503347352050781E+13
       coef(  497) =     0.30990282640794591904E+10
       coef(  498) =    -0.15874886530102294922E+12
       coef(  499) =    -0.33980099353177435303E+12
       coef(  500) =     0.27753526261916796875E+13
       coef(  501) =     0.29219799693882788086E+13
       coef(  502) =    -0.98791020312948378906E+13
       coef(  503) =    -0.60340879108706250000E+13
       coef(  504) =    -0.81789879048605334901E+04
       coef(  505) =    -0.34129843260929747485E+07
       coef(  506) =     0.20009294725250726938E+09
       coef(  507) =    -0.54559371139860868454E+08
       coef(  508) =    -0.51396651160566917419E+11
       coef(  509) =     0.28036292941877520752E+12
       coef(  510) =    -0.31500879483997625732E+12
       coef(  511) =     0.86319624591589778662E+08
       coef(  512) =    -0.99197109547646884918E+10
       coef(  513) =     0.25508215824355285645E+12
       coef(  514) =     0.45249031350627441406E+11
       coef(  515) =    -0.44148910793536660156E+13
       coef(  516) =     0.59608463928159873047E+13
       coef(  517) =     0.17048599186024890900E+11
       coef(  518) =    -0.93610945989102685547E+12
       coef(  519) =    -0.93725872104641835938E+13
       coef(  520) =     0.27324350329180890625E+14
       coef(  521) =    -0.21292614651992468750E+14
       coef(  522) =     0.66463472495654980469E+13
       coef(  523) =     0.10657065872326871094E+14
       coef(  524) =     0.97918137543150332031E+13
       coef(  525) =    -0.54959785226982468750E+14
       coef(  526) =     0.52787862849088811874E+09
       coef(  527) =    -0.28467830172042211914E+12
       coef(  528) =    -0.12084168988305273438E+13
       coef(  529) =     0.16915674837643117188E+14
       coef(  530) =    -0.32841459062889804688E+14
       coef(  531) =     0.25143535216811656250E+14
       coef(  532) =     0.22980464716644511719E+13
       coef(  533) =    -0.80655191540423603516E+13
       coef(  534) =     0.15405047017426644531E+14
       coef(  535) =    -0.86773713165191531250E+14
       coef(  536) =    -0.46227468105049593750E+14
       coef(  537) =     0.15083994886754471875E+15
       coef(  538) =    -0.75839451339420468750E+13
       coef(  539) =     0.74152867072952171875E+14
       coef(  540) =     0.16421027230773238281E+14

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
      maxorder=10

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

