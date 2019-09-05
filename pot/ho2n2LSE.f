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

       ncoef=  380
       coef(    1) =    -0.13441769310772675450E+01
       coef(    2) =    -0.50237575952314074357E+02
       coef(    3) =     0.24607395580368818628E+05
       coef(    4) =    -0.98089998333854530938E+06
       coef(    5) =     0.95705358385762795806E+07
       coef(    6) =    -0.11518299768164146371E+05
       coef(    7) =    -0.19637596671333210543E+07
       coef(    8) =     0.33512879794124297798E+08
       coef(    9) =    -0.61631727814743578434E+09
       coef(   10) =     0.41580594686354473233E+08
       coef(   11) =    -0.50480375149584794044E+09
       coef(   12) =     0.88171977224615454674E+09
       coef(   13) =     0.13446702902085533142E+10
       coef(   14) =     0.39117689056940414429E+11
       coef(   15) =    -0.13623091803255397034E+12
       coef(   16) =    -0.23861917596269184723E+07
       coef(   17) =     0.55974489739887773991E+09
       coef(   18) =    -0.39016430445658874512E+10
       coef(   19) =     0.34001650441722717285E+11
       coef(   20) =    -0.18174131648302520752E+11
       coef(   21) =     0.25374714166729226685E+12
       coef(   22) =    -0.83808507900056323242E+12
       coef(   23) =    -0.10236542913967197266E+13
       coef(   24) =     0.38699897168520185547E+13
       coef(   25) =    -0.11214749510786638184E+13
       coef(   26) =     0.82755211631453826904E+11
       coef(   27) =    -0.12958452443009724121E+13
       coef(   28) =     0.28092141899424877930E+13
       coef(   29) =     0.61268784687863320312E+13
       coef(   30) =    -0.10607491853772996094E+14
       coef(   31) =    -0.48682133942375849609E+13
       coef(   32) =    -0.36303060460264262695E+13
       coef(   33) =    -0.22354034686388500000E+14
       coef(   34) =    -0.33595676335325270891E+09
       coef(   35) =     0.17218439435529228210E+11
       coef(   36) =    -0.16623717852974566650E+12
       coef(   37) =     0.86107398944821868896E+11
       coef(   38) =    -0.19445901827595608521E+12
       coef(   39) =     0.43806114341653056641E+13
       coef(   40) =    -0.49191346745290683594E+13
       coef(   41) =    -0.14317218045811210938E+14
       coef(   42) =     0.62102375481799046875E+14
       coef(   43) =    -0.10144272717078759375E+15
       coef(   44) =    -0.30441641724236761475E+12
       coef(   45) =    -0.98301741692421503906E+13
       coef(   46) =    -0.12604771472010390625E+13
       coef(   47) =     0.79454996355444785156E+13
       coef(   48) =     0.33545682899028406250E+14
       coef(   49) =     0.89035477283680531250E+14
       coef(   50) =     0.18011885790046704102E+13
       coef(   51) =     0.26322921994537480469E+14
       coef(   52) =    -0.60552859216606156250E+14
       coef(   53) =    -0.94175469826962781250E+14
       coef(   54) =    -0.91995456221526443481E+11
       coef(   55) =    -0.18331986459138360596E+12
       coef(   56) =     0.10787642995802163696E+12
       coef(   57) =    -0.23564918263081577301E+11
       coef(   58) =     0.69034799034395736694E+11
       coef(   59) =    -0.55494222430280319214E+11
       coef(   60) =     0.42320250225190391541E+11
       coef(   61) =     0.38599544046957237244E+11
       coef(   62) =     0.96901567172672195435E+10
       coef(   63) =     0.17979779508467594147E+11
       coef(   64) =    -0.21359254003788108826E+10
       coef(   65) =     0.24221693261844573975E+12
       coef(   66) =    -0.27189794613951580811E+12
       coef(   67) =    -0.78168191813521652222E+11
       coef(   68) =     0.10846775458551979980E+13
       coef(   69) =    -0.36498498548867859840E+10
       coef(   70) =    -0.19596360783711185455E+10
       coef(   71) =     0.13577477350331845093E+12
       coef(   72) =     0.12026871181576846313E+12
       coef(   73) =    -0.75947052055313766479E+11
       coef(   74) =    -0.46541203816393134766E+13
       coef(   75) =     0.82750478295791748047E+13
       coef(   76) =     0.33935806997972019531E+14
       coef(   77) =    -0.14674527320351678125E+15
       coef(   78) =     0.15832998291478465625E+15
       coef(   79) =     0.29832202537540078125E+13
       coef(   80) =    -0.27821215460231921875E+14
       coef(   81) =     0.34258116445249703125E+14
       coef(   82) =     0.54872599881395664062E+14
       coef(   83) =    -0.44411038392128976562E+14
       coef(   84) =    -0.28475221169574136719E+14
       coef(   85) =    -0.47280343993577699661E+10
       coef(   86) =    -0.93425400464127655029E+10
       coef(   87) =    -0.20787029883217797852E+13
       coef(   88) =     0.28000448599364042969E+13
       coef(   89) =     0.16711313265589357910E+13
       coef(   90) =     0.32479113477404226562E+14
       coef(   91) =    -0.54552230153165625000E+14
       coef(   92) =    -0.13294697148198943750E+15
       coef(   93) =     0.33478290667942025000E+15
       coef(   94) =    -0.25731429443048816406E+14
       coef(   95) =     0.84404837834124468750E+14
       coef(   96) =     0.11644510648644267578E+14
       coef(   97) =    -0.99464103520580296875E+14
       coef(   98) =     0.34885984878096968750E+14
       coef(   99) =    -0.82608710951357218750E+14
       coef(  100) =     0.20136091164695566893E+09
       coef(  101) =    -0.62763787086045615375E+08
       coef(  102) =     0.37475372234113025665E+09
       coef(  103) =    -0.28763686986356263161E+10
       coef(  104) =     0.11867839699778327942E+11
       coef(  105) =    -0.30132169989650583267E+09
       coef(  106) =    -0.25440391244917593002E+10
       coef(  107) =    -0.64014989353076919556E+11
       coef(  108) =     0.24540964164603878784E+12
       coef(  109) =     0.27411007107614297485E+12
       coef(  110) =    -0.28963613987854682617E+13
       coef(  111) =     0.63740790745505578613E+12
       coef(  112) =     0.10218123199829587891E+14
       coef(  113) =    -0.14724597154068251953E+14
       coef(  114) =     0.32481359126769316406E+13
       coef(  115) =     0.11807823432342330933E+11
       coef(  116) =    -0.45625591437880761719E+12
       coef(  117) =     0.54181533852186699219E+13
       coef(  118) =    -0.12864821715334199219E+14
       coef(  119) =     0.44246017344869384766E+10
       coef(  120) =    -0.20784648605019753906E+14
       coef(  121) =     0.83890618844711968750E+14
       coef(  122) =     0.19348573198132949219E+13
       coef(  123) =    -0.92993379107278265625E+14
       coef(  124) =    -0.86711855036222265625E+12
       coef(  125) =     0.49767225002931367188E+14
       coef(  126) =    -0.74359402845422406250E+14
       coef(  127) =    -0.15047141292596248047E+14
       coef(  128) =    -0.29420172299751804352E+11
       coef(  129) =     0.15573572406842460938E+13
       coef(  130) =    -0.21601207591972078125E+14
       coef(  131) =     0.59333107526810687500E+14
       coef(  132) =    -0.10489729061716533203E+13
       coef(  133) =    -0.52990104448981421875E+14
       coef(  134) =    -0.22899322090502521875E+15
       coef(  135) =     0.10974847385351018750E+15
       coef(  136) =     0.24651543031511382812E+14
       coef(  137) =     0.54332420093458367188E+14
       coef(  138) =     0.17621472482677273438E+14
       coef(  139) =    -0.20685943432926083915E+06
       coef(  140) =     0.17244248185549184680E+08
       coef(  141) =    -0.29696180175754303932E+10
       coef(  142) =     0.43168509802950592041E+11
       coef(  143) =    -0.22805461403796969604E+12
       coef(  144) =     0.26553964826174201965E+10
       coef(  145) =     0.46867781100196502686E+11
       coef(  146) =     0.26933744889652880859E+12
       coef(  147) =     0.89889858906860729980E+12
       coef(  148) =    -0.31271914128706923828E+13
       coef(  149) =     0.32537364115702242188E+14
       coef(  150) =    -0.69374504289555869141E+13
       coef(  151) =    -0.11955517561900181250E+15
       coef(  152) =     0.15760423290286896875E+15
       coef(  153) =    -0.11128038543029281616E+12
       coef(  154) =     0.34472250482913115234E+13
       coef(  155) =    -0.47900431659684148438E+14
       coef(  156) =     0.58706620188156132812E+14
       coef(  157) =     0.22867024157250449219E+14
       coef(  158) =     0.79460727487548125000E+14
       coef(  159) =    -0.27921662341555459375E+15
       coef(  160) =     0.15921467798511193750E+15
       coef(  161) =    -0.76890790769311000000E+14
       coef(  162) =    -0.99926456436788781250E+14
       coef(  163) =     0.70949750926087243652E+12
       coef(  164) =    -0.24975733831547199219E+14
       coef(  165) =     0.13485825997516792188E+15
       coef(  166) =    -0.11022581480923396875E+15
       coef(  167) =     0.34495373723060710938E+14
       coef(  168) =     0.68035467668133031250E+14
       coef(  169) =    -0.12397025098705139062E+15
       coef(  170) =     0.63153221925071990117E+07
       coef(  171) =    -0.53731740130411076546E+09
       coef(  172) =     0.11745682059790689468E+11
       coef(  173) =    -0.14478706033070138550E+12
       coef(  174) =     0.65177591567008813477E+12
       coef(  175) =     0.62022176796288967133E+10
       coef(  176) =    -0.11532268859560745239E+12
       coef(  177) =     0.32445188404234570312E+12
       coef(  178) =    -0.76033582791867451172E+13
       coef(  179) =     0.28003604029933168945E+13
       coef(  180) =    -0.22132037639992308594E+14
       coef(  181) =    -0.41740329489709101562E+14
       coef(  182) =     0.10646808482510310938E+15
       coef(  183) =    -0.45450735740471939087E+11
       coef(  184) =    -0.21022046711928393555E+12
       coef(  185) =     0.51128466109249046875E+14
       coef(  186) =     0.26040589605806613281E+14
       coef(  187) =    -0.66722789807870203125E+14
       coef(  188) =    -0.78445058971249156250E+14
       coef(  189) =     0.21047968040086500000E+15
       coef(  190) =    -0.10813919326073872070E+13
       coef(  191) =     0.38001442295436156250E+14
       coef(  192) =    -0.20523411177890175000E+15
       coef(  193) =    -0.18808541842113109375E+14
       coef(  194) =     0.91995456219574295044E+11
       coef(  195) =     0.18331979485470101929E+12
       coef(  196) =    -0.10787664645777816772E+12
       coef(  197) =     0.23621717570360145569E+11
       coef(  198) =    -0.68544286422512298584E+11
       coef(  199) =     0.55509643262947868347E+11
       coef(  200) =    -0.43108482612986373901E+11
       coef(  201) =    -0.36194364115101310730E+11
       coef(  202) =    -0.78262135982206420898E+11
       coef(  203) =    -0.14836042310764045715E+11
       coef(  204) =     0.47877736312764358521E+11
       coef(  205) =     0.25331473213800970459E+12
       coef(  206) =    -0.30567986191299822998E+12
       coef(  207) =    -0.10633258622242489624E+12
       coef(  208) =     0.10739629869803725586E+13
       coef(  209) =     0.35853991687879872322E+10
       coef(  210) =     0.11991666318443634033E+11
       coef(  211) =     0.14462442314565411377E+12
       coef(  212) =     0.12775584737127972412E+12
       coef(  213) =    -0.10663756728846185303E+12
       coef(  214) =    -0.46437631535448115234E+13
       coef(  215) =     0.82666467583064697266E+13
       coef(  216) =     0.33919095626330535156E+14
       coef(  217) =    -0.14671351192611803125E+15
       coef(  218) =     0.15834414550924681250E+15
       coef(  219) =     0.29807416430788774414E+13
       coef(  220) =    -0.27823865057723339844E+14
       coef(  221) =     0.34263498095786789062E+14
       coef(  222) =     0.54866823089354695312E+14
       coef(  223) =    -0.44414376020976289062E+14
       coef(  224) =    -0.28479270004623226562E+14
       coef(  225) =     0.15739322369906470776E+10
       coef(  226) =    -0.77778045978924751282E+10
       coef(  227) =    -0.20846299527375520020E+13
       coef(  228) =     0.27979080370032612305E+13
       coef(  229) =     0.16772798153281809082E+13
       coef(  230) =     0.32477546843020671875E+14
       coef(  231) =    -0.54551616318693687500E+14
       coef(  232) =    -0.13294867797643850000E+15
       coef(  233) =     0.33478325307358168750E+15
       coef(  234) =    -0.25735614836975472656E+14
       coef(  235) =     0.84407491785969296875E+14
       coef(  236) =     0.11641820554263171875E+14
       coef(  237) =    -0.99466615928531750000E+14
       coef(  238) =     0.34883956046243531250E+14
       coef(  239) =    -0.82605844338554156250E+14
       coef(  240) =     0.38687178526586292719E+05
       coef(  241) =    -0.33381256503402497619E+07
       coef(  242) =    -0.23429782468096062541E+09
       coef(  243) =     0.88798990564641928673E+09
       coef(  244) =    -0.14923209630438783646E+11
       coef(  245) =     0.24871785054443985224E+09
       coef(  246) =     0.31861814391057456970E+11
       coef(  247) =    -0.10939022828599328613E+12
       coef(  248) =     0.56457384030903698730E+12
       coef(  249) =    -0.79733839191980603027E+12
       coef(  250) =     0.11160550770499710938E+14
       coef(  251) =    -0.19479393027212410156E+14
       coef(  252) =    -0.46386752319173203125E+14
       coef(  253) =     0.17526911818177909375E+15
       coef(  254) =    -0.18835794914201906250E+15
       coef(  255) =    -0.27193021267929210663E+11
       coef(  256) =     0.48070644296219543457E+12
       coef(  257) =    -0.56917011510414550781E+13
       coef(  258) =     0.40305372485425849609E+13
       coef(  259) =     0.10971703650977491455E+13
       coef(  260) =     0.67042284411024599609E+13
       coef(  261) =     0.15910572902894218750E+14
       coef(  262) =     0.18558878464698371094E+14
       coef(  263) =    -0.73170022277321875000E+13
       coef(  264) =     0.41850369697820864258E+13
       coef(  265) =    -0.41089761419933007812E+14
       coef(  266) =    -0.26381875811832320312E+14
       coef(  267) =    -0.83916220052779312500E+14
       coef(  268) =     0.25803000090327777100E+12
       coef(  269) =    -0.41831088851919282227E+13
       coef(  270) =     0.27559329647218867188E+14
       coef(  271) =    -0.42023052540250937500E+14
       coef(  272) =    -0.78867456322915439453E+13
       coef(  273) =     0.42266647493224593750E+14
       coef(  274) =    -0.97575412379729328125E+14
       coef(  275) =     0.20896274860422046875E+14
       coef(  276) =     0.24861771599140089844E+14
       coef(  277) =     0.18584640688022300781E+14
       coef(  278) =    -0.11866283866929917188E+15
       coef(  279) =    -0.12017987081542311353E+06
       coef(  280) =     0.11600143564538781345E+09
       coef(  281) =    -0.27679752040547595024E+10
       coef(  282) =     0.70250270541429779053E+11
       coef(  283) =    -0.12666610203222969055E+12
       coef(  284) =    -0.72247499508540191650E+10
       coef(  285) =     0.34417156416230819702E+11
       coef(  286) =    -0.31236111583693691406E+13
       coef(  287) =     0.66049778395263623047E+13
       coef(  288) =     0.22000586169251547852E+13
       coef(  289) =    -0.13801765516783761719E+14
       coef(  290) =    -0.38487650872521500000E+14
       coef(  291) =     0.11394373289859737500E+15
       coef(  292) =    -0.15475713566875818750E+15
       coef(  293) =     0.24026930739008856201E+12
       coef(  294) =    -0.17493530616675600586E+13
       coef(  295) =     0.42272813466795187500E+14
       coef(  296) =    -0.19027971777275765625E+14
       coef(  297) =    -0.19446214687437972656E+14
       coef(  298) =    -0.18103213937197165625E+15
       coef(  299) =     0.18406202219506353125E+15
       coef(  300) =     0.86223327388924687500E+13
       coef(  301) =     0.42491329606370562500E+14
       coef(  302) =     0.25718291242095250000E+15
       coef(  303) =    -0.35059543502374790039E+13
       coef(  304) =     0.49379384440434570312E+14
       coef(  305) =    -0.21835809879959725000E+15
       coef(  306) =     0.22058767049421950000E+15
       coef(  307) =    -0.48292388008629375000E+14
       coef(  308) =     0.17892217275724821875E+15
       coef(  309) =     0.43316207000575226562E+14
       coef(  310) =    -0.14564896870205471292E+08
       coef(  311) =     0.14912305552066645622E+10
       coef(  312) =    -0.11690790441230878830E+11
       coef(  313) =    -0.27610924382811181641E+12
       coef(  314) =     0.98386672063273181152E+12
       coef(  315) =    -0.43968623595625747681E+11
       coef(  316) =     0.15466718131118640137E+13
       coef(  317) =     0.54432367559544345703E+13
       coef(  318) =    -0.21941178234685132812E+14
       coef(  319) =    -0.94313079586048496094E+13
       coef(  320) =     0.12148786571297917969E+14
       coef(  321) =     0.14979516428692709375E+15
       coef(  322) =    -0.13546604790931809375E+15
       coef(  323) =    -0.22139013455360607910E+12
       coef(  324) =    -0.20835049810758847656E+14
       coef(  325) =     0.16083521535192458984E+14
       coef(  326) =    -0.12905745742358393750E+15
       coef(  327) =     0.16431655317806321875E+15
       coef(  328) =    -0.18555741446780355469E+14
       coef(  329) =    -0.23682728031293409375E+15
       coef(  330) =     0.13011229717310318359E+14
       coef(  331) =    -0.11436712954932409375E+15
       coef(  332) =     0.22298921718336015625E+15
       coef(  333) =     0.37358304534129820312E+14
       coef(  334) =    -0.30239169389934904873E+08
       coef(  335) =     0.18844675059867997169E+10
       coef(  336) =    -0.29704154018752819061E+11
       coef(  337) =     0.84800756281265405273E+12
       coef(  338) =    -0.29590299213457368164E+13
       coef(  339) =    -0.16985029454546087265E+11
       coef(  340) =    -0.16852416144078569336E+13
       coef(  341) =    -0.11206895102757658203E+14
       coef(  342) =     0.45633577044819078125E+14
       coef(  343) =     0.14239293009343630859E+14
       coef(  344) =    -0.21053617367207781250E+14
       coef(  345) =     0.14426451408055705566E+13
       coef(  346) =     0.15015738871204767578E+14
       coef(  347) =    -0.27882053500085097656E+14
       coef(  348) =    -0.10341654470002385938E+15
       coef(  349) =    -0.18937031076576117188E+14
       coef(  350) =     0.11580979756931146875E+15
       coef(  351) =    -0.20137460303407478333E+09
       coef(  352) =     0.62320766880137175322E+08
       coef(  353) =    -0.11431277054869353771E+08
       coef(  354) =    -0.16688850762346134186E+10
       coef(  355) =     0.11550738490590766907E+11
       coef(  356) =    -0.16472990118257844448E+09
       coef(  357) =    -0.25710786455498485565E+10
       coef(  358) =    -0.63659367195273147583E+11
       coef(  359) =     0.24453149536793521118E+12
       coef(  360) =     0.27371332287038302612E+12
       coef(  361) =    -0.28971098498020356445E+13
       coef(  362) =     0.63750463703386267090E+12
       coef(  363) =     0.10217737365239253906E+14
       coef(  364) =    -0.14724511725062001953E+14
       coef(  365) =     0.32475478017267812500E+13
       coef(  366) =     0.11933554908263921738E+11
       coef(  367) =    -0.45630473886651757812E+12
       coef(  368) =     0.54186804108844843750E+13
       coef(  369) =    -0.12864644412294660156E+14
       coef(  370) =     0.42317146694184570312E+10
       coef(  371) =    -0.20784252248650171875E+14
       coef(  372) =     0.83891062364290703125E+14
       coef(  373) =     0.19349631285385195312E+13
       coef(  374) =    -0.92993179201721968750E+14
       coef(  375) =    -0.86696713240406237793E+12
       coef(  376) =     0.49767256738549671875E+14
       coef(  377) =    -0.74359578217326406250E+14
       coef(  378) =    -0.15046979492951693359E+14
       coef(  379) =    -0.29423527490438888550E+11
       coef(  380) =     0.15569763938465585938E+13


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

      norder=4
      norder2=4
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

