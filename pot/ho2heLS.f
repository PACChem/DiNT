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
c        call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
        call lsbath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)
      endif
      if (natom2.ne.0) then
        call ho2(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom2,maxatom)
      endif

      tmpprint(1)=v

c      print *,v,v2,v3

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

      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      r(1)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oo
      r(1)=r(1)*autoang

      dx=x(2)-x(3)
      dy=y(2)-y(3)
      dz=z(2)-z(3)
      r(2)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oh1
      r(2)=r(2)*autoang

      dx=x(3)-x(1)
      dy=y(3)-y(1)
      dz=z(3)-z(1)
      r(3)=dsqrt(dx*dx+dy*dy+dz*dz)   ! oh2
      r(3)=r(3)*autoang

      call ho2fit(r,v)
      v=v/autocmi

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

      dimension r(3)

      do i=1,3
      r(i)=dsqrt((x(4-i)-x(4))**2+(y(4-i)-y(4))**2+(z(4-i)-z(4))**2)
      enddo

      call lsbathfit(r,v)

      resp=0.0001d0

      do i1=1,3
      if (i1.eq.1) i=3
      if (i1.eq.1) j=4
      if (i1.eq.2) i=2
      if (i1.eq.2) j=4
      if (i1.eq.3) i=1
      if (i1.eq.3) j=4
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i1)=r(i1)+resp
      call lsbathfit(r,vp)
      r(i1)=r(i1)-2.d0*resp
      call lsbathfit(r,vm)
      r(i1)=r(i1)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)
      dvdx(i) = dvdx(i) + dtmpdrr*dx/r(i1)
      dvdx(j) = dvdx(j) - dtmpdrr*dx/r(i1)
      dvdy(i) = dvdy(i) + dtmpdrr*dy/r(i1)
      dvdy(j) = dvdy(j) - dtmpdrr*dy/r(i1)
      dvdz(i) = dvdz(i) + dtmpdrr*dz/r(i1)
      dvdz(j) = dvdz(j) - dtmpdrr*dz/r(i1)
      enddo

      return
      end

      subroutine lsbathfit(r,v)

      implicit real*8(a-h,o-z)

      parameter(autoang=0.529177249d0)
      parameter(autocmi=219474.63067d0)

      dimension r(3)
      dimension coef(100),basis(100)

       ncoef=64   
       coef(    1) =    -0.464790203479406693887200000000D+00
       coef(    2) =     0.538926991220206659249920000000D+02
       coef(    3) =     0.652598414484149452619020000000D+04
       coef(    4) =    -0.834047773859128239564600000000D+06
       coef(    5) =     0.538913301571024945246790000000D+02
       coef(    6) =    -0.186110236584855447290470000000D+06
       coef(    7) =     0.559049876820417772978540000000D+07
       coef(    8) =     0.141676615031986888498070000000D+08
       coef(    9) =     0.652611566156440403574380000000D+04
       coef(   10) =     0.559048986467811465263370000000D+07
       coef(   11) =    -0.261916160114278435707090000000D+09
       coef(   12) =     0.126659279387620401382450000000D+10
       coef(   13) =    -0.834049736806231550872330000000D+06
       coef(   14) =     0.141678083196932636201380000000D+08
       coef(   15) =     0.126659083017548465728760000000D+10
       coef(   16) =    -0.890799947490126419067380000000D+10
       coef(   17) =    -0.267062116508042777240920000000D+03
       coef(   18) =     0.579551385265176795655860000000D+05
       coef(   19) =    -0.258316367364710196852680000000D+07
       coef(   20) =     0.435776294326732903718950000000D+08
       coef(   21) =     0.579550980521201272495090000000D+05
       coef(   22) =    -0.173766796102870628237720000000D+07
       coef(   23) =     0.318727649408190250396730000000D+08
       coef(   24) =    -0.834445874960983514785770000000D+09
       coef(   25) =    -0.258316235526929330080750000000D+07
       coef(   26) =     0.318728541806325465440750000000D+08
       coef(   27) =     0.173441500306848430633540000000D+10
       coef(   28) =    -0.125231725066485424041750000000D+11
       coef(   29) =     0.435776268485115766525270000000D+08
       coef(   30) =    -0.834447785969969868659970000000D+09
       coef(   31) =    -0.125231372328817520141600000000D+11
       coef(   32) =     0.139865097066250915527340000000D+12
       coef(   33) =    -0.229839068864324726746420000000D+05
       coef(   34) =     0.181085303980568336555730000000D+06
       coef(   35) =     0.154360450590725056827070000000D+08
       coef(   36) =    -0.524905102464531064033510000000D+09
       coef(   37) =     0.181085312773957703029740000000D+06
       coef(   38) =    -0.241355135309342965483670000000D+08
       coef(   39) =     0.440896135927837714552880000000D+08
       coef(   40) =     0.139220765105024585723880000000D+11
       coef(   41) =     0.154360007336929552257060000000D+08
       coef(   42) =     0.440905896914432421326640000000D+08
       coef(   43) =    -0.166215869643643531799320000000D+11
       coef(   44) =     0.255943912870476455688480000000D+11
       coef(   45) =    -0.524904493380205273628230000000D+09
       coef(   46) =     0.139220639862572097778320000000D+11
       coef(   47) =     0.255943944993595809936520000000D+11
       coef(   48) =    -0.649401774458572265625000000000D+12
       coef(   49) =     0.122298332474947441369300000000D+06
       coef(   50) =     0.241337962951077846810220000000D+07
       coef(   51) =    -0.277955712606859542429450000000D+08
       coef(   52) =     0.105647132056908738613130000000D+10
       coef(   53) =     0.241339773791181715205310000000D+07
       coef(   54) =    -0.419786168490824937820430000000D+09
       coef(   55) =     0.806425750421033763885500000000D+10
       coef(   56) =    -0.726945872434805450439450000000D+11
       coef(   57) =    -0.277958757646665535867210000000D+08
       coef(   58) =     0.806426278896270942687990000000D+10
       coef(   59) =    -0.155346462285521209716800000000D+12
       coef(   60) =     0.102540487817228417968750000000D+13
       coef(   61) =     0.105647167515809583663940000000D+10
       coef(   62) =    -0.726945917804177856445310000000D+11
       coef(   63) =     0.102540487815325610351560000000D+13
       coef(   64) =    -0.475664836059329785156250000000D+13

      r1=r(1)*autoang
      r2=r(2)*autoang
      r3=r(3)*autoang

      b1=1.d0
      b2=1.d0

      e1=dexp(-r1/b1)
      e2=dexp(-r2/b2)
      e3=dexp(-r3/b2)

      norder=3
      norder2=3

      ii=0
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

      v=0.d0
      do j=1,ncoef
        v=v+coef(j)*basis(j)
      enddo
      v=v/autocmi

      return

      end
