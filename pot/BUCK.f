      subroutine pot(x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c Lit HalP
      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom),c(3)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      parameter(autocmi=219474.63067d0)

c ADD THIS FOR ANT
      v1=0.d0
      do i=1,natom
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo
c END

      c(1)=8.038880581
      c(2)=0.289971618
      c(3)=9.735404523
      a1=10.**c(1)
      a2=c(2)
      a3=c(3)

      do 1 i=1,natom
      do 2 j=i+1,natom

      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)*autoang

      v1 = a1*exp(-rr/a2)-(a3/rr)**6
      v1 = v1/autocmi

c compute deriv of potential WRT rr
      dvdr = -a1*exp(-rr/a2)/a2+6.d0*a3**6/rr**7

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + dvdr*dx/rr/autocmi*autoang**2
      dvdx(j) = dvdx(j) - dvdr*dx/rr/autocmi*autoang**2
      dvdy(i) = dvdy(i) + dvdr*dy/rr/autocmi*autoang**2
      dvdy(j) = dvdy(j) - dvdr*dy/rr/autocmi*autoang**2
      dvdz(i) = dvdz(i) + dvdr*dz/rr/autocmi*autoang**2
      dvdz(j) = dvdz(j) - dvdr*dz/rr/autocmi*autoang**2

    2 continue
    1 continue
      v=v1
      return
      end

      subroutine prepot
      return
      end
