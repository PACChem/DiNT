      subroutine pot(x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c Lit HalP
      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
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

      rep=12.d0
     
      sig=3.68d0/autoang   ! CH4 + N2
      de=88.7d0/autocmi
      de=de/((6./rep)**(6./(rep-6.))-(6./rep)**(rep/(rep-6.)))

      do 1 i=1,natom
      do 2 j=i+1,natom

      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

      ven=(sig/rr)**6
      vem=(sig/rr)**int(rep)
      v1=v1+vem-ven

c compute deriv of potential WRT rr
      dvdr = 6.d0*ven/rr-rep*vem/rr

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + de*dvdr*dx/rr
      dvdx(j) = dvdx(j) - de*dvdr*dx/rr
      dvdy(i) = dvdy(i) + de*dvdr*dy/rr
      dvdy(j) = dvdy(j) - de*dvdr*dy/rr
      dvdz(i) = dvdz(i) + de*dvdr*dz/rr
      dvdz(j) = dvdz(j) - de*dvdr*dz/rr

    2 continue
    1 continue
      v=de*v1
      return
      end

      subroutine prepot
      return
      end
