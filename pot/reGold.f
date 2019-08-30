             subroutine pot(x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c ANT: update info
c   System:                        Al
c   Functional form:               Embedded atom
c   Common name:                   reGold
c   Number of derivatives:         1
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-1
c
c   Notes:  Many-body aluminum potential energy function.  The parameters 
c           for this PEF have been parameterized against a data set of
c           aluminum cluster energies and bulk data.  This routine includes
c           derivatives.
c
c  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c              and D. G. Truhlar, "Analytic Potentials Energy Functions
c              for Aluminum Clusters," Phys. Rev. B, accepted (2004).
c              H. Gollisch, Surface Science 166 (1986), 87.
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system (output)
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system (input)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximal number of atoms (input)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension d_u(3,maxatom),d_ud(3,maxatom)    ! ADD OTHER VARS AS NEEDED
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

c ADD THIS FOR ANT
      v1=0.d0
      do i=1,natom
      dvdx(i)=0.d0
      dvdy(i)=0.d0
      dvdz(i)=0.d0
      enddo
c END

c reparameterized
      aa = 641.335d0  / autoev
      bb = 12.2192d0  / autoev
      a  = 2.9913d0  * autoang
      b  = 1.364d0  * autoang
      en = 0.66258d0

      v=0.d0
      do i=1,natom
c initialize variables need for each i
      u=0.d0
      ud=0.d0    
      do k=1,3
      do l=1,natom
      d_u(k,l)=0.d0
      d_ud(k,l)=0.d0
      enddo
      enddo

      do j=1,natom
      if(i.ne.j) then
       dx=x(i)-x(j)
       dy=y(i)-y(j)
       dz=z(i)-z(j)
       rr=dsqrt(dx*dx+dy*dy+dz*dz)

       u=u+dexp(-rr*a)
       ud=ud+dexp(-rr*b)

c ANT:  deriv of parts of PEF WRT rr
       dudr = -a*dexp(-rr*a)
       duddr = -b*dexp(-rr*b)

c ANT:  derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      d_u(1,i) = d_u(1,i) + dudr*dx/rr
      d_u(1,j) = d_u(1,j) - dudr*dx/rr
      d_u(2,i) = d_u(2,i) + dudr*dy/rr
      d_u(2,j) = d_u(2,j) - dudr*dy/rr
      d_u(3,i) = d_u(3,i) + dudr*dz/rr
      d_u(3,j) = d_u(3,j) - dudr*dz/rr

      d_ud(1,i) = d_ud(1,i) + duddr*dx/rr
      d_ud(1,j) = d_ud(1,j) - duddr*dx/rr
      d_ud(2,i) = d_ud(2,i) + duddr*dy/rr
      d_ud(2,j) = d_ud(2,j) - duddr*dy/rr
      d_ud(3,i) = d_ud(3,i) + duddr*dz/rr
      d_ud(3,j) = d_ud(3,j) - duddr*dz/rr

      endif
      enddo

      v=v+aa*u-bb*ud**en
      do l=1,natom
      dvdx(l) = dvdx(l) + aa*d_u(1,l) - bb*d_ud(1,l)*en*ud**(en-1.d0)
      dvdy(l) = dvdy(l) + aa*d_u(2,l) - bb*d_ud(2,l)*en*ud**(en-1.d0)
      dvdz(l) = dvdz(l) + aa*d_u(3,l) - bb*d_ud(3,l)*en*ud**(en-1.d0)
      enddo

      enddo
      return
      end
