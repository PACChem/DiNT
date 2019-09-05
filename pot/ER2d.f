             subroutine pot(x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

C   System:                        Al
C   Functional form:               Extended Rydberg
C   Common name:                   ER2
C   Number of derivatives:         0
C   Number of bodies:              variable
C   Number of electronic surfaces: 1
C   Interface:                     HO-MM-0
C
c   Notes:  Many-body aluminum potential energy function.  The parameters
c           for this PEF have been parameterized against a data set of
c           aluminum cluster energies and bulk data.
c           This PEF has a cutoff at 6.5 angstroms.
c
c  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c              and D. G. Truhlar, "Analytic Potentials Energy Functions
c              for Aluminum Clusters," Phys. Rev. B, accepted (2004).
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c                    of the system
c  natom   --- number of atoms in the system
c  maxatom --- maximal number of atoms

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),xx(3,maxatom)
      dimension dvdy(maxatom)
      dimension dvdz(maxatom)
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

        de= 1.71013678553981441/autoev
        re= 5.08182706399609163
        a1= 1.24074255007327805
        a2= 0.551880801172447422
        a3= 0.129970688812896917
        au2= 0.143243771372740580
        bu2= 6.50000000000000000/autoang

c potential, includes cutoff
      v=0.d0
      do 1 i=1,natom
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rho = rr-re
      poly=1.d0+a1*rho+a2*rho**2+a3*rho**3
      ee=dexp(-a1*rho)
      co=0.d0
      if (rr.le.bu2) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
      tmp=de*ee*poly*co
      v=v-tmp

      dcodrr=0.d0
      if (rr.le.bu2) dcodrr=(-co/bu2)*au2/(1.d0-rr/bu2)**2
      dpolydrr=a1+2.d0*a2*rho+3.d0*a3*rho**2
      deedrr=-a1*ee
      dtmpdrr=-de*(ee*poly*dcodrr+ee*dpolydrr*co+deedrr*poly*co)
      dvdx(i) = dvdx(i) + dtmpdrr*dx/rr
      dvdx(j) = dvdx(j) - dtmpdrr*dx/rr
      dvdy(i) = dvdy(i) + dtmpdrr*dy/rr
      dvdy(j) = dvdy(j) - dtmpdrr*dy/rr
      dvdz(i) = dvdz(i) + dtmpdrr*dz/rr
      dvdz(j) = dvdz(j) - dtmpdrr*dz/rr
    1 continue

      return

      end
