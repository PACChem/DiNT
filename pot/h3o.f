      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      parameter(autokcal=627.5096d0)   ! Molpro's number
      dimension r(maxatom*(maxatom-1)/2)
      character*2 symb(maxatom)
      logical lh2

      lh2=.false.
      v=0.
      do i=1,natom   ! sometimes we call with fewer than 4 atoms
      x2(i)=x(i)
      y2(i)=y(i)
      z2(i)=z(i)
      enddo
      do i=1,4
      dvdx2(i)=0.d0
      dvdy2(i)=0.d0
      dvdz2(i)=0.d0
      enddo

      nn=natom

      if (natom.eq.1) return  ! H atom

      if (natom.eq.2.and.(symb(1).eq."O".or.symb(1).eq."o")) then   ! OH
      nn=4
      x2(3)=9.
      y2(3)=0.
      z2(3)=0.73707
      x2(4)=9.
      y2(4)=0.
      z2(4)=0.
      endif

      if (natom.eq.2.and.(symb(1).eq."H".or.symb(1).eq."h")) then   ! H2
      lh2=.true.
      nn=4
      x2(1)=9.
      y2(1)=0.
      z2(1)=0.97766
      x2(2)=9.
      y2(2)=0.
      z2(2)=0.
      x2(3)=x(1)
      y2(3)=y(1)
      z2(3)=z(1)
      x2(4)=x(2)
      y2(4)=y(2)
      z2(4)=z(2)
      endif
      
      if (natom.eq.3) then
      nn=4
      x2(4)=9.
      y2(4)=0.
      z2(4)=0.
      endif

      ii=0
      do i=1,nn
      do j=i+1,nn
      ii=ii+1
      dx=x2(i)-x2(j)
      dy=y2(i)-y2(j)
      dz=z2(i)-z2(j)
      r(ii)=dsqrt(dx*dx+dy*dy+dz*dz)*autoang
      enddo
      enddo

      call pip(r,v,maxatom)
      if (v.lt.-10.) print *,"HOLE IN THE PES? V=",v
      v=v/autokcal

      resp=0.001d0

      ii=0
      do i=1,nn
      do j=i+1,nn
      ii=ii+1
      dx=(x2(i)-x2(j))*autoang
      dy=(y2(i)-y2(j))*autoang
      dz=(z2(i)-z2(j))*autoang
      r(ii)=r(ii)+resp
      call pip(r,vp,maxatom)
      r(ii)=r(ii)-2.d0*resp
      call pip(r,vm,maxatom)
      r(ii)=r(ii)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)/autokcal*autoang
      if (natom.eq.3.and.(i.eq.4.or.j.eq.4)) dtmpdrr=0.d0
      if (natom.eq.2.and.(i.gt.2.or.j.gt.2).and..not.lh2) dtmpdrr=0.d0
      if (natom.eq.2.and.(i.lt.3.or.j.lt.3).and.lh2) dtmpdrr=0.d0
      dvdx2(i) = dvdx2(i) + dtmpdrr*dx/r(ii)
      dvdx2(j) = dvdx2(j) - dtmpdrr*dx/r(ii)
      dvdy2(i) = dvdy2(i) + dtmpdrr*dy/r(ii)
      dvdy2(j) = dvdy2(j) - dtmpdrr*dy/r(ii)
      dvdz2(i) = dvdz2(i) + dtmpdrr*dz/r(ii)
      dvdz2(j) = dvdz2(j) - dtmpdrr*dz/r(ii)
      enddo
      enddo

      if (.not.lh2) then
      do i=1,natom
      dvdx(i)=dvdx2(i)
      dvdy(i)=dvdy2(i)
      dvdz(i)=dvdz2(i)
      enddo
      else
      do i=1,natom
      dvdx(i)=dvdx2(i+2)
      dvdy(i)=dvdy2(i+2)
      dvdz(i)=dvdz2(i+2)
      enddo
      endif

      end

      subroutine pip(r,v,n)

      implicit double precision(a-h,o-z)
      dimension r(n*(n-1)/2)
      parameter(mcoef=1000)
      parameter(mterm=10000)
      parameter(mpairs=10)
      dimension coef(mcoef),ind(mterm,mpairs),ibasis(mterm)
      character dum
      save npairs,nterms,ibasis,coef,ind

      v=0.d0
      do i=1,nterms
      arg=1.d0
      do j=1,npairs
      arg=arg*(dexp(-r(j))**ind(i,j))
      enddo
      v=v+arg*coef(ibasis(i))
      enddo

      return 

      entry prepot   ! read basis and coefs once and save here
      print *,'prepot'
      open(55,file="basisX.dat")
      read(55,*)natom,npairs,ncoef,nterms
      read(55,*)
      do i=1,nterms
      read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
      enddo
      close(55)
      open(55,file="coefX.dat")
      do i=1,ncoef
      read(55,*)k,coef(i)
      enddo
      close(55)
      return

      end

