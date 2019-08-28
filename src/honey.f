      subroutine honey(itmp,xx,nclu,time)
c see Honeycut & Anderson J. Phys. Chem. 91, 4950 (1987).

      implicit none
      include 'param.f'
      include 'c_output.f'

      integer nmax
      parameter(nmax=10)
      integer i,j,k,ii,nclu,ind1,ind2,ind3,isav(nmax),
     & isav2(nclu),k1,k2,nhit,itmp
      double precision xx(3,mnat),rij,rjk,rik,rkk,rcut,time,
     & haparam(0:nmax,0:nmax,0:nmax),tmph(10)

c cutoff for first shell
      rcut = 3.5d0/autoang
      rcut = rcut**2

      nhit = 0
      do i=0,nmax
      do j=0,nmax
      do k=0,nmax
        haparam(i,j,k) = 0.d0
      enddo
      enddo
      enddo
      do i=1,nclu
        do j=i+1,nclu
          ind1 = 0
          ind2 = 0
          ind3 = 0
          do ii=1,nmax
            isav(ii)=0
          enddo
          do ii=1,nclu
            isav2(ii)=0
          enddo
          rij = 0.d0
          do ii=1,3
          rij = rij + (xx(ii,i)-xx(ii,j))**2
          enddo
          if (rij.le.rcut) then
            do 10 k=1,nclu
            if (k.eq.i.or.k.eq.j) go to 10
              rik = 0.d0
              rjk = 0.d0
              do ii=1,3
                rik = rik + (xx(ii,i)-xx(ii,k))**2
                rjk = rjk + (xx(ii,j)-xx(ii,k))**2
              enddo
              if (rik.le.rcut.and.rjk.le.rcut.and.ind1.le.nmax) then
                ind1 = ind1 + 1
                isav(ind1)=k
              endif
  10        continue
            do k1=1,ind1
            do k2=k1+1,ind1
              rkk = 0.d0
              do ii=1,3
                rkk = rkk + (xx(ii,isav(k1))-xx(ii,isav(k2)))**2
              enddo
              if (rkk.le.rcut) then
                ind2 = ind2 + 1
                isav2(k1) = isav2(k1) + 1
                isav2(k2) = isav2(k2) + 1
              endif
            enddo
            enddo
            ind3 = 1
            do k1=1,ind1
              if (isav2(k1).ge.2) ind3 = 2
            enddo
            nhit = nhit + 1
            if (ind1.ge.0.and.ind2.ge.0.and.ind3.ge.0) then
            if (ind1.le.nmax.and.ind2.le.nmax.and.ind3.le.nmax) then
            haparam(ind1,ind2,ind3) = haparam(ind1,ind2,ind3) + 1.d0
            endif
            endif
          endif
        enddo
      enddo

c      do i=0,nmax
c      do j=0,nmax
c      do k=0,nmax
c        if (haparam(i,j,k).gt.0.d0) write(45,100)i,j,k,
c     & haparam(i,j,k)/dble(nhit)
c      enddo
c      enddo
c      enddo
c      write(45,*)
c 100  format(3i5,f15.5)

       tmph(1) = 0.d0
       tmph(2) = 0.d0
       tmph(3) = 0.d0
       tmph(4) = 0.d0
       tmph(5) = 0.d0
       do i=0,nmax
       tmph(1) = tmph(1) + haparam(2,i,1) + haparam(2,i,2)
       tmph(2) = tmph(2) + haparam(3,i,1) + haparam(3,i,2)
       tmph(3) = tmph(3) + haparam(4,i,1) + haparam(4,i,2)
       tmph(4) = tmph(4) + haparam(5,i,1) + haparam(5,i,2)
       tmph(5) = tmph(5) + haparam(6,i,1) + haparam(6,i,2)
       enddo
       tmph(6) = haparam(4,2,1)
       tmph(7) = haparam(4,2,2)

       if (lwrite(42)) write(42,101)itmp,time*autofs,
     & (tmph(k)/dble(nhit),k=1,7)
 101   format(i7,100f10.3)

      return
 
      end
