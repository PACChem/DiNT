      subroutine getwellindex(wellindex1,wellindex2,xwell)

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'

      integer nlight,nheavy,ilight(mnat),iheavy(mnat),
     &        ncount(mnat),mcount(mnat),wellindex1,wellindex2
      integer k,k1,k2,kk1,kk2,kp,j,i,k3
      double precision rrr,rmin,r12,r13,scale,xwell
      character*20 species
      character*3 ctmp(3)

! wellindex
      wellindex1 = 0
      wellindex2 = 0

! count and index heavies and lights
      nheavy = 0
      nlight = 0
      do k=1,nat
      ncount(k) = 0
      mcount(k) = 0
        if(symbol(k).eq.'H' .or.symbol(k).eq.' H'.or.
     &     symbol(k).eq.'h' .or.symbol(k).eq.' h'.or.
     &     symbol(k).eq.' h'.or.symbol(k).eq.' H') then
          nlight = nlight + 1
          ilight(nlight) = k
        else
          nheavy = nheavy + 1
          iheavy(nheavy) = k
        endif
      enddo

      if (nheavy.eq.0) return

! assign lights to heavies
      do kk1=1,nlight
      k1 = ilight(kk1)
      rmin = 1.d20
      do kk2=1,nheavy
      k2 = iheavy(kk2)
        rrr = 0.d0
        do j=1,3
           rrr = rrr + (xx(j,k1)-xx(j,k2))**2
        enddo
        rrr=dsqrt(rrr)
        if (rrr.lt.rmin) then
          kp = kk2
          rmin = rrr
        endif
      enddo
      ncount(kp) = ncount(kp) + 1
      if (rmin.lt.(1.6d0/autoang)) mcount(kp) = mcount(kp) + 1
      enddo

      k2 = iheavy(2)   ! first C
      k1 = iheavy(1)   ! middle C
      k3 = iheavy(3)   ! O
      r12 = 0.d0
      r13 = 0.d0
      do j=1,3
         r12 = r12 + (xx(j,k1)-xx(j,k2))**2
         r13 = r13 + (xx(j,k1)-xx(j,k3))**2
      enddo
      r12=dsqrt(r12)
      r13=dsqrt(r13)
      scale = 1.d0
c      print *,"r",r12*autoang,r13*autoang
      if (r12.gt.(1.9d0/autoang)) scale = 100.d0
      if (r13.gt.(1.9d0/autoang)) scale = 10.d0
      if (r13.gt.(1.9d0/autoang).and.r12.gt.(1.9d0/autoang)) scale = -1.d0
      wellindex1 = 100*ncount(2)+10*ncount(1)+ncount(3)
      wellindex2 = 100*mcount(2)+10*mcount(1)+mcount(3)
      xwell = dble(wellindex2)/scale
c      if (xwell.gt.100.d0) then
c          write(6,'(i5)')wellindex2
c      else if (xwell.lt.0.d0) then
c          write(6,'(i5)')-wellindex2
c      else if (xwell.gt.9.999d0) then
c          write(6,'(f5.1)')xwell
c      else if (xwell.gt..9999d0) then
c          write(6,'(f5.2)')xwell
c      else
c          write(6,*)xwell
c      endif

c      print *,wellindex1,wellindex2,xwell

c      do i=1,3
c      ctmp(i) = "Hx"
c      if (mcount(i).eq.0) ctmp(i) = "H0"
c      if (mcount(i).eq.1) ctmp(i) = "H1"
c      if (mcount(i).eq.2) ctmp(i) = "H2"
c      if (mcount(i).eq.3) ctmp(i) = "H3"
c      if (mcount(i).eq.4) ctmp(i) = "H4"
c      enddo
c      species = symbol(2) // ctmp(2) // symbol(1) // ctmp(1)
c     & // symbol(3) // ctmp(3)
c      write(6,*)species

      return

      end

