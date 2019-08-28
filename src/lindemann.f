      subroutine lindemann(xx,nclu,step,time,arij,arij2,lind)

      implicit none
      include 'param.f'

      integer i,j,k,nclu
      double precision xx(3,mnat),lind,arij2(mnat,mnat),
     & arij(mnat,mnat),rij,tmp1,tmp2,step,time

      lind = 0.d0
      do i=1,nclu
      do j=i+1,nclu
        rij = 0.d0
        do k=1,3
        rij = rij + (xx(k,i)-xx(k,j))**2
        enddo
        arij(i,j) = arij(i,j) + dsqrt(rij)*step
        arij2(i,j) = arij2(i,j) + rij*step
        tmp1 = arij(i,j)  / time
        tmp2 = arij2(i,j) / time
        if (tmp2.gt.tmp1**2) lind = lind + dsqrt(tmp2-tmp1**2)/tmp1
      enddo
      enddo

      lind = lind*2.d0/dble(nclu*(nclu-1))

      return
 
      end
