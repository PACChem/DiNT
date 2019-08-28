      subroutine cn(xx,nclu,cntot,time)

      implicit none
      include 'param.f'

      integer i,j,k,nclu
      double precision xx(3,mnat),lind,rij,cnx,del,rnn,rcut,
     &  cntot,time,rcom

      rnn = (4.0217/dsqrt(2.d0))/autoang
      rcut = rnn*1.5d0
      del = rnn/2.d0

      cntot = 0.d0
      do i=1,nclu
        cnx = 0.d0
c compute distance of atom i from center of mass (assumed to be at the origin)
c        rcom = 0.d0
c        do k=1,3
c        rcom = rcom + xx(k,i)**2
c        enddo
c        rcom = dsqrt(rcom)
c compute coordination number for atom i
        do j=1,nclu
          if (i.ne.j) then
          rij = 0.d0
          do k=1,3
          rij = rij + (xx(k,i)-xx(k,j))**2
          enddo
          rij = dsqrt(rij)
          cnx = cnx + (1.d0-dtanh(del*(rij-rcut)))/2.d0
          endif
        enddo
c     compute average CN
      cntot = cntot + cnx
      enddo
      cntot = cntot/nclu

      return
 
      end
