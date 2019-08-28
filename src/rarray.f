      subroutine rarray(ithistraj,time,xx,nat)

      implicit none
      include 'param.f'

      integer i,j,k,nat,l,ithistraj
      double precision rr(mnat*(mnat-1)/2),rij,xx(3,mnat),time,
     & symbol(mnat)

      l=0
      do 10 i=1,nat
      do 10 j=i+1,nat
      rij=0.d0
      do k=1,3
      rij = rij + (xx(k,i)-xx(k,j))**2
      enddo
      rij = dsqrt(rij)
      l=l+1
      rr(l)=rij
 10   continue

      write(43,143)ithistraj,time*autofs,(rr(k)*autoang,k=1,l)
 143  format(i10,f20.8,100(f20.8))

      return
      end
