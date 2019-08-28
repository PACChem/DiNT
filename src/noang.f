      subroutine noang(xx,pp,mm,natom)

c Remove the component of PP that corresponds to overall rotation.

c see e.g., W. L. Hase, D. G. Buckowski, and K. N. Swamy, 
c           J. Phys. Chem., vol. 87, p. 2754, 1983.
c This subroutine was inspired by a subroutine from GenDyne.

      implicit none
      include 'param.f'

      integer natom
      double precision xx(3,mnat),pp(3,mnat),mm(mnat),
     &   bigj(3),bigjtot,xxp(3,mnat),ppp(3,mnat)

c local
      integer i,j,info
      double precision d1,d2,d3,d4,d5,d6,det,omeg(3),
     & c11,c12,c13,c21,c22,c23,c31,c32,c33,temp1,temp2,
     & mom(3,3),ap(6),eig(3),rot(3,3),work(9),temp3


c diatoms
c project out momenta ortogonal to the bond axis
c not debugged
      if (natom.eq.2) then
        temp1=0.d0
        temp2=0.d0
        temp3=0.d0
        do i=1,3
           temp3=temp3+(xx(i,1)-xx(i,2))**2
        enddo
        temp3=dsqrt(temp3)
        do i=1,3
           temp1=temp1+pp(i,1)*(xx(i,1)-xx(i,2))/temp3
           temp2=temp2+pp(i,2)*(xx(i,1)-xx(i,2))/temp3
        enddo
        do i=1,3
           pp(i,1)=temp1*(xx(i,1)-xx(i,2))/temp3
           pp(i,2)=temp2*(xx(i,1)-xx(i,2))/temp3
        enddo
       return
       endif

c compute moment of intertia matrix mom
      do i=1,3
      do j=1,3
          mom(i,j) = 0.d0
      enddo
      enddo

      do i=1,natom
         mom(1,1)=mom(1,1)+mm(i)*(xx(2,i)**2+xx(3,i)**2)
         mom(2,2)=mom(2,2)+mm(i)*(xx(1,i)**2+xx(3,i)**2)
         mom(3,3)=mom(3,3)+mm(i)*(xx(1,i)**2+xx(2,i)**2)
         mom(1,2)=mom(1,2)-mm(i)*(xx(1,i)*xx(2,i))
         mom(1,3)=mom(1,3)-mm(i)*(xx(1,i)*xx(3,i))
         mom(2,3)=mom(2,3)-mm(i)*(xx(2,i)*xx(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

c      print *,"mom"
c      print *,(mom(1,j),j=1,3)
c      print *,(mom(2,j),j=1,3)
c      print *,(mom(3,j),j=1,3)

c      call angmom(xx,pp,mm,natom,bigj,bigjtot)
c      print *,"j",bigj,bigjtot

c     diagonalize the mom matrix
      do i=1,3
      do j=i,3
        ap(i+(j-1)*j/2)=mom(i,j)
      enddo
      enddo
      call dspev( 'v','u',3,ap,eig,rot,3,work,info )

c      print *,"rot"
c      print *,(rot(1,j),j=1,3)
c      print *,(rot(2,j),j=1,3)
c      print *,(rot(3,j),j=1,3)
c      print *,"eig"
c      print *,(eig(j),j=1,3)

c rotate to diagonalize mom
      do i=1,natom
         temp1 = xx(1,i)
         temp2 = xx(2,i)
         temp3 = xx(3,i)
         xx(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         xx(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         xx(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
         temp1 = pp(1,i)
         temp2 = pp(2,i)
         temp3 = pp(3,i)
         pp(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         pp(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         pp(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
      enddo

      call angmom(xx,pp,mm,natom,bigj,bigjtot)
c      print *,bigj,bigjtot

c omega = inverse of moment of intertia matrix dotted into total angular momentum vector
c omega = I-1 . J
c I (mom) is diagonal now, so the inverse of I is the diagonal matrix of the inverses of the eigenvalues
      do i=1,3
          omeg(i) = bigj(i)/eig(i)
      enddo

c correction = r x omega
      do i=1,natom
         pp(1,i)=pp(1,i)-mm(i)*(omeg(2)*xx(3,i)-omeg(3)*xx(2,i))
         pp(2,i)=pp(2,i)-mm(i)*(omeg(3)*xx(1,i)-omeg(1)*xx(3,i))
         pp(3,i)=pp(3,i)-mm(i)*(omeg(1)*xx(2,i)-omeg(2)*xx(1,i))
      enddo

c invert the rot matrix
      d1 = rot(1,1)*rot(2,2)*rot(3,3)
      d2 = rot(1,2)*rot(2,3)*rot(3,1)
      d3 = rot(1,3)*rot(2,1)*rot(3,2)
      d4 = rot(1,3)*rot(2,2)*rot(3,1)
      d5 = rot(1,2)*rot(2,1)*rot(3,3)
      d6 = rot(1,1)*rot(2,3)*rot(3,2)
      det = d1 + d2 + d3 - d4 - d5 - d6
      c11 = rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2)
      c12 = rot(1,3)*rot(3,2)-rot(1,2)*rot(3,3)
      c13 = rot(1,2)*rot(2,3)-rot(1,3)*rot(2,2)
      c21 = rot(2,3)*rot(3,1)-rot(2,1)*rot(3,3)
      c22 = rot(1,1)*rot(3,3)-rot(1,3)*rot(3,1)
      c23 = rot(1,3)*rot(2,1)-rot(1,1)*rot(2,3)
      c31 = rot(2,1)*rot(3,2)-rot(2,2)*rot(3,1)
      c32 = rot(1,2)*rot(3,1)-rot(1,1)*rot(3,2)
      c33 = rot(1,1)*rot(2,2)-rot(1,2)*rot(2,1)
      rot(1,1) = c11/det
      rot(1,2) = c12/det
      rot(1,3) = c13/det
      rot(2,1) = c21/det
      rot(2,2) = c22/det
      rot(2,3) = c23/det
      rot(3,1) = c31/det
      rot(3,2) = c32/det
      rot(3,3) = c33/det

c rotate back to original coordinate system
      do i=1,natom
         temp1 = xx(1,i)
         temp2 = xx(2,i)
         temp3 = xx(3,i)
         xx(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         xx(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         xx(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
         temp1 = pp(1,i)
         temp2 = pp(2,i)
         temp3 = pp(3,i)
         pp(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         pp(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         pp(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
      enddo

      return
 
      end
