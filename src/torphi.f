      subroutine torphi(xx,nclu,a1,a2,phi,n1)

      implicit none
      include 'param.f'

      integer nclu,a1,a2,n1
      double precision xx(3,mnat),phi

c local
      integer i,j
      double precision axis(3),xx2(3,mnat),cp,op,sp,tmp,rot(3,3),
     & tmp1,tmp2,tmp3

c move to put one of the atoms at the origin
      do i=1,nclu
      do j=1,3
        xx2(j,i)=xx(j,i)-xx(j,a2)
      enddo
      enddo

c rotation axis
      tmp=0.d0
      do i=1,3
        axis(i)=xx2(i,a1)-xx2(i,a2)
        tmp=tmp+axis(i)**2
      enddo
      do i=1,3
      axis(i)=axis(i)/dsqrt(tmp)
      enddo

c rotation matrix
c      cp=dcos(phi/180.d0*pi)
c      sp=dsin(phi/180.d0*pi)
      cp=dcos(phi)
      sp=dsin(phi)
      op=1.d0-cp

      rot(1,1)=cp+axis(1)**2*op
      rot(2,2)=cp+axis(2)**2*op
      rot(3,3)=cp+axis(3)**2*op
      rot(1,2)=axis(1)*axis(2)*op-axis(3)*sp
      rot(2,1)=axis(1)*axis(2)*op+axis(3)*sp
      rot(1,3)=axis(1)*axis(3)*op+axis(2)*sp
      rot(3,1)=axis(1)*axis(3)*op-axis(2)*sp
      rot(2,3)=axis(2)*axis(3)*op-axis(1)*sp
      rot(3,2)=axis(2)*axis(3)*op+axis(1)*sp

c rotate
      do i=1,n1
         tmp1 = xx2(1,i)
         tmp2 = xx2(2,i)
         tmp3 = xx2(3,i)
         xx2(1,i)=tmp1*rot(1,1)+tmp2*rot(2,1)+tmp3*rot(3,1)
         xx2(2,i)=tmp1*rot(1,2)+tmp2*rot(2,2)+tmp3*rot(3,2)
         xx2(3,i)=tmp1*rot(1,3)+tmp2*rot(2,3)+tmp3*rot(3,3)
      enddo

c move back 
      do i=1,nclu
      do j=1,3
        xx2(j,i)=xx2(j,i)+xx(j,a2)
      enddo
      enddo
      do j=1,3
      do i=1,nclu
        xx(j,i)=xx2(j,i)
      enddo
      enddo

      return

      end
