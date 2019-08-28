      subroutine addrot(pp,xx,mmm,nclu,jj,kk,brot1,brot2,ek2,erot,dum)

      implicit none
      include 'param.f'
c      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu
      double precision pp(3,mnat),mmm(mnat),xx(3,mnat)
      double precision jj,kk,brot1,brot2,ek2

c local
      integer i,ix,iy,iz,j,info,iunique
      double precision r1,r2,x1,x2,w,y1,y2,tmp,mom(3,3),ap(6),
     & eig(3),rot(3,3),work(9),temp1,temp2,temp3,omeg(3),
     & d1,d2,d3,d4,d5,d6,det,c11,c12,c13,c21,c22,c23,c31,c32,c33,
     & theta,ia,ib,ic,tmp12,tmp23,tmp13,tmpx,tmpy,tmpz,
     & jx,jy,jz,jx0,jy0,jz0,erot,erot0,rtmp,ptmp,dum(3)

       integer ii,ij
       double precision phi,eta,snt,cst,snp,csp,sne,cse,
     &  tmp2(9),rr,xji,ppo,mi,mj

c compute moment of intertia matrix mom
      do i=1,3
      do j=1,3
          mom(i,j) = 0.d0
      enddo
      enddo

      do i=1,nclu
         mom(1,1)=mom(1,1)+mmm(i)*(xx(2,i)**2+xx(3,i)**2)
         mom(2,2)=mom(2,2)+mmm(i)*(xx(1,i)**2+xx(3,i)**2)
         mom(3,3)=mom(3,3)+mmm(i)*(xx(1,i)**2+xx(2,i)**2)
         mom(1,2)=mom(1,2)-mmm(i)*(xx(1,i)*xx(2,i))
         mom(1,3)=mom(1,3)-mmm(i)*(xx(1,i)*xx(3,i))
         mom(2,3)=mom(2,3)-mmm(i)*(xx(2,i)*xx(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

      do i=1,3
      do j=i,3
        ap(i+(j-1)*j/2)=mom(i,j)
      enddo
      enddo
      call dspev( 'v','u',3,ap,eig,rot,3,work,info )

      do i=1,3
      if (eig(i).le.0.d0) eig(i)=1.d-16
      enddo
      ia=0.5d0/eig(1)*autocmi
      ib=0.5d0/eig(2)*autocmi
      ic=0.5d0/eig(3)*autocmi

      write(6,*)"Calculated rotational constants for this structure: "
      write(6,*)ia,ib,ic," cm-1"
      dum(1)=ia
      dum(2)=ib
      dum(3)=ic

      if (brot2.lt.0.d0.and.nclu.ne.2) then
      else
      tmp12=dabs(ia-ib)
      tmp13=dabs(ia-ic)
      tmp23=dabs(ib-ic)
      if (tmp12.lt.tmp23.and.tmp12.lt.tmp13) then
         write(6,*)"Symmetrized to ",0.5d0*(ia+ib)," (x2) and "
     &                ,ic," cm-1"
         iunique = 3
      elseif (tmp23.lt.tmp12.and.tmp23.lt.tmp13) then
         write(6,*)"Symmetrized to ",0.5d0*(ic+ib)," (x2) and "
     &                ,ia," cm-1"
         iunique = 1
      else
         write(6,*)"Symmetrized to ",0.5d0*(ia+ic)," (x2) and "
     &                ,ib," cm-1"
         iunique = 2
      endif
      write(6,*)"Input rotational constants: "
      write(6,*)brot1*autocmi," (x2) ",brot2*autocmi," cm-1"
      write(6,*)
      endif 

c rotate to diagonalize mom
      do i=1,nclu
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

      if (brot2.lt.0.d0.and.nclu.ne.2) then
          write(6,*)"Choosing orientation of J evenly over a sphere."
          write(6,*)
 10       continue
          jx = sprng(rng_stream)*2.d0-1.d0
          jy = sprng(rng_stream)*2.d0-1.d0
          jz = sprng(rng_stream)*2.d0-1.d0
          tmp=dsqrt(jx**2+jy**2+jz**2)
          if (tmp.gt.1.d0) go to 10
          jx=jx/tmp*jj  ! normalize to a magnitude of JJ
          jy=jy/tmp*jj
          jz=jz/tmp*jj
          omeg(1) = jx/eig(1)
          omeg(2) = jy/eig(2)
          omeg(3) = jz/eig(3)
      else
          write(6,*)"Choosing orientation of J with K = ",kk
c          tmpx = jj/dsqrt(jj*(jj+1.d0))*kk      ! quantum
          tmpx = kk                             ! classical
          if (nclu.eq.2) kk=0.
          theta = sprng(rng_stream)*2.d0*pi
          tmpy=dsin(theta)*dsqrt(jj**2-tmpx**2)
          tmpz=dcos(theta)*dsqrt(jj**2-tmpx**2)
          if (iunique.eq.1) then
             jx=tmpx
             jy=tmpy
             jz=tmpz
             jx0=0.d0
             jy0=dsin(theta)*jj
             jz0=dcos(theta)*jj
          elseif (iunique.eq.2) then
             jx=tmpz
             jy=tmpx
             jz=tmpy
             jx0=dcos(theta)*jj
             jy0=0.d0
             jz0=dsin(theta)*jj
          else
             jx=tmpy
             jy=tmpz
             jz=tmpx
             jx0=dsin(theta)*jj
             jy0=dcos(theta)*jj
             jz0=0.d0
          endif

c         compute Erot for K=0
          omeg(1) = jx0/eig(1)
          omeg(2) = jy0/eig(2)
          omeg(3) = jz0/eig(3)
          tmp = 0.d0
          do i=1,nclu
           tmp = tmp + 0.5d0*mmm(i)*(omeg(2)*xx(3,i)-omeg(3)*xx(2,i))**2
           tmp = tmp + 0.5d0*mmm(i)*(omeg(3)*xx(1,i)-omeg(1)*xx(3,i))**2
           tmp = tmp + 0.5d0*mmm(i)*(omeg(1)*xx(2,i)-omeg(2)*xx(1,i))**2
          enddo
          erot0=tmp

c         compute Erot
          omeg(1) = jx/eig(1)
          omeg(2) = jy/eig(2)
          omeg(3) = jz/eig(3)
          tmp = 0.d0
          do i=1,nclu
           tmp = tmp + 0.5d0*mmm(i)*(omeg(2)*xx(3,i)-omeg(3)*xx(2,i))**2
           tmp = tmp + 0.5d0*mmm(i)*(omeg(3)*xx(1,i)-omeg(1)*xx(3,i))**2
           tmp = tmp + 0.5d0*mmm(i)*(omeg(1)*xx(2,i)-omeg(2)*xx(1,i))**2
          enddo
          erot=tmp
          ek2=erot-erot0

          print *,"Erot(K=0) = ",erot0*autoev," eV"
          print *,"Erot(K)   = ",erot*autoev," eV"
          print *,"EK2       = ",ek2*autoev," eV"
          print *

      endif

c Omeg = I-1 . J
c correction = r x omega
      do i=1,nclu
         pp(1,i)=pp(1,i)+mmm(i)*(omeg(2)*xx(3,i)-omeg(3)*xx(2,i))
         pp(2,i)=pp(2,i)+mmm(i)*(omeg(3)*xx(1,i)-omeg(1)*xx(3,i))
         pp(3,i)=pp(3,i)+mmm(i)*(omeg(1)*xx(2,i)-omeg(2)*xx(1,i))
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
      do i=1,nclu
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
