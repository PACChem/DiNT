      subroutine ranj(jmin,jmax,jtemp,jj,kk,brot1,brot2)

c picks J & K

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu
      double precision jj,jmax,jmin,jtemp,brot1,brot2,kk

c local
      integer i,ix,iy,iz,j,info,istep,istepmax,ij,ik,imin,imax
      parameter(istepmax=20000)
      double precision r1,r2,x1,x2,w,y1,y2,tmp,mom(3,3),ap(6),
     & eig(3),rot(3,3),work(9),temp1,temp2,temp3,omeg(3),
     & d1,d2,d3,d4,d5,d6,det,c11,c12,c13,c21,c22,c23,c31,c32,c33,
     & jstep,erot,pj,pjk,pjsum,pjksum,ji,ki,itmp,ptmp(istepmax)
      logical ljstyle

c HARD-CODED
c      ljstyle = .true. ! quantal J dist
      ljstyle = .false. ! classical J dist

c pick total J in one of two ways
      if (jtemp.lt.0.d0) then
c (1) pick evenly from jmin to jmax
      r1 = sprng(rng_stream)
      jj=jmin+r1*(jmax-jmin)
c     pick evenly from -j to j
      r1 = sprng(rng_stream)
      kk=-jj+r1*2.d0*jj
      write(6,*)"J = ",jj," hbar, selected evenly from ",
     &  jmin," to ",jmax," hbar"
      write(6,*)"K = ",kk," hbar, selected evenly from ",
     &  -jj," to ",jj," hbar"
      write(6,*)
      else
c (2) select from a thermally weighted dist
      jstep = 0.1d0  ! J step size
      istep = nint((jmax-jmin)/jstep)
      imin = nint(jmin/jstep)
      imax = nint(jmax/jstep)
      if (imax.gt.istepmax) then
       print *,istepmax,"<",imax
       stop
      endif
      pjsum=0.d0
      pjksum=0.d0
      do i=imin,imax   ! loop over J
        ji=dble(i)*jstep
        if (ljstyle) erot=(ji*(ji+1.d0))*brot1          ! rotational energy
        if (.not.ljstyle) erot=ji*ji*brot1          ! rotational energy
        pj=dexp(-erot/(kb*jtemp))*(2.d0*ji+1.d0)  ! Boltz weight including degeneracy
        pjsum=pjsum+pj
        ptmp(i) = 0.d0
      do j=-i,i       ! loop over K
        ki=dble(j)*jstep
        if (ljstyle) erot=(ji*(ji+1.d0))*brot1          ! rotational energy
        if (.not.ljstyle) erot=ji*ji*brot1          ! rotational energy
        erot=erot+ki*ki*(brot2-brot1)      ! K dependent term
        pjk=dexp(-erot/(kb*jtemp))*(2.d0*ji+1.d0)  ! Boltz weight including degeneracy
        pjksum=pjksum+pjk
        ptmp(i)=ptmp(i) + pjk
      enddo
      enddo

c      tmp =0.d0
c      do i=imin,imax   ! loop over J
c        ji=dble(i)*jstep
c        if (mod(i,50).eq.0) print *,i,ji,ptmp(i)/pjksum/0.1d0
c      enddo 
c      stop

c tmp
c      do i=1,istep
c        ji=dble(i-1)*jstep+jmin
c        print *,ji,pj(i)
c      enddo
c      do i=1,istep
c        ji=dble(i-1)*jstep+jmin
c        itmp=nint(ji/jstep)
c      do j=1,2*itmp+1
c        ki=dble(j-1)*jstep-ji
c        print *,ji,ki,pjk(i,j)
c      enddo
c      enddo
c      stop

      IF (brot2.lt.0.d0) then

c spherical top
      r1 = sprng(rng_stream)  ! random number
      tmp = 0.d0
      i=-1
      do while (tmp.lt.r1)
        i=i+1
        ji=dble(i)*jstep+jmin
        if (ljstyle) erot=(ji*(ji+1.d0))*brot1          ! rotational energy
        if (.not.ljstyle) erot=ji*ji*brot1          ! rotational energy
        pj=dexp(-erot/(kb*jtemp))*(2.d0*ji+1.d0)  ! Boltz weight including degeneracy
        tmp=tmp+pj/pjsum
c        print *,r1,tmp
      enddo 
      jj=ji
      kk=0
      write(6,*)"J = ",jj," hbar, selected from ",jmin," to ",
     &   jmax," hbar and ",jtemp," K thermal dist"
      write(6,*)

      ELSE

c symmetric top
      r1 = sprng(rng_stream)  ! random number
      tmp=0.d0
c      i=-1
c      j=-1
      i=imin-1
      j=imin-1
      do while (tmp.lt.r1.and.j.le.imax)
        if (j.eq.i) then
            i=i+1
            j=-i
        else
            j=j+1
        endif
        ji=dble(i)*jstep
        ki=dble(j)*jstep
        if (ljstyle) erot=(ji*(ji+1.d0))*brot1          ! rotational energy
        if (.not.ljstyle) erot=ji*ji*brot1          ! rotational energy
        erot=erot+ki*ki*(brot2-brot1)      ! K dependent term
        pjk=dexp(-erot/(kb*jtemp))*(2.d0*ji+1.d0)  ! Boltz weight including degeneracy
        tmp=tmp+pjk/pjksum
c        print *,i,j,ji,ki,tmp,r1
      enddo

      jj=ji
      kk=ki
      write(6,*)"J = ",jj," hbar and K = ",kk," hbar, selected from "
     & ,jmin," to ",jmax," hbar and ",jtemp," K thermal dist"
      write(6,*)

      ENDIF

      endif

      return

      end
