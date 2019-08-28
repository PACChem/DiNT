      subroutine ransamp(lbins,ifems,ntot,filex2,filep2,nn,xx,pp,molpe)


      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

c in/out
      integer ntot,nn
      character*10 filex2,filep2
      character*16 filexx,filepp
      character*2 dum
      double precision pp(3,mnat),xx(3,mnat),molpe
      logical lbins,ifems

c local
      integer n,i,j,nrcl,itmp1,itmp2
      double precision r
      character*2 ctmp(nn)

      r=sprng(rng_stream)
      r=r*dble(ntot)
      n=nint(r)

      filexx="../r0/" // filex2
      filepp="../r0/" // filep2
c      filexx=filex2
c      filepp=filep2

      write(6,*)"Reading sample ",n

      if (lbins) then
 1080   format(a2,3f15.5)
 1081   format(i5,i7,f15.5)
        nrcl=(3+5*nn)*10
        open(93,file=filexx,form='unformatted',
     &   status='unknown',access='direct',RECL=nrcl)
        read(93,REC=n)itmp1,itmp2,molpe,(ctmp(i),i=1,nn),
     &   ((xx(j,i),j=1,3),i=1,nn)
        write(6,1081)itmp1,itmp2,molpe
        do i=1,nn
          write(6,1080)ctmp(i),(xx(j,i),j=1,3)
          xx(1,i) = xx(1,i)/autoang
          xx(2,i) = xx(2,i)/autoang
          xx(3,i) = xx(3,i)/autoang
        enddo
      else
        open(unit=93,file=filexx)
        read(93,*)
        read(93,*)
        do i=1,n-1
          read(93,*)
          read(93,*)
          do j=1,nn
            read(93,*)
          enddo
        enddo
        read(93,*)
        read(93,*)
        do j=1,nn
          read(93,*)dum,xx(1,j),xx(2,j),xx(3,j)
          xx(1,j)=xx(1,j)/autoang
          xx(2,j)=xx(2,j)/autoang
          xx(3,j)=xx(3,j)/autoang
        enddo
      endif
      close(93)
      
      if (ifems) return
      if (lbins) then
        open(94,file=filepp,form='unformatted',
     &   status='unknown',access='direct',RECL=nrcl)
        read(94,REC=n)itmp1,itmp2,molpe,(ctmp(i),i=1,nn),
     &   ((pp(j,i),j=1,3),i=1,nn)
        write(6,1081)itmp1,itmp2,molpe
        do i=1,nn
          write(6,1080)ctmp(i),(pp(j,i),j=1,3)
        enddo
      else
        open(unit=94,file=filepp)
        read(94,*)
        read(94,*)
        do i=1,n-1
          read(94,*)
          read(94,*)
          do j=1,nn
            read(94,*)
          enddo
        enddo
        read(94,*)
        read(94,*)dum,molpe
        do j=1,nn
          read(94,*)dum,pp(1,j),pp(2,j),pp(3,j)
        enddo
      endif
      close(94)

      return

      end
