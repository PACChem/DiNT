      subroutine bootstrap(dd,wd,nd,md,id,perr)

      implicit none
c      include 'param.f'
c      include 'c_ran.f'
c#include <sprng_f.h>

      integer nd,md,id,iseed,nsamp,irand,j,i
      double precision dd(md),avgi(100000),x,avgavg,tot,wtot,
     & avg,std,perr,wd(md),rand
      integer*4 timeArray(3)

      call itime(timeArray)     ! Get the current time
      iseed = timeArray(1)+timeArray(2)+timeArray(3)
      x=rand(iseed)

      nsamp=1000

      avgavg=0.d0
      do j=1,nsamp
        tot=0.d0
        do i=1,id
          irand=int(rand()*dble(id))+1
          if (irand.le.nd) then
              tot=tot+dd(irand)*wd(irand)
          endif
        enddo
        avg = tot/dble(id)
        avgi(j)=avg
        avgavg=avgavg+avg
      enddo
      avgavg=avgavg/dble(nsamp)

      std=0.
      do j=1,nsamp
        std=std+(avgi(j)-avgavg)**2
      enddo
      std=std/dble(nsamp-1)
      std=dsqrt(std)

      perr=std/avg*100.d0
c      write(99,*)avg,std,std/avg*100.d0

      return
      end
