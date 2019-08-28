      subroutine nsphere(n,r,v)

c     Returns the volume v of a hypersphere of dimension n
c     with radius r, calculated using recurrence relations
c     Jasper 4/6/11

      implicit none
      integer n,i
      double precision r,v,pi

      pi=dacos(-1.d0)

      if (n.eq.1) then
        v=2.d0*r
      elseif (n.eq.2) then
        v=pi*r**2
      elseif (mod(n,2).eq.0) then
        v=pi*r**2
        do i=4,n
        if (mod(i,2).eq.0) v=v*2.d0*pi*r**2/dble(i)
        enddo
      else
        v=2.d0*r
        do i=3,n
        if (mod(i,2).eq.1) v=v*2.d0*pi*r**2/dble(i)
        enddo
      endif

      return
      end

