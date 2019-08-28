      subroutine rtox(xx,r,ido)

c Compute three atom-atom distances from Cartesian coordinates and
c vice-versa.  For three systems only.

      implicit none
      include 'param.f'

      integer i,j,ido
      double precision r(3),xx(3,mnat),costh,sinth

      if (ido.eq.1) then
c     transform xx to internals
        r(1)=0.d0
        r(2)=0.d0
        r(3)=0.d0
      do i=1,3
        r(1) = r(1) + (xx(i,1)-xx(i,2))**2
        r(2) = r(2) + (xx(i,2)-xx(i,3))**2
        r(3) = r(3) + (xx(i,1)-xx(i,3))**2
      enddo
      do i=1,3
        r(i)=dsqrt(r(i))
      enddo

      else
c     transform internals to xx 
      costh = ( -r(3)**2+r(1)**2+r(2)**2 )/(2.d0*r(1)*r(2))
      costh = max(-1.d0,costh)
      costh = min(1.d0,costh)
      sinth = dsqrt(1.d0-costh**2)
c coordinate system:  B at origin, A on x-axis x > 0, C in xy plane y > 0
c r(1) = rAB
c r(2) = rBC
c r(3) = rAC
      xx(1,1) = r(1)
      xx(2,1) = 0.d0
      xx(3,1) = 0.d0
      xx(1,2) = 0.d0
      xx(2,2) = 0.d0
      xx(3,2) = 0.d0
      xx(1,3) = r(2)*costh
      xx(2,3) = r(2)*sinth
      xx(3,3) = 0.d0
      endif

      return

      end
