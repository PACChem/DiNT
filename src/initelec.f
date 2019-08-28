      subroutine initelec

C PREPARE INITIAL ELECTRONIC VARIABLES
c Real and imaginary parts all set to 0 except for real part
c of the initially occupied state

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      integer i

      do i=1,nsurft
      crei(i) = 0.d0
      cimi(i) = 0.d0
      enddo
      crei(nsurf) = 1.d0

      do i=1,2*nsurft
      cre(i) = 0.d0
      cim(i) = 0.d0
      enddo
      cre(nsurf) = 1.d0
c     used for CSDM
      cre(nsurf+nsurft) = 1.d0
       
      return

      end
