      subroutine ejscale(im,jj,ejsc,etarget)

      implicit none
      include 'param.f'
      integer im
      double precision etarget,jj,ejsc(4,im)

c Vibrational dissociation energy 
      etarget=ejsc(1,im)+ejsc(2,im)*jj+ejsc(3,im)*jj**2
c     EJSC are input in eV and stay in eV

      etarget = etarget/autoev  ! eV --> au
      etarget = ejsc(4,im)*etarget ! scale by ejsc(4,im)

      return
      end
