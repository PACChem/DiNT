      subroutine getrho(cre,cim,rhor,rhoi,nsurft)

c Compute the density matrix from the electronic coefficients.
c Note:  Phases are handled separately, and RHOR and RHOI do
c not include phases.  These may need to be added to RHOR and RHOI,
c depending on how these quantities are used.

      implicit none
      include 'param.f'

      integer nsurft,i,j
      double precision cim(2*mnsurf),cre(2*mnsurf),
     & rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf)

      do i=1,nsurft
      do j=1,nsurft
        rhor(i,j) = cre(i)*cre(j)+cim(i)*cim(j)
        rhoi(i,j) = cre(i)*cim(j)-cre(j)*cim(i)
      enddo
      enddo

      return

      end
