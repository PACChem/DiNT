      subroutine radialdist(itmp,xx,nclu,step,time,
     &   nbinrad,raddist,iprint)

c really the pair correlation function

      implicit none
      include 'param.f'

      integer i,j,k,nclu,ibin,nbinrad,iprint,itmp
      double precision xx(3,mnat),rij,raddist(0:nbinrad+1),tmp,
     & step,rbin,time,xnorm,rmax,rmin

      rmin = 0.d0 / autoang
      rmax = 30.d0 / autoang
      rbin = (rmax-rmin)/dble(nbinrad)
    
      do i=1,nclu
      do j=i+1,nclu
        rij = 0.d0
        do k=1,3
        rij = rij + (xx(k,i)-xx(k,j))**2
        enddo
        rij = dsqrt(rij)
        if (rij.le.rmin) then
          ibin = 0
        elseif (rij.ge.rmax) then
          ibin = nbinrad + 1
        else
          rij = rij - rmin
          ibin = int(rij/rbin) + 1 ! determine bin number
          ibin = max(1,ibin)
          ibin = min(nbinrad,ibin)
        endif
        raddist(ibin) = raddist(ibin) + step
      enddo
      enddo

      if (iprint.eq.1) then
      tmp = time*dble(nclu)
      tmp = tmp*rbin*4.d0*pi
c      http://www.physics.emory.edu/~weeks/idl/gofr2.html
      write(41,10)itmp,time*autofs,raddist(0)/time,
     &  (raddist(k)/(tmp*((dble(k)-0.5d0)*rbin+rmin)**2),k=1,nbinrad),
     &   raddist(nbinrad+1)/time
      endif
 10   format(i7,100f12.4)

      return
 
      end
