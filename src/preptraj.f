      subroutine preptraj

c Put center of mass at origin

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
 
      double precision xtot,ytot,ztot,mtot

      integer i
      call gettemp(pp,mm,nat,temp,ke)
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)

      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
          xtot = xtot + mm(i)*xx(1,i)
          ytot = ytot + mm(i)*xx(2,i)
          ztot = ztot + mm(i)*xx(3,i)
          mtot = mtot + mm(i)
      enddo
      do i=1,nat
          xx(1,i) = xx(1,i) - xtot/mtot
          xx(2,i) = xx(2,i) - ytot/mtot
          xx(3,i) = xx(3,i) - ztot/mtot
      enddo

      call gettemp(pp,mm,nat,temp,ke)
      write(6,*)"Placed CoM at origin"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)

c     remove overall momenta
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
        xtot=xtot+pp(1,i)
        ytot=ytot+pp(2,i)
        ztot=ztot+pp(3,i)
        mtot = mtot + mm(i)
      enddo
      do i=1,nat
        pp(1,i)=pp(1,i)-xtot*mm(i)/mtot
        pp(2,i)=pp(2,i)-ytot*mm(i)/mtot
        pp(3,i)=pp(3,i)-ztot*mm(i)/mtot
      enddo
      call gettemp(pp,mm,nat,temp,ke)
      write(6,*)"Removed CoM motion"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)

      return

      end
