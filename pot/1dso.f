      subroutine pot(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)

      implicit none
      double precision autoang,autocmi,autoev,autokcal
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoev=27.2113961d0)
      parameter(autokcal=627.509d0)
      integer nat,mnat,nsurf,mnsurf
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     & gpemd(3,mnat,mnsurf,mnsurf),pemd(mnsurf,mnsurf)

      integer i,j
      double precision de,re,dx,dy,dz,rr,xe,dd,v,yy,eso,dvdr,yy0,re0,abc

      double precision bb,bb0

      double precision h12x,h12y
      common/pesh12/h12x,h12y

c masses should be 0.7 and 0.7
      eso = 48.d0 ! S0/T1 MSX CO2
      xe=-5.d0  !  MODEL 2
c      xe=56.d0  !  MODEL 1
      bb=1.d0
      bb0=1.d0
      de=50.d0  ! kcal/mol
      dd=45.d0  ! kcal/mol
      re=3.d0   ! A
      re0=3.d0   ! A

c masses should be 1.1 and 1.1
      eso = 78.77 ! S0/T1 collinear minimum
      xe=-140.d0  !  collinear 
      de=44.d0  ! kcal/mol
      dd=45.d0  ! kcal/mol
      re0=4.45d0   ! A
      re=4.d0  ! A
      bb0=1.d0
      bb=1.d0

      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      yy=rr*autoang*bb-re
      yy0=rr*autoang*bb0-re0

      v=de*dexp(-yy0)
      v=v/autokcal
      pemd(2,2)=v
      dvdr=-de*dexp(-yy0)*bb0
      dvdr=dvdr*autoang/autokcal
      gpemd(1,1,2,2) =  dvdr*dx/rr
      gpemd(1,2,2,2) = -dvdr*dx/rr
      gpemd(2,1,2,2) =  dvdr*dy/rr
      gpemd(2,2,2,2) = -dvdr*dy/rr
      gpemd(3,1,2,2) =  dvdr*dz/rr
      gpemd(3,2,2,2) = -dvdr*dz/rr

      v=(xe-dd)*dexp(-yy)+dd
      v=v/autokcal
      pemd(1,1)=v
      dvdr=-(xe-dd)*dexp(-yy)*bb
      dvdr=dvdr*autoang/autokcal
      gpemd(1,1,1,1) =  dvdr*dx/rr
      gpemd(1,2,1,1) = -dvdr*dx/rr
      gpemd(2,1,1,1) =  dvdr*dy/rr
      gpemd(2,2,1,1) = -dvdr*dy/rr
      gpemd(3,1,1,1) =  dvdr*dz/rr
      gpemd(3,2,1,1) = -dvdr*dz/rr

      pemd(1,2)=eso/autocmi
      pemd(2,1)=eso/autocmi
      h12x=pemd(1,2)
      do i=1,3
      do j=1,mnat
      gpemd(i,j,1,2)=0.d0
      gpemd(i,j,2,1)=0.d0
      enddo
      enddo

      return
      end

      subroutine prepot
      end
