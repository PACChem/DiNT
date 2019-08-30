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
      double precision de,re,dx,dy,dz,rr,xe,dd,v,yy,eso,dvdr

      double precision h12x,h12y
      common/pesh12/h12x,h12y

      eso = 48.d0 ! S0/T1 MSX CO2
      re=4.d0/autoang

      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      yy=rr*autoang-re

      a1=-50.d0
      a2=5.d0
      v=yy*a1/2.d0+dsqrt((yy*a1)**2/4.d0+a2**2)
      v=v/autokcal
      pemd(2,2)=v
      dvdr=a1/2.d0+0.5d0/dsqrt((yy*a1)**2/4.d0+a2**2)*2.d0*yy*a1**2/4.d0
      dvdr=dvdr*autoang/autokcal
      gpemd(1,1,2,2) =  dvdr*dx/rr
      gpemd(1,2,2,2) = -dvdr*dx/rr
      gpemd(2,1,2,2) =  dvdr*dy/rr
      gpemd(2,2,2,2) = -dvdr*dy/rr
      gpemd(3,1,2,2) =  dvdr*dz/rr
      gpemd(3,2,2,2) = -dvdr*dz/rr

      a1=-13.d0
      a2=3.d0
      v=yy*a1/2.d0+dsqrt((yy*a1)**2/4.d0+a2**2)
      v=v/autokcal
      pemd(1,1)=v
      dvdr=a1/2.d0+0.5d0/dsqrt((yy*a1)**2/4.d0+a2**2)*2.d0*yy*a1**2/4.d0
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
