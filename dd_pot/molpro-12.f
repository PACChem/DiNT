      subroutine pot(symbol,x,y,z,pemd,gpemd,nclu,mnclu,nsurf,mnsurf)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,idum,ithis(mnclu),nsurf,mnsurf
      character*2 symbol(mnclu),cdum
      character*80 string,tmpstr
      character*21 Estr
      character*18 Gstr
      character*18 Astr
      double precision x(mnclu),y(mnclu),z(mnclu),v
      double precision gpemd(3,mnclu,mnsurf,mnsurf),pemd(mnsurf,mnsurf)
      double precision xtmp,ytmp,ztmp,total,dum
      double precision dx(mnclu),dy(mnclu),dz(mnclu)
      double precision cfloat
      double precision autoang,so
      parameter(autoang=0.52917706d0)

      call pot1(symbol,x,y,z,v,dx,dy,dz,nclu,mnclu)
      pemd(1,1)=v
      do i=1,nclu
      gpemd(1,i,1,1)=dx(i)
      gpemd(2,i,1,1)=dy(i)
      gpemd(3,i,1,1)=dz(i)
      enddo

      call pot2(symbol,x,y,z,v,dx,dy,dz,nclu,mnclu)
      pemd(2,2)=v
      do i=1,nclu
      gpemd(1,i,2,2)=dx(i)
      gpemd(2,i,2,2)=dy(i)
      gpemd(3,i,2,2)=dz(i)
      enddo

c SO constant for now
      pemd(1,2)=0.0
      pemd(2,1)=0.0
      do i=1,nclu
      gpemd(1,i,2,1)=0.
      gpemd(2,i,2,1)=0.
      gpemd(3,i,2,1)=0.
      gpemd(1,i,1,2)=0.
      gpemd(2,i,1,2)=0.
      gpemd(3,i,1,2)=0.
      enddo

      return
      end
