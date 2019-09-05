      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom,ccc)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension x3(maxatom),y3(maxatom),z3(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      dimension dvdx3(maxatom),dvdy3(maxatom),dvdz3(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      character*2 symb(maxatom),symb2(maxatom),symb3(maxatom)
      integer at(maxatom),at2(maxatom),at3(maxatom)

      COMMON/USROCM/ PENGYGS,PENGYES(5),
     +               PENGYIJ(15),
     +               DGSCART(25,3),DESCART(25,3,5),
     +               DIJCART(25,3,15)
C
      COMMON/USRICM/ CART(25,3),ANUZERO,
     +               NULBL(25),NFLAG(20),
     +               NASURF(5+1,5+1),NDER
C
C READ INPUT POSITIONS,DERIVATIVES, EXCITED STATE, NFLAG1-5 INFO
      n=6
      ncart=18
      natoms=8
      nder=1
      nasurf=0
      nflag(1)=0
      nflag(2)=0
      nflag(3)=0
      nflag(4)=0
      nflag(5)=0
      do i = 1,n
        cart(i,1)=x(i)
        cart(i,2)=y(i)
        cart(i,3)=z(i)
      enddo

      call pot
 
      v=pengygs           
      do i=1,n
        dvdx(i)=dgscart(i,1)
        dvdy(i)=dgscart(i,2)
        dvdz(i)=dgscart(i,3)
      enddo

      return
      end




