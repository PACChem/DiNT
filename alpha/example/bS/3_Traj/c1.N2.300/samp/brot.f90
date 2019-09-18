program brot

implicit double precision (a-h,o-z)
integer seed
parameter(mnat=100)
character*3 s(mnat)
double precision x(3,mnat),m(mnat),dum(3),mom(3,3),ap(6),eig(3),rot(3,3),work(9)
parameter(autoang=0.52917706d0)
parameter(autocmi=219474.63067d0)
parameter(amutoau=1822.844987d0)

!call itime(time)
!seed=time(3)*24*60+time(2)*60+time(3)
call system_clock(seed)
call srand(seed) 
r1=rand()

open(10,file="fort.80")
read(10,*)
read(10,*)
ngeo=0
do while (.true.)
read(10,*,END=100)natom
read(10,*)
do i=1,natom
read(10,*)
enddo
ngeo=ngeo+1
if (ngeo.eq.10000) go to 100
enddo
100 continue

nsamp=1000
b1a=0.d0
b2a=0.d0
do is=1,nsamp
rewind(10)

j=int(rand()*real(ngeo))

read(10,*)
read(10,*)

k=0
do while (k.lt.j)
read(10,*)natom
read(10,*)
do i=1,natom
read(10,*)
enddo
k=k+1
enddo

read(10,*)natom
read(10,*)
do i=1,natom
read(10,*)s(i),x(1,i),x(2,i),x(3,i)
do l=1,3
x(l,i)=x(l,i)/autoang
enddo
if (s(i).eq."C") m(i)=12.d0*amutoau
if (s(i).eq."H") m(i)=1.007825d0*amutoau
if (s(i).eq."N") m(i)=14.003074d0*amutoau
enddo

      do i=1,3
      do j=1,3
          mom(i,j) = 0.d0
      enddo
      enddo

      do i=1,natom
         mom(1,1)=mom(1,1)+m(i)*(x(2,i)**2+x(3,i)**2)
         mom(2,2)=mom(2,2)+m(i)*(x(1,i)**2+x(3,i)**2)
         mom(3,3)=mom(3,3)+m(i)*(x(1,i)**2+x(2,i)**2)
         mom(1,2)=mom(1,2)-m(i)*(x(1,i)*x(2,i))
         mom(1,3)=mom(1,3)-m(i)*(x(1,i)*x(3,i))
         mom(2,3)=mom(2,3)-m(i)*(x(2,i)*x(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

      do i=1,3
      do j=i,3
        ap(i+(j-1)*j/2)=mom(i,j)
      enddo
      enddo
      call dspev( 'v','u',3,ap,eig,rot,3,work,info )

      do i=1,3
      if (eig(i).le.0.d0) eig(i)=1.d-16
      enddo
      a=0.5d0/eig(1)*autocmi
      b=0.5d0/eig(2)*autocmi
      c=0.5d0/eig(3)*autocmi

      if (dabs(a-b).lt.dabs(b-c)) then 
      b1=(a+b)/2.d0
      b2=c
      else
      b1=(c+b)/2.d0
      b2=a
      endif

!      write(6,*)is,a,b,c,b1,b2
      b1a=b1a+b1
      b2a=b2a+b2
enddo

open(11,file="dint.brot")
write(11,*)b1a/dble(nsamp),b2a/dble(nsamp)

end

