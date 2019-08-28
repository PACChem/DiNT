      subroutine period2(edia,jj,mm,rin,rout,peri,nsurf,symb)

c used by the atom-diatom initial conditions method
c computes the period in au of the diatomic vibration.

      implicit none

c input
      integer nsurf
      double precision rin,rout,mm(2),edia,xj,jj
      character*2 symb(2)

c output
      double precision peri,ptaudi(50)

c local
      integer ncheb,i,nstep
      double precision arg,x(12),rmid,rdif,r,term,v,sum,rmass,pi,step

      xj = jj
      pi=dacos(-1.d0)

      ncheb = 12
      do 5 i = 1,ncheb
        arg = dble(2*i-1)/dble(2*ncheb)*pi
        x(i) = dcos(arg)
    5 continue

      rdif = (rout-rin)/2.d0
      rmid = (rout+rin)/2.d0
      sum = 0.d0
      do 10 i = 1,ncheb
        r = rmid+rdif*x(i)
        call diapot2(r,v,xj,mm,nsurf,symb)
        term = sqrt((r-rin)*(rout-r)/(edia-v))
        sum = sum+term
   10 continue
      rmass = mm(1)*mm(2)/(mm(1)+mm(2))
c accurately integrated vibrational period
      peri = sqrt(2.d0*rmass)*sum*pi/dble(ncheb)

c approx integrate the vibrational period to store evenly spaced probabilities
      sum = 0.d0
      nstep=50
      step=(rout-rin)/dble(nstep)
      do i=0,nstep-1
        r=rin+step*(dble(i)+0.5d0)
        call diapot2(r,v,xj,mm,nsurf,symb)
        term = 1.d0/dsqrt(edia-v)*step
        ptaudi(i+1)=term
        sum = sum+term
      enddo
      rmass = mm(1)*mm(2)/(mm(1)+mm(2))
      print *,'period check in period2',sum*sqrt(2.d0*rmass),peri

      do i=1,nstep
        ptaudi(i)=ptaudi(i)/sum
      enddo

      return

      end

