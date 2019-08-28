      subroutine vwkb2(mm,edia,xj,xv,rin,rout,nsurf,symb)

c Computes the WKB vibrational action at the energy EDIA and
c for the rotational state XJ.  For the atom-diatom initial
c conditions method (INITx = 3) only.  From NAT8.1.

      implicit none
      include 'param.f'

      integer ncheb,i,nsurf
      double precision mm(mnat),edia,xv,summ,arg,x(12),wgt(12),
     &  v,rin,rout,rmass,term,r,rmid,rdif,xj
      character*2 symb(mnat)

c      print *,"aj in vwkb2"
      
      ncheb = 12
      do 5 i = 1,ncheb
        arg = dble(i)/dble(ncheb+1)*pi
        x(i) = dcos(arg)
        wgt(i) = pi/dble(ncheb+1)*(dsin(arg))**2
    5 continue

      call turn2(mm,edia,xj,rin,rout,nsurf,symb)
c      print *,'turning points',rin,rout
      if (rin .gt. rout) then                                  
        print*,'Outer turning point smaller than inner turning',
     &      ' point in vwkb, trajectory'
        stop
      else if ((rout - rin) .lt. 1.d-5) then
        xv = -0.5d0                    
      else                              
        rdif = 0.5d0*(rout-rin)
        rmid = 0.5d0*(rout+rin)
        summ = 0.d0
        do 10 i = 1,ncheb
          r = rmid+rdif*x(i)
          call diapot2(r,v,xj,mm,nsurf,symb)
c             print *,"aj",r,v
          if( (edia-v)/((r-rin)*(rout-r)) .lt. 0.0d0)then
            print*,'WARNING vwkb,  attempting ',
     &       'to take a sqrt of a negative number ',edia-v
            print *,'Change guess in diamin?'
            stop
          else
            term = sqrt((edia-v)/((r-rin)*(rout-r)))
          endif
          summ = summ + wgt(i)*term
   10   continue
      rmass = mm(1)*mm(2)/(mm(1)+mm(2))
      xv = sqrt(2.d0*rmass)*rdif**2*summ/pi-0.5d0
c      print *,xv
      endif

      return

      end
