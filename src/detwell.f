      subroutine detwell(pe,xx,welli)

      implicit none
      include 'param.f'

      integer i,welli,i341,i342,i431,i432,i4,i3,it
      double precision pe,xx(12),r1,r2,r3,r4,r5,r6,
     & pets,pevin,chts,ccts,hhdis,hccts,
     & hpi,a341,a342,a431,a432

      hpi = 0.5d0*pi
      pets=16382d0/autocmi
      pevin=15410d0/autocmi
      chts=1.3592d0/autoang
      ccts=1.2619d0/autoang
      hhdis=5.d0/autoang
      hccts=66.46161d0/180.d0
c     r1: H1-C3, r2: H2-C4, r3: C3-C4
c     r4: H1-C4, r5: H2-C3, r6: H1-H2

      r1=sqrt((xx(1)-xx(7))**2+(xx(2)-xx(8))**2+(xx(3)-xx(9))**2)
      r2=sqrt((xx(4)-xx(10))**2+(xx(5)-xx(11))**2+(xx(6)-xx(12))**2)
      r3=sqrt((xx(7)-xx(10))**2+(xx(8)-xx(11))**2+(xx(9)-xx(12))**2)
      r4=sqrt((xx(1)-xx(10))**2+(xx(2)-xx(11))**2+(xx(3)-xx(12))**2)
      r5=sqrt((xx(4)-xx(7))**2+(xx(5)-xx(8))**2+(xx(6)-xx(9))**2)
      r6=sqrt((xx(1)-xx(4))**2+(xx(2)-xx(5))**2+(xx(3)-xx(6))**2)

      a341 = dacos((r1**2+r3**2-r4**2)/(2.d0*r1*r3))
      a342 = dacos((r5**2+r3**2-r2**2)/(2.d0*r5*r3))

      a431 = dacos((r4**2+r3**2-r1**2)/(2.d0*r4*r3))
      a432 = dacos((r2**2+r3**2-r5**2)/(2.d0*r2*r3))

      i341=0
      i342=0
      i431=0
      i432=0
      if (a341.ge.hpi) i341=i341+1
      if (a342.ge.hpi) i342=i342+1    
      if (a431.ge.hpi) i431=i431+1
      if (a432.ge.hpi) i432=i432+1
      i4=i341+i342
      i3=i431+i432
      it=i4+i3

      welli=0
      if (it.eq.2) then
        if (i4.eq.1.and.i3.eq.1) then
          welli=1
        elseif (i4.eq.0.and.i3.eq.2) then
          welli=2
        elseif (i4.eq.2.and.i3.eq.0) then
          welli=2
        else
          write(6,*)"impossible geometry"
        endif
      elseif (it.eq.1) then
        if (i3.eq.1) then
          if (i431.eq.1.and.a432.ge.hccts) then
            welli=2
          elseif (i432.eq.1.and.a431.gt.hccts) then
            welli=2
          elseif (a341.ge.hccts.or.a342.ge.hccts) then
            welli=1
          else
            write(6,*)"impossible geometry"
          endif
        elseif (i4.eq.1) then
          if (i341.eq.1.and.a342.gt.hccts) then
            welli=2
          elseif (i432.eq.1.and.a431.gt.hccts) then
            welli=2
          elseif (a341.ge.hccts.or.a342.ge.hccts) then
            welli=1
          else
            write(6,*)"impossible geometry"
          endif
        else
          write(6,*)"impossible geometry"
        endif
      elseif(it.eq.0) then
         if (a431.ge.hccts.and.a342.ge.hccts) then
           welli=1
         elseif (a432.ge.hccts.and.a341.ge.hccts) then
           welli=1
         elseif (a431.ge.hccts.and.a432.ge.hccts) then
           if (a341.lt.hccts.and.a342.lt.hccts) then
             welli=2
           else
             write(6,*)"impossible geometry"
           endif
         elseif (a341.ge.hccts.and.a342.ge.hccts) then
           if (a431.lt.hccts.and.a432.lt.hccts) then
             welli=2
           else
             write(6,*)"impossible geometry"
           endif
         else
           write(6,*)"impossible geometry"
         endif
       else
         write(6,*)"impossible geometry"
       endif

       return
       
       end


























