      subroutine getplz(pp,mm,gpem,nat,plz,elz,h12,gtmp,pairy,n0)

      implicit none
      include 'param.f'

      integer i,j,nat,n0
      double precision gtmp,h12,mm(mnat),pp(3,mnat),gpem(3,mnat,mnsurf),
     & pdot,elz,plz,h12x,h12y,gtmp1,gtmp2,tmpf,e0,be,ex,xxx,pairy(2),
     & a1,a2,a3,a4,f1,f2

      common/pesh12/h12x,h12y

      gtmp=0.d0
      gtmp1=0.d0
      gtmp2=0.d0
      do i=1,nat
      do j=1,3
        gtmp=gtmp+((gpem(j,i,1)-gpem(j,i,2))**2)*mu/mm(i)
        gtmp1=gtmp1+(gpem(j,i,1)**2)*mu/mm(i)
        gtmp2=gtmp2+(gpem(j,i,2)**2)*mu/mm(i)
      enddo
      enddo
      gtmp=dsqrt(gtmp)
      gtmp1=dsqrt(gtmp1)
      gtmp2=dsqrt(gtmp2)

      pdot=0.d0
      f1=0.d0
      f2=0.d0
      do i=1,nat
      do j=1,3
        pdot=pdot+((gpem(j,i,1)-gpem(j,i,2)))*pp(j,i)/gtmp*mu/mm(i)
        f1=f1-((gpem(j,i,1)-gpem(j,i,2)))*gpem(j,i,1)/gtmp*mu/mm(i)
        f2=f2-((gpem(j,i,1)-gpem(j,i,2)))*gpem(j,i,2)/gtmp*mu/mm(i)
      enddo
      enddo

      elz=pdot**2/(2.d0*mu)

      h12=h12x

             plz=1.d0-dexp(-2.d0*pi*h12**2/gtmp
     &              *dsqrt(0.5d0*mu/elz))  ! LJ prob

c        write(6,'(a,10f13.5)')'getplz',h12,mu,elz,gtmp,f1,f2

! airy
       tmpf=dsqrt(dabs(f1*f2))
       e0=(tmpf**4/(2.d0*mu*gtmp**2))**(1.d0/3.d0)
       be=(2.d0*h12*tmpf/e0/gtmp)**(3.d0/2.d0)
       ex=elz*gtmp/(2.d0*h12*tmpf)
       xxx=-ex*(be**(2.d0/3.d0))
       call AIRYA(xxx,a1,a2,a3,a4)
       pairy(1)=pi**2*be**(4.d0/3.d0)*a1**2

       tmpf=dsqrt(dabs(f1*f2))
       if (n0.eq.1) tmpf=dabs(f1)
       if (n0.eq.2) tmpf=dabs(f2)
       e0=(tmpf**4/(2.d0*mu*gtmp**2))**(1.d0/3.d0)
       be=(2.d0*h12*tmpf/e0/gtmp)**(3.d0/2.d0)
       ex=elz*gtmp/(2.d0*h12*tmpf)
       xxx=-ex*(be**(2.d0/3.d0))
       call AIRYA(xxx,a1,a2,a3,a4)
       pairy(2)=pi**2*be**(4.d0/3.d0)*a1**2

!       be=4.d0*h12*dsqrt(mu*h12/tmpf/dabs(f1-f2))
!       xxx=-ex*(be**(2.d0/3.d0))
!       call AIRYA(xxx,a1,a2,a3,a4)
!       pairy=pi**2*be**(4.d0/3.d0)*a1**2

!       print *,elz*627.509,e0,be,ex,xxx,a1,pairy,tmpf,gtmp1,gtmp2,gtmp
!      print *,elz*autoev,gtmp,plz,pairy
!       print *,gtmp,gtmp1,gtmp2

!      stop

      return

      end 

