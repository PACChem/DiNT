      subroutine nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,d,
     & itype,mcpar,v)

c     returns the energy EE of a structure displaced from XX0 along NN
c     mass-scaled vector VEC by some distance DD

      implicit none
      include 'param.f'
      integer i,j,k,nclu,repflag,nn,itype(3*mnat),mcpar(2,3*mnat),
     & a1,a2,n1,nsurft,maxstep,im
      double precision xx(3,mnat),vec(3,mnat,3*mnat),
     & xx0(3,mnat),d(3*mnat),mm(mnat),pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),v(mnsurf),phi,gtmp,gtmp1,gapgrad,h12
      character*2 symb(mnat)

      common/lz/gapgrad,h12

      do j=1,3
      do k=1,nclu
        xx(j,k)=xx0(j,k)
      enddo
      enddo

c Do Cartesian normal mode displacements first
      do i=1,nn
        if (itype(i).eq.0) then ! Cartesian
        do j=1,3
        do k=1,nclu
          xx(j,k)=xx(j,k)+vec(j,k,i)*d(i)*dsqrt(mu/mm(k))
        enddo
        enddo
        endif
      enddo

c           write(80,*)nclu
c          write(80,1081)i,phi
c          do k=1,nclu
c            write(80,1080)symb(k),(xx(j,k)*autoang,j=1,3)
c          enddo

C Do torsions next
      do i=1,nn
        if (itype(i).eq.1) then ! torsion
        a1=mcpar(1,i)
        a2=mcpar(2,i)
        n1=a2-1
        phi=d(i)
           call torphi(xx,nclu,a1,a2,phi,n1)
c           write(80,*)nclu
c           write(80,1081)i,phi
c           do k=1,nclu
c             write(80,1080)symb(k),(xx(j,k)*autoang,j=1,3)
c           enddo
c 1080       format(a2,3f15.5)
c 1081       format(i7,50f15.5)
c          stop
        endif
      enddo

c get energy
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
      do i=1,nsurft
      if (repflag.eq.1) then
        v(i) = pemd(i,i)
      else
        v(i) = pema(i)
      endif
      enddo

c MINIMIZE GAP TO MSX  TEMP AJ AJ AJ
      go to 10 ! SKIP
      gtmp=0.d0
      do i=1,nclu
      do j=1,3
        gtmp=gtmp+((gpemd(j,i,1,1)-gpemd(j,i,2,2))**2)*mu/mm(i)
      enddo
      enddo
      gtmp=dsqrt(gtmp) 
      maxstep=50
      do im=1,maxstep
      if (dabs(v(2)-v(1))*autoev.lt.0.01d0) go to 10
      print *,'min',im,v(1)*autoev,v(2)*autoev,0.5d0*(v(1)+v(2))*autoev,
     &   (v(1)-v(2))*autoev 
      gtmp=0.d0
      gtmp1=0.d0
      do i=1,nclu
      do j=1,3
        gtmp=gtmp+((gpemd(j,i,1,1)-gpemd(j,i,2,2))**2)*mu/mm(i)
        gtmp1=gtmp1+((gpemd(j,i,1,1)-gpemd(j,i,2,2))**2)
      enddo
      enddo
      gtmp=dsqrt(gtmp) 
      gtmp1=dsqrt(gtmp1) 
      do i=1,nclu
      do j=1,3
        xx(j,i)=xx(j,i)-2.d0*(v(2)-v(1))
     *    *(gpemd(j,i,2,2)-gpemd(j,i,1,1))/gtmp1
      enddo
      enddo
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
      do i=1,nsurft
        v(i) = pemd(i,i)
      enddo
      enddo
c MINIMIZE GAP TO MSX  TEMP AJ AJ AJ

 10   continue

          write(80,*)nclu
          write(80,1081)0,v(1)*autoev,v(2)*autoev,gtmp,
     & pemd(1,2)*autocmi
          do i=1,nclu
            write(80,1080)symb(i),(xx(j,i)*autoang,j=1,3)
 1080       format(a2,3f15.5)
 1081       format(i7,50f15.5)
          enddo

      gapgrad=gtmp
      h12=pemd(1,2)

      return
      end

