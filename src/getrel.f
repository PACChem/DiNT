      subroutine getrel(rcom,eint)

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'

      integer i,j,k,arr,index,nprod
      double precision rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),tmp,tmpi,
     & rr(3),jp(9),jq(9),jm(9),xv,xj,rin,rout,evib,erotx,
     & ediatot,kerel,kedia,xjx,xjy,xjz,eint,mmtot,rcom

c     temp
      integer wellindex1,wellindex2
      double precision xwell,ptmp

c     for unit 31
      integer ifrag,nfrag(2),fragind(mnat,2),itmp,itmp1,itmp2
      double precision rr1,rr2,ppfrag(3,mnat),xxfrag(3,mnat),
     &  mmfrag(mnat),bigjag(3,2),bigjtotag(2),comxag(3,2),compag(3,2),
     &  relmu,kerelfrag(2),jfragcom(2),erotag(3,2),erottotag(2),
     &  jfragrot(3,2),jfragrottot,rf,comx(3),comp(3),lorb(3),lorbtot,
     &  eorb,erel,tmp3j(3),tmp3e(3),mmfragtot(2),tempfrag,kefrag(2),
     &  evibag(2),tdelay

      double precision tmpprint(50)
      common/tmp/tmpprint

C FOR DISSOCIATION REACTIONS
      nfrag(1) = natom(1)
      nfrag(2) = natom(2)
      mmtot = mmag(1)+mmag(2)
      mmfragtot(1) = mmag(1)
      mmfragtot(2) = mmag(2)
      do i=1,3
      comx(i)=0.d0
      comp(i)=0.d0
      comxag(i,1)=0.d0
      comxag(i,2)=0.d0
      compag(i,1)=0.d0
      compag(i,2)=0.d0
      enddo

      ifrag=1
      do i=1,nfrag(1)
        do k=1,3
          comx(k)=comx(k) + xx(k,i)*mm(i)                  ! center of mass
          comp(k)=comp(k) + pp(k,i)                        ! overall momentum
          comxag(k,ifrag)=comxag(k,ifrag) + xx(k,i)*mm(i)  ! center of mass for fragment
          compag(k,ifrag)=compag(k,ifrag) + pp(k,i)        ! overall momentum for fragment
        enddo
      enddo
      do k=1,3
        comx(k)=comx(k)/mmtot                                 ! overall center of mass
      enddo

      ifrag=2
      do i=nfrag(1)+1,nfrag(1)+nfrag(2)
        do k=1,3
          comx(k)=comx(k) + xx(k,i)*mm(i)                  ! center of mass
          comp(k)=comp(k) + pp(k,i)                        ! overall momentum
          comxag(k,ifrag)=comxag(k,ifrag) + xx(k,i)*mm(i)  ! center of mass for fragment
          compag(k,ifrag)=compag(k,ifrag) + pp(k,i)        ! overall momentum for fragment
        enddo
      enddo
      do k=1,3
        comx(k)=comx(k)/mmtot                                 ! overall center of mass
      enddo


      do j=1,2
      do k=1,3
        comxag(k,j)=comxag(k,j)/mmfragtot(j)-comx(k)          ! center of mass for fragment
        compag(k,j)=compag(k,j)-comp(k)/mmtot*mmfragtot(j)    ! overall momentum for fragment
      enddo
      enddo

c     relative reduced mass
      relmu = mmfragtot(1)*mmfragtot(2)/(mmfragtot(1)+mmfragtot(2))
      rf = 0.
      do i=1,3
        rf = rf + (comxag(i,1)-comxag(i,2))**2
      enddo
      rf = dsqrt(rf)  ! final fragment separation
      lorb(1) = (comxag(2,1)-comxag(2,2))*compag(3,1)
     &         -(comxag(3,1)-comxag(3,2))*compag(2,1)
      lorb(2) = (comxag(3,1)-comxag(3,2))*compag(1,1)
     &         -(comxag(1,1)-comxag(1,2))*compag(3,1)
      lorb(3) = (comxag(1,1)-comxag(1,2))*compag(2,1)
     &         -(comxag(2,1)-comxag(2,2))*compag(1,1)
      lorbtot = dsqrt(lorb(1)**2+lorb(2)**2+lorb(3)**2)
      eorb = lorbtot**2/(2*relmu*rf**2)
      eint = (compag(1,1)**2+compag(2,1)**2+compag(3,1)**2)/(2.d0*relmu)
      erel = eint - eorb

c      print *,"xxx",lorbtot,rf*autoang,eorb*autoev,eint*autoev,
c     &                erel*autoev

      tmpprint(2)=rf
      rcom=rf

      return
      end
