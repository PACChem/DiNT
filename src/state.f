      subroutine state

c Computes some quantities on the fly.
c Writes some output files.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'

      integer i,j,k,arr,index,nprod
      double precision rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),tmp,tmpi,
     & rr(3),jp(9),jq(9),jm(9),xv,xj,rin,rout,evib,erotx,theta,
     & ediatot,kerel,kedia,xjx,xjy,xjz,eint,mmtot

c     temp
      integer wellindex1,wellindex2
      double precision xwell,ptmp,etmp

c     for unit 31
      integer ifrag,nfrag(2),fragind(mnat,2),itmp,itmp1,itmp2
      double precision rr1,rr2,ppfrag(3,mnat),xxfrag(3,mnat),eig(3),
     &  mmfrag(mnat),bigjag(3,2),bigjtotag(2),comxag(3,2),compag(3,2),
     &  relmu,kerelfrag(2),jfragcom(2),erotag(3,2),erottotag(2),
     &  jfragrot(3,2),jfragrottot,rf,comx(3),comp(3),lorb(3),lorbtot,
     &  eorb,erel,tmp3j(3),tmp3e(3),mmfragtot(2),tempfrag,kefrag(2),
     &  evibag(2),tdelay

c     for fragment pes calls
      integer nf
      double precision pemaf(mnsurf),pemdf(mnsurf,mnsurf),
     & gpemaf(3,mnat,mnsurf),gpemdf(3,mnat,mnsurf,mnsurf),
     & df(3,mnat,mnsurf,mnsurf),peag(mnmol),rtmp
      character*2 symbf(mnat)

c GET TRAJECTORY INDEX
      if (tflag(2).eq.1) then
        index = trajlist(itraj)
      else
        index = itraj
      endif

      nprod=2
      nfrag(1) = 0
      nfrag(2) = 0
      mmtot = 0.d0
      mmfragtot(1) = 0.d0
      mmfragtot(2) = 0.d0
      do i=1,3
      comx(i)=0.d0
      comp(i)=0.d0
      comxag(i,1)=0.d0
      comxag(i,2)=0.d0
      compag(i,1)=0.d0
      compag(i,2)=0.d0
      enddo
      do i=1,nat
        if (nprod.eq.1) then       ! for undissociated complexes, assign all atoms to the first group
          ifrag = 1
        else
          ifrag=1
          if (i.gt.natom(1)) ifrag=2
        endif
c        write(6,*)"Atom #",i," ",symbol(i)," is in fragment",ifrag
        nfrag(ifrag) = nfrag(ifrag) + 1
        fragind(nfrag(ifrag),ifrag) = i
        mmfragtot(ifrag) = mmfragtot(ifrag) + mm(i)   ! total fragment mass
        mmtot = mmtot + mm(i)                         ! total mass
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
      do j=1,nprod
      do k=1,3
        comxag(k,j)=comxag(k,j)/mmfragtot(j)-comx(k)          ! center of mass for fragment
        compag(k,j)=compag(k,j)-comp(k)/mmtot*mmfragtot(j)    ! overall momentum for fragment
      enddo
      enddo
c      write(6,*)

c      do j=1,nprod
c      write(6,*)"Product group ",j," has ",nfrag(j),
c     & " atoms with a mass ",mmfragtot(j)/amutoau," amu",natom(1),nprod
c      enddo
c      write(6,*)

      if (nprod.eq.2) then
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
      tdelay = time - rf/dsqrt(erel*2.d0/relmu)
      tmp=dsqrt(compag(1,2)**2+compag(2,2)**2+compag(3,2)**2)
      theta=dacos(compag(3,2)/tmp)/dacos(-1.d0)*180.d0   ! Scattering angle. Assumes Pinitial = (0,0,Prel) for frag 2.
c      write(6,*)"RELATIVE"
c      write(6,*)"Final CoM separation           = ",rf*autoang," A"
c      write(6,*)"Final orbital angular momentum = ",lorbtot," hbar"
c      write(6,*)"Final relative orbital energy  = ",eorb*autoev," eV"
c      write(6,*)"Final relative trans energy    = ",erel*autoev," eV"
c      write(6,*)"Delay time                     = ",tdelay*autofs," fs"
c      write(6,*)"Scattering angle               = ",theta
c      write(6,*)
      else
      lorbtot = 0
      eorb = 0.d0
      erel = 0.d0
      bigjtotag(2) =0.d0
      erottotag(2) = 0.d0
      evibag(2) = 0.d0
      theta = 0.d0
      endif
c      write(6,*)"SUPERMOL"
c      write(6,*)"Final kinetic energy           = ",ke*autoev," eV"
c      write(6,*)"Final potential energy         = ",pe*autoev," eV"
c      write(6,*)"Final total energy             = ",te*autoev," eV"
c      write(6,*)

c      write(6,*)"FRAGMENTS"
      do ifrag=1,nprod
c        write(6,*)"Analyzing fragment ",ifrag
        if (nfrag(ifrag).eq.0) then
c          write(6,*)"No atoms in this fragment"
          peag(ifrag)=0.d0
        else if (nfrag(ifrag).eq.1) then
          itmp = fragind(1,ifrag)
c          write(6,*)" Atomic fragment"
c          write(6,*)" Atom # ",itmp," (",symbol(itmp),")"
          bigjtotag(ifrag)=0.d0
          erottotag(ifrag)=0.d0
          evibag(ifrag)=0.d0
          peag(ifrag)=0.d0
        else
c          write(6,*)"Polyatomic (N=>2) fragment"
          nf=nfrag(ifrag)
          do i=1,nf          ! reassign mass, xx, pp info for this AG to a single array
            itmp = fragind(i,ifrag)
            mmfrag(i) = mm(itmp)
            symbf(i)=symbol(itmp)
            do k=1,3
              ppfrag(k,i) = pp(k,itmp)
              xxfrag(k,i) = xx(k,itmp) - comxag(k,ifrag)
            enddo
          enddo
          if (ldofrag) then  ! TEMP
c          if (.false.) then
            call getpem(xxfrag,nf,pemaf,pemdf,gpemaf,gpemdf,df,symbf)
            if (repflag.eq.1) then
               peag(ifrag)=pemaf(nsurf)
            else
               peag(ifrag)=pemdf(nsurf,nsurf)
            endif
          else
            peag(ifrag)=0.d0
          endif
          peag(ifrag)=peag(ifrag)+ezero-ezeroim(ifrag)
c          write(6,*)"Final PE for this fragment = ",peag(ifrag)*autoev
          call ange(xxfrag,ppfrag,mmfrag,nfrag(ifrag),eig,
     &    tmp3j,bigjtotag(ifrag),
     &    tmp3e,erottotag(ifrag))
c          write(6,*)"Final angular momentum for frag ",ifrag,
c     &              " = ",bigjtotag(ifrag)
c          write(6,*)"Final rotational energy for frag ",ifrag,
c     &              " = ",erottotag(ifrag)*autoev
          do i=1,nfrag(ifrag)
          do k=1,3
            ppfrag(k,i)=ppfrag(k,i)
     &      -compag(k,ifrag)*mmfrag(i)/mmfragtot(ifrag)
          enddo
          enddo
          call gettemp(ppfrag,mmfrag,nfrag(ifrag),
     &                               tempfrag,kefrag(ifrag))
          evibag(ifrag)=kefrag(ifrag)-erottotag(ifrag)
c          write(6,*)"Final vib energy for frag ",ifrag,
c     &              " = ",evibag(ifrag)*autoev
c          write(6,*)"Final kinetic energy for frag ",ifrag,
c     &              " = ",kefrag(ifrag)*autoev
c          write(6,*)"Final total energy for frag ",ifrag,
c     &              " = ",(peag(ifrag)+kefrag(ifrag))*autoev
        endif
c        write(6,*)
      enddo


 88   continue

      if (lwrite(39)) then
        write(39,1031)index,nsurf,outcome,
     &  time*autofs,
     &  bqci*autoang,
     &  pe*autoev,
     &  lorbtot,eorb*autoev,erel*autoev,
     &  bigjtotag(1),erottotag(1)*autoev,evibag(1)*autoev,
     &  peag(1)*autoev,
     &  bigjtotag(2),erottotag(2)*autoev,evibag(2)*autoev,
     &  peag(2)*autoev,
     &  sampjtemp(1),theta
      endif
 1031 format(3i5,100f15.5)

c **********************************************************************

      return
      end
