      subroutine finalstate

c Computes some quantities after each trajectory finishes.
c Writes some output files.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'

      integer i,j,k,arr,index,nprod,ioutatom
      double precision rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),tmp,tmpi,
     & rr(3),jp(9),jq(9),jm(9),xv,xj,rin,rout,evib,erotx,theta,
     & ediatot,kerel,kedia,xjx,xjy,xjz,eint,mmtot

c     temp
      integer wellindex1,wellindex2
      double precision xwell,ptmp,etmp,etmp2,integrand

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

c GENERAL ANALYSES
c     check integration
      write(6,*)"Energy conservation"
      write(6,100)"Initial total energy = ",tei*autoev," eV"
      write(6,100)"Final total energy   = ",te*autoev," eV"
      write(6,101)"Converged to         = ",(te-tei)*autoev," eV"
      write(6,*)
 100  format(1x,a,f18.8,a)
 101  format(1x,a,e18.8,a)

      write(6,*)"Angular momentum conservation"
      write(6,100)"Initial angular mom. = ",bigjtoti," au"
      write(6,100)"Final angular mom.   = ",bigjtot," au"
      write(6,101)"Converged to         = ",(bigjtot-bigjtoti)," au"
      write(6,*)

      write(6,*)"Electronic state population conservation"
      call getrho(crei,cimi,rhor,rhoi,nsurft)
      tmpi = 0.d0
      do i=1,nsurft
       tmpi = tmpi + rhor(i,i)
      enddo
      write(6,102)"                       ","total",(i,i=1,nsurft)
      write(6,103)"Initial populations  = ",tmpi,(rhor(i,i),i=1,nsurft)
      call getrho(cre,cim,rhor,rhoi,nsurft)
      tmp = 0.d0
      do i=1,nsurft
       tmp = tmp + rhor(i,i)
      enddo
      write(6,103)"Final populations    = ",tmp,(rhor(i,i),i=1,nsurft)
      write(6,104)"Converged to         = ",(tmp-tmpi)
      write(6,*)
 102  format(1x,a,10x,a,100i15)
 103  format(1x,a,100f15.5)
 104  format(1x,a,e18.8)

      if (lwrite(36)) write(36,1036)index,nsurf,outcome,
     &   (te-tei)*autoev,(bigjtot-bigjtoti),(tmp-tmpi)
 1036 format(1x,3i5,10f15.5)

c WRITE FINAL STATE INFO
      write(6,*)"Initial electronic state  =",nsurf0
      write(6,*)"Final electronic state    =",nsurf
      write(6,*)

c ATOM-DIATOM INITIAL CONDITIONS, PERFORM THE FOLLOWING ANALYSES, WHICH FOLLOW NAT
      if (initx(1).eq.3) then  !!!!! TURN OFF FOR CO2
      if (outcome.ne.0) then
        call rtox(xx,rr,1)
        if (rr(1).lt.rr(2).and.rr(1).lt.rr(3)) arr=1
        if (rr(2).lt.rr(3).and.rr(2).lt.rr(1)) arr=2   ! bugfix
!        if (rr(2).lt.rr(3).and.rr(2).lt.rr(3)) arr=2
        if (rr(3).lt.rr(1).and.rr(3).lt.rr(2)) arr=3
        write(6,*)"Initial arrangement       =",arrad
        write(6,*)"Final arrangement         =",arr
        write(6,*)
       call carttojac(xx,pp,mm,jq,jp,jm,0,arr)
       kerel = 0.d0
       kedia = 0.d0
       do i=1,3
         kedia = kedia + jp(i)**2/(2.d0*jm(i)) ! this quantity oscillates as the diatom vibrates and rotates
         kerel = kerel + jp(i+3)**2/(2.d0*jm(i+3)) ! this quantity converges as the fragments separate
       enddo
       erottot = te - kedia - kerel - pe ! should be near zero for total angular momentum zero
       ediatot = pe + kedia ! total diatomic energy
c       print *,erottot*autoev,ediatot*autoev,"hi"
c       print *,'kerel=',kerel*autoev
!       do i=1,3
!       easym(1,i)=1.92539/autoev   ! IMLS CO2
!       easym(2,i)=-0.01126/autoev
!       easym(3,i)=-0.01626/autoev
!      enddo
       eint = ediatot - easym(nsurf,arr) ! set zero of energy at the classical equilibrium
       xjx = jq(2)*jp(3)-jq(3)*jp(2)
       xjy = jq(3)*jp(1)-jq(1)*jp(3)
       xjz = jq(1)*jp(2)-jq(2)*jp(1)
       xj = dsqrt(xjx*xjx+xjy*xjy+xjz*xjz) ! rotational quantum number is magnitude of RxP
       xj = xj-0.5d0 ! Bohr-Sommerfeld conditions
       xj = max(0.d0,xj)
       erotx = (xj+0.5d0)**2/(2.d0*rasym(nsurf,arr)**2*jm(1)) ! assume separability
!       rtmp = 1.1293/autoang  ! IMLS CO2
!       erotx = (xj+0.5d0)**2/(2.d0*rtmp**2*jm(1)) ! assume separability
       evib = eint - erotx
       call vwkb(arr,mm,ediatot,xj,xv,rin,rout,nsurf) ! compute vibrational action using WKB
c       print *,rin*autoang,rout*autoang
       eint = eint*autoev
       kerel = kerel*autoev
       erotx = erotx*autoev
       evib = evib*autoev
       etmp=easym(nsurf,arr)
      else
      te=0.d0
      arr=0
      etmp=0.d0
      kerel=0.d0
      eint=0.d0
      evib=0.d0
      erotx=0.d0
      xv=0.d0
      xj=0.d0
      endif
       if (lwrite(30)) write(30,115)index,nsurf,arr,time*autofs,istep,
     &            rhor(1,1),rhor(2,2),
     &            te*autoev,etmp*autoev,
     &            kerel,eint,evib,erotx,xv,xj,theta
       etmp2=eint
c NOTE:  XV is not accurate for SE method, because a pure state is used
c        instead of the mixed state in VWKB!
 115   format(3i5,f15.5,i10,100f10.5)
      endif   !!! turn off for CO2



c **********************************************************************
C FOR ASSOCIATION/DISSOCIATION REACTIONS
      if (termflag.eq.3) then
c **********************************************************************

      write(6,*)"Final state analysis for TERMFLAG = 3"
      write(6,*)"Outcome = ",outcome

      if (outcome.eq.-1) then
c     somehow the trajectory failed
c     don't analyze
      write(6,*)"Outcome = -1 indicates a problem with this trajectory"
      write(6,*)"Skipping analysis"
      write(6,*)
      tdelay = 0.
      go to 88
      endif

      if (outcome.eq.-2) then
c     this trajectory could not be photoexcited
c     don't analyze
      write(6,*)"Outcome = -2 "
      write(6,*)"Skipping analysis"
      write(6,*)
      tdelay = 0.
      go to 88
      endif

      if (outcome.eq.0) then
c     we hit the max steps before we dissociated
        write(6,*)"Hit the maximum number of steps. ",
     &  "Products will be analyzed as a single molecule."
        nprod = 1
      else
        if (t_r(outcome).gt.0) then
        write(6,*)"Dissociation"
        write(6,*)"Broke a ",t_symb(outcome,1),"-",t_symb(outcome,2),
     &  " bond with atom indices",aind(outcome,1),
     &  " and ",aind(outcome,2)
        write(6,*)"Assigning fragment # based on shortest distance to ",
     &  "the atoms in the broken bond"
        nprod = 2
        else
        write(6,*)"Association"
        write(6,*)"Formed a ",t_symb(outcome,1),"-",t_symb(outcome,2),
     &  " bond with atom indices",aind(outcome,1),
     &  " and ",aind(outcome,2)
        write(6,*)"Assigning all atoms to the same product fragment"
        nprod = 1
        endif
      endif
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
          if(t_symb(outcome,1).eq."cm".and.
     &       t_symb(outcome,2).eq."cm") then
          ifrag=1
          if (i.gt.natom(1)) ifrag=2
          else
          rr1 = 0.d0               ! for dissociated fragments, assign atoms to group based on shortest distance to atoms in broken bond
          rr2 = 0.d0
          do j=1,3
            rr1 = rr1 + (xx(j,i)-xx(j,aind(outcome,1)))**2
            rr2 = rr2 + (xx(j,i)-xx(j,aind(outcome,2)))**2
          enddo
          ifrag = 1
          if (rr2.lt.rr1) ifrag = 2
          endif
        endif
c H2O+2AR HACK TEMP AJ
c        if (i.le.3) ifrag=1
c        if (i.gt.3) ifrag=2
c H2O+2AR HACK TEMP AJ
        write(6,*)"Atom #",i," ",symbol(i)," is in fragment",ifrag
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
      write(6,*)

      do j=1,nprod
      write(6,*)"Product group ",j," has ",nfrag(j),
     &          " atoms with a mass ",mmfragtot(j)/amutoau," amu"
      if (nfrag(j).eq.1.and.mmfragtot(j).eq.15.9949d0) print *,"REACT"
      enddo
      write(6,*)

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
      ioutatom=2
c      ioutatom=2  ! hack for CO2
c      if(outcome.eq.2) ioutatom=1  ! hack for CO2
      tmp=dsqrt(compag(1,ioutatom)**2
     & +compag(2,ioutatom)**2+compag(3,ioutatom)**2)
      theta=dacos(compag(3,ioutatom)/tmp)/dacos(-1.d0)*180.d0   ! Scattering angle. Assumes Pinitial = (0,0,Prel) for frag 2.
      write(6,*)"RELATIVE"
      write(6,*)"Final CoM separation           = ",rf*autoang," A"
      write(6,*)"Final orbital angular momentum = ",lorbtot," hbar"
      write(6,*)"Final relative orbital energy  = ",eorb*autoev," eV"
      write(6,*)"Final relative trans energy    = ",erel*autoev," eV"
      write(6,*)"Delay time                     = ",tdelay*autofs," fs"
      write(6,*)"Scattering angle               = ",theta
      write(6,*)
c uncomment for CO2
c       if (lwrite(30)) write(30,115)index,nsurf,arr,time*autofs,istep,
c     &            rhor(1,1),rhor(2,2),
c     &            te*autoev,etmp*autoev,
c     &            kerel,etmp2,evib,erotx,xv,xj,theta
      else
      lorbtot = 0
      eorb = 0.d0
      erel = 0.d0
      bigjtotag(2) =0.d0
      erottotag(2) = 0.d0
      evibag(2) = 0.d0
      theta = 0.d0
      endif
      write(6,*)"SUPERMOL"
      write(6,*)"Final kinetic energy           = ",ke*autoev," eV"
      write(6,*)"Final potential energy         = ",pe*autoev," eV"
      write(6,*)"Final total energy             = ",te*autoev," eV"
      write(6,*)

      write(6,*)"FRAGMENTS"
      do ifrag=1,nprod
        write(6,*)"Analyzing fragment ",ifrag
        if (nfrag(ifrag).eq.0) then
          write(6,*)"No atoms in this fragment"
          peag(ifrag)=0.d0
        else if (nfrag(ifrag).eq.1) then
          itmp = fragind(1,ifrag)
          write(6,*)" Atomic fragment"
          write(6,*)" Atom # ",itmp," (",symbol(itmp),")"
          bigjtotag(ifrag)=0.d0
          erottotag(ifrag)=0.d0
          evibag(ifrag)=0.d0
          peag(ifrag)=0.d0
c        else if (nfrag(ifrag).eq.2) then
c maybe include special analysis for diatoms, but for now use poly analyses
        else if (.false.) then
          itmp1 = fragind(1,ifrag)
          itmp2 = fragind(2,ifrag)
          write(6,*)" Diatomic fragment"
          write(6,*)" Atom #s ",itmp1," (",symbol(itmp1),") and ",itmp2
     & ," (",symbol(itmp2),")"
        else
          write(6,*)"Polyatomic (N=>2) fragment"
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
          write(6,*)"Final PE for this fragment = ",peag(ifrag)*autoev
          call ange(xxfrag,ppfrag,mmfrag,nfrag(ifrag),eig,
     &    tmp3j,bigjtotag(ifrag),
     &    tmp3e,erottotag(ifrag))
          write(6,*)"Final angular momentum for frag ",ifrag,
     &              " = ",bigjtotag(ifrag)
          write(6,*)"Final rotational energy for frag ",ifrag,
     &              " = ",erottotag(ifrag)*autoev
          do i=1,nfrag(ifrag)
          do k=1,3
            ppfrag(k,i)=ppfrag(k,i)
     &      -compag(k,ifrag)*mmfrag(i)/mmfragtot(ifrag)
          enddo
          enddo
          call gettemp(ppfrag,mmfrag,nfrag(ifrag),
     &                               tempfrag,kefrag(ifrag))
          evibag(ifrag)=kefrag(ifrag)-erottotag(ifrag)
          write(6,*)"Final vib energy for frag ",ifrag,
     &              " = ",evibag(ifrag)*autoev
          write(6,*)"Final kinetic energy for frag ",ifrag,
     &              " = ",kefrag(ifrag)*autoev
          write(6,*)"Final total energy for frag ",ifrag,
     &              " = ",(peag(ifrag)+kefrag(ifrag))*autoev
        endif
        write(6,*)
      enddo


 88   continue

c hack delay time
      if (lwrite(32)) write(32,1032)index,nsurf,outcome,
     &   rhor(1,1),rhor(2,2),tdelay*autofs
 1032 format(3i5,3f15.5)


      if (lwrite(31)) then
        write(31,1031)index,nsurf,outcome,
     &  time*autofs,
     &  atemp/time,stemp,
     &  bqci*autoang,lqci,eorbqci*autoev,erelqci*autoev,
     &  pei*autoev,
     &  bigjtotagi(1),erottotagi(1)*autoev,evibagi(1)*autoev,
     &  peagi(1)*autoev,
     &  bigjtotagi(2),erottotagi(2)*autoev,evibagi(2)*autoev,
     &  peagi(2)*autoev,
     &  pe*autoev,
     &  lorbtot,eorb*autoev,erel*autoev,
     &  bigjtotag(1),erottotag(1)*autoev,evibag(1)*autoev,
     &  peag(1)*autoev,
     &  bigjtotag(2),erottotag(2)*autoev,evibag(2)*autoev,
     &  peag(2)*autoev,
     &  sampjtemp(1),theta
      endif

      if (lwrite(38)) then
        write(38,1038)index,bqci*autoang,(eorbqci+erelqci)*autoev,
     &  (eorb+erel)*autoev,pei*autoev,pe*autoev,theta
      endif

 1031 format(3i5,100f15.5)
 1038 format(i5,100d15.5)

      endif
c **********************************************************************

c      integrand=bqci*erelqci/tempqc/kb
c     &   *(1.d0-dcos(theta/180.d0*dacos(-1.d0)))
c      write(6,'(10f15.5)')bqci,erelqci/tempqc/kb,theta,integrand

      return
      end
