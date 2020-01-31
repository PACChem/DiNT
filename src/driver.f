      subroutine driver

c Computes a full trajectory from start to finish.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
c hack
      common/tmp/tmpprint
      double precision tmpprint(50)
      double precision taup(3),taup2(3)
      common/tauprint/taup,taup2
      double precision h12x,h12y
      common/pesh12/h12x,h12y

c     number of bins for binning the radial distribution function
      integer nbinrad,im
      parameter(nbinrad=60)

      integer i,j,k,ii,iturn,i1,i2,ithistraj,frusflag(mnsurf),
     &  newsurf,iramp,ncross
      double precision rminmax,ravg,ce,hrad,dot,peim
      double precision rr(3),rhor(mnsurf,mnsurf),eig(3),dum3(3),
     &  rhoi(mnsurf,mnsurf),step,timeprev,phop(mnsurf),tmp1,tmp2,tmp3,
     &  dmag,dmag1,dmag2,gvmag,rrr,rmax(mnoutcome),rmin(mnoutcome),
     &  dvec(3,mnat,mnsurf,mnsurf),gpem(3,mnat,mnsurf),
     &  arij(mnat,mnat),arij2(mnat,mnat),lind,
     &  raddist(0:nbinrad+1),egap,rtmp,rcom,tmp4,ecom,
     &  lastgap0,lastgap,plz,plz0,elz,plzx,elz0,elzx,lastcross,
     &  hlz,hlz0,hlzx,glz,glzx,glz0,r1,r10,r1x,r2,r20,r2x, 
     7  a123,a1230,a123x,r3,vlz,vlz0,vlzx,pairy(2),pairy0(2),pairyx(2)
      logical lhop

c used for FSTU
      double precision tu_xx(3,mnat,mnsurf),tu_pp(3,mnat,mnsurf),
     & tu_cre(mnsurf,mnsurf),tu_cim(mnsurf,mnsurf),tu_hstep(mnsurf),
     & tu_time(mnsurf),tu_phase(mnsurf,mnsurf),tu_maxt(mnsurf),
     & deltforw,deltback,tu_maxtime
      integer tu_istep(mnsurf)
      logical lfrust

c Andersen thermostat
      logical lhit

c stochastic decoherence
      integer nhop
      double precision stodecotime,stodecotau,electau

      integer wellindex1,wellindex2,n1,n2,ijk
      double precision xwell,rvdw,evdw
      logical lzflag
      lzflag=.false.
      if (tflag(4).lt.0) lzflag=.true.
c      rvdw=0.d0
      evdw=1000.d0

c initialize for this trajectory
      lastgap0 =0.d0 ! TEMP AJ
      lastgap =0.d0 ! TEMP AJ
      lastcross =0.d0 ! TEMP AJ
      ncross =0 ! TEMP AJ
      nhop = 0
      time = 0.d0
      stodecotime = -1.d0
      stodecotau = 4.d0/autofs
      iturn = 0
      dmag1 = 0.d0
      dmag2 = 0.d0
      do i=1,nsurft
        phase(i) = 0.d0
        tu_time(i) = -1.d0
      enddo
      lfrust = .false.
      do i=1,mnat
      do j=1,mnat
        arij(i,j) = 0.d0
        arij2(i,j) = 0.d0
      enddo
      enddo
      atemp = 0.d0
      atemp2 = 0.d0
      stemp = 0.d0
      iramp = 0
      do i=0,nbinrad+1
        raddist(i) = 0.d0
      enddo
      outcome = -1   ! OUTCOME should get reassigned when the traj finishes. A final value of -1 indicates a failure of some kind.
      istep = 0
      istepw = 0
      step = 0.d0

c PIMD/RPMD add initialization

c calculate and save initial values of things for later analysis
c     initial coordinates and momenta
      do i=1,nat
      do j=1,3
         xxi(j,i)=xx(j,i)
         ppi(j,i)=pp(j,i)
      enddo
      enddo
      call gettemp(pp,mm,nat,tempi,kei)
      call getgrad(xx,pp,nsurf,pei,cre,cim,gv,gcre,gcim,nat,
     & phop,dmag,dvec,pem,gpem,phase)
      lastgap0 = lastgap  ! TEMP AJ
      lastgap = pem(1)-pem(2)  ! TEMP AJ
      call getplz(pp,mm,gpem,nat,plz,elz,hlz,glz,pairy,nsurf)  ! TEMP AJ

c if initx=2 or 5 and scale0im>0 scale momenta to etot=scale0im
c      do im=1,nmol
c      if (ldofrag) then
c         peim=peagi(im)
c      else
c         peim=pei
c      endif
c      if (initx(im).eq.2.or.initx(im).eq.5) then 
c        if (scale0im(im).gt.0.d0.and.pei.lt.scale0im(im)) then
c          do i=1,natom(im)
c          ii=i+iatom(im)
c          do j=1,3
c             pp(j,ii)=pp(j,ii)-compp(j,im)*mm(ii)/mmag(im)
c             pp(j,ii)=pp(j,ii)*dsqrt((scale0im(im)-peim)/(kei-erelqc))
c             pp(j,ii)=pp(j,ii)+compp(j,im)*mm(ii)/mmag(im)
c          enddo
c          enddo
c          call gettemp(pp,mm,nat,tempi,kei)
c          call getgrad(xx,pp,nsurf,pei,cre,cim,gv,gcre,gcim,nat,
c     &     phop,dmag,dvec,pem,gpem,phase)
c      else if (scale0im(im).gt.0.d0.and.pei.gt.scale0im(im)) then
c        ldriver=0
c        go to 1001
c      endif
c      endif
c      enddo

      tei=kei+pei
      call ange(xx,pp,mm,nat,eig,bigji,bigjtoti,eroti,erottoti)
      write(6,106)"Initial overall ang momentum = ",bigjtoti," au"
      write(6,106)"Initial overall rot energy   = ",erottoti*autoev,
     & " eV"
      write(6,106)"Initial kinetic energy       = ",kei*autoev," eV"
      write(6,106)"Initial potential energy     = ",pei*autoev," eV"
      write(6,106)"Initial total energy         = ",tei*autoev," eV"
      write(6,106)"Initial temp                 = ",tempi," K"
 106  format(1x,a,f15.5,1x,a)
 107  format(1x,a,3f12.5,1x,a)
      write(6,*)

      if(nsurft.gt.1) then
      do i=1,nsurft
      write(6,*)"Potential for surf ",i," = ",pem(i)*autoev," eV"
      enddo
      endif
      write(6,*)

c photoexcite
      if (tflag(3).eq.1) then
        egap = dabs(pem(ntarget)-pem(nsurf))
        if (dabs(egap-ephoton).lt.wphoton) then
c         successful excitation
          print *,"Exciting from state ",nsurf," to state ",ntarget,
     &               " with a photon energy of ",egap*autoev
          nsurf=ntarget
          pei = pem(nsurf)
          tei = kei + pei
          call initelec
        else
          print *,"Unsuccessful excitation"
          print *,"Egap = ",egap*autoev," eV"
          print *,"Photon energy = ",ephoton*autoev," +/- ",
     &             wphoton*autoev," eV"
          outcome = -2
          go to 999
        endif
      endif

c get trajectory index
      if (tflag(2).eq.1) then
        ithistraj = trajlist(itraj)
      else
        ithistraj = itraj
      endif

c write initial coordinates and momenta to output
      if (lwrite(10)) then
      open(10)
      write(10,*)ithistraj,pei*autoev
      do i=1,nat
      write(10,1010)symbol(i),mm(i)/amutoau,(xx(j,i)*autoang,j=1,3)
      enddo
      endif
      if (lwrite(11)) then
      open(11)
      write(11,*)ithistraj,kei*autoev
      do i=1,nat
      write(11,1010)symbol(i),mm(i)/amutoau,(pp(j,i),j=1,3)
      enddo
      endif

 1010 format(5x,a2,4f20.8)

c trajectory output
      if (termflag.eq.2) then
      write(6,'(20a15)')"          step",
     &                  "      time(fs)",
     &                  "        PE(eV)",
     &                  " |gradV|(eV/A)"
      else
      if (nsurft.eq.1.and.lwell.le.0) then
      write(6,'(20a15)')"           step",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "      Lindemann",
     &                  "         <T(K)>",
     &                  "   st.dev(T(K))"
      else if (nsurft.eq.1.and.lwell.gt.0) then
      write(6,'(20a15)')"           step",
     &                  "          stepw",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "      Lindemann",
     &                  "         <T(K)>",
     &                  "   st.dev(T(K))"
      else
      write(6,'(20a15)')"           step",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "          Probs"
      endif
      endif

c ######################################################################
c INTEGRATE TRAJECTORY
c ######################################################################
      do 10 while (.true.)
      istep = istep + 1      

c      call state
 
c check for trajectory termination here
c when termination occurs we skip out of the infinite 
c "do 10 while" loop with "go to 999"

c     quit if we exceed max number of steps
      if (lwell.le.0.and.istep.gt.t_nstep.and.t_nstep.ne.-1) then
        outcome = 0
        go to 999
      else if (lwell.gt.0.and.istepw.ge.t_nstep.and.t_nstep.ne.-1) then
        outcome = 0
        go to 999
      endif

c     quit if we've integrated long enough (TERMFLAG = 1)
      if (termflag.eq.1.and.time.gt.t_stime) then
        outcome = 0
        go to 999
      endif

c     converge gradient (TERMFLAG = 2)
      if (termflag.eq.2) then
        gvmag = 0.d0
        do i=1,nat
        do j=1,3
        gvmag = gvmag + gv(j,i)**2
        enddo
        enddo
        if (gvmag.le.t_gradmag**2) then
          write(6,*)"Gradient converged to ",
     &               dsqrt(gvmag)*autoev/autoang,"eV/A"
          write(6,*)"At the geometry "
          do i=1,nat
            write(*,1090)symbol(i),mm(i)/amutoau,(xx(j,i)*autoang,j=1,3)
 1090       format(a2,f15.6,3f21.8)
          enddo
          tmp1=-1.d0
          tmp2=-1.d0
          call addrot(pp,xx,mm,nat,tmp1,tmp2,sampjbrot1,sampjbrot2,
     &       tmp3,tmp4,dum3)

          if (lwrite(77)) then  ! special output file
            open(77,file="dint.geo") 
            write(77,*)pem(1)*autoev,
     &" ! Optimized energy (eV), rot constants (cm-1), and geometry (A)"
            write(77,*)dum3
            do i=1,nat
           write(77,1090)symbol(i),mm(i)/amutoau,(xx(j,i)*autoang,j=1,3)
            enddo
          endif

          go to 999
        endif
      endif

c     monitor bond dissociation (TERMFLAG = 3)
      if (termflag.eq.3) then
c       find maximum and minimum bond distance corresponding to each outcome
        do k=1,t_noutcome
        rmax(k) = 0.d0
        rmin(k) = 1000.d0
        enddo
        do i1=1,nat
        do i2=i1+1,nat
          do k=1,t_noutcome
            if ( ( t_symb(k,1).eq.symbol(i1).and.
     &             t_symb(k,2).eq.symbol(i2)      ) .or.
     &           ( t_symb(k,1).eq.symbol(i2).and.
     &             t_symb(k,2).eq.symbol(i1)      ) ) then
              rrr = 0.d0
              do j=1,3
              rrr = rrr + (xx(j,i1)-xx(j,i2))**2
              enddo
              rrr = dsqrt(rrr)
c              if (rmax(k).lt.rrr.and.t_r(k).gt.0.d0) then
              if (rmax(k).lt.rrr) then
                rmax(k) = rrr
                if (t_r(k).gt.0.d0) then
                aind(k,1) = i1
                aind(k,2) = i2
                endif
              endif
c              if (rmin(k).gt.rrr.and.t_r(k).lt.0.d0) then
              if (rmin(k).gt.rrr) then
                rmin(k) = rrr
                if (t_r(k).lt.0.d0) then
                aind(k,1) = i1
                aind(k,2) = i2
                endif
              endif
c              print *,i1,i2,rrr*autoang,rmin(k)*autoang,rmax(k)*autoang
            endif
            if ( ( t_symb(k,1).eq.'cm'.and.
     &             t_symb(k,2).eq.'cm'      ) ) then
              call getrel(rcom,ecom)
              rmin(k)=rcom
              rmax(k)=rcom
            endif
          enddo
        enddo
        enddo

c       ITURN starts each trajectory at 0
c       ITURN is set to 1 when all the distances are less than 95% of the
c             dissociation termination distances
c       This prevents the trajectories from "dissociating" before the collision
        if (iturn.eq.0) then
          iturn = 1
c         print *,rmax(1)*autoang,t_r(1)*autoang,tmpprint(2)*autoang
          do k=1,t_noutcome
            if (rmax(k).gt.0.95d0*t_r(k).and.t_r(k).gt.0.d0) iturn = 0
          enddo
        endif
c       monitor dissociation
c          iturn = 1   ! TEMP HACK AJ
        if (iturn.eq.1) then
          do k=1,t_noutcome
            rtmp=rmax(k)
            if (t_symb(k,2).eq.'x') rtmp=rmin(k)
c            if (k.eq.1) rtmp=rmin(k)  ! TEMP HACK AJ
c       print *,k,rmin(k)*autoang,rmax(k)*autoang,iturn,t_r(k)*autoang
            if (rtmp.gt.t_r(k).and.t_r(k).gt.0.d0) then
c           intercept FSTU trajectories that are still looking to unfrustrate
            if (lfrust) then
            else
              outcome = k
              go to 999
            endif
            endif
          enddo
        endif
c       monitor association
        do k=1,t_noutcome
          if (rmin(k).lt.dabs(t_r(k)).and.t_r(k).lt.0.d0) then
c TMP AJ for CO2 association
c            call getrho(cre,cim,rhor,rhoi,nsurft)
c            if (rhor(1,1).gt.0.99d0) then   ! CSDM or FS/SD
c            if (nsurf.eq.1) then    ! FS
c TMP AJ
            outcome = k
            go to 999
c            endif
          endif
        enddo
      endif

c Andersen thermostat
      lhit=.false.
      if (tflag(1).eq.3) call andersen(pp,mm,nat,step,lhit)
      if (lhit.and.scandth.gt.0.d0) then  
        call gettemp(pp,mm,nat,temp,ke)
        do i=1,nat
        do j=1,3
          pp(j,i)=pp(j,i)*dsqrt((scandth+pei-pe)/ke)
        enddo
        enddo
        call gettemp(pp,mm,nat,temp,ke)
      endif

c     update information if we changed momentum
      if (lhit) call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     & phop,dmag,dvec,pem,gpem,phase)

c take step
      timeprev = time
      call takestep
      step = time - timeprev
c      if (step*autofs.lt.1.d-2) then
c         outcome=0
c         go to 999
c      endif
      if (lwrite(85)) write(85,*)
      if (lwrite(85)) write(85,*)ithistraj,istep,time*autofs

c transport
c      rtmp=xx(1,1)**2+xx(2,1)**2+xx(3,1)**2
c      rtmp=rtmp/6.d0
c      rtrans=(rtrans*(time-step)+rtmp*step)/time
c      if (mod(istep,1000).eq.0) print *,"D = ",rtrans/time

c ramp temperature
      if (tflag(1).eq.2) then
c hard-coded parameters
c      ramptime = 2.d4 ! in fs
c      nramp = 10 ! number of ramping cycles
c      rampfact = 1.02 ! increase temp by this factor each ramp
      if (time > ramptime) then
      call lindemann(xx,nat,step,time,arij,arij2,lind)
      call radialdist(ithistraj,xx,nat,step,time,nbinrad,raddist,1)
      if (lwrite(42)) call honey(ithistraj,xx,nat,time)
      iramp = iramp + 1
      if (lwrite(40)) write(40,*)ithistraj,iramp,atemp/time,stemp,lind
 1040 format(2i7,3f15.5)
      if (iramp.eq.nramp) go to 999
      do i=1,mnat
      do j=1,mnat
        arij(i,j) = 0.d0
        arij2(i,j) = 0.d0
      enddo
      enddo
      do i=0,nbinrad+1
        raddist(i) = 0.d0
      enddo
      atemp = 0.d0
      atemp2 = 0.d0
      stemp = 0.d0
      time = step
      do i=1,nat
      do j=1,3
        pp(j,i)=pp(j,i)*rampfact
      enddo
      enddo
      endif
      endif

c stochastic decoherence
      if (stodecoflag) then
      call stodeco(ithistraj,istep,cre,cim,pem,nsurf,time,
     &   stodecotime,stodecotau,step)
      endif

c update info at new geometry
      call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     &   phop,dmag,dvec,pem,gpem,phase)
      if (lzflag) then
      r10=r1
      r20=r2
      a1230=a123
      r1=0.d0
      r2=0.d0
      r3=0.d0
      do i=1,3
      r1=r1+(xx(i,1)-xx(i,2))**2
      r2=r2+(xx(i,3)-xx(i,2))**2
      r3=r3+(xx(i,3)-xx(i,1))**2
      enddo
      a123=(r1+r2-r3)/dsqrt(4.d0*r1*r2)
      r1=dsqrt(r1)*autoang
      r2=dsqrt(r2)*autoang
      a123=dacos(a123)/dacos(-1.d0)*180.d0
      lastgap0 = lastgap  ! TEMP AJ
      lastgap = pem(1)-pem(2)  ! TEMP AJ
      vlz0=vlz
      vlz=(pem(1)+pem(2))/2.d0
      pairy0(1)=pairy(1)
      pairy0(2)=pairy(2)
      plz0=plz
      elz0=elz
      hlz0=hlz
      vlz0=vlz
      glz0=glz
      call getplz(pp,mm,gpem,nat,plz,elz,hlz,glz,pairy,nsurf)  ! TEMP AJ
      if (lastgap0*lastgap.lt.0.d0) then
      ncross=ncross+1
      if (ncross.eq.abs(tflag(4))) go to 999
c     crude interpolation
      r1x=(dabs(lastgap0)*r1+dabs(lastgap)*r10)
     &         /(dabs(lastgap)+dabs(lastgap0))
      r2x=(dabs(lastgap0)*r2+dabs(lastgap)*r20)
     &         /(dabs(lastgap)+dabs(lastgap0))
      a123x=(dabs(lastgap0)*a123+dabs(lastgap)*a1230)
     &         /(dabs(lastgap)+dabs(lastgap0))
      plzx=(dabs(lastgap0)*plz+dabs(lastgap)*plz0)
     &         /(dabs(lastgap)+dabs(lastgap0))
      pairyx(1)=(dabs(lastgap0)*pairy(1)+dabs(lastgap)*pairy0(1))
     &         /(dabs(lastgap)+dabs(lastgap0))
      pairyx(2)=(dabs(lastgap0)*pairy(2)+dabs(lastgap)*pairy0(2))
     &         /(dabs(lastgap)+dabs(lastgap0))
      elzx=(dabs(lastgap0)*elz+dabs(lastgap)*elz0)
     &         /(dabs(lastgap)+dabs(lastgap0))
      hlzx=(dabs(lastgap0)*hlz+dabs(lastgap)*hlz0)
     &         /(dabs(lastgap)+dabs(lastgap0))
      vlzx=(dabs(lastgap0)*vlz+dabs(lastgap)*vlz0)
     &         /(dabs(lastgap)+dabs(lastgap0))
      glzx=(dabs(lastgap0)*glz+dabs(lastgap)*glz0)
     &         /(dabs(lastgap)+dabs(lastgap0))
        n1 = 2  ! current state
        n2 = 1  ! target state
        call elecdeco(ithistraj,istep,time,nat,n1,n2,
     &   pem,gpem,dvec,xx,pp,stodecotau)   ! returns STODECOTAU, which is the 
        call getrho(cre,cim,rhor,rhoi,nsurft)
      if (lwrite(16))
     &  write(16,'(2f13.5,2e18.5,3e18.5,i5,100e20.5)')
     &  time*autofs,dmag,phop(1)*step,rhor(1,1),
c     &  (taup(i)*autofs,i=1,3),
     &  (taup2(i)*autofs,i=1,3),
     &  ncross,
     &  (time-lastcross)*autofs,
     &  plzx,elzx*autoev,hlzx*autocmi,glzx,
     &  vlzx*autoev,pairyx(1),pairyx(2)
c     &  r1x,r2x,a123x,vlzx*autoev,pairyx(1),pairyx(2)
!      print *,"CROSS",plzx,pairyx,pairy,pairy0
      lastcross=time
c      ijk=1
      else
        n1 = 2  ! current state
        n2 = 1  ! target state
        call elecdeco(ithistraj,istep,time,nat,n1,n2,
     &   pem,gpem,dvec,xx,pp,stodecotau)   ! returns STODECOTAU, which is the 
      if (lwrite(16)) 
     & write(16,'(2f13.5,2e18.5,3e18.5)')time*autofs,dmag,
     & phop(1)*step,rhor(1,1),
     & (taup2(i)*autofs,i=1,3)
      endif
      endif

      call getrel(rcom,ecom)
      if (lwrite(90)) 
     &  write(90,'(i5,10f15.5)')ithistraj,time*autofs,
     & tmpprint(2)*autoang,
     &  tmpprint(1)*autocmi,
     &  ecom*autocmi,(tmpprint(1)+ecom)*autocmi
      if (tmpprint(1).lt.evdw) then
         evdw=tmpprint(1)
         rvdw=tmpprint(2)
      endif

c Stop at the first DMAG minimum
      if (termflag.eq.4) then
      if( (dmag1-dmag2).lt.0.d0 .AND. (dmag-dmag1).gt.0.d0 ) go to 999
      dmag2 = dmag1
      dmag1 = dmag
      endif

c reinitialize CSDM
      if (methflag.eq.4) then
      if (istep.gt.2) then
        if( (dmag1-dmag2).lt.0.d0 .AND. (dmag-dmag1).gt.0.d0 ) then
          write(6,*)"Electronic reinitialization!"
          do i=1,nsurft
            i2 = i + nsurft
            cre(i2)=cre(i)
            cim(i2)=cim(i)
            gcre(i2)=gcre(i)
            gcim(i2)=gcim(i)
          enddo
c         update information (check this)
          call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     &     phop,dmag,dvec,pem,gpem,phase)
        endif
      endif
      dmag2 = dmag1
      dmag1 = dmag
      endif

c multiply hopping probability by stepsize for both surface hopping
c and DM methods
      phop(nsurf) = 1.d0
      do i=1,nsurft
        if (i.ne.nsurf) then
          phop(i)=phop(i)*step
          phop(i)=max(phop(i),0.d0)
c         Note:  phop can be > 1 if the probabilities are changing too
c         quickly.  But this shouldn't happen too often.  For better
c         results we'd need to integrate the hopping probability so
c         that the integrator takes smaller steps when necessary.
          phop(i)=min(phop(i),1.d0)
c         prob of staying in NSURF = 1 - sum of all hopping probs
          phop(nsurf) = phop(nsurf) - phop(i)
        endif
      enddo

c FSTU, check for frustrated hop
      if (methflag.eq.5) then
        call checkfrus(nat,nsurf,pp,pem,dvec,frusflag,tu_maxt)
        do k=1,nsurft
          if ((.not.lfrust).and.frusflag(k).eq.0) then
c           at this time a hop is allowed, save information
c           for k = nsurf, we save the current info
c           when lfrust = true, we are looking for a hopping
c           location, so we do not save info
            do i=1,nat
            do j=1,3
              tu_xx(j,i,k) = xx(j,i)
              tu_pp(j,i,k) = pp(j,i)
            enddo
            enddo
            do i=1,nsurft
              tu_cre(i,k) = cre(i)
              tu_cim(i,k) = cim(i)
              tu_phase(i,k) = phase(i)
            enddo
            tu_hstep(k) = hstep
            tu_istep(k) = istep
            tu_time(k) = time
          endif
        enddo
      endif

c TFS and FSTU check for a hop
      if (methflag.eq.1.or.(methflag.eq.5.and.(.not.lfrust))) then
        lhop = .false.
        call checkhop(phop,nsurf,lhop,newsurf)
      endif

c intercept a frustrated hop for FSTU
      if (methflag.eq.5.and.lhop) then
        if (frusflag(newsurf).ne.0) then
c       this hop would be frustrated
        lfrust = .true.
        lhop = .false.
        tu_maxtime = tu_maxt(newsurf)
        endif
      endif

c handle frustrated hop with FSTU
      if (methflag.eq.5.and.lfrust) then
        deltback = tu_time(nsurf) - tu_time(newsurf)
        deltforw = time - tu_time(nsurf)
        if (tu_time(newsurf).ge.0.d0.and.deltback.lt.deltforw.and.
     &                                   deltback.le.tu_maxtime) then
c         do a TU hop backward
          write(6,*)"BACKWARD HOP!"
          do i=1,nat
          do j=1,3
            xx(j,i) = tu_xx(j,i,newsurf)
            pp(j,i) = tu_pp(j,i,newsurf)
          enddo
          enddo
          do i=1,nsurft
            cre(i) = tu_cre(i,newsurf)
            cim(i) = tu_cim(i,newsurf)
            phase(i) = tu_phase(i,newsurf)
          enddo
          hstep = tu_hstep(newsurf)
          istep = tu_istep(newsurf)
          time = tu_time(newsurf)
          lhop = .true.
        else if (frusflag(newsurf).eq.0.and.deltforw.le.tu_maxtime) then
c         do a TU hop forward
          write(6,*)"FORWARD HOP!"
          lhop = .true.
        else if (deltforw.gt.tu_maxtime) then
c         do a frustrated hop at the original hopping location
          do i=1,nat
          do j=1,3
            xx(j,i) = tu_xx(j,i,nsurf)
            pp(j,i) = tu_pp(j,i,nsurf)
          enddo
          enddo
          do i=1,nsurft
            cre(i) = tu_cre(i,nsurf)
            cim(i) = tu_cim(i,nsurf)
            phase(i) = tu_phase(i,nsurf)
          enddo
          hstep = tu_hstep(nsurf)
          istep = tu_istep(nsurf)
          time = tu_time(nsurf)
          lhop = .true.
        endif
        if (lhop) then
        lfrust = .false.
        call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     &   phop,dmag,dvec,pem,gpem,phase)
          do i=1,nsurft
c           reset all backward hops
            tu_time(i) = -1.d0
          enddo
        endif
      endif

c hop or do a frustrated hop
      if ((methflag.eq.1.or.methflag.eq.5).and.lhop) then
        n1 = nsurf  ! current state
        n2 = newsurf  ! target state
        if (stodecoflag) then
        stodecotime=time ! save the time of this hopping event in STODECOTIME
        call elecdeco(ithistraj,istep,time,nat,n1,n2,
     &   pem,gpem,dvec,xx,pp,stodecotau)   ! returns STODECOTAU, which is the 
                                           ! characteristic decoherence time computed 
                                           ! using the JCP paper formula
!        write(6,*)"Setting stodecotau to ",stodecotau*autofs," fs"
        endif

        nhop = nhop + 1
c       switch surfaces or reflect/ignore frustrated hop
        call hop(ithistraj,istep,time,nat,nsurf,newsurf,
     &   pem,gpem,dvec,xx,pp,stodecotau)
        call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     &   phop,dmag,dvec,pem,gpem,phase)
c AJAJAJAJ
c Switch commented out parts to stop allowing TU hops up at the last hop down
            do i=1,nat
            do j=1,3
              tu_pp(j,i,n1) = pp(j,i)
              tu_pp(j,i,n2) = pp(j,i)
            enddo
            enddo
c        do i=1,nsurft
c         reset all backward hops
c          tu_time(i) = -1.d0
c        enddo
c AJAJAJAJ
      endif

c for Semiclassical Ehrenfest, change NSURF according to RHOR
      if (methflag.eq.2) then
        call getrho(cre,cim,rhor,rhoi,nsurft)
        tmp1 = 0.d0
        do i=1,nsurft
          if (rhor(i,i).gt.tmp1) then
            tmp1 = rhor(i,i)
            nsurf = i
          endif
        enddo
      endif

c SCDM and CSDM, check for a decoherence state switch
      if (methflag.eq.3.or.methflag.eq.4) then
        lhop = .false.
        call decocheck(phop,nsurf,lhop)
        if (lhop) then
        hstep=hstep0/5.d0
        call getgrad(xx,pp,nsurf,pe,cre,cim,gv,gcre,gcim,nat,
     &    phop,dmag,dvec,pem,gpem,phase)
        endif
      endif

c for steepest descent (TFLAG(1) = 1 or 4)
c set PP to zero at every step
      if (tflag(1).eq.1.or.tflag(1).eq.4) then
        do i=1,nat
        do j=1,3
        pp(j,i) = 0.d0
        enddo
        enddo
      endif

c on-the-fly analyses
      ce = -pe/dble(nat)
      call gettemp(pp,mm,nat,temp,ke)
      call angmom(xx,pp,mm,nat,bigj,bigjtot)
      te = ke + pe
      call lindemann(xx,nat,step,time,arij,arij2,lind)
      if ( (mod(istep,nprint).eq.0.or.istep.eq.1) .and. lwrite(43) ) 
     &                           call rarray(ithistraj,time,xx,nat)
      if (lwrite(41)) then
        if (mod(istep,nprint).eq.0.or.istep.eq.1) then
          call radialdist(ithistraj,xx,nat,step,time,nbinrad,raddist,1)
        else
          call radialdist(ithistraj,xx,nat,step,time,nbinrad,raddist,0)
        endif
      endif
      atemp = atemp + temp*step
      atemp2 = atemp2 + (temp**2)*step
      stemp = dsqrt(atemp2/time - (atemp/time)**2)

      if (rhor(1,1).ge.0.2d0.and.rhor(1,1).le.0.8d0.and.lwrite(35)) then
      write(35,1035)ithistraj,istep,time*autofs,step*autofs,
     &     taup(1)*autofs,taup(2)*autofs,
     &     taup(3)*autofs,rhor(1,1),rhor(2,2)
      endif
 1035  format(2i8,10f15.5)

c write info every nprint steps
      call getrho(cre,cim,rhor,rhoi,nsurft)
      if (mod(istep,nprint).eq.0.or.istep.eq.1) then
c molden format
       nistep = nistep + 1
       if (lwrite(81)) then
          write(81,*)nat
          write(81,1081)istep,pe
          do i=1,nat
            write(81,1080)symbol(i),(pp(j,i),j=1,3)
          enddo
       endif
       if (lwrite(80)) then
          write(80,*)nat
          write(80,1081)istep,pe
          do i=1,nat
            write(80,1080)symbol(i),(xx(j,i)*autoang,j=1,3)
 1080       format(a2,3f15.5)
 1081       format(i7,50f15.5)
          enddo
       endif
       if (lwrite(82)) write(82,REC=nistep)nat,istep,pe,
     &  (symbol(i),i=1,nat),((xx(j,i)*autoang,j=1,3),i=1,nat)
       if (lwrite(83)) write(83,REC=nistep)nat,istep,pe,
     &  (symbol(i),i=1,nat),((pp(j,i),j=1,3),i=1,nat)
       if (lwell.gt.0) then
         call detwell(pe,xx,welli)
         if (welli.eq.lwell) then
           istepw=istepw+1
           nistepw=nistepw+1
           if (lwrite(86)) then
             write(86,*)nat
             write(86,1081)istepw,pe
             do i=1,nat
               write(86,1080)symbol(i),(xx(j,i)*autoang,j=1,3)
             enddo
           endif
           if (lwrite(87)) then
             write(87,*)nat
             write(87,1081)istepw,pe
             do i=1,nat
               write(87,1080)symbol(i),(pp(j,i),j=1,3)
             enddo
           endif
           if (lwrite(88)) write(88,REC=nistepw)nat,istepw,pe,
     &      (symbol(i),i=1,nat),((xx(j,i)*autoang,j=1,3),i=1,nat)
           if (lwrite(89)) write(89,REC=nistepw)nat,istepw,pe,
     &      (symbol(i),i=1,nat),((pp(j,i),j=1,3),i=1,nat)
         else if (welli.eq.3) then 
           print *,"welli=",welli," not supported"
         endif
       endif
       if (lwrite(42)) call honey(ithistraj,xx,nat,time)
       if (termflag.eq.2) then
          if(nsurft.eq.1) then
            write(6,101)istep,time*autofs,pe*autoev,
     &    dsqrt(gvmag)*autoev/autoang
          else
            write(6,101)istep,time*autofs,(pem(j)*autoev,j=1,nsurft),
     &    dsqrt(gvmag)*autoev/autoang
          endif
       else
        if (nsurft.eq.1.and.lwell.eq.0) then
c          write(6,101)istep,time*autofs,temp,
c     &    ke*autoev,pe*autoev,te*autoev,
c     &    lind,atemp/time,stemp
          write(6,103)istep,time*autofs,temp,
     &    ke*autoev,pe*autoev,te*autoev,
     &    lind,atemp/time,stemp
c     &    ,wellindex1,wellindex2,xwell
        else if (nsurft.eq.1.and.lwell.gt.0) then
          write(6,104)istep,istepw,time*autofs,temp,
     &    ke*autoev,pe*autoev,te*autoev,
     &    lind,atemp/time,stemp
        else
c          write(6,101)istep,time*autofs,temp,
c     &    ke*autoev,pe*autoev,te*autoev,
c     &    (rhor(k,k),k=1,nsurft),pem(1)*autoev,pem(2)*autoev
c temp
          write(6,102)istep,time*autofs,temp,
     &    ke*autoev,pe*autoev,te*autoev,
     &    (rhor(k,k),k=1,nsurft),(pem(k)*autoev,k=1,nsurft)
     &     ,dmag,h12x*autocmi
c     &     dmag,h12x*autocmi,h12y*autocmi
c     &    ,wellindex1,wellindex2,xwell
        endif
       endif
      endif
 101  format(i15,f15.3,20f21.9)
 102  format(i15,f15.3,100f15.5)
 103  format(i15,f15.3,7f15.5,2i5,f8.2)
 104  format(i15,i15,f15.3,7f15.5,2i5,f8.2)
 105  format(i15,i15,f15.3,20f15.5)

 10   continue

c ######################################################################
c ######################################################################

 999  continue

      if (lwrite(16)) 
     & write(16,'(2f13.5,2e18.5,3e18.5)')time*autofs,1.d4,
     & phop(1)*step,rhor(1,1),
     & (taup2(i)*autofs,i=1,3)

c trajectory has ended

c write final coordinates and momenta to output
      IF (outcome.ge.0) THEN
      if (termflag.eq.2) then
          write(6,101)istep,time*autofs,
     &    pe*autoev,dsqrt(gvmag)*autoev/autoang
      else
        if (nsurft.eq.1.and.lwell.eq.0) then
          write(6,101)istep,time*autofs,temp,
     &    ke*autoev,pe*autoev,te*autoev,
     &    lind,atemp/time,stemp
        else if (nsurft.eq.1.and.lwell.gt.0) then
          write(6,105)istep,istepw,time*autofs,temp,
     &    ke*autoev,pe*autoev,te*autoev,
     &    lind,atemp/time,stemp
        else
          write(6,101)istep,time*autofs,temp,
     &     ke*autoev,pe*autoev,te*autoev,
     &     (rhor(k,k),k=1,nsurft)
        endif
      endif
      if (lwrite(41)) then
      call radialdist(ithistraj,xx,nat,step,time,nbinrad,raddist,1)
      endif
      if (lwrite(42)) call honey(ithistraj,xx,nat,time)
      ENDIF

      if (lwrite(20)) then
      write(20,*)ithistraj,pe*autoev
      do i=1,nat
      write(20,1010)symbol(i),mm(i)/amutoau,(xx(j,i)*autoang,j=1,3)
      enddo
      endif

      if (lwrite(21)) then
      write(21,*)ithistraj,ke*autoev
      do i=1,nat
      write(21,1010)symbol(i),mm(i)/amutoau,(pp(j,i),j=1,3)
      enddo
      endif

      if (lwrite(91)) write(91,*)ithistraj,rvdw*autoang,evdw*autocmi

 1001 return
      end

