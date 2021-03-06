      subroutine mc(im)

c Monte Carlo evaluations of classical N(E) and rho(E) or Q(T)

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer i,nbound,im,ii,j,nclu,ido,isamp,hit,hitho,neval,ibin,
     & k,nho,nfr,samptype,icalc,a1,a2,n1,nn,itype(3*mnat),negrid
      double precision xvec(mnat*3),mmm(mnat),xxm(3,mnat),ppm(3,mnat),
     & rin(3*mnat),rout(3*mnat),vec(3,mnat),etot,rtmp1,rtmp2,
     & disp(3*mnat),vec2(3,mnat,3*mnat),esamp,esampho,ne,neho,
     & nemc,nemcho,pp,pref,hscal,neana,tmp,ceval,nsav(mntraj),perr,
     & rhoana,etmp,nehom,nehop,wsav(mntraj),wwtot,eetmp(mnsurf),
     & rhoho,estep,nep,nem,rhomc,hrad,rrr,zpe,qpg,qpgtilde,
     & ptest,tmp1,tmp2,cevalmax,irot(mnat*3),x10,x21,x00,sinth,dd,
     & mom(2),rad,radmax,ww,rad0,alpha,kt,qscal,qana,qanaq,
     & qvib,qho,qturn,phi,em,ep,ftor,ee,rholz,rholzi,eex,eey,
     & plz,h12,gapgrad,neim,neip,rhoi,est
      character*2 symb(mnat)

      common/lz/gapgrad,h12

c ********************
c PARAMETERS
c      icalc=0 ! N(E) and r(E)
      icalc=1 ! Q(T)
c NUMERICAL PARAMETERS
      ceval=3.d0 ! skip evaluations at HO estimated energies > CEVAL*ETOT
c      samptype=0 ! dumb sampling, uniform in a hypercube
      samptype=1 ! uniform in d (size of hypercube) WORKS MUCH BETTER THAN ^
c      samptype=2 ! gaussian in d DOESN'T SEEM TO HELP AT ALL
c      samptype=3 ! uniform in r (generalized hyperradius) DOESN'T WORK YET!
      rad0=.0d0 ! center of gaussian
      alpha=5.d0 ! 1 sigma width of gaussian
      estep=0.5d0/autocmi ! numerically differentiate N(E) with stepsize ESTEP
      qturn=3.d0 ! use turning points for Q(T) calc for E=QTURN*kT
c ********************

c Convert to local variables for convenience
      do i=1,natom(im)
       ii=i+iatom(im)
       mmm(i) = mm(ii)
       symb(i) = symbol(ii)
       do j=1,3
        xxm(j,i) = xx0(j,ii)
       enddo
      enddo
      nclu=natom(im)
      etot=scale0im(im)
      kt=scale0im(im)*kb*autoev  ! scale0im is converted to eV in READIN...
      if (icalc.eq.1) etot=kt*qturn


c NBOUND is the dimension of the sampled space
      if (nmtype(im).eq.0) then
        nbound=3*nclu-6
      else
        nbound=3*nclu-7
      endif
      if (nclu.eq.2) then
        nbound=1
      endif

      do i=1,nbound
        itype(i)=0
        do j=1,mccurv
          if (mcmode(j).eq.i) itype(i)=mctype(j)
        enddo
      enddo

c Precompute things for the curvilinear modes
      do i=1,mccurv
        if (mctype(i).eq.1) then
c         calculate Brot
          ii=0
          nn=mcmode(i)
          mom(1)=0.d0
          mom(2)=0.d0
          do j=1,nclu
            if (j.eq.mcpar(1,nn)) then
            ii=1
            elseif (j.eq.mcpar(2,nn)) then
            ii=2
            else
            x21=0.d0
            x10=0.d0
            x00=0.d0
            do k=1,3
            x21=x21+(xxm(k,mcpar(2,nn))-xxm(k,mcpar(1,nn)))**2
            x10=x10+(xxm(k,mcpar(1,nn))-xxm(k,j))**2
            x00=x00+(xxm(k,mcpar(1,nn))-xxm(k,j))*
     &            (xxm(k,mcpar(2,nn))-xxm(k,mcpar(1,nn)))
            enddo
            x21=dsqrt(x21)
            x10=dsqrt(x10)
            sinth=dsqrt(1.d0-(x00/(x21*x10))**2)
            dd=x10*sinth  ! the distance of atom j from the vector defining the torsion
            if (ii.ne.0) mom(ii)=mom(ii)+dd*mmm(j)
            endif
          enddo
          irot(mcmode(i))=mom(1)*mom(2)/(mom(1)+mom(2))
          write(6,1066)"For torsion ",mcmode(i),"   B1 = ",
     &        autocmi/(2.d0*mom(1))," cm-1"
          write(6,1067)"                    B2 = ",
     &        autocmi/(2.d0*mom(2))
          write(6,1067)"                  Beff = ",
     &        autocmi/(2.d0*irot(mcmode(i)))
          write(6,*)
c numerical freq of torsion
          do j=1,nbound
            disp(j)=0.d0
          enddo
          call nmpot(symb,xx0,mm,nclu,nbound,
     &             vec2,repflag,nsurft,disp,itype,mcpar,eetmp)
          ee=eetmp(nsurf0)
          disp(nn)=0.01d0
          call nmpot(symb,xx0,mm,nclu,nbound,
     &             vec2,repflag,nsurft,disp,itype,mcpar,eetmp)
          ep=eetmp(nsurf0)
          disp(nn)=-0.01d0
          call nmpot(symb,xx0,mm,nclu,nbound,
     &             vec2,repflag,nsurft,disp,itype,mcpar,eetmp)
          em=eetmp(nsurf0)
          ftor=(ep-2.d0*ee+em)/0.01d0**2   ! 2nd deriv
          ftor=dsqrt(ftor/irot(mcmode(i))) ! freq
        print *,"q* (",ftor*autocmi,") = ",
     & dexp(-0.5d0*ftor/kt)/(1.d0-dexp(-ftor/kt))
        print *,"q (",freq(nn,im)*autocmi,") = ",
     & dexp(-0.5d0*freq(nn,im)/kt)/(1.d0-dexp(-ftor/kt))
        print *,"q*/q = ",dexp(-0.5d0*ftor/kt)/(1.d0-dexp(-ftor/kt))/
     & (dexp(-0.5d0*freq(nn,im)/kt)/(1.d0-dexp(-ftor/kt)))

 1066   format(a,i5,a,f10.4,a)
 1067   format(a,f10.4)
        endif
      enddo
 
c find turning points on the real PES for each normal mode
c calculate prefactor HSCAL
      hscal=1.d0
      qscal=1.d0
      rrr=1.d0
      if (lreadhess) write(6,*)"Reading turning points from unit 70"
      write(6,*)"Eturn= ",etot*autoev," eV"
      write(6,'(2a)')" Mode    Freq(cm-1)    Rin(m.s.au)   Rout(m.s.au)"
     &         ,"   RturnHO(m.s.au)    VHO(Rout)/E     VHO(Rout)/E"
      do ido=1,nbound
        do i=1,nclu
        do j=1,3
          vec(j,i)=nmvec(j,i,ido,im) ! mass-scaled normal mode vectors
          vec2(j,i,ido)=nmvec(j,i,ido,im)
        enddo
        enddo

        if (itype(ido).eq.0) then
c       Cartesian displacement

        if (lreadhess) then ! read from fort.70
c        if (.false.) then
          read(70,*)ii,rin(ido),rout(ido)
        else
c         calculate turning points
          call nmturn(symb,xxm,mmm,nclu,vec,etot,rtmp1,rtmp2,
     &                repflag,nsurf0,nsurft) ! returns turning points
          rin(ido)=rtmp1
          rout(ido)=rtmp2
          write(70,*)ido,rin(ido),rout(ido) ! save restart info to fort.70
        endif

        rtmp1 = dsqrt(2.d0*etot/(freq(ido,im)**2*mu)) ! harmonic turning point
        tmp1=0.5d0*mu*freq(ido,im)**2*rin(ido)**2/etot ! harmonic PES at real turning point
        tmp2=0.5d0*mu*freq(ido,im)**2*rout(ido)**2/etot
        write(6,1010)ido,freq(ido,im)*autocmi,rin(ido),rout(ido),rtmp1,
     &       tmp1,tmp2
        hscal=hscal*dabs(rout(ido)-rin(ido))/(2.d0*pi) ! hypervolume of box/h**nbound
        qscal=qscal*dabs(rout(ido)-rin(ido))*dsqrt(mu*kt/(2.d0*pi))

        elseif (itype(ido).eq.1) then
c       torsion

        rin(ido)=-pi*dsqrt(irot(ido)/mu) ! mass scale the torsion
        rout(ido)=pi*dsqrt(irot(ido)/mu) ! mass scale the torsion
        rtmp1=0.d0
        tmp1=0.d0
        tmp2=0.d0
        write(6,1010)ido,freq(ido,im)*autocmi,rin(ido),rout(ido),rtmp1,
     &       tmp1,tmp2," TORSION"
        hscal=hscal*dabs(rout(ido)-rin(ido))/(2.d0*pi) ! hypervolume of box/h**nbound
        qscal=qscal*dabs(rout(ido)-rin(ido))*dsqrt(mu*kt/(2.d0*pi))

        elseif (itype(ido).eq.-1) then
c skip
        rtmp1=0.d0
        tmp1=0.d0
        tmp2=0.d0
        write(6,1010)ido,freq(ido,im)*autocmi,0.d0,0.d0,rtmp1,
     &       tmp1,tmp2," SKIPPING"
        endif

      enddo
      write(6,*)

c analytic HO solution
c NOT ACCURATE FOR N(E) and r(E) for NFR != 0
      neana=1.d0
      qana=1.d0
      qanaq=1.d0
      nho=0
      nfr=0
      zpe=0.d0
      do i=1,nbound
        zpe=zpe+freq(i,im)/2.d0
        if (itype(i).eq.0) then
          nho=nho+1
          neana=neana*(etot/freq(i,im))
          qana=qana*(kt/freq(i,im))
          qanaq=qanaq*dexp(-0.5d0*freq(i,im)/kt)/
     &                (1.d0-dexp(-freq(i,im)/kt))
        elseif (itype(i).eq.1) then
          neana=neana*dsqrt(etot*2.d0*irot(i))   ! root(E/B) = root(E 2 I)
          nfr=nfr+1
          qana=qana*dsqrt((kt*irot(i)*(2.d0*pi)))
          qanaq=qanaq*dsqrt((kt*irot(i)*(2.d0*pi)))
        endif
      enddo
      if (nfr.eq.0) then
c     only HOs
        do i=1,nho
          neana=neana/dble(i)
        enddo
        rhoana=neana/etot*dble(nho)
      elseif (nfr.eq.1) then
c     one torsion
        neana=neana*2.d0**(nho+1)
        ii=2*nho+1
        do i=1,ii,2
        neana=neana/dble(i)
        enddo
        rhoana=neana/etot*dble(2*nho+1)/2.d0
      endif

c normalize gaussian weights 
      wwtot=0.d0
      do i=1,1000
        rad=dble(i)/1000.d0
        wwtot=wwtot+dexp(-(rad-rad0)**2/(2.d0*alpha**2))/1000.d0
      enddo
c      do i=1,100
c        rad=dble(i)/100.d0
c        print *,rad,dexp(-(rad-rad0)**2/(2.d0*alpha**2))/wwtot
c      enddo

c MC SAMPLING LOOP
      cevalmax=1.d0
      tmp=0.d0
      ne=0.d0
      neho=0.d0
      nehop=0.d0
      nehom=0.d0
      hit=0
      hitho=0
      isamp=0
      neval=0
      qvib=0.d0
      qho=0.d0
      rholz=0.d0

      if (icalc.eq.0) then
      write(6,'(4a)')"      Nhit     Nsamp     Neval   %hit   ",
     & "    W        %err     r,1/cm-1    W/W(HO,an) ",
     & " W(HO)/W(HO,an) r/r(HO,an) r(HO)/r(HO,an)",
     & "    Cmax "
      elseif (icalc.eq.1) then
      write(6,'(6a)')"      Nhit     Nsamp     Neval   %hit   ",
     & "    Q        %err     Q(HO)       Q(HO,an) ",
     & "     Q(HO,q)      Q/Q(HO,an) Q(HO)/Q(HO,an)",
     & "      QPG        Q~PG   "
      endif

      DO WHILE(hit.lt.ntraj)

        isamp=isamp+1 ! total number of trials

        rad=1.d0
        if (samptype.eq.1) call ranno(rad,0.d0,1.d0)
        if (samptype.eq.2) call rangaussian(rad,rad0,alpha)
        if (samptype.eq.3) call ranno(rad,0.d0,1.d0)
c        rad=0.2

c       get geometry (expressed here as a vector of displacements)
  10    radmax=0.d0
        hrad=0.d0
        do i=1,nbound
        rtmp1=rin(i)
        rtmp2=rout(i)
        if (samptype.eq.3) rtmp1=-1.d0
        if (samptype.eq.3) rtmp2=1.d0
        if (itype(i).eq.0) then
        call ranno(disp(i),rtmp1,rtmp2) ! random, uniform displacements
        hrad=hrad+disp(i)**2
        radmax=max(radmax,disp(i)/rtmp1)
        radmax=max(radmax,disp(i)/rtmp2)
c        print *,i,rin(i),rout(i),radmax
        else
        disp(i)=0.d0
        endif
        enddo
        hrad=dsqrt(hrad)
        if (samptype.eq.3) then
        if (hrad.gt.1.d0) go to 10
        do i=1,nbound
        disp(i)=disp(i)/hrad
        enddo
        endif

c SCAN TMP
c        if (isamp.eq.52) stop
c        do i=1,nbound
c        if (itype(i).eq.0) then
c        rad=dble(isamp-1)/50.d0
c        disp(i)=rin(i)+(rout(i)-rin(i))*rad
c        endif
c        enddo

        esampho=0.d0
        do i=1,nbound
        if (itype(i).eq.0) then
        if (samptype.eq.1.or.samptype.eq.2) disp(i)=disp(i)*rad/radmax
        if (samptype.eq.3) disp(i)=disp(i)*max(-rin(i),rout(i))*rad
        esampho=esampho+0.5d0*mu*freq(i,im)**2*disp(i)**2
        endif
        enddo

        ww=1.d0
        if (samptype.eq.2) 
     &    ww=dexp(-(rad-rad0)**2/(2.d0*alpha**2))/wwtot
c     &    ww=dexp(-(rad-rad0)**2/(2.d0*alpha**2))/sqrt(2.d0*pi*alpha**2)
        if (samptype.ne.0.and.samptype.ne.3.and.nho.ne.0) 
     &   ww=dble(nho)*rad**(nho-1)/ww
c     &   ww=dble(nbound)*rad**(nbound-1)/ww
        if (samptype.eq.3) 
     &   ww=dble(nbound)*rad**(nbound-1)/ww

c sample torsion
        do i=1,nbound
        if (itype(i).eq.1) then
        rtmp1=rin(i)
        rtmp2=rout(i)
        call ranno(disp(i),rtmp1,rtmp2) ! random, uniform displacements
        endif
        enddo



c N(E) and r(E)
        IF (icalc.eq.0) THEN
c       evaluate potential (or skip if HO estimate is too high)
        if (esampho.lt.ceval*etot.or.ceval.lt.0.d0) then
        neval=neval+1
c       make displacements & return energy
        call nmpot(symb,xx0,mm,nclu,nbound,
     &             vec2,repflag,nsurft,disp,itype,mcpar,eetmp)
        esamp=eetmp(nsurf0)
c        esamp=(eetmp(1)+eetmp(2))/2.d0
        write(71,'(20f15.5)')rad,(disp(k),k=1,nbound),
     7             eetmp(1)*autoev,eetmp(2)*autoev
c        print *,isamp,esamp*autoev,esampho*autoev
        else
        esamp=esampho
        endif

c       evaluate MC sampled harmonic N(E)
        if (esampho.lt.etot) then
          hitho=hitho+1
          pp=dsqrt(2.d0*mu*(etot-esampho))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          neho=neho+pref*ww
          nemcho=neho/dble(isamp)*hscal
        endif

c       numerically differentiate harmonic N(E) 
        etmp=etot-estep
        if (esampho.lt.etmp) then
          pp=dsqrt(2.d0*mu*(etmp-esampho))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          nehom=nehom+pref*hscal*ww
c        endif

        etmp=etot+estep
c        if (esampho.lt.etmp) then
          pp=dsqrt(2.d0*mu*(etmp-esampho))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          nehop=nehop+pref*hscal*ww
        rhoho=(nehop-nehom)/(2.d0*estep*dble(isamp))
        endif

c       numerically differentiate true N(E) 
        etmp=etot-estep
        if (esamp.lt.etmp) then
          pp=dsqrt(2.d0*mu*(etmp-esamp))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          nem=nem+pref*hscal*ww
c        endif

        etmp=etot+estep
c        if (esamp.lt.etmp) then
          pp=dsqrt(2.d0*mu*(etmp-esamp))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          nep=nep+pref*hscal*ww
          rhomc=(nep-nem)/(2.d0*estep*dble(isamp))
        endif

c       evaluate MC sampled true N(E)
        if (esamp.lt.etot) then

c LZ EVAL
          negrid=1000
c          h12=75.d0/autocmi
c          gapgrad=0.0529193396    ! gmp2dz mag at MSX
c          print *,gapgrad,h12*autocmi
          rholzi=0.d0
          est=(etot-esamp)/dble(negrid)
          do i=0,negrid-1
             eex=(0.5d0+dble(i))*est
             eey=etot-esamp-eex

c             plz=1.d0-dexp(-2.d0*pi*h12**2/gapgrad
c     &              *dsqrt(0.5d0*mu/eex))  ! LJ prob
             plz=1.d0

             etmp=eey+est/2.d0
             pp=dsqrt(2.d0*mu*(etmp))
             call nsphere(nho,pp,pref)
             neip=pref*hscal*ww

             etmp=eey-est/2.d0
             if (etmp.lt.0.d0) etmp=0.d0
             pp=dsqrt(2.d0*mu*(etmp))
             call nsphere(nho,pp,pref)
             neim=pref*hscal*ww

             rhoi=(neip-neim)/est

             rholzi=rholzi+plz*rhoi*est
c             print *,eex,eey,plz,rhoi
          enddo
          rholz=rholz+rholzi

          if (esampho/etot.gt.cevalmax) cevalmax=esampho/etot
          hit=hit+1
          pp=dsqrt(2.d0*mu*(etot-esamp))
c          call nsphere(nbound,pp,pref)
          call nsphere(nho,pp,pref)
          nsav(hit)=pref*hscal ! save for bootstrap analysis
          wsav(hit)=ww
          ne=ne+nsav(hit)*wsav(hit) ! treat sampling over momentum analytically
c          write(99,*)rad,hit,nsav(hit),wsav(hit),nsav(hit)*wsav(hit)
          nemc=ne/dble(isamp)
          perr=0.d0
          if (hit.gt.nprint*10) nprint=nprint*10
          if (mod(hit,nprint).eq.0) then
            call bootstrap(nsav,wsav,mntraj,hit,isamp,perr) ! compute boostrap uncertainty
            write(6,1060)hit,isamp,neval,
     &            dble(hit)/dble(neval)*100.d0,nemc,perr,rhomc/autocmi,
     &            nemc/neana,nemcho/neana,
     &            rhomc/rhoana,rhoho/rhoana,cevalmax,nemc/nemcho,
     &            rholz/dble(isamp),rholz/(dble(isamp)*nemc)
          endif
c          if ((esampho/esamp).gt.tmp) print *,tmp
c          if ((esampho/esamp).gt.tmp) tmp=esampho/esamp
c          write(99,*)isamp,rad,hit,nsav(hit),ww,nsav(hit)*wsav(hit)
        else
c          write(99,*)isamp,rad,0,0.,wsav(hit),0.
        endif

      ELSEIF (icalc.eq.1) THEN ! Q(T)
          qho=qho+dexp(-esampho/kt)*ww*qscal
          neval=neval+1
          hit=hit+1
          call nmpot(symb,xx0,mm,nclu,nbound,
     &             vec2,repflag,nsurft,disp,itype,mcpar,eetmp)
c          esamp=(eetmp(1)+eetmp(2))/2.d0
          esamp=eetmp(nsurf0)
          qvib=qvib+dexp(-esamp/kt)*ww*qscal
          nsav(hit)=dexp(-esamp/kt)*qscal ! save for bootstrap analysis
          wsav(hit)=ww
          if (hit.gt.nprint*10) nprint=nprint*10
c          if (.true.) print *,esamp*autoev,
c     & dexp(-esamp/kt),ww,qscal,rad,dexp(-esamp/kt)*ww*qscal
c          if ((qvib/(dble(isamp)*qana)).gt.10.d0) print *,esamp*autoev,
c     & dexp(-esamp/kt),ww,qscal,qvib/dble(isamp),qana
          if ((qvib/(dble(isamp)*qana)).gt.10.d0) stop
          if (mod(hit,nprint).eq.0) then
            call bootstrap(nsav,wsav,mntraj,hit,isamp,perr) ! compute boostrap uncertainty
          qpg=qvib/dble(isamp)*qanaq/qana
          qpgtilde=qpg/dexp(-zpe/kt)
            write(6,1070)hit,isamp,neval,
     &            dble(hit)/dble(neval)*100.d0,qvib/dble(isamp),perr,
     &            qho/dble(isamp),qana,qanaq,qvib/(dble(isamp)*qana),
     &            qho/(dble(isamp)*qana),qpg,qpgtilde
          endif
      ENDIF

      ENDDO
c     END SAMPLING LOOP

      stop
 1010 format(i5,6f15.5,1a)
 1011 format(i10,10f15.5)
 1060 format(3i10,0pf7.1,1pe13.5,0pf7.1,0p20f13.5)
 1070 format(3i10,0pf7.1,1pe13.5,0pf7.1,1p3e13.5,0p2f13.5,1p2e13.5)
      return
      end
