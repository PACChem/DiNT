      program moments
      implicit double precision(a-h,o-z)

      parameter(mdat=100000)
      parameter(mobs=40)
      character*5 m1,m2
      integer bsample
      integer traj,state,outcome,idat(mobs)
      double precision li,j1i,j2i,lf,j1f,j2f,kbev,kbcmi,
     & maxediffi,maxedifff,dat(mobs),maxd,mind,sd(mdat)

      dimension ejmom(-1:1,0:3,-1:1,0:3,0:3)

      read (5,*)ntraj,bsample,bmax,m1,m2,ttt
      if (ntraj.gt.mdat) then
         print *,"Increase MDAT to ",ntraj
         stop
      endif

c orders
       ndblej=2  ! dble J'
       npwre=3   ! max order of E
       npwrj=3   ! max order of J

c constants 
      evtocmi = 8065.7208
      evtoau = 1./27.211
      kbev = 8.61692d-05
      kbcmi = 0.695016d0
      pi=dacos(-1.d0)
      dsmall=2.d0 ! 2 cm-1
      djsmall=.02d0 ! 0.02 hbar

c output
      open(7,file='dint.mom')

c print
      print *,m1," + ",m2," collisions "
      print *,ntraj," trajectories"
      if(bsample.eq.1) print *,"b was sampled evenly (0,",bmax,")"
      if(bsample.eq.2) print *,"b**2 was sampled evenly (0,",bmax,"**2)"
      print *

c compute Z/Z
      print *,"   T,K    ZHS/ZLJ     ZHS cm3/s    ZLJ cm3/s"
      call zonz(ttt,bmax,m1,m2,zz,zhs,zlj)
      print *," ",ttt," ",zz,zhs,zlj

c initialize
       do i=-1,1
       do j=0,npwre
       do k=-1,1
       do l=0,npwrj
       do m=0,ndblej
       ejmom(i,j,k,l,m)=0.d0
       enddo
       enddo
       enddo
       enddo
       enddo
       etotdmax = -1.d0
       tt=-1.d0
       ttmin=1d10
       ttmax=-1d10
       ttav=0.d0

       etriav=0.d0
       etotiav=0.d0
       etotdmax=-1.d0
       etotdmue=0.d0
       etotdrms=0.d0
       eeiavg=0.d0
       xjavg=0.d0

       maxd=0.d0
       mind=0.d0
       maxediffi=0.d0
       maxedifff=0.d0

       bav=0.d0

       nt=0
       md=0
       wd=0
       wu=0
       wt=0
       mm=0
       mu=0
       mjd=0
       mju=0
       mdd=0
       mdu=0
       mud=0
       muu=0

       dnd=0
       dnw=0

       nsmall=0
       njsmall=0
       nsmallr=0
       nsmallj=0

       tdel=0.d0
       ndel=0

c loop over trajectories
       nobs=32
       if (nbos.gt.mobs) then
         print *,"Increase NOBS to ",nobs
         stop
       endif
       do it=1,ntraj
       read (5,*)(idat(j),j=1,3)

        traj = idat(1)
        state = idat(2)
        outcome = idat(3)

       if (outcome.ne.0) then
       backspace(5)        
       read (5,*)(idat(j),j=1,3),(dat(j),j=4,nobs)

        time = dat(4)
        temp = dat(5)
        temps = dat(6)
        bi = dat(7)
        li = dat(8)
        eorbi = dat(9)
        ereli = dat(10)
        vi = dat(11)
        j1i = dat(12)
        erot1i = dat(13)
        evib1i = dat(14)
        v1i = dat(15)
        j2i = dat(16)
        erot2i = dat(17)
        evib2i = dat(18)
        v2i = dat(19)
        vf = dat(20)
        lf = dat(21)
        eorbf = dat(22)
        erelf = dat(23)
        j1f = dat(24)
        erot1f = dat(25)
        evib1f = dat(26)
        v1f = dat(27)
        j2f = dat(28)
        erot2f = dat(29)
        evib2f = dat(30)
        v2f = dat(31)
        tt = dat(32) 
        k1i=0.
       endif

! Check energy error due to finite fragment separation
        ediffi = vi-v1i-v2i
        edifff = vf-v1f-v2f
        if (dabs(edifff)*evtocmi.gt.1.d3) then
c           print *,"Vf = ",edifff*evtocmi," (traj ",it,")"
           outcome = 3 
c           write(11,*)bi,traj
        endif

! Total energy
        eei = vi+erot1i+evib1i+erot2i+evib2i
        etoti = eorbi+ereli+vi+erot1i+evib1i+erot2i+evib2i
        etotf = eorbf+erelf+vf+erot1f+evib1f+erot2f+evib2f
        etotd = etotf-etoti
        if (dabs(etotd)*evtocmi.gt.1000.) then
           print *,"Delta E = ",etotd*evtocmi," (traj ",it,")"
           outcome = 3
        endif

! check for bad trajectories
        iskip=0
        if (outcome.eq.0.or.outcome.ge.2) then
c                print *,"Outcome = ",outcome," (traj ",it,")"
                iskip=1
        endif

        IF (iskip.eq.0) THEN
        nt=nt+1

! sampling temperature; compute max, min, avg
        if (tt.gt.ttmax) ttmax = tt
        if (tt.lt.ttmin) ttmin = tt
        ttav = ttav+tt

! Check energy conservation, avg, mue, rms, max
        etotiav = etotiav + etoti
        etotdmue = etotdmue + dabs(etotd)
        etotdrms = etotdrms + (etotd)**2
        if (dabs(etotd).gt.etotdmax) itetotdmax = it
        if (dabs(etotd).gt.etotdmax) etotdmax = dabs(etotd)

! Check energy error due to finite fragment separation
        ediffi = vi-v1i-v2i
        edifff = vf-v1f-v2f
        if (dabs(ediffi).gt.maxediffi) maxediffitraj = it
        if (dabs(edifff).gt.maxedifff) maxediffftraj = it
        maxediffi=max(dabs(ediffi),maxediffi)
        maxedifff=max(dabs(edifff),maxedifff)
        avgediffi=avgediffi+dabs(ediffi)
        avgedifff=avgedifff+dabs(edifff)

! Initial translational energy (orbital + relative) & average
        etri = eorbi+ereli
        etriav = etriav + etri

! Impact parameter
        bav = bav + bi                                  ! compute average impact parameter b

! Delta J
        delj = j1f - j1i
        if (dabs(delj).lt.djsmall) then
           delj=djsmall/2.d0
           njsmall=njsmall+1
           if (mod(njsmall,2).eq.0) delj=-djsmall/2.d0
        endif

! Delta E
        delr = erelf+eorbf-ereli-eorbi                 ! change in trans energy
        delf1 = erot1f+evib1f+v1f-(erot1i+evib1i+v1i)      ! change in internal energy
        delf2 = erot2f+evib2f+v2f-(erot2i+evib2i+v2i)      ! change in internal energy
        if (v1i.eq.0d0.and.v2i.eq.0d0) then 
!       if LDOFRAG=.false.
        d = -delr                                         ! delta E 
        else
!       if LDOFRAG=.true., e.g., for a molecular collider
        d = delf1                                        ! delta E 
        endif
        d=d*evtocmi
        if (dabs(d).lt.dsmall) then
           d=dsmall/2.d0
           nsmall=nsmall+1
           if (mod(nsmall,2).eq.0) d=-dsmall/2.d0
        endif
        if (d.gt.maxd) then
                 maxd=d
                 maxdit=it
        endif
        if (d.lt.mind) then
                 mind=d
                 mindit=it
        endif
        eeiavg=eeiavg+eei
        xjavg=xjavg+j1i

c        print *,bi,j1i,d,delj

!       Weights
        w = 1.d0
        if (bsample.eq.1) w = w*2.d0*bi/bmax    ! w=impact parameter b scaling

        if (d.lt.0.d0) then
          md=md+1
          sd(md)=d*w*zz
        endif

        do kj=0,ndblej
        xkj=dble(j1i**kj)
        do ie=0,npwre
        do ij=0,npwrj
        ejmom(0,ie,0,ij,kj)=ejmom(0,ie,0,ij,kj)+d**ie*delj**ij*w*xkj
        iesign=nint(d/dabs(d))
        ijsign=nint(delj/dabs(delj))
        ejmom(iesign,ie,0,ij,kj)=ejmom(iesign,ie,0,ij,kj)
     &     +d**ie*delj**ij*w*xkj
        ejmom(0,ie,ijsign,ij,kj)=ejmom(0,ie,ijsign,ij,kj)
     &     +d**ie*delj**ij*w*xkj
        ejmom(iesign,ie,ijsign,ij,kj)=ejmom(iesign,ie,ijsign,ij,kj)
     &     +d**ie*delj**ij*w*xkj
c        if (d.lt.0.d0.and.kj.eq.0) then
c        ejmom(-1,ie,0,ij,kj)=ejmom(-1,ie,0,ij,kj)+d**ie*delj**ij*w*xkj
c        else
c        ejmom(1,ie,0,ij,kj)=ejmom(1,ie,0,ij,kj)+d**ie*delj**ij*w*xkj
c        endif
c        if (delj.lt.0.d0) then
c        ejmom(0,ie,-1,ij,kj)=ejmom(0,ie,-1,ij,kj)+d**ie*delj**ij*w*xkj
c        else
c        ejmom(0,ie,1,ij,kj)=ejmom(0,ie,1,ij,kj)+d**ie*delj**ij*w*xkj
c        endif
        enddo
        enddo
        enddo

       ENDIF
      enddo

      print *
      print *,"N traj = ",nt
      print *,"N E down = ",md," (",
     &    (dble(md)/dble(nt))*100.," %)",md,mu,md+mu,wd,wu,wd+wu
      print *,"N J down = ",mjd," (",
     &    (dble(mjd)/dble(nt))*100.," %)",mjd,mju,mjd+mju

      dn=dble(nt)
      etriav = etriav/dn
      ttav = ttav/dn
      etotiav = etotiav/dn
      etotdmue = etotdmue/dn
      etotdrms = dsqrt(etotdrms/dn)
      avgediffi = avgediffi/dn
      avgedifff = avgedifff/dn
      eeiavg=eeiavg/dn
      xjavg=xjavg/dn

      bav = bav/dn

      print *
      print *,"2x<bi> = ",(2.d0*bav)," (",bmax,")"
      print *
      print *," <Etot>i = ",etotiav," eV"
      print *," <Etr>i  = ",etriav," eV"
      print *," <Eint>i = ",(etotiav-etriav)," eV"
      print *
      print *," E conservation (cm-1):"
      print *,"    MUE ",(etotdmue*evtocmi)
      print *,"    RMS ",(etotdrms*evtocmi)
      print *,"    max ",(etotdmax*evtocmi)," traj ",itetotdmax
      print *
      print *," Vi diff (cm-1): max =",maxediffi*evtocmi,
     & ", mue =",avgediffi*evtocmi," traj #",maxediffitraj
      print *," Vf diff (cm-1): max =",maxedifff*evtocmi,
     & ", mue =",avgedifff*evtocmi," traj #",maxediffftraj
      print *
      print *," Max DE (cm-1) = ",(maxd)," (traj # ",maxdit,")"
      print *," Min DE (cm-1) = ",(mind)," (traj # ",mindit,")"
      print *

        call boot(md,mdat,sd,avg,sig)
        write(6,'(a10,3f15.5)')'<DEd>',avg,sig,(sig/avg)*100.d0
        write(7,*)"! alpha = ",-avg,sig,dabs(sig/avg)*100.d0,"%"

 99    format('      ejmomq(',4(i2,','),i2,') = ',f25.10)
        print *,ejmom(0,0,0,0,0),dn
        print *
        print *,'! ZLJ = ',zlj
        print *,'! alpha = ',avg,sig

        do kj=0,ndblej
        print *,"! J' order ",kj
        write(7,*)"! J' order ",kj

        do ie=1,3
        print *,"! Pure moments order ",ie
        write(7,*)"! Pure moments order ",ie

! energy moments
        do i=0,-2,-1
        ix=i
        if (ix.eq.-2) ix=1
        tmp=ejmom(ix,ie,0,0,kj)/ejmom(ix,0,0,0,kj)
        tmp=sign(dabs(tmp*zz)**(1.d0/dble(ie)),tmp)
        write(6,99)ix,ie,0,0,kj,tmp
        write(7,99)ix,ie,0,0,kj,tmp
        enddo
! j moments
        do i=0,-2,-1
        ix=i
        if (ix.eq.-2) ix=1
        tmp=ejmom(0,0,ix,ie,kj)/ejmom(0,0,ix,0,kj)
        tmp=sign(dabs(tmp*zz)**(1.d0/dble(ie)),tmp)
        write(6,99)0,0,ix,ie,kj,tmp
        write(7,99)0,0,ix,ie,kj,tmp
        enddo

        enddo

! ej moments
        do k=1,3
        if (k.eq.1) then
        print *,"! Cross moments order",2
        write(7,*)"! Cross moments order",2
        ia=1
        ib=1
        elseif (k.eq.2) then
        print *,"! Cross moments order",3
        write(7,*)"! Cross moments order",3
        ia=2
        ib=1
        else
        ia=1
        ib=2
        endif
        do i=0,-2,-1
        ix=i
        if (ix.eq.-2) ix=1
        do j=0,-2,-1
        iy=j
        if (iy.eq.-2) iy=1
        tmp=ejmom(ix,ia,iy,ib,kj)/ejmom(ix,0,iy,0,kj)*zz
        write(6,99)ix,ia,iy,ib,kj,tmp
        write(7,99)ix,ia,iy,ib,kj,tmp
        enddo
        enddo
        enddo

        enddo



      end

      subroutine boot(nd,md,dd,avg,sig)

      implicit double precision(a-h,o-z)
      dimension dd(md),avgi(100000)
      integer*4 timeArray(3)

      call itime(timeArray)     ! Get the current time
      iseed = timeArray(1)+timeArray(2)+timeArray(3)
      x=rand(iseed)

      nsamp=100000

      avgavg=0.
      do j=1,nsamp
        tot=0
        do i=1,nd
          irand=int(rand()*dble(nd))+1;
          v=dd(irand)
          tot=tot+v
        enddo
        avg = tot/dble(nd)
        avgi(j)=avg
        avgavg=avgavg+avg
      enddo
      avgavg=avgavg/dble(nsamp)

      std=0.
      do j=1,nsamp
        std=std+(avgi(j)-avgavg)**2
      enddo
      std=std/dble(nsamp-1)
      std=dsqrt(std)

      avg=avgavg
      sig=std

      return
      end

      subroutine binw(nd,md,dd,ww)

      implicit double precision(a-h,o-z)
      parameter(maxbin=20000)
      dimension dd(md),ww(md),nbin(maxbin),wbin(maxbin)

      evtocmi = 8065.7208

        demax=dd(1)
        demin=dd(1)
        do i=1,nd
          demax=max(dd(i),demax)
          demin=min(dd(i),demin)
        enddo

        do i=1,maxbin
         nbin(i)=0
         wbin(i)=0.d0
        enddo
        binsize = (demax-demin)/dble(nb)

        demin=-4020
        demax=4020
        demin=-8000
        demax=8000
        binsize=50.
        nb=(demax-demin)/binsize
        if (nb.gt.maxbin) print *,nb,maxbin
        if (nb.gt.maxbin) stop
        demin=demin-binsize

        wtot=0.d0
        do i=1,nd
          ib = int((dd(i)-demin)/binsize)
          if (ib.lt.nb.and.ib.gt.0) then
          nbin(ib)=nbin(ib)+1
          wbin(ib)=wbin(ib)+ww(i)
          wtot=wtot+ww(i)
c          print *,dble(ib)*binsize+demin,
c     &      dble(ib+1)*binsize+demin,dd(i)
          endif
        enddo

        eup=0.d0
        eup2=0.d0
        wup=0.d0
        edown=0.d0
        edown2=0.d0
        wdown=0.d0
        e0=0.d0
        e02=0.d0
        w0=0.d0
        do i=1,nb
          a=dble(i)*binsize+demin
          b=dble(i+1)*binsize+demin
          c=(a+b)/2.d0
c          print  *,i,c,dble(nbin(i))/dble(nd),nbin(i),nd,wbin(i)/wtot
          print  *,i,c,wbin(i)/wtot
          if (c.gt.0) then
               eup=eup+wbin(i)/wtot*c
               eup2=eup2+wbin(i)/wtot*c**2
               wup=wup+wbin(i)/wtot
          elseif (c.lt.0) then
               edown=edown+wbin(i)/wtot*c
               edown2=edown2+wbin(i)/wtot*c**2
               wdown=wdown+wbin(i)/wtot
          else
               e0=e0+wbin(i)/wtot*c
               e02=e02+wbin(i)/wtot*c**2
               w0=w0+wbin(i)/wtot
          endif
        enddo

        print *,"tot down up down2 up2 ",
     &                 (e0+edown+eup)/(w0+wdown+wup),
     &                  edown/wdown,eup/wup,
     &                  edown2/wdown,eup2/wup
        print *,"tot down up down2 up2 ",
     &                 (e0+edown+eup)/(w0+wdown+wup),
     &                  edown/(wdown+.5*w0),eup/(wup+.5*w0),
     &                  edown2/(wdown+.5*w0),eup2/(wup+.5*w0)

c         stop

        return
        end

