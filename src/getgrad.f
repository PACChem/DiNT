      subroutine getgrad(xx,pp,nsurf,v,cre,cim,gv,gcre,gcim,nclu,
     &  phop,dmag,dvec,pem,gpem,phase)

c Calls GETPEM once, and then computes several important variables:
c V = Potential energy
c GV = Gradient of the potential energy
c DMAG = Magnitude of the nonadiabatic coupling vector (for CSDM)
c GCRE and GCIM = Time derivatives of the real and imaginary parts
c   of the electronic variables, including the DM terms for the
c   SCDM and CSDM methods.
c PHOP = Hopping probability (for surface hoppping) or switching
c   probability (for SCDM and CSDM) divided by the stepsize.
c PEM = Array of adiabatic or diabatic potential energies.
c GPEM = Gradients of the adiabatic or diabatic potential energies.

      implicit none
      include 'param.f'
      include 'c_sys.f'

c temp
      double precision taup(3),taup2(3)
      common/tauprint/taup,taup2

      integer i,j,k,l,m,nsurf,ii,nclu,i2,j2
      double precision xx(3,mnat),pp(3,mnat),gv(3,mnat),v,
     & gcim(2*mnsurf),gcre(2*mnsurf),cre(2*mnsurf),cim(2*mnsurf)
      double precision pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),gpem(3,mnat,mnsurf),
     & rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),vdotd(mnsurf,mnsurf),
     & phop(mnsurf),bb,pem(mnsurf),phase(mnsurf)

      double precision ppvib(3,mnat),svec(3,mnat,mnsurf),
     & Cparam,E0param,es,taui(mnsurf),ps(mnsurf),tmp,dum,dvecmag,pdotd,
     & u2b2(2,2),gu2b2(3,mnat,2,2),dvec2b2(3,mnat),gcred(mnsurf),
     & gcimd(mnsurf),vd(mnsurf),pdeco(3,mnat),smagms,
     & rhorc(mnsurf,mnsurf),rhoic(mnsurf,mnsurf),tmpr,
     & tmpi,dsum2,dmag,sntmp,cstmp,rhotmp,
     & rhorc2(mnsurf,mnsurf),rhoic2(mnsurf,mnsurf)

c msx search variables
       integer nupper
       double precision x1mag,x1(3,mnat),gupp(3,mnat),
     &   dot1,egap,gf(3,mnat),gg(3,mnat),msx_grad(3,mnat)

c      print *,((xx(i,j),i=1,3),j=1,3),((pp(i,j),i=1,3),j=1,3)
c ZERO
      do i=1,nsurft
        phop(i) = 0.d0
      enddo

c GET ENERGIES
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symbol)

c POTENTIAL ENERGIES AND GRADIENTS
      IF (METHFLAG.EQ.0.OR.METHFLAG.EQ.1.OR.METHFLAG.EQ.5) THEN
c     single surface propagation or surface hopping calculation

      if (repflag.eq.0) then
c       adiabatic
        v = pema(nsurf)
        do i=1,3
        do j=1,nclu
          gv(i,j)=gpema(i,j,nsurf)
        enddo
        enddo
      else if (repflag.eq.1) then
c       diabatic
        v = pemd(nsurf,nsurf)
        do i=1,3
        do j=1,nclu
          gv(i,j)=gpemd(i,j,nsurf,nsurf)
        enddo
        enddo
      else
        write(6,*)"REPFLAG = ",repflag," in GETPOT"
        stop
      endif

      ELSE IF (METHFLAG.EQ.2.OR.METHFLAG.EQ.3.OR.METHFLAG.EQ.4) THEN
c     Semiclassical Ehrenfest and DM methods

c     convert electronic variables to density matrix
      call getrho(cre,cim,rhor,rhoi,nsurft)

      if (repflag.eq.1) then

c       diabatic representation
        v = 0.d0
        do k=1,nsurft
        do l=1,nsurft
c integrate phase angle separately
        tmp = phase(l)-phase(k)
        sntmp = dsin(tmp)
        cstmp = dcos(tmp)
        rhotmp = (rhor(k,l)*cstmp+rhoi(k,l)*sntmp)
c       imaginary terms cancel exactly
        v=v+rhotmp*pemd(k,l)
c end
c integrate whole coefficient
c        v = v + rhor(k,l)*pemd(k,l)
c end
        enddo
        enddo
        do 15 i=1,3
        do 15 j=1,nclu
        gv(i,j) = 0.d0
        do 15 k=1,nsurft
        do 15 l=1,nsurft
c integrate phase angle separately
        tmp = phase(l)-phase(k)
        sntmp = dsin(tmp)
        cstmp = dcos(tmp)
        rhotmp = (rhor(k,l)*cstmp+rhoi(k,l)*sntmp)
        gv(i,j)=gv(i,j)+rhotmp*gpemd(i,j,k,l)
c end
c integrate whole coefficient
c        gv(i,j) = gv(i,j) + rhor(k,l)*gpemd(i,j,k,l)
c end
  15    continue

      else if (repflag.eq.0) then

c       adiabatic representation
        v = 0.d0
        do k=1,nsurft
        v = v + rhor(k,k)*pema(k)
        enddo
        do 20 i=1,3
        do 20 j=1,nclu
        gv(i,j) = 0.d0
        do 20 k=1,nsurft
        gv(i,j) = gv(i,j) + rhor(k,k)*gpema(i,j,k)
        do 20 l=1,nsurft
c integrate whole coefficient
c        gv(i,j) = gv(i,j) - 2.d0*rhor(k,l)*dvec(i,j,k,l)*pema(k)
c end
c integrate phase angle separately
        tmp = phase(k)-phase(l)
        sntmp = dsin(tmp)
        cstmp = dcos(tmp)
        rhotmp = (rhor(k,l)*cstmp-rhoi(k,l)*sntmp)
        gv(i,j) = gv(i,j) - 2.d0*rhotmp*dvec(i,j,k,l)*pema(k)
c end
  20    continue

      else
        write(6,*)"REPFLAG = ",repflag," in GETGRAD"
        stop
      endif

      ELSE

      write(6,*)"METHFLAG = ",methflag," is not allowed in GETPOT"
      stop

      ENDIF

c FOR USE BY CSDM, MAGNITUDE OF D
c      if (methflag.eq.4) then
c       compute sum of magnitude of coupled DVECs
        dmag = 0.d0
        do k=1,nsurft
          dsum2 = 0.d0
          if (k.ne.nsurf) then
            do i=1,3
            do j=1,nclu
              dsum2 = dsum2 + dvec(i,j,k,nsurf)**2
            enddo
            enddo
          endif
          dsum2 = max(0.d0,dsum2)
          dsum2 = dsqrt(dsum2)
          dmag = dmag + dsum2
        enddo
c      endif

c TIME-DERIVATIVES OF THE ELECTRONIC COORDINATES
c     All methods, DM methods add terms later
      if (nsurft.eq.1) then
        gcre(1) = 0.d0
        gcim(1) = 0.d0
      else
        if (repflag.eq.1) then
c         diabatic rep
          do i=1,nsurft
            gcim(i) = 0.d0
            gcre(i) = 0.d0
            do j=1,nsurft
c integrate whole coefficient
c              gcre(i) = gcre(i) - cim(j)*pemd(i,j)
c              gcim(i) = gcim(i) + cre(j)*pemd(i,j)
c end
c integrate phase angle separately
              if (i.ne.j) then
                tmp = phase(j)-phase(i)
                sntmp = dsin(tmp)
                cstmp = dcos(tmp)
                tmpr = -(sntmp*cre(j)-cstmp*cim(j))*pemd(i,j)
                tmpi = -(sntmp*cim(j)+cstmp*cre(j))*pemd(i,j)
                gcre(i) = gcre(i) + tmpr
                gcim(i) = gcim(i) + tmpi
                if (i.eq.nsurf) then
                  phop(j) = -2.d0*(tmpr*cre(i)+tmpi*cim(i))
     &               /(cre(i)**2+cim(i)**2)
                endif 
              endif 
c end
            enddo
          enddo
        elseif (repflag.eq.0) then
c         adiabatic rep
c         compute velocity-dot-dvec
          do 10 k=1,nsurft
          do 10 l=1,nsurft
          vdotd(k,l) = 0.d0
          do 10 i=1,3
          do 10 j=1,nclu
            vdotd(k,l) = vdotd(k,l) + dvec(i,j,k,l)*pp(i,j)/mm(j)
  10      continue
          do i=1,nsurft
c integrate whole coefficient
c            gcre(i) =  cim(i)*pema(i)
c            gcim(i) = -cre(i)*pema(i)
c end
c integrate phase angle separately
            gcre(i) = 0.d0
            gcim(i) = 0.d0
            do j=1,nsurft
c end
c integrate whole coefficient
c              gcre(i) =  gcre(i) - cre(j)*vdotd(i,j)
c              gcim(i) =  gcim(i) - cim(j)*vdotd(i,j)
c end
c integrate phase angle separately
              if (i.ne.j) then
                tmp = phase(j)-phase(i)
                sntmp = dsin(tmp)
                cstmp = dcos(tmp)
                tmpr = -(cstmp*cre(j)+sntmp*cim(j))*vdotd(i,j)
                tmpi = -(cstmp*cim(j)-sntmp*cre(j))*vdotd(i,j)
                gcre(i) = gcre(i) + tmpr
                gcim(i) = gcim(i) + tmpi
                if (i.eq.nsurf) then
                  phop(j) = -2.d0*(tmpr*cre(i)+tmpi*cim(i))
     &               /(cre(i)**2+cim(i)**2)
                endif
              endif
c end
            enddo
          enddo
        else
          write(6,*)"REPFLAG = ",repflag," in GETGRAD"
          stop
        endif
      endif

      if (methflag.eq.4) then
c       CSDM
c       propagate coherent part of the electronic coordinates
c       put these quantities in CRE and CIM after real coefficients
c       the integrator will integrate the whole thing
c       Also compute phop using the coherent variables.  Phop was computed
c       above using the decoherent variables, here we overwrite it for CSDM.
        if (repflag.eq.1) then
c         diabatic rep
          do i=1,nsurft
            i2 = i + nsurft
            gcim(i2) = 0.d0
            gcre(i2) = 0.d0
            do j=1,nsurft
              j2 = j + nsurft
c integrate whole coefficient
c              gcre(i2) = gcre(i2) - cim(j2)*pemd(i,j)
c              gcim(i2) = gcim(i2) + cre(j2)*pemd(i,j)
c end
c integrate phase angle separately
              if (i.ne.j) then
                tmp = phase(j)-phase(i)
                sntmp = dsin(tmp)
                cstmp = dcos(tmp)
                tmpr=-(sntmp*cre(j2)-cstmp*cim(j2))*pemd(i,j)
                tmpi=-(sntmp*cim(j2)+cstmp*cre(j2))*pemd(i,j)
                gcre(i2)=gcre(i2)+tmpr
                gcim(i2)=gcim(i2)+tmpi
                if (i.eq.nsurf) then
                  phop(j) = -2.d0*(tmpr*cre(i2)+tmpi*cim(i2))
     &               /(cre(i2)**2+cim(i2)**2)
                endif
              endif
c end
            enddo
          enddo
        elseif (repflag.eq.0) then
          do i=1,nsurft
            i2 = i + nsurft
c integrate whole coefficient
c            gcre(i2) =  cim(i2)*pema(i)
c            gcim(i2) = -cre(i2)*pema(i)
c end
c integrate phase angle separately
            gcim(i2) = 0.d0
            gcre(i2) = 0.d0
c end
            do j=1,nsurft
              j2 = j + nsurft
c integrate whole coefficient
c              gcre(i2) =  gcre(i2) - cre(j2)*vdotd(i,j)
c              gcim(i2) =  gcim(i2) - cim(j2)*vdotd(i,j)
c end
c integrate phase angle separately
              if (i.ne.j) then
                tmp = phase(j)-phase(i)
                sntmp = dsin(tmp)
                cstmp = dcos(tmp)
                tmpr=-(cstmp*cre(j2)+sntmp*cim(j2))*vdotd(i,j)
                tmpi=-(cstmp*cim(j2)-sntmp*cre(j2))*vdotd(i,j)
                gcre(i2)=gcre(i2)+tmpr
                gcim(i2)=gcim(i2)+tmpi
                if (i.eq.nsurf) then
                  phop(j) = -2.d0*(tmpr*cre(i2)+tmpi*cim(i2))
     &               /(cre(i2)**2+cim(i2)**2)
                endif
              endif
c end
            enddo
          enddo
        else
          write(6,*)"REPFLAG = ",repflag," in GETGRAD"
          stop
        endif
      endif

c SPECIAL DM TERMS
c tmp
c tmp
c tmp
      if (methflag.eq.3.or.methflag.eq.4) then
c tmp
c tmp
c tmp
c tmp
c       DM methods
c       compute pvib
c        if (nclu.eq.3) then
c        call pjsplit(xx,pp,mm,ppvib)
c        else
c       THIS PART (FOR MORE THAN 3 ATOMS) HASN'T BEEN TESTED
        do i=1,3
        do j=1,nclu
          ppvib(i,j) = pp(i,j)
        enddo
        enddo
        call noang(xx,ppvib,mm,nclu)
c        endif

        do k=1,nsurft
        if (k.eq.nsurf) then
c       skip
        else
c         compute "reduced" nonadiabatic coupling for states k and l
          if (repflag.eq.1) then
c           diabatic
            u2b2(1,1) = pemd(k,k)
            u2b2(1,2) = pemd(k,nsurf)
            u2b2(2,1) = pemd(nsurf,k)
            u2b2(2,2) = pemd(nsurf,nsurf)
            do i=1,3
            do j=1,nclu
              gu2b2(i,j,1,1) = gpemd(i,j,k,k)
              gu2b2(i,j,1,2) = gpemd(i,j,k,nsurf)
              gu2b2(i,j,2,1) = gpemd(i,j,nsurf,k)
              gu2b2(i,j,2,2) = gpemd(i,j,nsurf,nsurf)
            enddo
            enddo
            call getdvec2(nclu,u2b2,gu2b2,dvec2b2)
          else
c           adiabatic rep
c           just use real DVEC
            do i=1,3
            do j=1,nclu
              dvec2b2(i,j) = dvec(i,j,k,nsurf)
            enddo
            enddo
          endif

c         compute PP dot DVEC2b2 / |DVEC2b2|
          pdotd = 0.d0
          dvecmag = 0.d0
          do i=1,3
          do j=1,nclu
            pdotd = pdotd + dvec2b2(i,j)*pp(i,j)/mm(j)
            dvecmag = dvecmag + (dvec2b2(i,j)**2)/mm(j)
          enddo
          enddo
          if (dvecmag.lt.1.d-40) then
          pdotd = 0.d0
          else
          dvecmag = dsqrt(dvecmag)
          pdotd = pdotd/dvecmag
          endif

c         compute SVEC(i,j,k) = DVEC2b2*PP_DVEC+PPVIB
          smagms = 0.d0
          do i=1,3
          do j=1,nclu
            svec(i,j,k)=dvec2b2(i,j)*pdotd+ppvib(i,j)
            smagms = smagms + (svec(i,j,k)**2)/mm(j)
          enddo
          enddo
          if (smagms.le.0.d0) then
c         problem
          go to 99
          endif
          smagms = dsqrt(smagms)
c         Note:  SVEC is not mass-scaled
c                SMAGMS is the magnitude of the mass-scaled SVEC

c         compute energy along s vector (ES)
          ps(k) = 0.d0
          do i=1,3
          do j=1,nclu
c           use mass-scaling, SMAGMS normalizes the mass-scaled SVEC
c           when computing ES
            ps(k) = ps(k) + svec(i,j,k)*pp(i,j)/mm(j)
          enddo
          enddo
          Es = (ps(k)**2)/(2.d0*smagms**2)
   
c         compute tau matrix
          Cparam = 1.d0
c          E0param = 0.1d0
          E0param = 0.01d0
          if (repflag.eq.0) tmp = dabs(pema(k)-pema(nsurf))
          if (repflag.eq.1) tmp = dabs(pemd(k,k)-pemd(nsurf,nsurf))
c         inverse of tau
          taui(k) = tmp/(Cparam+E0param/Es)
c hack
          taup(1) = Cparam/tmp
          taup(2) = E0param/(Es*tmp)
          taup(3) = 1.d0/taui(k)
        endif
        taui(nsurf) = 0.d0
        enddo

c tmp
c tmp
c tmp
c        if (.false.) then
c tmp
c tmp
c tmp
c       compute decoherent time-derivs of electronic coords
        gcred(nsurf) = 0.d0
        gcimd(nsurf) = 0.d0
        do k=1,nsurft
          if (k.ne.nsurf) then
          gcred(nsurf) = gcred(nsurf) + rhor(k,k)*taui(k)
          gcimd(nsurf) = gcimd(nsurf) + rhor(k,k)*taui(k)
          gcred(k) = -0.5d0*taui(k)*cre(k)
          gcimd(k) = -0.5d0*taui(k)*cim(k)
          endif
        enddo
        if (rhor(nsurf,nsurf).ne.0.d0) then
        gcred(nsurf) = gcred(nsurf)*0.5d0*cre(nsurf)/rhor(nsurf,nsurf)
        gcimd(nsurf) = gcimd(nsurf)*0.5d0*cim(nsurf)/rhor(nsurf,nsurf)
        else
        go to 99
        endif

c       compute time-deriv of decoherent potential
        do k=1,nsurft
        vd(k) = 0.d0
        enddo
        if (repflag.eq.1) then
c         diabatic
          do k=1,nsurft
c           diagonal terms
            if (k.ne.nsurf) then
              vd(k) = vd(k) - rhor(k,k)*taui(k)*pemd(k,k)
            elseif (k.eq.nsurf) then
              do l=1,nsurft
                if (l.ne.k) then
                  vd(l) = vd(l) + rhor(l,l)*taui(l)*pemd(k,k)
                endif
              enddo
            else
              write(6,*)"Strange things in GETGRAD (1)..."
              stop
            endif
            do l=1,nsurft
c             off-diagonal terms
              tmp = phase(l)-phase(k)
              sntmp = dsin(tmp)
              cstmp = dcos(tmp)
              tmpr=rhor(k,l)*cstmp+rhoi(k,l)*sntmp
              if (k.ne.nsurf.and.l.ne.nsurf.and.k.ne.l) then
                vd(k)=vd(k)-0.5d0*(taui(k)+taui(l))*tmpr*pemd(k,l)
              elseif (k.eq.nsurf.and.l.ne.nsurf) then
                vd(l) = vd(l)-0.5d0*taui(l)*tmpr*pemd(k,l)
                do m=1,nsurft
                  if (m.ne.nsurf) then
                    vd(m) = vd(m) + 0.5d0*rhor(m,m)*taui(m)/
     &                      rhor(nsurf,nsurf)*tmpr*pemd(k,l)
                  endif
                enddo
              elseif (l.eq.nsurf.and.k.ne.nsurf) then
                vd(k) = vd(k)-0.5d0*taui(k)*tmpr*pemd(k,l)
                do m=1,nsurft
                  if (m.ne.nsurf) then
                    vd(m) = vd(m) + 0.5d0*rhor(m,m)*taui(m)/
     &                      rhor(nsurf,nsurf)*tmpr*pemd(k,l)
                  endif
                enddo
              endif
            enddo
          enddo
        else
c         adiabatic
          do k=1,nsurft
c           diagonal terms
            if (k.ne.nsurf) then
              vd(k) = vd(k) - rhor(k,k)*taui(k)*pema(k)
            elseif (k.eq.nsurf) then
              do l=1,nsurft
                if (l.ne.k) then
                  vd(l) = vd(l) + rhor(l,l)*taui(l)*pema(k)
                endif
              enddo
            else
              write(6,*)"Strange things in GETGRAD (1)..."
              stop
            endif
          enddo
        endif

c       compute decoherent force
        do i=1,3
        do j=1,nclu
        pdeco(i,j) = 0.d0
        do k=1,nsurft
        if (k.ne.nsurf) then
c       magnitude of svec cancels magnitude of svec in ps(k)
        pdeco(i,j)=pdeco(i,j)-vd(k)*svec(i,j,k)/ps(k)
        endif
        enddo
        enddo
        enddo

c       add decoherent terms to electronic probs
        do k=1,nsurft
        gcre(k) = gcre(k) + gcred(k)
        gcim(k) = gcim(k) + gcimd(k)
        enddo

c       add force to grad v
c       add the negative so it is like a gradient, its negative will be taken
c       in XPTOY

        do i=1,3
        do j=1,nclu
        gv(i,j) = gv(i,j) - pdeco(i,j)
        enddo
        enddo
      endif

   99 continue

c save PEM as the diagonal elements of whichever represenation we are using
      do i=1,nsurft
        if (repflag.eq.0) then
        pem(i) = pema(i)
        do j=1,3
        do k=1,nclu
          gpem(j,k,i) = gpema(j,k,i)
        enddo
        enddo
        else
        pem(i) = pemd(i,i)
        do j=1,3
        do k=1,nclu
          gpem(j,k,i) = gpemd(j,k,i,i)
        enddo
        enddo
        endif
      enddo

c overwrite gradients if we are doing a MSX search
      if (tflag(1).eq.4) then

c     get upper state index
      nupper=2
      if (pem(1).gt.pem(2)) nupper=1
      egap=pem(1)-pem(2)

c     get x1 and gupper

      x1mag = 0.d0
      do i=1,nclu
      do j=1,3
      x1(j,i)   = gpem(j,i,1)-gpem(j,i,2)    ! grad(E1-E2)
      x1mag     = x1mag+x1(j,i)**2
      gupp(j,i) = gpem(j,i,nupper)           ! grad of upper surface
      enddo
      enddo
      x1mag = dsqrt(x1mag)                   ! magnitude of grad(E1-E2)

c     for zero nonadiabatic coupling...
c     compute dot prod
      dot1=0.d0
      do i=1,nclu
      do j=1,3
        dot1=dot1+x1(j,i)*gupp(j,i)/x1mag    ! grad upper . x1/x1mag
      enddo
      enddo

c     project out componenets in the direction of x1
      do i=1,nclu
      do j=1,3
      gf(j,i) = 2.d0*egap*x1(j,i)/x1mag      ! f = 2*(E1-E2)*x1/x1mag
      gg(j,i) = gupp(j,i)-dot1*x1(j,i)/x1mag ! g = grad upper projected onto plane perp to f
      msx_grad(j,i) = gg(j,i)+gf(j,i)            ! total gradient
      enddo
      enddo

c     replace gradient with msx_grad
      do i=1,nclu
      do j=1,3
        gv(j,i)=msx_grad(j,i)
      enddo
      enddo

      endif

      return
      end
