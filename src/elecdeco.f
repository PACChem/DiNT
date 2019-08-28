      subroutine elecdeco(ithistraj,istep,time,nclu,nsurf,newsurf,
     &    pem,gpem,dvec,xx,pp,tau)

c compute the decoherence time according to J. Chem. Phys. 123, 064103 (2005)
      implicit none
      include 'param.f'
      include 'c_sys.f'

c input
      integer nclu,nsurf,newsurf,ithistraj,istep
      double precision time,xx(3,mnat)
      double precision pp(3,mnat),dvec(3,mnat,mnsurf,mnsurf)

c output
      double precision tau

c local
      integer i,j
      double precision pem(mnsurf),gpem(3,mnat,mnsurf),
     & sectrm,a,b,a2,ab,ra,radi,scr,tmp,deltae,dot1,dot2
      double precision pdotd,dvecmag,dg,deltaf,deltap,pavg,
     & tau_f,tau_p,phase

      double precision taup(3),taup2(3)
      common/tauprint/taup,taup2

c compute p along d = pdotd = the momentum along the nonadiabatic coupling vector
      pdotd = 0.d0
      dvecmag = 0.d0
      do i=1,3
      do j=1,nclu
        pdotd = pdotd + dvec(i,j,nsurf,newsurf)*pp(i,j)/mm(j)
        dvecmag = dvecmag + (dvec(i,j,nsurf,newsurf)**2)/mm(j)
      enddo
      enddo
      if (dvecmag.lt.1.d-40) then
        pdotd = 0.d0
      else
        dvecmag = dsqrt(dvecmag)
        pdotd = pdotd/dvecmag
      endif

c     d (the nonadiabatic coupling vector) has an arbitrary phase
c     here we define the phase of d to make p.d positive
      phase = 1.d0
      if (pdotd.lt.0.d0) then 
        phase = -1.d0
        pdotd = -pdotd
      endif

c delta e
      deltae = pem(newsurf) - pem(nsurf)

c delta F along d, with attention to mass scaling
      deltaf = 0.d0
      do i=1,3
      do j=1,nclu
        dg = (gpem(i,j,nsurf) - gpem(i,j,newsurf))/dsqrt(mm(j))
        deltaf = deltaf - phase*dg*dvec(i,j,nsurf,newsurf)/dsqrt(mm(j))
      enddo
      enddo
      deltaf = deltaf/dvecmag

c compute deltap and pavg
c deltap is the change in momentum that would occur if we hopped to the other state
c both of these are defined to be positive
        a=0.0d0
        b=0.0d0
        do i=1,3
        do j=1,nclu
          scr=dvec(i,j,nsurf,newsurf)/mm(j)
          a=a+scr*pp(i,j)
          b=b+scr*dvec(i,j,nsurf,newsurf)
        enddo
        enddo
        a2=a*a
        ab=a/b
c       sectrm = deltae/kinetic energy along d
        sectrm = 2.0d0*deltae*b/a2

        radi = 1.0d0-sectrm
c       if radi > 0 then delta > kinetic energy along d, so we can hop
c       note that for a hop down, sectrm is negative, so we can always hop
        if(radi.ge.0.0d0) then
c         a successful hop could be possible
          radi=sqrt(radi)
          ra=ab*(1.0d0-radi)
          deltap = dabs(ra)*dsqrt(b)
          if (pem(nsurf).gt.pem(newsurf)) then
c           we are in the excited state, so the avg p is current p minus 1/2 the gap
            pavg=pdotd+deltap/2.d0
          else 
c           we are in the lower energy state, so the avg p is current p plus 1/2 the gap
            pavg=pdotd-deltap/2.d0
          endif
        else
c         a successful hop is not possible
c         in this case, we set the momentum for the other surface to zero
c         so delta p is just p for the current surface
          deltap = pdotd
c         and the average p is the current p/2
          pavg = pdotd/2.d0
        endif

c         print *,pavg,pdotd,deltap

c tau_p (the delta momentum term) and tau_f (the delta force term) are 
c computed separately for diagnostic purposes
        tau_p = 2.d0*pi*dsqrt(pavg/deltap)/dabs(deltae)   ! eq 19 
        tau_f = 2.d0*pavg/(pi*deltaf)                     ! reciprical of first term of eq 16
        if (pem(nsurf).gt.pem(newsurf)) tau_f = -tau_f     ! an additional rule is required to fix
                                                         ! the arbitrary sign of delta f. We choose
                                                         ! to define things such that the imaginary
                                                         ! packet on the lower surface is "advanced"
                                                         ! relative to the packet on the upper surface.
                                                         ! There was some sensitivity to this choice for
                                                         ! the collinear NaFH system.
c PROBLEM?: When tau_p is small, you get 1/tau = 1/tau_f + sqrt(1/tau_f^2), which can be ~0 for tau_f < 0.
c Instead: When tau_p is small, then all the decoherence should arise from the DF term. There is
c                no reason that this term will cancel out. Therefore, we should use |tau_f|.
c                This change is not well tested.
        tau_f=dabs(tau_f)

c put it together
        tau = 1.d0/tau_f + dsqrt(1.d0/tau_f**2 + 1.d0/tau_p**2)  ! eq 17
        tau = 1.d0/tau                                        ! take inverse to convert from rate to time

        taup2(1)=tau_p
        taup2(2)=tau_f
        taup2(3)=tau

      return
      end
