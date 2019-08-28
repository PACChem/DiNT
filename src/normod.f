      subroutine normod(xx0,mm,nclu,nmvec,freq,ewell,
     &   repflag,nsurf,rturn,nmqn,lreadhess,ldofrag,symb,nmtype)

c Compute the normal modes using a numerical Hessian.
c The code assumes that XX0 is a minimum.
c The code returns FREQ = Frequency in 1/a.u.time for each mode
c                  EWELL = Potential energy of the well.
c                  RTURN = HO approximation of the turning point 
c                          distance for each mode (mass-scaled).

      implicit none
      include 'param.f'

c     input
      logical lreadhess,ldofrag
      integer nclu,repflag,nsurf,nmtype
      double precision xx0(3,mnat),mm(mnat)
      character*2 symb(mnat)

c     output
      double precision nmvec(3,mnat,3*mnat),freq(3*mnat)

c     local
      integer lwork,nbound
      parameter(lwork=9*mnat-1)
      integer i,j,ij,k,i1,i2,l,kl,nmax,ndim,info
      double precision xx(3,mnat),hh,ewell,evib,nmqn(3*mnat)
      double precision pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),gv1(3,mnat),gv2(3,mnat),
     & nmvec1(3*mnat,3*mnat),work(lwork),hess(3*mnat,3*mnat),
     & gvsum,rturn(3*mnat),etotvib,temp(3*mnat),tmp,rr

      if (nclu.lt.2) then
        write(6,*)"cant do normal modes for ",nclu," atoms"
        stop
      endif

c Initialize
      do i=1,3*mnat
      do j=1,3*mnat
      hess(i,j) = 0.d0
      enddo
      enddo

c Check geometry
      if (ldofrag) then
      do k=1,3
      do l=1,nclu
        xx(k,l) = xx0(k,l)
      enddo
      enddo
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
      if (repflag.eq.0) then
c       adiabatic
        ewell = pema(nsurf)
        do k=1,3
        do l=1,nclu
          gv1(k,l) = gpema(k,l,nsurf)
        enddo
        enddo
      elseif (repflag.eq.1) then
c       diabatic
        ewell = pemd(nsurf,nsurf)
        do k=1,3
        do l=1,nclu
          gv1(k,l) = gpemd(k,l,nsurf,nsurf)
        enddo
        enddo
      else
        write(6,*)"REPFLAG = ",repflag," in NORMOD"
        stop
      endif
      write(6,'(a,f19.10)')"Energy at this geometry = ",ewell*autoev
      write(6,*)
      write(6,*)"      Atom      x,y,z      Gradient (eV/A)"
      gvsum = 0.d0
      do l=1,nclu
      do k=1,3
        write(6,100)l,k,gv1(k,l)*autoev/autoang
        gvsum = gvsum + gv1(k,l)**2
      enddo
      enddo
      gvsum = dsqrt(gvsum)
      write(6,*)"Magnitude of gradient = ",gvsum*autoev/autoang," eV/A"
      write(6,*)
 100  format(2i10,f20.10)
      endif
      write(6,*)"Computing and diagonalizing Hessian..."
      write(6,*)

      IF (lreadhess) THEN
c     read in Hessin from fort.70
      write(6,*)"Reading Hessian from unit 70"
      open(70)
      do i=1,3
      do j=1,nclu
        ij = (i-1)*nclu + j
        do k=1,3
        do l=1,nclu
          kl = (k-1)*nclu + l
             read(70,*)i1,i2,hess(i1,i2)
             hess(i1,i2) = hess(i1,i2)*mu/dsqrt(mm(j)*mm(l))
             if (i1.ne.ij.or.i2.ne.kl) then
               write(6,*)"Index mismatch reading Hessian from unit 70"
               stop
             endif
        enddo
        enddo
      enddo
      enddo
      ELSE
c Compute Hessian matrix from gradients
c     stepsize
      hh = 0.01d0
c      hh = 0.001d0
      do i=1,3
      do j=1,nclu
        print *,'doing step ',i,',',j,' out of 3,',nclu
        ij = (i-1)*nclu + j
        do k=1,3
        do l=1,nclu
          xx(k,l) = xx0(k,l)
        enddo
        enddo
        xx(i,j) = xx0(i,j) + hh
        call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
        if (repflag.eq.0) then
c         adiabatic
          do k=1,3
          do l=1,nclu
            gv1(k,l) = gpema(k,l,nsurf)
          enddo
          enddo
        elseif (repflag.eq.1) then
c         diabatic
          do k=1,3
          do l=1,nclu
            gv1(k,l) = gpemd(k,l,nsurf,nsurf)
          enddo
          enddo
        else
          write(6,*)"REPFLAG = ",repflag," in NORMOD"
          stop
        endif
        xx(i,j) = xx0(i,j) - hh
        call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
        if (repflag.eq.0) then
c         adiabatic
          do k=1,3
          do l=1,nclu
            gv2(k,l) = gpema(k,l,nsurf)
          enddo
          enddo
        elseif (repflag.eq.1) then
c         diabatic
          do k=1,3
          do l=1,nclu
            gv2(k,l) = gpemd(k,l,nsurf,nsurf)
          enddo
          enddo
        else
          write(6,*)"REPFLAG = ",repflag," in NORMOD"
          stop
        endif
        open(70)
        do k=1,3
        do l=1,nclu
          kl = (k-1)*nclu + l
c         Hessian matrix, 3NCLU X 3NCLU matrix
c         data ordered (x1,x2,...,y1,...,z1,...,zNCLU)
          hess(ij,kl) = (gv1(k,l) - gv2(k,l))/(2.d0*hh)
          write(70,770)ij,kl,hess(ij,kl)
 770      format(2i10,d20.8)
c         mass-scale
          hess(ij,kl) = hess(ij,kl)*mu/dsqrt(mm(j)*mm(l))
        enddo
        enddo
      enddo
      enddo
      write(6,*)"Hessian written to unit 70"
      ENDIF

c diagonlize the hessian
      ndim = 3*nclu
      nmax = 3*mnat
      call dsyev( 'v','u',ndim,hess,nmax,freq,work,lwork,info )
c     FREQ is always ordered from smallest value to largest

      write(6,*)"   Index  Force Const (mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)
 150  format(i10,e15.5,f15.5)
      write(6,*)"Keeping 3N-6 modes with largest magnitude frequencies."
      write(6,*)

c repackage (order from biggest to smallest)
      do k=1,ndim
        temp(k) = freq(k)
      enddo
      do k=1,ndim
        freq(k) = temp(ndim-k+1)
      enddo
      do i=1,3
      do j=1,nclu
        ij = (i-1)*nclu + j
        do k=1,ndim
          nmvec(i,j,ndim-k+1) = hess(ij,k)
        enddo
      enddo
      enddo

      nbound=ndim-6   ! well
      if (nmtype.eq.1) nbound=ndim-7  ! saddle point
      if (nclu.eq.2) nbound=1 ! diatomic
      do k=1,nbound
        if (freq(k).gt.0.d0) then
c       convert force constant to units of 1/time in atomic units
        freq(k) = dsqrt(freq(k)/mu)
        else
        write(6,*)"Negative frequency for mode ",k," !"
        write(6,*)"Cant do unbound modes!"
        stop
        endif
      enddo
      if (nmtype.eq.1) then
        if (freq(ndim).lt.0.d0) then
          freq(ndim-6) = -dsqrt(-freq(ndim)/mu)
        else
          write(6,*)"Positive frequency for mode ",k," !"
          write(6,*)"Cant do bound modes!"
          stop
        endif
        do i=1,3
        do j=1,nclu
          ij = (i-1)*nclu + j
          nmvec(i,j,ndim-6) = hess(ij,1)
        enddo
        enddo
      endif

      write(6,*)"     Mode   Freq (1/fs)    ZPE (eV) Freq (cm-1) "
c      do k=1,ndim-6
      do k=1,nbound
        write(6,200)k,freq(k)/autofs,freq(k)/2.d0*autoev,freq(k)*autocmi
      enddo
      write(6,*)

      write(6,*)"     Mode      Mass-scaled normal",
     & " mode eigenvectors(Ax,Ay,Az,Bx,...)"
      do k=1,ndim-6
        write(6,200)k,((nmvec(i,j,k),i=1,3),j=1,nclu)
      enddo
      write(6,*)
 200  format(i10,100f13.5)

c compute HO turning points for each mode
      write(6,*)"Harmonic turning point info"
      write(6,*)"     Mode     Q. Number  Energy (eV) ",
     &  " Distance (mass-scaled bohr)"
      etotvib = 0.d0
      do k=1,nbound
        evib = freq(k)*(nmqn(k)+0.5d0)
        etotvib = etotvib + evib
        rturn(k) = dsqrt(2.d0*evib/(freq(k)**2*mu)) 
        write(6,300)k,nmqn(k),evib*autoev,rturn(k)
 300    format(i10,f8.2,3f12.5)
      enddo

      write(6,*)"HO approximation to the internal energy = ",
     &   etotvib*autoev,"eV"
      write(6,*)"HO approximation to the total energy    = ",
     &   (ewell+etotvib)*autoev,"eV"
      write(6,*)
      end
