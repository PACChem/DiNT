      program lsfit
c
      implicit double precision (a-h,o-z)
c
      include 'param.inc'
      parameter(npp=maxdata)
      parameter(mpp=maxdata)
      parameter(autocmi=219474.63067d0)
      dimension u(mpp,npp),v(npp,npp),w(npp)
      dimension coef(maxterm),sig(maxdata)
      dimension basis(maxterm)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension vv2(maxdata),rrr2(maxdata,maxpair)
      dimension rcom2(maxdata),xprint(50,20)
      dimension x(maxatom),y(maxatom),z(maxatom)

      character*2 dum

      common/foox/rrr,nncoef

        write(6,*)'initializing data'
        open(7,file='ai.dat')
        read(7,*)ndat2,natom1
        write(6,*)'Reading ',ndat2,' data '
        if (ndat2.gt.maxdata) then
          write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata,
     &          ". Change in LSS-PIP."
          stop
        endif

        vvmin=1.d10
        do i=1,ndat2
           read(7,*)natom
           read(7,*)dum,dum,vv2(i)
           if (vv2(i).lt.vvmin) vvmin=vv2(i)
           do j=1,natom
             read(7,*)dum,x(j),y(j),z(j)
           enddo
           ii=0
           do j=1,natom1  ! frag 1
           do k=natom1+1,natom  ! He
           ii=ii+1
           ij=ii
          rrr2(i,ij)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
           enddo  
           enddo  
        enddo  

c fitting
       open(56,file="coef.dat")
       cut0=8000
       cut1=4000.
       cut2=-0.5
       epsilon=700.d0
       read(5,*)cut0,cut1,cut2
       read(5,*)epsilon

        ndat=0
        do i=1,ndat2
        if (vv2(i).lt.cut0) then
        ndat=ndat+1
        vv(ndat)=vv2(i)
        sig(ndat)=1.d0/(epsilon/(vv(ndat)-vvmin+epsilon))
        do k=1,natom1*(natom-natom1)
        rrr(ndat,k)=rrr2(i,k)
        enddo
        write(20,'(i10,8f18.8)')ndat,(rrr(ndat,k),k=1,6),
     &    vv(ndat),1.d0/sig(ndat)
        endif
        enddo

        write(6,*)'Using ',ndat,' data '

        call prepot
        ncoef=nncoef
        print *,ncoef," coefficients"

        call svdfit(vv,sig,ndat,coef,ncoef,u,v,w,mpp,npp,chisq)

        err=0.d0
        err2=0.d0 
        errx=0.d0
        errx2=0.d0
        erry=0.d0 
        erry2=0.d0 
        errz=0.d0 
        errz2=0.d0 
        nn1=0
        nn2=0
        nn3=0
        wn=0.d0
        wnx=0.d0
        wny=0.d0
        vvxm=10.d0
        vvim=10.d0
        do i=1,ndat
           call funcs1(i,basis,ncoef) 
           vvx=0.d0
           do j=1,ncoef
              vvx=vvx+coef(j)*basis(j)
         if (i.eq.1) write(6,'(a,i5,a,e20.10)')
     &    '       coef(',j,') = ',coef(j)
         if (i.eq.1) write(56,'(i5,e20.10)')j,coef(j)
c           write(21,'(2i10,5e18.8)')i,j,basis(j)
           enddo
c         if (i.le.1000) write(6,'(i7,99e20.10)')i,vvx,vv(i),1.d0/sig(i)
c         write(22,'(i10,10f18.8)')i,(rrr(i,k),k=1,6),
c     &    vvx,vv(i),1.d0/sig(i)
         err=err+dabs(vvx-vv(i))
         err2=err2+(vvx-vv(i))**2/sig(i)**2
         wn=wn+1.d0/sig(i)**2
         if (vv(i).lt.cut1) then
         errx=errx+dabs(vvx-vv(i))
         errx2=errx2+(vvx-vv(i))**2/sig(i)**2
         wnx=wnx+1.d0/sig(i)**2
         nn1=nn1+1
         endif
         if (vv(i).lt.cut2) then
         erry=erry+dabs(vvx-vv(i))
         erry2=erry2+(vvx-vv(i))**2/sig(i)**2
         wny=wny+1.d0/sig(i)**2
         nn2=nn2+1
         endif
         if (vv(i).lt.vvim) then
            vvim2=vvx
            vvim=vv(i)
         endif
         if (vvx.lt.vvxm) then
            vvxm=vvx
            vvxm2=vv(i)
         endif
        enddo
        close(56)
c        err=err/dble(ndat)
c        err2=dsqrt(err2/dble(ndat))
c        errx=errx/dble(nn1)
c        errx2=dsqrt(errx2/dble(nn1))
c        erry=erry/dble(nn2)
c        erry2=dsqrt(erry2/dble(nn2))
        err=err/wn
        err2=dsqrt(err2/wn)
        errx=errx/wnx
        errx2=dsqrt(errx2/wnx)
        erry=erry/wny
        erry2=dsqrt(erry2/wny)
         print *,'ERRORS'
         print *,' < ',cut0,ndat,err2
         print *,' < ',cut1,nn1,errx2
         print *,' < ',cut2,nn2,erry2
         print *,"         fit  ai"
         print *,"fit min",vvxm,vvxm2,vvxm-vvxm2
         print *,"ai  min",vvim2,vvim,vvim2-vvim

c test
        open(11,file='fit.dat')
        close(7)
        open(7,file='ai2.dat')
        rewind(7)
        read(7,*)ndat2,natom1
        write(6,*)'Reading ',ndat2,' data '
        if (ndat2.gt.maxdata) then
          write(6,*)"ndat = ",ndat2," while maxdat = ",maxdata,
     &          ". Change in FUNC."
          stop
        endif

        do i=1,ndat2
           read(7,*)natom
           read(7,*)dum,rcom2(i),vv2(i)
           sig(i)=1.d0/(epsilon/(vv2(i)-vvmin+epsilon))
           do j=1,natom
             read(7,*)dum,x(j),y(j),z(j)
           enddo
           ii=0
           do j=1,natom1
           do k=natom1+1,natom
           ii=ii+1
           ij=ii
          rrr(i,ij)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
           enddo
           enddo
        enddo

        rlast=0.d0
        k=2
        j=0
        err2=0.d0
        ndat=0
        wn=0.d0
        do i=1,ndat2
           j=j+1
           call funcs1(i,basis,ncoef)
           vvx=0.d0
           do l=1,ncoef
              vvx=vvx+coef(l)*basis(l)
           enddo
         if (vv2(i).lt.cut1) err2=err2+(vvx-vv2(i))**2/sig(i)**2
         if (vv2(i).lt.cut1) wn=wn+1.d0/sig(i)**2
         if (vv2(i).lt.cut1) ndat=ndat+1
        if (rcom2(i).lt.rlast) k=k+2
        if (rcom2(i).lt.rlast) j=1
        if (rcom2(i).ne.xprint(j,1).and.i.ne.j) 
     &     print *,"com distance mismatch"
        rlast=rcom2(i)
        xprint(j,1)=rcom2(i)
        xprint(j,k)=vv2(i)
        xprint(j,k+1)=vvx
        enddo
c        err2=dsqrt(err2/dble(ndat))
        err2=dsqrt(err2/wn)
        print *,' test set < ',cut1,ndat,err2

         do i=1,j
           write(11,'(100f15.6)')(xprint(i,k),k=1,13)
         enddo

         print *,"errors",errx2,err2
 
          end
