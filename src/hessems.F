      subroutine hessems(xxmi,ppmi,vm,mmmi,symbi,nhtot,imm)

      implicit none      
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
#include <sprng_f.h>

      integer i,j,k,ii,ibrea,imm,nems,inems,nhtot,rnwell,nfq
      double precision bgj(3),bgjt,er(3),ert,tmp,tempi0,kem,
     &  wt,ratio,rn,vm,vmt,evb
      double precision xxmi(3,mnat),ppmi(3,mnat),mmmi(mnat),
     &  xxmt(3,mnat),ppmt(3,mnat),egv(3),xcm(3),pcm(3)
      character*2 symbi(mnat)
      logical laccpt
      double precision pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf)


 1082 format(a2,3f15.5)
 1083 format(i5,50f15.5)
      tempi0=1/kb

      inems=0
      nems=0
      do while (nems.lt.nemstot(imm))
        ibrea=0
        do while (ibrea.lt.nhtot+1)
          inems = inems+1
          write(6,*)
          write(6,*)"EMS TRIAL GEOMETRY ",inems
c         pick random quantum levels for each mode from a thermal dist
          call vibwells(xxmt,ppmt,emicr(imm),imm)
c         populate normal modes and select trial coordinates
          call popnorm(xxmt,ppmt,mmmi,natom(imm),imm)
c         place CoM at origin
          do j=1,3
            xcm(j)=0.d0
          enddo
          do i=1,natom(imm)
          do j=1,3
            xcm(j) = xcm(j)+mmmi(i)*xxmt(j,i)
          enddo
          enddo
          do i=1,natom(imm)
          do j=1,3
            xxmt(j,i) = xxmt(j,i)-xcm(j)/mmag(imm)
          enddo
          enddo
c          call pot(symbi,xt,yt,zt,vmt,emin,dxt,dyt,dzt,
c     &          natom(imm),natom(imm))
          call getpem(xxmt,natom(imm),pema,pemd,gpema,gpemd,dvec,symbi)
          if (repflag.eq.0) then
c          adiabatic
           vmt = pema(nsurf)
          elseif (repflag.eq.1) then
c          diabatic
           vmt = pemd(nsurf,nsurf)
          else
           write(6,*)"REPFLAG = ",repflag," in NORMOD"
           stop
          endif
c         correct for this atom group's zero
          write (6,*)"vmt ezero ezeroim",vmt,ezero,ezeroim(imm)
          vmt = vmt + ezero - ezeroim(imm)
c          call detwell(vmt,xxmt,welli)
c          if (welli.ne.rnwell) then
c            write(6,*)"randomwell=",rnwell," welli=",welli
c            do i=1,natom(imm)
c              write(6,*)(xxmt(k,i)*autoang,k=1,3)
c            enddo
c            go to 201
c          endif
          write(6,*)"Trial PE ",vmt*autoev," eV"
          if (vmt.gt.emicr(imm)) then
             write(6,*)"Potential energy (",vmt*autoev," eV) > E (",
     &         emicr(imm)*autoev," eV)"
             laccpt = .false.
             write(6,*)"Trial geometry REJECTED"
             write(6,*)"Ensemble size = ",nems
             go to 11
          elseif (vmt.lt.0.d0) then
            write(6,*)"Potential energy (",vmt*autoev,
     &        " eV) < zero of energy"
             laccpt = .false.
             write(6,*)"Trial geometry REJECTED"
             write(6,*)"Ensemble size = ",nems
             go to 11
          endif
c     calculate mom of inertia
          call ange(xxmt,ppmt,mmmi,natom(imm),egv,bgj,bgjt,er,ert)
c     calculate weighting factor (J=0 ensemble)
          wt = (emicr(imm)-vmt)**((3.d0*dble(natom(imm))-8.d0)/2.d0)
          wt = wt/dsqrt(egv(1)*egv(2)*egv(3))
c     accept structure or not
          if (wt.gt.emsw(imm)) then
            laccpt = .true.
            write(6,*)"Wt(i) = ",wt
            write(6,*)"Wt(i-1) = ",emsw(imm)
            write(6,*)"Trial geometry ACCEPTED"
          else
            rn=sprng(rng_stream)
            ratio=wt/emsw(imm)
            write(6,*)"Wt(i) = ",wt
            write(6,*)"Wt(i-1) = ",emsw(imm)
            write(6,*)"ratio   = ",ratio
            if (ratio.ge.rn) then
              laccpt = .true.
              write(6,*)"Trial geometry ACCEPTED"
            else
              laccpt = .false.
              write(6,*)"Trial geometry REJECTED"
              write(6,*)"Ensemble size = ",nems
            endif
          endif
c     output accepted coordinates
 11       if (laccpt) then
            ibrea = ibrea+1
            nems = nems+1
            write(6,*)"Ensemble size = ",nems
            vm=vmt
            emsw(imm)=wt
            do i=1,natom(imm)
            do j=1,3
              xxmi(j,i) = xxmt(j,i)
            enddo
            enddo
          endif
        enddo

c     generate random momenta for accepted coordinates
        call rantherm(ppmi,mmmi,natom(imm),tempi0)
c     remove angular momentum
        call noang(xxmi,ppmi,mmmi,natom(imm))
c     remove CoM momentum
        do j=1,3
          pcm(j)=0.d0
        enddo
        do i=1,natom(imm)
        do j=1,3
          pcm(j) = pcm(j)+ppmi(j,i)
        enddo
        enddo
        do i=1,natom(imm)
        do j=1,3
          ppmi(j,i) = ppmi(j,i)-mmmi(i)*pcm(j)/mmag(imm)
        enddo
        enddo
c     scale momenta to KE=Etot-PE
        call gettemp(ppmi,mmmi,natom(imm),tmp,kem)
        do i=1,natom(imm)
        do j=1,3
          ppmi(j,i) = ppmi(j,i)*dsqrt((emicr(imm)-vm)/kem)
        enddo
        enddo
        if (lwrite(86)) then
          write(86,*)natom(imm)
          write(86,1083)nems,vm
          do i=1,natom(imm)
            write(86,1082)symbi(i),(xxmi(j,i)*autoang,j=1,3)
          enddo
        endif
        if (lwrite(87)) then
          write(87,*)natom(imm)
          write(87,1083)nems,vm
          do i=1,natom(imm)
            write(87,1082)symbi(i),(ppmi(j,i),j=1,3)
          enddo
        endif
        if (lwrite(88)) write(88,REC=nems)natom(imm),nems,vm,
     &   (symbi(i),i=1,natom(imm)),
     &   ((xxmi(j,i)*autoang,j=1,3),i=1,natom(imm))
        if (lwrite(89)) write(89,REC=nems)natom(imm),nems,vm,
     &   (symbi(i),i=1,natom(imm)),
     &   ((ppmi(j,i),j=1,3),i=1,natom(imm))
      enddo

      return

      end

