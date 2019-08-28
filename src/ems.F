      subroutine ems(xxmi,ppmi,vm,mmmi,symbi,nhtot,imm)

c Efficient Microcanonical Sampling (EMS)
c See H. W. Schranz, S. Nordhlm, and G. Nyman, J Chem Phys 94, 1487 (1991)
c     G. Nyman, S. Nordholm, and H. W. Schranz, J Chem Phys 93, 6767 (1990)

      implicit none      
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
#include <sprng_f.h>

      integer i,j,ii,ibrea,imm,nems,irn,inems,nhtot,rnwell
      double precision molpe,bgj(3),bgjt,er(3),ert,tmp
      double precision wt,ratio,rn,vm,vmt,tempi0,kem
      double precision xxmi(3,mnat),ppmi(3,mnat),rrmi(6),
     & mmmi(mnat),xxmt(3,mnat),ppmt(3,mnat),
     & egv(3),xcm(3),pcm(3),
     & xt(mnat),yt(mnat),zt(mnat),dxt(mnat),dyt(mnat),dzt(mnat)
      character*2 symbi(mnat)
      logical laccpt
      
 1082 format(a2,3f15.5)
 1083 format(i5,50f15.5)
      tempi0=1/kb

      inems=0
      nems=0
      do while (nems.lt.nemstot(imm))
        ibrea = 0
        do while (ibrea.lt.nhtot+1) 
          inems = inems+1
          write(6,*)
          write(6,*)"TRIAL GEOMETRY ",inems
c     read random coordinates
          if (sampwell(imm).eq.2) then
            call ranwell(rnwell,relwell1(imm))
            write(6,*)"Randomly reading initial states from well ",
     &       rnwell
            if (rnwell.eq.1) then
              call ransamp(lbinsamp(imm),lems(imm),samptot(imm),
     &         sampfilexx(imm),sampfilepp(imm),
     &         natom(imm),xxmt,ppmt,vmt)
            elseif (rnwell.eq.2) then
              call ransamp(lbinsamp2(imm),lems(imm),samptot2(imm),
     &         sampfilexx2(imm),sampfilepp2(imm),
     &         natom(imm),xxmt,ppmt,vmt)
            endif
          else
            write(6,*)"Randomly reading initial states"
            call ransamp(lbinsamp(imm),lems(imm),samptot(imm),
     &       sampfilexx(imm),sampfilepp(imm),
     &       natom(imm),xxmt,ppmt,vmt)
          endif
c     place CoM at origin
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
c     calculate PE (read, do not call pot again, expensive for dd)
          write(6,*)"Trial PE ",vmt," au"
          if (vmt.gt.emicr(imm)) then
             write(6,*)"Potential energy > Etot"
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
c     remove angular momentum
        call noang(xxmi,ppmi,mmmi,natom(imm))
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
c     output coordinates and momenta
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
