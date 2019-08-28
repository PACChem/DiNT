      subroutine vibwells(xxmt,ppmt,sceng,imm)

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
#include <sprng_f.h>

      integer i,j,ii,imm,rnwell,nfq
      double precision xxmt(3,mnat),ppmt(3,mnat),evb,sceng,
     & r1(mnat*3),tmp

      if (temp0im(imm).gt.0.d0) then
c     random thermal
      if (sampwell(imm).eq.2) then
        call ranwell(rnwell,relwell1(imm))
        write(6,*)"Randomly assigning initial states from well ",
     &   rnwell," for each normal mode for a target temp ",
     &   temp0im(imm)
      else
        rnwell=1
        write(6,*)"Picking initial vibrational states",
     &   " for ",temp0im(imm)," K"
      endif
      write(6,*)"Mode #   Quanta"
      nfq=3*natom(imm)-6
      if (natom(imm).eq.2) nfq=1  ! diatom
      if (rnwell.eq.1) then
c 201    evb=0.d0
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx0(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
        do i=1,nfq
          call ranstate(freq(i,imm),nmqn(i,imm),temp0im(imm))
          write(6,*)i,nmqn(i,imm)
c          evb=evb+freq(i,imm)*(0.5d0+nmqn(i,imm))
        enddo
        write(6,*)
      else if (rnwell.eq.2) then
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx02(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
c 202    evb=emin2(imm)
        do i=1,nfq
          call ranstate(freq2(i,imm),nmqn(i,imm),temp0im(imm))
          write(6,*)i,nmqn(i,imm)
c          evb=evb+freq2(i,imm)*(0.5d0+nmqn(i,imm))
        enddo
      endif

      else
c     random microcanonical
      if (sampwell(imm).eq.2) then 
        call ranwell(rnwell,relwell1(imm))
      else
         rnwell=1
      endif
      nfq=3*natom(imm)-6
      if (rnwell.eq.1) then
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx0(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
        tmp = 0.d0
        do i=1,nfq
          r1(i) = sprng(rng_stream)
          tmp = tmp + r1(i)
        enddo
        write(6,*)"Mode QN and/or energy in unbound mode (au)"
        do i=1,nfq
          r1(i)=r1(i)/tmp*sceng   ! randomly distrubuted harmonic energy in mode i
          if (nmtype(imm).eq.1.and.i.eq.nfq) then
          nmqn(i,imm)=r1(i) ! Energy in unbound mode
          else
          nmqn(i,imm)=r1(i)/freq(i,imm)-0.5d0 ! effective quantum number for that much HO energy
          endif
          write(6,*)i,nmqn(i,imm)
        enddo
        write(6,*)"For a total harmonic energy of ",
     &  sceng*autoev," eV"
      else if (rnwell.eq.2) then 
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx02(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
        tmp = 0.d0
        do i=1,nfq
          r1(i) = sprng(rng_stream)
          tmp = tmp + r1(i)
        enddo
        write(6,*)"Mode QN (PE well 2)"
        do i=1,nfq
          r1(i)=r1(i)/tmp*(sceng-emin2(imm))   ! randomly distrubuted harmonic energy in mode i
          nmqn(i,imm)=r1(i)/freq2(i,imm)-0.5d0 ! effective quantum number for that much HO energy
          write(6,*)i,nmqn(i,imm)
        enddo
        write(6,*)"For a total harmonic energy of ",
     &  (sceng-emin2(imm))*autoev," eV"
      endif

      endif
      return
          
      end
