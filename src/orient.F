      subroutine orient(im)

c ORIENT AGs WITH RESPECT TO EACH OTHER

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>
 
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat),murel,prel,
     & bimpact,ran,pref,tmp,tempqc2
      integer im,i,j,k,ii

      write(6,*)"Orienting AG ",im
      write(6,*)"------------------"

c     store coords for this molecule in their own array for easier manipulation
      do i=1,natom(im)
       ii=i+iatom(im)
       mmm(i)=mm(ii)
       do j=1,3
        xxm(j,i) = xx(j,ii)
        ppm(j,i) = pp(j,ii)
       enddo
      enddo

C ************* IORIENT = 0 ****************
      if (iorient.eq.0) then
C ******************************************
c     assign CoM coorinates and momenta to fixed values which were read by READIN
      write(6,100)"Setting initial CoM momenta for this AG     to ",
     &             (compp(i,im),i=1,3)," au"
 100  format(1x,a,3f12.5,1x,a)
      write(6,100)"Setting initial CoM coordinates for this AG to ",
     &             (comxx(i,im)*autoang,i=1,3)," A"
      write(6,*)
      do i=1,natom(im)
        do k=1,3
        ppm(k,i)=ppm(k,i)+compp(k,im)*mmm(i)/mmag(im)
        xxm(k,i)=xxm(k,i)+comxx(k,im) 
        enddo
      enddo

C ************* IORIENT = 1 ****************
      elseif (iorient.eq.1) then
C ******************************************
        if (im.eq.1) then
          if (tempqc.gt.0.d0) then
c         random thermal dist for erelqc
          write(6,*)"Sampling from a thermal distribution at T = ",
     &    tempqc," K"
          call rancol(tempqc,erelqc)         
          tempqc2=tempqc
c          erelqc=erelqc+57.d0/627.509d0
          elseif (tempqc.eq.0.d0) then
c         use random sampled temp from AG 1
          tempqc2=sampjtemp(1)
          write(6,*)"Sampling from a thermal distribution at T = ",
     &    tempqc2," K"
          call rancol(tempqc2,erelqc)         
          else
c         constant relative energy
          erelqc=dabs(tempqc)
c          call no2arbeam(erelqc)
          endif
c         compute relative reduced mass
          murel=mmag(1)*mmag(2)/(mmag(1)+mmag(2))
          prel=dsqrt(2.d0*murel*erelqc)
          comxx(1,1)=0.d0
          comxx(2,1)=0.d0
          comxx(3,1)=0.d0
          compp(1,1)=0.d0
          compp(2,1)=0.d0
          compp(3,1)=-prel
      write(6,100)"Setting initial CoM momenta for this AG     to ",
     &             (compp(i,im),i=1,3)," au"
      write(6,100)"Setting initial CoM coordinates for this AG to ",
     &             (comxx(i,im)*autoang,i=1,3)," A"
      write(6,*)
      do i=1,natom(im)
        do k=1,3
        ppm(k,i)=ppm(k,i)+compp(k,im)*mmm(i)/mmag(im)
        xxm(k,i)=xxm(k,i)+comxx(k,im)
        enddo
      enddo
        elseif (im.eq.2) then
c         compute relative reduced mass
          murel=mmag(1)*mmag(2)/(mmag(1)+mmag(2))
          prel=dsqrt(2.d0*murel*erelqc)
          ran = sprng(rng_stream)

          if (bmaxqc.lt.0.d0) then
c select b**2 evenly from bmin to bmax (this is the correct distribution if bmin=0...)
          bimpact = bminqc**2+(bmaxqc**2-bminqc**2)*ran
          bimpact = dsqrt(max(bimpact,0.d0))
          else
c select b evenly from bmin to bmax (this is not the correct distribution, 
c may be used for importance sampling)
          bimpact = bminqc+(bmaxqc-bminqc)*ran
          endif

          bqci = bimpact
          lqci = bqci*prel
c          eorbqci = lqci**2/(2.d0*murel*(rel0qc**2+bqci**2))
          eorbqci = lqci**2/(2.d0*murel*(rel0qc**2))
          erelqci = erelqc - eorbqci
          write(6,*)"Impact parameter = ",bimpact*autoang," A"
          write(6,*)"Orbital angular momentum = ",bimpact*prel," au"
          write(6,*)"Orbital energy  = ",eorbqci*autoev," eV"
          write(6,*)"Translational energy = ",erelqci*autoev," eV"
          write(6,*)"Total relative energy = ",erelqc*autoev," eV"
          if (tempqc.gt.0d0) then
          write(6,*)"Total relative temp = ",(2.d0/3.d0)*erelqc/kb,
     &      " K (target is ",tempqc2," K)"
          pref = dsqrt(8.d0*kb*tempqc2/(pi*murel))*pi*bmaxqc**2   ! au
          pref = pref*(autoang)**3/autofs                        ! Ang^3/fs
          pref = pref*(1.d-8**3)*1.d15                           ! cm^3/s
          write(6,*)"k(T) prefactor = ",pref," cm3/s"
          endif
          write(6,*)
          comxx(1,2)=0.d0
          comxx(2,2)=bimpact
          if (bimpact.gt.rel0qc) then
            write(6,*)"Impact parameter (",bimpact*autoang," A) greater"
     &      ," than initial separation (",rel0qc*autoang," A)" 
            stop
          endif
          tmp = dsqrt(rel0qc**2-bimpact**2)
          comxx(3,2)=-tmp
          compp(1,2)=0.d0
          compp(2,2)=0.d0
          compp(3,2)=prel
      write(6,100)"Setting initial CoM momenta for this AG     to ",
     &             (compp(i,im),i=1,3)," au"
      write(6,100)"Setting initial CoM coordinates for this AG to ",
     &             (comxx(i,im)*autoang,i=1,3)," A"
      write(6,*)
      do i=1,natom(im)
        do k=1,3
        ppm(k,i)=ppm(k,i)+compp(k,im)*mmm(i)/mmag(im)
        xxm(k,i)=xxm(k,i)+comxx(k,im)
        enddo
      enddo

        else
          write(6,*)"IM = ",im," in ORIENT"
          stop
        endif
C ******************************************
      else
C ******************************************
        write(6,*)"IORIENT = ",iorient," in ORIENT"
        stop
C ******************************************
      endif
C ******************************************

c write manipulated data back to xx and pp arrays
      do i=1,natom(im)
       ii=i+iatom(im)
       do j=1,3
        xx(j,ii) = xxm(j,i)
        pp(j,ii) = ppm(j,i)
       enddo
      enddo

      return

      end
