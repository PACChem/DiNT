      program dint
c     DiNT:  Adiabatix and Nonadiabatic Trajectories
c     Main program.  This routine calls HEADER and READIN once, then
c     calls PREMOL once per AG to precalculate some constants.
c     DRIVER is called NTRAJ times, once for each trajectory.
c     FINALSTATE is called after each trajectory to write some data.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
      include 'c_traj.f'
#include <sprng_f.h>

      integer im,initrng,minitrng,i,ia
      double precision t1,t2,t3,cputime,t3old,r

      call header
 
      write(6,'("...entering READIN...",/)')
      call readin
      write(6,'("...exiting READIN...",/)')
 
 
 
      call timing(t1)
      write(6,'("...entering PREMOL...",/)')
c     prepare each molecule
c     compute quantities that are the same for all trajectories
c     initialize electronic coordinates
      call initelec
      do im = 1 , nmol
      call premol(im)
      enddo
      write(6,'("...exiting PREMOL...",/)')
      call timing(t2)
      write(6,'(" CPU time in PREMOL is ",f10.5," s",/)')(t2-t1)
    
c     redirect for N(E)
      if (methflag.eq.10) then
      initrng=1
      rng_stream = init_sprng(initrng,maxtraj,ranseed,SPRNG_DEFAULT)
        im=1
        call mc(im)
        go to 999
      endif 
      
      t3 = 0.d0
      cputime = 0.d0
c     run only on head process

c     loop over trajectories
      nistep=0
      nistepw=0

      do itraj=1,ntraj

      write(6,*)"***TRAJECTORY ",itraj," OUT OF ",ntraj,"***"
      write(6,*)

      if (tflag(2).eq.1) then
          write(6,*)"Restarting trajectory #",trajlist(itraj)
          write(6,*)
          initrng = trajlist(itraj)-1
      else
          initrng = itraj-1
      endif

c     initialize this trajectory
      nsurf = nsurf0
      hstep = hstep0
c     initialize electronic coordinates
      call initelec
c     initialize the randum number generator
      rng_stream = init_sprng(initrng,maxtraj,ranseed,SPRNG_DEFAULT)

      write(6,'("...entering INITMOL...",/)')
      do im = 1 , nmol
c     generate specific initial coordinates for each molecule
      call initmol(im)
      enddo
      write(6,'("...exiting INITMOL...",/)')

      write(6,'("...entering ORIENT...",/)')
      do im = 1 , nmol
c     orient AGs with respect to each other
      call orient(im)
      enddo
      write(6,'("...exiting ORIENT...",/)')

      write(6,'("...entering PREPTRAJ...",/)')
      call preptraj
      write(6,'("...exiting PREPTRAJ...",/)')

c     propagate this trajectory
      write(6,'("...entering DRIVER...",/)')
      call driver
      write(6,'("...exiting DRIVER...",/)')
      t3old = t3
      call timing(t3)
      cputime = cputime + (t3-t3old)
      write(6,'(" CPU time for this trajectory is ",f10.5," s",/)')
     & (t3-t3old)

c     analyze trajectory
      write(6,'("...entering FINALSTATE...",/)')
      call finalstate
      write(6,'("...exiting FINALSTATE...",/)')

c     end loop over trajectories
      enddo

c     from ransamp.F

      write(6,'(" CPU time for all trajectories is ",f10.5," s or ",
     & f10.5," s per trajectory.",/)')
     & cputime,cputime/dble(ntraj)

 999  write(6,'("Well done!")')
      end
