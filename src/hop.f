      subroutine hop(ithistraj,istep,time,nclu,nsurf,newsurf,
     &    pem,gpem,dvec,xx,pp,tau)

c Perform a successful or a frustrated hop from NSURF to NEWSURF.
c For a successful hop, NSURF is set to NEWSURF and PP is adjusted
c   to conserve energy.
c For a frustrated hop, NSURF is unchanged, and PP is adjusted for
c   a reflected-and-hop-back.  PP is unchanged for an ingored
c   frustrated hop.

      implicit none
      include 'param.f'
      include 'c_sys.f'

c input
      integer nclu,nsurf,newsurf,ithistraj,istep
      double precision time,xx(3,mnat),tau

c input/output
      double precision pp(3,mnat)

c local
      integer frusflag
      integer i,j
      double precision dxi
      double precision pem(mnsurf),gpem(3,mnat,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),r(3),
     & sectrm,a,b,a2,ab,ra,radi,scr,tmp,deltae,dot1,dot2

c new
      double precision radi2,c,sectrm2,ppvib(3,mnat),scr2,deltae2

c compute distance matrix for output
      if (nclu.eq.3) call rtox(xx,r,1) 

c check for forbidden hops
        deltae = pem(newsurf) - pem(nsurf)

        do i=1,3
        do j=1,nclu
          ppvib(i,j)=pp(i,j)
        enddo
        enddo
        call noang(xx,ppvib,mm,nclu)

        a=0.0d0
        b=0.0d0
        c=0.0d0
        do i=1,3
        do j=1,nclu
          scr=dvec(i,j,nsurf,newsurf)/mm(j)
c         a = momentum.dot.dvec/mass
c         b = dvec.dot.dvec/mass
c         c = pvib.dot.pvib/mass
          a=a+scr*pp(i,j)
          b=b+scr*dvec(i,j,nsurf,newsurf)
          c=c+ppvib(i,j)*ppvib(i,j)/mm(j)
        enddo
        enddo
        a2=a*a
        ab=a/b
c       sectrm = deltae/kinetic energy along d
        sectrm = 2.0d0*deltae*b/a2
        deltae2=deltae-a2/(2.d0*b)
        sectrm2 = deltae2/(c/2.d0-a2/(2.d0*b))
        radi = 1.0d0-sectrm
        radi2 = 1.0d0-sectrm2
c       if radi > 0 then delta > kinetic energy along d, so we can hop
c       note that for a hop down, sectrm is negative, so we can always hop
        if(radi.ge.0.0d0) then
c         successful hop
          write(6,*)"HOP!",nsurf," ---> ",newsurf
          write(33,1033)ithistraj,istep,time*autofs,"H",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
 1033     format(i5,i10,f12.3,1x,a,2i5,5f12.3)
          radi=sqrt(radi)
          ra=ab*(1.0d0-radi)
          do i=1,3
          do j=1,nclu
           pp(i,j) = pp(i,j) - dvec(i,j,nsurf,newsurf)*ra
          enddo
          enddo
          nsurf = newsurf



c NEW NEW NEW
c currently set up to skip this part
c        elseif(radi2.ge.0.0d0) then
        elseif(.false.) then
c modified successful hop
          write(6,*)"NEW HOP!",nsurf," ---> ",newsurf
          write(33,1033)ithistraj,istep,time*autofs,"V",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
c         take out all the KE along d first
          do i=1,3
          do j=1,nclu
           pp(i,j) = pp(i,j) - dvec(i,j,nsurf,newsurf)*ab
          enddo
          enddo
c         now take out energy along ppvib
          do i=1,3
          do j=1,nclu
            ppvib(i,j)=pp(i,j)
          enddo
          enddo
          call noang(xx,ppvib,mm,nclu)
        a=0.0d0
        b=0.0d0
        do i=1,3
        do j=1,nclu
          scr=ppvib(i,j)/mm(j)
c         a = momentum.dot.dvec/mass
c         b = dvec.dot.dvec/mass
          a=a+scr*pp(i,j)
          b=b+scr*ppvib(i,j)
        enddo
        enddo
        a2=a*a
        ab=a/b
c       sectrm = deltae/kinetic energy along d
        sectrm = 2.0d0*deltae2*b/a2
        radi = 1.0d0-sectrm
          radi=sqrt(radi)
          ra=ab*(1.0d0-radi)
          do i=1,3
          do j=1,nclu
           pp(i,j) = pp(i,j) - ppvib(i,j)*ra
          enddo
          enddo

          nsurf = newsurf
c END NEW NEW NEW





        else
c         frustrated hop
c NOTE:   FRUSFLAG IS HARDCODED
c         FRUSFLAG = 0 --->  IGNORE ALL FRUSTRATED HOPS (I.E., "+")
c         FRUSFLAG = 1 --->  REFLECT ALL FRUSTRATED HOPS (I.E., "-")
c         FRUSFLAG = 2 --->  GRADV FOR ALL FRUSTRATED HOPS (I.E., "gV")
          if (methflag.eq.1) then
c           for FS we use "+"
            frusflag = 0
          else if (methflag.eq.5) then
c           for FSTU we use "gradV"
            frusflag = 2
          else
            print *,"METHFLAG != 1 or 5 in HOP"
            stop
          endif
          if (frusflag.eq.0) then
c           do nothing
            write(6,*)"FRUST! IGNORE"
          write(33,1033)ithistraj,istep,time*autofs,"I",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
          elseif (frusflag.eq.1) then
c           frustrated hop, reflect
            write(6,*)"FRUST! REFLECT"
          write(33,1033)ithistraj,istep,time*autofs,"R",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
            ra=2.0d0*ab
            do i=1,3
            do j=1,nclu
              pp(i,j) = pp(i,j) - dvec(i,j,nsurf,newsurf)*ra
            enddo
            enddo
          elseif (frusflag.eq.2) then
c           gradV method
c           dot1 is the gradient of the NEWSURF in the direction of dvec
            dot1 = 0.d0
            do i=1,3
            do j=1,nclu
              dot1 = dot1 + gpem(i,j,newsurf)*dvec(i,j,nsurf,newsurf)
            enddo
            enddo
c           dot2 is the momentum in the direction of dvec
            dot2 = 0.d0
            do i=1,3
            do j=1,nclu
              dot2 = dot2 + pp(i,j)*dvec(i,j,nsurf,newsurf)
            enddo
            enddo
            if (dot1*dot2.gt.0.d0) then
c           reflect
            write(6,*)"FRUST! REFLECT"
          write(33,1033)ithistraj,istep,time*autofs,"R",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
            ra=2.0d0*ab
            do i=1,3
            do j=1,nclu
              pp(i,j) = pp(i,j) - dvec(i,j,nsurf,newsurf)*ra
            enddo
            enddo
            else
            write(6,*)"FRUST! IGNORE"
          write(33,1033)ithistraj,istep,time*autofs,"I",nsurf,
     &      newsurf,pem(nsurf)*autoev,pem(newsurf)*autoev,tau*autofs
c           ignore
            endif
          else
            write(6,*)"FRUSFLAG = ",frusflag," in HOPCHECK"
            stop
          endif
        endif

      return
      end
