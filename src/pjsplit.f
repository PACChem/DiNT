      subroutine pjsplit(xx,pp,mm,ppvib)

c Removes overall angular motion from PP and places the resulting
c momenta in PPVIB.  For three atoms only.
c From NAT8.1.  This subroutine calculates
c the internal/J coordinates (Jint) and projects them into 
c the vibrational subspace in the jacobi rep (pvib) and the
c rotational subspace in the jacobi rep (prot)
c All of the angular momentum for the p vector will be
c in the last three coords of jint, and in the vector prot

c xk is the angular momentum of the system, evib is the energy
c in internal motion, erot is the energy in rotational motion

      implicit none
      include 'param.f'

      double precision xx(3,mnat),pp(3,mnat),mm(mnat)
      double precision p(9),q(9),jm(9)
      double precision pr(6),r(6),np(3),nq(6),Jint(6)
      double precision pvib(6),prot(6),xdum(3,3),ppvib(3,3)

      double precision nj(6),nk(6),njnorm,nknorm,rr(3)

c q and p are the original jacobi coordinates, jm is the reduced mass

      double precision xk,qnorm,Qnorm2, npnorm
      double precision nqnorm, nqnorm2,erot,evib
      double precision coschi
      integer i,arr

c     get Jacobis from cartesians
c     get internals
      call rtox(xx,rr,1)

c     pick which jacobis to use
      if (rr(1).lt.rr(2).and.rr(1).lt.rr(3)) arr=1
      if (rr(2).lt.rr(3).and.rr(2).lt.rr(3)) arr=2
      if (rr(3).lt.rr(1).and.rr(3).lt.rr(2)) arr=3
      arr = 2
      call carttojac(xx,pp,mm,q,p,jm,0,arr)

c  transform q's and p's to mass weighted coordinates 
      do i = 1,6 
        p(i) = p(i)/dsqrt(jm(i))
        q(i) = q(i)*dsqrt(jm(i))
      enddo

      qnorm =  dsqrt(q(1)**2 + q(2)**2 +q(3)**2)
      Qnorm2 = dsqrt(Q(4)**2 + Q(5)**2 +Q(6)**2)
c calculate angle between q and Q . . .
      coschi = (q(1)*Q(4)+ q(2)*Q(5)+q(3)*Q(6))/(qnorm*Qnorm2)

      if (dabs(coschi) .ne. 1.0d0) then

c now we calculate vectors we need for the transformation
c np is a vector normal to the Q,q plane
      np(1) = (q(2)*Q(6) - q(3)*Q(5))
      np(2) = (q(3)*Q(4) - q(1)*Q(6))
      np(3) = (q(1)*Q(5) - q(2)*Q(4))
      npnorm = dsqrt(np(1)**2+np(2)**2+np(3)**2)

      do i =1,3
        np(i) = np(i)/npnorm
      enddo

c nq(1-3) is a vector normal to q and np
      nq(1) = (np(2)*q(3)-np(3)*q(2))
      nq(2) = (np(3)*q(1)-np(1)*q(3))
      nq(3) = (np(1)*q(2)-np(2)*q(1))
      nqnorm = dsqrt(nq(1)**2+nq(2)**2+nq(3)**2) 
      do i =1,3
        nq(i) = nq(i)/nqnorm
      enddo

c nQ(4-6) is a vector normal to Q and np
      nQ(4) = (np(2)*Q(6)-np(3)*Q(5))
      nQ(5) = (np(3)*Q(4)-np(1)*Q(6))
      nQ(6) = (np(1)*Q(5)-np(2)*Q(4))
      nqnorm2 = dsqrt(nq(4)**2+nq(5)**2+nq(6)**2) 
      do i =4,6
        nQ(i) = nQ(i)/nqnorm2
      enddo

c we project p and P onto the 6 new vectors (q,nq,np,Q,nQ,np)

      pr(1) = (p(1)*q(1)+p(2)*q(2)+p(3)*q(3))
     &       / qnorm
      pr(2) = (p(1)*nq(1) + p(2)*nq(2)+ p(3)*nq(3))
      pr(3) = (p(1)*np(1) + p(2)*np(2) + p(3)*np(3))
      pr(4) = (P(4)*Q(4)+P(5)*Q(5)+P(6)*Q(6))
     &      / Qnorm2
      pr(5) = (P(4)*nQ(4) + P(5)*nQ(5)+ P(6)*nQ(6))
      pr(6) = (P(4)*np(1) + P(5)*np(2) + P(6)*np(3))

c now lets transform into internal/J coords . . 
c pr(1) and pr(4) are projections on r and R,
c pr(3) and pr(6) are out of plane motion.
c pr(2) and pr(5) are in plane motion . . .

c coords 1-3 have no J components
      Jint(1) = pr(1)
      Jint(2) = pr(4)
      Jint(3) = (Qnorm2*pr(2) - qnorm*pr(5))
     &        /dsqrt(qnorm**2 + Qnorm2**2)

c coords 4-6 have no internal components 
      Jint(4) = (qnorm*pr(2) + Qnorm2*pr(5))
     &        /dsqrt(qnorm**2 + Qnorm2**2)
      Jint(5) = pr(3)
      Jint(6) = pr(6)

c calculatate J and energies (here's how)

      evib = 0.0d0
      erot = 0.0d0
      do i=1,3
        evib = evib + 0.5d0*Jint(i)**2 
        erot = erot + 0.5d0*Jint(i+3)**2 
      enddo

      xk = 0.0d0
      do i =1,3
        xk = xk 
     &    + (np(i)*Jint(4)*dsqrt(qnorm**2+Qnorm2**2)
     &    - nq(i)*qnorm*Jint(5)
     &    - nQ(i+3)*Qnorm2*Jint(6))**2
      enddo
      xk = dsqrt(xk)


c now transform back into prot and pvib

c 1st set Jint(1-3) to zero and transform 
c to get prot components . . .
      
      pr(1) = 0.0d0
      pr(2) = (qnorm*Jint(4))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = Jint(5)
      pr(4) = 0.0d0
      pr(5) = (Qnorm2*Jint(4))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = Jint(6)

      do i=1,3
        prot(i) = pr(1)*q(i)/qnorm + pr(2)*nq(i) + pr(3)*np(i)
      end do

      do i=4,6
        prot(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nq(i) + pr(6)*np(i-3)
      end do


c 2nd set Jint(4-6) to zero and transform 
c to get pvib components . . .

      pr(1) = Jint(1)
      pr(2) = (Qnorm2*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = 0.0d0
      pr(4) = Jint(2)
      pr(5) = (-qnorm*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = 0.0d0

      do i=1,3
        pvib(i) = pr(1)*q(i)/qnorm + pr(2)*nq(i) + pr(3)*np(i)
      end do

      do i=4,6
        pvib(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nq(i) + pr(6)*np(i-3)
      end do


c 3rd transform all Jint's back . . . should obtain original p's
      pr(1) = Jint(1)
      pr(2) = (qnorm*Jint(4)+Qnorm2*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = Jint(5)
      pr(4) = Jint(2)
      pr(5) = (Qnorm2*Jint(4)-qnorm*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = Jint(6)

      do i=1,3
        p(i) = pr(1)*q(i)/qnorm + pr(2)*nq(i) + pr(3)*np(i)
      end do

      do i=4,6
        p(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nq(i) + pr(6)*np(i-3)
      end do

c finally transform back to non mass scaled coordinates
      do i = 1,6 
        q(i) = q(i)/dsqrt(jm(i))
        p(i) = p(i)*dsqrt(jm(i))
        prot(i) = prot(i)*dsqrt(jm(i))
        pvib(i) = pvib(i)*dsqrt(jm(i))
      enddo

      else

c here the vectors are collinear, which means
c a different transformation procedure . . .

      if ((q(2) .eq. 0.0d0) .and. (q(3) .eq. 0.0d0)) then
c q points along 1,0,0 . . . so 
        nj(1) = 0.0d0
        nj(2) = 1.0d0
        nj(3) = 0.0d0

        nk(1) = 0.0d0
        nk(2) = 0.0d0
        nk(3) = 1.0d0

      else
c q does not point along 1,0,0.  So calculate nk = 1,0,0 x q
        nj(1) = 0.0d0
        nj(2) = -q(3)
        nj(3) = q(2)
        njnorm = dsqrt(q(2)**2+q(3)**2)
        do i =1,3
          nj(i) = nj(i)/njnorm
        enddo

c now calculate a vector nj . . .

        nk(1) = (q(2)*nj(3) - q(3)*nj(2))
        nk(2) = (q(3)*nj(1) - q(1)*nj(3))
        nk(3) = (q(1)*nj(2) - q(2)*nj(1))
    
        nknorm = dsqrt(nk(1)**2+nk(2)**2+nk(3)**2)
        do i =1,3
          nk(i) = nk(i)/nknorm
        enddo
      endif

c we project p and P onto the 3 new vectors (q,nj,nk,Q,nj,nk)

      pr(1) = (p(1)*q(1)+p(2)*q(2)+p(3)*q(3))
     &       / qnorm
      pr(2) = (p(1)*nj(1) + p(2)*nj(2)+ p(3)*nj(3))
      pr(3) = (p(1)*nk(1) + p(2)*nk(2) + p(3)*nk(3))
      pr(4) = (P(4)*Q(4)+P(5)*Q(5)+P(6)*Q(6))
     &      / Qnorm2
      pr(5) = (P(4)*nj(1) + P(5)*nj(2)+ P(6)*nj(3))
      pr(6) = (P(4)*nk(1) + P(5)*nk(2) + P(6)*nk(3))

c now lets transform into internal/J coords . .

c coords 1-4 have no J components
      Jint(1) = pr(1)
      Jint(2) = pr(4)
      Jint(3) = (Qnorm2*pr(2) - qnorm*coschi*pr(5))
     &        /dsqrt(qnorm**2 + Qnorm2**2)
      Jint(4) = (Qnorm2*pr(3) - qnorm*coschi*pr(6))
     &        /dsqrt(qnorm**2 + Qnorm2**2)

c coords 5-6 have no internal components
      Jint(5) = (qnorm*coschi*pr(2) + Qnorm2*pr(5))
     &        /dsqrt(qnorm**2 + Qnorm2**2)
      Jint(6) = (qnorm*coschi*pr(3) + Qnorm2*pr(6))
     &        /dsqrt(qnorm**2 + Qnorm2**2)

c calculatate J and energies (here's how)

      evib = 0.0d0
      erot = 0.0d0
      do i=1,4
        evib = evib + 0.5d0*Jint(i)**2
      enddo
      do i=5,6
        erot = erot + 0.5d0*Jint(i)**2
      enddo

      xk = (Jint(5)**2+Jint(6)**2)*(qnorm**2+Qnorm2**2)
      xk = dsqrt(xk)

c now transform back to prot and pvib

c 1st set Jint(1-4) to zero and transform
c to get prot components . . .

      pr(1) = 0.0d0
      pr(2) = (qnorm*Jint(5)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = (qnorm*Jint(6)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(4) = 0.0d0

      pr(5) = (Qnorm2*Jint(5))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = (Qnorm2*Jint(6))
     &               /dsqrt(qnorm**2+Qnorm2**2)

      do i=1,3
        prot(i) = pr(1)*q(i)/qnorm + pr(2)*nj(i) + pr(3)*nk(i)
      end do

      do i=4,6
        prot(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nj(i-3) + pr(6)*nk(i-3)
      end do

c 2nd set Jint(4-6) to zero and transform
c to get pvib components . . .

      pr(1) = Jint(1)
      pr(2) = (Qnorm2*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = (Qnorm2*Jint(4))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(4) = Jint(2)

      pr(5) = (-qnorm*Jint(3)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = (-qnorm*Jint(4)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)

      do i=1,3
        pvib(i) = pr(1)*q(i)/qnorm + pr(2)*nj(i) + pr(3)*nk(i)
      end do

      do i=4,6
        pvib(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nj(i-3) + pr(6)*nk(i-3)
      end do


c 3rd transform all Jint's back . . . should obtain original p's
      pr(1) = Jint(1)
      pr(2) = (qnorm*Jint(5)*coschi+Qnorm2*Jint(3))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(3) = (qnorm*Jint(6)*coschi+Qnorm2*Jint(4))
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(4) = Jint(2)
 
      pr(5) = (Qnorm2*Jint(5)-qnorm*Jint(3)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)
      pr(6) = (Qnorm2*Jint(6)-qnorm*Jint(4)*coschi)
     &               /dsqrt(qnorm**2+Qnorm2**2)

      do i=1,3
        p(i) = pr(1)*q(i)/qnorm + pr(2)*nj(i) + pr(3)*nk(i)
      end do

      do i=4,6
        p(i) = pr(4)*q(i)/Qnorm2 + pr(5)*nj(i-3) + pr(6)*nk(i-3)
      end do

c finally transform back to non mass scaled coordinates
      do i = 1,6
        q(i) = q(i)/dsqrt(jm(i))
        p(i) = p(i)*dsqrt(jm(i))
        prot(i) = prot(i)*dsqrt(jm(i))
        pvib(i) = pvib(i)*dsqrt(jm(i))
      enddo

      endif

c     get cartesians from Jacobis
      call carttojac(xdum,ppvib,mm,q,pvib,jm,1,arr)

      end


