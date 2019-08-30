c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between Al and Al.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine al_al_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,i
      double precision zetas,zetap,zetasp
      double precision sss,dsss,sps,dsps,pps,dpps,ppp,dppp,pss,dpss
      double precision r,kwh
      double precision x,kss,ksp,kpp
      double precision tmp1,tmp2,tmp3,tmp4,temp,intgr,dintgr

      common /coulom/ vsip
      common /distan/ r
c
c STO parameters for Al
c
      DATA zetas/1.3724/ ! Al s exponent in overlap parameter in au
      DATA zetap/1.3552/ ! Al p exponent in overlap parameter in au
      zetasp = (zetas+zetap)/2.d0
c
c Wolfberg-Holmholz constant for Al
c
      DATA kwh /0.39280/
c
c Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(9,1)+vsip(9,1))/2.0
      tmp2=(vsip(9,1)+vsip(9,2))/2.0
      tmp3=(vsip(9,2)+vsip(9,2))/2.0

      kss = kwh*tmp1
      ksp = kwh*tmp2
      kpp = kwh*tmp3
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
      call ovlcon(3,0,zetas,3,0,0,zetas,intgr,dintgr)
      dsss=dintgr*kss
       sss= intgr*kss

c s p sigma bonding overlap integral in z-axis
      call ovlcon(3,0,zetasp,3,1,0,zetasp,intgr,dintgr)
      dsps=dintgr*ksp
       sps= intgr*ksp

c p p sigma bonding overlap integral in z-axis
      call ovlcon(3,1,zetap,3,1,0,zetap,intgr,dintgr)
      dpps=dintgr*kpp
       pps= intgr*kpp

c p p pi bonding overlap integral in z-axis
      call ovlcon(3,1,zetap,3,1,1,zetap,intgr,dintgr)
      dppp=dintgr*kpp
       ppp= intgr*kpp

c symmetry
      pss = -sps
      dpss = -dsps

c Project the two-center integral along the basis set.
c

      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)
c      print *,"AJ ham 11",t

      call sxyz('x',sps,dsps,t,dt)       !s-px
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)
c      print *,"AJ ham 12",t

      call sxyz('y',sps,dsps,t,dt)       !s-py
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)
c      print *,"AJ ham 13",t

      call sxyz('z',sps,dsps,t,dt)       !s-pz
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)
c      print *,"AJ ham 14",t

      call sxyz('x',pss,dpss,t,dt)       !px-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)
c      print *,"AJ ham 21",t

      call xx('x',pps,dpps,ppp,dppp,t,dt)   !px-px
      hamil(2,2)=t
      dhamil(2,2,1)=dt(1)
      dhamil(2,2,2)=dt(2)
      dhamil(2,2,3)=dt(3)
c      print *,"AJ ham 22",t

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !px-py
      hamil(2,3)=t
      dhamil(2,3,1)=dt(1)
      dhamil(2,3,2)=dt(2)
      dhamil(2,3,3)=dt(3)
c      print *,"AJ ham 23",t

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !px-pz
      hamil(2,4)=t
      dhamil(2,4,1)=dt(1)
      dhamil(2,4,2)=dt(2)
      dhamil(2,4,3)=dt(3)
c      print *,"AJ ham 24",t

      call sxyz('y',pss,dpss,t,dt)          !py-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)
c      print *,"AJ ham 31",t

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !py-px
      hamil(3,2)=t
      dhamil(3,2,1)=dt(1)
      dhamil(3,2,2)=dt(2)
      dhamil(3,2,3)=dt(3)
c      print *,"AJ ham 32",t

      call xx('y',pps,dpps,ppp,dppp,t,dt)   !py-py
      hamil(3,3)=t
      dhamil(3,3,1)=dt(1)
      dhamil(3,3,2)=dt(2)
      dhamil(3,3,3)=dt(3)
c      print *,"AJ ham 33",t

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !py-pz
      hamil(3,4)=t
      dhamil(3,4,1)=dt(1)
      dhamil(3,4,2)=dt(2)
      dhamil(3,4,3)=dt(3)
c      print *,"AJ ham 34",t

      call sxyz('z',pss,dpss,t,dt)          !pz-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)
c      print *,"AJ ham 41",t

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !pz-px
      hamil(4,2)=t
      dhamil(4,2,1)=dt(1)
      dhamil(4,2,2)=dt(2)
      dhamil(4,2,3)=dt(3)
c      print *,"AJ ham 42",t

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !pz-py
      hamil(4,3)=t
      dhamil(4,3,1)=dt(1)
      dhamil(4,3,2)=dt(2)
      dhamil(4,3,3)=dt(3)
c      print *,"AJ ham 43",t

      call xx('z',pps,dpps,ppp,dppp,t,dt)   !pz-pz
      hamil(4,4)=t
      dhamil(4,4,1)=dt(1)
      dhamil(4,4,2)=dt(2)
      dhamil(4,4,3)=dt(3)
c      print *,"AJ ham 44",t

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Al and Al.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine al_al_rp (repul,drepul)
c
      implicit none
      double precision repul,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision aa,bb,mu,r,dt
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /distan/ r
c
c Repulsive parameters for carbon-carbon
c
      DATA aa/1142.d0/
      DATA bb/2.921d0/
      DATA mu/0.02936d0/

      repul=aa/(r**mu)*dexp(-bb*r)
      dt=-bb*aa/(r**mu)*dexp(-bb*r)+aa/(r**(mu-1.d0))*dexp(-bb*r)*(-mu)
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between Al and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine al_h_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,i
      double precision zetas,zetap,zetasp
      double precision sss,dsss,sps,dsps,pps,dpps,ppp,dppp,pss,dpss
      double precision r,kwh
      double precision x,kss,ksp,kpp
      double precision tmp1,tmp2,tmp3,tmp4,temp,intgr,dintgr
      double precision szetapm

      common /coulom/ vsip
      common /distan/ r
c
c STP parameters for H
c
      DATA szetapm /1.2/
c
c STO parameters for Al
c
      DATA zetas/1.3724/ ! Al s exponent in overlap parameter in au
      DATA zetap/1.3552/ ! Al p exponent in overlap parameter in au
      zetasp = (zetas+zetap)/2.d0
c
c Wolfberg-Holmholz constant for Al&H
c
      DATA kwh /0.39280/
c
c Wolfsberg-Helmholtz relation for Al(3s)/Al(3p) and H(1s) VSIPs.
c
      tmp1=(vsip(9,1)+vsip(1,1))/2.0
      tmp2=(vsip(9,2)+vsip(1,1))/2.0

      kss = kwh*tmp1
      ksp = kwh*tmp2
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
      call ovlcon(1,0,szetapm,3,0,0,zetas,intgr,dintgr)
      dsss=dintgr*kss
       sss= intgr*kss

c s p sigma bonding overlap integral in z-axis
      call ovlcon(1,0,szetapm,3,1,0,zetasp,intgr,dintgr)
      dsps=dintgr*ksp
       sps= intgr*ksp

      dpss=-dsps
       pss=- sps

c Project the two-center integral along the basis set.
c


      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxyz('x',pss,dpss,t,dt)       !px-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)

      call sxyz('y',pss,dpss,t,dt)       !py-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)

      call sxyz('z',pss,dpss,t,dt)       !pz-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)

      return
      end
CC Note: The following section gives a fortran program to compute
CC coefficients by removing "C "(C and a space).

CC This program is made to obtain the C matrix  for a given l,N,L,M.
CC The formula for C matrix is given in reference:
CC H. W. Jones and C.A. Weatherford, Int. J. Quantum Chem. Symp. 12, 483 (1978).
CC In Page 486 they wrote: We build up a C matrix Cl(NA,NR) using eight nested 
CC DO loops to sort and collect the sums of all coefficients of all like powers
CC of r and a.
CC
CC This program only gives all C matrix coefficients, the modulus of 
CC N-L+2p+2q+2qprime-k-kprime to 2, the exponential index of a, and that of r.
CC After this program is compiled with f77, you can obtain a series of data
CC for C matrix.
C
C      program gen
C
C      implicit none
C      integer sl,N,L,M,p,q,v,pp,qp,vp,t,k,kp
C      integer result2,result3,result4,tmp
C      double precision result,alpha,beta,efunc
C       
C      print *, "Please input l, N, L, M"
C      read *,sl,N,L,M
C      if (sl .lt. M) stop "Irrational pair: l and M"
C      print 98,sl,N,L,M
C      do 1 p=0,(L+M)/2
C        do 2 q=0,(L+M-2*p)
C          do 3 v=0,(L+M-2*p-q)
C            do 4 pp=0,((sl-M)/2)
C              do 5 qp=0,sl-M-2*pp
C                do 6 vp=0,sl-M-2*pp-qp
C                  t=N-L+2*p+2*q+2*qp
C                  do 7 k=0,t
C                    do 8 kp=0,t-k
C                      result=(-1)**(L+v+qp)
C                      result=result/4.0**(L+sl-p-pp)
C                      tmp=L-M
C                      result=result*alpha(L,tmp,p)
C                      tmp=L+M-2*p
C                      result=result*beta(tmp,q)
C                      tmp=L+M-2*p-q
C                      result=result*beta(tmp,v)
C                      tmp=sl+M
C                      result=result*alpha(sl,tmp,pp)
C                      tmp=sl-M-2*pp
C                      result=result*beta(tmp,qp)
C                      tmp=sl-M-2*pp-qp
C                      result=result*beta(tmp,vp)
C                      tmp=N-L+2*p+2*q+2*qp
C                      result=result*efunc(tmp,k)
C                      tmp=N-L+2*p+2*q+2*qp-k
C                      result=result*beta(tmp,kp)
C                 
C                      result2=N-L+2*p+2*q+2*qp-k-kp
C                      result2=MOD(result2,2)
C
C                      result3=N+L+2*sl-2*pp-2*vp-2*v-k-kp
C
C                      result4=2*pp+2*v+2*vp+kp
C                     
C                      print 99, result,result2,result3,result4
C
C 8                  continue
C 7                continue
C 6              continue
C 5            continue
C 4          continue
C 3        continue
C 2      continue
C 1    continue
C
C 98   format(5x,"tttt ",4(I4,2x))
C 99   format(5x,f20.8,3(I4,2x))
C
C      end
C       
C
C      double precision function alpha(sl,n,p)
C      implicit none
C      integer sl,n,p,tmpa,tmpb,tmpc
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=2*sl-2*p
C      tmpb=2*sl-2*p-n
C      tmpc=sl-p
C
C      alpha=(-1)**p*fact(tmpa)/(fact(tmpb)*fact(tmpc)*fact(p))
C
C      return
C      end
C
C      double precision function beta(n,k)
C      implicit none
C      integer n,k,tmpa
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=n-k
C
C      beta=fact(n)/(fact(tmpa)*fact(k))
C
C      return
C      end
C
C      double precision function efunc(t,k)
C      implicit none
C      integer t,k,tmpa
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=t-k
C
C      efunc=fact(t)/fact(tmpa)
C
C      return
C      end
C
C      BLOCK DATA
C      double precision fact(0:14)
C      common /factorial/ fact
C      data fact/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
C     &          362880.0,3628800.0,39916800.0,479001600.0,622020800.0,
C     &          87178291200.0/
C      end
C
CC Suppose the executable file of this program is a.out and you are in a
CC UNIX system, the following commands let you sort and collect the sums of 
CC all coefficients of all like powers C of r and a.
CC     a.out | sort +2 -4 | awk -f awk.prog 
CC awk.prog is a file to do some special job for unix command awk.
CC File awk.prog begins.
CC /ttt/ { sl=$2;n=$3;l=$4;m=$5}
CC {pop[$3,$4] += $1}
CC END { printf("\nC Those are data for Lprime,N,L,M %d%d%d%d.\n\n",sl,n,l,m);
CC   for (i=0;i<=(sl+n+l);i++) 
CC     for (j=0;j<=(sl+n);j++)  
CC       printf("      DATA cmatrx(%2d,%2d,%2d,%2d,%2d,%2d)/%20.8f/\n",sl,n,l,
CC         m,i,j,pop[i,j]);
CC }
CC File awk.prog ends.
CC This is how the all above data are produced.

C C This is how the all following data are produced.

C Main data begins here.

       BLOCK DATA set1

      implicit none

      double precision cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)

      common /cmatrix/ cmatrx

C This is all the data block for the C matrix elements
C   l (lprime): 0-2
C   N         : 1-6
C   L         : 0-2
C   M         : 0-2
C If you want to calculate the overlap integral involving principal quantum
C number 7 or lagger, go ahead made your extra C matrix.
C I have some little program and notes followed those data.
C Advice: Do not trust me too much.

C Those are data for Lprime,N,L,M 0100.

      DATA cmatrx( 0, 1, 0, 0, 0, 0)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 0, 1)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 1, 0)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 1, 1)/          0.00000000/

C Those are data for Lprime,N,L,M 0200.

      DATA cmatrx( 0, 2, 0, 0, 0, 0)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 0, 1)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 0, 2)/          1.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 0)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 1)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 2)/          0.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 0)/          1.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 1)/          0.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 0210.

      DATA cmatrx( 0, 2, 1, 0, 0, 0)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 0, 1)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 0, 2)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 0)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 1)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 2)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 0)/         -2.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 1)/         -2.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 2)/          0.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 0)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 1)/          0.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 0300.

      DATA cmatrx( 0, 3, 0, 0, 0, 0)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 1)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 3)/          1.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 0)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 1)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 0)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 1)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 0)/          1.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0310.

      DATA cmatrx( 0, 3, 1, 0, 0, 0)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 1)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 2)/         -5.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 3)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 0)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 1)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 2)/         -5.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 3)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 0)/         -7.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 1)/         -7.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 2)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 0)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 1)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 0)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0320.

      DATA cmatrx( 0, 3, 2, 0, 0, 0)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 1)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 2)/         18.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 3)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 0)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 1)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 2)/         18.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 3)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 0)/         24.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 1)/         24.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 2)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 0)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 1)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 0)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 1)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0400.

      DATA cmatrx( 0, 4, 0, 0, 0, 0)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 1)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 2)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 4)/          1.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 0)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 1)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 2)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 0)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 1)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 2)/          6.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 0)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 1)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 0)/          1.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0410.

      DATA cmatrx( 0, 4, 1, 0, 0, 0)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 1)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 2)/        -27.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 3)/         -7.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 4)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 0)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 1)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 2)/        -27.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 3)/         -7.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 4)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 0)/        -33.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 1)/        -33.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 2)/        -15.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 3)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 0)/        -13.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 1)/        -13.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 2)/         -6.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 0)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 1)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 0)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0420.

      DATA cmatrx( 0, 4, 2, 0, 0, 0)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 1)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 2)/        117.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 3)/         27.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 4)/          3.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 0)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 1)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 2)/        117.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 3)/         27.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 4)/          3.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 0)/        141.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 1)/        141.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 2)/         60.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 3)/         13.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 0)/         51.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 1)/         51.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 2)/         21.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 0)/         15.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 1)/         15.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 2)/          6.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 0)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 1)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0500.

      DATA cmatrx( 0, 5, 0, 0, 0, 0)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 1)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 2)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 3)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 5)/          1.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 0)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 1)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 2)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 3)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 0)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 1)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 2)/         30.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 3)/         10.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 0)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 1)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 2)/         10.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 0)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 1)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0510.

      DATA cmatrx( 0, 5, 1, 0, 0, 0)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 1)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 2)/       -168.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 3)/        -48.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 5)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 0)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 1)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 2)/       -168.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 3)/        -48.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 5)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 0)/       -192.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 1)/       -192.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 2)/        -90.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 3)/        -26.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 4)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 0)/        -72.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 1)/        -72.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 2)/        -34.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 3)/        -10.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 0)/        -21.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 1)/        -21.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 2)/        -10.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 0)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 1)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0520.

      DATA cmatrx( 0, 5, 2, 0, 0, 0)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 1)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 2)/        855.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 3)/        225.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 4)/         36.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 5)/          3.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 0)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 1)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 2)/        855.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 3)/        225.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 4)/         36.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 5)/          3.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 0)/        975.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 1)/        975.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 2)/        438.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 3)/        113.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 4)/         17.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 0)/        345.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 1)/        345.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 2)/        153.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 3)/         38.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 0)/         96.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 1)/         96.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 2)/         42.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 3)/         10.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 0)/         23.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 1)/         23.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 2)/         10.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 0)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 1)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 0)/          1.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0600.

      DATA cmatrx( 0, 6, 0, 0, 0, 0)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 1)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 2)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 3)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 4)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 6)/          1.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 0)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 1)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 2)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 3)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 4)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 0)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 1)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 2)/        180.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 3)/         60.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 4)/         15.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 0)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 1)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 2)/         60.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 3)/         20.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 0)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 1)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 2)/         15.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 0)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 1)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 0610.

      DATA cmatrx( 0, 6, 1, 0, 0, 0)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 1)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 2)/      -1200.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 3)/       -360.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 4)/        -75.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 5)/        -11.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 6)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 0)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 1)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 2)/      -1200.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 3)/       -360.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 4)/        -75.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 5)/        -11.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 6)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 0)/      -1320.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 1)/      -1320.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 2)/       -630.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 3)/       -190.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 4)/        -40.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 5)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 0)/       -480.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 1)/       -480.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 2)/       -230.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 3)/        -70.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 4)/        -15.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 0)/       -135.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 1)/       -135.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 2)/        -65.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 3)/        -20.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 0)/        -31.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 1)/        -31.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 2)/        -15.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 0)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 1)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 0620.

      DATA cmatrx( 0, 6, 2, 0, 0, 0)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 1)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 2)/       7020.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 3)/       1980.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 4)/        369.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 5)/         45.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 6)/          3.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 0)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 1)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 2)/       7020.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 3)/       1980.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 4)/        369.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 5)/         45.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 6)/          3.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 0)/       7740.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 1)/       7740.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 2)/       3582.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 3)/       1002.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 4)/        183.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 5)/         21.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 0)/       2700.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 1)/       2700.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 2)/       1242.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 3)/        342.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 4)/         60.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 0)/        729.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 1)/        729.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 2)/        333.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 3)/         90.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 4)/         15.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 0)/        165.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 1)/        165.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 2)/         75.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 3)/         20.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 0)/         33.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 1)/         33.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 2)/         15.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 0)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 1)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 0)/          1.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1100.

      DATA cmatrx( 1, 1, 0, 0, 0, 0)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 0, 1)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 0, 2)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 0)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 1)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 2)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 0)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 1)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 1200.

      DATA cmatrx( 1, 2, 0, 0, 0, 0)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 1)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 2)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 0)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 1)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 2)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 0)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 1)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 2)/         -2.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 0)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 1)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1210.

      DATA cmatrx( 1, 2, 1, 0, 0, 0)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 1)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 2)/         12.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 3)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 0)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 1)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 2)/         12.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 3)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 0)/         15.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 1)/         15.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 2)/          6.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 0)/          5.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 1)/          5.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 2)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 0)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 1)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1211.

      DATA cmatrx( 1, 2, 1, 1, 0, 0)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 1)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 2)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 0)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 1)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 2)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 0)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 1)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 2)/         -2.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 0)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 1)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 0)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 1)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1300.

      DATA cmatrx( 1, 3, 0, 0, 0, 0)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 1)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 2)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 3)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 0)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 1)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 2)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 3)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 0)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 1)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 2)/        -12.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 0)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 1)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1310.

      DATA cmatrx( 1, 3, 1, 0, 0, 0)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 1)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 2)/         78.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 3)/         18.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 4)/          2.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 0)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 1)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 2)/         78.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 3)/         18.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 4)/          2.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 0)/         90.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 1)/         90.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 2)/         39.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 3)/          9.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 0)/         30.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 1)/         30.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 2)/         13.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 3)/          3.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 0)/          7.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 1)/          7.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 2)/          3.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 1)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1311.

      DATA cmatrx( 1, 3, 1, 1, 0, 0)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 1)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 2)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 3)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 0)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 1)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 2)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 3)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 0)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 1)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 2)/        -16.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 0)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 1)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 0)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 1)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1320.

      DATA cmatrx( 1, 3, 2, 0, 0, 0)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 1)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 2)/       -405.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 3)/        -90.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 0)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 1)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 2)/       -405.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 3)/        -90.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 0)/       -465.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 1)/       -465.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 2)/       -198.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 3)/        -43.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 0)/       -150.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 1)/       -150.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 2)/        -63.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 3)/        -13.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 0)/        -36.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 1)/        -36.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 2)/        -15.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 0)/         -7.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 1)/         -7.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1321.

      DATA cmatrx( 1, 3, 2, 1, 0, 0)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 1)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 2)/        135.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 3)/         30.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 4)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 0)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 1)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 2)/        135.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 3)/         30.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 4)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 0)/        150.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 1)/        150.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 2)/         63.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 3)/         13.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 4)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 0)/         45.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 1)/         45.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 2)/         18.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 3)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 0)/          9.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 1)/          9.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 2)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 0)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 1)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1400.

      DATA cmatrx( 1, 4, 0, 0, 0, 0)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 1)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 2)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 3)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 0)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 1)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 2)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 3)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 0)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 1)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 2)/        -78.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 3)/        -22.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 0)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 1)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 2)/        -22.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 0)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 1)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1410.

      DATA cmatrx( 1, 4, 1, 0, 0, 0)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 1)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 2)/        570.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 3)/        150.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 4)/         24.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 5)/          2.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 0)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 1)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 2)/        570.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 3)/        150.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 4)/         24.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 5)/          2.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 0)/        630.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 1)/        630.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 2)/        285.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 3)/         75.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 4)/         12.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 0)/        210.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 1)/        210.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 2)/         95.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 3)/         25.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 4)/          4.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 0)/         51.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 1)/         51.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 2)/         23.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 3)/          6.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 0)/          9.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 1)/          9.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 2)/          4.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 1)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1411.

      DATA cmatrx( 1, 4, 1, 1, 0, 0)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 1)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 2)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 3)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 4)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 0)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 1)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 2)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 3)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 4)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 0)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 1)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 2)/       -126.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 3)/        -31.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 0)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 1)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 2)/        -31.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 0)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 1)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1420.

      DATA cmatrx( 1, 4, 2, 0, 0, 0)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 1)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 2)/      -3375.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 3)/       -855.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 4)/       -126.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 5)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 0)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 1)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 2)/      -3375.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 3)/       -855.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 4)/       -126.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 5)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 0)/      -3735.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 1)/      -3735.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 2)/      -1662.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 3)/       -417.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 4)/        -60.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 5)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 0)/      -1215.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 1)/      -1215.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 2)/       -537.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 3)/       -132.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 4)/        -18.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 0)/       -294.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 1)/       -294.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 2)/       -129.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 3)/        -31.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 0)/        -57.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 1)/        -57.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 2)/        -25.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 0)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 1)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1421.

      DATA cmatrx( 1, 4, 2, 1, 0, 0)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 1)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 2)/       1125.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 3)/        285.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 4)/         42.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 5)/          3.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 0)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 1)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 2)/       1125.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 3)/        285.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 4)/         42.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 5)/          3.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 0)/       1215.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 1)/       1215.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 2)/        537.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 3)/        132.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 4)/         18.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 5)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 0)/        375.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 1)/        375.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 2)/        162.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 3)/         37.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 4)/          4.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 0)/         81.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 1)/         81.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 2)/         33.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 3)/          6.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 0)/         12.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 1)/         12.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 2)/          4.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 0)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 1)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1500.

      DATA cmatrx( 1, 5, 0, 0, 0, 0)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 1)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 2)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 3)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 4)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 5)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 0)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 1)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 2)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 3)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 4)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 5)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 0)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 1)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 2)/       -570.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 3)/       -170.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 4)/        -35.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 0)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 1)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 2)/       -170.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 3)/        -50.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 0)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 1)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 2)/        -35.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 0)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 1)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1510.

      DATA cmatrx( 1, 5, 1, 0, 0, 0)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 1)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 2)/       4680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 3)/       1320.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 4)/        246.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 5)/         30.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 6)/          2.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 0)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 1)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 2)/       4680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 3)/       1320.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 4)/        246.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 5)/         30.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 6)/          2.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 0)/       5040.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 1)/       5040.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 2)/       2340.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 3)/        660.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 4)/        123.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 5)/         15.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 0)/       1680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 1)/       1680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 2)/        780.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 3)/        220.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 4)/         41.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 5)/          5.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 0)/        414.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 1)/        414.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 2)/        192.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 3)/         54.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 4)/         10.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 0)/         78.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 1)/         78.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 2)/         36.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 3)/         10.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 0)/         11.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 1)/         11.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 2)/          5.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 0)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 1)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1511.

      DATA cmatrx( 1, 5, 1, 1, 0, 0)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 1)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 2)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 3)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 4)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 5)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 0)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 1)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 2)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 3)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 4)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 5)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 0)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 1)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 2)/      -1074.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 3)/       -294.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 4)/        -51.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 0)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 1)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 2)/       -294.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 3)/        -74.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 0)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 1)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 2)/        -51.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 0)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 1)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1520.

      DATA cmatrx( 1, 5, 2, 0, 0, 0)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 1)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 2)/     -31185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 3)/      -8505.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 4)/      -1485.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 5)/       -162.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 6)/         -9.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 0)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 1)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 2)/     -31185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 3)/      -8505.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 4)/      -1485.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 5)/       -162.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 6)/         -9.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 0)/     -33705.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 1)/     -33705.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 2)/     -15420.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 3)/      -4185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 4)/       -723.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 5)/        -77.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 6)/         -4.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 0)/     -11025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 1)/     -11025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 2)/      -5025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 3)/      -1350.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 4)/       -228.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 5)/        -23.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 0)/      -2685.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 1)/      -2685.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 2)/      -1218.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 3)/       -323.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 4)/        -53.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 0)/       -522.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 1)/       -522.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 2)/       -236.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 3)/        -62.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 0)/        -84.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 1)/        -84.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 2)/        -38.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 0)/        -11.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 1)/        -11.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1521.

      DATA cmatrx( 1, 5, 2, 1, 0, 0)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 1)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 2)/      10395.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 3)/       2835.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 4)/        495.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 5)/         54.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 6)/          3.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 0)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 1)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 2)/      10395.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 3)/       2835.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 4)/        495.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 5)/         54.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 6)/          3.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 0)/      11025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 1)/      11025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 2)/       5025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 3)/       1350.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 4)/        228.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 5)/         23.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 6)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 0)/       3465.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 1)/       3465.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 2)/       1560.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 3)/        405.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 4)/         63.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 5)/          5.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 0)/        780.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 1)/        780.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 2)/        342.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 3)/         82.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 4)/         10.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 0)/        129.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 1)/        129.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 2)/         53.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 3)/         10.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 0)/         15.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 1)/         15.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 2)/          5.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 0)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 1)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1600.

      DATA cmatrx( 1, 6, 0, 0, 0, 0)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 1)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 2)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 3)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 4)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 5)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 6)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 0)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 1)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 2)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 3)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 4)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 5)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 6)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 0)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 1)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 2)/      -4680.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 3)/      -1440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 4)/       -315.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 5)/        -51.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 0)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 1)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 2)/      -1440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 3)/       -440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 4)/        -95.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 0)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 1)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 2)/       -315.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 3)/        -95.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 0)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 1)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 2)/        -51.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 0)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 1)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1610.

      DATA cmatrx( 1, 6, 1, 0, 0, 0)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 1)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 2)/      42840.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 3)/      12600.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 4)/       2550.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 5)/        366.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 6)/         36.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 7)/          2.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 0)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 1)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 2)/      42840.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 3)/      12600.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 4)/       2550.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 5)/        366.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 6)/         36.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 7)/          2.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 0)/      45360.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 1)/      45360.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 2)/      21420.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 3)/       6300.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 4)/       1275.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 5)/        183.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 6)/         18.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 7)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 0)/      15120.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 1)/      15120.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 2)/       7140.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 3)/       2100.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 4)/        425.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 5)/         61.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 6)/          6.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 0)/       3750.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 1)/       3750.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 2)/       1770.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 3)/        520.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 4)/        105.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 5)/         15.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 0)/        726.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 1)/        726.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 2)/        342.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 3)/        100.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 4)/         20.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 0)/        111.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 1)/        111.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 2)/         52.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 3)/         15.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 0)/         13.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 1)/         13.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 2)/          6.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 0)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 1)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1611.

      DATA cmatrx( 1, 6, 1, 1, 0, 0)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 1)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 2)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 3)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 4)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 5)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 6)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 0)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 1)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 2)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 3)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 4)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 5)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 6)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 0)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 1)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 2)/     -10050.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 3)/      -2910.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 4)/       -570.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 5)/        -76.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 0)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 1)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 2)/      -2910.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 3)/       -810.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 4)/       -145.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 0)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 1)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 2)/       -570.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 3)/       -145.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 0)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 1)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 2)/        -76.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 0)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 1)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1620.

      DATA cmatrx( 1, 6, 2, 0, 0, 0)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 1)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 2)/    -317520.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 3)/     -90720.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 4)/     -17415.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 5)/      -2295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 6)/       -198.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 7)/         -9.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 0)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 1)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 2)/    -317520.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 3)/     -90720.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 4)/     -17415.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 5)/      -2295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 6)/       -198.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 7)/         -9.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 0)/    -337680.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 1)/    -337680.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 2)/    -157410.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 3)/     -44850.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 4)/      -8562.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 5)/      -1116.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 6)/        -94.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 7)/         -4.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 0)/    -110880.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 1)/    -110880.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 2)/     -51570.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 3)/     -14610.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 4)/      -2757.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 5)/       -351.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 6)/        -28.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 0)/     -27135.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 1)/     -27135.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 2)/     -12582.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 3)/      -3537.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 4)/       -657.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 5)/        -81.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 0)/      -5295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 1)/      -5295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 2)/      -2448.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 3)/       -683.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 4)/       -125.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 0)/       -858.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 1)/       -858.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 2)/       -396.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 3)/       -110.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 0)/       -117.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 1)/       -117.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 2)/        -54.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 0)/        -13.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 1)/        -13.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1621.

      DATA cmatrx( 1, 6, 2, 1, 0, 0)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 1)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 2)/     105840.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 3)/      30240.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 4)/       5805.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 5)/        765.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 6)/         66.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 7)/          3.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 0)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 1)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 2)/     105840.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 3)/      30240.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 4)/       5805.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 5)/        765.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 6)/         66.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 7)/          3.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 0)/     110880.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 1)/     110880.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 2)/      51570.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 3)/      14610.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 4)/       2757.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 5)/        351.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 6)/         28.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 7)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 0)/      35280.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 1)/      35280.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 2)/      16290.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 3)/       4530.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 4)/        822.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 5)/         96.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 6)/          6.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 0)/       8145.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 1)/       8145.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 2)/       3708.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 3)/        993.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 4)/        165.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 5)/         15.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 0)/       1425.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 1)/       1425.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 2)/        630.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 3)/        155.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 4)/         20.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 0)/        189.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 1)/        189.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 2)/         78.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 3)/         15.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 0)/         18.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 1)/         18.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 2)/          6.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 0)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 1)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2100.

      DATA cmatrx( 2, 1, 0, 0, 0, 0)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 1)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 2)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 3)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 0)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 1)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 2)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 3)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 0)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 1)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 2)/          7.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 0)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 1)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 2)/          1.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 2200.

      DATA cmatrx( 2, 2, 0, 0, 0, 0)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 1)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 2)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 3)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 0)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 1)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 2)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 3)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 0)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 1)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 2)/         50.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 3)/         11.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 0)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 1)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 2)/         11.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 3)/          2.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 0)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 1)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 2)/          1.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2210.

      DATA cmatrx( 2, 2, 1, 0, 0, 0)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 1)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 2)/       -405.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 3)/        -90.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 0)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 1)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 2)/       -405.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 3)/        -90.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 0)/       -450.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 1)/       -450.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 2)/       -192.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 3)/        -42.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 0)/       -135.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 1)/       -135.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 2)/        -57.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 3)/        -12.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 4)/         -1.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 0)/        -27.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 1)/        -27.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 2)/        -11.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 3)/         -2.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 0)/         -3.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 1)/         -3.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 2)/         -1.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2211.

      DATA cmatrx( 2, 2, 1, 1, 0, 0)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 1)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 2)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 3)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 4)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 0)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 1)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 2)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 3)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 4)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 0)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 1)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 2)/        171.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 3)/         36.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 0)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 1)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 2)/         36.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 3)/          6.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 0)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 1)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 2)/          3.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 0)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 1)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2300.

      DATA cmatrx( 2, 3, 0, 0, 0, 0)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 1)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 2)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 3)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 4)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 0)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 1)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 2)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 3)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 4)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 0)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 1)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 2)/        384.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 3)/         99.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 4)/         15.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 0)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 1)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 2)/         99.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 3)/         24.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 4)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 0)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 1)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 2)/         15.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 3)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 0)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 1)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 2)/          1.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2310.

      DATA cmatrx( 2, 3, 1, 0, 0, 0)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 1)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 2)/      -3375.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 3)/       -855.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 4)/       -126.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 5)/         -9.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 0)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 1)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 2)/      -3375.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 3)/       -855.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 4)/       -126.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 5)/         -9.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 0)/      -3645.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 1)/      -3645.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 2)/      -1623.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 3)/       -408.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 4)/        -59.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 5)/         -4.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 0)/      -1125.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 1)/      -1125.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 2)/       -498.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 3)/       -123.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 4)/        -17.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 5)/         -1.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 0)/       -243.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 1)/       -243.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 2)/       -106.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 3)/        -25.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 4)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 0)/        -36.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 1)/        -36.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 2)/        -15.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 3)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 0)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 1)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 2)/         -1.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2311.

      DATA cmatrx( 2, 3, 1, 1, 0, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 2)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 3)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 4)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 2)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 3)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 4)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 0)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 1)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 2)/       1494.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 3)/        369.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 4)/         51.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 0)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 1)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 2)/        369.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 3)/         84.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 4)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 0)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 1)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 2)/         51.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 3)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 0)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 1)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 2)/          3.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2320.

      DATA cmatrx( 2, 3, 2, 0, 0, 0)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 1)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 2)/      22680.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 3)/       5670.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 4)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 5)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 0)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 1)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 2)/      22680.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 3)/       5670.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 4)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 5)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 0)/      24570.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 1)/      24570.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 2)/      10890.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 3)/       2700.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 4)/        378.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 5)/         24.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 2)/       3330.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 3)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 4)/        108.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 5)/          6.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 0)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 1)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 2)/        726.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 3)/        171.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 4)/         21.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 5)/          1.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 0)/        279.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 1)/        279.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 2)/        120.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 3)/         27.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 4)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 0)/         36.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 1)/         36.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 2)/         15.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 3)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 0)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 1)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 2)/          1.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2321.

      DATA cmatrx( 2, 3, 2, 1, 0, 0)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 1)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 2)/     -15120.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 3)/      -3780.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 4)/       -540.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 5)/        -36.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 0)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 1)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 2)/     -15120.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 3)/      -3780.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 4)/       -540.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 5)/        -36.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 0)/     -16065.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 1)/     -16065.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 2)/      -7110.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 3)/      -1755.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 4)/       -243.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 5)/        -15.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 0)/      -4725.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 1)/      -4725.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 2)/      -2070.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 3)/       -495.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 4)/        -63.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 5)/         -3.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 0)/       -945.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 1)/       -945.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 2)/       -405.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 3)/        -90.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 4)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 0)/       -126.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 1)/       -126.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 2)/        -51.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 3)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 0)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 1)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 2)/         -3.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2322.

      DATA cmatrx( 2, 3, 2, 2, 0, 0)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 1)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 2)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 3)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 4)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 0)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 1)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 2)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 3)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 4)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 0)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 1)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 2)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 3)/        405.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 4)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 0)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 1)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 2)/        405.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 3)/         90.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 4)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 0)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 1)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 2)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 3)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 0)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 1)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 2)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2400.

      DATA cmatrx( 2, 4, 0, 0, 0, 0)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 1)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 2)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 3)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 4)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 5)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 0)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 1)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 2)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 3)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 4)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 5)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 0)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 1)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 2)/       3246.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 3)/        906.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 4)/        165.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 5)/         19.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 0)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 1)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 2)/        906.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 3)/        246.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 4)/         42.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 5)/          4.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 0)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 1)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 2)/        165.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 3)/         42.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 4)/          6.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 0)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 1)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 2)/         19.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 3)/          4.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 0)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 1)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 2)/          1.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2410.

      DATA cmatrx( 2, 4, 1, 0, 0, 0)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 1)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 2)/     -31185.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 3)/      -8505.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 4)/      -1485.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 5)/       -162.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 6)/         -9.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 0)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 1)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 2)/     -31185.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 3)/      -8505.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 4)/      -1485.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 5)/       -162.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 6)/         -9.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 0)/     -33075.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 1)/     -33075.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 2)/     -15135.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 3)/      -4110.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 4)/       -711.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 5)/        -76.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 6)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 0)/     -10395.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 1)/     -10395.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 2)/      -4740.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 3)/      -1275.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 4)/       -216.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 5)/        -22.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 6)/         -1.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 0)/      -2340.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 1)/      -2340.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 2)/      -1059.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 3)/       -279.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 4)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 5)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 0)/       -387.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 1)/       -387.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 2)/       -172.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 3)/        -43.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 4)/         -6.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 0)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 1)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 2)/        -19.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 3)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 0)/         -3.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 1)/         -3.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 2)/         -1.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2411.

      DATA cmatrx( 2, 4, 1, 1, 0, 0)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 1)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 2)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 3)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 4)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 5)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 0)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 1)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 2)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 3)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 4)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 5)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 0)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 1)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 2)/      14220.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 3)/       3825.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 4)/        648.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 5)/         66.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 0)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 1)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 2)/       3825.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 3)/        990.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 4)/        153.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 5)/         12.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 0)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 1)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 2)/        648.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 3)/        153.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 4)/         18.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 0)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 1)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 2)/         66.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 3)/         12.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 0)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 1)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 2)/          3.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2420.

      DATA cmatrx( 2, 4, 2, 0, 0, 0)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 1)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 2)/     232470.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 3)/      62370.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 4)/      10530.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 5)/       1080.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 6)/         54.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 0)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 1)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 2)/     232470.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 3)/      62370.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 4)/      10530.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 5)/       1080.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 6)/         54.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 0)/     247590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 1)/     247590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 2)/     112590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 3)/      30060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 4)/       5022.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 5)/        504.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 6)/         24.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 0)/      77490.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 1)/      77490.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 2)/      35100.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 3)/       9270.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 4)/       1512.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 5)/        144.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 6)/          6.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 0)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 1)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 2)/       7899.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 3)/       2049.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 4)/        321.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 5)/         28.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 6)/          1.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 0)/       3060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 1)/       3060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 2)/       1365.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 3)/        345.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 4)/         51.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 5)/          4.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 0)/        423.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 1)/        423.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 2)/        186.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 3)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 4)/          6.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 0)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 1)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 2)/         19.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 3)/          4.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 0)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 1)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 2)/          1.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2421.

      DATA cmatrx( 2, 4, 2, 1, 0, 0)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 1)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 2)/    -154980.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 3)/     -41580.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 4)/      -7020.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 5)/       -720.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 6)/        -36.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 0)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 1)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 2)/    -154980.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 3)/     -41580.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 4)/      -7020.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 5)/       -720.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 6)/        -36.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 0)/    -162540.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 1)/    -162540.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 2)/     -73845.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 3)/     -19665.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 4)/      -3267.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 5)/       -324.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 6)/        -15.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 0)/     -49140.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 1)/     -49140.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 2)/     -22185.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 3)/      -5805.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 4)/       -927.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 5)/        -84.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 6)/         -3.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 0)/     -10395.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 1)/     -10395.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 2)/      -4635.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 3)/      -1170.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 4)/       -171.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 5)/        -12.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 0)/      -1575.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 1)/      -1575.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 2)/       -684.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 3)/       -159.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 4)/        -18.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 0)/       -162.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 1)/       -162.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 2)/        -66.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 3)/        -12.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 0)/         -9.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 1)/         -9.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 2)/         -3.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2422.

      DATA cmatrx( 2, 4, 2, 2, 0, 0)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 1)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 2)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 3)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 4)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 5)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 0)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 1)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 2)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 3)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 4)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 5)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 0)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 1)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 2)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 3)/       4635.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 4)/        756.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 5)/         72.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 0)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 1)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 2)/       4635.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 3)/       1170.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 4)/        171.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 5)/         12.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 0)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 1)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 2)/        756.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 3)/        171.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 4)/         18.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 0)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 1)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 2)/         72.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 3)/         12.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 0)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 1)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 2)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2500.

      DATA cmatrx( 2, 5, 0, 0, 0, 0)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 1)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 2)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 3)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 4)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 5)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 6)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 0)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 1)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 2)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 3)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 4)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 5)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 6)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 0)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 1)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 2)/      30270.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 3)/       8850.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 4)/       1770.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 5)/        248.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 6)/         23.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 7)/          1.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 0)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 1)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 2)/       8850.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 3)/       2550.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 4)/        495.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 5)/         65.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 6)/          5.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 0)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 1)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 2)/       1770.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 3)/        495.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 4)/         90.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 5)/         10.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 0)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 1)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 2)/        248.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 3)/         65.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 4)/         10.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 0)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 1)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 2)/         23.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 3)/          5.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 0)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 1)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 2)/          1.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2510.

      DATA cmatrx( 2, 5, 1, 0, 0, 0)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 1)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 2)/    -317520.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 3)/     -90720.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 4)/     -17415.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 5)/      -2295.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 6)/       -198.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 7)/         -9.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 0)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 1)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 2)/    -317520.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 3)/     -90720.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 4)/     -17415.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 5)/      -2295.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 6)/       -198.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 7)/         -9.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 0)/    -332640.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 1)/    -332640.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 2)/    -155070.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 3)/     -44190.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 4)/      -8439.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 5)/      -1101.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 6)/        -93.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 7)/         -4.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 0)/    -105840.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 1)/    -105840.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 2)/     -49230.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 3)/     -13950.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 4)/      -2634.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 5)/       -336.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 6)/        -27.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 7)/         -1.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 0)/     -24435.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 1)/     -24435.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 2)/     -11316.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 3)/      -3171.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 4)/       -585.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 5)/        -71.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 6)/         -5.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 0)/      -4275.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 1)/      -4275.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 2)/      -1962.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 3)/       -537.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 4)/        -94.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 5)/        -10.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 0)/       -567.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 1)/       -567.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 2)/       -255.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 3)/        -66.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 4)/        -10.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 0)/        -54.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 1)/        -54.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 2)/        -23.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 3)/         -5.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 0)/         -3.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 1)/         -3.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 2)/         -1.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2511.

      DATA cmatrx( 2, 5, 1, 1, 0, 0)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 1)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 2)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 3)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 4)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 5)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 6)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 0)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 1)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 2)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 3)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 4)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 5)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 6)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 0)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 1)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 2)/     147690.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 3)/      41850.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 4)/       7902.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 5)/       1008.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 6)/         81.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 0)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 1)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 2)/      41850.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 3)/      11610.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 4)/       2097.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 5)/        243.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 6)/         15.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 0)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 1)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 2)/       7902.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 3)/       2097.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 4)/        342.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 5)/         30.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 0)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 1)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 2)/       1008.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 3)/        243.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 4)/         30.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 0)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 1)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 2)/         81.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 3)/         15.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 0)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 1)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 2)/          3.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2520.

      DATA cmatrx( 2, 5, 2, 0, 0, 0)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 1)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 2)/    2602530.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 3)/     731430.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 4)/     136080.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 5)/      17010.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 6)/       1350.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 7)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 0)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 1)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 2)/    2602530.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 3)/     731430.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 4)/     136080.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 5)/      17010.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 6)/       1350.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 7)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 0)/    2738610.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 1)/    2738610.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 2)/    1268190.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 3)/     355320.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 4)/      65700.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 5)/       8118.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 6)/        630.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 7)/         24.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 0)/     867510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 1)/     867510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 2)/     400680.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 3)/     111510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 4)/      20340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 5)/       2448.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 6)/        180.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 7)/          6.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 0)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 1)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 2)/      92145.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 3)/      25365.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 4)/       4524.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 5)/        521.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 6)/         35.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 7)/          1.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 0)/      35910.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 1)/      35910.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 2)/      16419.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 3)/       4449.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 4)/        768.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 5)/         83.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 6)/          5.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 0)/       5175.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 1)/       5175.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 2)/       2346.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 3)/        621.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 4)/        102.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 5)/         10.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 0)/        603.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 1)/        603.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 2)/        269.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 3)/         68.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 4)/         10.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 0)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 1)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 2)/         23.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 3)/          5.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 0)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 1)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 2)/          1.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2521.

      DATA cmatrx( 2, 5, 2, 1, 0, 0)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 0, 1)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 0, 2)/   -1735019.99999999/
      DATA cmatrx( 2, 5, 2, 1, 0, 3)/    -487619.99999999/
      DATA cmatrx( 2, 5, 2, 1, 0, 4)/     -90720.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 5)/     -11340.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 6)/       -900.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 7)/        -36.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 0)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 1, 1)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 1, 2)/   -1735019.99999998/
      DATA cmatrx( 2, 5, 2, 1, 1, 3)/    -487619.99999999/
      DATA cmatrx( 2, 5, 2, 1, 1, 4)/     -90720.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 5)/     -11340.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 6)/       -900.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 7)/        -36.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 0)/   -1803059.99999999/
      DATA cmatrx( 2, 5, 2, 1, 2, 1)/   -1803059.99999998/
      DATA cmatrx( 2, 5, 2, 1, 2, 2)/    -834434.99999999/
      DATA cmatrx( 2, 5, 2, 1, 2, 3)/    -233415.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 4)/     -43020.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 5)/      -5283.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 6)/       -405.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 7)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 0)/    -555659.99999999/
      DATA cmatrx( 2, 5, 2, 1, 3, 1)/    -555659.99999999/
      DATA cmatrx( 2, 5, 2, 1, 3, 2)/    -256095.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 3)/     -70875.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 4)/     -12780.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 5)/      -1503.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 6)/       -105.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 7)/         -3.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 0)/    -121905.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 1)/    -121905.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 2)/     -55755.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 3)/     -15120.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 4)/      -2610.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 5)/       -279.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 6)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 0)/     -19845.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 1)/     -19845.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 2)/      -8946.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 3)/      -2331.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 4)/       -366.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 5)/        -30.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 0)/      -2385.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 1)/      -2385.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 2)/      -1044.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 3)/       -249.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 4)/        -30.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 0)/       -198.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 1)/       -198.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 2)/        -81.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 3)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 0)/         -9.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 1)/         -9.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 2)/         -3.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2522.

      DATA cmatrx( 2, 5, 2, 2, 0, 0)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 0, 1)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 0, 2)/     433755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 3)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 4)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 5)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 6)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 0)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 1)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 2)/     433754.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 3)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 4)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 5)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 6)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 0)/     433755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 1)/     433754.99999999/
      DATA cmatrx( 2, 5, 2, 2, 2, 2)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 3)/      55755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 4)/      10170.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 5)/       1224.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 6)/         90.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 0)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 1)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 2)/      55755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 3)/      15120.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 4)/       2610.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 5)/        279.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 6)/         15.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 0)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 1)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 2)/      10170.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 3)/       2610.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 4)/        396.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 5)/         30.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 0)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 1)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 2)/       1224.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 3)/        279.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 4)/         30.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 0)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 1)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 2)/         90.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 3)/         15.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 0)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 1)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 2)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2600.

      DATA cmatrx( 2, 6, 0, 0, 0, 0)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 1)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 2)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 3)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 4)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 5)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 6)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 7)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 0)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 1)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 2)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 3)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 4)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 5)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 6)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 7)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 0)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 1)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 2)/     310140.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 3)/      93420.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 4)/      19755.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 5)/       3063.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 6)/        348.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 7)/         27.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 8)/          1.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 0)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 1)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 2)/      93420.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 3)/      27900.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 4)/       5805.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 5)/        873.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 6)/         93.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 7)/          6.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 0)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 1)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 2)/      19755.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 3)/       5805.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 4)/       1170.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 5)/        165.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 6)/         15.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 0)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 1)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 2)/       3063.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 3)/        873.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 4)/        165.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 5)/         20.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 0)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 1)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 2)/        348.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 3)/         93.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 4)/         15.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 0)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 1)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 2)/         27.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 3)/          6.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 0)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 1)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 2)/          1.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2610.

      DATA cmatrx( 2, 6, 1, 0, 0, 0)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 1)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 2)/   -3538080.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 3)/   -1043279.99999999/
      DATA cmatrx( 2, 6, 1, 0, 0, 4)/    -212625.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 5)/     -31185.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 6)/      -3285.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 7)/       -234.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 8)/         -9.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 0)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 1)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 2)/   -3538080.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 3)/   -1043279.99999999/
      DATA cmatrx( 2, 6, 1, 0, 1, 4)/    -212625.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 5)/     -31185.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 6)/      -3285.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 7)/       -234.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 8)/         -9.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 0)/   -3674160.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 1)/   -3674160.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 2)/   -1735650.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 3)/    -510930.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 4)/    -103800.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 5)/     -15138.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 6)/      -1578.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 7)/       -110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 8)/         -4.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 0)/   -1179359.99999999/
      DATA cmatrx( 2, 6, 1, 0, 3, 1)/   -1179359.99999999/
      DATA cmatrx( 2, 6, 1, 0, 3, 2)/    -556290.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 3)/    -163170.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 4)/     -32925.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 5)/      -4743.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 6)/       -483.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 7)/        -32.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 8)/         -1.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 0)/    -276885.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 1)/    -276885.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 2)/    -130245.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 3)/     -37950.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 4)/      -7560.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 5)/      -1063.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 6)/       -103.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 7)/         -6.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 0)/     -50085.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 1)/     -50085.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 2)/     -23439.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 3)/      -6744.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 4)/      -1310.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 5)/       -175.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 6)/        -15.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 0)/      -7110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 1)/      -7110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 2)/      -3294.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 3)/       -924.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 4)/       -170.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 5)/        -20.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 0)/       -783.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 1)/       -783.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 2)/       -355.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 3)/        -94.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 4)/        -15.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 0)/        -63.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 1)/        -63.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 2)/        -27.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 3)/         -6.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 0)/         -3.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 1)/         -3.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 2)/         -1.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2611.

      DATA cmatrx( 2, 6, 1, 1, 0, 0)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 1)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 2)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 3)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 0, 4)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 5)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 6)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 7)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 0)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 1)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 2)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 3)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 1, 4)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 5)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 6)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 7)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 0)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 1)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 2)/    1668870.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 3)/     489510.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 4)/      98775.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 5)/      14229.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 6)/       1449.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 7)/         96.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 0)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 3, 1)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 3, 2)/     489510.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 3)/     141750.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 4)/      27900.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 5)/       3834.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 6)/        354.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 7)/         18.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 0)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 1)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 2)/      98775.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 3)/      27900.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 4)/       5220.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 5)/        645.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 6)/         45.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 0)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 1)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 2)/      14229.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 3)/       3834.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 4)/        645.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 5)/         60.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 0)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 1)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 2)/       1449.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 3)/        354.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 4)/         45.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 0)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 1)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 2)/         96.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 3)/         18.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 0)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 1)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 2)/          3.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2620.

      DATA cmatrx( 2, 6, 2, 0, 0, 0)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 1)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 2)/   31638600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 3)/    9185399.99999988/
      DATA cmatrx( 2, 6, 2, 0, 0, 4)/    1820069.99999997/
      DATA cmatrx( 2, 6, 2, 0, 0, 5)/     255149.99999999/
      DATA cmatrx( 2, 6, 2, 0, 0, 6)/      25110.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 7)/       1620.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 8)/         54.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 0)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 1)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 2)/   31638599.99999970/
      DATA cmatrx( 2, 6, 2, 0, 1, 3)/    9185399.99999988/
      DATA cmatrx( 2, 6, 2, 0, 1, 4)/    1820069.99999997/
      DATA cmatrx( 2, 6, 2, 0, 1, 5)/     255150.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 6)/      25110.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 7)/       1620.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 8)/         54.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 0)/   32999400.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 1)/   32999399.99999970/
      DATA cmatrx( 2, 6, 2, 0, 2, 2)/   15486659.99999980/
      DATA cmatrx( 2, 6, 2, 0, 2, 3)/    4486859.99999994/
      DATA cmatrx( 2, 6, 2, 0, 2, 4)/     885599.99999999/
      DATA cmatrx( 2, 6, 2, 0, 2, 5)/     123300.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 6)/      11988.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 7)/        756.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 8)/         24.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 0)/   10546199.99999980/
      DATA cmatrx( 2, 6, 2, 0, 3, 1)/   10546199.99999980/
      DATA cmatrx( 2, 6, 2, 0, 3, 2)/    4940459.99999994/
      DATA cmatrx( 2, 6, 2, 0, 3, 3)/    1425059.99999998/
      DATA cmatrx( 2, 6, 2, 0, 3, 4)/     278909.99999999/
      DATA cmatrx( 2, 6, 2, 0, 3, 5)/      38250.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 6)/       3618.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 7)/        216.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 8)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 0)/    2470229.99999997/
      DATA cmatrx( 2, 6, 2, 0, 4, 1)/    2470229.99999997/
      DATA cmatrx( 2, 6, 2, 0, 4, 2)/    1153889.99999999/
      DATA cmatrx( 2, 6, 2, 0, 4, 3)/     330479.99999999/
      DATA cmatrx( 2, 6, 2, 0, 4, 4)/      63801.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 5)/       8535.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 6)/        771.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 7)/         42.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 8)/          1.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 0)/     451709.99999999/
      DATA cmatrx( 2, 6, 2, 0, 5, 1)/     451710.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 2)/     210150.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 3)/      59580.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 4)/      11277.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 5)/       1455.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 6)/        123.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 7)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 0)/      66960.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 1)/      66960.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 2)/      30978.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 3)/       8658.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 4)/       1593.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 5)/        195.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 6)/         15.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 0)/       8190.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 1)/       8190.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 2)/       3756.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 3)/       1026.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 4)/        180.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 5)/         20.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 0)/        819.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 1)/        819.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 2)/        369.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 3)/         96.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 4)/         15.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 0)/         63.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 1)/         63.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 2)/         27.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 3)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 0)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 1)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 2)/          1.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2621.

      DATA cmatrx( 2, 6, 2, 1, 0, 0)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 0, 1)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 0, 2)/  -21092399.99999970/
      DATA cmatrx( 2, 6, 2, 1, 0, 3)/   -6123599.99999985/
      DATA cmatrx( 2, 6, 2, 1, 0, 4)/   -1213379.99999996/
      DATA cmatrx( 2, 6, 2, 1, 0, 5)/    -170099.99999999/
      DATA cmatrx( 2, 6, 2, 1, 0, 6)/     -16740.00000000/
      DATA cmatrx( 2, 6, 2, 1, 0, 7)/      -1080.00000000/
      DATA cmatrx( 2, 6, 2, 1, 0, 8)/        -36.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 0)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 1, 1)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 1, 2)/  -21092399.99999960/
      DATA cmatrx( 2, 6, 2, 1, 1, 3)/   -6123599.99999985/
      DATA cmatrx( 2, 6, 2, 1, 1, 4)/   -1213379.99999998/
      DATA cmatrx( 2, 6, 2, 1, 1, 5)/    -170099.99999999/
      DATA cmatrx( 2, 6, 2, 1, 1, 6)/     -16740.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 7)/      -1080.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 8)/        -36.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 0)/  -21772799.99999970/
      DATA cmatrx( 2, 6, 2, 1, 2, 1)/  -21772799.99999960/
      DATA cmatrx( 2, 6, 2, 1, 2, 2)/  -10213559.99999980/
      DATA cmatrx( 2, 6, 2, 1, 2, 3)/   -2955959.99999994/
      DATA cmatrx( 2, 6, 2, 1, 2, 4)/    -582254.99999998/
      DATA cmatrx( 2, 6, 2, 1, 2, 5)/     -80775.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 6)/      -7803.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 7)/       -486.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 8)/        -15.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 0)/   -6803999.99999985/
      DATA cmatrx( 2, 6, 2, 1, 3, 1)/   -6803999.99999985/
      DATA cmatrx( 2, 6, 2, 1, 3, 2)/   -3182759.99999991/
      DATA cmatrx( 2, 6, 2, 1, 3, 3)/    -914759.99999997/
      DATA cmatrx( 2, 6, 2, 1, 3, 4)/    -177794.99999999/
      DATA cmatrx( 2, 6, 2, 1, 3, 5)/     -24075.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 6)/      -2223.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 7)/       -126.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 8)/         -3.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 0)/   -1530899.99999996/
      DATA cmatrx( 2, 6, 2, 1, 4, 1)/   -1530899.99999996/
      DATA cmatrx( 2, 6, 2, 1, 4, 2)/    -712529.99999998/
      DATA cmatrx( 2, 6, 2, 1, 4, 3)/    -202229.99999999/
      DATA cmatrx( 2, 6, 2, 1, 4, 4)/     -38340.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 5)/      -4950.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 6)/       -414.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 7)/        -18.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 0)/    -260819.99999999/
      DATA cmatrx( 2, 6, 2, 1, 5, 1)/    -260819.99999999/
      DATA cmatrx( 2, 6, 2, 1, 5, 2)/    -120330.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 3)/     -33390.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 4)/      -6039.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 5)/       -705.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 6)/        -45.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 0)/     -34155.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 1)/     -34155.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 2)/     -15507.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 3)/      -4122.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 4)/       -675.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 5)/        -60.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 0)/      -3375.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 1)/      -3375.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 2)/      -1485.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 3)/       -360.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 4)/        -45.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 0)/       -234.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 1)/       -234.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 2)/        -96.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 3)/        -18.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 0)/         -9.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 1)/         -9.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 2)/         -3.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2622.

      DATA cmatrx( 2, 6, 2, 2, 0, 0)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 0, 1)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 0, 2)/    5273099.99999994/
      DATA cmatrx( 2, 6, 2, 2, 0, 3)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 0, 4)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 0, 5)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 6)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 7)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 0)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 1, 1)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 1, 2)/    5273099.99999991/
      DATA cmatrx( 2, 6, 2, 2, 1, 3)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 1, 4)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 1, 5)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 6)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 7)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 0)/    5273099.99999994/
      DATA cmatrx( 2, 6, 2, 2, 2, 1)/    5273099.99999991/
      DATA cmatrx( 2, 6, 2, 2, 2, 2)/    2470229.99999996/
      DATA cmatrx( 2, 6, 2, 2, 2, 3)/     712529.99999998/
      DATA cmatrx( 2, 6, 2, 2, 2, 4)/     139455.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 5)/      19125.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 6)/       1809.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 7)/        108.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 0)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 3, 1)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 3, 2)/     712529.99999998/
      DATA cmatrx( 2, 6, 2, 2, 3, 3)/     202229.99999999/
      DATA cmatrx( 2, 6, 2, 2, 3, 4)/      38340.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 5)/       4950.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 6)/        414.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 7)/         18.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 0)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 4, 1)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 4, 2)/     139455.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 3)/      38340.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 4)/       6804.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 5)/        765.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 6)/         45.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 0)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 1)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 2)/      19125.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 3)/       4950.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 4)/        765.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 5)/         60.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 0)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 1)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 2)/       1809.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 3)/        414.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 4)/         45.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 0)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 1)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 2)/        108.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 3)/         18.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 0)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 1)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 2)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 8)/          0.00000000/

      END
      BLOCK DATA set2

      implicit none

      integer i
     
      double precision fact(0:12)
      double precision vsip(9,2)
      
      common /factorial/ fact
      common /coulom/ vsip

      data fact/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
     *          362880.0,3628800.0,39916800.0,479001600.0/

      data (vsip(1,i),i=1,2) /- 8.83, 0.00/
      data (vsip(2,i),i=1,2) /-18.62,-8.33/
      data (vsip(3,i),i=1,2) /-8.200,-9.10/
      data (vsip(4,i),i=1,2) /-8.65 ,-9.35 /
      data (vsip(5,i),i=1,2) /-8.150,-9.150/
      data (vsip(6,i),i=1,2) /-8.650,-9.450/
      data (vsip(7,i),i=1,2) /-8.05,-9.25/
      data (vsip(8,i),i=1,2) /-8.650,-9.650/
      data (vsip(9,i),i=1,2) /-10.620,-5.986/

      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine c_c_hd (hamil,dhamil)
c
      implicit none
      character*5 tc
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss
      double precision sps0,   spsna,spsnb,   spsnc,sps,dsps
      double precision pss0,   pssna,pssnb,   pssnc,pss,dpss
      double precision pps0,   ppsna,ppsnb,   ppsnc,pps,dpps
      double precision ppp0,   pppna,pppnb,   pppnc,ppp,dppp
c
c Hopping parameters for carbon-hydrogen from Wang and Mak.
c
      DATA r0 /1.312d0/
      DATA rt /2.d0/
c s s sigma
c      DATA sss0 /-7.948d0/
      DATA sss0 /-8.42256d0/
      DATA sssna /1.29827d0/
      DATA sssnb /1.d0/
      DATA sssnc /5.d0/
c s p sigma
c      DATA sps0 /7.413d0/
      DATA sps0 /8.08162d0/
      DATA spsna /0.99055d0/
      DATA spsnb /1.d0/
      DATA spsnc /5.d0/
c p s sigma
c      DATA pss0 /-7.413d0/
      DATA pss0 /-8.08162d0/
      DATA pssna /0.99055d0/
      DATA pssnb /1.d0/
      DATA pssnc /5.d0/
c p p sigma
c      DATA pps0 /6.224d0/
      DATA pps0 /7.75792d0/
      DATA ppsna / 1.01545d0/
      DATA ppsnb /1.d0/
      DATA ppsnc /5.d0/
c p p pi
c      DATA ppp0 /-2.980d0/
      DATA ppp0 /-3.67510d0/
      DATA pppna / 1.8246d0/
      DATA pppnb /1.d0/
      DATA pppnc /5.d0/
c  
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      call wangma(sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss)
      call wangma(sps0,r0,spsna,spsnb,rt,spsnc,sps,dsps)
      call wangma(pss0,r0,pssna,pssnb,rt,pssnc,pss,dpss)
      call wangma(pps0,r0,ppsna,ppsnb,rt,ppsnc,pps,dpps)
      call wangma(ppp0,r0,pppna,pppnb,rt,pppnc,ppp,dppp)
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxyz('x',sps,dsps,t,dt)       !s-px
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)

      call sxyz('y',sps,dsps,t,dt)       !s-py
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)

      call sxyz('z',sps,dsps,t,dt)       !s-pz
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)

      call sxyz('x',pss,dpss,t,dt)       !px-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)

      call xx('x',pps,dpps,ppp,dppp,t,dt)   !px-px
      hamil(2,2)=t
      dhamil(2,2,1)=dt(1)
      dhamil(2,2,2)=dt(2)
      dhamil(2,2,3)=dt(3)

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !px-py
      hamil(2,3)=t
      dhamil(2,3,1)=dt(1)
      dhamil(2,3,2)=dt(2)
      dhamil(2,3,3)=dt(3)

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !px-pz
      hamil(2,4)=t
      dhamil(2,4,1)=dt(1)
      dhamil(2,4,2)=dt(2)
      dhamil(2,4,3)=dt(3)

      call sxyz('y',pss,dpss,t,dt)          !py-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !py-px
      hamil(3,2)=t
      dhamil(3,2,1)=dt(1)
      dhamil(3,2,2)=dt(2)
      dhamil(3,2,3)=dt(3)

      call xx('y',pps,dpps,ppp,dppp,t,dt)   !py-py
      hamil(3,3)=t
      dhamil(3,3,1)=dt(1)
      dhamil(3,3,2)=dt(2)
      dhamil(3,3,3)=dt(3)

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !py-pz
      hamil(3,4)=t
      dhamil(3,4,1)=dt(1)
      dhamil(3,4,2)=dt(2)
      dhamil(3,4,3)=dt(3)

      call sxyz('z',pss,dpss,t,dt)          !pz-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !pz-px
      hamil(4,2)=t
      dhamil(4,2,1)=dt(1)
      dhamil(4,2,2)=dt(2)
      dhamil(4,2,3)=dt(3)

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !pz-py
      hamil(4,3)=t
      dhamil(4,3,1)=dt(1)
      dhamil(4,3,2)=dt(2)
      dhamil(4,3,3)=dt(3)

      call xx('z',pps,dpps,ppp,dppp,t,dt)   !pz-pz
      hamil(4,4)=t
      dhamil(4,4,1)=dt(1)
      dhamil(4,4,2)=dt(2)
      dhamil(4,4,3)=dt(3)

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between C and C.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine c_c_rp (repul,drepul)
c
      implicit none
      double precision repul,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision ecore,r0,ma,mb,rc,mc,t,dt
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for carbon-carbon
c
      DATA r0/1.312d0/
      DATA rc/1.9d0/
c      DATA ecore/20.356d0/
      DATA ecore/22.68939d0/
      DATA ma/2.72405d0/
      DATA mb/1.d0/
      DATA mc/7.d0/
c
      call wangma(ecore,r0,ma,mb,rc,mc,t,dt)

      repul=t
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine c_h_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss
      double precision pss0,   pssna,pssnb,   pssnc,pss,dpss
c
c Hopping parameters for carbon-hydrogen from Wang ang Mak.
c
      DATA r0 /1.09d0/
      DATA rt /2.d0/
c s s sigma
c      DATA sss0 /-6.921d0/
      DATA sss0 /-6.9986d0/
c      DATA sssna /1.9000/
c      DATA sssnb /1.9000/
      DATA sssna /1.9700/
      DATA sssnb /1.9700/
      DATA sssnc /9.000/
c p s sigma
c      DATA pss0 /-5.594/
      DATA pss0 /-7.39d0/
      DATA pssna /1.6030/
      DATA pssnb /1.6030/
      DATA pssnc /9.000/
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      call wangma(sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss)
      call wangma(pss0,r0,pssna,pssnb,rt,pssnc,pss,dpss)
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxyz('x',pss,dpss,t,dt)       !px-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)

      call sxyz('y',pss,dpss,t,dt)       !py-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)

      call sxyz('z',pss,dpss,t,dt)       !pz-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and Ni.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine c_trhd (ip,hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,jp,i
      integer sN(6),dN(6)
      double precision szetapm(6),pzetapm(6),szeta(6)
      double precision sss,dsss,sds,dsds,pss,dpss,pds,dpds,pdp,dpdp
      double precision coef1(6),coef2(6),dzeta1(6),dzeta2(6)
      double precision r,k(6),r0(6),mu0(6)
      double precision dkllp,kllp,x,param(20)
      double precision tmp1,tmp2,tmp3,tmp4,temp,intgr,dintgr

      common /coulom/ vsip
      common /distan/ r
c
c STO parameters for carbon element
c
      DATA (szetapm(i),i=1,6)/1.56068,1.13260,1.59663,1.52369,
     &                        1.57727,1.25328/
      DATA (pzetapm(i),i=1,6)/1.56068,1.13260,1.59663,1.52369,
     &                        1.57727,1.25328/
c
c STO parameters for transition element
c
      DATA (sN(i),i=1,6) /4,4,5,5,6,6/
      DATA (dN(i),i=1,6) /3,3,4,4,5,5/
      DATA (szeta(i),i=1,6) /1.46952,2.14462,1.75195,1.81756,
     &                       2.04336,2.53186/
      DATA (coef1(i),i=1,6) /0.58579,0.59552,0.54226,0.56019,
     &                       0.65678,0.64695/
      DATA (coef2(i),i=1,6) /0.64856,0.57655,0.65649,0.55617,
     &                       0.57165,0.53788/
      DATA (dzeta1(i),i=1,6) /5.39452,5.89462,5.54495 ,4.90056,
     &                       5.50236,5.89462/
      DATA (dzeta2(i),i=1,6) /1.64452,2.24462,2.17495,2.49556,
     &                        2.18536,2.24462/
c
c Wolfberg-Holmholz constants
c
      DATA (  k(i),i=1,6) /0.26441,0.18125,0.29889,0.20342,
     &                     0.34261,0.30277/
      DATA ( r0(i),i=1,6) /1.9928 ,1.2870 ,2.9265,1.5743,
     &                     2.0176 ,2.9275/
      DATA (mu0(i),i=1,6) /0.01756,1.27014,0.04703,1.37253,
     &                     0.01994,0.67298/
c
      jp=ip-2
c
      temp=k(jp)*(r/r0(jp))**mu0(jp)
      kllp=k(jp)+temp
      dkllp=mu0(jp)*temp/r
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c Actually, the final result is the product of overlap integral and constant
c K (distance dependent).
c
c s s sigma bonding overlap integral in z-axis
c
      call ovlcon(2,0,szetapm(jp),sN(jp),0,0,szeta(jp),intgr,dintgr)
      dsss=dintgr*kllp+intgr*dkllp
       sss= intgr*kllp
c
c p s sigma bonding overlap integral in z-axis
c
      call ovlcon(2,1,pzetapm(jp),sN(jp),0,0,szeta(jp),intgr,dintgr)
      dpss=dintgr*kllp+intgr*dkllp
       pss=intgr*kllp
c
c s d sigma bonding overlap integral in z-axis
c
      call ovlcon(2,0,szetapm(jp),dN(jp),2,0,dzeta1(jp),intgr,dintgr)
       sds= intgr*coef1(jp)
      dsds=dintgr*coef1(jp)

      call ovlcon(2,0,szetapm(jp),dN(jp),2,0,dzeta2(jp),intgr,dintgr)
       sds= intgr*coef2(jp)+ sds
      dsds=dintgr*coef2(jp)+dsds
c
      dsds=dsds*kllp+sds*dkllp
       sds= sds*kllp
c
c p d sigma bonding overlap integral in z-axis
c
      call ovlcon(2,1,pzetapm(jp),dN(jp),2,0,dzeta1(jp),intgr,dintgr)
       pds= intgr*coef1(jp)
      dpds=dintgr*coef1(jp)

      call ovlcon(2,1,pzetapm(jp),dN(jp),2,0,dzeta2(jp),intgr,dintgr)
       pds= intgr*coef2(jp)+ pds
      dpds=dintgr*coef2(jp)+dpds
c
      dpds=dpds*kllp+pds*dkllp
       pds= pds*kllp
c
c p d pi bonding overlap integral in z-axis
c
      call ovlcon(2,1,pzetapm(jp),dN(jp),2,1,dzeta1(jp),intgr,dintgr)
       pdp= intgr*coef1(jp)
      dpdp=dintgr*coef1(jp)

      call ovlcon(2,1,pzetapm(jp),dN(jp),2,1,dzeta2(jp),intgr,dintgr)
       pdp= intgr*coef2(jp)+ pdp
      dpdp=dintgr*coef2(jp)+dpdp
c
      dpdp=dpdp*kllp+pdp*dkllp
       pdp= pdp*kllp
c
c Multiply kmn with (Hm(Ni)+Hn(C))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(2,1)+vsip(ip,1))/2.0
      tmp2=(vsip(2,1)+vsip(ip,2))/2.0
      tmp3=(vsip(2,2)+vsip(ip,1))/2.0
      tmp4=(vsip(2,2)+vsip(ip,2))/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t*tmp1
      dhamil(1,1,1)=dt(1)*tmp1
      dhamil(1,1,2)=dt(2)*tmp1
      dhamil(1,1,3)=dt(3)*tmp1

      call sxy('xy',sds,dsds,t,dt)       !s-dxy
      hamil(1,2)=t*tmp2
      dhamil(1,2,1)=dt(1)*tmp2
      dhamil(1,2,2)=dt(2)*tmp2
      dhamil(1,2,3)=dt(3)*tmp2

      call sxy('yz',sds,dsds,t,dt)       !s-dyz
      hamil(1,3)=t*tmp2
      dhamil(1,3,1)=dt(1)*tmp2
      dhamil(1,3,2)=dt(2)*tmp2
      dhamil(1,3,3)=dt(3)*tmp2

      call sxy('zx',sds,dsds,t,dt)       !s-dzx
      hamil(1,4)=t*tmp2
      dhamil(1,4,1)=dt(1)*tmp2
      dhamil(1,4,2)=dt(2)*tmp2
      dhamil(1,4,3)=dt(3)*tmp2

      call sxxyy(sds,dsds,t,dt)          !s-dx2_y2
      hamil(1,5)=t*tmp2
      dhamil(1,5,1)=dt(1)*tmp2
      dhamil(1,5,2)=dt(2)*tmp2
      dhamil(1,5,3)=dt(3)*tmp2

      call szz(sds,dsds,t,dt)            !s-dz2
      hamil(1,6)=t*tmp2
      dhamil(1,6,1)=dt(1)*tmp2
      dhamil(1,6,2)=dt(2)*tmp2
      dhamil(1,6,3)=dt(3)*tmp2

      call sxyz('x',pss,dpss,t,dt)       !px-s 
      hamil(2,1)=t*tmp3
      dhamil(2,1,1)=dt(1)*tmp3
      dhamil(2,1,2)=dt(2)*tmp3
      dhamil(2,1,3)=dt(3)*tmp3

      call xxy('xy',pds,dpds,pdp,dpdp,t,dt)         !px-dxy
      hamil(2,2)=t*tmp4
      dhamil(2,2,1)=dt(1)*tmp4
      dhamil(2,2,2)=dt(2)*tmp4
      dhamil(2,2,3)=dt(3)*tmp4

      call xyz(pds,dpds,pdp,dpdp,t,dt)              !px-dyz
      hamil(2,3)=t*tmp4
      dhamil(2,3,1)=dt(1)*tmp4
      dhamil(2,3,2)=dt(2)*tmp4
      dhamil(2,3,3)=dt(3)*tmp4

      call xxy('xz',pds,dpds,pdp,dpdp,t,dt)         !px-dzx
      hamil(2,4)=t*tmp4
      dhamil(2,4,1)=dt(1)*tmp4
      dhamil(2,4,2)=dt(2)*tmp4
      dhamil(2,4,3)=dt(3)*tmp4

      call xxxyy(pds,dpds,pdp,dpdp,t,dt)            !px-dx2_y2
      hamil(2,5)=t*tmp4
      dhamil(2,5,1)=dt(1)*tmp4
      dhamil(2,5,2)=dt(2)*tmp4
      dhamil(2,5,3)=dt(3)*tmp4

      call xzz('x',pds,dpds,pdp,dpdp,t,dt)              !px-dz2
      hamil(2,6)=t*tmp4
      dhamil(2,6,1)=dt(1)*tmp4
      dhamil(2,6,2)=dt(2)*tmp4
      dhamil(2,6,3)=dt(3)*tmp4

      call sxyz('y',pss,dpss,t,dt)       !py-s 
      hamil(3,1)=t*tmp3
      dhamil(3,1,1)=dt(1)*tmp3
      dhamil(3,1,2)=dt(2)*tmp3
      dhamil(3,1,3)=dt(3)*tmp3

      call xxy('yx',pds,dpds,pdp,dpdp,t,dt)         !py-dxy
      hamil(3,2)=t*tmp4
      dhamil(3,2,1)=dt(1)*tmp4
      dhamil(3,2,2)=dt(2)*tmp4
      dhamil(3,2,3)=dt(3)*tmp4

      call xxy('yz',pds,dpds,pdp,dpdp,t,dt)         !py-dyz
      hamil(3,3)=t*tmp4
      dhamil(3,3,1)=dt(1)*tmp4
      dhamil(3,3,2)=dt(2)*tmp4
      dhamil(3,3,3)=dt(3)*tmp4

      call xyz(pds,dpds,pdp,dpdp,t,dt)              !py-dzx
      hamil(3,4)=t*tmp4
      dhamil(3,4,1)=dt(1)*tmp4
      dhamil(3,4,2)=dt(2)*tmp4
      dhamil(3,4,3)=dt(3)*tmp4

      call yxxyy(pds,dpds,pdp,dpdp,t,dt)            !py-dx2_y2
      hamil(3,5)=t*tmp4
      dhamil(3,5,1)=dt(1)*tmp4
      dhamil(3,5,2)=dt(2)*tmp4
      dhamil(3,5,3)=dt(3)*tmp4

      call xzz('y',pds,dpds,pdp,dpdp,t,dt)              !py-dz2
      hamil(3,6)=t*tmp4
      dhamil(3,6,1)=dt(1)*tmp4
      dhamil(3,6,2)=dt(2)*tmp4
      dhamil(3,6,3)=dt(3)*tmp4

      call sxyz('z',pss,dpss,t,dt)       !pz-s 
      hamil(4,1)=t*tmp3
      dhamil(4,1,1)=dt(1)*tmp3
      dhamil(4,1,2)=dt(2)*tmp3
      dhamil(4,1,3)=dt(3)*tmp3

      call xyz(pds,dpds,pdp,dpdp,t,dt)               !pz-dxy
      hamil(4,2)=t*tmp4
      dhamil(4,2,1)=dt(1)*tmp4
      dhamil(4,2,2)=dt(2)*tmp4
      dhamil(4,2,3)=dt(3)*tmp4

      call xxy('zy',pds,dpds,pdp,dpdp,t,dt)          !pz-dyz
      hamil(4,3)=t*tmp4
      dhamil(4,3,1)=dt(1)*tmp4
      dhamil(4,3,2)=dt(2)*tmp4
      dhamil(4,3,3)=dt(3)*tmp4

      call xxy('zx',pds,dpds,pdp,dpdp,t,dt)          !pz-dzx
      hamil(4,4)=t*tmp4
      dhamil(4,4,1)=dt(1)*tmp4
      dhamil(4,4,2)=dt(2)*tmp4
      dhamil(4,4,3)=dt(3)*tmp4

      call zxxyy(pds,dpds,pdp,dpdp,t,dt)             !pz-dx2_y2
      hamil(4,5)=t*tmp4
      dhamil(4,5,1)=dt(1)*tmp4
      dhamil(4,5,2)=dt(2)*tmp4
      dhamil(4,5,3)=dt(3)*tmp4

      call zzz(pds,dpds,pdp,dpdp,t,dt)               !pz-dz2
      hamil(4,6)=t*tmp4
      dhamil(4,6,1)=dt(1)*tmp4
      dhamil(4,6,2)=dt(2)*tmp4
      dhamil(4,6,3)=dt(3)*tmp4

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between C and H, or H and C.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine chhcrp (repul,drepul)
c
      implicit none
      double precision repul,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision ecore,r0,ma,mb,rc,mc,t,dt
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for carbon-hydrogen
c
      DATA rc/1.9d0/
      DATA r0/1.09d0/
c      DATA ecore/9.086d0/
      DATA ecore/10.8647d0/
      DATA ma/3.1d0/
      DATA mb/3.1d0/
      DATA mc/10.0d0/
 
      call wangma(ecore,r0,ma,mb,rc,mc,t,dt)

      repul=t
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Ni and C or C and Ni.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine ctrcrp (ip,repul,drepul)
c
      implicit none
      integer i,ip,jp
      double precision repul,dt,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn,r
      double precision alphap(6),betap(6),gammap(6)
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for transition metal-carbon
c
      DATA (alphap(i),i=1,6) /1846.09,3760.37,2884.69,2145.03,
     &                        2725.73,3413.9 /
      DATA (betap(i), i=1,6) /4.3979 ,5.4333 ,4.2078 ,4.3155 ,
     &                        4.5494,4.5461/
      DATA (gammap(i),i=1,6) /0.28764,1.20464,0.20619,0.73337,
     &                        1.07947 ,0.79440/
c
      jp=ip-2
c
      call grepfm (alphap(jp),gammap(jp),betap(jp),repul,dt)
      
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
C The program is coded for Equation (8) in reference:
C   H.W. Jones, Int. J. Quantum Chem. 18, 709 (1980).
C The program basically can work for any combinations of the quantum numbers,
C but the C matrix coefficientd provided in the data block preclude that kind
C of job with the principal quantum number bigger than 4.

      subroutine dovlap(Nprime,Lprime,zetapm,N,L,M,zeta,S,dS)

      implicit none
      double precision autoang,autoev
      parameter (autoang=0.529177249d0)
      parameter (autoev = 27.2113961d0)

      integer Nprime, Lprime,N,L,M,ns,k,i,j,tmp              !ns lower n
      double precision fact(0:12)
      double precision a,S,Keba,cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)
      double precision zeta,zetapm,tmpa,tmpb,ratioz,azeta,azetap
      double precision expzet,expzpm,azeadd,azminu,Kbrati
      double precision dS,dtempa,dtempb,r

      common /factorial/ fact
      common /cmatrix/ cmatrx
      common /distan/ r

      ratioz=(zetapm/zeta)**(Nprime+0.50)

      a=r/autoang

      if (r.gt.10.0) then
        S=0.0
        dS=0.0
        return
      endif

      azetap=a*zetapm
      azeta=a*zeta

      expzet=dexp(-azeta)
      expzpm=dexp(-azetap)

      azeadd=azetap+azeta
      azminu=azetap-azeta

      Keba=2**(Nprime+N)*(-1)**M

      tmpa=(2*L+1)*(2*Lprime+1)
      tmp=L+M
      tmpa=tmpa*fact(tmp)
      tmp=Lprime-M
      tmpa=tmpa*fact(tmp)
      tmp=2*Nprime
      tmpa=tmpa/fact(tmp)
      tmp=2*N
      tmpa=tmpa/fact(tmp)
      tmp=Lprime+M
      tmpa=tmpa/fact(tmp)
      tmp=L-M
      tmpa=tmpa/fact(tmp)

      keba=Keba*dsqrt(tmpa)

      S=0.0
      dS=0.0

      do 1 i=0,(N+L+Lprime)
        do 2 j=0,(N+Lprime)
          ns=Nprime-Lprime+j
          tmpa=0.0
          dtempa=0.0
          do 3 k=0,ns
            tmpb=(-1)**j/azminu**(k+1)
            tmpb=(-1)**i/azeadd**(k+1)-tmpb 
            tmp=ns-k
            tmpb=tmpb/fact(tmp)
C Prepare to compute the first derivative
            dtempb=-(k+1)*tmpb/a
            tmpa=tmpa+tmpb
            dtempa=dtempa+dtempb
 3        continue   
          tmpa=tmpa*expzpm
          dtempa=dtempa*expzpm-zetapm*tmpa
C
          tmpb=1.0/azeadd**(ns+1)
          tmpb=(-1)**j/azminu**(ns+1)-tmpb
          dtempb=-(ns+1)*tmpb/a
          tmpb=tmpb*expzet
          tmpa=tmpa+tmpb
          dtempb=dtempb*expzet-zeta*tmpb
          dtempa=dtempa+dtempb

          tmp=Nprime-2*Lprime-L+i+j
          tmpb=fact(ns)*cmatrx(Lprime,N,L,M,i,j)*azeta**tmp
          dtempb=tmp*tmpb/a

          S=S+tmpa*tmpb
          dS=dS+tmpa*dtempb+dtempa*tmpb

 2      continue
 1    continue
 
      Kbrati=Keba*ratioz

      S=S*Keba*ratioz*(-1)**M
      dS=dS*Keba*ratioz*(-1)**M
c
c Convert from bohr -1 to A -1
c
      dS=dS/autoang
c     if (M.eq.2) then
c     endif

      return
      end
c
c Arrange all the valence electrons to the necessary molecular orbitals 
c so that the potential energy of the system is minimized.
c
      subroutine filloc(ndim,nElect,eva,etb,occ,natom)
c
      implicit none
c
      integer ndim,nElect,i,j,nbottom,ntop,k,natom
      double precision eva(natom*6),u,etb,tmp,upenalty
      integer occ(natom*6)
c
c nElect --- number of valence electrons in the system.
c ndim --- dimension of the Hamiltonian matrix.
c
c AJ Penalty energy
c      u=0.64d0 ! used in Tiqing's parameterizations for H, C, and transition metals
c      u=0.07d0 ! used by Grazyna's parameteritaztions for Al
c      u = 0.d0 ! used in newer Al and Al&H parameterizations by Jasper

      common /penalty/ upenalty

      u = upenalty

c
c Initialize occupation number set
c
      do 10 i=1,ndim
         occ(i)=0
   10 continue
c
      if (nElect .gt. 2*ndim) stop "There is not enough MOs"
c
c Arrange electrons simply from lowest orbitals to highest orbitals
c
      nbottom=nElect/2
      if (nElect .eq. 1) then
        occ(1)=1
        etb=eva(1)
        return
      else
        do 20 i=1, nElect/2
          occ(i)=2
   20   continue
        if (2*(nElect/2) .ne. nElect) then
          occ(i)=1
          ntop=nbottom+2
        else
          ntop=nbottom+1
        endif
      endif
c
c Rearrange the occupation according to penalty energy to obtain a set of 
c occupancy to give minimum energy.
c At first, the orbitals are occupied from lowest to higher orbital with two
c electrons if there are enough electrons, which is called the starting point. 
c The rearrangement of the the molecular orbitals begins here.  The highest 
c doubly occupied orbital is called nbottom, and the lowest unoccupied MO is
c called ntop.  If the energy difference between ntop and nbottom is smaller
c than the penalty energy (u or 0.80 eV in this case), the electron from 
c nbottom jump to ntop.  Otherwise, the rearrangement will stop. 
c
      if (nbottom .eq. 0) goto 102
      if (ntop .gt. ndim) goto 102
c
 101  if (eva(nbottom)+u .gt. eva(ntop)) then 
c
c If the energy of the top orbital is higher by u than that of the bottom 
c orbital, the bottom becomes top and the number of the bottom orbital is 
c decreased by 1.
c
        occ(ntop)=1
        ntop=ntop+1
        occ(nbottom)=1
        nbottom=nbottom-1
c
        if (ntop .gt. ndim) goto 102
        if (nbottom .eq. 0) goto 102
c
        goto 101
c
      endif
c
c Calculate the energy of resulting structure.
c
 102  continue
      etb=0.0
      do k=1,ndim
        if (occ(k).eq.2) then
          etb=etb+u+eva(k)*2
        else
          etb=etb+eva(k)*occ(k)
        endif
      enddo
c
      return
      end
      subroutine genrhd(i,j,hamil,dhamil)

      implicit none
      integer i,j
      double precision hamil(6,6),dhamil(6,6,3)
c
c Hydrogen -- other elements.
c
      if (i.eq.1) then
        if (j.eq.1) then
          call h_h_hd(hamil,dhamil)
        else if (j.eq.2) then
          call h_c_hd(hamil,dhamil)
        else if (j.eq.9) then
          call h_al_hd(hamil,dhamil)
        else 
          call h_trhd(j,hamil,dhamil)
        endif
c
c Carbon -- other elements.
c
      else if (i.eq.2) then
        if (j.eq.1) then
          call c_h_hd(hamil,dhamil)
        else if (j.eq.2) then
          call c_c_hd(hamil,dhamil)
        else if (j.eq.9) then
          write(6,*)"Cant do Al and C together!"
          stop
        else
          call c_trhd(j,hamil,dhamil)
        endif
c
c Aluminum -- other elements.
c
      else if (i.eq.9) then
        if (j.eq.9) then
          call al_al_hd(hamil,dhamil)
        else if (j.eq.1) then
          call al_h_hd(hamil,dhamil)
        else
          stop "The parameters for the two elements are unavailable."
        endif
c
c Ni, Cu, Pd, Ag, Pt, and Au --- C, H, and itself.
c
      else 
        if (j.eq.1) then
          call trh_hd(i,hamil,dhamil)
        else if (j.eq.2) then
          call trc_hd(i,hamil,dhamil)
        else
          if (i.eq.j) then
            call trtrhd(i,hamil,dhamil)
          else
            stop "The parameters for the two elements are unavailable."
          endif
        endif
      endif

      return
      end
      subroutine genrrp(i,j,repul,drepul)

      implicit none
      integer i,j
      double precision repul,drepul(3)
c
c Repulsion: Hydrogen --- other elements.
c
      if (i.eq.1) then
        if (j.eq.1) then
          call h_h_rp(repul,drepul)
        else if (j.eq.2) then
          call chhcrp(repul,drepul)
        else if (j.eq.9) then
          call h_al_rp(repul,drepul)
        else 
          call htrhrp(j,repul,drepul)
        endif
c
c Repulsion: Carbon   --- other elements.
c
      else if (i.eq.2) then
        if (j.eq.1) then
          call chhcrp(repul,drepul)
        else if (j.eq.2) then
          call c_c_rp(repul,drepul)
        else if (j.eq.9) then
          write(6,*)"Cant do Al and C"
          stop
        else
          call ctrcrp(j,repul,drepul)
        endif
c
c
c Repulsion: Aluminum   --- other elements.
c
      else if (i.eq.9) then
        if (j.eq.9) then
          call al_al_rp(repul,drepul)
        else if (j.eq.1) then
          call h_al_rp(repul,drepul)
        else
          stop "The parameters for the two elements are unavailable."
        endif
c
c Repulsion: Ni, Cu, Pd, Ag, Pt, and Au  --- C, H, and itself.
c
      else
        if (j.eq.1) then
          call htrhrp(i,repul,drepul)
        else if (j.eq.2) then
          call ctrcrp(i,repul,drepul)
        else
          if (i.eq.j) then
            call trtrrp(i,repul,drepul)
          else
            stop "The parameters for the two elements are unavailable."
          endif
        endif
      endif

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c The following function is called by tr_c, tr_h, and tr-tr.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine grepfm (alpha,gamma,beta,t,dt)
c
      implicit none
      double precision alpha,gamma,beta,t,dt,r,tempa

      common /distan/ r

      tempa=alpha*exp(-beta*r)

      t=(1.0/r+gamma)*tempa

      dt=-tempa/r/r-beta*t

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between H and Al.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_al_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,i
      double precision zetas,zetap,zetasp
      double precision sss,dsss,sps,dsps,pps,dpps,ppp,dppp,pss,dpss
      double precision r,kwh
      double precision x,kss,ksp,kpp
      double precision tmp1,tmp2,tmp3,tmp4,temp,intgr,dintgr
      double precision szetapm

      common /coulom/ vsip
      common /distan/ r
c
c STP parameters for H
c
      DATA szetapm /1.2/
c
c STO parameters for Al
c
      DATA zetas/1.3724/ ! Al s exponent in overlap parameter in au
      DATA zetap/1.3552/ ! Al p exponent in overlap parameter in au
      zetasp = (zetas+zetap)/2.d0
c
c Wolfberg-Holmholz constant for Al&H
c
      DATA kwh /0.39280/
c
c Wolfsberg-Helmholtz relation for Al(3s)/Al(3p) and H(1s) VSIPs.
c
      tmp1=(vsip(9,1)+vsip(1,1))/2.0
      tmp2=(vsip(9,2)+vsip(1,1))/2.0

      kss = kwh*tmp1
      ksp = kwh*tmp2
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
      call ovlcon(1,0,szetapm,3,0,0,zetas,intgr,dintgr)
      dsss=dintgr*kss
       sss= intgr*kss

c s p sigma bonding overlap integral in z-axis
      call ovlcon(1,0,szetapm,3,1,0,zetasp,intgr,dintgr)
      dsps=dintgr*ksp
       sps= intgr*ksp

c Project the two-center integral along the basis set.
c

      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxyz('x',sps,dsps,t,dt)       !s-px
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)

      call sxyz('y',sps,dsps,t,dt)       !s-py
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)

      call sxyz('z',sps,dsps,t,dt)       !s-pz
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Al and Al.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_al_rp (repul,drepul)
c
      implicit none
      double precision repul,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision aa,bb,mu,r,dt
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /distan/ r
c
c Repulsive parameters for carbon-carbon
c
      DATA aa/1142.d0/
      DATA bb/2.921d0/
      DATA mu/0.02936d0/

      repul=aa/(r**mu)*dexp(-bb*r)
      dt=-bb*aa/(r**mu)*dexp(-bb*r)+aa/(r**(mu-1.d0))*dexp(-bb*r)*(-mu)
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_c_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss
      double precision sps0,   spsna,spsnb,   spsnc,sps,dsps
c
c Hopping parameters for hydrogen-carbon from Wang and Mak.
c
      DATA r0 /1.09d0/
      DATA rt /2.d0/
c s s sigma
c      DATA sss0 /-6.921d0/
      DATA sss0 /-6.9986d0/
c      DATA sssna /1.9000/
c      DATA sssnb /1.9000/
      DATA sssna /1.970/
      DATA sssnb /1.970/
      DATA sssnc /9.000/
c s p sigma
c      DATA sps0 /5.594d0/
      DATA sps0 /7.39d0/
      DATA spsna /1.603d0/
      DATA spsnb /1.603d0/
      DATA spsnc /9.0/
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      call wangma(sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss)
      call wangma(sps0,r0,spsna,spsnb,rt,spsnc,sps,dsps)
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxyz('x',sps,dsps,t,dt)       !s-px
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)

      call sxyz('y',sps,dsps,t,dt)       !s-py
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)

      call sxyz('z',sps,dsps,t,dt)       !s-pz
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between H and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_h_hd (hamil,dhamil)
c
      implicit none
      integer j
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss
c
c Hopping parameters for hydrogen-hydrogen
c
      DATA r0 /0.7414d0/
      DATA rt /4.176d0/
      DATA sss0 /-4.842d0/
      DATA sssna /2.049D0/
      DATA sssnb /1.762D0/
      DATA sssnc /5.609d0/
c
      call wangma(sss0,r0,sssna,sssnb,rt,sssnc,sss,dsss)
c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hamil(1,1)=0.d0
      do j=1,3
      dhamil(1,1,j)=0.d0
      enddo

      return

      call ss(sss,dsss,t,dt)     !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between H and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_h_rp (repul,drepul)
c
      implicit none
      double precision repul,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision ecore,r0,ma,mb,rc,mc,t,dt
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for hydrogen-hydrogen
c
      DATA r0/0.7414d0/
      DATA rc/1.162d0/
      DATA ecore/4.279D0/
      DATA ma/2.436/
      DATA mb/3.309/
      DATA mc/3.653/
c
      call wangma(ecore,r0,ma,mb,rc,mc,t,dt)

      repul=t
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n

      repul = 0.d0
      drepul(1)=0.d0
      drepul(2)=0.d0
      drepul(3)=0.d0
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between H and Ni.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine h_trhd (ip,hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,jp,i,j
      integer sN(6),dN(6)
      double precision szetapm(6),szeta(6)
      double precision sss,dsss,sds,dsds
      double precision kss(6),ksd(6)
      double precision tmp1,tmp2,intgr,dintgr
      double precision coef1(6),coef2(6),dzeta1(6),dzeta2(6)

      common /coulom/ vsip
c
c STO parameters for hydrogen element
c
      DATA (szetapm(i),i=1,6) /1.22576,1.26397,1.21011,1.23966,
     &                         1.25792,1.12836/
c
c STO parameters for transition elements
c
      DATA (sN(i),i=1,6) /4,4,5,5,6,6/
      DATA (dN(i),i=1,6) /3,3,4,4,5,5/
      DATA (szeta(i),i=1,6) /1.76329,1.81529,2.03423,1.74711,
     &                       2.36819,2.12506/
      DATA (coef1(i),i=1,6) /0.57107,0.61029,0.53160,0.56239,
     &                       0.64119,0.66519/
      DATA (coef2(i),i=1,6) /0.63227,0.59085,0.64359,0.55836,
     &                       0.55808,0.55304/
      DATA (dzeta1(i),i=1,6) /5.68829,5.56529,5.82723,4.83011,
     &                       5.82719,5.68606/
      DATA (dzeta2(i),i=1,6) /1.93829,1.91529,2.45723,2.42511,
     &                       2.51109,2.31706/
c
c Wolfberg-Holmholz constants
c
      DATA (kss(i),i=1,6) /0.75730 ,0.80843,0.67048,0.57090,
     &                     0.91655 ,0.61968/
      DATA (ksd(i),i=1,6) /0.18030 ,0.02514,0.03081,0.05065,
     &                     0.31296 ,0.003524/
c
      jp=ip-2
c     AJ:  jp = 1-6 for Ni, Cu, Pd, Ag, Pt, Au?
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
c
      call ovlcon(1,0,szetapm(jp),sN(jp),0,0,szeta(jp),intgr,dintgr)
       sss=intgr
      dsss=dintgr

c s d sigma bonding overlap integral in z-axis
c        
      call ovlcon(1,0,szetapm(jp),dN(jp),2,0,dzeta1(jp),intgr,dintgr)
       sds= intgr*coef1(jp)
      dsds=dintgr*coef1(jp)

      call ovlcon(1,0,szetapm(jp),dN(jp),2,0,dzeta2(jp),intgr,dintgr)
       sds= intgr*coef2(jp)+ sds
      dsds=dintgr*coef2(jp)+dsds
c
c Multiply kmn with (Hm(Ni)+Hn(H))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(ip,1)+vsip(1,1))*kss(jp)/2.0
      tmp2=(vsip(ip,2)+vsip(1,1))*ksd(jp)/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t*tmp1
      dhamil(1,1,1)=dt(1)*tmp1
      dhamil(1,1,2)=dt(2)*tmp1
      dhamil(1,1,3)=dt(3)*tmp1

      call sxy('xy',sds,dsds,t,dt)       !s-dxy
      hamil(1,2)=t*tmp2
      dhamil(1,2,1)=dt(1)*tmp2
      dhamil(1,2,2)=dt(2)*tmp2
      dhamil(1,2,3)=dt(3)*tmp2

      call sxy('yz',sds,dsds,t,dt)       !s-dyz
      hamil(1,3)=t*tmp2
      dhamil(1,3,1)=dt(1)*tmp2
      dhamil(1,3,2)=dt(2)*tmp2
      dhamil(1,3,3)=dt(3)*tmp2

      call sxy('zx',sds,dsds,t,dt)       !s-dzx
      hamil(1,4)=t*tmp2
      dhamil(1,4,1)=dt(1)*tmp2
      dhamil(1,4,2)=dt(2)*tmp2
      dhamil(1,4,3)=dt(3)*tmp2

      call sxxyy(sds,dsds,t,dt)          !s-dx2_y2
      hamil(1,5)=t*tmp2
      dhamil(1,5,1)=dt(1)*tmp2
      dhamil(1,5,2)=dt(2)*tmp2
      dhamil(1,5,3)=dt(3)*tmp2

      call szz(sds,dsds,t,dt)            !s-dz2
      hamil(1,6)=t*tmp2
      dhamil(1,6,1)=dt(1)*tmp2
      dhamil(1,6,2)=dt(2)*tmp2
      dhamil(1,6,3)=dt(3)*tmp2

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Ni and H or H and Ni.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine htrhrp (ip,repul,drepul)
c
      implicit none
      integer i,ip,jp
      double precision repul,dt,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision alpha(6),gamma(6),beta(6)
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for nickel-hydrogen
c
      DATA (alpha(i),i=1,6) /161.153,283.417,300.257,570.769,
     &                       255.462,514.39/
      DATA (gamma(i),i=1,6) /0.26570,0.11325,0.18122,1.02695,
     &                       0.30071,0.50950/
      DATA (beta(i), i=1,6) /3.3225 ,3.6555 ,3.5970,4.5567,
     &                       3.3099 ,4.3302/
c
      jp=ip-2
c
      call grepfm (alpha(jp),gamma(jp),beta(jp),repul,dt)
      
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c  This formula is taken from:
c  N.N. Lathiotakis, A.N. Andriotis, M. Menon, and J. Connolly, J. Chem. 
c  Phys. 104, 992-1003 (1996).
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine lathio (eta,rd,d,tau,alpha,V,dV)
c
      implicit none
      double precision eta,rd,d,tau,alpha,V,dV,r

      common /distan/ r

      V=eta*7.62d0*rd**tau/d**(tau+2.0)*exp(-alpha*(r-d))
      dV=-alpha*V

      return
      end
c This fragment calculate directional cosines between atom i and atom j in 
c the Hamiltonian matrix.

      subroutine lmndlm(i,j,ratom,natom)

      implicit none
      integer i,j,natom
      double precision ratom(3*natom)
      double precision dx,dy,dz,r
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      common /distan/ r

      dx=ratom(3*(i-1)+1)-ratom(3*(j-1)+1)
      dy=ratom(3*(i-1)+2)-ratom(3*(j-1)+2)
      dz=ratom(3*(i-1)+3)-ratom(3*(j-1)+3)
      
      r=dsqrt(dx*dx+dy*dy+dz*dz)

      l=dx/r
      m=dy/r
      n=dz/r

      lladmm=l*l+m*m
      llmimm=l*l-m*m
      nn=n*n

      dlx=(1.0-l*l)/r
      dly=-l*m/r
      dlz=-l*n/r
      dmx=dly
      dmy=(1.0-m*m)/r
      dmz=-m*n/r
      dnx=dlz
      dny=dmz
      dnz=(1.0-n*n)/r

      return
      end
C The program is coded for Equation (8) in reference:
C   H.W. Jones, Int. J. Quantum Chem. 19, 567 (1981).
C The program basically can work for any combinations of the quantum numbers,
C but the C matrix coefficients provided in the data block preclude that kind
C of job with the principal quantum number bigger than 4.
C For the same zeta.

      subroutine ovlapd(CN,CL,zeta,n,l,M,S,dS)

      implicit none

      double precision autoang,autoev
      parameter (autoang=0.529177249d0)
      parameter (autoev = 27.2113961d0)
      integer CN,CL,n,l,M,i,j,t,u,v,k,w,tmp
      double precision fact(0:12)
      double precision zeta,S,dS,r,a,p,Keba
      double precision cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)
      double precision tmpa,tmpb,tmpc,tmp1,tmp2,tmp3

      double precision AJtmp(0:10)

      common /factorial/ fact
      common /cmatrix/ cmatrx
      common /distan/ r

c convert a to au, zeta is in au, p is therefore unitless
      a=r/autoang
      p=zeta*a

      if (r.gt.10.0) then
        S=0.0
        dS=0.0
        return
      endif

      Keba=2**(CN+n)*(-1)**M 

      tmpa=(2*CL+1)*(2*l+1)
      tmp=l+M
      tmpa=tmpa*fact(tmp)
      tmp=CL-M
      tmpa=tmpa*fact(tmp)
      tmp=2*CN
      tmpa=tmpa/fact(tmp)
      tmp=2*n
      tmpa=tmpa/fact(tmp)
      tmp=CL+M
      tmpa=tmpa/fact(tmp)
      tmp=l-M
      tmpa=tmpa/fact(tmp)

      Keba=Keba*sqrt(tmpa)
c      print *,"keba",keba

      S=0.0
      dS=0.0

c AJ
      do i=0,10
        AJtmp(i) = 0.d0
      enddo

      do 1 i=0,n+CL+l
        do 2 j=0,n+CL
          t=j+CN-CL
          u=j+CN-2*CL-l+i
          v=i-CL-l-1
c AJ
          tmp1 = cmatrx(CL,n,l,M,i,j)*(-1)**j/(t+1)
          tmp2 = cmatrx(CL,n,l,M,i,j)*(-fact(t)/2**(t+1))
          if (u.ge.0) AJtmp(u)=AJtmp(u) + tmp1
          if (v.ge.0) AJtmp(v)=AJtmp(v) + tmp2
          S=S+cmatrx(CL,n,l,M,i,j)*
     &       ((-1)**j*p**u/(t+1)-fact(t)*p**v/2**(t+1))
          dS=dS+cmatrx(CL,n,l,M,i,j)*zeta*
     &       (u*(-1)**j*p**(u-1)/(t+1)-v*fact(t)*p**(v-1)/2**(t+1))
          tmpa=0.0
          tmpb=0.0
          do 3 k=0,t
            tmp=t-k
            w=t-k+i-l-CL-1
c           print 99,t,u,v,w
          tmp3 = fact(t)/fact(tmp)/2**(k+1)
     &        *cmatrx(CL,n,l,M,i,j)*(-1)**i
          if (w.ge.0) AJtmp(w)=AJtmp(w) + tmp3

c99    format("t,u,v,and w ",4(I2,1x))
            tmpa=tmpa+fact(t)/fact(tmp)*p**w/2**(k+1)
            tmpb=tmpb+w*fact(t)*zeta/fact(tmp)*p**(w-1)/2**(k+1)
 3        continue
c         print *,"tmpa tmpb ",tmpa,tmpb
          S=S+cmatrx(CL,n,l,M,i,j)*(-1)**i*tmpa
          dS=dS+cmatrx(CL,n,l,M,i,j)*(-1)**i*tmpb
 2      continue
 1    continue

      dS=Keba*(exp(-p)*dS-zeta*exp(-p)*S)*(-1)**M
      S=Keba*exp(-p)*S*(-1)**M
c      print *,"AJ w ep",S/exp(-p),exp(-p),S,p
      do i=0,10
c      print *,i,Keba*AJtmp(i),Keba*AJtmp(i)*p**i
      enddo

c
c Convert from bohr -1 to A -1
c
      dS=dS/autoang

      return
      end
C This subroutine is an interface to call two subroutines to compute the 
C overlap integrals (one has the same zeta values, the other has different
C zeta values.
      subroutine ovlcon(Np,Lp,zetap,N,L,M,zeta,intgr,dintgr)
c
      implicit none
c
      integer Np,Lp,N,L,M
c
      double precision zetap,zeta,intgr,dintgr
c
      if (Np.le.Lp .or. N.le.L) then
         intgr=0.0 
        dintgr=0.0 
        return
      endif
c
      if (abs(zetap-zeta) .lt. 0.001) then
        zeta=(zeta+zetap)/2.0
c
        call ovlapd(Np,Lp,zeta,N,L,M,intgr,dintgr)
      else
c
        call dovlap(Np,Lp,zetap,N,L,M,zeta,intgr,dintgr)
      endif
c
      return
      end
      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c      subroutine pot(symb,x,y,z,v,natom,maxatom)
      implicit none
      integer natom,maxatom,i
      character*2 symb(maxatom)
      double precision x(maxatom),y(maxatom),z(maxatom),
     & dvdx(maxatom),dvdy(maxatom),dvdz(maxatom),v,
     & dxx(maxatom*3),xx(maxatom*3)

      do i=1,natom
        xx(3*(i-1)+1) = x(i)
        xx(3*(i-1)+2) = y(i)
        xx(3*(i-1)+3) = z(i)
      enddo

c      call surf(symb,v,xx,dxx,natom,natom)
      call surf(symb,v,xx,dxx,natom)

      do i=1,natom
        dvdx(i) = dxx(3*(i-1)+1)
        dvdy(i) = dxx(3*(i-1)+2)
        dvdz(i) = dxx(3*(i-1)+3)
      enddo

      return
      end

      subroutine prepot
      end

c
c This is a subroutine to compute the overlap integral between s orbital and
c s orbital when the directional cosines of the two atoms is l, m, and n.
c
c In sss the first s means s orbital (angular quantum number is 0), the second
c s means s orbital (angular quantum number is 0), and the third s means sigma
c bonding (the magnetic quantum number for the two orbitals is 0, otherwise
c the sss is 0).
c

      subroutine ss(sss,dsss,hamil,dhamil)

      implicit none
      double precision sss,dsss,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some 
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hamil=sss

      dhamil(1)=dsss*l
      dhamil(2)=dsss*m
      dhamil(3)=dsss*n

      return
      end
c AJ 12/04
c removed MAXDIM.H header that set MAXATM, MAX3, and MAX6, replaced with NATOM passed
c changed call to SURF
c      SUBROUTINE SURF(V,X,DX,N3TM)
      subroutine surf(symb,v,x,dx,natom)
c
c The base orbitals will be ns + np for principal elements, (n-1)d + ns for
c transition elements.
c
      IMPLICIT NONE

      double precision autoang,autoev
      parameter (autoang=0.529177249d0)
      parameter (autoev = 27.2113961d0)
     
      character*75 space, hyphen
      parameter (space= "                                  ")
      parameter (hyphen="-----------------------------------------------
     *---------------------------")
c
c
c i,j,ip,jp,ii,jj,k --- temporary integer variables.
c
      integer i,j,ip,jp,ii,jj,k
c
c ntype --- number of atom types
      integer ntype,natom
      parameter(ntype=9)
c
c symb --- an array of atom labels (H, C, Ni, etc...)
c
      character*2 symb(natom)
c
c idtype --- identification for each element in the system.
c 	H   Hydrogen  identification number: 1
c 	C   Hydrogen  identification number: 2
c 	Ni  Nickel    identification number: 3
c 	Cu  Copper    identification number: 4
c 	Pd  Palladium identification number: 5
c 	Ag  Silver    identification number: 6
c 	Pt  Platium   identification number: 7
c 	Au  Gold      identification number: 8
c 	Al  Aluminum  identification number: 9
c ncharge --- charge of the system
c norbit --- holding the number of orbitals for each element.
c natom --- number of atoms in the system.
c nelem --- number of valence electrons for each atom.
c nElect --- number of valence electrons in the system.
c ndim --- dimension of the Hamiltonian matrix.
c
      integer ndim,nElect,ncharge
      integer idtype(natom),norbit(natom),nelem(natom)
c 
c nrdim,ncdim,temp --- temporary variables.
c
      integer nrdim,ncdim,temp
c
c X --- Cartesian coordinates, passed in bohr, converted later to A.
c vsip --- valence state ionization potential.
c
      double precision X(natom*3),vsip(ntype,2)
c
c hamil  --- interaction submatrix for two atoms.
c hamilt --- final Hamiltonian matrix for the system.
c
      double precision hamil(6,6),hamilt(natom*6,natom*6)
c
c ap  --- array passing Hamiltonian matrix to blas.
c work --- work space for blas
c evec --- molecular orbital coefficients
c eva  --- eigenvalues.
c occ  --- occupation number.
c

      double precision ap(natom*6*(natom*6+1)/2),work(3*natom*6)

      double precision evec(natom*6,natom*6),eva(natom*6)
      integer info,occ(natom*6)
c  
c dx --- first derivative of the system.
c dhamil --- First derivative of hamil.
c dhamlt --- First derivatives of dhamil.
c forcex,forcey,forcez,tmp --- temporary variables.
c
      double precision dx(natom*3),dhamlt(natom*6,natom*6,3),
     & dhamil(6,6,3)
      double precision forcex,forcey,forcez,tmp
c
c repul  --- pairwise repulsive energy.
c drepul --- first derivative of repul.
c eerep  --- total repulsive energy.
c
      double precision repul,drepul(3),eerep
c
c etb    --- valence band energy.
c V      --- potential energy.
c
      double precision etb,V
c
c q --- Mulliken net AO population
c qcharg --- gross AO populations for one atom.
c
      double precision q(natom*6),qcharg(natom)
      double precision upenalty
c
c
c Tzero --- Total zero energy.
c
      double precision Tzero
c
      common /coulom/ vsip
      common /penalty/ upenalty

c
c Main procedure begins here.
c

      ndim=0
      Tzero=0.0
      ncharge = 0 ! assume zero charge
      nElect = -ncharge ! number of "extra" electrons = -total charge
      upenalty = 3.d0 ! different atom combinations have different 
c                         default penalty energies
      do 20 i=1, natom
        if (symb(i).eq.'H '.or.symb(i).eq.' H') then
          nElect=nElect+1
          nelem(i)=1
          norbit(i)=1
          idtype(i)=1
          ndim=ndim+1
c          Tzero=Tzero-8.83
          Tzero=Tzero+vsip(1,1)
        else if (symb(i).eq.'C '.or.symb(i).eq.' C') then
          nElect=nElect+4
          nelem(i)=4
          norbit(i)=4
          idtype(i)=2
          ndim=ndim+4
c          Tzero=Tzero-53.26
          Tzero=Tzero+vsip(2,1)*2.d0+vsip(2,2)*2.d0+upenalty
        else if ( symb(i) .eq. 'Ni' .or. symb(i) .eq. 'NI') then
          nElect=nElect+10
          nelem(i)=10
          norbit(i)=6
          idtype(i)=3
          ndim=ndim+6
          Tzero=Tzero-87.80
        else if ( symb(i) .eq. 'Cu' .or. symb(i) .eq. 'CU') then
          nElect=nElect+11
          nelem(i)=11
          norbit(i)=6
          idtype(i)=4
          ndim=ndim+6
          Tzero=Tzero-93.50-8.65+3.20
        else if ( symb(i) .eq. 'Pd' .or. symb(i) .eq. 'PD') then
          nElect=nElect+10
          nelem(i)=10
          norbit(i)=6
          idtype(i)=5
          ndim=ndim+6
          Tzero=Tzero-91.50+3.20
        else if ( symb(i) .eq. 'Ag' .or. symb(i) .eq. 'AG') then
          nElect=nElect+11
          nelem(i)=11
          norbit(i)=6
          idtype(i)=6
          ndim=ndim+6
          Tzero=Tzero-94.50-8.65+3.20
        else if ( symb(i) .eq. 'Pt' .or. symb(i) .eq. 'PT') then
          nElect=nElect+10
          nelem(i)=10
          norbit(i)=6
          idtype(i)=7
          ndim=ndim+6
          Tzero=Tzero-92.50+3.20
        else if ( symb(i) .eq. 'Au' .or. symb(i) .eq. 'AU') then
          nElect=nElect+11
          nelem(i)=11
          norbit(i)=6
          idtype(i)=8
          ndim=ndim+6
          Tzero=Tzero-96.50-8.65+3.20
        else if ( symb(i) .eq. 'Al' .or. symb(i) .eq. 'AL') then
          nElect=nElect+3
          nelem(i)=3
          norbit(i)=4
          idtype(i)=9
          ndim=ndim+4
          Tzero=Tzero+(2.d0*vsip(9,1)+vsip(9,2))+upenalty
c
c USER SECTION BEGINS
c   Place the information of a new added element here.
c
c       else if ( symb(i) .eq. 'Ae' .or. symb(i) .eq. 'AE') then
c         nElect=nElect+x
c         nelem(i)=y
c         norbit(i)=z
c         idtype(i)=9
c         ndim=ndim+z
c
c USER SECTION ENDS
c
        else
          stop "I don't have the parameters to handle the element"
        endif
 20   continue
c
c Print the header information
c
c      write (30,*) "   Tight Binding Theory   "
c      write (30,*) space
      do i = 1, natom
c         write (30,95) symb(i),(X(3*(i-1)+j),j=1,3)
      enddo
c      write (30,*) hyphen
c      write (30,*) space
c
c Convert the unit of the coordinate from bohr to angstrom.
c
      do 11 i = 1, 3*natom
         X(i)=X(i)*autoang
   11 continue
c
c      write (30,*) "Following is the geometry file that was input in ang
c     *stroms"
c      write (30,*) hyphen
      do i = 1, natom
c         write (30,95) symb(i),(X(3*(i-1)+j),j=1,3)
      enddo
c      write (30,*) hyphen
c      write (30,*) space
c
      nrdim=0

c
c form matrix
      do 30 i=1,natom
        ip=idtype(i)
        ncdim=0
        do 31 j=1,natom
          jp=idtype(j)
          if (j.eq.i) then
            do 32 ii=1,norbit(i)
              do 33 jj=1,norbit(j)
                dhamlt(nrdim+ii,ncdim+jj,1)=0.0
                dhamlt(nrdim+ii,ncdim+jj,2)=0.0
                dhamlt(nrdim+ii,ncdim+jj,3)=0.0
                if (ii.eq.jj) then
                  if (ii .eq. 1) then
                    hamilt(nrdim+ii,ncdim+jj)=vsip(ip,1)
                  else
                    hamilt(nrdim+ii,ncdim+jj)=vsip(ip,2)
                  endif
                else 
                  hamilt(nrdim+ii,ncdim+jj)=0.0
                endif
 33           continue
 32         continue
            goto 39
          else
            call lmndlm(i,j,X,natom)
            call genrhd(ip,jp,hamil,dhamil)
            do 34 ii=1,norbit(i)
              do 35 jj=1,norbit(j)
                hamilt(nrdim+ii,ncdim+jj)=hamil(ii,jj)
                dhamlt(nrdim+ii,ncdim+jj,1)=dhamil(ii,jj,1)
                dhamlt(nrdim+ii,ncdim+jj,2)=dhamil(ii,jj,2)
                dhamlt(nrdim+ii,ncdim+jj,3)=dhamil(ii,jj,3)
 35           continue
 34         continue
          endif
 39       ncdim=ncdim+norbit(j)
 31     continue
        nrdim=nrdim+norbit(i)
 30   continue
c 
c Set up the uppermatrix to do eigenvalue calculations.
c
      do 40 i =1,ndim 
        do 41 j=i,ndim
          ap(i+(j-1)*j/2)=hamilt(i,j)
 41     continue 
 40   continue
c
      call dspev( 'v','u',ndim,ap,eva,evec,6*natom,work,info ) 
c      call dspev(21,ap,eva,evec,ndim,ndim,work,12*natom)
c
c Call a subroutine to arrange all the valence electrons to the necessary
c molecular orbitals so that the potential energy of the system is minimized.
c
      call filloc(ndim,nElect,eva,etb,occ,natom)
c
c Calculate Mulliken net AO Population for each atomic orbital
c Please refer to J. P. Lowe, Quantum Chemistry, p336 (Second edition,
c Acedemic press, 1993).  Note: Mulliken overlap population is zero since
c the overlap matrix here is assumed as unit matrix.
c
      do j = 1, ndim
        q(j)=0.0
        do i = 1, ndim
          q(j) = q(j) + occ(i)*evec(j,i)*evec(j,i)
        enddo
      enddo
c
      nrdim=0
      do i=1,natom
        qcharg(i)=0.0
        do j=1,norbit(i)
          qcharg(i)=qcharg(i)+q(nrdim+j)
c         do k=1,ndim
c           if (k.ne.(nrdim+j)) then
c             qcharg(i)=qcharg(i)+p(nrdim+j,k)/2.0
c           endif
c         enddo
        enddo
        nrdim=nrdim+norbit(i)
      enddo
c
c Calculate HOMO orbitals
c
c      write (30,*) "Following is information about HOMO and LUMO"
c      write (30,*) hyphen
c      write (30,97) " MO No.", "Eigenvalue(eV)","Description"
c      write (30,*) hyphen
      do i = 1, ndim
        if (occ(i).eq.2) then
          if (occ(i+1).eq.1) then
c            write(30,98) i,eva(i),' This is highest doubly occupied orbi
c     &tal.'
c            write(30,98) i+1,eva(i+1),' This is lowest singly occupied.'
          else if (occ(i+1).eq.0) then
c            write(30,98) i,eva(i),' This is highest doubly occupied orbi
c     &tal.'
c            write(30,*) "There is no singly-occupied orbital."
c            write(30,98) i,eva(i),' This is HOMO.'
            if ((i+1) .gt. ndim) then
c              write (30,*) 'There is no LUMO'
            else
c              write(30,98) i+1,eva(i+1),' This is LUMO. '
            endif
            goto 44
          endif
        else if (occ(i).eq.1) then
          if (occ(i+1).eq.0) then
c            write(30,98) i,eva(i),' This is HOMO.'
            if ((i+1) .gt. ndim) then
c              write (30,*) 'There is no LUMO'
            else
c              write(30,98) i+1,eva(i+1),' This is LUMO. '
            endif
            goto 44
          endif
        endif
      enddo
c
c Print Out eva, occ, q
c
 44   continue
c      write (30,*) hyphen
c      write (30,*) space
c
      do i=0,ndim/8
        if (i.lt.ndim/8) then
c          write (30,*) "Eigenvalues (E in eV), occupation number(n), and
c     * partial charge (q)"
c          write (30,*) " followed by Hamiltonian matrix H (in eV)"
c          write (30,*) hyphen
c          write (30,99) "  E   ",(eva(i*8+j),j=1,8)
c          write (30,103) "  n   ",(occ(i*8+j),j=1,8)
c          write (30,99) "  q   ",(q(i*8+j),j=1,8)
c          write (30,*) hyphen
          do k=1,ndim
c            write (30,100) "  H",k,(hamilt(i*8+j,k),j=1,8)
          enddo
c          write (30,*) hyphen
        else if (mod(ndim,8).ne.0) then
c          write (30,*) "Eigenvalues (E in eV), occupation number(n), and
c     * partial charge (q)"
c          write (30,*) "followed by overlap matrix S and Hamiltonian mat
c     *rix H(in eV)"
c          write (30,*) hyphen
c          write (30,99) "  E   ",(eva(j),j=8*i+1,ndim)
c          write (30,103) "  n   ",(occ(j),j=8*i+1,ndim)
c          write (30,99) "  q   ",(q(j),j=8*i+1,ndim)
c          write (30,*) hyphen
          do k=1,ndim
c            write (30,100) "  H",k,(hamilt(j,k),j=8*i+1,ndim)
          enddo
c          write (30,*) hyphen
        endif
      enddo
c      write (30,*) space
c
c Calculate the first derivative analytically.
c
        eerep=0.0
        nrdim=0
c
        do 51 i=1,natom
          ip=idtype(i)
          ncdim=0

          dx(3*(i-1)+1)=0.0
          dx(3*(i-1)+2)=0.0
          dx(3*(i-1)+3)=0.0

          do 52 j=1,natom
            jp=idtype(j)

            if (j.eq.i) goto 59
              
            call lmndlm(i,j,X,natom)
            call genrrp(ip,jp,repul,drepul)
            
            eerep=eerep+repul/2.0

            dx(3*(i-1)+1)=dx(3*(i-1)+1)+drepul(1)
            dx(3*(i-1)+2)=dx(3*(i-1)+2)+drepul(2)
            dx(3*(i-1)+3)=dx(3*(i-1)+3)+drepul(3)

            forcex=0.d0
            forcey=0.d0
            forcez=0.d0

            do 53 k=1,ndim
              do 54 ii=1,norbit(i)
                do 55 jj=1,norbit(j)
                  tmp=evec(nrdim+ii,k)*evec(ncdim+jj,k)*occ(k)
                  forcex=forcex+tmp*dhamlt(nrdim+ii,ncdim+jj,1)
                  forcey=forcey+tmp*dhamlt(nrdim+ii,ncdim+jj,2)
                  forcez=forcez+tmp*dhamlt(nrdim+ii,ncdim+jj,3)
 55             continue
 54           continue
 53         continue
            dx(3*(i-1)+1)=dx(3*(i-1)+1)+2.d0*forcex
            dx(3*(i-1)+2)=dx(3*(i-1)+2)+2.d0*forcey
            dx(3*(i-1)+3)=dx(3*(i-1)+3)+2.d0*forcez
 59         ncdim=ncdim+norbit(j)
 52       continue
          nrdim=nrdim+norbit(i)
c
 51     continue
c
c Convert Cartesian coordinates from Angstrom to bohr, first derivative from
c eV/(Angstrom) to Hartree/bohr.
c
c The binding energy is given in unit of eV.
c
      V=etb+eerep-Tzero
      V=V/autoev
      do 60 i = 1, 3*natom
        X(i)=X(i)/autoang
        DX(i)=DX(i)*autoang/autoev
   60 continue
c
c print out some information about the configuration.
c
c      write (30,102) "The energy of the system is ", V, " hartrees = ",
c     &V*autoev,"eV = ", V*627.50, " kcal/mol."
c      write (30,*) hyphen
c      write (30,*) space
c      write (30,*) "Element    Electrons  Charges     dx          dy
c     *   dz(a.u./bohr)"
c      write (30,*) hyphen
      do i=1,natom
c        write (30,101)symb(i),qcharg(i),nelem(i)-qcharg(i),
c     *         (dx(3*(i-1)+j),j=1,3)
      enddo
c      write (30,*) hyphen
c      write (30,*) "                                   THE END"
 95   format (T5,A,T10,3f10.5)
 97   format (T1,A,T13,A,T39,A)
 98   format (T1,I3,T9,f12.4,T30,A)
 99   format(T1,A8,8(f9.4))
 100  format(T1,A3,I3,8(f9.4))
 101  format(T3,A,T10,5(f10.5))
 102  format(T1,A,T26,f8.5,T34,A,T46,f10.3,T57,A,T61,f10.3,T71,A)
 103  format(T1,A8,8(I9))
c
      return
      END 
c
c This is a subroutine to compute the overlap integral between s
c orbital and dxxyy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In sds the first s means s orbital (angular quantum number is 0),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c

      subroutine sxxyy(sds,dsds,hamil,dhamil)

      implicit none
      double precision sds,dsds,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you want to compute the integral between dxxyy and s, just replace 
c sds and dsds with dss and dsss.
c

      hsqrt3=sqrt(3.0)/2.0

      hamil=hsqrt3*llmimm*sds

      dhamil(1)=hsqrt3*(llmimm*dsds*l+2.0*(l*dlx-m*dmx)*sds)
      dhamil(2)=hsqrt3*(llmimm*dsds*m+2.0*(l*dly-m*dmy)*sds)
      dhamil(3)=hsqrt3*(llmimm*dsds*n+2.0*(l*dlz-m*dmz)*sds)

      return
      end

c
c This is a subroutine to compute the overlap integral between s
c orbital and dxy (dyz or dzx) orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c sds is the two-center integral between the s orbital and d orbital when the
c line connecting the two atoms is the z axis, and dpps is the first derivative
c of sds with respect to the distance between the two atoms.
c
c In sds the first s means s orbital (angular quantum number is 0),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c For the input the character c is xy, yz, and zx reprenting the overlap
c integral between s orbital and dxy, or dyz, or dzx respectively.
c

      subroutine sxy(c,sds,dsds,hamil,dhamil)

      implicit none
      character*2 c
      double precision sds,dsds,el,em,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,sqrt3
      double precision delx,dely,delz,demx,demy,demz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      sqrt3=sqrt(3.0)

      if(c(1:1).eq.'x') then
        if(c(2:2).eq.'y') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
        endif
      else if(c(1:1).eq.'y')then
        if(c(2:2).eq.'z') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
      else if(c(1:1).eq.'z') then
        if(c(2:2).eq.'x') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hamil=sqrt3*el*em*sds

      dhamil(1)=sqrt3*((delx*em+el*demx)*sds+el*em*dsds*l)
      dhamil(2)=sqrt3*((dely*em+el*demy)*sds+el*em*dsds*m)
      dhamil(3)=sqrt3*((delz*em+el*demz)*sds+el*em*dsds*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between s orbital and
c px (py or pz) orbital when the directional cosines of the two atoms is l, m,
c and n.
c
c In sps the first s means s orbital (angular quantum number is 0), the 
c second p also means p orbital (angular quantum number is 1), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c For the input the character c is x, y, and z reprenting the overlap integral
c between s and px, py, and pz respectively.
c

      subroutine sxyz(c,sps,dsps,hamil,dhamil)

      implicit none
      character*1 c
      double precision sps,dsps,el,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,dx,dy,dz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if(c.eq.'x') then
        el=l
        dx=dlx
        dy=dly
        dz=dlz
      else if(c.eq.'y')then
        el=m
        dx=dmx
        dy=dmy
        dz=dmz
      else if(c.eq.'z')then
        el=n
        dx=dnx
        dy=dny
        dz=dnz
      else
        stop "Wrong in subroutine sxyz."
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c 
c If you want to get the overlap integral between px (py or pz) and s
c chang sps and dsps into pss and dpss, and reverse the sign of  the 
c following terms.
c

      hamil=el*sps

      dhamil(1)=el*dsps*l+dx*sps
      dhamil(2)=el*dsps*m+dy*sps
      dhamil(3)=el*dsps*n+dz*sps

      return
      end
c
c This is a subroutine to compute the overlap integral between s
c orbital and dxxyy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In sds the first s means s orbital (angular quantum number is 0),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c

      subroutine szz (sds,dsds,hamil,dhamil)

      implicit none
      double precision sds,dsds,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you want to compute the integral between dxxyy and s, just replace
c sds and dsds with dss and dsss.
c

      hamil=(1.5*nn-0.5)*sds

      dhamil(1)=(1.5*nn-0.5)*dsds*l+3.0*n*dnx*sds
      dhamil(2)=(1.5*nn-0.5)*dsds*m+3.0*n*dny*sds
      dhamil(3)=(1.5*nn-0.5)*dsds*n+3.0*n*dnz*sds

      return
      end

c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between Ni and C.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trc_hd (ip,hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,jp,i
      integer sNp(6),dNp(6)
      double precision szetapm(6),szeta(6),pzeta(6)
      double precision sss,dsss,dss,ddss,sps,dsps,dps,ddps,dpp,ddpp
      double precision coef1(6),coef2(6),dzeta1(6),dzeta2(6)
      double precision r,k(6),r0(6),mu0(6)
      double precision dkllp,kllp,x
      double precision tmp1,tmp2,tmp3,tmp4,temp,intgr,dintgr
c
      common /coulom/ vsip
      common /distan/ r
c
c STO parameters for transition element
c
      DATA (sNp(i),i=1,6) /4,4,5,5,6,6/
      DATA (dNp(i),i=1,6) /3,3,4,4,5,5/
      DATA (szetapm(i),i=1,6) /1.46952,2.14462,1.75195,1.81756,
     &                         2.04336,2.53186/
      DATA (coef1(i),i=1,6) /0.58579,0.59552,0.54226,0.56019,
     &                       0.65678,0.64695/
      DATA (coef2(i),i=1,6) /0.64856,0.57655,0.65649,0.55617,
     &                       0.57165,0.53788/
      DATA (dzeta1(i),i=1,6) /5.39452,5.89462,5.54495,4.90056,
     &                       5.50236,6.09286/
      DATA (dzeta2(i),i=1,6) /1.64452,2.24462,2.17495,2.49556,
     &                       2.18536,2.72386/
c
c STO parameters for carbon element
c
      DATA (szeta(i),i=1,6) /1.56068,1.13260,1.59663,1.52369,
     &                       1.57727,1.25328/
      DATA (pzeta(i),i=1,6) /1.56068,1.13260,1.59663,1.52369,
     &                       1.57727,1.25328/
c
c Wolfberg-Holmholz constants
c
      DATA (  k(i),i=1,6) /0.26441,0.18125,0.29889,0.20342,
     &                     0.34261,0.30277/
      DATA ( r0(i),i=1,6) /1.9928 ,1.2870,2.9265,1.5743,
     &                     2.0176 ,2.9275/
      DATA (mu0(i),i=1,6) /0.01756,1.27014,0.04703,1.37253,
     &                     0.01994,0.67298/
c
      jp=ip-2
c
      temp=k(jp)*(r/r0(jp))**mu0(jp)
      kllp=k(jp)+temp
      dkllp=mu0(jp)*temp/r
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c Actually, the final result is the product of overlap integral and constant
c K (distance dependent).
c
c s s sigma bonding overlap integral in z-axis
c
      call ovlcon(sNp(jp),0,szetapm(jp),2,0,0,szeta(jp),intgr,dintgr)
      dsss=dintgr*kllp+intgr*dkllp
       sss= intgr*kllp

c d s sigma bonding overlap integral in z-axis
c
      call ovlcon(dNp(jp),2,dzeta1(jp),2,0,0,szeta(jp),intgr,dintgr)
       dss= intgr*coef1(jp)
      ddss=dintgr*coef1(jp)

      call ovlcon(dNp(jp),2,dzeta2(jp),2,0,0,szeta(jp),intgr,dintgr)
       dss= intgr*coef2(jp)+ dss
      ddss=dintgr*coef2(jp)+ddss
c
      ddss=ddss*kllp+dss*dkllp
       dss= dss*kllp
c
c s p sigma bonding overlap integral in z-axis
c
      call ovlcon(sNp(jp),0,szetapm(jp),2,1,0,pzeta(jp),intgr,dintgr)
      dsps=dintgr*kllp+intgr*dkllp
       sps=intgr*kllp

c d p sigma bonding overlap integral in z-axis
c
      call ovlcon(dNp(jp),2,dzeta1(jp),2,1,0,pzeta(jp),intgr,dintgr)
       dps= intgr*coef1(jp)
      ddps=dintgr*coef1(jp)

      call ovlcon(dNp(jp),2,dzeta2(jp),2,1,0,pzeta(jp),intgr,dintgr)
       dps= intgr*coef2(jp)+ dps
      ddps=dintgr*coef2(jp)+ddps
c
      ddps=ddps*kllp+dps*dkllp
       dps= dps*kllp

c d p pi bonding overlap integral in z-axis
c
      call ovlcon(dNp(jp),2,dzeta1(jp),2,1,1,pzeta(jp),intgr,dintgr)
       dpp= intgr*coef1(jp)
      ddpp=dintgr*coef1(jp)

      call ovlcon(dNp(jp),2,dzeta2(jp),2,1,1,pzeta(jp),intgr,dintgr)
       dpp= intgr*coef2(jp)+ dpp
      ddpp=dintgr*coef2(jp)+ddpp
c
      ddpp=ddpp*kllp+dpp*dkllp
       dpp= dpp*kllp
c
c Multiply kmn with (Hm(Ni)+Hn(C))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(ip,1)+vsip(2,1))/2.0
      tmp2=(vsip(ip,1)+vsip(2,2))/2.0
      tmp3=(vsip(ip,2)+vsip(2,1))/2.0
      tmp4=(vsip(ip,2)+vsip(2,2))/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t*tmp1
      dhamil(1,1,1)=dt(1)*tmp1
      dhamil(1,1,2)=dt(2)*tmp1
      dhamil(1,1,3)=dt(3)*tmp1

      call sxyz('x',sps,dsps,t,dt)       !s-px 
      hamil(1,2)=t*tmp2
      dhamil(1,2,1)=dt(1)*tmp2
      dhamil(1,2,2)=dt(2)*tmp2
      dhamil(1,2,3)=dt(3)*tmp2

      call sxyz('y',sps,dsps,t,dt)       !s-py 
      hamil(1,3)=t*tmp2
      dhamil(1,3,1)=dt(1)*tmp2
      dhamil(1,3,2)=dt(2)*tmp2
      dhamil(1,3,3)=dt(3)*tmp2

      call sxyz('z',sps,dsps,t,dt)       !s-pz 
      hamil(1,4)=t*tmp2
      dhamil(1,4,1)=dt(1)*tmp2
      dhamil(1,4,2)=dt(2)*tmp2
      dhamil(1,4,3)=dt(3)*tmp2

      call sxy('xy',dss,ddss,t,dt)       !dxy-s
      hamil(2,1)=t*tmp3
      dhamil(2,1,1)=dt(1)*tmp3
      dhamil(2,1,2)=dt(2)*tmp3
      dhamil(2,1,3)=dt(3)*tmp3

      call xxy('xy',dps,ddps,dpp,ddpp,t,dt)         !dxy-px
      hamil(2,2)=t*tmp4
      dhamil(2,2,1)=dt(1)*tmp4
      dhamil(2,2,2)=dt(2)*tmp4
      dhamil(2,2,3)=dt(3)*tmp4

      call xxy('yx',dps,ddps,dpp,ddpp,t,dt)         !dxy-py
      hamil(2,3)=t*tmp4
      dhamil(2,3,1)=dt(1)*tmp4
      dhamil(2,3,2)=dt(2)*tmp4
      dhamil(2,3,3)=dt(3)*tmp4

      call xyz(dps,ddps,dpp,ddpp,t,dt)               !dxy-pz
      hamil(2,4)=t*tmp4
      dhamil(2,4,1)=dt(1)*tmp4
      dhamil(2,4,2)=dt(2)*tmp4
      dhamil(2,4,3)=dt(3)*tmp4

      call sxy('yz',dss,ddss,t,dt)       !dyz-s
      hamil(3,1)=t*tmp3
      dhamil(3,1,1)=dt(1)*tmp3
      dhamil(3,1,2)=dt(2)*tmp3
      dhamil(3,1,3)=dt(3)*tmp3

      call xyz(dps,ddps,dpp,ddpp,t,dt)              !dyz-px
      hamil(3,2)=t*tmp4
      dhamil(3,2,1)=dt(1)*tmp4
      dhamil(3,2,2)=dt(2)*tmp4
      dhamil(3,2,3)=dt(3)*tmp4

      call xxy('yz',dps,ddps,dpp,ddpp,t,dt)         !dyz-py
      hamil(3,3)=t*tmp4
      dhamil(3,3,1)=dt(1)*tmp4
      dhamil(3,3,2)=dt(2)*tmp4
      dhamil(3,3,3)=dt(3)*tmp4

      call xxy('zy',dps,ddps,dpp,ddpp,t,dt)          !dyz-pz
      hamil(3,4)=t*tmp4
      dhamil(3,4,1)=dt(1)*tmp4
      dhamil(3,4,2)=dt(2)*tmp4
      dhamil(3,4,3)=dt(3)*tmp4

      call sxy('zx',dss,ddss,t,dt)       !dzx-s
      hamil(4,1)=t*tmp3
      dhamil(4,1,1)=dt(1)*tmp3
      dhamil(4,1,2)=dt(2)*tmp3
      dhamil(4,1,3)=dt(3)*tmp3

      call xxy('xz',dps,ddps,dpp,ddpp,t,dt)         !dzx-px
      hamil(4,2)=t*tmp4
      dhamil(4,2,1)=dt(1)*tmp4
      dhamil(4,2,2)=dt(2)*tmp4
      dhamil(4,2,3)=dt(3)*tmp4

      call xyz(dps,ddps,dpp,ddpp,t,dt)              !dzx-py
      hamil(4,3)=t*tmp4
      dhamil(4,3,1)=dt(1)*tmp4
      dhamil(4,3,2)=dt(2)*tmp4
      dhamil(4,3,3)=dt(3)*tmp4

      call xxy('zx',dps,ddps,dpp,ddpp,t,dt)          !dzx-pz
      hamil(4,4)=t*tmp4
      dhamil(4,4,1)=dt(1)*tmp4
      dhamil(4,4,2)=dt(2)*tmp4
      dhamil(4,4,3)=dt(3)*tmp4

      call sxxyy(dss,ddss,t,dt)          !dx2_y2-s
      hamil(5,1)=t*tmp3
      dhamil(5,1,1)=dt(1)*tmp3
      dhamil(5,1,2)=dt(2)*tmp3
      dhamil(5,1,3)=dt(3)*tmp3

      call xxxyy(dps,ddps,dpp,ddpp,t,dt)            !dx2_y2-px
      hamil(5,2)=t*tmp4
      dhamil(5,2,1)=dt(1)*tmp4
      dhamil(5,2,2)=dt(2)*tmp4
      dhamil(5,2,3)=dt(3)*tmp4

      call yxxyy(dps,ddps,dpp,ddpp,t,dt)            !dx2_y2-py
      hamil(5,3)=t*tmp4
      dhamil(5,3,1)=dt(1)*tmp4
      dhamil(5,3,2)=dt(2)*tmp4
      dhamil(5,3,3)=dt(3)*tmp4

      call zxxyy(dps,ddps,dpp,ddpp,t,dt)             !dx2_y2-pz
      hamil(5,4)=t*tmp4
      dhamil(5,4,1)=dt(1)*tmp4
      dhamil(5,4,2)=dt(2)*tmp4
      dhamil(5,4,3)=dt(3)*tmp4

      call szz(dss,ddss,t,dt)            !dz2-s
      hamil(6,1)=t*tmp3
      dhamil(6,1,1)=dt(1)*tmp3
      dhamil(6,1,2)=dt(2)*tmp3
      dhamil(6,1,3)=dt(3)*tmp3

      call xzz('x',dps,ddps,dpp,ddpp,t,dt)              !dz2-px
      hamil(6,2)=t*tmp4
      dhamil(6,2,1)=dt(1)*tmp4
      dhamil(6,2,2)=dt(2)*tmp4
      dhamil(6,2,3)=dt(3)*tmp4

      call xzz('y',dps,ddps,dpp,ddpp,t,dt)              !dz2-py
      hamil(6,3)=t*tmp4
      dhamil(6,3,1)=dt(1)*tmp4
      dhamil(6,3,2)=dt(2)*tmp4
      dhamil(6,3,3)=dt(3)*tmp4

      call zzz(dps,ddps,dpp,ddpp,t,dt)               !dz2-pz
      hamil(6,4)=t*tmp4
      dhamil(6,4,1)=dt(1)*tmp4
      dhamil(6,4,2)=dt(2)*tmp4
      dhamil(6,4,3)=dt(3)*tmp4

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trh_hd (ip,hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3),vsip(9,2)
      integer ip,jp,i
      integer sNp(6),dNp(6)
      double precision szetapm(6),szeta(6)
      double precision sss,dsss,dss,ddss
      double precision kss(6),kds(6)
      double precision tmp1,tmp2,intgr,dintgr
      double precision coef1(6),coef2(6),dzeta1(6),dzeta2(6)

      common /coulom/ vsip
c
c STO parameters for transition element
c
      DATA (sNp(i),i=1,6) /4,4,5,5,6,6/
      DATA (dNp(i),i=1,6) /3,3,4,4,5,5/
      DATA (szetapm(i),i=1,6) /1.76329,1.81529,2.03423,1.74711,
     &                         2.36819,2.12506/
      DATA (coef1(i),i=1,6) /0.57107,0.61029,0.53160,0.56239,
     &                       0.64119,0.66519/
      DATA (coef2(i),i=1,6) /0.63227,0.59085,0.64359,0.55836,
     &                       0.55808,0.55304/
      DATA (dzeta1(i),i=1,6) /5.68829,5.56529,5.82723,4.83011,
     &                       5.82719,5.68606/
      DATA (dzeta2(i),i=1,6) /1.93829,1.91529,2.45723,2.42511,
     &                       2.51019,2.31706/
c
c STO parameters for hydrogen element
c
      DATA (szeta(i),i=1,6) /1.22576,1.26397,1.21011,1.23966,
     &                       1.25792,1.12836/
c
c Wolfberg-Hol relation k constant.
c
      DATA (kss(i),i=1,6) /0.75730 ,0.80843,0.67048,0.57090,
     &                     0.91655,0.61968/
      DATA (kds(i),i=1,6) /0.18030 ,0.02514,0.03081,0.05065,
     &                     0.31296,0.003524/
c
      jp=ip-2
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
c
      call ovlcon(sNp(jp),0,szetapm(jp),1,0,0,szeta(jp),intgr,dintgr)
       sss=intgr
      dsss=dintgr

c d s sigma bonding overlap integral in z-axis
c
      call ovlcon(dNp(jp),2,dzeta1(jp),1,0,0,szeta(jp),intgr,dintgr)
       dss= intgr*coef1(jp)
      ddss=dintgr*coef1(jp)

      call ovlcon(dNp(jp),2,dzeta2(jp),1,0,0,szeta(jp),intgr,dintgr)
       dss= intgr*coef2(jp)+ dss
      ddss=dintgr*coef2(jp)+ddss
c
c Multiply kmn with (Hm(Ni)+Hn(H))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(ip,1)+vsip(1,1))*kss(jp)/2.0
      tmp2=(vsip(ip,2)+vsip(1,1))*kds(jp)/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)                          !s-s
      hamil(1,1)=t*tmp1
      dhamil(1,1,1)=dt(1)*tmp1
      dhamil(1,1,2)=dt(2)*tmp1
      dhamil(1,1,3)=dt(3)*tmp1

      call sxy('xy',dss,ddss,t,dt)                    !dxy-s
      hamil(2,1)=t*tmp2
      dhamil(2,1,1)=dt(1)*tmp2
      dhamil(2,1,2)=dt(2)*tmp2
      dhamil(2,1,3)=dt(3)*tmp2

      call sxy('yz',dss,ddss,t,dt)                    !dyz-s
      hamil(3,1)=t*tmp2
      dhamil(3,1,1)=dt(1)*tmp2
      dhamil(3,1,2)=dt(2)*tmp2
      dhamil(3,1,3)=dt(3)*tmp2

      call sxy('zx',dss,ddss,t,dt)                    !dzx-s
      hamil(4,1)=t*tmp2
      dhamil(4,1,1)=dt(1)*tmp2
      dhamil(4,1,2)=dt(2)*tmp2
      dhamil(4,1,3)=dt(3)*tmp2

      call sxxyy(dss,ddss,t,dt)                       !dx2_y2-s
      hamil(5,1)=t*tmp2
      dhamil(5,1,1)=dt(1)*tmp2
      dhamil(5,1,2)=dt(2)*tmp2
      dhamil(5,1,3)=dt(3)*tmp2

      call szz(dss,ddss,t,dt)                         !dz2-s
      hamil(6,1)=t*tmp2
      dhamil(6,1,1)=dt(1)*tmp2
      dhamil(6,1,2)=dt(2)*tmp2
      dhamil(6,1,3)=dt(3)*tmp2

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between two same transition
c metal atoms.
c If you want to some calculation between two different metal atoms, find
c the parameters first.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trtrhd (ip,hamil,dhamil)
c
      implicit none
      integer i,ip,jp
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision ssseta0(6),rd(6),d(6),ssstau,sssalpha(6),sss,dsss
      double precision sdseta0,             sdstau,sdsalpha(6),sds,dsds
      double precision dsseta0,             dsstau,dssalpha(6),dss,ddss
      double precision ddseta0,             ddstau,ddsalpha(6),dds,ddds
      double precision ddpeta0,             ddptau,ddpalpha(6),ddp,dddp
      double precision dddeta0,             dddtau,dddalpha(6),ddd,dddd
c
c Hopping parameters for metal and metal
c
      DATA ( d(i),i=1,6) /2.4761,1.9606,2.8120,2.4307,2.5862,2.3083/
      DATA (rd(i),i=1,6) /0.71,0.67,0.94,0.89,1.04,1.01/
c s s sigma
      DATA (ssseta0(i),i=1,6) /-0.4530,-0.5447,-0.5652,-0.4340,
     &                         -0.6322,-0.3817/
      DATA ssstau /0.d0/
      DATA (sssalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c s d sigma
      DATA sdseta0 /-3.16d0/
      DATA sdstau /1.5d0/
      DATA (sdsalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c d s sigma
      DATA dsseta0 /-3.16d0/
      DATA dsstau /1.5d0/
      DATA (dssalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c d d sigma
      DATA ddseta0 /-16.2d0/
      DATA ddstau /3.d0/
      DATA (ddsalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c d d pi
      DATA ddpeta0 /8.75d0/
      DATA ddptau /3.d0/
      DATA (ddpalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c d d delta
      DATA dddeta0 /0.0d0/
      DATA dddtau /3.d0/
      DATA (dddalpha(i),i=1,6) /0.41132,0.60676,0.41411,0.50593,
     &                          0.51571,0.58339/
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      jp=ip-2
c
      call lathio(ssseta0(jp),rd(jp),d(jp),ssstau,sssalpha(jp),sss,dsss)
      call lathio(sdseta0,    rd(jp),d(jp),sdstau,sdsalpha(jp),sds,dsds)
      call lathio(dsseta0,    rd(jp),d(jp),dsstau,dssalpha(jp),dss,ddss)
      call lathio(ddseta0,    rd(jp),d(jp),ddstau,ddsalpha(jp),dds,ddds)
      call lathio(ddpeta0,    rd(jp),d(jp),ddptau,ddpalpha(jp),ddp,dddp)
      call lathio(dddeta0,    rd(jp),d(jp),dddtau,dddalpha(jp),ddd,dddd)
c 
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxy('xy',sds,dsds,t,dt)       !s-dxy
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)

      call sxy('yz',sds,dsds,t,dt)       !s-dyz
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)

      call sxy('zx',sds,dsds,t,dt)       !s-dzx
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)

      call sxxyy(sds,dsds,t,dt)          !s-dx2_y2
      hamil(1,5)=t
      dhamil(1,5,1)=dt(1)
      dhamil(1,5,2)=dt(2)
      dhamil(1,5,3)=dt(3)

      call szz(sds,dsds,t,dt)            !s-dz2
      hamil(1,6)=t
      dhamil(1,6,1)=dt(1)
      dhamil(1,6,2)=dt(2)
      dhamil(1,6,3)=dt(3)

      call sxy('xy',dss,ddss,t,dt)             !dxy-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)

      call xyxy('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dxy
      hamil(2,2)=t
      dhamil(2,2,1)=dt(1)
      dhamil(2,2,2)=dt(2)
      dhamil(2,2,3)=dt(3)

      call xyyz('xz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dyz
      hamil(2,3)=t
      dhamil(2,3,1)=dt(1)
      dhamil(2,3,2)=dt(2)
      dhamil(2,3,3)=dt(3)

      call xyyz('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dzx
      hamil(2,4)=t
      dhamil(2,4,1)=dt(1)
      dhamil(2,4,2)=dt(2)
      dhamil(2,4,3)=dt(3)

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dxy-dx2_y2
      hamil(2,5)=t
      dhamil(2,5,1)=dt(1)
      dhamil(2,5,2)=dt(2)
      dhamil(2,5,3)=dt(3)

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      hamil(2,6)=t
      dhamil(2,6,1)=dt(1)
      dhamil(2,6,2)=dt(2)
      dhamil(2,6,3)=dt(3)

      call sxy('yz',dss,ddss,t,dt)                     !dyz-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)

      call xyyz('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dxy
      hamil(3,2)=t
      dhamil(3,2,1)=dt(1)
      dhamil(3,2,2)=dt(2)
      dhamil(3,2,3)=dt(3)

      call xyxy('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dyz
      hamil(3,3)=t
      dhamil(3,3,1)=dt(1)
      dhamil(3,3,2)=dt(2)
      dhamil(3,3,3)=dt(3)

      call xyyz('yx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dzx
      hamil(3,4)=t
      dhamil(3,4,1)=dt(1)
      dhamil(3,4,2)=dt(2)
      dhamil(3,4,3)=dt(3)

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dxy-dx2_y2
      hamil(3,5)=t
      dhamil(3,5,1)=dt(1)
      dhamil(3,5,2)=dt(2)
      dhamil(3,5,3)=dt(3)

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      hamil(3,6)=t
      dhamil(3,6,1)=dt(1)
      dhamil(3,6,2)=dt(2)
      dhamil(3,6,3)=dt(3)

      call sxy('zx',dss,ddss,t,dt)                     !dzx-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)

      call xyyz('zy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dxy
      hamil(4,2)=t
      dhamil(4,2,1)=dt(1)
      dhamil(4,2,2)=dt(2)
      dhamil(4,2,3)=dt(3)

      call xyyz('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dyz
      hamil(4,3)=t
      dhamil(4,3,1)=dt(1)
      dhamil(4,3,2)=dt(2)
      dhamil(4,3,3)=dt(3)

      call xyxy('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dzx
      hamil(4,4)=t
      dhamil(4,4,1)=dt(1)
      dhamil(4,4,2)=dt(2)
      dhamil(4,4,3)=dt(3)

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      hamil(4,5)=t
      dhamil(4,5,1)=dt(1)
      dhamil(4,5,2)=dt(2)
      dhamil(4,5,3)=dt(3)

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dzx-dz2
      hamil(4,6)=t
      dhamil(4,6,1)=dt(1)
      dhamil(4,6,2)=dt(2)
      dhamil(4,6,3)=dt(3)

      call sxxyy(sds,dsds,t,dt)                       !dx2_y2-s
      hamil(5,1)=t
      dhamil(5,1,1)=dt(1)
      dhamil(5,1,2)=dt(2)
      dhamil(5,1,3)=dt(3)

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      hamil(5,2)=t
      dhamil(5,2,1)=dt(1)
      dhamil(5,2,2)=dt(2)
      dhamil(5,2,3)=dt(3)

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      hamil(5,3)=t
      dhamil(5,3,1)=dt(1)
      dhamil(5,3,2)=dt(2)
      dhamil(5,3,3)=dt(3)

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dzx
      hamil(5,4)=t
      dhamil(5,4,1)=dt(1)
      dhamil(5,4,2)=dt(2)
      dhamil(5,4,3)=dt(3)

      call xxyy_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dx2_y2
      hamil(5,5)=t
      dhamil(5,5,1)=dt(1)
      dhamil(5,5,2)=dt(2)
      dhamil(5,5,3)=dt(3)

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      hamil(5,6)=t
      dhamil(5,6,1)=dt(1)
      dhamil(5,6,2)=dt(2)
      dhamil(5,6,3)=dt(3)

      call szz(sds,dsds,t,dt)                          !dz2-s
      hamil(6,1)=t
      dhamil(6,1,1)=dt(1)
      dhamil(6,1,2)=dt(2)
      dhamil(6,1,3)=dt(3)

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      hamil(6,2)=t
      dhamil(6,2,1)=dt(1)
      dhamil(6,2,2)=dt(2)
      dhamil(6,2,3)=dt(3)

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      hamil(6,3)=t
      dhamil(6,3,1)=dt(1)
      dhamil(6,3,2)=dt(2)
      dhamil(6,3,3)=dt(3)

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dzx
      hamil(6,4)=t
      dhamil(6,4,1)=dt(1)
      dhamil(6,4,2)=dt(2)
      dhamil(6,4,3)=dt(3)

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dz2-dx2_y2
      hamil(6,5)=t
      dhamil(6,5,1)=dt(1)
      dhamil(6,5,2)=dt(2)
      dhamil(6,5,3)=dt(3)

      call zz_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dz2
      hamil(6,6)=t
      dhamil(6,6,1)=dt(1)
      dhamil(6,6,2)=dt(2)
      dhamil(6,6,3)=dt(3)

      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Ni and Ni
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trtrrp (ip,repul,drepul)
c
      implicit none
      integer i,ip,jp
      double precision repul,dt,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision alpha(6),gamma(6),beta(6)
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for two same transition metal atom.
c
      DATA (alpha(i),i=1,6) /4593.4 ,3921.32,10488.1,11040.6,
     &                       9929.2 ,9594.20/
      DATA (gamma(i),i=1,6) /0.94507,1.13089,0.50079,0.41131,
     &                       0.98402,0.64168/
      DATA (beta(i), i=1,6) /4.3988 ,4.2130,3.9416 ,3.8614 ,
     &                       3.8869,3.6789/
c
      jp=ip-2
c
      call grepfm (alpha(jp),gamma(jp),beta(jp),repul,dt)
      
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
c 88888888888888888888888888888888888888888888888888888888888888888888888
c The following function is taken from
c Y. Wang and C.H. Mak, Chem. Phys. Lett. 235, 37-46 (1995).
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine wangma (t0,r0,na,nb,rt,nc,t,dt)
c
      implicit none
      double precision t0,r0,na,nb,rt,nc,t,dt,r,tmp

      common /distan/ r

      tmp=(r/rt)**nc

      t=t0*(r0/r)**na*(exp(-tmp+(r0/rt)**nc))**nb
      dt=(-na-nb*nc*tmp)*t/r

      return
      end
c
c This is a subroutine to compute the overlap integral between px (py, pz)
c orbital and px (py or pz) orbital when the directional cosines of the two 
c atoms is l, m, and n. 
c
c In pps the first p means p orbital (angular quantum number is 1), 
c the second p means p orbital (angular quantum number is 1), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ppp the first p means p orbital (angular quantum number is 1), 
c the second p means p orbital (angular quantum number is 1), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1 
c or -1, otherwise the ppp is 0).
c
c For the input the character c is x, y, and z reprenting the overlap integral
c between px and px, py and py, or pz and pz respectively.
c

      subroutine xx(c,pps,dpps,ppp,dppp,hamil,dhamil)

      implicit none
      character*1 c
      double precision pps,dpps,ppp,dppp,el,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,temp,dtemp
      double precision delx,dely,delz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if(c(1:1).eq.'x') then
        el=l
        delx=dlx
        dely=dly
        delz=dlz
      else if(c(1:1).eq.'y')then
        el=m
        delx=dmx
        dely=dmy
        delz=dmz
      else if(c(1:1).eq.'z') then
        el=n
        delx=dnx
        dely=dny
        delz=dnz
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      temp=pps-ppp
      dtemp=dpps-dppp

      hamil=el*el*temp+ppp

      dhamil(1)=2.0*el*delx*temp+(el*el*dtemp+dppp)*l
      dhamil(2)=2.0*el*dely*temp+(el*el*dtemp+dppp)*m
      dhamil(3)=2.0*el*delz*temp+(el*el*dtemp+dppp)*n

      return
      end
c
c This is a subroutine to compute the overlap integral between px
c orbital and dx2-y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c

      subroutine xxxyy(pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      double precision pds,dpds,pdp,dpdp,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3,tmp,dtmp
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      hsqrt3=sqrt(3.0)/2.0

      tmp=hsqrt3*pds-pdp
      dtmp=hsqrt3*dpds-dpdp

      hamil=l*(llmimm*tmp+pdp)

      dhamil(1)=dlx*(llmimm*tmp+pdp)+l*(2.0*(l*dlx-m*dmx)*tmp
     +    +(llmimm*dtmp+dpdp)*l)

      dhamil(2)=dly*(llmimm*tmp+pdp)+l*(2.0*(l*dly-m*dmy)*tmp
     +    +(llmimm*dtmp+dpdp)*m)

      dhamil(3)=dlz*(llmimm*tmp+pdp)+l*(2.0*(l*dlz-m*dmz)*tmp
     +    +(llmimm*dtmp+dpdp)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between px (py, pz)
c orbital and dxy(dyz or dzx) orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c For the input the character c is xy, xz, yz, yx, zy, and zx reprenting the 
c overlap integral between px and dxy, px and dzx, py and dxy, py and dyz, 
c pz and dyz, and pz and dzx respectively.
c

      subroutine xxy(c,pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      character*2 c
      double precision pds,dpds,pdp,dpdp,el,em,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,sqrt3,tmp,dtmp
      double precision delx,dely,delz,demx,demy,demz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      sqrt3=sqrt(3.0)

      if(c(1:1).eq.'x') then
        if(c(2:2).eq.'y') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
        else if(c(2:2).eq.'z') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
      else if(c(1:1).eq.'y')then
        if(c(2:2).eq.'z') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
        if(c(2:2).eq.'x') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=l
          demx=dlx
          demy=dly
          demz=dlz
        endif
      else if(c(1:1).eq.'z')then
        if(c(2:2).eq.'y') then
          el=n
          delx=dnx
          dely=dny
          delz=dnz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
        endif
        if(c(2:2).eq.'x') then
          el=n
          delx=dnx
          dely=dny
          delz=dnz
          em=l
          demx=dlx
          demy=dly
          demz=dlz
        endif
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dxy and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      tmp=sqrt3*pds-2.0*pdp
      dtmp=sqrt3*dpds-2.0*dpdp

      hamil=em*(el*el*tmp+pdp)

      dhamil(1)=demx*(el*el*tmp+pdp)+em*(2.0*el*delx*tmp
     +    +(el*el*dtmp+dpdp)*l)    

      dhamil(2)=demy*(el*el*tmp+pdp)+em*(2.0*el*dely*tmp
     +    +(el*el*dtmp+dpdp)*m)    

      dhamil(3)=demz*(el*el*tmp+pdp)+em*(2.0*el*delz*tmp
     +    +(el*el*dtmp+dpdp)*n)    

      return
      end
c
c This is a subroutine to compute the overlap integral between dx2-y2
c orbital and dx2-y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine xxyy_2 (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz


      tmpa=(3.0*dds-4.0*ddp+ddd)/4.0
      dtmpa=(3.0*ddds-4.0*dddp+dddd)/4.0
      tmpb=ddp-ddd
      dtmpb=dddp-dddd

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      hamil=llmimm*llmimm*tmpa+lladmm*tmpb+ddd

      dhamil(1)=4.0*llmimm*(l*dlx-m*dmx)*tmpa+llmimm*llmimm*dtmpa*l
     +         +2.0*(l*dlx+m*dmx)*tmpb+(lladmm*dtmpb+dddd)*l

      dhamil(2)=4.0*llmimm*(l*dly-m*dmy)*tmpa+llmimm*llmimm*dtmpa*m
     +         +2.0*(l*dly+m*dmy)*tmpb+(lladmm*dtmpb+dddd)*m

      dhamil(3)=4.0*llmimm*(l*dlz-m*dmz)*tmpa+llmimm*llmimm*dtmpa*n
     +         +2.0*(l*dlz+m*dmz)*tmpb+(lladmm*dtmpb+dddd)*n

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dxy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine xxyyzz (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,qsqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      qsqrt3=sqrt(3.0)/4.0

      tmpa=qsqrt3*(3.0*dds-4.0*ddp+ddd)
      dtmpa=qsqrt3*(3.0*ddds-4.0*dddp+dddd)
      tmpb=-qsqrt3*(dds-ddd)
      dtmpb=-qsqrt3*(ddds-dddd)

      hamil=llmimm*(nn*tmpa+tmpb)

      dhamil(1)=2.0*(l*dlx-m*dmx)*(nn*tmpa+tmpb)
     +    +llmimm*(2.0*n*dnx*tmpa+(nn*dtmpa+dtmpb)*l)

      dhamil(2)=2.0*(l*dly-m*dmy)*(nn*tmpa+tmpb)
     +    +llmimm*(2.0*n*dny*tmpa+(nn*dtmpa+dtmpb)*m)

      dhamil(3)=2.0*(l*dlz-m*dmz)*(nn*tmpa+tmpb)
     +    +llmimm*(2.0*n*dnz*tmpa+(nn*dtmpa+dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between px (py, pz)
c orbital and py (pz or px) orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pps the first p means p orbital (angular quantum number is 1),
c the second p means p orbital (angular quantum number is 1), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ppp the first p means p orbital (angular quantum number is 1),
c the second p means p orbital (angular quantum number is 1), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c For the input the character c is xy, xz, and yz reprenting the overlap 
c integral between px and py or py and px, px and pz or pz and px, and py and
c pz or pz and py respectively.
c

      subroutine xy(c,pps,dpps,ppp,dppp,hamil,dhamil)

      implicit none
      character*2 c
      double precision pps,dpps,ppp,dppp,el,em,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,temp,dtemp
      double precision delx,dely,delz,demx,demy,demz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if(c(1:1).eq.'x') then
        if(c(2:2).eq.'y') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
        else if(c(2:2).eq.'z') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
      else if(c(1:1).eq.'y')then
        if(c(2:2).eq.'z') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
        endif
      endif

      temp=pps-ppp
      dtemp=dpps-dppp

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hamil=el*em*temp

      dhamil(1)=(delx*em+el*demx)*temp+el*em*dtemp*l
      dhamil(2)=(dely*em+el*demy)*temp+el*em*dtemp*m
      dhamil(3)=(delz*em+el*demz)*temp+el*em*dtemp*n

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dx2_y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine xyxxyy (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      tmpa=1.5*dds-2.0*ddp+0.5*ddd
      dtmpa=1.5*ddds-2.0*dddp+0.5*dddd

      hamil=l*m*llmimm*tmpa

      dhamil(1)=(dlx*m*llmimm+l*dmx*llmimm+l*m*2.0*(l*dlx-m*dmx))*tmpa
     +    +l*m*llmimm*dtmpa*l

      dhamil(2)=(dly*m*llmimm+l*dmy*llmimm+l*m*2.0*(l*dly-m*dmy))*tmpa
     +    +l*m*llmimm*dtmpa*m

      dhamil(3)=(dlz*m*llmimm+l*dmz*llmimm+l*m*2.0*(l*dlz-m*dmz))*tmpa
     +    +l*m*llmimm*dtmpa*n

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dxy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c
c Character c xy representing dxy and dxy, yz dyz and dyz, zx dzx and dzx.
c

      subroutine xyxy(c,dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      character*2 c
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,el,em,en
      double precision delx,dely,delz,demx,demy,demz,denx,deny,denz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,tmpb,dtmpa,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if (c(1:2).eq.'xy') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
          en=n
          denx=dnx
          deny=dny
          denz=dnz
      else if (c(1:2).eq.'yz') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
          en=l
          denx=dlx
          deny=dly
          denz=dlz
      else if (c(1:2).eq.'zx') then
          el=n
          delx=dnx
          dely=dny
          delz=dnz
          em=l
          demx=dlx
          demy=dly
          demz=dlz
          en=m
          denx=dmx
          deny=dmy
          denz=dmz
      else
        stop "Wrong for input to subroutine xyxy."
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      tmpa=3.0*dds-4.0*ddp+ddd
      dtmpa=3.0*ddds-4.0*dddp+dddd
      tmpb=ddp-ddd
      dtmpb=dddp-dddd

      hamil=el*el*em*em*tmpa+(el*el+em*em)*tmpb+ddd

      dhamil(1)=2.0*(el*delx*em*em+el*el*em*demx)*tmpa
     +         +(el*el*em*em*dtmpa+dddd+(el*el+em*em)*dtmpb)*l
     +         +2.0*(el*delx+em*demx)*tmpb

      dhamil(2)=2.0*(el*dely*em*em+el*el*em*demy)*tmpa
     +         +(el*el*em*em*dtmpa+dddd+(el*el+em*em)*dtmpb)*m
     +         +2.0*(el*dely+em*demy)*tmpb

      dhamil(3)=2.0*(el*delz*em*em+el*el*em*demz)*tmpa
     +         +(el*el*em*em*dtmpa+dddd+(el*el+em*em)*dtmpb)*n
     +         +2.0*(el*delz+em*demz)*tmpb

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dyz orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c
c Character c xz representing dxy and dyz or dyz and dxy, yz dxy and dzx or 
c dzx and dxy, xy dyz and dzx or dzx and dyz.
c

      subroutine xyyz (c,dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      character*2 c
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,el,em,en
      double precision delx,dely,delz,demx,demy,demz,denx,deny,denz
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,tmpb,dtmpa,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if (c(1:2).eq.'xz'.or.c(1:2).eq.'zx') then
          el=l
          delx=dlx
          dely=dly
          delz=dlz
          em=m
          demx=dmx
          demy=dmy
          demz=dmz
          en=n
          denx=dnx
          deny=dny
          denz=dnz
      else if (c(1:2).eq.'zy'.or.c(1:2).eq.'yz') then
          el=n
          delx=dnx
          dely=dny
          delz=dnz
          em=l
          demx=dlx
          demy=dly
          demz=dlz
          en=m
          denx=dmx
          deny=dmy
          denz=dmz
      else if (c(1:2).eq.'yx'.or.c(1:2).eq.'xy') then
          el=m
          delx=dmx
          dely=dmy
          delz=dmz
          em=n
          demx=dnx
          demy=dny
          demz=dnz
          en=l
          denx=dlx
          deny=dly
          denz=dlz
      else
        stop "Wrong for input to subroutine xyyz."
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      tmpa=3.0*dds-4.0*ddp+ddd
      dtmpa=3.0*ddds-4.0*dddp+dddd
      tmpb=ddp-ddd
      dtmpb=dddp-dddd

      hamil=el*en*(em*em*tmpa+tmpb)

      dhamil(1)=(delx*en+el*denx)*(em*em*tmpa+tmpb)
     +    +el*en*(2.0*em*demx*tmpa+(em*em*dtmpa+dtmpb)*l)

      dhamil(2)=(dely*en+el*deny)*(em*em*tmpa+tmpb)
     +    +el*en*(2.0*em*demy*tmpa+(em*em*dtmpa+dtmpb)*m)

      dhamil(3)=(delz*en+el*denz)*(em*em*tmpa+tmpb)
     +    +el*en*(2.0*em*demz*tmpa+(em*em*dtmpa+dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between px (py, pz)
c orbital and px (py or pz) orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c This is for px and dyz, py and dzx, or pz and dxy.
c

      subroutine xyz(pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      double precision pds,dpds,pdp,dpdp,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,sqrt3,temp,dtemp
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dyz and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      sqrt3=sqrt(3.0)

      temp=sqrt3*pds-2.0*pdp
      dtemp=sqrt3*dpds-2.0*dpdp

      hamil=l*m*n*temp

      dhamil(1)=(dlx*m*n+l*dmx*n+l*m*dnx)*temp+l*m*n*dtemp*l
      dhamil(2)=(dly*m*n+l*dmy*n+l*m*dny)*temp+l*m*n*dtemp*m
      dhamil(3)=(dlz*m*n+l*dmz*n+l*m*dnz)*temp+l*m*n*dtemp*n

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dz2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine xyzz (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hsqrt3=sqrt(3.0)/2.0

      tmpa=hsqrt3*(3.0*dds-4.0*ddp+ddd)
      dtmpa=hsqrt3*(3.0*ddds-4.0*dddp+dddd)
      tmpb=hsqrt3*(dds-ddd)
      dtmpb=hsqrt3*(ddds-dddd)

      hamil=l*m*(nn*tmpa-tmpb)

      dhamil(1)=(dlx*m+l*dmx)*(nn*tmpa-tmpb)
     +  +l*m*(2.0*n*dnx*tmpa+(nn*dtmpa-dtmpb)*l)

      dhamil(2)=(dly*m+l*dmy)*(nn*tmpa-tmpb)
     +  +l*m*(2.0*n*dny*tmpa+(nn*dtmpa-dtmpb)*m)

      dhamil(3)=(dlz*m+l*dmz)*(nn*tmpa-tmpb)
     +  +l*m*(2.0*n*dnz*tmpa+(nn*dtmpa-dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between px (or py)
c orbital and dz2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c The input character c with x, y representing px, py with dz2.
c

      subroutine xzz(c,pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      character*1 c
      double precision pds,dpds,pdp,dpdp,el,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,delx,dely,delz,sqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz,tmp,dtmp

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      if (c.eq.'x') then
        el=l
        delx=dlx
        dely=dly
        delz=dlz
      else if (c.eq.'y') then
        el=m
        delx=dmx
        dely=dmy
        delz=dmz
      else
        stop "Wrong for input to subroutine xzz."
      endif

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dz2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      sqrt3=sqrt(3.0)

      tmp=1.5*pds-sqrt3*pdp
      dtmp=1.5*dpds-sqrt3*dpdp

      hamil=el*(nn*tmp-0.5*pds)

      dhamil(1)=delx*(nn*tmp-0.5*pds)+el*(2.0*n*dnx*tmp
     +    +(nn*dtmp-0.5*dpds)*l)

      dhamil(2)=dely*(nn*tmp-0.5*pds)+el*(2.0*n*dny*tmp
     +    +(nn*dtmp-0.5*dpds)*m)

      dhamil(3)=delz*(nn*tmp-0.5*pds)+el*(2.0*n*dnz*tmp
     +    +(nn*dtmp-0.5*dpds)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between py
c orbital and dx2-y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c

      subroutine yxxyy(pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      double precision pds,dpds,pdp,dpdp,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3,tmp,dtmp
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and py replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      hsqrt3=sqrt(3.0)/2.0

      tmp=hsqrt3*pds-pdp
      dtmp=hsqrt3*dpds-dpdp

      hamil=m*(llmimm*tmp-pdp)

      dhamil(1)=dmx*(llmimm*tmp-pdp)+m*(2.0*(l*dlx-m*dmx)*tmp
     +    +(llmimm*dtmp-dpdp)*l)

      dhamil(2)=dmy*(llmimm*tmp-pdp)+m*(2.0*(l*dly-m*dmy)*tmp
     +    +(llmimm*dtmp-dpdp)*m)

      dhamil(3)=dmz*(llmimm*tmp-pdp)+m*(2.0*(l*dlz-m*dmz)*tmp
     +    +(llmimm*dtmp-dpdp)*n)

      return
      end

c
c This is a subroutine to compute the overlap integral between dyz
c orbital and dx2_y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine yzxxyy (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      tmpa=1.5*dds-2.0*ddp+0.5*ddd
      dtmpa=1.5*ddds-2.0*dddp+0.5*dddd
      tmpb=ddp-ddd
      dtmpb=dddp-dddd

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hamil=m*n*(llmimm*tmpa-tmpb)

      dhamil(1)=(dmx*n+m*dnx)*(llmimm*tmpa-tmpb)
     +    +m*n*(2.0*(l*dlx-m*dmx)*tmpa+(llmimm*dtmpa-dtmpb)*l)

      dhamil(2)=(dmy*n+m*dny)*(llmimm*tmpa-tmpb)
     +    +m*n*(2.0*(l*dly-m*dmy)*tmpa+(llmimm*dtmpa-dtmpb)*m)

      dhamil(3)=(dmz*n+m*dnz)*(llmimm*tmpa-tmpb)
     +    +m*n*(2.0*(l*dlz-m*dmz)*tmpa+(llmimm*dtmpa-dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between dyz
c orbital and dz2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine yzzz (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hsqrt3=sqrt(3.0)/2.0

      tmpa=-hsqrt3*(3.0*dds-4.0*ddp+ddd)
      dtmpa=-hsqrt3*(3.0*ddds-4.0*dddp+dddd)
      tmpb=hsqrt3*(dds-ddp)*2.0
      dtmpb=hsqrt3*(ddds-dddp)*2.0

      hamil=m*n*(lladmm*tmpa+tmpb)

      dhamil(1)=(dmx*n+m*dnx)*(lladmm*tmpa+tmpb)
     +  +m*n*(2.0*(l*dlx+m*dmx)*tmpa+(lladmm*dtmpa+dtmpb)*l)

      dhamil(2)=(dmy*n+m*dny)*(lladmm*tmpa+tmpb)
     +  +m*n*(2.0*(l*dly+m*dmy)*tmpa+(lladmm*dtmpa+dtmpb)*m)

      dhamil(3)=(dmz*n+m*dnz)*(lladmm*tmpa+tmpb)
     +  +m*n*(2.0*(l*dlz+m*dmz)*tmpa+(lladmm*dtmpa+dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dxy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine zxxxyy (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      tmpa=1.5*dds-2.0*ddp+0.5*ddd
      dtmpa=1.5*ddds-2.0*dddp+0.5*dddd
      tmpb=ddp-ddd
      dtmpb=dddp-dddd

      hamil=n*l*(llmimm*tmpa+tmpb)

      dhamil(1)=(dnx*l+n*dlx)*(llmimm*tmpa+tmpb)
     +    +n*l*(2.0*(l*dlx-m*dmx)*tmpa+(llmimm*dtmpa+dtmpb)*l)

      dhamil(2)=(dny*l+n*dly)*(llmimm*tmpa+tmpb)
     +    +n*l*(2.0*(l*dly-m*dmy)*tmpa+(llmimm*dtmpa+dtmpb)*m)

      dhamil(3)=(dnz*l+n*dlz)*(llmimm*tmpa+tmpb)
     +    +n*l*(2.0*(l*dlz-m*dmz)*tmpa+(llmimm*dtmpa+dtmpb)*n)

      return
      end

c
c This is a subroutine to compute the overlap integral between pz
c orbital and dx2-y2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c

      subroutine zxxyy(pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      double precision pds,dpds,pdp,dpdp,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3,temp,dtemp
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and pz replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      hsqrt3=sqrt(3.0)/2.0

      temp=hsqrt3*pds-pdp
      dtemp=hsqrt3*dpds-dpdp

      hamil=n*llmimm*temp

      dhamil(1)=(dnx*llmimm+n*2.0*(l*dlx-m*dmx))*temp+n*llmimm*dtemp*l
      dhamil(2)=(dny*llmimm+n*2.0*(l*dly-m*dmy))*temp+n*llmimm*dtemp*m
      dhamil(3)=(dnz*llmimm+n*2.0*(l*dlz-m*dmz))*temp+n*llmimm*dtemp*n

      return
      end
c
c This is a subroutine to compute the overlap integral between dzx
c orbital and dz2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine zxzz (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,hsqrt3
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      hsqrt3=sqrt(3.0)/2.0

      tmpa=-hsqrt3*(3.0*dds-4.0*ddp+ddd)
      dtmpa=-hsqrt3*(3.0*ddds-4.0*dddp+dddd)
      tmpb=hsqrt3*(dds-ddp)*2.0
      dtmpb=hsqrt3*(ddds-dddp)*2.0

      hamil=l*n*(lladmm*tmpa+tmpb)

      dhamil(1)=(dlx*n+l*dnx)*(lladmm*tmpa+tmpb)
     +  +l*n*(2.0*(l*dlx+m*dmx)*tmpa+(lladmm*dtmpa+dtmpb)*l)

      dhamil(2)=(dly*n+l*dny)*(lladmm*tmpa+tmpb)
     +  +l*n*(2.0*(l*dly+m*dmy)*tmpa+(lladmm*dtmpa+dtmpb)*m)

      dhamil(3)=(dlz*n+l*dnz)*(lladmm*tmpa+tmpb)
     +  +l*n*(2.0*(l*dlz+m*dmz)*tmpa+(lladmm*dtmpa+dtmpb)*n)

      return
      end
c
c This is a subroutine to compute the overlap integral between dxy
c orbital and dxy orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In dds the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In ddp the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c
c In ddd the first d means d orbital (angular quantum number is 2),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 2
c or -2, otherwise the ppp is 0).
c

      subroutine zz_2 (dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil)

      implicit none
      double precision dds,ddds,ddp,dddp,ddd,dddd,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      double precision tmpa,dtmpa,tmpb,dtmpb

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c

      tmpa=0.75*(3.0*dds-4.0*ddp+ddd)
      dtmpa=0.75*(3.0*ddds-4.0*dddp+dddd)
      tmpb=3.0*(dds-ddp)
      dtmpb=3.0*(ddds-dddp)

      hamil=lladmm*(lladmm*tmpa-tmpb)+dds

      dhamil(1)=2.0*(l*dlx+m*dmx)*(lladmm*tmpa-tmpb)+ddds*l
     +    +lladmm*(2.0*(l*dlx+m*dmx)*tmpa+(lladmm*dtmpa-dtmpb)*l)

      dhamil(2)=2.0*(l*dly+m*dmy)*(lladmm*tmpa-tmpb)+ddds*m
     +    +lladmm*(2.0*(l*dly+m*dmy)*tmpa+(lladmm*dtmpa-dtmpb)*m)

      dhamil(3)=2.0*(l*dlz+m*dmz)*(lladmm*tmpa-tmpb)+ddds*n
     +    +lladmm*(2.0*(l*dlz+m*dmz)*tmpa+(lladmm*dtmpa-dtmpb)*n)

      return
      end

c
c This is a subroutine to compute the overlap integral between pz
c orbital and dz2 orbital when the directional cosines of the two
c atoms is l, m, and n.
c
c In pds the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third s
c means sigma bonding (the magnetic quantum number for the two orbitals is 0,
c otherwise the sps is 0).
c
c In pdp the first p means p orbital (angular quantum number is 1),
c the second d means d orbital (angular quantum number is 2), and the third p
c means sigma bonding (the magnetic quantum number for the two orbitals is 1
c or -1, otherwise the ppp is 0).
c

      subroutine zzz(pds,dpds,pdp,dpdp,hamil,dhamil)

      implicit none
      double precision pds,dpds,pdp,dpdp,hamil,dhamil(3)
      double precision l,m,n,lladmm,llmimm,nn,sqrt3,tmp,dtmp
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following common blocks give the directional cosines and some
c combinations of the directional cosines, and the first derivatives of the
c directional cosines, which are generated by a subroutine lmndlm.
c

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

c
c The following codes uses Slater-Koster scheme to build Hamiltonian
c Interaction elements.  The reference:  J.C. Slater and G.F. Koster,
c Phys. Rev. 94, 1498-1524 (1954).  Especially Table 1, pp. 1503.
c
c If you have to calculate the overlap integrel between dx2-y2 and px replace
c pds and pdp with dps and dpp and then reverse the sign the results.
c

      sqrt3=sqrt(3.0)

      tmp=-0.5*pds+sqrt3*pdp
      dtmp=-0.5*dpds+sqrt3*dpdp

      hamil=n*(nn*pds+lladmm*tmp)

      dhamil(1)=dnx*(nn*pds+lladmm*tmp)+n*((nn*dpds+lladmm*dtmp)*l
     +    +2.0*(n*dnx*pds+(l*dlx+m*dmx)*tmp))

      dhamil(2)=dny*(nn*pds+lladmm*tmp)+n*((nn*dpds+lladmm*dtmp)*m
     +    +2.0*(n*dny*pds+(l*dly+m*dmy)*tmp))

      dhamil(3)=dnz*(nn*pds+lladmm*tmp)+n*((nn*dpds+lladmm*dtmp)*n
     +    +2.0*(n*dnz*pds+(l*dlz+m*dmz)*tmp))

      return
      end
