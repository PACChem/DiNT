       SUBROUTINE svdfit(y,sig,ndata,a,ma,u,v,w,mp,np,chisq)
       IMPLICIT NONE
       INTEGER ma,mp,ndata,np,NMAX,MMAX,i,j,ma2,idum(ndata)
       DOUBLE PRECISION chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np),
     * x(ndata),y(ndata),TOL,yy(ndata),dum(ndata)
c       EXTERNAL funcs1
c       PARAMETER (NMAX=100000,MMAX=2000,TOL=1.e-14)
       include 'param.inc'
c         NMAX is the maximum expected value of ndata; MMAX the maximum expected for ma; the default TOL val
       DOUBLE PRECISION sum,thresh,tmp,wmax,afunc(maxterm),b(maxdata)
       PARAMETER (TOL=1.e-14)
       do 12 i=1,ndata 					! Accumulate coefficients of the fitting ma
       call funcs1(i,afunc,ma)
c        if (i.eq.2) print *,i,y(i)
       tmp=1./sig(i)
       do 11 j=1,ma
       u(i,j)=afunc(j)*tmp
   11  enddo
       b(i)=y(i)*tmp
   12  enddo
       call svdcmp(u,ndata,ma,mp,np,w,v)                 ! Singular value decomposition.
       wmax=0.                                           ! Edit the singular values, given TOL from the
       do 13 j=1,ma                                      ! parameter statement, between here ...
       if(w(j).gt.wmax)wmax=w(j)
   13  enddo
       thresh=TOL*wmax
       do 14 j=1,ma
       if(w(j).lt.thresh)w(j)=0.
   14  enddo                                             ! ...and here.
       call svbksb(u,w,v,ndata,ma,mp,np,b,a)
       chisq=0.                                          ! Evaluate chi-square.
       do 16 i=1,ndata
       call funcs1(i,afunc,ma)
       sum=0.
       do 15 j=1,ma
       sum=sum+a(j)*afunc(j)
   15  enddo
       chisq=chisq+((y(i)-sum)/sig(i))**2
   16  enddo
       return
       END

