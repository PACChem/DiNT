      subroutine diamin2(xmin,val,xj,mm,nsurf,symb,rguess)

c Borrowed with modifications from NAT8.1.

c MNBRAK and BRENT are taken from Numerical Recipes.
c This routine is a combination of MNBRAK and BRENT.
c It minimizes the diatomic energy (obtained by calling
c DIAPOT) with respect to bond distance.  The value
c of the miminum-energy bond distance that is returned
c is XMIN, and the potential at the value is VAL.

c NOTE:  The code will automatically check for minima
c for all surfaces and for all three molecular 
c arrangements.  If one of these combinations is
c unbound, the code may return an error message
c and non-physical values.  This is OK because
c the code does not use these values for molecular
c arrangements that are not energetically accessible.

      implicit none
      include 'param.f'
c mnbrak variables
      integer im,nsurf
      double precision ax,bx,cx,fa,fb,fc,func,gold,glimit,tiny,guess
      parameter (gold=1.618034d0, glimit=100.0d0,tiny=1.d-20)
      double precision dum,fu,q,r,u,ulim,xj,rguess

c brent variables
      integer itmax
      double precision brent,tol,CGOLD,ZEPS
      parameter (itmax=100,CGOLD=.3819660,ZEPS=1.0e-10)
      integer iter
      double precision a,b,d,e,etemp,fv,fw,fx,p,mm(mnat),
     &                 tol1,tol2,v,w,x,xm,val,xmin
      character*2 symb(mnat)

c hard-coded
c      guess = 2.d0
      guess = rguess
      tol = 1.d-6

c      print *,"aj in diamin2"

c ***********MNBRAK***********
c initial two points are guess, and 80% of guess
      ax = guess*0.8d0
      bx = guess

c func returns y(x), where y is the function to be minimized
      call diapot2(ax,fa,xj,mm,nsurf,symb)
      call diapot2(bx,fb,xj,mm,nsurf,symb)

c make sure that fa > fb
      if (fb.gt.fa) then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif

c get third point by golden ration
      cx=bx+gold*(bx-ax)
      call diapot2(cx,fc,xj,mm,nsurf,symb)

c bracket minimum so that ax < bx < cx and fa > fb and fc > fb
 11   if (fb.ge.fc) then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/
     &     (2.0d0*dsign(max(dabs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        if ((bx-u)*(u-cx).gt.0.0d0) then
      call diapot2(u,fu,xj,mm,nsurf,symb)
          if (fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          elseif (fu.gt.fb) then
            cx=u
            fc=fu
            return
          endif
          u=cx+gold*(cx-bx)
      call diapot2(u,fu,xj,mm,nsurf,symb)
        elseif ((cx-u)*(u-ulim).gt.0.0d0) then
      call diapot2(u,fu,xj,mm,nsurf,symb)
          if (fu.lt.fc) then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            fb=fc
            fc=fu
      call diapot2(u,fu,xj,mm,nsurf,symb)
          endif
        elseif ((u-ulim)*(ulim-cx).ge.0.0d0) then
          u=ulim
      call diapot2(u,fu,xj,mm,nsurf,symb)
        else
          u=cx+gold*(cx-bx)
      call diapot2(u,fu,xj,mm,nsurf,symb)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 11
      endif
 
c uses ax,bx,cx, and the tolerance tol
c ********** BRENT ***********

      a=min(ax,cx)
      b=max(ax,cx)
      if (cx.lt.ax) a=cx
      if (cx.lt.ax) b=ax

      v=bx
      w=v
      x=v
      e=0.d0
      call diapot2(x,fx,xj,mm,nsurf,symb)
      fv=fx
      fw=fx
      do iter=1,itmax
         xm=0.5d0*(a+b)
         tol1=tol*dabs(x)+ZEPS
         tol2=2.d0*tol1
         if(dabs(x-xm).le.(tol2-.5d0*(b-a))) go to 3
         if(dabs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            if(q.gt.0.d0) p=-p
            q=dabs(q)
            etemp=e
            e=d
            if(dabs(p).ge.dabs(.5d0*q*etemp).or.p.le.q*(a-x).or.
     *           p.ge.q*(b-x)) go to 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=dsign(tol1,xm-x)
            go to 2
         endif
    1    if(x.ge.xm) then
            e=a-x
         else
            e=b-x
         endif
         d=CGOLD*e
    2    if(dabs(d).ge.tol1) then
            u=x+d
         else
            u=x+dsign(tol1,d)
         endif
      call diapot2(u,fu,xj,mm,nsurf,symb)
         if (fu.le.fx) then
            if (u.ge.x) then
               a=x
            else
               b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         else
            if(u.lt.x) then
               a=u
            else
               b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
               v=w
               fv=fw
               w=u
               fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
               v=u
               fv=fu
            endif
         endif
      enddo
      print *,"brent exceed maximum iterations"
      x = 1.d10
      fx = 1.d10

    3 xmin=x
      val=fx

      return
      END
