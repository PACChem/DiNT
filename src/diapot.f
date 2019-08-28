      subroutine diapot(x,im,v,xj,mmm,nsurf)

c Compute the diatomic energy for molecular arrangement
c number IM. Specialized for atom-diatom initial conditions.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      integer im,j,nsurf,i
      double precision v,xx(3,mnat),r(3),xj,erot,rmass,x,
     &  mmm(mnat)

      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)

c      r(1) = 100.0d0
c      r(2) = 100.0d0
c      r(3) = 100.0d0
      r(1) = 7.0d0
      r(2) = 7.0d0
      r(3) = 7.0d0
c      r(1) = 4.0d0   ! IMLS
c      r(2) = 4.0d0
c      r(3) = 4.0d0
      r(im) = x
      call rtox(xx,r,0)
      call getpem(xx,3,pema,pemd,gpema,gpemd,dvec,symbol)
      if (repflag.eq.1) then
        v = pemd(nsurf,nsurf)
      else
        v = pema(nsurf)
      endif

c      print *,x*autoang,v*autoev

c add rotation if j > 0
      if (xj.ge.0) then
      if (im.eq.1) rmass = mmm(1)*mmm(2)/(mmm(1)+mmm(2))
      if (im.eq.2) rmass = mmm(2)*mmm(3)/(mmm(2)+mmm(3))
      if (im.eq.3) rmass = mmm(3)*mmm(1)/(mmm(3)+mmm(1))
      erot = 0.5d0*(xj+0.5d0)**2/(rmass*x**2)
      v = v + erot
      endif

      return
      end

