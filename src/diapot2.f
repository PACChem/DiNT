      subroutine diapot2(x,v,xj,mmm,nsurf,symb)

c Compute the diatomic energy.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      integer j,nsurf,i
      double precision v,xx(3,mnat),xj,erot,rmass,x,
     &  mmm(2)

      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)

      character*2 symb(mnat)

      xx(1,1) = 0.d0
      xx(2,1) = 0.d0
      xx(3,1) = 0.d0
      xx(1,2) = 0.d0
      xx(2,2) = 0.d0
      xx(3,2) = x
      call getpem(xx,2,pema,pemd,gpema,gpemd,dvec,symb)
      if (repflag.eq.1) then
        v = pemd(nsurf,nsurf)
      else
        v = pema(nsurf)
      endif

c      print *,x*autoang,v

c add rotation if j > 0
      if (xj.ge.0) then
       rmass = mmm(1)*mmm(2)/(mmm(1)+mmm(2))
       erot = 0.5d0*(xj+0.5d0)**2/(rmass*x**2)
       v = v + erot
      endif

      return
      end

