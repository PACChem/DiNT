      subroutine derivs(x,yn,yout,nv,nsurf)

c This subroutine (along with XPTOY) puts together an single array 
c and its derivatives for integration.  It is called by the
c integrators whenever they need derivative information.  YN
c is the array of variables to be integrated, and YOUT is the
c derivatives.

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer nsurf,nv
      double precision xx(3,mnat),pp(3,mnat),gv(3,mnat),v,
     & cre(2*mnsurf),cim(2*mnsurf),gcre(2*mnsurf),gcim(2*mnsurf),
     & pem(mnsurf),phase(mnsurf),dvec(3,mnat,mnsurf,mnsurf)
      double precision yn(mnyarray),yout(mnyarray),phop(mnsurf)
      double precision x,dmag,gpem(3,mnat,mnsurf)

c     get xx from yn
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,yn,yout,nat,nsurft,
     & 1,methflag,pem,phase)
c     get gradients
      call getgrad(xx,pp,nsurf,v,cre,cim,gv,gcre,gcim,nat,phop,dmag,
     & dvec,pem,gpem,phase)
c     transform to 1-D
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,yn,yout,nat,nsurft,
     & 0,methflag,pem,phase)

      return
      end
