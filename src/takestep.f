      subroutine takestep

c Take a step with the integrator.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'

      integer i,k,nv
      double precision y(mnyarray),dydx(mnyarray),yscal(mnyarray)
      double precision hnext,hdid
c bugfix, next line, Rosendo 11/29/07
      save hnext

 10   continue
      nv = 6*nat+3*nsurft
      if (methflag.eq.4) nv = nv + 2*nsurft
c put xx and pp and their derivs into a 1-D array for integration
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,y,dydx,nat,nsurft,
     & 0,methflag,pem,phase)

c Bulirsch-Stoer
      if (intflag.eq.0) then
      do i=1,mnyarray
      yscal(i)=1.d0
      enddo 
      call bsstep(y,dydx,nv,time,hstep,eps,yscal,hdid,hnext,nsurf)
      hstep=hnext

c Runge-Kutte 4th order
      elseif (intflag.eq.1) then
      call rk4(y,dydx,nv,hstep,nsurf)
      time=time+hstep

c Other
      else
      write(6,*)"Dont know intflag = ",intflag
      stop
      endif

c transform 1-D array into xx and pp
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,y,dydx,nat,nsurft,
     & 1,methflag,pem,phase)

      return

      end
