      subroutine getdvec2(nclu,u2b2,gu2b2,dvec2b2)

c Computes an effective nonadiabatic coupling vector c
c from a 2x2 diabatic matrix.

      implicit none
      include 'param.f'

      integer i1,i2,k,l,nclu
      double precision u2b2(2,2),gu2b2(3,mnat,2,2),dvec2b2(3,mnat)
      double precision u11,u12,u22,u21,cc(2,2),v1,v2,sn1,cs1,tmp

c     compute adiabatic info
c     diagonalize the 2x2
      u11 = u2b2(1,1)
      u12 = u2b2(1,2)
      u21 = u2b2(2,1)
      u22 = u2b2(2,2)
      call dlaev2(u11,u12,u22,v2,v1,cs1,sn1)
c     phase convention (this part depends on how DLAEV2 works)
      if (v2.ge.v1) then
        cc(1,1) = -sn1
        cc(2,1) = cs1
        cc(1,2) = cs1
        cc(2,2) = sn1
      else
        tmp = v1
        v1 = v2
        v2 = tmp
        cc(1,1) = cs1
        cc(2,1) = sn1
        cc(1,2) = -sn1
        cc(2,2) = cs1
      endif

c compute d
      do i1 = 1,3
      do i2 = 1,nclu
        dvec2b2(i1,i2) = 0.d0
        do k = 1, 2
        do l = 1, 2
          dvec2b2(i1,i2) = dvec2b2(i1,i2)
     &                    +cc(k,2)*cc(l,1)*gu2b2(i1,i2,k,l)
        enddo
        enddo
        if ((v1-v2) .ne. 0.0d0) then
          dvec2b2(i1,i2) = dvec2b2(i1,i2)/(v2-v1)
        else
          dvec2b2(i1,i2) = 0.0d0
        endif
      enddo
      enddo

      return

      end
