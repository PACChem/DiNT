      subroutine carttojac(xx,pp,mm,jq,jp,jm,ido,arr)

c if IDO = 0, transform Cartesian XX and PP to Jacobi coordinates JQ and JP
c if IDO = 1, transform Jacobi JQ and JP to Cartesian coordinates XX and PP

      implicit none
      include 'param.f'

c input
      integer ido
      double precision xx(3,mnat),pp(3,mnat),mm(mnat)

c output
      double precision jq(9),jp(9),jm(9)

c local
      integer arr,i,j,irel,ii,ij
      double precision m1,m2

c jq(1-3) are xyz for diatom
c jq(4-6) are xyz for relative translation
c jq(7-9) are xyz for overall translation
c jp are conjugate momenta

      if (arr.eq.1) then
c       C + AB
        irel = 3
        ii = 1
        ij = 2
      endif
      if (arr.eq.2) then
c       A + BC
        irel = 1
        ii = 2
        ij = 3
      endif
      if (arr.eq.3) then
c       B + CA
        irel = 2
        ii = 3
        ij = 1
      endif

      m1 = mm(irel)+mm(ii)+mm(ij)
      m2 = mm(ii)+mm(ij)

      if (ido.eq.0) then
c     get Jac from XX and PP
      do i=1,3
       jq(i) = xx(i,ij) - xx(i,ii)
       jq(i+3) = xx(i,irel) - (mm(ii)*xx(i,ii)+mm(ij)*xx(i,ij))/m2
       jq(i+6) = (mm(irel)*xx(i,irel)+mm(ij)*xx(i,ij)+mm(ii)*xx(i,ii))
     &                                       /m1
       jp(i) = -mm(ij)/m2*pp(i,ii)+mm(ii)/m2*pp(i,ij)
       jp(i+3) = m2/m1*pp(i,irel)
     &           -mm(irel)*(pp(i,ii)+pp(i,ij))/m1
       jp(i+6) = pp(i,irel)+pp(i,ii)+pp(i,ij)
       jm(i) = mm(ii)*mm(ij)/m2
       jm(i+3) = mm(irel)*m2/m1
       jm(i+6) = m1
      enddo

      else
c     get cartesians from JQ and JP
      do i=1,3
       xx(i,1) =                        m2/m1*jq(i+3)+jq(i+6)
       xx(i,2) = -mm(ij)/m2*jq(i)-mm(irel)/m1*jq(i+3)+jq(i+6)
       xx(i,3) =  mm(ii)/m2*jq(i)-mm(irel)/m1*jq(i+3)+jq(i+6)
       pp(i,1) =                  jp(i+3)+mm(irel)/m1*jp(i+6)
       pp(i,2) = -jp(i)-mm(ii)/m2*jp(i+3)+mm(ii)  /m1*jp(i+6)
       pp(i,3) =  jp(i)-mm(ij)/m2*jp(i+3)+mm(ij)  /m1*jp(i+6)
      enddo

      endif

      return
      end
