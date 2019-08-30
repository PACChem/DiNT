      SUBROUTINE pot(symb,x,y,z,pemd,gpemd,nat,mnat,nsurft,mnsurf)

c Wrapper to interface the NH3 potential with ant using the HE-MM-1 protocol.
c Jasper 6/2/08

      implicit none
      integer i,j,mnat,nat,nsurft,mnsurf
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     &     pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     &     xcart(12),u11,u22,u12,v1,v2,gu11(12),gu22(12),gu12(12),
     &     gv1(12),gv2(12),countn,counth

      if (nat.ne.4) then
        write(6,*)"This potential requires 4 atoms"
        stop
      endif

      countn=0
      counth=0
      do i=1,nat
        if (symb(i).eq."H".or.symb(i).eq."h") then
          counth=counth+1
          xcart(1+counth*3)=x(i)
          xcart(2+counth*3)=y(i)
          xcart(3+counth*3)=z(i)
        elseif (symb(i).eq."N".or.symb(i).eq."n") then
          xcart(1)=x(i)
          xcart(2)=y(i)
          xcart(3)=z(i)
          countn=countn+1
        else
          write(6,*)"This potential can't handle atom type: ",symb(i)
          stop
        endif
      enddo

      if (countn.eq.1.and.counth.eq.3) then
      else
        write(6,*)"This potential requires 1 N atom and 3 H atoms"
        stop
      endif

      call nh3(xcart,U11,U22,U12,V1,V2,
     1           gU11,gU22,gU12,
     2           gV1,gV2)

      pemd(1,1)=u11
      pemd(1,2)=u12
      pemd(2,1)=u12
      pemd(2,2)=u22
 
      counth=0
      do i=1,nat
        if (symb(i).eq."H".or.symb(i).eq."h") then
          counth=counth+1
          do j=1,3
          gpemd(j,i,1,1)=gu11(j+counth*3)
          gpemd(j,i,1,2)=gu12(j+counth*3)
          gpemd(j,i,2,1)=gu12(j+counth*3)
          gpemd(j,i,2,2)=gu22(j+counth*3)
          enddo
        elseif (symb(i).eq."N".or.symb(i).eq."n") then
          do j=1,3
          gpemd(j,i,1,1)=gu11(j)
          gpemd(j,i,1,2)=gu12(j)
          gpemd(j,i,2,1)=gu12(j)
          gpemd(j,i,2,2)=gu22(j)
          enddo
        endif
      enddo

      return

      end

