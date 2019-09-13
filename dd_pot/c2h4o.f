c **********************************************************************
c **********************************************************************
c     DD: A common interface for direct dynamics calls to the Gaussian & 
c     Molpro packages. Jasper, Dec 2007
c
c     The subroutine is specialized for diabatic surfaces
c **********************************************************************
c **********************************************************************

      subroutine pot(symb,x,y,z,pemd,gpemd,nat,mnat,nsurf,mnsurf)

      implicit none
      double precision autoang,autocmi
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)

c in/out
      integer nat,mnat,nsurf,mnsurf
      character*2 symb(mnat)
      double precision x(mnat),y(mnat),z(mnat),
     & gpemd(3,mnat,mnsurf,mnsurf),pemd(mnsurf,mnsurf)

c local
      integer mnc,ifail
      parameter(mnc=10) ! max number of QC calls per geometry
      integer ic,nc,ns(mnc),np(mnc),ntmp,ictot,i,j,k,l,ii,jj
      double precision gptmp(3,mnat,mnsurf,mnsurf),ptmp(mnsurf,mnsurf)
      integer it1,it2
      character*8 tname(mnc),tnameB(mnc),tnam

c SO coupling
      integer icc,iff
      double precision x1,x2,somax,so,del,rcent,rcom,arg
      double precision dx,dy,dz,darg,ddel,dso

      double precision xso(4),yso(4),zso(4),eso,a,b,g,eps

      double precision h12x,h12y
      common/pesh12/h12x,h12y

c -------------------
c SET THESE
c     NC = number of separate QC calls per geom
c     NS(i) = number of states for call #i
c     NP(i) = QC package to be used for call #i
c           = 1 for G03
c           = 2 for Molpro 2006
      nc    = 2
      ns(1) = 1
      np(1) = 1
      ns(2) = 1
      np(2) = 1
      tname(1) = "qc.1"
      tname(2) = "qc.3"
c -------------------

      if (nc.gt.mnc) then
          print *,"nc (",nc,") is bigger than mnc (",mnc,")" 
          stop
      endif

c zero
      do i=1,mnsurf
      do j=1,mnsurf
        pemd(i,j)=0.d0
        do k=1,3
        do l=1,mnat
          gpemd(k,l,i,j)=0.d0
        enddo
        enddo
      enddo
      enddo

c make calls
      ifail=0
      ictot=0
      do ic = 1,nc
        ntmp = ns(ic)
        tnam = tname(ic)
        if (np(ic).eq.1) then
c         call G03
          call dd_g03(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
        elseif (np(ic).eq.2) then
c         call Molpro 2006
          call system_clock(it1)
          call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
          call system_clock(it2)
c          write(83,*)dble(it2-it1)/1.d2,ifail
c         handle failures
          if (ifail.eq.0) then
c           everything OK
          elseif (ifail.eq.3) then
c           retry with numerical gradient in backup template
            ifail=0
            do i=1,mnsurf
            do j=1,mnsurf
              pemd(i,j)=0.d0
              do k=1,3
              do l=1,mnat
                gpemd(k,l,i,j)=0.d0
              enddo
              enddo
            enddo
            enddo
            tnam=tnameB(ic)
            ifail=0
            call system_clock(it1)
            call dd_m06(symb,tnam,x,y,z,ptmp,gptmp,nat,
     &                                     mnat,ntmp,mnsurf,ifail)
            call system_clock(it2)
c            write(83,*)dble(it2-it1)/1.d2,ifail
            if (ifail.ne.0) then
c             failed again
              write(6,*)"QC backup file failure"
              stop
            endif
          else
c           can't fix
            write(6,*)"QC failure"
            stop
          endif
        else
          print *,"np(",ic,") = ",np(ic),", which is not allowed"
          stop
        endif
        do i=1,ntmp
        do j=1,ntmp
          ii=i+ictot
          jj=j+ictot
          pemd(ii,jj)=ptmp(i,j)
          do k=1,3
          do l=1,nat
            gpemd(k,l,ii,jj)=gptmp(k,l,i,j)
          enddo
          enddo
        enddo
        enddo
        ictot=ictot+ntmp
      enddo

c ------
c Set diabatic coupling energy and gradient here
c ------

      pemd(1,2)=29.d0/autocmi
      pemd(2,1)=29.d0/autocmi
      h12x=pemd(1,2)
      return

      end
c **********************************************************************
c **********************************************************************





c **********************************************************************
c **********************************************************************
      subroutine prepot
      return
      end
c **********************************************************************
c **********************************************************************




c **********************************************************************
c **********************************************************************
      subroutine dd_g03(symbol,tnam,x,y,z,pemd,gpemd,nclu,
     $                                      mnclu,nsurf,mnsurf,ifail)


c INPUT
c
c SYMBOL(MNCLU) : Array of atomic symbols (H, C, etc...)
c X,Y,Z(MNCLU) :  Arrays of cartesian coordinates in bohr
c NCLU :          Number of atoms
c MNCLU :         Max number of atoms for declaring arrays
c 
c OUTPUT:
c
c V:               Potential energy in hartree
c DX,DY,DZ(MNCLU): Gradients in hartree/bohr
c IFAIL :          Flag giving info about QC failures
c                  (Not yet implemented for Gaussian)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,nsurf,mnsurf,ifail
      character*2 symbol(mnclu)
      character*80 string,tmpstr
      character*12 Estr
      character*7 Gstr
      double precision x(mnclu),y(mnclu),z(mnclu),xtmp(mnclu*3),v
      double precision dx(mnclu),dy(mnclu),dz(mnclu)
      double precision cfloat
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnclu,mnsurf,mnsurf)
      character*8 tnam

c read template & write QC input file
      open(unit=7,file=tnam)       ! template file
      open(unit=10,file='qc.in')   ! temporary input file
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
          do j=1,nclu
            write(10,fmt='(a2,2x,3f20.10)')symbol(j),x(j),y(j),z(j)
          enddo
        else
          write(10,fmt='(a80)')string
        endif
      enddo
      write(6,*)"QC TEMPLATE IS TOO LONG (> 100 lines)"
      stop
 100  close(7)
      close(10)

c do gaussian calculation
      call system('./g.x ')

c read the formatted checkpoint file
      open (8,file='Test.FChk')

c     get energy
      Estr='Total Energy'
      lstr=len(Estr)
 200  read (8,fmt='(a80)',end=220) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Estr) then
          tmpstr=string(46:80)
          v=cfloat(tmpstr)
          goto 299
        endif
      enddo
      goto 200

 220  write(6,*)"Couldn't find string '",estr,"' in Test.FChk"
      stop

 299  continue

c     get gradients
      Gstr='Cartesian Gradient'
      lstr=len(Gstr)
 300  read (8,fmt='(a80)',end=320) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Gstr) then
          j=0
          do while (j.lt.nclu*3)
            if (nclu*3-j.gt.5) then
              read(8,*)(xtmp(j+k),k=1,5)
              j=j+5
            else
c I think it goes x1,y1,z1,x2,...,zN
              read(8,*)(xtmp(j+k),k=1,nclu*3-j)
              j=nclu*3
            endif
          enddo
          goto 399
        endif
      enddo
      goto 300

 320  write(6,*)"Couldn't find string '",gstr,"' in Test.FChk"
      stop

 399  continue
      close(8)

      do i=1,nclu
        j=(i-1)*3
        dx(i)=xtmp(j+1)
        dy(i)=xtmp(j+2)
        dz(i)=xtmp(j+3)
      enddo

c organize things
      do i=1,nsurf
        pemd(i,i)=v
        do j=1,nclu
          gpemd(1,j,i,i)=dx(j)
          gpemd(2,j,i,i)=dy(j)
          gpemd(3,j,i,i)=dz(j)
        enddo
      enddo

      end
c **********************************************************************
c **********************************************************************



c **********************************************************************
c **********************************************************************
      subroutine dd_m06(symbol,tnam,x,y,z,pemd,gpemd,nclu,
     &                                      mnclu,nsurf,mnsurf,ifail)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,lstr2,idum,ithis(mnclu),nsurf,
     &   mnsurf,isurf
      integer ifail
      character*2 symbol(mnclu),cdum
      character*80 string,tmpstr
      character*21 Estr
      character*18 Gstr,Gstr2
      character*18 Astr
      double precision x(mnclu),y(mnclu),z(mnclu),v(mnsurf)
      double precision xk,dx1,dy1,dz1
      double precision gpemd(3,mnclu,mnsurf,mnsurf),pemd(mnsurf,mnsurf)
      double precision xtmp,ytmp,ztmp,total,dum
      double precision dx(mnclu,mnsurf),dy(mnclu,mnsurf),
     &                                  dz(mnclu,mnsurf)
      double precision cfloat
      double precision autoang,so,autoev
      character*8 tnam
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.211d0)

c read template & write QC input file
      open(unit=7,file=tnam)
      open(unit=10,file='qc.in')
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
c          write(10,*)"geomtyp=xyz"
c          write(10,*)"geometry"
c          write(10,*)"nosym"
c          write(10,*)"noorient"
          write(10,*)nclu
          write(10,*)"ANT Direct Dynamics Calculation"
          do j=1,nclu
            write(10,fmt='(a2,2x,3f20.10)')symbol(j),
     &       x(j)*autoang,y(j)*autoang,z(j)*autoang   ! default is Angstroms for geomtyp=xyz
          enddo
c          write(10,*)"end"
        else
          write(10,fmt='(a80)')string
        endif
      enddo
      write(6,*)"QC TEMPLATE IS TOO LONG (> 100 lines)"
      stop
 100  close(7)
      close(10)

c do molpro calculation
      call system('./m.x ')

c read the formatted checkpoint file
      open (8,file='qc.out')

c     molpro reorders the atoms, and I can't figure out how to make it stop, so...
      Astr='ATOMIC COORDINATES'
      lstr=len(Astr)
 150  read (8,fmt='(a80)',end=170) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Astr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            ithis(j)=-1
          enddo
          do j=1,nclu
            read(8,*)idum,cdum,dum,xtmp,ytmp,ztmp  ! read where Molpro reprints the reordered atoms
            do k=1,nclu
             total=dabs(xtmp-x(k))+dabs(ytmp-y(k))+dabs(ztmp-z(k))  ! check each (x,y,z) against input (both are in bohr)
             if (total.le.1.d-4) ithis(j)=k ! if the difference is small, then this is it 
c                                             (obviously this will fail for weird geometries with very small distances)
c                                             if everything works, then line number J in the molpro output corresponds
c                                             to ANT atom K
            enddo
          enddo
          goto 199
        endif
      enddo
      goto 150

 170  write(6,*)"Couldn't find string '",astr,"' in qc.out"
      ifail=1
      go to 999

 199  continue
      do j=1,nclu
        if (ithis(j).eq.-1) then
          write(6,*)"Problem with atom ordering in Molpro output file"
          ifail=2
          go to 999
        endif
      enddo

      DO ISURF=1,NSURF

c     get energy
      Estr='SETTING MOLPRO_ENERGY'
      lstr=len(Estr)
 200  read (8,fmt='(a80)',end=220) string
      do i=1,3
        if (string(i:i+lstr-1).eq.Estr) then
          tmpstr=string(28:46)
          v(isurf)=cfloat(tmpstr)
          goto 299
        endif
      enddo
      goto 200

 220  write(6,*)"Couldn't find string '",estr,"' in qc.out"
      ifail=3
      go to 999

 299  continue

c     get gradients
      Gstr='MOLGRAD'
      lstr=len(Gstr)
      Gstr2='Total Energy'
      lstr2=len(Gstr2)
 300  read (8,fmt='(a80)',end=320) string
      do i=1,10
        if (string(i:i+lstr-1).eq.Gstr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)xk,dx1,dy1,dz1
            k=nint(xk)
            dx(ithis(k),isurf)=dx1
            dy(ithis(k),isurf)=dy1
            dz(ithis(k),isurf)=dz1
          enddo
          goto 399
        elseif (string(i:i+lstr2-1).eq.Gstr2) then
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)xk,dx1,dy1,dz1
            k=nint(xk)
            dx(ithis(k),isurf)=dx1
            dy(ithis(k),isurf)=dy1
            dz(ithis(k),isurf)=dz1
          enddo
          goto 399
        endif
      enddo
      goto 300

 320  continue
      ifail=4
      go to 999

 399  continue

      ENDDO


c organize things
      do i=1,nsurf
        pemd(i,i)=v(i)
        do j=1,nclu
          gpemd(1,j,i,i)=dx(j,i)
          gpemd(2,j,i,i)=dy(j,i)
          gpemd(3,j,i,i)=dz(j,i)
        enddo
      enddo

c off-diagonal elements (set to zero here)
      do i=1,nsurf
      do j=i+1,nsurf
        pemd(i,j)=0.d0/autoev
        pemd(j,i)=pemd(i,j)
        do k=1,nclu
          gpemd(1,k,i,j)=0.d0
          gpemd(2,k,i,j)=0.d0
          gpemd(3,k,i,j)=0.d0
          gpemd(1,k,j,i)=0.d0
          gpemd(2,k,j,i)=0.d0
          gpemd(3,k,j,i)=0.d0
        enddo
      enddo
      enddo

 999  continue
      close(8)
      return

      end
C**********************************************************************
C**********************************************************************




C**********************************************************************
C**********************************************************************
C CFLOAT
C
      double precision function cfloat(string)
C
      implicit double precision(a-h,o-z)
      character*80 string,numbe
      character ch
      logical lexp,ldec

c AJ
      integer fu6
      fu6=6
C

      LEXP = .FALSE.
      LDEC = .FALSE.
      LENGTH = LEN(STRING)
      IF (LENGTH .EQ. 0) THEN
         CFLOAT = 0.0D0
         RETURN
      ENDIF
C
C     Find the first nonblank character
C
      I = 1
10    IF (STRING(I:I) .EQ. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 10
      ENDIF
C
C     If it is a blank string set function to zero
C
      IF (I .GT. LENGTH) THEN
         CFLOAT = 0.0D0
         RETURN
      ENDIF
      IBEG = I
C
C     Find the first blank character after the number
C
      I = IBEG+1
20    IF (STRING(I:I) .NE. ' ' .AND. I .LE. LENGTH) THEN
         I = I + 1
         GOTO 20
      ENDIF
      IEND = I-1
C
C     Stripe the blanks before and after the number
C
      NUMBE = STRING(IBEG:IEND)
      LENGTH = IEND - IBEG + 1
C   
C     Make sure there is no blank left
C
      IF (INDEX(NUMBE,' ') .LE. LENGTH) THEN
         WRITE(FU6,1000) STRING
         STOP 'CFLOAT 1'
      ENDIF
C
C     Find the decimal point
C
      IDEC = INDEX(NUMBE,'.')
      IF (IDEC .NE. 0) LDEC = .TRUE.
C
C     Find the exponential symbol
C
      IUE = INDEX(NUMBE,'E')
      ILE = INDEX(NUMBE,'e')
      IUD = INDEX(NUMBE,'D')
      ILD = INDEX(NUMBE,'d')
      ISUM = IUE + ILE + IUD + ILD
      IEXP = MAX0(IUE,ILE,IUD,ILD)
      IF (ISUM .GT. IEXP) THEN
         WRITE(FU6,1000) STRING
         STOP 'CFLOAT 2'
      ENDIF
      IF (IEXP .NE. 0) THEN
         LEXP = .TRUE.
      ELSE
         IEXP = LENGTH + 1
      ENDIF
C
      IF (.NOT. LDEC) IDEC = IEXP
C
C     Get the number before decimal
C
      IBEG = 2
      IF (NUMBE(1:1) .EQ. '+') THEN
         SIGN = 1.0D0
      ELSEIF(NUMBE(1:1) .EQ. '-') THEN
         SIGN = -1.0D0
      ELSE
         SIGN = 1.0D0
         IBEG = 1
      ENDIF
      IF (IBEG .EQ. IEXP) THEN
         F1 = 1.0D0
      ELSE
         F1 = 0.0D0
      ENDIF
      DO 50 I = IBEG,IDEC-1
         CH = NUMBE(I:I)
         IF (CH .GE. '0' .AND. CH .LE. '9') THEN
            N = ICHAR(CH) - ICHAR('0')
            F1 = F1 * 10.0D0 + DBLE(N)
         ELSE
            WRITE(FU6,1000) STRING
            STOP 'CFLOAT 3'
         ENDIF
50    CONTINUE
C
C     Get the number after decimal 
C
      F2 = 0.0D0
      IF (LDEC) THEN
         J = 0
         DO 60 I = IDEC+1,IEXP-1
            CH = NUMBE(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F2 = F2 * 10.0D0 + DBLE(N)
               J = J + 1
            ELSE
               WRITE(FU6,1000) STRING
               STOP 'CFLOAT 4'
            ENDIF
60       CONTINUE
         F2 = F2 / 10.0D0 ** DBLE(J)
      ENDIF
C
C    Get the exponent
C
      ESIGN = 1.0D0
      F3 = 0.0D0
      IF (LEXP) THEN 
         IBEG = IEXP + 2
         IF (NUMBE(IEXP+1:IEXP+1) .EQ. '+') THEN
            ESIGN = 1.0D0
         ELSEIF(NUMBE(IEXP+1:IEXP+1) .EQ. '-') THEN
            ESIGN = -1.0D0
         ELSE
            ESIGN = 1.0D0
            IBEG = IEXP + 1
         ENDIF
         DO 70 I = IBEG,LENGTH
            CH = NUMBE(I:I)
            IF (CH .GE. '0' .AND. CH .LE. '9') THEN
               N = ICHAR(CH) - ICHAR('0')
               F3 = F3 * 10.0D0 + DBLE(N)
            ELSE
               WRITE(FU6,1000) STRING
               STOP 'CFLOAT 5'
            ENDIF
70       CONTINUE
      ENDIF 
C
      CFLOAT = (SIGN * (F1 + F2)) * 10.0D0 ** (ESIGN*F3)
C
      RETURN
C
1000  FORMAT(/1X,'Illegal number: ',A80)
C
      END

C**********************************************************************
C**********************************************************************

      subroutine hncoso(x,y,z,eso)

      implicit double precision(a-h,o-z)
      parameter (maxcoef=600)
      dimension coef(maxcoef)

      double precision rcn,rco,rhn,rhc,rno,ahnc,anco,
     & cn(3),co(3),hn(3),dot,xco,xhn,eso,x(4),y(4),z(4)

c 1 O
c 2 C
c 3 N
c 4 H

          rcn=(z(2)-z(3))**2+(y(2)-y(3))**2+(x(2)-x(3))**2
          rcn=dsqrt(rcn)
          rco=(z(2)-z(1))**2+(y(2)-y(1))**2+(x(2)-x(1))**2
          rco=dsqrt(rco)
          rhn=(z(4)-z(3))**2+(y(4)-y(3))**2+(x(4)-x(3))**2
          rhn=dsqrt(rhn)
          rhc=(z(4)-z(2))**2+(y(4)-y(2))**2+(x(4)-x(2))**2
          rhc=dsqrt(rhc)
          rno=(z(3)-z(1))**2+(y(3)-y(1))**2+(x(3)-x(1))**2
          rno=dsqrt(rno)

          ahnc=(rhc**2-rcn**2-rhn**2)/(-2.d0*rcn*rhn)
          ahnc=dacos(ahnc)/dacos(-1.d0)*180.d0

          anco=(rno**2-rco**2-rcn**2)/(-2.d0*rcn*rco)
          anco=dacos(anco)/dacos(-1.d0)*180.d0

          cn(1)=(x(2)-x(3))/rcn  ! normalized axis
          cn(2)=(y(2)-y(3))/rcn
          cn(3)=(z(2)-z(3))/rcn

          co(1)=(x(2)-x(1))
          co(2)=(y(2)-y(1))
          co(3)=(z(2)-z(1))
          dot=co(1)*cn(1)+co(2)*cn(2)+co(3)*cn(3)
          co(1)=co(1)-dot*cn(1)
          co(2)=co(2)-dot*cn(2)
          co(3)=co(3)-dot*cn(3)
          xco=co(1)**2+co(2)**2+co(3)**2
          xco=dsqrt(xco)
          

          hn(1)=(x(3)-x(4))
          hn(2)=(y(3)-y(4))
          hn(3)=(z(3)-z(4))
          dot=hn(1)*cn(1)+hn(2)*cn(2)+hn(3)*cn(3)
          hn(1)=hn(1)-dot*cn(1)
          hn(2)=hn(2)-dot*cn(2)
          hn(3)=hn(3)-dot*cn(3)
          xhn=hn(1)**2+hn(2)**2+hn(3)**2
          xhn=dsqrt(xhn)
         
          dot=(hn(1)*co(1)+hn(2)*co(2)+hn(3)*co(3))/(xhn*xco)
          dot=dacos(dot)/dacos(-1.d0)*180.d0

          r1=rcn
          a1=ahnc
          a2=dot
          a3=anco

       coef(    1) =     0.227591360451142609235830000000D+05
       coef(    2) =    -0.822407524609554791823030000000D+06
       coef(    3) =     0.117524232108611594885590000000D+08
       coef(    4) =    -0.107103989917339354753490000000D+09
       coef(    5) =     0.776243793396884918212890000000D+09
       coef(    6) =    -0.394613286441402864456180000000D+10
       coef(    7) =     0.868247792791468620300290000000D+10
       coef(    8) =    -0.674781082654259662376720000000D+05
       coef(    9) =     0.243968848358457861468200000000D+07
       coef(   10) =    -0.396364506895566284656520000000D+08
       coef(   11) =     0.309750277340195059776310000000D+09
       coef(   12) =    -0.124354236738145971298220000000D+10
       coef(   13) =     0.397439674537977695465090000000D+10
       coef(   14) =    -0.100582409262498874664310000000D+11
       coef(   15) =     0.716704528497670544311400000000D+05
       coef(   16) =    -0.353895348454024037346240000000D+07
       coef(   17) =     0.709138295749869346618650000000D+08
       coef(   18) =    -0.432033619259702146053310000000D+09
       coef(   19) =    -0.668318205405892014503480000000D+09
       coef(   20) =     0.120886786521938629150390000000D+11
       coef(   21) =    -0.228136599310415763854980000000D+11
       coef(   22) =    -0.330062617549845526809800000000D+05
       coef(   23) =     0.202486594202605611644690000000D+07
       coef(   24) =    -0.440036458025800213217740000000D+08
       coef(   25) =     0.305583930454960227012630000000D+09
       coef(   26) =    -0.429478099847694993019100000000D+09
       coef(   27) =    -0.155510461820331764221190000000D+10
       coef(   28) =     0.139754922596545410156250000000D+10
       coef(   29) =    -0.536714155257538426667450000000D+06
       coef(   30) =     0.184168531137855201959610000000D+08
       coef(   31) =    -0.273729555694134533405300000000D+09
       coef(   32) =     0.225227360266963195800780000000D+10
       coef(   33) =    -0.100801301721240615844730000000D+11
       coef(   34) =     0.251651000685326614379880000000D+11
       coef(   35) =    -0.361182336730548248291020000000D+11
       coef(   36) =     0.893535493459939258173110000000D+06
       coef(   37) =    -0.133307182665102146565910000000D+08
       coef(   38) =     0.183436141194377899169920000000D+09
       coef(   39) =    -0.156850998074770426750180000000D+10
       coef(   40) =     0.174966922964447879791260000000D+10
       coef(   41) =     0.304102545058543968200680000000D+11
       coef(   42) =    -0.690567713248935241699220000000D+11
       coef(   43) =    -0.135720426768955599982290000000D+06
       coef(   44) =    -0.347432258743930421769620000000D+07
       coef(   45) =    -0.341834564229797601699830000000D+09
       coef(   46) =     0.452100101356552219390870000000D+10
       coef(   47) =    -0.854099543510501861572270000000D+10
       coef(   48) =    -0.760444160894592590332030000000D+11
       coef(   49) =     0.226359165236269897460940000000D+12
       coef(   50) =    -0.317317682318152219522740000000D+04
       coef(   51) =    -0.931635935909743420779710000000D+07
       coef(   52) =     0.546782763507752060890200000000D+09
       coef(   53) =    -0.667919127937339496612550000000D+10
       coef(   54) =     0.353048032706913833618160000000D+11
       coef(   55) =    -0.931661828647721862792970000000D+11
       coef(   56) =     0.118714054783224273681640000000D+12
       coef(   57) =     0.264531760259628435596820000000D+07
       coef(   58) =    -0.660238415669600665569310000000D+08
       coef(   59) =     0.834412453352769017219540000000D+09
       coef(   60) =    -0.702579317134916305541990000000D+10
       coef(   61) =     0.325578562532474822998050000000D+11
       coef(   62) =    -0.651957981211770172119140000000D+11
       coef(   63) =     0.557555687954662399291990000000D+11
       coef(   64) =    -0.659560226111940108239650000000D+07
       coef(   65) =     0.796303988425320982933040000000D+08
       coef(   66) =    -0.134001971429140710830690000000D+10
       coef(   67) =     0.146875543220876369476320000000D+11
       coef(   68) =    -0.536825355096459045410160000000D+11
       coef(   69) =    -0.752644521340117797851560000000D+11
       coef(   70) =     0.435625657016280273437500000000D+12
       coef(   71) =     0.222467334809276391752060000000D+06
       coef(   72) =     0.822946692347428947687150000000D+08
       coef(   73) =     0.207415067185633039474490000000D+10
       coef(   74) =    -0.381361155030282211303710000000D+11
       coef(   75) =     0.179405381762135314941410000000D+12
       coef(   76) =    -0.803488348406454315185550000000D+11
       coef(   77) =    -0.644649983381622680664060000000D+12
       coef(   78) =     0.185161069293389981612560000000D+07
       coef(   79) =    -0.240299829214095696806910000000D+08
       coef(   80) =    -0.266103561503818941116330000000D+10
       coef(   81) =     0.405191017215658111572270000000D+11
       coef(   82) =    -0.236500773567798156738280000000D+12
       coef(   83) =     0.598960526447369873046880000000D+12
       coef(   84) =    -0.546878868110156616210940000000D+12
       coef(   85) =    -0.327922792878083745017650000000D+07
       coef(   86) =     0.213611758828648701310160000000D+08
       coef(   87) =     0.431514151640766859054570000000D+08
       coef(   88) =     0.193561924620626926422120000000D+10
       coef(   89) =    -0.168154737274260826110840000000D+11
       coef(   90) =    -0.219718844116054840087890000000D+11
       coef(   91) =     0.193614049397848876953120000000D+12
       coef(   92) =     0.109559651699653659015890000000D+08
       coef(   93) =     0.221075211880658566951750000000D+08
       coef(   94) =     0.509263669702546596527100000000D+09
       coef(   95) =    -0.194966227767165679931640000000D+11
       coef(   96) =     0.974796071572872924804690000000D+11
       coef(   97) =     0.162804169671104888916020000000D+12
       coef(   98) =    -0.104985347177404687500000000000D+13
       coef(   99) =     0.450797851780133880674840000000D+07
       coef(  100) =    -0.477002881780247926712040000000D+09
       coef(  101) =    -0.944554270607900619506840000000D+09
       coef(  102) =     0.685264385218427581787110000000D+11
       coef(  103) =    -0.426338466374163818359380000000D+12
       coef(  104) =     0.573364610762022949218750000000D+12
       coef(  105) =     0.695554097641152099609380000000D+12
       coef(  106) =    -0.951000998927538469433780000000D+07
       coef(  107) =     0.357574825279772400856020000000D+09
       coef(  108) =     0.112042854565874099731450000000D+10
       coef(  109) =    -0.568676681413460693359380000000D+11
       coef(  110) =     0.423097089417847412109380000000D+12
       coef(  111) =    -0.118622502922962841796880000000D+13
       coef(  112) =     0.106113871127973791503910000000D+13
       coef(  113) =    -0.420466643780098005663600000000D+06
       coef(  114) =     0.139876705991836842149500000000D+08
       coef(  115) =    -0.159754655523943901062010000000D+09
       coef(  116) =     0.101348991697055387496950000000D+10
       coef(  117) =    -0.663350882041236305236820000000D+10
       coef(  118) =     0.423053412024098587036130000000D+11
       coef(  119) =    -0.111644818799077972412110000000D+12
       coef(  120) =     0.121426365579350292682650000000D+07
       coef(  121) =    -0.406135823381338268518450000000D+08
       coef(  122) =     0.583682155200336813926700000000D+09
       coef(  123) =    -0.394093660005711603164670000000D+10
       coef(  124) =     0.135041074554679241180420000000D+11
       coef(  125) =    -0.476254024054843597412110000000D+11
       coef(  126) =     0.150617971247638946533200000000D+12
       coef(  127) =    -0.120409144152997853234410000000D+07
       coef(  128) =     0.553762867328761592507360000000D+08
       coef(  129) =    -0.106764006883262324333190000000D+10
       coef(  130) =     0.623687505559921741485600000000D+10
       coef(  131) =     0.107743077797811622619630000000D+11
       coef(  132) =    -0.173817225962027618408200000000D+12
       coef(  133) =     0.308905469686999816894530000000D+12
       coef(  134) =     0.469044634133662097156050000000D+06
       coef(  135) =    -0.255180181745448298752310000000D+08
       coef(  136) =     0.471661732359754979610440000000D+09
       coef(  137) =    -0.139435115832022023200990000000D+10
       coef(  138) =    -0.177710911555648918151860000000D+11
       coef(  139) =     0.106271080754078018188480000000D+12
       coef(  140) =    -0.123830964486406539916990000000D+12
       coef(  141) =     0.100440004886206984519960000000D+08
       coef(  142) =    -0.327735271756812632083890000000D+09
       coef(  143) =     0.436662110639355373382570000000D+10
       coef(  144) =    -0.310390122775062294006350000000D+11
       coef(  145) =     0.119404096927939270019530000000D+12
       coef(  146) =    -0.274688611205116943359380000000D+12
       coef(  147) =     0.458529173304999450683590000000D+12
       coef(  148) =    -0.161003005168480202555660000000D+08
       coef(  149) =     0.177808984877083986997600000000D+09
       coef(  150) =    -0.122655133571495366096500000000D+10
       coef(  151) =     0.113956504223659553527830000000D+11
       coef(  152) =    -0.613884030233282623291020000000D+11
       coef(  153) =     0.133207007950458099365230000000D+12
       coef(  154) =    -0.426958448315176574707030000000D+12
       coef(  155) =    -0.464914643090819241479040000000D+06
       coef(  156) =     0.321641661437448024749760000000D+09
       coef(  157) =    -0.470578487537751078605650000000D+09
       coef(  158) =    -0.226886588337274932861330000000D+11
       coef(  159) =     0.726775120641153869628910000000D+11
       coef(  160) =     0.396301907919979187011720000000D+12
       coef(  161) =    -0.112277945187365576171880000000D+13
       coef(  162) =     0.283076949197583133354780000000D+07
       coef(  163) =    -0.699338374157409369945530000000D+08
       coef(  164) =    -0.296713131806315946578980000000D+10
       coef(  165) =     0.334138702392853584289550000000D+11
       coef(  166) =    -0.129340754128040893554690000000D+12
       coef(  167) =     0.332130721636885681152340000000D+12
       coef(  168) =    -0.823429019112961425781250000000D+12
       coef(  169) =    -0.505547268569631800055500000000D+08
       coef(  170) =     0.119244253844632577896120000000D+10
       coef(  171) =    -0.121995406941366939544680000000D+11
       coef(  172) =     0.791596395828834838867190000000D+11
       coef(  173) =    -0.333508566794323425292970000000D+12
       coef(  174) =     0.798816119324198486328120000000D+12
       coef(  175) =    -0.128231488388710009765620000000D+13
       coef(  176) =     0.121352011183028459548950000000D+09
       coef(  177) =    -0.103555119182126235961910000000D+10
       coef(  178) =     0.669337627015719127655030000000D+10
       coef(  179) =    -0.657250512705591583251950000000D+11
       coef(  180) =     0.456231196159002075195310000000D+12
       coef(  181) =    -0.885255033012258911132810000000D+12
       coef(  182) =     0.980430991067901367187500000000D+12
       coef(  183) =     0.109381784490384086966510000000D+08
       coef(  184) =    -0.321678922942690753936770000000D+10
       coef(  185) =     0.166433447282898559570310000000D+11
       coef(  186) =     0.578320492986435394287110000000D+11
       coef(  187) =    -0.596716832118312500000000000000D+12
       coef(  188) =     0.171429116245019653320310000000D+12
       coef(  189) =     0.206037271897242919921880000000D+13
       coef(  190) =    -0.465305116409961357712750000000D+08
       coef(  191) =     0.176561493899216222763060000000D+10
       coef(  192) =     0.476281187275918006896970000000D+10
       coef(  193) =    -0.141586778784539764404300000000D+12
       coef(  194) =     0.675515357752444335937500000000D+12
       coef(  195) =    -0.118048146455915258789060000000D+13
       coef(  196) =     0.972179364375206787109380000000D+12
       coef(  197) =     0.657266461485663205385210000000D+08
       coef(  198) =    -0.449542469673758029937740000000D+09
       coef(  199) =    -0.450023016903790092468260000000D+10
       coef(  200) =     0.227030168227515640258790000000D+11
       coef(  201) =     0.183805588929474639892580000000D+12
       coef(  202) =    -0.779550683757423217773440000000D+12
       coef(  203) =     0.637640836436785156250000000000D+12
       coef(  204) =    -0.209980411776942163705830000000D+09
       coef(  205) =    -0.950912043150328159332280000000D+09
       coef(  206) =     0.253746861028542900085450000000D+11
       coef(  207) =    -0.895665906122153320312500000000D+11
       coef(  208) =    -0.409557956613270751953120000000D+12
       coef(  209) =     0.121965879997559838867190000000D+13
       coef(  210) =    -0.258280418820494750976560000000D+12
       coef(  211) =    -0.923228355912820398807530000000D+08
       coef(  212) =     0.113402646397943572998050000000D+11
       coef(  213) =    -0.883807534548449707031250000000D+11
       coef(  214) =     0.167574789515394683837890000000D+12
       coef(  215) =     0.708924405354067260742190000000D+12
       coef(  216) =    -0.106759310819517333984380000000D+13
       coef(  217) =    -0.208951359701558691406250000000D+13
       coef(  218) =     0.183250600307854920625690000000D+09
       coef(  219) =    -0.838385182154230976104740000000D+10
       coef(  220) =     0.557441039953632125854490000000D+11
       coef(  221) =    -0.114856757170365463256840000000D+12
       coef(  222) =    -0.306843578847260681152340000000D+12
       coef(  223) =     0.114433247512769799804690000000D+13
       coef(  224) =    -0.155472808526854003906250000000D+12
       coef(  225) =     0.288611828659026557579640000000D+07
       coef(  226) =    -0.906371181327782869338990000000D+08
       coef(  227) =     0.830897481384882330894470000000D+09
       coef(  228) =    -0.232851006135761260986330000000D+10
       coef(  229) =     0.719586701379900550842290000000D+10
       coef(  230) =    -0.130439660625315597534180000000D+12
       coef(  231) =     0.496076887340369567871090000000D+12
       coef(  232) =    -0.831875202398184873163700000000D+07
       coef(  233) =     0.268437545957465946674350000000D+09
       coef(  234) =    -0.355528124535090541839600000000D+10
       coef(  235) =     0.215630368979825706481930000000D+11
       coef(  236) =    -0.649326830201865844726560000000D+11
       coef(  237) =     0.237158367118463073730470000000D+12
       coef(  238) =    -0.842741882845394165039060000000D+12
       coef(  239) =     0.839087630767405033111570000000D+07
       coef(  240) =    -0.384124429173964500427250000000D+09
       coef(  241) =     0.762448711957059478759770000000D+10
       coef(  242) =    -0.556397445948636627197270000000D+11
       coef(  243) =     0.946807487836897888183590000000D+11
       coef(  244) =     0.370362572697788635253910000000D+12
       coef(  245) =    -0.795268936979929565429690000000D+12
       coef(  246) =    -0.317527138294442417100070000000D+07
       coef(  247) =     0.167653183590920537710190000000D+09
       coef(  248) =    -0.310032079036917543411250000000D+10
       coef(  249) =     0.117776611918938827514650000000D+11
       coef(  250) =     0.711675977007880706787110000000D+11
       coef(  251) =    -0.455352687276838684082030000000D+12
       coef(  252) =     0.411885723877481689453120000000D+12
       coef(  253) =    -0.697533184108883738517760000000D+08
       coef(  254) =     0.222340609114949607849120000000D+10
       coef(  255) =    -0.272362579304672584533690000000D+11
       coef(  256) =     0.161431823421851623535160000000D+12
       coef(  257) =    -0.424099257722075805664060000000D+12
       coef(  258) =     0.421944576301376586914060000000D+12
       coef(  259) =    -0.892538355661591186523440000000D+12
       coef(  260) =     0.115223280610326334834100000000D+09
       coef(  261) =    -0.132579839263484287261960000000D+10
       coef(  262) =     0.460192034714242935180660000000D+10
       coef(  263) =     0.867446946878662109375000000000D+08
       coef(  264) =     0.943257434373796844482420000000D+10
       coef(  265) =    -0.374114576620662612915040000000D+11
       coef(  266) =     0.172691044109432446289060000000D+13
       coef(  267) =    -0.219916290589901246130470000000D+07
       coef(  268) =    -0.221255515785266447067260000000D+10
       coef(  269) =     0.114056946325572605133060000000D+11
       coef(  270) =     0.594264130642853469848630000000D+11
       coef(  271) =    -0.503654541072154724121090000000D+12
       coef(  272) =    -0.171746100710094726562500000000D+12
       coef(  273) =     0.207661652271763330078120000000D+13
       coef(  274) =    -0.185738640401843190193180000000D+08
       coef(  275) =     0.674221426448468923568730000000D+09
       coef(  276) =     0.101726417877162609100340000000D+11
       coef(  277) =    -0.834914031263602294921880000000D+11
       coef(  278) =     0.639942469064036636352540000000D+11
       coef(  279) =    -0.686622465398386077880860000000D+11
       coef(  280) =     0.273241671138644824218750000000D+13
       coef(  281) =     0.363358821271631836891170000000D+09
       coef(  282) =    -0.872566707162492752075200000000D+10
       coef(  283) =     0.809279914302359771728520000000D+11
       coef(  284) =    -0.365528221959974243164060000000D+12
       coef(  285) =     0.706570504866792114257810000000D+12
       coef(  286) =     0.274150783924787536621090000000D+12
       coef(  287) =     0.588149471852296264648440000000D+12
       coef(  288) =    -0.905712465868008494377140000000D+09
       coef(  289) =     0.102476707858646564483640000000D+11
       coef(  290) =    -0.594300322597367477416990000000D+11
       coef(  291) =     0.137816889706669494628910000000D+12
       coef(  292) =    -0.505223293926401184082030000000D+12
       coef(  293) =    -0.520488653084554687500000000000D+12
       coef(  294) =    -0.503411448372449462890620000000D+12
       coef(  295) =     0.426663254240191429853440000000D+08
       coef(  296) =     0.177638460592382583618160000000D+11
       coef(  297) =    -0.113058784247958740234380000000D+12
       coef(  298) =     0.244916881002108489990230000000D+12
       coef(  299) =     0.580987122581054687500000000000D+12
       coef(  300) =    -0.141266155828030670166020000000D+12
       coef(  301) =    -0.825246373970281250000000000000D+12
       coef(  302) =     0.254370180178004413843150000000D+09
       coef(  303) =    -0.100444210923315353393550000000D+11
       coef(  304) =    -0.244632911938973960876460000000D+11
       coef(  305) =     0.446685933815136230468750000000D+12
       coef(  306) =    -0.968386050954097900390620000000D+12
       coef(  307) =    -0.419484028202063232421880000000D+12
       coef(  308) =    -0.727933531311538330078120000000D+12
       coef(  309) =    -0.502247551731739044189450000000D+09
       coef(  310) =     0.521640687592238426208500000000D+10
       coef(  311) =     0.167678955695976638793950000000D+11
       coef(  312) =    -0.309078025389781005859380000000D+12
       coef(  313) =     0.515933697167967712402340000000D+12
       coef(  314) =    -0.208156183765656433105470000000D+12
       coef(  315) =     0.106440734174480065917970000000D+13
       coef(  316) =     0.167854465310113167762760000000D+10
       coef(  317) =    -0.479107296416229343414310000000D+10
       coef(  318) =    -0.685503956260692520141600000000D+11
       coef(  319) =     0.744491810703077270507810000000D+12
       coef(  320) =    -0.940476937456233276367190000000D+12
       coef(  321) =     0.673620387511822631835940000000D+12
       coef(  322) =    -0.446100409649051879882810000000D+12
       coef(  323) =     0.190581070211094796657560000000D+09
       coef(  324) =    -0.577095990405402984619140000000D+11
       coef(  325) =     0.419715794360809326171880000000D+12
       coef(  326) =    -0.133974953065164721679690000000D+13
       coef(  327) =     0.625397164142054687500000000000D+12
       coef(  328) =     0.730268604870091552734380000000D+11
       coef(  329) =    -0.159227518918631689453120000000D+13
       coef(  330) =    -0.985140522243207812309270000000D+09
       coef(  331) =     0.453381844799718017578120000000D+11
       coef(  332) =    -0.257030069420148773193360000000D+12
       coef(  333) =     0.769170860027774658203120000000D+12
       coef(  334) =    -0.130892138279145361328120000000D+13
       coef(  335) =     0.105333129367731494140620000000D+13
       coef(  336) =    -0.551596546769862182617190000000D+12
       coef(  337) =    -0.848305274708347581326960000000D+07
       coef(  338) =     0.251468252374846994876860000000D+09
       coef(  339) =    -0.173675823619745302200320000000D+10
       coef(  340) =    -0.507961224753224945068360000000D+10
       coef(  341) =     0.754594908111085968017580000000D+11
       coef(  342) =    -0.140557705786494255065920000000D+11
       coef(  343) =    -0.790248064773350708007810000000D+12
       coef(  344) =     0.242150814095731973648070000000D+08
       coef(  345) =    -0.734669983546177983284000000000D+09
       coef(  346) =     0.835750260607081031799320000000D+10
       coef(  347) =    -0.350179262176918411254880000000D+11
       coef(  348) =     0.697364947854302787780760000000D+10
       coef(  349) =    -0.202468825220369567871090000000D+11
       coef(  350) =     0.128476554400535034179690000000D+13
       coef(  351) =    -0.248903333912154696881770000000D+08
       coef(  352) =     0.111331018789930343627930000000D+10
       coef(  353) =    -0.221371579327320327758790000000D+11
       coef(  354) =     0.178199930978403686523440000000D+12
       coef(  355) =    -0.547806550708060607910160000000D+12
       coef(  356) =     0.307093551838130859375000000000D+12
       coef(  357) =     0.151085723770688842773440000000D+12
       coef(  358) =     0.953410896866091340780260000000D+07
       coef(  359) =    -0.493937171277942180633540000000D+09
       coef(  360) =     0.934946310983517265319820000000D+10
       coef(  361) =    -0.481754418140493240356450000000D+11
       coef(  362) =    -0.282669605659199218750000000000D+11
       coef(  363) =     0.492127458752572387695310000000D+12
       coef(  364) =     0.431921721853302383422850000000D+11
       coef(  365) =     0.209216832049640595912930000000D+09
       coef(  366) =    -0.661365920367007732391360000000D+10
       coef(  367) =     0.775232849596382141113280000000D+11
       coef(  368) =    -0.413759554594708129882810000000D+12
       coef(  369) =     0.819562154741207275390620000000D+12
       coef(  370) =     0.184864798390042266845700000000D+12
       coef(  371) =     0.368573756322556884765620000000D+12
       coef(  372) =    -0.356297762272236883640290000000D+09
       coef(  373) =     0.434560010781825733184810000000D+10
       coef(  374) =    -0.705734627370188522338870000000D+10
       coef(  375) =    -0.160110214129443847656250000000D+12
       coef(  376) =     0.770512978090154907226560000000D+12
       coef(  377) =    -0.157147003858765820312500000000D+13
       coef(  378) =    -0.189119605030888378906250000000D+13
       coef(  379) =     0.291490672176066040992740000000D+08
       coef(  380) =     0.616714481886519718170170000000D+10
       coef(  381) =    -0.487525196942003326416020000000D+11
       coef(  382) =     0.551643520994449844360350000000D+11
       coef(  383) =     0.103157951825925109863280000000D+13
       coef(  384) =    -0.165930483618140747070310000000D+13
       coef(  385) =    -0.245315228945315869140620000000D+13
       coef(  386) =     0.451756415105110555887220000000D+08
       coef(  387) =    -0.198029963178930616378780000000D+10
       coef(  388) =    -0.186048581616471366882320000000D+11
       coef(  389) =     0.702936226444275970458980000000D+11
       coef(  390) =     0.694192755232234863281250000000D+12
       coef(  391) =    -0.232503849813204736328120000000D+13
       coef(  392) =    -0.256744431389774609375000000000D+13
       coef(  393) =    -0.112758216215065956115720000000D+10
       coef(  394) =     0.282687046679627761840820000000D+11
       coef(  395) =    -0.266361488293147430419920000000D+12
       coef(  396) =     0.113410871730770898437500000000D+13
       coef(  397) =    -0.183480135479166430664060000000D+13
       coef(  398) =    -0.183224393956406054687500000000D+13
       coef(  399) =     0.338004254740234130859380000000D+12
       coef(  400) =     0.291510809384830856323240000000D+10
       coef(  401) =    -0.417912413063760833740230000000D+11
       coef(  402) =     0.298635800524513488769530000000D+12
       coef(  403) =    -0.689412597934034667968750000000D+12
       coef(  404) =     0.129252642201776464843750000000D+13
       coef(  405) =     0.105080181927923950195310000000D+13
       coef(  406) =    -0.558911373627964721679690000000D+12
       coef(  407) =    -0.523517104840029120445250000000D+09
       coef(  408) =    -0.352529972369605560302730000000D+11
       coef(  409) =     0.181247237454072296142580000000D+12
       coef(  410) =    -0.790649624390348266601560000000D+12
       coef(  411) =     0.682838703636981567382810000000D+12
       coef(  412) =     0.477362470755019897460940000000D+12
       coef(  413) =    -0.139919712089200537109380000000D+13
       coef(  414) =    -0.537870163395693421363830000000D+09
       coef(  415) =     0.208181873337961044311520000000D+11
       coef(  416) =     0.156708836352592803955080000000D+12
       coef(  417) =    -0.137437727173534887695310000000D+13
       coef(  418) =     0.194918769922311547851560000000D+13
       coef(  419) =     0.203981520143190991210940000000D+13
       coef(  420) =    -0.630825565378871826171880000000D+12
       coef(  421) =     0.163831770294090938568120000000D+10
       coef(  422) =    -0.229488673743924255371090000000D+11
       coef(  423) =     0.467790251138009262084960000000D+11
       coef(  424) =     0.506579998533538085937500000000D+12
       coef(  425) =    -0.869082674046520019531250000000D+12
       coef(  426) =    -0.101613296900951660156250000000D+13
       coef(  427) =     0.374602584795738952636720000000D+12
       coef(  428) =    -0.567359495063316154479980000000D+10
       coef(  429) =     0.508143170087616729736330000000D+11
       coef(  430) =    -0.265709412787918548583980000000D+12
       coef(  431) =    -0.236167724936382629394530000000D+12
       coef(  432) =     0.550432956592742614746090000000D+11
       coef(  433) =     0.857115611384540771484380000000D+12
       coef(  434) =    -0.172245848667363403320310000000D+12
       coef(  435) =     0.745814626153419375419620000000D+09
       coef(  436) =     0.106009638754919647216800000000D+12
       coef(  437) =    -0.420976406083545654296880000000D+12
       coef(  438) =     0.674400992848071533203120000000D+12
       coef(  439) =     0.129242010152720922851560000000D+13
       coef(  440) =     0.799060792280504882812500000000D+12
       coef(  441) =    -0.735336921725433837890620000000D+12
       coef(  442) =     0.216366217680346012115480000000D+10
       coef(  443) =    -0.980028331401249084472660000000D+11
       coef(  444) =     0.280151468488865844726560000000D+12
       coef(  445) =    -0.831576241558233032226560000000D+11
       coef(  446) =     0.648670838291647338867190000000D+11
       coef(  447) =     0.141156357638805834960940000000D+13
       coef(  448) =    -0.165983247774429199218750000000D+12
       coef(  449) =     0.890503340161270089447500000000D+07
       coef(  450) =    -0.243778339033547520637510000000D+09
       coef(  451) =     0.914627591172404527664180000000D+09
       coef(  452) =     0.202455123856228065490720000000D+11
       coef(  453) =    -0.197189978598620880126950000000D+12
       coef(  454) =     0.484604078042398315429690000000D+12
       coef(  455) =     0.707853112052564544677730000000D+11
       coef(  456) =    -0.250783863950135223567490000000D+08
       coef(  457) =     0.681752331268677473068240000000D+09
       coef(  458) =    -0.536920612623529243469240000000D+10
       coef(  459) =    -0.126504745345291080474850000000D+11
       coef(  460) =     0.346443288284323852539060000000D+12
       coef(  461) =    -0.132335673568719580078120000000D+13
       coef(  462) =     0.852080071740771972656250000000D+12
       coef(  463) =     0.259764572568356245756150000000D+08
       coef(  464) =    -0.109484649650756096839900000000D+10
       coef(  465) =     0.208482222858830871582030000000D+11
       coef(  466) =    -0.165362497868251159667970000000D+12
       coef(  467) =     0.528798221553710571289060000000D+12
       coef(  468) =    -0.395598002773572692871090000000D+12
       coef(  469) =    -0.822206818467578430175780000000D+11
       coef(  470) =    -0.101138109413729608058930000000D+08
       coef(  471) =     0.503287002091175138950350000000D+09
       coef(  472) =    -0.943456768499755287170410000000D+10
       coef(  473) =     0.527591118514281539916990000000D+11
       coef(  474) =    -0.480266027519365921020510000000D+11
       coef(  475) =    -0.111380919952261138916020000000D+12
       coef(  476) =    -0.591606516888461669921880000000D+12
       coef(  477) =    -0.226143162546206325292590000000D+09
       coef(  478) =     0.708257962555968570709230000000D+10
       coef(  479) =    -0.799683849244004516601560000000D+11
       coef(  480) =     0.393325974100023315429690000000D+12
       coef(  481) =    -0.620045940847812866210940000000D+12
       coef(  482) =    -0.799804140000675170898440000000D+12
       coef(  483) =     0.516595068199702514648440000000D+12
       coef(  484) =     0.395622323175564587116240000000D+09
       coef(  485) =    -0.486039851719141769409180000000D+10
       coef(  486) =    -0.304327804732994508743290000000D+10
       coef(  487) =     0.345589680578906188964840000000D+12
       coef(  488) =    -0.149276483183504541015620000000D+13
       coef(  489) =     0.237568155441847070312500000000D+13
       coef(  490) =     0.518182709384265869140620000000D+12
       coef(  491) =    -0.507285589459194093942640000000D+08
       coef(  492) =    -0.656032677755815982818600000000D+10
       coef(  493) =     0.731276712924544067382810000000D+11
       coef(  494) =    -0.323317728861003479003910000000D+12
       coef(  495) =    -0.566243189660601318359380000000D+12
       coef(  496) =     0.354163395243913525390620000000D+13
       coef(  497) =     0.481722623635561706542970000000D+12
       coef(  498) =    -0.401811062035230547189710000000D+08
       coef(  499) =     0.219805055402766609191890000000D+10
       coef(  500) =     0.941045193862219619750980000000D+10
       coef(  501) =     0.441981599985169754028320000000D+11
       coef(  502) =    -0.103791282240916748046880000000D+13
       coef(  503) =     0.275317244506686425781250000000D+13
       coef(  504) =     0.846544226264884796142580000000D+11
       coef(  505) =     0.125158300737453794479370000000D+10
       coef(  506) =    -0.322608116777000579833980000000D+11
       coef(  507) =     0.306538686488593933105470000000D+12
       coef(  508) =    -0.127683142531066552734380000000D+13
       coef(  509) =     0.216192613307625903320310000000D+13
       coef(  510) =     0.704021308952226562500000000000D+12
       coef(  511) =     0.970665829699770751953120000000D+12
       coef(  512) =    -0.333215637489406776428220000000D+10
       coef(  513) =     0.539275639818861160278320000000D+11
       coef(  514) =    -0.401223278055239807128910000000D+12
       coef(  515) =     0.772561823571442749023440000000D+12
       coef(  516) =    -0.924821211566028442382810000000D+12
       coef(  517) =     0.102339911651577563476560000000D+13
       coef(  518) =     0.373880515577537078857420000000D+11
       coef(  519) =     0.896806152629122972488400000000D+09
       coef(  520) =     0.246038172313322868347170000000D+11
       coef(  521) =    -0.105088403481351943969730000000D+12
       coef(  522) =     0.123825670457690795898440000000D+13
       coef(  523) =    -0.346308418862348583984380000000D+13
       coef(  524) =    -0.279081791933069641113280000000D+12
       coef(  525) =    -0.781969099113280395507810000000D+12
       coef(  526) =     0.412279719129535079002380000000D+09
       coef(  527) =    -0.163788933991197338104250000000D+11
       coef(  528) =    -0.208658689040943054199220000000D+12
       coef(  529) =     0.135978979192592382812500000000D+13
       coef(  530) =    -0.175259555948422558593750000000D+13
       coef(  531) =     0.107509286727178369140620000000D+13
       coef(  532) =    -0.212679444603247863769530000000D+12
       coef(  533) =    -0.187737768614117932319640000000D+10
       coef(  534) =     0.303791153002554168701170000000D+11
       coef(  535) =    -0.117314501879527633666990000000D+12
       coef(  536) =    -0.346308111962231506347660000000D+12
       coef(  537) =     0.111461193343622924804690000000D+13
       coef(  538) =     0.285206290999530944824220000000D+12
       coef(  539) =     0.468606110420266113281250000000D+12
       coef(  540) =     0.666705105631488609313960000000D+10
       coef(  541) =    -0.827727384128704681396480000000D+11
       coef(  542) =     0.603421910728106079101560000000D+12
       coef(  543) =    -0.608570338811336547851560000000D+12
       coef(  544) =    -0.846377825622058837890620000000D+12
       coef(  545) =     0.327992464219600097656250000000D+12
       coef(  546) =    -0.443675983818184814453120000000D+11
       coef(  547) =    -0.184048674801085448265080000000D+10
       coef(  548) =    -0.649759230459321899414060000000D+11
       coef(  549) =    -0.106178470410757522583010000000D+12
       coef(  550) =     0.890042815624353271484380000000D+12
       coef(  551) =    -0.957476084068588256835940000000D+12
       coef(  552) =    -0.402367979388060607910160000000D+11
       coef(  553) =    -0.420316173532886474609380000000D+12
       coef(  554) =    -0.178205200140485405921940000000D+10
       coef(  555) =     0.814424165951492614746090000000D+11
       coef(  556) =    -0.421360519367182693481450000000D+11
       coef(  557) =    -0.445199696144614257812500000000D+12
       coef(  558) =    -0.162013900280427832031250000000D+13
       coef(  559) =     0.261951693298928985595700000000D+12
       coef(  560) =    -0.157614380519611450195310000000D+12

      norder1=4 ! nc
      norder2=3 ! hnc
      norder3=3 ! dih
      norder4=6 ! nco

      r1=r1
      r2=a1
      r3=a2
      r4=a3

      b1=1.d0
      b2=50.d0
      b3=100.d0

      e1=dexp(-r1/b1)
      e2=dexp(-r2/b2)
      e3=dexp(-r3/b3)
      e4=dexp(-r4/b2)

c      print "e",e1,e2,e3,e4

      ii=0
      eso=0.d0
      do i1=0,norder1
      x1=e1**i1
      do i2=0,norder2
      x2=x1*(e2**i2)
      do i3=0,norder3
      x3=x2*(e3**i3)
      do i4=0,norder4
      x4=x3*(e4**i4)
      ii=ii+1
      basis=1.d0
      if (i1+i2+i3+i4.ne.0) basis=x4
      eso=eso+basis*coef(ii)
c      print *,ii,eso,basis,coef(ii)
      enddo
      enddo
      enddo
      enddo

c      print *,r1,a1,a2,a3,eso

      return 
      end

***************************************************
