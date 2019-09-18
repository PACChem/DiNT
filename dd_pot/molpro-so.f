      subroutine pot(symbol,x,y,z,pemd,gpemd,nclu,mnclu,nsurf,mnsurf)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,idum,ithis(mnclu),nsurf,mnsurf,isurf
      character*2 symbol(mnclu),cdum
      character*80 string,tmpstr
      character*21 Estr
      character*18 Gstr
      character*18 Astr
      character*9 Sstr
      double precision x(mnclu),y(mnclu),z(mnclu),v(mnsurf)
      double precision gpemd(3,mnclu,mnsurf,mnsurf),pemd(mnsurf,mnsurf)
      double precision xtmp,ytmp,ztmp,total,dum
      double precision dx(mnclu,mnsurf),dy(mnclu,mnsurf),
     &                                  dz(mnclu,mnsurf)
      double precision cfloat
      double precision autoang,so,autoev
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.211d0)

c hack
      common/tmp/tmpprint,ntmpprint
      double precision tmpprint(50)
      integer ntmpprint

c read template & write QC input file
      open(unit=7,file='qc.tmol')
      open(unit=10,file='qc.in')
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
          write(10,*)"geomtyp=xyz"
          write(10,*)"geometry"
          write(10,*)"nosym"
          write(10,*)"noorient"
          write(10,*)nclu
          write(10,*)"ANT Direct Dynamics Calculation"
          do j=1,nclu
            write(10,fmt='(a2,2x,3f20.10)')symbol(j),
     &       x(j)*autoang,y(j)*autoang,z(j)*autoang   ! default is Angstroms for geomtyp=xyz
          enddo
          write(10,*)"end"
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
            xtmp = xtmp*0.529177249/autoang        ! remove error due to different au to A conversions in Molpro and ANT
            ytmp = ytmp*0.529177249/autoang
            ztmp = ztmp*0.529177249/autoang
            do k=1,nclu
             total=dabs(xtmp-x(k))+dabs(ytmp-y(k))+dabs(ztmp-z(k))  ! check each (x,y,z) against input (both are in bohr)
             if (total.le.1.d-5) ithis(j)=k ! if the difference is small, then this is it 
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
      stop

 199  continue
      do j=1,nclu
        if (ithis(j).eq.-1) then
          write(6,*)"Problem with atom ordering in Molpro output file"
          stop
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
      stop

 299  continue

c     get gradients
      Gstr='GRADIENT FOR STATE'
      lstr=len(Gstr)
 300  read (8,fmt='(a80)',end=320) string
      do i=1,10
        if (string(i:i+lstr-1).eq.Gstr) then
          read(8,*)
          read(8,*)
          read(8,*)
          do j=1,nclu
            read(8,*)k,dx(ithis(k),isurf),
     &                 dy(ithis(k),isurf),dz(ithis(k),isurf) ! Molpro 2002.6 prints Hartrees/bohr
          enddo
          goto 399
        endif
      enddo
      goto 300

 320  write(6,*)"Couldn't find string '",gstr,"' in qc.out"
      stop

 399  continue

      ENDDO

c **************************************
c get extra data
      Sstr='LOOK HERE'
      lstr=len(Sstr)
 400  read (8,fmt='(a80)',end=420) string
      do i=1,10
        if (string(i:i+lstr-1).eq.Sstr) then
          read(8,fmt='(a80)')string
          tmpstr=string(30:45)
          ntmpprint=nint(cfloat(tmpstr))
c          print *,"AJ1",ntmpprint
          do j=1,ntmpprint
            read(8,fmt='(a80)')string
            tmpstr=string(30:45)
            tmpprint(j)=cfloat(tmpstr)
c            print *,"AJ",j,tmpprint(j)
          enddo
          goto 499
        endif
      enddo
      goto 400

 420  write(6,*)"Couldn't find string '",sstr,"' in qc.out"
      stop

 499  continue
c **************************************

      close(8)

c organize things
      do i=1,nsurf
        pemd(i,i)=v(i)
        do j=1,nclu
          gpemd(1,j,i,i)=dx(j,i)
          gpemd(2,j,i,i)=dy(j,i)
          gpemd(3,j,i,i)=dz(j,i)
        enddo
      enddo
      do i=1,nsurf
      do j=i+1,nsurf
        pemd(i,j)=5.d-5/autoev
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

      return

      end

C**********************************************************************
C CFLOAT
C**********************************************************************
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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepot
      return
      end
