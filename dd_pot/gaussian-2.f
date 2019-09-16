      subroutine pot2(symbol,x,y,z,v,dx,dy,dz,nclu,mnclu)

      implicit none
      integer i,j,k,nclu,mnclu,lstr,mult
      character*2 symbol(mnclu)
      character*80 string,tmpstr
      character*12 Estr
      character*29 Gstr
      double precision x(mnclu),y(mnclu),z(mnclu),xtmp(mnclu*3),v
      double precision dx(mnclu),dy(mnclu),dz(mnclu)
      double precision cfloat

c different multiplicities for reactants and total system
      mult=0
      do i=1,nclu
        if(symbol(i).eq.'H'.or.symbol(i).eq.' H'.or.
     &     symbol(i).eq.'h'.or.symbol(i).eq.' h'.or.
     &     symbol(i).eq.' h'.or.symbol(i).eq.' H') mult=mult+1
        if(symbol(i).eq.'C'.or.symbol(i).eq.' C'.or.
     &     symbol(i).eq.'c'.or.symbol(i).eq.' c'.or.
     &     symbol(i).eq.' c'.or.symbol(i).eq.' C') mult=mult+6
        if(symbol(i).eq.'O'.or.symbol(i).eq.' O'.or.
     &     symbol(i).eq.'o'.or.symbol(i).eq.' o'.or.
     &     symbol(i).eq.' o'.or.symbol(i).eq.' O') mult=mult+8
        if(symbol(i).eq.'He'.or.symbol(i).eq.'HE') mult=mult+2
      enddo
      mult = mod(mult,2) + 1

c read template & write QC input file
      open(unit=7,file='qc2.tg98')
      open(unit=10,file='qc.in')
      read(7,*)mult
      do i=1,100
        read(unit=7,end=100,fmt='(a80)') string
        if (string.eq."GEOMETRY") then
          write(10,*)" 0 ",mult
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
