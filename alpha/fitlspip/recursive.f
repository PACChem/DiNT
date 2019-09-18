      program heapppp

      integer i,ia(1000),n
      integer iatom(1000),npermute

      n=5
      npermute=0
      do i=1,n
      ia(i)=i
      iatom(i)=0
      enddo
      
      call heapp(ia,n,n,iatom,npermute)
      
      end
      
      recursive subroutine heapp(ia,size,n,iia,ii)
      
      integer i,ia(1000),n,size,ii,iia(1000)

      if (size.eq.1) then
	 ii=ii+1
         print *,ii,(ia(i),i=1,n)
         return
      endif
      
      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif
      
      enddo
      
      end
      
