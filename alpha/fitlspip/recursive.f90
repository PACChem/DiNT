program heapppp

integer i,ia(100),n

n=4
do i=1,n
ia(i)=i
enddo

call heapp(ia,n,n)

end

recursive subroutine heapp(ia,size,n)

integer i,ia(100),n,size,ii

if (size.eq.1) then
   print *,(ia(i),i=1,n)
   return
endif

do i=1,size
  call heapp(ia,size-1,n)
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

