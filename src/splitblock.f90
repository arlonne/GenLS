subroutine splitblock(ablock1,nblock,iblock)
  implicit none
  character(100),intent(IN) :: ablock1
  integer,intent(OUT) :: nblock
  character(10),intent(OUT)  :: iblock(100)
  character(100)       :: ablock
  integer             :: i,toend,nblank

  nblock=1
  ablock=trim(adjustl(ablock1))
  toend=len_trim(adjustl(ablock1))
! write(*,*)
! write(*,"(A,I4)") "Length of whole words read in: ",toend
  if(toend==0) stop 'Blank string!'
  iblock=' '

  do i=1,toend
    if(i<toend) then
      if(ablock(i:i)==' ' .and. ablock(i+1:i+1)/= ' ') then
        nblock=nblock+1
      endif
    endif
    if(ablock(i:i)/=' ') then
      iblock(nblock)=trim(iblock(nblock))//ablock(i:i)
    endif
  enddo

!  write(*,*)
!  write(*,"(' Whole input were splited into',I4,' words, which are:')") nblock
!  write(*,"(100A10)") iblock(1:nblock)

  return
end subroutine
