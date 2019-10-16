Module partpermu_mod
implicit none

contains

!decrease sort list(a,n) 
subroutine sorta(array,aind,array_size)
implicit none
integer,allocatable,intent(inout) :: array(:)
integer,intent(in) :: array_size,aind
integer :: i,j,t

do i=1,aind-2          ! each sort will be repeat aind-2 times
  do j=1,aind-2        ! sort numbers that before aind
    if(array(j)<array(j+1)) then
      t=array(j)
      array(j)=array(j+1)
      array(j+1)=t
    endif
  enddo
enddo

return
end subroutine

subroutine swapa(array,aind,bind)
implicit none
integer,allocatable,intent(inout) :: array(:)
integer,intent(in) :: aind,bind
integer :: t

t=array(aind)
array(aind)=array(bind)
array(bind)=t

return
end subroutine

subroutine permunext(array,array_size,bsize)
implicit none
integer,allocatable,intent(inout) :: array(:)
integer,intent(in) :: array_size
integer,intent(inout) :: bsize
integer :: a,b,i,j
! from left to right, find the first 1 0 combination
! index a
a=0
do i=1,array_size-1
  if(array(i)==1 .and. array(i+1)==0) then
    a=i
    b=i+1
    exit
  endif
enddo

! swap 1 0 to 0 1 and decrease sort left part
if (a>=1) then
!do j=array_size,a,-1
!  if(array(j)>array(a)) then
!    b=j
!    exit
!  endif 
!enddo
 bsize=bsize+1 
 call swapa(array,a,b)
 call sorta(array,a,array_size)
endif

return
end subroutine

! print function 
subroutine aprint(array, array_size,filehd)
implicit none
integer,allocatable,intent(in) :: array(:)
integer,intent(in) :: array_size,filehd
character(100) :: itxt,forma

write(itxt,*) array_size
forma='('//trim(adjustl(itxt))//'I6)'
write(filehd,forma) array(1:array_size)

return
end subroutine

! print function 
!subroutine readallpermu(apermu,asize,bsize)
!implicit none
!integer,allocatable,intent(inout) :: apermu(:,:)
!integer,intent(in) :: asize,bsize
!integer :: i
!logical :: alive
!
!inquire(file='Allpermu.out',exist=alive)
!if(.not.alive) stop "'Allpermu.out' is not Found! Run Again!"
!open(35,file='Allpermu.out',status='old')
!do i=1,bsize
!  read(35,*) apermu(1:asize,i)
!enddo
!read(35,*) i
!close(35)
!if(i/=bsize) then
!  write(*,"(A)") "Number of Permutations Is not Agree with 'Allpermu.out'!"
!  stop "'Allpermu.out' is Broken! Run Again!"
!endif

!return
!end subroutine

!* partly permutation
subroutine partpermu(array,array_size,bsize)
implicit none
integer,allocatable,intent(inout) :: array(:)
integer,intent(in) :: array_size
integer,intent(out) :: bsize
logical :: lreachmax
integer :: a,i,n,filehd

!output the inital sequence
filehd=55
open(filehd,file='Allpermu.out',status='replace')
call aprint(array, array_size,filehd);
n=1
lreachmax=.false.
do while(.not. lreachmax)
  call permunext(array,array_size,n)  
  call aprint(array,array_size,filehd)
  a=0
  do i=1,array_size-1
    if(array(i)==1 .and. array(i+1)==0) then
      a=i
      exit
    endif
  enddo
  if(a<1) then
    lreachmax=.true.
  endif
enddo
bsize=n
write(filehd,"(I10)") bsize
close(filehd)

return
end subroutine

end module
