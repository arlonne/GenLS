module modrdpermu
implicit none

contains

subroutine readallpermu(apermu,asize,bsize)
implicit none
integer,allocatable,intent(inout) :: apermu(:,:)
integer,intent(in) :: asize,bsize
integer :: i
logical :: alive

inquire(file='Allpermu.out',exist=alive)
if(.not.alive) stop "'Allpermu.out' is not Found! Run Again!"
open(35,file='Allpermu.out',status='old')
do i=1,bsize
  read(35,*) apermu(1:asize,i)
enddo
read(35,*) i
close(35)
if(i/=bsize) then
  write(*,"(A)") "Number of Permutations Is not Agree with 'Allpermu.out'!"
  stop "'Allpermu.out' is Broken! Run Again!"
endif

return
end subroutine

subroutine readfinpermu(apermu,asize,bsize)
implicit none
integer,allocatable,intent(inout) :: apermu(:,:)
integer,intent(in) :: asize,bsize
integer :: i,j
logical :: alive

inquire(file='Finalpermu.out',exist=alive)
if(.not.alive) stop "'Finalpermu.out' is not Found! Run Again!"
open(35,file='Finalpermu.out',status='old')
read(35,*) j
do i=1,bsize
  read(35,*) apermu(1:asize,i)
enddo
close(35)
if(j/=bsize) then
  write(*,"(A)") "Number of Permutations Is not Agree with 'Finalpermu.out'!"
  stop "'Finalpermu.out' is Broken! Run Again!"
endif

return
end subroutine

end module
