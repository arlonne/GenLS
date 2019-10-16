
module construct_mod
implicit none
  real(8) :: scale, latt_sub(3,3)
  integer,allocatable :: iontp_sub(:),iontp_ly(:)
  real(8),allocatable :: pos_sub(:,:),pos_ly(:,:)
  character(2),allocatable :: elemname_sub(:),elemname_ly(:)
  character(20) :: title_sub
  character(1),allocatable :: ismove_sub(:,:)
  character(2),allocatable :: ischg(:)
  integer :: ntype_ly,ntype_sub,nions_sub,nions_ly

contains

subroutine structbuild(buildarray,in_sub)
implicit none
character(100),intent(in) :: in_sub
integer,allocatable,intent(in) :: buildarray(:)
integer :: natom
  
  natom=size(buildarray)
  call readpos(in_sub,natom)
  call chgpos(buildarray)
  call writepos

deallocate(iontp_sub,iontp_ly,pos_sub,pos_ly,ischg)
deallocate(elemname_sub,elemname_ly,ismove_sub)

return
end subroutine

subroutine readpos(in_sub,nchgatom)
implicit none
  integer :: i
  character(100) :: in_sub,nomean
  character(300) :: posblock
  integer :: nblock,stat,nchgatom
  character(10) :: iblock(100)
  logical :: alive

!--read input structure
  ntype_sub=0
  inquire(file=in_sub,exist=alive)
  if(.not.alive) then
    write(*,"(A)") "Structure File: ",trim(adjustl(in_sub))," is not exist!"
    stop
  endif
  open(33,file=in_sub,status='old')
  read(33,'(A)') title_sub
  read(33,*) scale
  do i=1,3
    read(33,*,iostat=stat) latt_sub(1:3,i) 
  enddo
  read(33,"(A100)",iostat=stat) nomean
call splitblock(nomean,nblock,iblock)
 ntype_sub=nblock
 if(allocated(iontp_sub)) deallocate(iontp_sub)
 if(allocated(elemname_sub)) deallocate(elemname_sub)
 allocate(iontp_sub(ntype_sub),elemname_sub(ntype_sub))
 Backspace(33)
! Backspace(33)
  read(33,*) elemname_sub(1:ntype_sub)
  read(33,*) iontp_sub(1:ntype_sub)
  read(33,*) nomean
! write(*,*) "From Substrate POSCAR, I found :"
! write(*,*) "   ",ntype_sub," type ions: ",elemname_sub(1:ntype_sub)
! write(*,*) "have ",iontp_sub(1:ntype_sub)," repectively"
! write(*,*)
  nions_sub=sum(iontp_sub)
  if(allocated(pos_sub)) deallocate(pos_sub)
  allocate(pos_sub(3,nions_sub))
  if(allocated(ismove_sub)) deallocate(ismove_sub)
  allocate(ismove_sub(3,nions_sub))
  if(allocated(ischg)) deallocate(ischg)
  allocate(ischg(nions_sub))
  ischg='M'
  if(nomean(1:1)=='s' .or. nomean(1:1)=='S') then
    read(33,*)
    do i=1,nions_sub
        read(33,*) pos_sub(1:3,i),ismove_sub(1:3,i),ischg(i)
    enddo
  else
    do i=1,nions_sub
      read(33,*) pos_sub(1:3,i),ischg(i)
    enddo
    ismove_sub='F'
  endif
  close(33)
!-----------------------
!re-exam
nions_ly=0
do i=1,nions_sub
  if(ischg(i)/='M') then
    nions_ly=nions_ly+1
  endif
enddo
if(nions_ly/=nchgatom) stop 'Number of Specified Atoms in POSCAR Error!'

return
end subroutine

subroutine chgpos(pattern)
implicit none
real(8),allocatable :: pos_all(:,:),pos_tmp(:,:)
integer,allocatable,intent(in) :: pattern(:)
integer,allocatable :: iontp_tmp(:)
character(2),allocatable :: elem_all(:),elemname_tmp(:)
character(2),allocatable :: elemname_tmp2(:)
character(2),allocatable :: ischg_tmp(:)
integer :: i,j,k,m,n,natom,nspec

! all pos and elem
if(allocated(pos_all)) deallocate(pos_all)
allocate(pos_all(3,nions_sub))
if(allocated(elem_all)) deallocate(elem_all)
allocate(elem_all(nions_sub))
pos_all=pos_sub

do i=1,ntype_sub
  n=iontp_sub(i)
  if(i==1) then
    k=0
  else
    k=iontp_sub(i-1)+k
  endif
    elem_all(k+1:k+n)=elemname_sub(i)
enddo

! re-calculate the fix part of structure
n=0
do i=1,nions_sub
  if(ischg(i)=='M') then
    n=n+1
  endif
enddo
if(allocated(pos_sub)) deallocate(pos_sub)
allocate(pos_sub(3,n))
if(allocated(elemname_tmp)) deallocate(elemname_tmp)
allocate(elemname_tmp(n))
if(allocated(elemname_tmp2)) deallocate(elemname_tmp2)
allocate(elemname_tmp2(n))
j=0
do i=1,nions_sub
  if(ischg(i)=='M') then
    j=j+1
    pos_sub(:,j)=pos_all(:,i)
    elemname_tmp(j)=elem_all(i)
  endif
enddo
elemname_tmp2=elemname_tmp
natom=j
!iontp_sub() elemname_sub()
! count ion types
do i=1,natom-1
  if(elemname_tmp(i)/='no') then
    do j=i+1,natom
      if(elemname_tmp(i) == elemname_tmp(j)) then
         elemname_tmp(j)='no'
      endif
    enddo
  endif
enddo
nspec=0
do i=1,natom
  if(elemname_tmp(i)/='no') then
    nspec=nspec+1
  endif
enddo
if(allocated(elemname_sub)) deallocate(elemname_sub)
allocate(elemname_sub(nspec))
n=0
do i=1,natom
  if(elemname_tmp(i)/='no') then
    n=n+1
    elemname_sub(n)=elemname_tmp(i)
  endif
enddo
deallocate(elemname_tmp)
do i=1,nspec
  iontp_sub(i)=0
  do j=1,natom
    if(elemname_sub(i)==elemname_tmp2(j)) then
      iontp_sub(i)=iontp_sub(i)+1
    endif
  enddo
enddo
deallocate(elemname_tmp2)

!calculate the unfix part
n=0
do i=1,nions_sub
  if(ischg(i)/='M') then
    n=n+1
  endif
enddo
if(allocated(pos_tmp)) deallocate(pos_tmp)
if(allocated(elemname_tmp)) deallocate(elemname_tmp)
if(allocated(ischg_tmp)) deallocate(ischg_tmp)
allocate(pos_tmp(3,n),elemname_tmp(n))
allocate(ischg_tmp(n))
k=0
do i=1,nions_sub
  if(ischg(i)/='M') then
    k=k+1
    pos_tmp(:,k)=pos_all(:,i)
    elemname_tmp(k)=elem_all(i)
    ischg_tmp(k)=ischg(i)
  endif
enddo
if(allocated(elemname_ly)) deallocate(elemname_ly)
allocate(elemname_ly(nions_ly))
!change ly according permutation
do i=1,nions_ly
  if(pattern(i)==0) then    ! '0' means un-selected/deleted part
    elemname_ly(i)=elemname_tmp(i)
  else                      ! '1' means selected/deleted ions
    elemname_ly(i)=ischg_tmp(i)
  endif
enddo
if(allocated(elemname_tmp)) deallocate(elemname_tmp)
allocate(elemname_tmp(nions_ly))
if(allocated(elemname_tmp2)) deallocate(elemname_tmp2)
allocate(elemname_tmp2(nions_ly))
elemname_tmp=elemname_ly
elemname_tmp2=elemname_ly
do i=1,nions_ly-1
  if(elemname_tmp(i)/='no') then
    do j=i+1,nions_ly 
      if(elemname_tmp(i)==elemname_tmp(j)) then
        elemname_tmp(j)='no'
      endif
    enddo
  endif
enddo
n=0     !ntype of ly (may incluce X)
do i=1,nions_ly
  if(elemname_tmp(i)/='no') then
    n=n+1
  endif
enddo
ntype_ly=n
if(allocated(elemname_ly)) deallocate(elemname_ly)
allocate(elemname_ly(ntype_ly))
if(allocated(iontp_ly)) deallocate(iontp_ly)
allocate(iontp_ly(ntype_ly))
if(allocated(pos_ly)) deallocate(pos_ly)
allocate(pos_ly(3,nions_ly))
n=0
do i=1,nions_ly
  if(elemname_tmp(i)/='no') then
    n=n+1
    elemname_ly(n)=elemname_tmp(i)
  endif
enddo
k=0
do i=1,ntype_ly
  iontp_ly(i)=0
  do j=1,nions_ly
    if(elemname_ly(i)==elemname_tmp2(j)) then
      iontp_ly(i)=iontp_ly(i)+1
      k=k+1
      pos_ly(:,k)=pos_tmp(:,j)
    endif
  enddo
enddo

!refine pos_ly to remove 'X' part
if(allocated(iontp_tmp)) deallocate(iontp_tmp)
allocate(iontp_tmp(ntype_ly))
if(allocated(pos_tmp)) deallocate(pos_tmp)
allocate(pos_tmp(3,nions_ly))
if(allocated(elemname_tmp)) deallocate(elemname_tmp)
allocate(elemname_tmp(ntype_ly))
iontp_tmp=iontp_ly
pos_tmp=pos_ly
elemname_tmp=elemname_ly

n=0
do i=1,ntype_ly
  if(elemname_tmp(i)/='X') then
    n=n+1
  endif
enddo
if(allocated(iontp_ly)) deallocate(iontp_ly)
allocate(iontp_ly(n))
if(allocated(elemname_ly)) deallocate(elemname_ly)
allocate(elemname_ly(n))

k=0
n=0
do i=1,ntype_ly
    if(i==1) then
      j=0
    else
      j=j+iontp_tmp(i-1)
    endif
  if(elemname_tmp(i)/='X') then
    k=k+1
    iontp_ly(k)=iontp_tmp(i)
    elemname_ly(k)=elemname_tmp(i)
    do m=j+1,j+iontp_tmp(i)
      n=n+1
      pos_ly(:,n)=pos_tmp(:,m)
    enddo
  endif
enddo
deallocate(elemname_tmp,elemname_tmp2,pos_tmp,pos_all)
deallocate(iontp_tmp)

!atoms of fix and unfixed part is changed
nions_sub=natom
ntype_sub=nspec
nions_ly=n
ntype_ly=k

return
end subroutine

subroutine writepos
integer :: i
character(200) :: ntype_txt,forma
real(8) :: torl
character(1),allocatable :: moveornot(:,:)
integer :: nspec,natom

torl=1.0d-5
! merge elem and position
open(33,file='POSCAR_ran',status='replace')
!write out poscar
write(33,"(A)") 'MergedCell'
write(33,*) scale
do i=1,3
  write(33,"(3F18.11)") latt_sub(1:3,i)
enddo
nspec=ntype_sub+ntype_ly
write(ntype_txt,*) nspec 
forma='('//trim(adjustl(ntype_txt))//'A4)'
write(33,forma) elemname_sub(1:ntype_sub),elemname_ly(1:ntype_ly)
forma='('//trim(adjustl(ntype_txt))//'I6)'
write(33,forma) iontp_sub(1:ntype_sub),iontp_ly(1:ntype_ly)
write(33,"('Selective')")
write(33,"('Direct')")
natom=nions_sub+nions_ly
if(allocated(moveornot)) deallocate(moveornot)
allocate(moveornot(3,natom))
moveornot='T'
do i=1,nions_sub
  if(abs(pos_sub(1,i))<=torl) then
    moveornot(1,i)='F'
  endif
  if(abs(pos_sub(2,i))<=torl) then
    moveornot(2,i)='F'
  endif
  if(abs(pos_sub(3,i))<=torl) then
    moveornot(3,i)='F'
  endif
  write(33,"(3F18.11,X,3(X,A))") pos_sub(1:3,i),moveornot(1:3,i)
enddo
moveornot='T'
do i=1,nions_ly
  if(abs(pos_ly(1,i))<=torl) then
    moveornot(1,i)='F'
  endif
  if(abs(pos_ly(2,i))<=torl) then
    moveornot(2,i)='F'
  endif
  if(abs(pos_ly(3,i))<=torl) then
    moveornot(3,i)='F'
  endif
  write(33,"(3F18.11,X,3(X,A))") pos_ly(1:3,i),moveornot(1:3,i)
enddo
close(33)
deallocate(moveornot)

return
end subroutine

end module
