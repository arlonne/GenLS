subroutine gendefect
use partpermu_mod
use construct_mod
use modrdpermu
!use mod_interface
implicit none
integer,allocatable :: a(:),repindex(:)
integer,allocatable :: allpermu(:,:),finpermu(:,:),finpermu_tmp(:,:)
!integer,allocatable :: permub(:),sym(:,:,:)
integer,allocatable :: permub(:)
integer :: natom,nselecatom,npermu,effnpermu
integer :: i,j,t,stat
logical :: alive,alive2
character(100) :: subunit,addunit,itxt,jtxt,tofile
INTEGER*4 :: rename
character(100) :: forma,substr,permufile1

!substr='POSCAR.all'
!natom=2
!nselecatom=1

inquire(file='inpara.in',exist=alive)
if(.not. alive) then
  stop "No 'inpara.in' Found"
else
  call readinpara(natom,nselecatom,substr)
endif

!find Allpermu.out/Finalpermu.out then Goto write function
permufile1="Allpermu.out"
inquire(file=permufile1,exist=alive)
if(alive) then
  write(*,*) "Find '",trim(permufile1),"'! Permutation Step is Skipped."
  write(*,*) "Enter Structures Constructing Steps..."
  write(*,*) 
  open(35,file=permufile1,status='old',position='append')
  backspace(35)
  read(35,*,iostat=stat) npermu
  if(stat/=0) then
    write(*,*) "Reading ",permufile1,"Error!"
    stop
  endif
  close(35)
  Goto 1333
endif

!initial index array to 1:m=1 others 0: e.g., 11100000
allocate(a(natom))
a=0
a(1:nselecatom)=1

!calculate permutations
call partpermu(a,natom,npermu)
deallocate(a)

1333 continue
!remove repeat permutations: symmetry and periodict
if(allocated(allpermu)) deallocate(allpermu)
allocate(allpermu(natom,npermu))
!if(allocated(repindex)) deallocate(repindex)
!allocate(repindex(npermu))
call readallpermu(allpermu,natom,npermu)

!forall(i=1:npermu,i>0) repindex(i)=i
!call rmreppermu(allpermu,repindex,a_size,npermu,sym)
!write(*,*) 'repindex= ',repindex(1:npermu)
!effective permuation
!allocate(finpermu_tmp(a_size,npermu))
!effnpermu=0
!do i=1,npermu
!  if(repindex(i)/=0) then
!    effnpermu=effnpermu+1
!    finpermu_tmp(:,effnpermu)=allpermu(:,repindex(i))
!  endif
!enddo
effnpermu=npermu
!allocate(finpermu(natom,effnpermu))
!finpermu=allpermu
!if(effnpermu==npermu) then
!  finpermu=allpermu
!else
!  do i=1,effnpermu
!    finpermu(:,i)=finpermu_tmp(:,i)
!  enddo
!endif
!deallocate(a,allpermu)

!open(55,file='Finalpermu.out',status='replace')
!write(55,"(I10)") effnpermu
!write(itxt,*) natom
!forma='('//trim(adjustl(itxt))//'I6)'
!do j=1,effnpermu
!  write(55,forma) finpermu(:,j)
!enddo
!close(55)

if(npermu<10) then
forma='(A,I1,A)'
elseif(npermu>=10 .and. npermu<100) then
forma='(A,I2,A)'
elseif(npermu>=100 .and. npermu<1000) then
forma='(A,I3,A)'
elseif(npermu>=1000 .and. npermu<10000) then
forma='(A,I4,A)'
elseif(npermu>=10000 .and. npermu<100000) then
forma='(A,I5,A)'
elseif(npermu>=100000 .and. npermu<1000000) then
forma='(A,I6,A)'
elseif(npermu>=1000000 .and. npermu<10000000) then
forma='(A,I7,A)'
elseif(npermu>=10000000 .and. npermu<100000000) then
forma='(A,I8,A)'
else
forma='(A,I10,A)'
endif

write(*,forma) "I Have Found ",npermu," Permutations. See in 'Allpermu.out'"
if(effnpermu<npermu) then
  write(*,forma) "However, Only ",effnpermu," Permutations is Unique Due to &
                   Symmetrys. See in 'Finalpermu.out'"
endif
write(*,*)

if(effnpermu<10) then
forma='(A,I1,A)'
elseif(effnpermu>=10 .and. effnpermu<100) then
forma='(A,I2,A)'
elseif(effnpermu>=100 .and. effnpermu<1000) then
forma='(A,I3,A)'
elseif(effnpermu>=1000 .and. effnpermu<10000) then
forma='(A,I4,A)'
elseif(effnpermu>=10000 .and. effnpermu<100000) then
forma='(A,I5,A)'
elseif(effnpermu>=100000 .and. effnpermu<1000000) then
forma='(A,I6,A)'
elseif(effnpermu>=1000000 .and. effnpermu<10000000) then
forma='(A,I7,A)'
elseif(effnpermu>=10000000 .and. effnpermu<100000000) then
forma='(A,I8,A)'
else
forma='(A,I10,A)'
endif

! enter structure construct function
if(allocated(permub)) deallocate(permub)
allocate(permub(natom))
do i=1,effnpermu          !npermu -> effnpermu
! permub=finpermu(:,i)
  permub=allpermu(:,i)
  write(*,forma) 'Enter permutation: ',i,' ...'
  write(itxt,*) i
!     write(jtxt,*) finpermu(1,i)  ! allpermu -> finpermu
!     subunit='POSCAR.'//trim(adjustl(jtxt))
!     write(jtxt,*) finpermu(2,i)
!     addunit='POSCAR.'//trim(adjustl(jtxt))
  call structbuild(permub,substr)
  subunit='POSCAR_ran'
  tofile='POSCAR_permu'//trim(adjustl(itxt))
  stat = rename(subunit,tofile)
  if ( stat .ne. 0 ) then
    write(*,forma) 'Permutation: ',i,' Structure is UN-successful!'
    stop 'Rename: Error! Check POSCAR_alloy and Run Angain!'
  endif
  write(*,forma) 'Permutation: ',i,' Structure is Builted Successfully!'
  write(*,*)
enddo

!deallocate(finpermu,permub)
deallocate(allpermu,permub)
return
end subroutine

subroutine readinpara(natom,nselecatom,substr)
implicit none
integer, INTENT(OUT) :: natom,nselecatom
character(100),INTENT(OUT) :: substr
integer :: stat
logical :: alive

open(33,file='inpara.in',status='old')
read(33,*,iostat=stat) natom
if(stat/=0) stop 'Error in reading inpara.in: &
 Total numuber of atoms that may be Substituted/Deleted.'
read(33,*,iostat=stat) nselecatom
if(stat/=0) stop 'Error in reading inpara.in: numuber of defects.'
if(natom==1) stop 'Only one atom is specified, no permutation!'
if(nselecatom==natom) stop 'Substituted/Deleted atoms are All atoms, no permutation!'
read(33,*,iostat=stat) substr
if(stat/=0) stop 'Error in reading inpara.in: structure file is not specified!'
inquire(file=substr,exist=alive)
if(.not. alive) then
  write(*,*) adjustl(trim(substr))," is not Found!"
  stop
endif
!read sysmmtry operations
!read(33,*) nsym
!allocate(sym(3,3,nsym))
!do i=1,nsym
!  do j=1,3
!    read(33,*) sym(:,j,i)
!  enddo
!enddo
close(33)
return
end subroutine
