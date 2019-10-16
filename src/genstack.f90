subroutine genstack
use permu_mod
use addmove_mod
use mod_interface
use modrdpermu

implicit none
integer,allocatable :: a(:),repindex(:)
integer,allocatable :: allpermu(:,:),finpermu(:,:),finpermu_tmp(:,:)
integer :: a_size,npermu,effnpermu
!integer :: ntype,a_size,npermu,effnpermu
!integer,allocatable :: nintype(:)
!integer :: i,j,t,stat
integer :: i,j,stat
logical :: lp,alive,alive2
!logical :: alive,lp,lptxterr,alive2
!character(100) :: forma,ni,initfile,permufile1,permufile2
character(100) :: subunit,addunit,itxt,jtxt,tofile
INTEGER*4 :: rename
character(100) :: forma,permufile1,permufile2
!character(1) :: lptxt

inquire(file='inpara.in',exist=alive)
if(.not. alive) then
  stop "'inpara.in' is not Exist!"
else
  call rdinpara(a,a_size,lp)
!write(*,*) a,a_size,lp
endif

!find Allpermu.out/Finalpermu.out then Goto write function
permufile1="Allpermu.out"
permufile2="Finalpermu.out"
inquire(file=permufile1,exist=alive)
inquire(file=permufile2,exist=alive2)
if(alive2) then
  write(*,*) "Find '",trim(permufile2),"'! Permutation Step is Skipped."
  write(*,*) "Enter Structures Constructing Steps..."
  write(*,*) 
  open(35,file=permufile2,status='old')
  read(35,*,iostat=stat) effnpermu
  if(stat/=0) then
    write(*,*) "Reading ",trim(permufile2)," Error!"
    stop
  endif
! write(*,*) "Effective P: ",effnpermu
! write(*,*) "n_size: ",a_size
  close(35)
  if(allocated(finpermu)) deallocate(finpermu)
  allocate(finpermu(a_size,effnpermu))
  call readfinpermu(finpermu,a_size,effnpermu)
  Goto 1330
endif

!calculate permutations
call dictpermu(a,a_size,npermu)

!remove repeat permutations: symmetry and periodict
if(allocated(allpermu)) deallocate(allpermu)
allocate(allpermu(a_size,npermu))
if(allocated(repindex)) deallocate(repindex)
allocate(repindex(npermu))
call readallpermu(allpermu,a_size,npermu)
forall(i=1:npermu,i>0) repindex(i)=i
call rmreppermu(allpermu,repindex,a_size,npermu,lp)
!write(*,*) 'repindex= ',repindex(1:npermu)
!effective permuation
allocate(finpermu_tmp(a_size,npermu))
effnpermu=0
do i=1,npermu
  if(repindex(i)/=0) then
    effnpermu=effnpermu+1
    finpermu_tmp(:,effnpermu)=allpermu(:,repindex(i))
  endif
enddo
allocate(finpermu(a_size,effnpermu))
if(effnpermu==npermu) then
  finpermu=allpermu
else
  do i=1,effnpermu
    finpermu(:,i)=finpermu_tmp(:,i)
  enddo
endif
deallocate(a,allpermu,finpermu_tmp,repindex)

open(55,file='Finalpermu.out',status='replace')
write(55,"(I10)") effnpermu
write(itxt,*) a_size
forma='('//trim(adjustl(itxt))//'I6)'
do j=1,effnpermu
  write(55,forma) finpermu(:,j)
enddo
close(55)

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

1330 continue

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

! enter addmove function
do i=1,effnpermu          !npermu -> effnpermu
  write(*,forma) 'Enter permutation: ',i,' ...'
  write(itxt,*) i
  do j=2,a_size
    if(j==2) then
      write(jtxt,*) finpermu(1,i)  ! allpermu -> finpermu
      subunit='POSCAR.'//trim(adjustl(jtxt))
      write(jtxt,*) finpermu(2,i)
      addunit='POSCAR.'//trim(adjustl(jtxt))
      call addmove(subunit,addunit)
      subunit='POSCAR_sub_ly'
    else
      write(jtxt,*) finpermu(j,i)
      subunit='POSCAR_sub_ly'
      addunit='POSCAR.'//trim(adjustl(jtxt))
      call addmove(subunit,addunit)
    endif
  enddo
  tofile='POSCAR_permu'//trim(adjustl(itxt))
  stat = rename(subunit,tofile)
  if ( stat .ne. 0 ) then
    write(*,forma) 'Permutation: ',i,' Structure is UN-successful!'
    stop 'Rename: Error! Check POSCAR_sub_ly and Run Angain!'
  endif
  write(*,forma) 'Permutation: ',i,' Structure is Builted Successfully!'
  write(*,*)
enddo

deallocate(finpermu)
return
end subroutine

subroutine rdinpara(a,a_size,lp)
implicit none
integer,allocatable,intent(out) :: a(:)
logical,intent(out) :: lp
integer,intent(out) :: a_size
integer :: ntype
integer,allocatable :: nintype(:)
integer :: i,j,t,stat
logical :: alive,lptxterr
character(100) :: ni,initfile
character(1) :: lptxt

lptxterr=.True.
  open(33,file='inpara.in',status='old')
  read(33,*,iostat=stat) lptxt
  if(lptxt=='T' .or. 't') then
    lp=.true.
    lptxterr=.false.
  endif
  if(lptxt=='F' .or. 'f') then
    lp=.false.
    lptxterr=.false.
  endif
  if(stat/=0 .or. lptxterr) stop 'Error in reading inpara.in: Periodic &
&                            strucuture? T/F'
  read(33,*,iostat=stat) ntype
  if(stat/=0) stop 'Error in reading inpara.in: numuber of build block types'
  allocate(nintype(ntype))
  read(33,*,iostat=stat) nintype(1:ntype)
  if(stat/=0) stop 'Error in reading inpara.in: numuber of blocks in each type'
  a_size=0
  do i=1,ntype
    do j=1,nintype(i)
      a_size=a_size+1
    enddo
  enddo
  if(a_size==1) stop 'Only one block is specified, no permutation!'
  if(allocated(a)) deallocate(a)
  allocate(a(a_size))
  t=0
  do i=1,ntype
    do j=1,nintype(i)
      t=t+1
      a(t)=i
    enddo
  enddo
if(allocated(nintype)) deallocate(nintype)
close(33)

!inquire starting blocks structure
do i=1,ntype
  write(ni,*) i
  initfile="POSCAR."//trim(adjustl(ni))
  inquire(file=initfile,exist=alive)
  if(.not. alive) then
     write(*,*) trim(adjustl(initfile))," is not exist!"
     stop
  endif
enddo

return
end subroutine
