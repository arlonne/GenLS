!first version done in 2016.Jul 14
!by arlonne@ujs
! plan:
!second version: include .cif output
!third version: sort elements in output files
!fourth version: combine functions: writepos and writepos_move
!fifth version: change logical of INPUT and movement
!sixth version: add rotation
!seventh version : use head and tail atoms to join two structures
!eigth version:    add terminator H
!The biggest Drawback is: the difference of lattice(x,y) 
!between two structures should not be too much
! 

module addmove_mod
implicit none
  real(8) :: scale, latt_sub(3,3),latt_ly(3,3)
  real(8) :: vaccum,dist_req,lattz_req
  integer,allocatable :: iontp_sub(:),iontp_ly(:)
  real(8),allocatable :: pos_sub(:,:),pos_ly(:,:)
  character(2),allocatable :: elemname_sub(:),elemname_ly(:)
  character(20) :: title_sub, title_ly, prefix_outfile
  logical :: lmove(3)
  character(1),allocatable :: ismove_sub(:,:),ismove_ly(:,:)
  integer :: ntype_ly,ntype_sub,nions_sub,nions_ly
  real(8),allocatable :: operation_x(:),operation_y(:),operation_z(:)
  real(8) :: step
  integer :: num_operation(3),movemode
  character(1) :: dir_x,dir_y,dir_z
  logical :: isr(2)
  real(8) :: rthe(2)
  real(8) :: rcenter(3)

contains

subroutine addmove(in_sub,in_ly)
implicit none
character(100) :: in_sub,in_ly
integer :: i,n
logical :: alive
!,isr(2)
!real(8) :: rthe(2)
!real(8) :: rcenter(3)
character(100) :: operation_txt

inquire(file=in_sub,exist=alive)
if(.not. alive) stop 'Specified file1 is not exist!'
inquire(file=in_ly,exist=alive)
if(.not. alive) stop 'Specified file2 is not exist!'

inquire(file='addmovecontrol.in',exist=alive)
if(alive) then
  call readinput
else
  write(*,"(A)") "'addmovecontrol.in' not Found! Default Parameters Will be Used:"
  write(*,"(A)") "1.5   !default vaccum"
  write(*,"(A)") "1.5   !default intermolecule distance"
  write(*,"(A)") "F F F !atoms do not move in-plane"
  write(*,"(A)") "F F   !sub-ly do not rotation"
  vaccum=1.5    !default vaccum
  dist_req=1.5 !default interlayer distance
  lmove=.false.
  isr=.false.
  rthe=0.0
  rcenter=0.0
endif
  call readpos(in_sub,in_ly)
! call normalcell
!version_6:
  if(isr(1) .or. isr(2)) then
    call rotation(isr,rthe,rcenter)
  endif
!----------------
  call mergcell
  call chgdist
  call centrapos
  call writepos
  if(lmove(1) .or. lmove(2) .or. lmove(3)) then
    stop 'Movements is Ignored At Present'
    call movepos
  endif

deallocate(iontp_sub,iontp_ly,pos_sub,pos_ly)
deallocate(elemname_sub,elemname_ly,ismove_sub,ismove_ly)
if(allocated(operation_x)) deallocate(operation_x)
if(allocated(operation_y)) deallocate(operation_y)
if(allocated(operation_z)) deallocate(operation_z)

return
end subroutine

subroutine rotation(isr,rthe,rcenter)
implicit none
logical :: isr(2)
real(8) :: rthe(2)
real(8) :: rcenter(3)
real(8),allocatable :: pos_new(:,:)
real(8),parameter :: pi=3.14159265358979

!original cart (x1,y1) rotation with the center(x2,y2)
!x_new= (x1 - x2)*cos(pi/rthe) - (y1 - y2)*sin(pi/rhte) + x2 ;
!y_new= (x1 - x2)*sin(pi/rthe) + (y1 - y2)*cos(pi/rhte) + y2 ;
!write(*,*) isr(1:2),rthe(1:2),rcenter(1:3)
if(isr(1)) then
if(allocated(pos_new)) deallocate(pos_new)
allocate(pos_new(3,nions_sub))
pos_new(1,:)=(pos_sub(1,:)-rcenter(1))*cos(pi/180*rthe(1))-(pos_sub(2,:)-rcenter(2))*sin(pi/180*rthe(1))+rcenter(1)
pos_new(2,:)=(pos_sub(1,:)-rcenter(1))*sin(pi/180*rthe(1))+(pos_sub(2,:)-rcenter(2))*cos(pi/180*rthe(1))+rcenter(2)
pos_sub(1,:)=pos_new(1,:)
pos_sub(2,:)=pos_new(2,:)
endif

if(isr(2)) then
if(allocated(pos_new)) deallocate(pos_new)
allocate(pos_new(3,nions_ly))
pos_new(1,:)=(pos_ly(1,:)-rcenter(1))*cos(pi/180*rthe(2))-(pos_ly(2,:)-rcenter(2))*sin(pi/180*rthe(2))+rcenter(1)
pos_new(2,:)=(pos_ly(1,:)-rcenter(1))*sin(pi/180*rthe(2))+(pos_ly(2,:)-rcenter(2))*cos(pi/180*rthe(2))+rcenter(2)
pos_ly(1,:)=pos_new(1,:)
pos_ly(2,:)=pos_new(2,:)
endif

deallocate(pos_new)
 
return
end subroutine

subroutine readinput
implicit none
integer :: i,stat
character(1) :: lmove_txt(3),isr_txt(2)
!logical :: isr(2)
!real(8) :: rthe(2),rcenter(3)

dir_x='x'
dir_y='y'
dir_z='n'
lmove_txt='F'
isr_txt='F'
lmove=.FALSE.
isr=.false.
num_operation=0
step=0
vaccum=1.5   !default vaccum
dist_req=1.5 !default intermolecule distance
rthe=0.0
rcenter=0.0

open(33,file='addmovecontrol.in',status='old')
read(33,*,iostat=stat) vaccum       ! specify vaccum length (Angstrom) 
if(stat/=0) stop 'vaccum Reading Error!'
read(33,*,iostat=stat) dist_req     ! specify required intermolecule distance
if(stat/=0) stop 'dist_req Reading Error!'
read(33,*,iostat=stat) lmove_txt(1:3)        ! if move or not
if(stat/=0) lmove=.False.
!stop 'lmove Reading Error!'
do i=1,3
  if(lmove_txt(i)=='t' .or. lmove_txt(i)=='T') lmove(i)=.TRUE.
enddo
!if(lmove(1) .or. lmove(2) .or. lmove(3)) then
!  read(33,*,iostat=stat) prefix_outfile ! prefix of output file
!  if(stat/=0) stop 'prefix_outfile Reading Error!'
!endif
if(lmove(3)) then
  if(lmove(1) .or. lmove(2)) then
    write(*,"(A)") "Movement z is not Compatible with x or y!"
    write(*,"(A)") "Only Four Movement Models is Supported:"
    write(*,"(A)") "T F F //x"
    write(*,"(A)") "F T F //y"
    write(*,"(A)") "T T F //xy"
    write(*,"(A)") "F F T //z"
    write(*,"(A)") "Please Check Your 'addmovecontrol.in' File!"
    write(*,*)
  endif
endif
!x-axis
if(lmove(1)) then
  read(33,*,iostat=stat) dir_x,num_operation(1),movemode
  if(stat/=0) stop 'movements along x Reading Error!'
  if(dir_x=='x' .or. dir_x=='X') then
    if(allocated(operation_x)) deallocate(operation_x)
    allocate(operation_x(num_operation(1)))
    if(movemode==1) then
      read(33,*,iostat=stat) operation_x(1:num_operation(1))
      if(stat/=0) stop 'operation_x Reading Error!'
    elseif(movemode==2) then
      read(33,*,iostat=stat) step         !move step
      if(stat/=0) stop 'step_x Reading Error!'
      do i=1,num_operation(1)
        operation_x(i)=i*step
      enddo
    else
      stop 'The select movement mode is wrong'
    endif
  else
    write(*,"(A)") " Movements on x-direction is set to True, &
                   But I cannot found 'x or X' mark!"
    stop
  endif 
endif 

!y-axis
if(lmove(2)) then
  read(33,*,iostat=stat) dir_y,num_operation(2),movemode
  if(stat/=0) stop 'movements along y Reading Error!'
  if(dir_y=='y' .or. dir_y=='Y') then
    if(allocated(operation_y)) deallocate(operation_y)
    allocate(operation_y(num_operation(2)))
    if(movemode==1) then
      read(33,*,iostat=stat) operation_y(1:num_operation(2))
      if(stat/=0) stop 'operation_y Reading Error!'
    elseif(movemode==2) then
      read(33,*,iostat=stat) step         !move step
      if(stat/=0) stop 'step_y Reading Error!'
      do i=1,num_operation(2)
        operation_y(i)=i*step
      enddo
    else
      stop 'The select movement mode_y is wrong'
    endif
  else
    write(*,"(A)") " Movements on y-direction is set to True, &
                   But I cannot found 'y or Y' mark!"
    stop
  endif
endif

if(lmove(1) .and. lmove(2)) then
  if(num_operation(1)/=num_operation(2)) then
    write(*,"(A)") "Movements Along Both x and y-direction is Found!"
    write(*,"(A)") "However, Number of movements is not Equal with Each Other!"
    stop
  endif
endif

!z-axis
if(lmove(3)) then
  read(33,*,iostat=stat) dir_z,num_operation(3),movemode
  if(stat/=0) stop 'movements along z Reading Error!'
  if(dir_z=='z' .or. dir_z=='Z') then
    if(allocated(operation_z)) deallocate(operation_z)
    allocate(operation_z(num_operation(3)))
    if(movemode==1) then
      read(33,*,iostat=stat) operation_z(1:num_operation(3))
      if(stat/=0) stop 'operation_z Reading Error!'
    elseif(movemode==2) then
      read(33,*,iostat=stat) step         !move step
      if(stat/=0) stop 'step_z Reading Error!'
      do i=1,num_operation(3)
        operation_z(i)=i*step
      enddo
    else
      stop 'The select movement mode_z is wrong'
    endif
  else
    write(*,"(A)") " Movements on z-direction is set to True, &
                   But I cannot found 'z or Z' mark!"
    stop
  endif
endif 

read(33,*,iostat=stat) isr_txt(1:2)        ! if rotation or not
if(stat/=0) isr=.False.
!stop 'rotation Reading Error!'
do i=1,2
  if(isr_txt(i)=='t' .or. isr_txt(i)=='T') isr(i)=.TRUE.
enddo
if(isr(1) .or. isr(2)) then
  read(33,*,iostat=stat) rthe(1:2)
  if(stat/=0) stop 'rotation angel Reading Error!'
  read(33,*,iostat=stat) rcenter(1:3)
  if(stat/=0) stop 'rotation center Reading Error!'
endif

close(33)

return
end subroutine

subroutine readpos(in_sub,in_layer)
implicit none
  integer :: i,tmp_iontp(60)
  character(100) :: in_sub,in_layer,nomean
  integer :: nblock
  character(10) :: iblock(100)

!--read substrate
  tmp_iontp=0
  ntype_sub=0
  open(33,file=in_sub,status='old')
  read(33,'(A)') title_sub
  read(33,*) scale
  do i=1,3
    read(33,*) latt_sub(1:3,i) 
  enddo
  read(33,"(A100)") nomean

call splitblock(nomean,nblock,iblock)

 ntype_sub=nblock
! read(33,*) tmp_iontp(1:60)
! do i=1,60
!   if(tmp_iontp(i)>0) ntype_sub=ntype_sub+1
! enddo
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
  if(nomean(1:1)=='s' .or. nomean(1:1)=='S') then
    read(33,*)
    do i=1,nions_sub
      read(33,*) pos_sub(1:3,i),ismove_sub(1:3,i)
    enddo
  else
    do i=1,nions_sub
      read(33,*) pos_sub(1:3,i)
    enddo
    ismove_sub='F'
  endif
  close(33)
!-----------------------
!read adsorber layer
  tmp_iontp=0
  ntype_ly=0
  open(33,file=in_layer,status='old')
  read(33,'(A)') title_ly
  read(33,*) scale
  do i=1,3
    read(33,*) latt_ly(1:3,i)
  enddo
  read(33,"(A100)") nomean
 call splitblock(nomean,nblock,iblock)

 ntype_ly=nblock
! read(33,*) tmp_iontp(1:60)
! do i=1,60
!   if(tmp_iontp(i)>0) ntype_ly=ntype_ly+1
! enddo
  if(allocated(iontp_ly)) deallocate(iontp_ly)
  if(allocated(elemname_ly)) deallocate(elemname_ly)
  allocate(iontp_ly(ntype_ly),elemname_ly(ntype_ly))
! Backspace(33)
  Backspace(33)
  read(33,*) elemname_ly(1:ntype_ly)
  read(33,*) iontp_ly(1:ntype_ly)
  read(33,*) nomean
! write(*,*) "From Substrate POSCAR, I found :"
! write(*,*) "   ",ntype_ly," type ions: ",elemname_ly(1:ntype_ly)
! write(*,*) "have ",iontp_ly(1:ntype_ly)," repectively"
! write(*,*)
  nions_ly=sum(iontp_ly)
  if(allocated(pos_ly)) deallocate(pos_ly)
  allocate(pos_ly(3,nions_ly))
  if(allocated(ismove_ly)) deallocate(ismove_ly)
  allocate(ismove_ly(3,nions_ly))
  if(nomean(1:1)=='s' .or. nomean(1:1)=='S') then
    read(33,*) 
    do i=1,nions_ly
      read(33,*) pos_ly(1:3,i),ismove_ly(1:3,i)
    enddo
  else
    do i=1,nions_ly
      read(33,*) pos_ly(1:3,i)
    enddo
    ismove_ly='T'
  endif
  close(33)
!-----------------------
!re-exam

if(latt_ly(1,1) /= latt_sub(1,1) .or. latt_ly(2,2) /= latt_sub(2,2)) then
write(*,*) 
write(*,*) "Causion! The lattice (x or y dimension) of adsorber is not mathched &
&           that of thesubstrate!!"
write(*,*) "The lattice (x and y) of substrate will be used!"
write(*,*) 
endif
return
end subroutine

!----------------------
subroutine mergcell
implicit none
integer :: new_lattz

new_lattz=latt_sub(3,3)+latt_ly(3,3)  !Assuming z parperd x and y
pos_sub(3,:)=pos_sub(3,:)*latt_sub(3,3)/new_lattz
pos_ly(3,:)=pos_ly(3,:)*latt_ly(3,3)/new_lattz

return
end subroutine
!------------------------

!---------------------------
subroutine chgdist
implicit none
integer :: i,new_lattz
real(8) :: zmax_sub,zmax_ly,zmin_sub,zmin_ly,new_zmin_ly
real(8) :: delta_dist,dist,delta_z
real(8) :: thick_sub,thick_ly

new_lattz=latt_sub(3,3)+latt_ly(3,3)  !Assuming z parperd x and y

zmax_sub=-1000d0
zmin_sub=1000d0
do i=1,nions_sub
  if(abs(pos_sub(3,i))>zmax_sub) then
    zmax_sub=pos_sub(3,i)
  endif
  if(abs(pos_sub(3,i))<zmin_sub) then
    zmin_sub=pos_sub(3,i)
  endif
enddo

zmax_ly=-1000
zmin_ly=1000
do i=1,nions_ly
  if(abs(pos_ly(3,i))>zmax_ly) then
    zmax_ly=pos_ly(3,i)
  endif
  if(abs(pos_ly(3,i))<zmin_ly) then
    zmin_ly=pos_ly(3,i)
  endif
enddo

thick_sub=(zmax_sub-zmin_sub)*new_lattz
thick_ly=(zmax_ly-zmin_ly)*new_lattz
!dist=(zmin_ly-zmax_sub)*new_lattz

!delta_dist=dist_req-dist

new_zmin_ly=dist_req/new_lattz+zmax_sub

delta_z=new_zmin_ly-zmin_ly
pos_ly(3,:)=pos_ly(3,:)+delta_z

!change lattice with disired vaccum
lattz_req=vaccum+thick_sub+thick_ly+dist_req

pos_sub(3,:)=pos_sub(3,:)*new_lattz/lattz_req
pos_ly(3,:)=pos_ly(3,:)*new_lattz/lattz_req

return
end subroutine
!-------------------------

!------------------
subroutine centrapos
implicit none
integer :: i
real(8) :: zmax_sub,zmin_sub,new_zmin_sub,delta_z

!zmax_sub=-1000d0
zmin_sub=1000d0
do i=1,nions_sub
! if(pos_sub(3,i)>zmax_sub) then
!   zmax_sub=pos_sub(3,i)
! endif
  if(pos_sub(3,i)<zmin_sub) then
    zmin_sub=pos_sub(3,i)
  endif
enddo

new_zmin_sub=vaccum/2/lattz_req
delta_z=new_zmin_sub-zmin_sub
pos_sub(3,:)=pos_sub(3,:)+delta_z
pos_ly(3,:)=pos_ly(3,:)+delta_z

return
end subroutine
!-----------------

subroutine writepos(operation_txt)
integer :: i,j,k,n
character(200) :: ntype_txt,forma
real(8) :: torl
character(2),allocatable :: elemname_all(:),elemname_tmp(:)
character(2),allocatable :: spsymb(:)
character(1),allocatable :: moveornot(:,:,:)
integer :: nspec,natom
real(8),allocatable :: pos(:,:),wpos(:,:,:)
integer,allocatable :: nwpos(:)
character(100),intent(IN),optional :: operation_txt

torl=1.0d-5
! merge elem and position
natom=nions_sub+nions_ly
if(allocated(elemname_all)) deallocate(elemname_all)
allocate(elemname_all(natom))
if(allocated(elemname_tmp)) deallocate(elemname_tmp)
allocate(elemname_tmp(natom))
do i=1,ntype_sub
  n=iontp_sub(i)
  if(i==1) then
    k=0
  else
    k=iontp_sub(i-1)+k
  endif
    elemname_all(k+1:k+n)=elemname_sub(i)
enddo
!write(*,*) "elem_sub: ",elemname_all(1:nions_sub)
do i=1,ntype_ly
  n=iontp_ly(i)
  if(i==1) then
    k=0
  else
    k=iontp_ly(i-1)+k
  endif
  elemname_all(k+1+nions_sub:k+n+nions_sub)=elemname_ly(i)
enddo
!write(*,*) "elem_ly: ",elemname_all(nions_sub+1:natom)
elemname_tmp=elemname_all
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
if(allocated(spsymb)) deallocate(spsymb)
allocate(spsymb(ntype_sub+ntype_ly))
nspec=0
!spsymb='no'
do i=1,natom
  if(elemname_tmp(i)/='no') then
    nspec=nspec+1
    spsymb(nspec)=elemname_tmp(i)
  endif
enddo
deallocate(elemname_tmp)

if(allocated(nwpos)) deallocate(nwpos)
if(allocated(wpos)) deallocate(wpos)
if(allocated(pos)) deallocate(pos)
allocate(nwpos(nspec),wpos(3,natom,nspec))
allocate(pos(3,natom))
do i=1,nions_sub
  pos(:,i)=pos_sub(:,i)
enddo
do i=1,nions_ly
  pos(:,i+nions_sub)=pos_ly(:,i)
enddo

do i=1,nspec
  nwpos(i)=0
  do j=1,natom
    if(spsymb(i)==elemname_all(j)) then
      nwpos(i)=nwpos(i)+1
      wpos(:,nwpos(i),i)=pos(:,j)
    endif
  enddo
enddo

if(present(operation_txt)) then
  open(33,file=trim(adjustl(operation_txt)),status='replace')
else
  open(33,file='POSCAR_sub_ly',status='replace')
endif
!write out poscar
write(33,"(A)") 'MergedCell'
write(33,*) scale
do j=1,2
  write(33,"(3F18.11)") latt_sub(1:3,j)
enddo
write(33,"(3F18.11)") latt_sub(1:2,3),lattz_req
write(ntype_txt,*) nspec
forma='('//trim(adjustl(ntype_txt))//'A4)'
write(33,forma) spsymb(1:nspec)
forma='('//trim(adjustl(ntype_txt))//'I6)'
write(33,forma) nwpos(1:nspec)
write(33,"('Selective')")
write(33,"('Direct')")

if(allocated(moveornot)) deallocate(moveornot)
allocate(moveornot(3,natom,nspec))
moveornot='T'
  do i=1,nspec
    do j=1,nwpos(i)
      do k=1,3
        if(wpos(k,j,i)>0.999999999999d0) then
          wpos(k,j,i)=abs(wpos(k,j,i)-1.0d0)
        else if(wpos(k,j,i)<-0.999999999999d0) then
          wpos(k,j,i)=wpos(k,j,i)+1.0d0
        endif
      enddo
      if(abs(wpos(1,j,i))<=torl) then
        moveornot(1,j,i)='F'
      endif
      if(abs(wpos(2,j,i))<=torl) then
        moveornot(2,j,i)='F'
      endif
      if(abs(wpos(3,j,i))<=torl) then
        moveornot(3,j,i)='F'
      endif
      write(33,"(3F18.11,X,3(X,A))") wpos(1:3,j,i),moveornot(1:3,j,i)
    enddo
enddo
close(33)

deallocate(moveornot,wpos,nwpos,pos,elemname_all,spsymb)
return
end subroutine

!subroutine writepos_move(operation_txt)
!use commvar
!implicit none
!integer :: i,j,k
!character(100) :: ntype_txt,forma
!character(100),intent(IN) :: operation_txt
!real(8) :: torl,torl2
!
!torl=1.0d-4
!torl2=1.0d-5
!
!  open(33,file='POSCAR_'//trim(adjustl(operation_txt)),status='replace')
!write out poscar
!    write(33,"(A)") trim(adjustl(title_sub))//'_'//trim(adjustl(title_ly))
!    write(33,*) scale
!    do j=1,2
!     write(33,"(3F18.11)") latt_sub(1:3,j)
!   enddo
!   write(33,"(3F18.11)") latt_sub(1:2,3),lattz_req
!   write(ntype_txt,*) ntype_sub+ntype_ly
!   forma='('//trim(adjustl(ntype_txt))//'A4)'
!   write(33,forma) elemname_sub(:),elemname_ly(:)
!   forma='('//trim(adjustl(ntype_txt))//'I6)'
!   write(33,forma) iontp_sub(:),iontp_ly(:)
!   write(33,"('Selective')")
!   write(33,"('Direct')")
!   do j=1,nions_sub
!     do i=1,3
!       if(abs(pos_sub(i,j))<=torl) then
!         ismove_sub(i,j)='F'
!       else
!         ismove_sub(i,j)='T'
!       endif
!     enddo
!     write(33,"(3F18.11,2X,3(A,X))") pos_sub(1:3,j),ismove_sub(1:3,j)
!   enddo
!   do j=1,nions_ly
!   do i=1,3
!     if(pos_ly(i,j)>0.999999999999d0) then
!       pos_ly(i,j)=abs(pos_ly(i,j)-1.0d0)
!     else if(pos_ly(i,j)<-0.999999999999d0) then
!       pos_ly(i,j)=pos_ly(i,j)+1.0d0
!     endif
!     
!     if(abs(pos_ly(i,j))<=torl) then
!       ismove_ly(i,j)='F'
!     else
!       ismove_ly(i,j)='T'
!     endif
!   enddo
!     write(33,"(3F18.11,2X,3(A,X))") pos_ly(1:3,j),ismove_ly(1:3,j)
!   enddo
!   close(33)
!
!return
!end subroutine

!------------------
subroutine movepos
implicit none
integer :: i,j
character(100) :: itxt,optxt
character(100) :: operation_txt
real(8),allocatable :: pos_ly_old(:,:)

  if(allocated(pos_ly_old)) deallocate(pos_ly_old)
  allocate(pos_ly_old(3,nions_ly))
  pos_ly_old=pos_ly

! T F F
if(lmove(1) .and. (.not. lmove(2))) then
  optxt='POS_x_'
  do i=1,num_operation(1)
    write(itxt,*) i
    operation_txt=trim(adjustl(optxt))//trim(adjustl(itxt))
    pos_ly(1,:)=pos_ly_old(1,:)+operation_x(i)
    call writepos(operation_txt)
  enddo
endif

! F T F
if(lmove(2) .and. (.not. lmove(1))) then
  optxt='POS_y_'
  do i=1,num_operation(2)
     write(itxt,*) i
     operation_txt=trim(adjustl(optxt))//trim(adjustl(itxt))
     pos_ly(2,:)=pos_ly_old(2,:)+operation_y(i)
     call writepos(operation_txt)
  enddo
endif

! T T F
if(lmove(1) .and. lmove(2)) then
  optxt='POS_xy_'
  do i=1,num_operation(1)
    write(itxt,*) i
    operation_txt=trim(adjustl(optxt))//trim(adjustl(itxt))
    pos_ly(1,:)=pos_ly_old(1,:)+operation_x(i)
    pos_ly(2,:)=pos_ly_old(2,:)+operation_y(i)
    call writepos(operation_txt)
  enddo
endif

! F F T
if(lmove(3)) then
  if((.not. lmove(1)) .and. (.not. lmove(2))) then
    operation_txt='POS_z_'
    do i=1,num_operation(3)
      write(itxt,*) i
      pos_ly(3,:)=pos_ly_old(3,:)+operation_z(i)
      call writepos(operation_txt)
    enddo
  endif
endif

deallocate(pos_ly_old)
return
end subroutine

end module
