program GenLS
implicit none
integer :: n,mode,iargc
character(100) :: smode

! Direct Way
n=iargc()
if(n>=1) then
  call getarg(1,smode)
endif
smode=trim(adjustl(smode))
if(smode(1:2) == '-1' ) then
  write(*,*)
  write(*,'(A)') 'Entering Vertical/Lateral Stacking Strucutures Constructing...'
  call genstack

else if (smode(1:2) == '-2' ) then
  write(*,*)
  write(*,'(A)') 'Entering Point Defects Strucutures Constructing...'
  call gendefect

else if (smode(1:2) == '-h' ) then
  call genlshelp

else if (smode(1:3) == '-g2' ) then
  write(*,*)  
  write(*,'(A)') "Preparing Point Defects Mode File (Example):"
  write(*,'(A)') "'inpara.in'"
  open(20,file='inpara.in',status='replace')
  write(20,'(A)') "8       ! Number of atoms that may be substituted/vacant"
  write(20,'(A)') "2       ! Number of defects"
  write(20,'(A)') "POSCAR.subst    ! POSCAR name"
  close(20)

else if (smode(1:3) == '-g1' ) then
  write(*,*) 
  write(*,'(A)') "Preparing Stacking Mode Files (Example):"
  write(*,'(A)') "'inpara.in' and 'addmovecontrol.in'"
  open(20,file='inpara.in',status='replace')
  write(20,'(A)') "F       ! periodic condition? T/F"
  write(20,'(A)') "2       ! how many types of build blocks?"
  write(20,'(A)') "2 2     ! number of build blocks in each type!"
  close(20)
  open(20,file='addmovecontrol.in',status='replace')
  write(20,'(A)') "15      ! specify vaccum length (Angstrom)" 
  write(20,'(A)') "1.45    ! specify required interlayer distance"
  close(20)

! Indirect Way 'using menu'
else
n=0
do while (.true.)
  n=n+1
  call genlsmenu(smode)
  read(smode,*,ERR=900) mode
  if(mode .eq. 1) then
    write(*,*)
    write(*,'(A)') 'Entering Vertical/Lateral Stacking Strucutures Constructing...'
    call genstack
    exit
  else if(mode .eq. 2) then
    write(*,*)
    write(*,'(A)') 'Entering Point Defects Strucutures Constructing...'
    call gendefect
    exit
  else if(mode .eq. 3) then
    write(*,*)
    write(*,'(A)') "Preparing Stacking Mode Files (Example):"
    write(*,'(A)') "'inpara.in' and 'addmovecontrol.in'"
    open(20,file='inpara.in',status='replace')
    write(20,'(A)') "F       ! periodic condition? T/F"
    write(20,'(A)') "2       ! how many types of build blocks?"
    write(20,'(A)') "2 2     ! number of build blocks in each type!"
    close(20)
    open(20,file='addmovecontrol.in',status='replace')
    write(20,'(A)') "15      ! specify vaccum length (Angstrom)"
    write(20,'(A)') "1.45    ! specify required interlayer distance"
    close(20)
    exit

  else if(mode .eq. 4) then
    write(*,*)
    write(*,'(A)') "Preparing Point Defects Mode File (Example):"
    write(*,'(A)') "'inpara.in'"
    open(20,file='inpara.in',status='replace')
    write(20,'(A)') "8       ! Number of atoms that may be substituted/vacant"
    write(20,'(A)') "2       ! Number of defects"
    write(20,'(A)') "POSCAR.subst    ! POSCAR name"
    close(20)
    exit
  
  else if(mode .eq. 5) then
    write(*,*)
    stop '! Exit GenLS !'
  else
    if(n<3) then
      write(*,*)
      write(*,'(I2,A)') mode,' is not Supported! Please Choose the Correct Mode Again!'
    else
      write(*,*)
      write(*,'(A)') 'Serious? Too Much Errors in Modes Selecting, Exit!'
      stop
    endif
  endif

900 write(*,*) 'Only Integer Number is Supported, Choose Again!'

enddo
endif

stop
end program

subroutine genlsmenu(smode)
implicit none
character(100) :: smode

write(*,*)
write(*,'(A)') "#-------------------------------------------#"
write(*,'(A)') "#---------------    GenLS   ----------------#"
write(*,'(A)') "#--   v1.1                                --#"
write(*,'(A)') "#--  by arlonne                           --#"
write(*,'(A)') "#--  longhuali@ujs.edu.cn                 --#"
write(*,'(A)') "#-------------------------------------------#"
write(*,*)
write(*,'(A)') "Please Select the Generation Mode:"
write(*,'(A)') "   1 : Stacking Mode (For Vertical/Lateral Superlattice)"
write(*,'(A)') "   2 : Defect Mode (For Point Defects Strucutures)"
write(*,'(A)') "   3 : Generate Input Files for Stacking Mode"
write(*,'(A)') "   4 : Generate Input File  for Defect Mode"
write(*,'(A)') "   5 : Do Nothing and Exit"
write(*,'(A5)',advance='no') "-->: "
read(*,'(A20)') smode
return
end subroutine

subroutine genlshelp()
implicit none
write(*,'(A)') "USEAGE: GenLS.ne [OPTIONS]"
write(*,*)
write(*,'(A)') "OPTIONS: "
write(*,'(A)') "   -h  : show this help message and exit"
write(*,*)
write(*,'(A)') "  -g1  : generate an example file for Stacking Mode"
write(*,*)
write(*,'(A)') "  -g2  : generate an example file for Defect Mode"
write(*,*)
write(*,'(A)') "   -1  : Enter Stacking Mode"
write(*,'(A)') "         You should have Correct 'inpara.in' and 'addmovecontrol.in'"
write(*,'(A)') "         As well as building blocks, such as 'POSCAR.1''POSCAR.2'..."
write(*,*)
write(*,'(A)') "   -2  : Enter Defect Mode"
write(*,'(A)') "         You should have Correct 'inpara.in'"
write(*,'(A)') "         As well as protype structure, such as 'POSCAR.vacan'"
write(*,*)
return
end subroutine
