module mod_interface
implicit none

interface

subroutine rmreppermu(aarray,barray,asize,bsize,lp)
implicit none
integer,allocatable,intent(in) :: aarray(:,:)
integer,intent(in) :: asize
integer,intent(in) :: bsize
logical,intent(in) :: lp
integer,allocatable,intent(inout) :: barray(:)
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
end subroutine

end interface

end module
