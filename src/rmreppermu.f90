subroutine rmreppermu(aarray,barray,asize,bsize,lp)
implicit none
integer,allocatable,intent(in) :: aarray(:,:)
integer,intent(in) :: asize
integer,intent(in) :: bsize
logical,intent(in) :: lp
integer,allocatable,intent(inout) :: barray(:) !recoding repeat position of aarray
integer,allocatable :: tmparray(:),tmparray_rev(:)
integer,allocatable :: tmparray_fwd(:),tmparray_rev2(:)
integer :: i,j,k,l,m,n,o

allocate(tmparray(asize),tmparray_rev(asize))
if(lp) then
  allocate(tmparray_fwd(asize),tmparray_rev2(asize))
endif
!check periodic condition
do i=1,bsize-1
  if(barray(i)/=0) then
   tmparray(:)=aarray(:,i)
   if(lp) then
     do k=1,asize
       tmparray_rev(k)=tmparray(asize-k+1)
     enddo
     do l=1,asize-1
  ! move forword
       forall(j=1:asize,j>1) tmparray_fwd(j-1)=tmparray(j)
       tmparray_fwd(asize)=tmparray(1)
       tmparray=tmparray_fwd
  !reverse tmparray_fwd
       do k=1,asize
         tmparray_rev2(k)=tmparray_fwd(asize-k+1)
       enddo
  !compare with aarray()
       do k=i+1,bsize
         n=0
         m=0
         o=0
         do j=1,asize
           if(tmparray_fwd(j)==aarray(j,k)) n=n+1
           if(tmparray_rev2(j)==aarray(j,k)) m=m+1
           if(tmparray_rev(j)==aarray(j,k)) o=o+1
         enddo
         if(n==asize .or. m==asize .or. o==asize) then
           barray(k)=0
         endif
       enddo 
     enddo
   else
     do k=1,asize
       tmparray_rev(k)=tmparray(asize-k+1)
     enddo
  !compare with aarray()
     do k=i+1,bsize
       m=0
       do j=1,asize
         if(tmparray_rev(j)==aarray(j,k)) m=m+1
       enddo
       if(m==asize) then
         barray(k)=0
       endif
     enddo
   endif
 endif
enddo

!

deallocate(tmparray,tmparray_rev)
if(allocated(tmparray_fwd)) deallocate(tmparray_fwd)
if(allocated(tmparray_rev2)) deallocate(tmparray_rev2)
return
end subroutine
