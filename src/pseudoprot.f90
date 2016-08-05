!
! Subroutine pseudoprot: Computes pseudoprotein A for initial point
!                        generation
!

subroutine pseudoprot(prota,pseudoa,na)

  use sizes
  implicit none
  integer :: i, na
  double precision :: prota(maxatom,3), pseudoa(maxatom,3), dist

  if(na.le.5) then
    return
  end if

  do i = 1, na - 4
    pseudoa(i,1) = dist(prota,i,i+2)
    pseudoa(i,2) = dist(prota,i,i+3)
    pseudoa(i,3) = dist(prota,i,i+4)
  end do
  pseudoa(na-3,1) = dist(prota,na-3,na-1) 
  pseudoa(na-3,2) = dist(prota,na-3,na)
  pseudoa(na-3,3) = dist(prota,na-2,na)

end subroutine pseudoprot
