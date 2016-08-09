!        
! Subroutine writenonbije: Writes the correspondence of a non-bijective 
!                          score
!                         

subroutine writenonbije(bije,nbij)

  use sizes
  implicit none
  integer :: i, bije(maxatom,2), nbij

  write(*,"(/,'  ',25('-'),' CORRESPONDENCE ',26('-'))")
  do i = 1, nbij
    write(*,"( i5,' ->', i5 )") bije(i,1), bije(i,2)
  end do
  write(*,*)

end subroutine writenonbije
