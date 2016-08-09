!
! Subroutine that performs a random initial alignment
!

subroutine randomini(na,nb,prota,protb,nbij,bije)

  use sizes
  integer :: ia, ib, na, nb, i, nbij, bije(maxatom,2)
  double precision :: prota(maxatom,3), protb(maxatom,3)
  double precision :: random

  call random_number(random)
  ia = int(random*(na-4))+1
  call random_number(random)
  ib = int(random*(nb-4))+1
  nbij = 4
  do i = 1, 4
    bije(i,1) = ia + i - 1
    bije(i,2) = ib + i - 1
  end do
  call procrustes(nbij,na,bije,prota,protb)

end subroutine randomini
