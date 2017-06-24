!
! Subroutine that performs a random initial alignment
!

subroutine randomini(na,nb,prota,protb,nbij,bije)

  use sizes
  use inputpars
  use bijetype
  integer :: ia, ib, na, nb, i, nbij, bije(maxatom,2)
  double precision :: prota(maxatom,3), protb(maxatom,3)
  double precision :: random

  nbij = max(4,min(na/10,nb/10))
  call random_number(random)
  ia = int(random*(na-nbij))+1

  if ( seqtype > 0 ) then
    do i = 1, nbij
      bije(i,1) = fixbije(ia+i-1,1)
      bije(i,2) = fixbije(ia+i-1,2)
    end do
  else
    call random_number(random)
    if ( random > 0.5d0 ) then
      call random_number(random)
      ib = int(random*(nb-nbij))+1
    else
      ib = ia
    end if
    do i = 1, nbij
      bije(i,1) = ia + i - 1
      bije(i,2) = ib + i - 1
    end do
  end if
  call procrustes(nbij,na,bije,prota,protb)

end subroutine randomini
