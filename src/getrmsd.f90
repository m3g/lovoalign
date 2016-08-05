!
! Subroutine getrmsd: Computes the rmsd given the bijection
!

subroutine getrmsd(prota,protb,bije,nbij,rmsd)

  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2)
  double precision :: rmsd, dist, prota(maxatom,3), protb(maxatom,3)

  rmsd = 0.d0
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    rmsd = rmsd + dist
  end do
  rmsd = dsqrt(rmsd/dfloat(nbij))

end subroutine getrmsd
