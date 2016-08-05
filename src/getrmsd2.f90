!
! Subroutine getrmsd2: Computes the rmsd given the bijection,
!                      only for atoms which are closer than
!                      some tolerance deffined by the dtri parameter
!

subroutine getrmsd2(prota,protb,bije,nbij,rmsd,nbij_dtri,dtri)

  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2), nbij_dtri
  double precision :: rmsd, dist, prota(maxatom,3), protb(maxatom,3),&
                      dtri, dtri2

  dtri2 = dtri*dtri
  rmsd = 0.d0
  nbij_dtri = nbij
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    if(dist.le.dtri2) then
      rmsd = rmsd + dist
    else
      nbij_dtri = nbij_dtri - 1
    end if         
  end do
  if(nbij_dtri.gt.0) then
    rmsd = dsqrt(rmsd/dfloat(nbij_dtri))
  else
    rmsd = 0.d0
  end if

end subroutine getrmsd2
