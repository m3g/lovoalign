!
! subroutine tocm: Put the center of mass of protein A
!                  together with the center of mass of protein B,
!                  used as a simple initial point if the 
!                  pseudoprotein initial point is not wanted
!

subroutine tocm(prota,protb,na,nb)

  use sizes
  implicit none
  integer :: na, nb, i
  double precision :: prota(maxatom,3), protb(maxatom,3), cma(3), cmb(3)

  cma(1) = 0.d0
  cma(2) = 0.d0
  cma(3) = 0.d0
  cmb(1) = 0.d0
  cmb(2) = 0.d0
  cmb(3) = 0.d0

  do i = 1, na
    cma(1) = cma(1) + prota(i,1)
    cma(2) = cma(2) + prota(i,2)
    cma(3) = cma(3) + prota(i,3)
  end do
  cma(1) = cma(1) / dfloat(na)
  cma(2) = cma(2) / dfloat(na)
  cma(3) = cma(3) / dfloat(na)

  do i = 1, nb
    cmb(1) = cmb(1) + protb(i,1)
    cmb(2) = cmb(2) + protb(i,2)
    cmb(3) = cmb(3) + protb(i,3)
  end do
  cmb(1) = cmb(1) / dfloat(nb)
  cmb(2) = cmb(2) / dfloat(nb)
  cmb(3) = cmb(3) / dfloat(nb)

  do i = 1, na 
    prota(i,1) = prota(i,1) + ( cmb(1) - cma(1) )
    prota(i,2) = prota(i,2) + ( cmb(2) - cma(2) )
    prota(i,3) = prota(i,3) + ( cmb(3) - cma(3) )
  end do

end subroutine tocm
