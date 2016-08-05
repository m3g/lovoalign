!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                                      !
! Subroutine orprot: Computes the ordered square distance matrices for !
! the fast evaluation of correspondences                               !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine orprot(nprot,proteina,disord,indisord)

  use sizes
  implicit none
  
  ! disprot(nb,nb) contains the internal square distances of protein B
  ! disord(nb,nb-1) contains, for each atom of B, the distances to the other
  ! atoms, in increasing order indisord(nb,nb-1) contains the indices of the
  ! atoms corresponding to disord.
 
  double precision :: proteina(maxatom, 3), disprot(maxatom, maxatom),&
                      disord(maxatom-1, maxatom), z
  real :: aux(maxatom)
  integer :: i, j, k, nprot
  integer :: indisord(maxatom-1, maxatom), lflash(maxatom), mflash, indflash(maxatom) 

  !  Compute disprot, squared internal distances

  mflash = 1 + int(float(nprot)/10.)
  do i = 1, nprot
    disprot(i, i) = 0.d0
  end do

  do j = 1, nprot-1
    do i = j+1, nprot
      z = 0.d0
      do k = 1, 3
        z = z + (proteina(i, k) - proteina(j, k))**2
      end do
      disprot(i, j) = z
      disprot(j, i) = disprot(i, j)
    end do
  end do

  ! Sort the distances to atom i of the protein

  do j = 1, nprot
    do i = 1, nprot
      aux(i) = real(disprot(i,j))
    end do
    call flash1(aux, nprot, lflash, mflash, indflash)
    do i = 1, nprot-1
      disord(i,j) = aux(i+1)
      indisord(i,j) = indflash(i+1)
    end do
  end do

end subroutine orprot
