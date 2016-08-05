!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                                   !
! Subroutine initial: using pseudoproteins compute the bijection    !
! that maximizes the score using dynamic programming and return     !
! protein A ready for alignment                                     !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial(pseudoa,pseudob,na,nb,prota,protb,seqfix)

  use sizes
  implicit none
  integer :: na, nb, nbij, ngaps, bije(maxatom,2), i
  double precision :: pseudoa(maxatom,3), pseudob(maxatom,3), gap,&
                      dzero2, bijscore(maxatom), score,&
                      prota(maxatom,3), protb(maxatom,3)
  logical :: seqfix

  if ( seqfix ) return

  if(min(na,nb).le.5) then
    write(*,*) ' Too few atoms. Ignoring pseudoprot initial point.'
    return
  end if

  ! Parameters for scoring internal distances

      dzero2 = 100.
      gap = 1.

  ! Initialization based on internal distances and dynamic programming

  call structal(pseudoa,pseudob,na-3,nb-3,dzero2,gap,bije,nbij,&
                bijscore,ngaps,score,seqfix)

  do i = 1, nbij
    bije(i,1) = bije(i,1) + 1
    bije(i,2) = bije(i,2) + 1
  end do

  call procrustes(nbij,na,bije,prota,protb)

end subroutine initial
