!
! Subroutine nonbscore: Subroutine that computes the non-bijective
!                       triangular score.
! 
! On input:
!   na: Number of atoms of protein A
!   nb: Number of atoms of protein B
!   prota: Coordinates of protein A
!   protb: Coordinates of protein B
!   dtri2: Squared threshold distance (above which atoms will
!          not be considered in the correspondence)
!   gap: penalization for gaps
!   disord: Matrix containing internal distances of protein B
!   indisord: Matrix containing the indexes of internal distances
!             of protein B, from the closer to the farther atom,
!             for each atom
!   it: Current iteration of the algorithms
!   pair: The atom o B associated with each atom of A
!
! On return:
!   score: The non-bijective triangular score.
!   bije: The correspondence between atoms obtained.
!   nbij: The number of pairs of the correspondence.
!   pair: The new atom to atom pair association (to be used
!         in the next iteration of the algorithm).
!   ngaps: Number of gaps of the correspondence
!

subroutine nonbscore(na, nb, prota, protb, dtri2, gap,&
                     disord, indisord, it, pair,& 
                     score, bije, nbij, ngaps)

  use sizes
  implicit none
  double precision :: prota(maxatom,3), protb(maxatom,3),&
                      disord(maxatom-1, maxatom),&
                      distpair(maxatom), score, dtri2, dmin, gap
  integer :: na, nb, nbij, bije(maxatom,2), pair(maxatom),&
             indisord(maxatom-1, maxatom),&
             i, nearest, it, ngaps

  ! If this is the first iteration, the first guess is not very good 

  if(it.eq.0) then

    ! First guess for first atom

    nearest = 1
    call get_nearest(protb,nb,indisord,disord,&
                     prota(1,1), prota(1,2), prota(1,3),&
                     nearest,dmin)
    pair(1) = nearest
    distpair(1) = dmin

    ! Guesses for other atoms (use the result for first atom to
    ! improve guess

    do i = 2, na
      nearest = pair(i-1)
      call get_nearest(protb,nb,indisord,disord,&
                       prota(i,1), prota(i,2), prota(i,3),&
                       nearest,dmin)
      pair(i) = nearest
      distpair(i) = dmin
    end do

    ! For the second iteration on, use the previous assignment as the 
    ! guess for the next one

  else
    do i = 1, na
      nearest = pair(i)
      call get_nearest(protb,nb,indisord,disord,&
                       prota(i,1), prota(i,2), prota(i,3),&
                       nearest,dmin)
      pair(i) = nearest
      distpair(i) = dmin
    end do
  end if

  ! Include in the correspondence only pairs for which the distance
  ! is smaller than dtri, and compute TRIANGULAR SCORE

  nbij = 0
  score = 0.d0
  do i = 1, na 
    if( distpair(i) .lt. dtri2 ) then
      nbij = nbij + 1
      bije(nbij,1) = i
      bije(nbij,2) = pair(i)
      score = score + 20.d0*( 1.d0 - distpair(i) / dtri2 )
    end if
  end do

  ! Penalization for gaps (not reasonable to use for this type of alignment)

  if( gap .gt. 1.d-10 ) then
    ngaps = 0
    do i = 1, nbij - 1
      if ( bije(i+1,1) .ne. bije(i,1) + 1 .or.&
           bije(i+1,2) .ne. bije(i,2) + 1 ) then
        ngaps = ngaps + 1
      end if
    end do
    score = score - gap*ngaps
  end if

end subroutine nonbscore
