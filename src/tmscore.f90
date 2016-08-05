!
! Subroutine tmscore: Computes the TM-SCORE if the
!                     bijection needs to be determined 
!                     using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dzero2: square of the structal normalization factor
!           gap: penalization for gaps (usually 0 for TM-SCORE)
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           ngaps: number of gaps of the bijection
!           score: the TM-SCORE (maximum for current orientation)
!

subroutine tmscore(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                   bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2)
  double precision prota(maxatom,3), protb(maxatom,3), dzero2,&
                   dist, score, scorin(maxatom,maxatom), gap,&
                   bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    score = 0.d0
    do i = 1, na
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      score = score + 1.d0 / ( 1.d0 + dist / dzero2 )
      bije(i,1) = i
      bije(i,2) = i
    end do
    score = score / dfloat(na)
    nbij = na
    ngaps = 0
    return
  end if

  ! Computes individual scores for all pairs

  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = 1.d0 / ( 1.d0 + dist / dzero2 )
    end do
  end do

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,score)
  score = score / dfloat(nb)

end subroutine tmscore
