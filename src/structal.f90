!
! Subroutine structal: Computes the structal score if the
!                      bijection needs to be determined 
!                      using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dzero2: square of the structal normalization factor
!           gap: penalization for gaps
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           bijscore: individual scores for each pair of the bijection
!           ngaps: number of gaps of the bijection
!           score: the structal score (maximum for current orientation)
!

subroutine structal(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                    bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2), nmin
  double precision :: prota(maxatom,3), protb(maxatom,3), dzero2,&
                      dist, score, scorin(maxatom,maxatom), gap,&
                      bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    nmin = min(na,nb)
    score = 0.d0
    do i = 1, nmin
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      bijscore(i) = 20.d0 / ( 1.d0 + dist / dzero2 )
      score = score + bijscore(i)
      bije(i,1) = i
      bije(i,2) = i
    end do
    nbij = nmin
    ngaps = 0
    return
  end if

  ! Computes individual scores for all pairs

  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = 20.d0 / ( 1.d0 + dist / dzero2 )
    end do
  end do

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,score)

end subroutine structal
