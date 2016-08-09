!
! Subroutine triang: Computes the structal score if the
!                    bijection needs to be determined 
!                    using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dtri2: square of the triangular distance 
!           gap: penalization for gaps
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           bijscore: individual scores for each pair of the bijection
!           ngaps: number of gaps of the bijection
!           score: the triangular score (maximum for current orientation)
!

subroutine triang(prota,protb,na,nb,dtri2,gap,bije,nbij,&
                  bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2), npos,&
            nbij_temp, bije_temp(maxatom,2)
  double precision :: prota(maxatom,3), protb(maxatom,3), dtri2,&
                      dist, score, scorin(maxatom,maxatom), gap,&
                      bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    score = 0.d0
    nbij = 0
    do i = 1, na
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      if ( dist < dtri2 ) then 
        nbij = nbij + 1
        bije(nbij,1) = i 
        bije(nbij,2) = i 
        bijscore(nbij) = 1.d0 - dist / dtri2
        score = score + bijscore(nbij)
      end if
    end do
    ngaps = 0
    if ( nbij == 0 ) then
      nbij = na
      do i = 1, na
        bije(i,1) = i
        bije(i,2) = i
        bijscore(i) = 0.d0
      end do
    end if
    return
  end if

  ! Computes individual scores for all pairs

  npos = 0
  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = dmax1(0.d0, 1.d0 - dist / dtri2)
      if(scorin(i,j).gt.0.d0) npos = npos + 1
    end do
  end do

  ! Test if the number of atoms in the bijection is greater than zero

  if( npos == 0 ) then
    nbij_temp = 1
    bije_temp(1,1) = 1
    bije_temp(1,2) = 1
    bijscore(1) = 1.d0
    score = 1.d0
    ngaps = 0
    return
  end if

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije_temp,nbij_temp,&
              bijscore,score)

  ! Include ontly atoms for which the score is positive in the bijection

  nbij = 0
  do i = 1, nbij_temp
    if(bijscore(i).gt.0.d0) then
      nbij = nbij + 1
      bije(nbij,1) = bije_temp(i,1)
      bije(nbij,2) = bije_temp(i,2)
    end if
  end do

end subroutine triang
