!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
! Subroutine prodin: Performs Dynamic Programming in a matrix of    !
!                    individual scores and returns the bijection    !
!                    that globally maximizes the score and the      !
!                    corresponding score                            !
!                                                                   !
! On input:                                                         !
!          na: number of rows of the individual score matrix        !
!          nb: number of columns of the individual score matrix     !
!          scorin: the individual score matrix                      !
!          gap: Penalization for gaps                               !
!                                                                   !
! On return:                                                        !
!          bije: the bijection that maximizes the score             !
!          nbij: number of correspondences in the bijection         !
!          ngaps: number of gaps of the bijection                   !
!                                                                   !
! This subroutine takes most of the time of the calculations.       !
! If you have a better way to perform this dynamic programming      !
! it would have an important impact on the total alignment time.    !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,scomax)

  use sizes
  implicit none
  integer :: i, j, na, nb, bije(maxatom,2), nbij, ngaps, nl, nc, maxj,&
             maxi(maxatom), seg(maxatom,maxatom,2), imax, jmax 
  double precision scorin(maxatom,maxatom), aux(maxatom,maxatom),&
                   gap, gaf(maxatom,maxatom), bijscore(maxatom),&
                   scomax, gamax

  ! As currently implemented, the cost of this procedure depends on
  ! which protein is the largest one. Therefore, we always set the
  ! the largest protein to be the first

  if(na.ge.nb) then
    do j = 1, nb
      do i = 1, na
        aux(i,j) = scorin(i,j)
      end do
    end do
    nc = na
    nl = nb
  else 
    do j = 1, nb
      do i = 1, na
        aux(j,i) = scorin(i,j)
      end do
    end do
    nc = nb
    nl = na
  end if

  ! Here begin the dynamic programming procedure

  do i = 1, nc
    gaf(i,nl) = aux(i,nl)
    maxi(i) = nl
  end do

  ! Dynamic programming in O(n^2)

  do j = nl - 1, 1, -1
    gaf(nc,j) = aux(nc,j)
    maxj = nc
    do i = nc - 1, 1, -1
      gamax = gaf(i+1,j+1)
      seg(i,j,1) = i + 1
      seg(i,j,2) = j + 1
      if(gaf(i+1,maxi(i+1))-gap.gt.gamax) then
        gamax = gaf(i+1, maxi(i+1)) - gap
        seg(i,j,1) = i+1
        seg(i,j,2) = maxi(i+1)
      end if
      if(gaf(maxj, j+1)-gap.gt.gamax) then
        gamax = gaf(maxj,j+1) - gap
        seg(i,j,1) = maxj
        seg(i,j,2) = j + 1
      end if
      if(gaf(i,j+1).ge.gaf(maxj,j+1)) maxj = i
      gaf(i,j) = aux(i,j) + gamax
    end do
    do i = 1, nc
      if(gaf(i,j).ge.gaf(i, maxi(i))) maxi(i)=j
    end do
  end do
      
  ! Now will compute the score to the end

  imax = 0 ! avoid warning
  jmax = 0 ! avoid warning 
  scomax = 0.d0
  do j = 1, nl
    if(gaf(1,j).ge.scomax) then
      imax = 1
      jmax = j
      scomax = gaf(1,j)
    end if
  end do
  do i = 1, nc
    if(gaf(i,1).ge.scomax) then
      imax = i
      jmax = 1
      scomax = gaf(i,1)
    endif
  end do

  ! Save maximum score obtained and the corredponding bijection      

  bije(1,1) = imax
  bije(1,2) = jmax
  nbij = 1
  i = 2
  do i = 2, nc
    if(bije(i-1,1).eq.nc.or.bije(i-1,2).eq.nl) exit
    bije(i,1) = seg(bije(i-1,1),bije(i-1,2),1)
    bije(i,2) = seg(bije(i-1,1),bije(i-1,2),2)
    nbij = i
  end do

  ngaps = 0
  do i = 1, nbij-1
    if(bije(i+1,1).ne.bije(i,1)+1.or.&
       bije(i+1,2).ne.bije(i,2)+1) then
      ngaps = ngaps+1
    endif
  end do

  do i = 1, nbij
    bijscore(i) = aux(bije(i,1),bije(i,2))
  end do

  ! If the largest protein was protein B, restore

  if(na.lt.nb) then
    do i = 1, nbij
      j = bije(i,1)
      bije(i,1) = bije(i,2)
      bije(i,2) = j
    end do
  end if

end subroutine prodin
