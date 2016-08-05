!
! Subroutine get_nearest: gets the nearest atom of protein B to point 
! X using the fast algorithm that uses the ordered internal distances
! of protein B
!
! On input:
!   protb: The coordinates of proteinB
!   nb: The number of atoms of proteinB
!   indisord: Matrix containing the indexes of of the atoms of 
!             proteinB such that, for each atom, the other 
!             atoms are ordered from closer to farther. 
!   disord: Matrix containing the internal distances of proteinB
!   x, y, z: coordinates of current atom (of proteinA)
!   nearest: First guess of what atom is the closest one
!
! On output:
!   nearest: The index of the closest atom.
!   dmin: The (squared) distance of the nearest atom to x.
!
      
subroutine get_nearest(protb,nb,indisord,disord,&
                       x, y, z,&
                       nearest, dmin)

  use sizes
  implicit none
  double precision :: protb(maxatom,3), disord(maxatom-1,maxatom),&
                      x, y, z, dmin, d, dbase
  integer :: nearest, i, nb, indisord(maxatom-1,maxatom), ibase

  d = ( x - protb(nearest,1) )**2 + &
      ( y - protb(nearest,2) )**2 + &
      ( z - protb(nearest,3) )**2
  dmin = d
  ibase = nearest
  dbase = dmin

  do i = 1, nb - 1

    if( disord(i,ibase) .ge. 4.d0*dbase ) return

    d = ( x - protb(indisord(i,ibase),1) )**2 + &
        ( y - protb(indisord(i,ibase),2) )**2 + &
        ( z - protb(indisord(i,ibase),3) )**2

    if( d .lt. dmin ) then
      dmin = d
      nearest = indisord(i,ibase)
    end if

  end do

end subroutine get_nearest

