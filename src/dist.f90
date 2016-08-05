! Function that computes the distance between two atoms

double precision function dist(prota,i,j)

  use sizes
  implicit none
  integer :: i, j
  double precision :: prota(maxatom,3)

  dist = dsqrt((prota(i,1)-prota(j,1))**2 + &
               (prota(i,2)-prota(j,2))**2 + &
               (prota(i,3)-prota(j,3))**2)
  
end function dist
