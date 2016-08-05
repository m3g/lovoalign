!
! Subroutine moveprot: Moves a protein relative to its initial
!                      orientation given the displacement and
!                      the euler angles of rotation
!
! On input: 
!          x(1)...x(3): Displacement of the center of mass of
!                       the protein to be applied
!          x(4)...x(6): Euler angles of the rotation to be
!                       be applied
!          na: Number of atoms of the protein
!          prota: Original coordinates of the protein
!
! On return:
!          prota: The coordinates rotated and translated
!

subroutine moveprot(x,na,prota)

  use sizes
  implicit none
  integer :: i, na 
  double precision :: x(6), prota(maxatom,3), ca, sa, cb, sb, cg, sg,&
                      xt, yt, zt

  ca = dcos(x(4))
  sa = dsin(x(4))
  cb = dcos(x(5))
  sb = dsin(x(5))
  cg = dcos(x(6))
  sg = dsin(x(6))

  do i = 1, na
    xt = prota(i,1) 
    yt = prota(i,2) 
    zt = prota(i,3) 
    prota(i,1) = x(1) + cb*ca*xt - cb*sa*yt - sb*zt
    prota(i,2) = x(2) + ( cg*sa - sg*sb*ca )*xt &
                      + ( cg*ca + sg*sb*sa )*yt &
                      - sg*cb*zt
    prota(i,3) = x(3) + ( sg*sa + cg*sb*ca )*xt &
                      + ( sg*ca - cg*sb*sa )*yt &
                      + cg*cb*zt
  end do

end subroutine moveprot
