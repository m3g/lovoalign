!
! Module to pass some warnings from some routines to the printing routine
!
module warnings

  integer, parameter :: maxwarn = 100
  integer :: nwarn
  character(len=200) :: warn(maxwarn)

end module warnings
