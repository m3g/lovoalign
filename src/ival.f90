!
! Function that gets the integer value from a variable string
!

integer function ival(string)

  implicit none
  integer :: ioerr
  character(len=200) :: string

  read(string,*,iostat=ioerr) ival
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read integer value from some'
    write(*,*) '        of the parameters. Some parameter with'
    write(*,*) '        a expected integer value was not set '
    write(*,*) '        using -keyword [integer]'
    stop
  end if

end function ival
