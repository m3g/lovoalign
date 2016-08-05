!
! Function that gets the double precision value from a variable string
!

double precision function dval(string)

  implicit none
  integer :: ioerr
  character(len=200) :: string

  read(string,*,iostat=ioerr) dval
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read real value from some'
    write(*,*) '        of the parameters. Some parameter with'
    write(*,*) '        a expected real value was not set '
    write(*,*) '        using -keyword [real]'
    stop
  end if

end function dval
