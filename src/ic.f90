!
! Function that gets the first character of the file
! name without the path
!

integer function ic(string)

  implicit none
  integer :: length
  character(len=200) :: string

  ic = length(string)   
  do while(string(ic:ic).ne.'/'.and.string(ic:ic).ne.'\\')
    ic = ic - 1
    if ( ic == 0 ) exit
  end do
  ic = ic + 1

end function ic
