!
! Function that sets the length of a string
!

integer function length(string)

  implicit none
  character(len=200) :: string

  length = 200
  do while(string(length:length).le.' ')
    length = length - 1
  end do

  return
end function length
