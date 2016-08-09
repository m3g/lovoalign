!
! Subroutine that returns a random number given a the seed
!

subroutine init_random_number(iseed)
  integer :: i, seed(12), iseed
  do i = 1, 12
    seed(i) = i*iseed
  end do
  call random_seed(put=seed)
  return
end subroutine init_random_number

!
! Subroutine that uses the date to create a random seed
! 

subroutine seed_from_time(seed)

  implicit none
  integer :: seed, value(8)
  character(len=10) :: b(3)
  call date_and_time( b(1), b(2), b(3), value )
  seed = value(1)+value(2)+value(3)+value(4)+value(5)+value(6)+value(7)+value(8)
  seed = seed + value(1)+value(2)+value(3)+value(4)+value(5)/100+value(6)*100+value(7)/10+value(8)*10

end subroutine seed_from_time

