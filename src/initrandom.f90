module initrandom

  implicit none
  integer :: seed

  contains

    !
    ! Subroutine that returns a random number given a the seed
    !
    
    subroutine init_random_number(seed)
      implicit none
      integer :: i, seed, iseed(12)
      do i = 1, 12
        iseed(i) = i*seed
      end do
      call random_seed(put=iseed)
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

end module initrandom
