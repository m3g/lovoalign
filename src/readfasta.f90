!
! Subroutine that reads the sequence alignment in fasta
! format
!
subroutine readfasta()

  use sizes
  use inputpars
  use bijetype
  implicit none
  integer :: i, j, k, ioerr
  character(len=200) :: record
  integer :: length
  logical :: empty_char

  open(10,file=fastafile,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open fasta file: ', fastafile(1:length(fastafile))
    stop
  end if

  k = 0
  i = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( record(1:1) == ">" ) then
      i = i + 1
      k = 0
      cycle
    end if
    do j = 1, length(record)
      if ( .not. empty_char(record(j:j)) ) then
        k = k + 1
        fasta(i)%seq(k:k) = record(j:j)
      end if
    end do
  end do

  !write(*,*) ' FASTA alignment: '
  !write(*,*) trim(adjustl(fasta(1)%seq))
  !write(*,*) trim(adjustl(fasta(2)%seq))
  !stop

end subroutine readfasta

!
! Function that determines if a character is empty (empty, space, or tab)
! (nice suggestion from Ian Harvey -IanH0073- at github)
!

function empty_char(ch)
  character :: ch
  logical empty_char
  empty_char = .false.
  if ( ch == '' .or. &
       ch == achar(9) .or. &
       ch == achar(32) ) then
    empty_char = .true.
  end if
end function empty_char

