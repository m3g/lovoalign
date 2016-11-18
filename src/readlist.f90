!
! Subroutine readlist: Reads names from the pdblist for running
! modes 1 and 2
!

subroutine readlist(pdblist,pdbfiles,nfiles)

  use sizes
  use ioformat
  implicit none
  integer :: nfiles, length, ioerr, namesize, ic
  character(len=200) :: pdblist, pdbfiles(maxfiles), record

  open(10,file=pdblist,status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', pdblist(1:length(pdblist))
    stop
  end if

  max_filename_size = 0
  nfiles = 0
  do while(.true.)
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(record.gt.' ') then
      nfiles = nfiles + 1
      namesize = length(record)
      pdbfiles(nfiles) = record(1:namesize)
      max_filename_size = max(max_filename_size,namesize-ic(record)+1)
    end if
  end do
  close(10)

  if(nfiles.gt.maxfiles) then
    write(*,*) ' Number of files in list greater than MAXFILES. '
    write(*,*) ' Increase the maxfiles parameter. '
    stop
  end if

  call listformats()

end subroutine readlist
