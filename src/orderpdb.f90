!
! Subroutine orderpdb: Determines the number of atoms of all proteins in
! a list of pdb files and order the list from the largest to the
! smallest protein. Used for all-on-all database comparison using 
! methods with fast-nonbijective correspondences.
!

subroutine orderpdb(pdbfiles,nfiles)

  use sizes
  implicit none

  real :: natoms(maxfiles)
  integer :: nfiles, i, indflash(maxfiles), lflash(maxfiles), mflash, length, ioerr
  character(len=200) :: record, pdbfiles(maxfiles), aux(maxfiles)

  ! Opening files and reading the number of CA atoms

  do i = 1, nfiles
    record = pdbfiles(i)
    aux(i) = record
    natoms(i) = 0.
    open(10,file=record(1:length(record)),status='old',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    do while(.true.)
      read(10,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit 
      if(natoms(i).ge.1.and.record(1:3).eq.'END') exit
      if(record(1:4).eq.'ATOM'.and.record(14:15).eq.'CA') natoms(i) = natoms(i) + 1.
    end do
    close(10)
  end do

  ! Ordering the file names according to protein size 

  mflash = 1 + int(float(nfiles)/10.)
  call flash1(natoms, nfiles, lflash, mflash, indflash)
  do i = 1, nfiles
    pdbfiles(i) = aux(indflash(nfiles-i+1))
  end do

end subroutine orderpdb
