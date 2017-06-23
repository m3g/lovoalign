!
! Subroutine that reads the parameters from the command line
! 

subroutine getpars()

  use sizes
  use bijetype
  use inputpars
  use ioformat
  implicit none

  integer :: narg, length, i, ival, iargc, ioerr
  double precision :: dval
  character(len=200) :: keyword, value

  ! Reading the command line specifications

  narg = iargc() 
  i = 1
  gap = 1.d30
  skip = .false.
  do while(i.le.narg)
    call getarg(i,keyword)
    if(keyword(1:length(keyword)).eq.'-p1') then
      call getarg(i+1,value)
      protea = value
    else if(keyword(1:length(keyword)).eq.'-p2') then 
      call getarg(i+1,value)
      proteb = value
    else if(keyword(1:length(keyword)).eq.'-c1') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) chaina
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read chain A from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-c2') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) chainb
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read chain A from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-beta1') then
      beta1 = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-beta2') then
      beta2 = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-ocup1') then
      ocup1 = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-ocup2') then
      ocup2 = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-rmin1') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) rmin1
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read rmin1 from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-rmax1') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) rmax1
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read rmax1 from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-rmin2') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) rmin2
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read rmin2 from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-rmax2') then
      call getarg(i+1,value)
      read(value,*,iostat=ioerr) rmax2
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not read rmax2 from command line. '
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-pdblist') then
      call getarg(i+1,value)
      pdblist = value
    else if(keyword(1:length(keyword)).eq.'-skip') then
      skip = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-m') then
      call getarg(i+1,value)
      method = ival(value)
      if(method.ne.1.and.&
         method.ne.2.and.&
         method.ne.3.and.&
         method.ne.4) then
        write(*,*) ' ERROR: Wrong method specification (1, 2, 3 or 4)'
        stop
      end if
    else if(keyword(1:length(keyword)).eq.'-g') then
      call getarg(i+1,value)
      gap = dval(value)
    else if(keyword(1:length(keyword)).eq.'-dtri') then
      call getarg(i+1,value)
      dtri = dval(value)
    else if(keyword(1:length(keyword)).eq.'-gdt_threshold') then
      call getarg(i+1,value)
      gdt_threshold = dval(value)
    else if(keyword(1:length(keyword)).eq.'-maxit') then
      call getarg(i+1,value)
      maxit = ival(value)
    else if(keyword(1:length(keyword)).eq.'-print') then 
      call getarg(i+1,value)
      iprint = ival(value)
    else if(keyword(1:length(keyword)).eq.'-o') then
      call getarg(i+1,value)
      output = .true.
      pdbout = value(1:length(value))
    else if(keyword(1:length(keyword)).eq.'-rmsf') then
      call getarg(i+1,value)
      rmsf = .true.
      rmsfout = value(1:length(value))
    else if(keyword(1:length(keyword)).eq.'-rmsftrend') then
      call getarg(i+1,value)
      rmsftrend = .true.
      rmsftrendout = value(1:length(value))
    else if(keyword(1:length(keyword)).eq.'-all') then
      all = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-seqoff') then
      seqoff = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-nglobal') then
      call getarg(i+1,value)
      nglobal = ival(value)
    else if(keyword(1:length(keyword)).eq.'-maxtrial') then
      call getarg(i+1,value)
      maxtrial = ival(value)
    else if(keyword(1:length(keyword)).eq.'-seqfix') then
      seqtype = 1
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-seqnum') then
      seqtype = 2
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-fasta') then
      seqtype = 3
      call getarg(i+1,value)
      fastafile = value(1:length(value))
      call readfasta()
    else if(keyword(1:length(keyword)).eq.'-noini') then
      useini = .false.
      i = i - 1
    else
      write(*,*) ' Unrecognized command line argument:',&
                 keyword(1:length(keyword))
      stop
    end if
    i = i + 2
  end do

  if(gap.gt.1.d20) then
    if(method.eq.1) gap = 10. 
    if(method.eq.2) gap = 0. 
    if(method.eq.3) gap = 0.
    if(method.eq.4) gap = 0.
  end if 

  if(iprint.eq.1) call title()

end subroutine getpars
