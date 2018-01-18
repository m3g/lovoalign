!
! Subroutine that reads protein coordinates both in the pdb format or in
! simple clean coordinates file used for tests or other purposes
!

subroutine readfile(protea,prota,chaina,beta1,ocup1,rmin1,rmax1,na,resa,numa,all,error)

  use sizes
  implicit none
  integer :: na, length, ic, numa(maxatom),ioerr
  integer :: rmin1, rmax1, resnum
  real :: beta, occupancy
  double precision :: prota(maxatom,3)
  logical :: all, error, beta1, ocup1
  character(len=1) :: chaina, resa(maxatom), letter
  character(len=3) :: resid
  character(len=200) :: protea, record

  error = .false. 

  !
  ! Trying to read coordinates in pdb format
  !

  open(10,file=protea(1:length(protea)),status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file:',protea(ic(protea):length(protea))
    error=.true.
    return
  end if

  na = 0
  do while(.true.)
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(na.ge.1.and.record(1:3).eq.'END') exit

    ! Reading coordinates for all atoms (not only CAs), if specified

    if(all) then
      if((record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM').and.&
         (record(22:22).eq.chaina.or.chaina.eq.'#')) then
       
        if(ocup1) read(record(55:60),*,iostat=ioerr) occupancy
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Tried to read occupancy from file: ',&
                     protea(ic(protea):length(protea)),&
                     ' and failed. '
          write(*,*) '        However, the Occupancy option was set. '
          error = .true.
          return
        end if

        if(beta1) read(record(61:66),*,iostat=ioerr) beta
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Tried to read beta factor from file: ',&
                     protea(ic(protea):length(protea)),&
                     ' and failed. '
          write(*,*) '        However, the Beta option was set. '
          error = .true.
          return
        end if

        read(record(23:26),*,iostat=ioerr) resnum 
        if ( ioerr == 0 ) then
          if ( resnum < rmin1 .or. resnum > rmax1 ) cycle
        end if

        if( ( .not.ocup1 .and. .not.beta1 ) .or.&
            ( ocup1 .and. occupancy.gt.0. ) .or.&
            ( beta1 .and. beta.gt.0. ) ) then 
          na = na + 1
          if(na.gt.maxatom) then 
            write(*,*) ' ERROR: ',protea(ic(protea):length(protea)),&
                       ' atoms exceed MAXATOM.'
            error=.true.
            return
          end if
          read(record(31:38),*,iostat=ioerr) prota(na,1)
          if ( ioerr /= 0 ) error = .true.
          read(record(39:46),*,iostat=ioerr) prota(na,2)
          if ( ioerr /= 0 ) error = .true.
          read(record(47:54),*,iostat=ioerr) prota(na,3)
          if ( ioerr /= 0 ) error = .true.
          if ( error ) then
            write(*,*) ' ERROR: Failed reading coordinates in file ',&
                       protea(ic(protea):length(protea))
            return
          end if

          ! Reading residue name and getting one-letter code
      
          read(record(18:20),*,iostat=ioerr) resid
          if ( ioerr /= 0 ) resid = 'XXX'
          resa(na) = letter(resid)

          ! Reading the residue number

          numa(na) = resnum
          read(record(23:26),*,iostat=ioerr) numa(na)
          if ( ioerr /= 0 ) numa(na) = na

        end if
      end if

      ! Default: Reading coordinates for CA atoms only

    else 
      if(record(1:4).eq.'ATOM'.and.&
         (record(14:15).eq.'CA'.or.record(13:14).eq.'CA').and.&
         (record(22:22).eq.chaina.or.chaina.eq.'#')) then

        if(ocup1) read(record(55:60),*,iostat=ioerr) occupancy
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Tried to read occupancy from file: ',&
                     protea(ic(protea):length(protea)),&
                     ' and failed. '
          write(*,*) '        However, the Occupancy option was set. '
          error = .true.
          return
        end if

        if(beta1) read(record(61:66),*,iostat=ioerr) beta
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Tried to read beta factor from file: ',&
                     protea(ic(protea):length(protea)),&
                     ' and failed. '
          write(*,*) '        However, the Beta option was set. '
          error = .true.
          return
        end if

        read(record(23:26),*,iostat=ioerr) resnum 
        if ( ioerr == 0 ) then
          if ( resnum < rmin1 .or. resnum > rmax1 ) cycle
        end if

        if( ( .not.ocup1 .and. .not.beta1 ) .or.&
            ( ocup1 .and. occupancy.gt.0. ) .or.&
            ( beta1 .and. beta.gt.0. ) ) then 
          na = na + 1
          if(na.gt.maxatom) then 
            write(*,*) ' ERROR: ',protea(ic(protea):length(protea)),&
                       ' Calpha atoms exceed MAXATOM.'
            error=.true.
            return
          end if                              
          read(record(31:38),*,iostat=ioerr) prota(na,1)
          if ( ioerr /= 0 ) error = .true.
          read(record(39:46),*,iostat=ioerr) prota(na,2)
          if ( ioerr /= 0 ) error = .true.
          read(record(47:54),*,iostat=ioerr) prota(na,3)
          if ( ioerr /= 0 ) error = .true.
          if ( error ) then
            write(*,*) ' ERROR: Failed reading coordinates in file ',&
                       protea(ic(protea):length(protea))
            return
          end if
        end if

        ! Reading residue name and getting one-letter code
      
        read(record(18:20),*,iostat=ioerr) resid
        if ( ioerr /= 0 ) resid = 'XXX'
        resa(na) = letter(resid)

        ! Reading the residue number

        read(record(23:26),*,iostat=ioerr) numa(na)
        if ( ioerr /= 0 ) numa(na) = na

      end if
    end if
  end do
  close(10)

  ! If coordinates where not found

  if(na.eq.0) then
    write(*,*) ' ERROR: Could not read coordinates from',&
               ' file: ', protea(ic(protea):length(protea)),','
    write(*,*) '        or some selection has no atoms.' 
    error=.true.
    return
  end if

  ! With less than tree CA atoms we are not talking about a protein

  if(na.lt.3) then
    write(*,*) ' ERROR: Protein with less than three atoms:',&
               protea(ic(protea):length(protea))
    error=.true.
    return
  end if

! If the number of atoms is greater than maxatom, report error

  if(na.gt.maxatom) then
    write(*,*) ' ERROR: Number of atoms to be read in file',&
               protea(ic(protea):length(protea)),&
               ' is greater than maxatom. '
    error = .true.
  end if

end subroutine readfile
