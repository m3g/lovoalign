!
! Program LOVOALIGN: Low Order Value Optimization Methods for Protein
!                    Alignment.
!
!         Copyright:
!         Leandro Martinez, Roberto Andreani, Jose Mario Martinez
!
!         Department of Applied Mathematics, IMECC, University of
!         Campinas and Institute of Chemistry, University of
!         Campinas, Brazil. 
!
! Home-Page: http://www.ime.unicamp.br/~martinez/lovoalign
!            Contains the software, instructions, regular updates.
!
! Primary reference, please cite this work when using LOVOALIGN:
!
!            L. Martinez, R. Andreani, J. M. Martinez. 
!            Convergent algorithms for protein structural alignment
!            BMC Bioinformatics, 8:306 2007
!
! Related work:
!
!            R. Andreani, J. M. Martinez, L. Martinez, F. Yano.
!            Continuous Optimization Methods for Structural Alignment.
!            Mathematical Pogramming, 112:93-124, 2008
!
!            R. Andreani, J. M. Martinez, L. Martinez
!            Trust region superposition methods for protein alignment.
!            IMA Journal on Numerical Analysis, 28, 690-710, 2008.
!
!            R. Andreani, J. M. Martinez, L. Martinez, F. Yano.
!            Low Order Value Optimization and Applications.
!            Journal of Global Optimization, 43, 1-22, 2009.
!
! Dimensions: maxatom:  Maximum number of CA atoms of each protein
!             maxfiles: Maximum number of files in the pdb file list
!                       for database comparisons
!
! Options:
!    -p1 [filename]        Set file of first protein
!    -p2 [filename]        Set file of target protein
!    -c1 [character]       Chain to be considered on first protein
!    -c2 [character]       Chain to be considered on second protein
!    -beta1                Consider atoms with beta > 0 on first protein
!    -beta2                Consider atoms with beta > 0 on second protein
!    -ocup1                Consider atoms with occupancy > 0 on first protein
!    -ocup2                Consider atoms with occupancy > 0 on second protein
!    -rmin1 [integer]      Consider only residues of first protein with number greater than rmin1
!    -rmax1 [integer]      Consider only residues of first protein with number smaller than rmax1
!    -rmin2 [integer]      Consider only residues of second protein with number greater than rmin1
!    -rmax2 [integer]      Consider only residues of second protein with number smaller than rmax1
!    -m 1                  Maximize STRUCTAL score
!    -m 2                  Maximize TM-SCORE
!    -m 3                  Optimize TRIANGULAR score
!    -m 4                  Optimize the NON-BIJECTIVE TRIANGULAR score
!    -pdblist [filename]   Sets file containing protein list
!    -noini                Do not use pseudoprotein initial point
!    -maxit [integer]      Maximum number of outer iterations
!    -g [real]             Penalization for gaps
!    -print 0 or 1         Concise or extensive output
!    -all                  Consider all atoms (not only CAs)
!    -o [filename]         Sets the name for output file 
!    -dtri [real]          Cutoff distance for triangular score (default: 3.d0)
!                          and for small RMSD cutoff output for other methods.
!    -gdt_threshold [real] Threshold distance for GDT scores (default: 4.d0)
!    -seqoff               Do not write sequence alignment
!    -seqfix               Use a fixed sequence alignment (1-1,2-2,...)
!    -rmsf [filename]      Write rmsf profile to file.
!    -rmsftrend [filename] Write rmsf profile trend to file (fraction of pairs
!                          with rmsf smaller than threshold)
!    -nglobal              Number of times the best point must be found to be
!                          accepted to be the global minimizer (default=3)
!

program lovoalign

  use sizes
  use inputpars
  use ioformat
  use file_operations
  use initrandom
  implicit none
  integer :: i, j, narg, mode, &
             na, nb, iargc, &
             nfiles, numa(maxatom), numb(maxatom), &
             indisord(maxatom-1,maxatom)
  double precision :: prota(maxatom,3), protb(maxatom,3), &
            pseudoa(maxatom,3), pseudob(maxatom,3), &
            disord(maxatom-1,maxatom) 
  real :: etime, tarray(2), time0 
  logical :: error
  character(len=1) :: resa(maxatom), resb(maxatom)
  character(len=200) :: record, pdbfiles(maxfiles), file1, file2

  ! Computing running time

  time0 = etime(tarray)
       
  ! The LOVOALIGN program can be run in three modes: 
  ! Mode 0: Align two protein structures.
  ! Mode 1: Align a protein structure to a database.
  ! Mode 2: Perform an all-on-all structural alignment in a database.

  mode = 0
  narg = iargc()
 
  ! If there are no command line arguments, print instructions

  if(narg.eq.0) call help
                 
  ! If the '-pdblist' was found in the command line, the program will run in
  ! modes 1 or 2. If no specific coordinate file was set with '-p1' this
  ! will be a all-on-all comparison. If '-p1' and '-pdblist' were set this
  ! will be a comparison of one protein to a database
                                  
  do i = 1, narg
    call getarg(i,record)
    if(record(1:length(record)).eq.'-pdblist') mode = 2
    if(record(1:length(record)).eq.'-h'.or. &
       record(1:length(record)).eq.'-help'.or. &
       record(1:length(record)).eq.'--help') call help
  end do
  if(mode.eq.2) then
    do i = 1, narg
      call getarg(i,record)
      if(record(1:length(record)).eq.'-p1') mode = 1
    end do
  end if
  
  ! Default parameters

  maxit = 10000
  iprint = 0
  output = .false.
  chaina = '#'
  chainb = '#'
  beta1 = .false.
  beta2 = .false.
  ocup1 = .false.
  ocup2 = .false.
  rmin1 = 0
  rmax1 = 1000000
  rmin2 = 0
  rmax2 = 1000000
  all = .false.
  seqoff = .false.
  seqfix = .false.
  useini = .true.
  dtri = 3.d0
  gdt_threshold = 4.d0
  nglobal = 3
  if ( mode == 0 ) then
    maxtrial = 1000
  else
    maxtrial = 4
  end if
 
  ! Default method

  method = 2

  ! Initialize random number generator
  seed = 1234567
  call init_random_number(seed)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                    !
  ! Mode 0: Align two protein structures (A and B)     !
  !                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(mode.eq.0) then

    ! The default here is to use extensive printing
      
    iprint = 1

    ! Read command line parameters

    if(narg.gt.0) call getpars()
    max_filename_size = max(length(remove_path(protea)),length(remove_path(proteb)))
    call outputformats()
    if(iprint.eq.0) then
      write(*,header_list) trim(adjustl(basename(protea))), &
                           trim(adjustl(basename(proteb))), "none", 1, dtri
    end if

    ! Read protein coordinates

    call readfile(protea,prota,chaina,beta1,ocup1,&
                  rmin1,rmax1,na,resa,numa,all,error)
    call readfile(proteb,protb,chainb,beta2,ocup2,&
                  rmin2,rmax2,nb,resb,numb,all,error)
    if(error) stop

    ! Printing all data from this run:
      
    if(iprint.ge.1) call printdata(protea,proteb,na,nb,chaina,&
                                   chainb,method,gap,maxit,dtri,gdt_threshold,&
                                   useini)

    ! Compute pseudoproteins A and B, d1 is the distance between atom i and
    ! atom i+2

    if(useini) call pseudoprot(prota,pseudoa,na)
    if(useini) call pseudoprot(protb,pseudob,nb)

    ! Do dynamic programming with these pseudoproteins and return the
    ! proteins ready for structural alignment

    if(useini) then
      call initial(pseudoa,pseudob,na,nb,prota,protb,seqfix)
    else
      call tocm(prota,protb,na,nb)
    end if

    ! If using non-bijective method, sort internal distances of protein B

    if(method.eq.4) call orprot(nb,protb,disord,indisord)

    ! Performing the structural alignment
                                        
    call protall(prota,protb,na,nb,disord,indisord,resa,resb,numa,numb)

    ! If required, print output file with protein A aligned to protein B

    if(output) then
      call writepdb(pdbout,protea,prota,chaina,beta1,ocup1,&
                    rmin1,rmax1,na,resa,numa,proteb,all)
      if(iprint.ge.1) then
        write(*,*) ' Wrote file: ', trim(adjustl(pdbout))
        write(*,dash_line)
      end if
    end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                    !
  ! Mode 1: Align a protein to a database              !
  !                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else if(mode.eq.1) then

    ! Get parameters from the command line

    call getpars()

    ! Read file of protein B (the specified protein will be allways B)

    call readfile(protea,protb,chainb,beta1,ocup1,&
                  rmin1,rmax1,nb,resb,numb,all,error)
    proteb = protea
    if(error) stop

    ! Compute pseudoprotein B for initial point

    if(useini) call pseudoprot(protb,pseudob,nb)

    ! Compute ordered distance matrix for protein B for non-bijective alignments

    if(method.eq.4) call orprot(nb,protb,disord,indisord)

    ! Read file list and align proteins

    call readlist(pdblist,pdbfiles,nfiles)
    max_filename_size = max(max_filename_size,length(remove_path(proteb)))
    call outputformats()
    if(iprint.eq.0) write(*,header_list) "from list", &
                                         trim(adjustl(proteb)), &
                                         trim(adjustl(pdblist)), &
                                         nfiles, dtri

    ! Start alignments

    file1 = trim(adjustl(proteb))
    do i = 1, nfiles 

      ! Read file for protein A

      protea = pdbfiles(i)

      ! Skip this alignment if skip was set and the file was
      ! not reached, if reached, will compute from the next one

      if ( skip ) then
        file2 = pdbfiles(i)
        file2 = trim(adjustl(file2))
        if ( file2 == file1 ) then
          skip = .false. 
        end if
        cycle
      end if

      call readfile(protea,prota,chaina,beta2,ocup2,&
                    rmin2,rmax2,na,resa,numa,all,error)

      if(.not.error) then
 
        ! If specificaly required, print data for this problem
      
        if(iprint.ge.1) call printdata(protea,proteb,na,nb,chaina,&
                                       chainb,method,gap,maxit,dtri,gdt_threshold,&
                                       useini)

        ! Compute pseudoprotein A for initial point
       
        if(useini) call pseudoprot(prota,pseudoa,na)

        ! Compute initial point

        if(useini) then
          call initial(pseudoa,pseudob,na,nb,prota,protb,seqfix)
        else
          call tocm(prota,protb,na,nb)
        end if

        ! Performing the structural alignment
                                        
        call protall(prota,protb,na,nb,disord,indisord,resa,resb,numa,numb)

      end if

    end do 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                        !
  ! Mode 2: Perform an all-on-all comparison in a database !
  !                                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else if(mode.eq.2) then 

    ! Get parameters from the command line

    call getpars()

    ! Read the list of pdb files

    call readlist(pdblist,pdbfiles,nfiles)
    call outputformats()
    if(iprint.eq.0) write(*,header_list) "from list", &
                                         "from list", &
                                         trim(adjustl(pdblist)), &
                                         nfiles, dtri

    ! For non-bijective methods, get the number of atoms of all proteins
    ! and order them from biggest to smallest. 

    if(method.eq.4) call orderpdb(pdbfiles,nfiles)

    ! Perform the alignment

    do i = 1, nfiles - 1

      ! Read protein file B

      proteb = pdbfiles(i)
      call readfile(proteb,protb,chainb,beta1,ocup1,&
                    rmin1,rmax1,nb,resb,numb,all,error)

      if(.not.error) then

        ! Compute pseudoprotein B for initial point

        if(useini) call pseudoprot(protb,pseudob,nb)

        ! Computing ordered distance matrix for protein B for non-bijective method

        if(method.eq.4) call orprot(nb,protb,disord,indisord)

        ! Running over all other proteins
     
        do j = i + 1, nfiles

          ! Read file for protein A

          protea = pdbfiles(j)
          call readfile(protea,prota,chaina,beta2,ocup2,&
                        rmin2,rmax2,na,resa,numa,all,error)

          if(.not.error) then

            ! If specificaly required, print data for this problem
      
            if(iprint.ge.1) call printdata(protea,proteb,na,nb,chaina,&
                                           chainb,method,gap,maxit,dtri,gdt_threshold,&
                                           useini)
 
            ! Compute pseudoprotein A for obtaining initial point

            if(useini) call pseudoprot(prota,pseudoa,na)

            ! Compute initial point
            
            if(useini) then 
              call initial(pseudoa,pseudob,na,nb,prota,protb,seqfix)
            else
              call tocm(prota,protb,na,nb)
            end if

            ! Performing the structural alignment
                                        
            call protall(prota,protb,na,nb,disord,indisord,resa,resb,numa,numb)

          end if
        end do
      end if
    end do
  end if
 
  ! Compute total running time

  time0 = etime(tarray) - time0
  if(iprint.ge.1) then
    write(*,*) ' TOTAL RUNNING TIME: ', time0
    write(*,dash_line)
  else
    write(*,"(a)") trim(hash_line)
    write(*,"('# TOTAL RUNNING TIME:',i10,' h ',i6,' min ',f12.4,' s')") &
            int(time0/3600),&
            int(mod(time0,3600.)/60),&
            mod(time0,60.)
    write(*,"(a)") trim(hash_line)
  end if
 
end program lovoalign
