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
!

! Module that set maximum problem size

module sizes

  integer, parameter :: maxatom = 4500
  integer, parameter :: maxfiles = 3000

end module sizes

! Module containing ahestetic formats

module ioformat

  implicit none
  character(len=*), parameter :: dash_line = "('  ',71('-'))"

end module ioformat

!
! Main program
!

program lovoalign

  use sizes
  use ioformat
  implicit none
  integer :: i, j, narg, mode, method, maxit, iprint, &
             na, nb, length, ic, iargc, &
             nfiles, numa(maxatom), numb(maxatom), &
             indisord(maxatom-1,maxatom)
  double precision :: gap, prota(maxatom,3), protb(maxatom,3), &
            pseudoa(maxatom,3), pseudob(maxatom,3), dtri, &
            disord(maxatom-1,maxatom), gdt_threshold
  real :: etime, tarray(2), time0 
  logical :: all, error, output, useini, seqoff, seqfix
  character(len=1) :: chaina, chainb, resa(maxatom), resb(maxatom)
  character(len=200) :: record, protea, proteb, pdblist, pdbfiles(maxfiles), &
                        pdbout
  character(len=1000) :: header_list
  logical :: beta1, beta2, ocup1, ocup2

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
  all = .false.
  seqoff = .false.
  seqfix = .false.
  useini = .true.
  dtri = 3.d0
  gdt_threshold = 4.d0
 
  ! Default method

  method = 2

  ! Header for concise output printing

  header_list = "(111('#'),/,&
           & '# LOVOALIGN ',/&
           &,'# http://www.ime.unicamp.br/~martinez/lovoalign',/,&
           &111('#'),/,&
           &'# Prot A: Variable protein.',/,&
           &'# Prot B: Target (fixed) protein.',/,&
           &'# SCORE: Best score obtained.',/,&
           &'# COV: Coverage (number of corresponding atoms).',/,&
           &'# RMSD: Root mean square deviation of COV atoms.',/,&
           &'# COV2: Number of atoms closer than ',f8.3,' Angstroms.',/,&
           &'# RMSD2: Root mean square deviation of COV2 atoms.',/,&
           &'# GDT_TM: Global Distance Test (GDT) score.',/,&
           &'# GDT_HA: High-accuracy GDT score.',/,&
           &'# TIME: Time used in this alignment.',/,&
           &111('#'),/,&
           &'# Prot A',t16,'Prot B',t35,'SCORE',t46,'COV',t58,'RMSD',&
           &t64,'COV2',t76,'RMSD2',t84,'GDT_TM',t93,'GDT_HA',t108,'TIME')"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                    !
  ! Mode 0: Align two protein structures (A and B)     !
  !                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(mode.eq.0) then

    ! The default here is to use extensive printing
      
    iprint = 1

    ! Read command line parameters

    if(narg.gt.0) call getpars(method,gap,maxit,dtri,gdt_threshold,&
                               protea,proteb,chaina,chainb,iprint,&
                               all,seqoff,seqfix,pdblist,output,pdbout,&
                               useini,beta1,beta2,ocup1,ocup2)
    if(iprint.eq.0) write(*,header_list) dtri

    ! Read protein coordinates

    call readfile(protea,prota,chaina,beta1,ocup1,&
                 na,resa,numa,all,error)
    call readfile(proteb,protb,chainb,beta2,ocup2,&
                  nb,resb,numb,all,error)
    if(error) stop

    ! Printing all data from this run:
      
    if(iprint.eq.1) call printdata(protea,proteb,na,nb,chaina,&
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
                                        
    call protall(prota,protb,na,nb,method,gap,maxit,dtri,gdt_threshold,&
                 iprint,disord,indisord,&
                 protea,proteb,resa,resb,numa,numb,seqoff,seqfix)

    ! If required, print output file with protein A aligned to protein B

    if(output) then
      call writepdb(pdbout,protea,prota,chaina,beta1,ocup1,&
                    na,resa,numa,proteb,all)
      if(iprint.eq.1) then
        write(*,*) ' Wrote file: ', pdbout(ic(pdbout):length(pdbout))
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

    call getpars(method,gap,maxit,dtri,gdt_threshold,&
                 protea,proteb,chaina,chainb,iprint,&
                 all,seqoff,seqfix,&
                 pdblist,output,pdbout,useini,&
                 beta1,beta2,ocup2,ocup2)
    if(iprint.eq.0) write(*,header_list) dtri

    ! Read file of protein B (the specified protein will be allways B)

    call readfile(protea,protb,chainb,beta1,ocup1,&
                  nb,resb,numb,all,error)
    proteb = protea
    if(error) stop

    ! Compute pseudoprotein B for initial point

    if(useini) call pseudoprot(protb,pseudob,nb)

    ! Compute ordered distance matrix for protein B for non-bijective alignments

    if(method.eq.4) call orprot(nb,protb,disord,indisord)

    ! Read file list and align proteins

    call readlist(pdblist,pdbfiles,nfiles)

    ! Start alignments

    do i = 1, nfiles 

      ! Read file for protein A

      protea = pdbfiles(i)
      call readfile(protea,prota,chaina,beta2,ocup2,&
                    na,resa,numa,all,error)

      if(.not.error) then
 
        ! If specificaly required, print data for this problem
      
        if(iprint.eq.1) call printdata(protea,proteb,na,nb,chaina,&
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
                                        
        call protall(prota,protb,na,nb,method,&
                     gap,maxit,dtri,gdt_threshold,&
                     iprint,disord,indisord,&
                     protea,proteb,resa,resb,numa,numb,&
                     seqoff,seqfix)

      end if

    end do 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                        !
  ! Mode 2: Perform an all-on-all comparison in a database !
  !                                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else if(mode.eq.2) then 

    ! Get parameters from the command line

    call getpars(method,gap,maxit,dtri,gdt_threshold,&
                 protea,proteb,chaina,chainb,iprint,&
                 all,seqoff,seqfix,&
                 pdblist,output,pdbout,useini,&
                 beta1,beta2,ocup1,ocup2)
    if(iprint.eq.0) write(*,header_list) dtri

    ! Read the list of pdb files

    call readlist(pdblist,pdbfiles,nfiles)

    ! For non-bijective methods, get the number of atoms of all proteins
    ! and order them from biggest to smallest. 

    if(method.eq.4) call orderpdb(pdbfiles,nfiles)

    ! Perform the alignment

    do i = 1, nfiles - 1

      ! Read protein file B

      proteb = pdbfiles(i)
      call readfile(proteb,protb,chainb,beta1,ocup1,&
                    nb,resb,numb,all,error)

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
                        na,resa,numa,all,error)

          if(.not.error) then

            ! If specificaly required, print data for this problem
      
            if(iprint.eq.1) call printdata(protea,proteb,na,nb,chaina,&
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
                                        
            call protall(prota,protb,na,nb,method,&
                         gap,maxit,dtri,gdt_threshold,&
                         iprint,disord,indisord,protea,proteb,&
                         resa,resb,numa,numb,seqoff,seqfix)

          end if
        end do
      end if
    end do
  end if
 
! Compute total running time

  time0 = etime(tarray) - time0
  if(iprint.eq.1) then
    write(*,*) ' TOTAL RUNNING TIME: ', time0
    write(*,dash_line)
  else
    write(*,"(111('#'),/,&
              &'# TOTAL RUNNING TIME:',i10,' h ',i6,' min ',f12.4,' s',/,&
              &111('#'))") int(time0/3600),&
                          int(mod(time0,3600.)/60),&
                          mod(time0,60.)
  end if
 
end program lovoalign

!!!!!!!!!!!!!!!!!!!!!!!!
!                      !
! END OF MAIN PROGRAM  !
!                      !
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
! SUBROUTINE PROTALL: This is the main subroutine of the program. The !
! actual subroutine that performs the alignment between two proteins. !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
subroutine protall(prota,protb,na,nb,method,&
                   gap,maxit,dtri,gdt_threshold,&
                   iprint,disord,indisord,&
                   protea,proteb,resa,resb,numa,numb,&
                   seqoff,seqfix)

  use sizes
  use ioformat
  implicit none

  double precision :: dzero, gap, prota(maxatom,3), protb(maxatom,3),&
                      bijscore(maxatom), gnor,&
                      score, dzero2, tol, scale,&
                      prevscore, rmsd, rmsd2, dtri, dtri2,&
                      disord(maxatom-1,maxatom), gdt_threshold,&
                      gdt_tm, gdt_ha
  real :: etime, tarray(2), time1
  integer :: maxit, na, nb,&
             bije(maxatom,2),&
             ngaps, nbij, method, nef, nbij_dtri,&
             iprint, length, ic, numa(maxatom),&
             numb(maxatom), it,&
             indisord(maxatom-1,maxatom), pair(maxatom)
  character(len=1) :: resa(maxatom), resb(maxatom)
  character(len=200) :: protea, proteb
  character(len=200) :: title_format, data_format 
  logical :: seqoff, seqfix
  external :: structal, tmscore

  title_format = "(t3,'ITER',t20,'SCORE',t30,'GRADIENT NORM',&
                  &t45,'COVERAGE',t56,'GAPS',t64,'NEF')"
  data_format = "(i6,tr1,e17.10,tr1,e17.10,tr4,i6,tr1,i6,tr1,i6)"

  ! Time used in this alignment is computed from here

  time1 = etime(tarray)

  ! This is the relative precision for score convergence

  tol = 1.d-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                !
  ! Method 1: Maximize the STRUCTAL score of Gerstein and Levitt   !
  !                                                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  if(method.eq.1) then

    ! Writes the titles for regular iteration printing

    if(iprint.eq.1) write(*,title_format)

    ! Normalization of STRUCTAL score
  
    dzero = 2.24d0
    dzero2 = dzero*dzero
    scale = 20.d0
 
    ! Number iterations, functional evaluations and DP calls

    nef = 0 
    it = 0
    prevscore = 0.d0
 
    ! Compute the DP bijection and the score at initial point          

    call structal(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                  bijscore,ngaps,score,seqfix)

    nef = nef + 1
    if(iprint.eq.1) write(*,data_format) it, score, 0.d0, nbij, ngaps,nef
          
    ! Here begin the iteration loop

    do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

      it = it + 1
      prevscore = score

      ! Perform a newton step to maximize the score given the bijection

      call newton(structal,na,nb,prota,protb,score,bije,bijscore,&
                  dzero2,scale,nbij,gap,ngaps,nef,gnor,seqfix)

      ! Output regular iteration data
          
      if(iprint.eq.1) write(*,data_format) it, score, gnor, nbij, ngaps,nef

    end do
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                !
  ! Method 2: Maximize the TM-SCORE of Zhang and Skolnick          !
  !                                                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(method.eq.2) then

    ! Checks if protein B has more than 15 atoms, otherwise stop (because
    ! of the deffinition of dzero, below)

    if(nb.le.15) then
      write(*,*) ' ERROR: For using the TM-SCORE, the number of '
      write(*,*) '        atoms of protein B must be at least 15.'
      stop
    end if

    ! Writes the titles for regular iteration printing

    if(iprint.eq.1) write(*,title_format)

    ! Normalization of the TM-SCORE score

    dzero = 1.24d0 * (nb-15.d0)**(1.d0/3.d0) - 1.8d0
    dzero2 = dzero*dzero
    scale = 1.d0 / dfloat(nb)

    ! Number iterations, functional evaluations and DP calls

    nef = 0 
    it = 0
    prevscore = 0.d0
 
    ! Compute the DP bijection and the score at initial point          

    call tmscore(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                 bijscore,ngaps,score,seqfix)
    nef = nef + 1
    if(iprint.eq.1) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
          
    ! Here begin the iteration loop

    do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

      it = it + 1
      prevscore = score

      ! Perform a newton step to maximize the score given the bijection

      call newton(tmscore,na,nb,prota,protb,score,bije,bijscore,&
                  dzero2,scale,nbij,gap,ngaps,nef,gnor,seqfix)

      ! Output regular iteration data
          
      if(iprint.eq.1) write(*,data_format) it, score, gnor, nbij, ngaps, nef

    end do
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                !
  ! Method 3: Maximize the TRIANGULAR score                        !
  !                                                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(method.eq.3) then

    ! Square of the minimum distance
  
    dtri2 = dtri*dtri
 
    ! Writes the titles for regular iteration printing

    if(iprint.eq.1) write(*,title_format)

    ! Number iterations, functional evaluations and DP calls

    nef = 0 
    it = 0
    prevscore = 0.d0
 
    ! Compute the correspondence and the score at initial point          

    call triang(prota,protb,na,nb,dtri2,gap,bije,nbij,&
                bijscore,ngaps,score,seqfix)
    nef = nef + 1
    if(iprint.eq.1) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
          
    ! Here begin the iteration loop

    do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

      it = it + 1
      prevscore = score
 
      ! Given the correspondence, perform a Procrustes RMSD alignment

      call procrustes(nbij,na,bije,prota,protb)

      ! Compute the DP bijection and the score at the new orientation

      call triang(prota,protb,na,nb,dtri2,gap,bije,nbij,bijscore,ngaps,&
                  score,seqfix)

      ! Output regular iteration data
          
      if(iprint.eq.1) write(*,data_format) it, score, score-prevscore, nbij, ngaps, nef

    end do

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
! Method 4: Maximize the NON-BIJECTIVE TRIANGULAR score          !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(method.eq.4) then

    ! Square of the minimum distance
  
    dtri2 = dtri*dtri
 
    ! Writes the titles for regular iteration printing

    if(iprint.eq.1) write(*,title_format)

    ! Number iterations and functional evaluations

    nef = 0 
    it = 0
    prevscore = 0.d0 
 
    ! Compute the correspondence and the score at initial point          

    call nonbscore(na, nb, prota, protb, dtri2, gap,&
                   disord, indisord, it, pair,&
                   score, bije, nbij, ngaps)
    nef = nef + 1
    if(iprint.eq.1) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
          
    ! Here begin the iteration loop

    do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

      it = it + 1
      prevscore = score
 
      ! Given the correspondence, perform a Procrustes RMSD alignment

      call procrustes(nbij,na,bije,prota,protb)

      ! Compute non-bijective correspondence and the score at the new orientation

      call nonbscore(na, nb, prota, protb, dtri2, gap,&
                     disord, indisord, it, pair,& 
                     score, bije, nbij, ngaps)
      nef = nef + 1

      ! Output regular iteration data
          
      if(iprint.eq.1) write(*,data_format) it, score, score-prevscore, nbij, ngaps, nef

    end do

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                               !
  ! POST ANALYSIS AND REPORT                      !
  !                                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ! Writting the final bijection obtained

  if(iprint.eq.1 .and. .not. seqoff) call writebije(na,nb,resa,resb,numa,numb,bije,nbij)
        
  ! Computing the RMSD of aligned residues at the solution

  call getrmsd(prota,protb,bije,nbij,rmsd)

  ! Computing the GDT scores at the solution

  call computegdt(na,nb,prota,protb,bije,nbij,gdt_threshold,gdt_tm,gdt_ha)
 
  ! Printing the final score

  if(iprint.eq.1) then
    write(*,dash_line)
    write(*,"(a14,tr1,f12.6,tr1,a10,tr1,i5,tr1,a6,tr1,f10.6,tr1,a6,i4)")&
            '  FINAL SCORE:', score,' COVERAGE:', nbij,' RMSD:', rmsd,' GAPS:', ngaps
  endif
 
  ! Compute rmsd for atoms which are closer than some tolerance

  call getrmsd2(prota,protb,bije,nbij,rmsd2,nbij_dtri,dtri)
 
  ! Printing the final score

  if(iprint.eq.1) then
    write(*,dash_line)
    write(*,"(a,f8.4,a,f10.6,a,i6)")&
          '  ATOMS CLOSER THAN ',dtri,' Ang: RMSD: ',rmsd2,' COVERAGE: ', nbij_dtri
    write(*,"(a,f8.3,t34,a,f8.3)")&
          '  GDT_TM SCORE: ', gdt_tm, ' GDT_HA SCORE: ', gdt_ha

  endif

  ! Alignment time

  time1 = etime(tarray) - time1
  
  ! Printing concise output for database comparisons

  if(iprint.eq.0) then
    if(length(protea)-ic(protea)+1.le.10.and.&
       length(proteb)-ic(proteb)+1.le.10) then
      write(*,"(t1,a,t12,a,tr1,f12.6,2(tr1,i5,tr1,f12.6),2(tr1,f8.3),tr1,f12.6)")&
              protea(ic(protea):length(protea)),&
              proteb(ic(proteb):length(proteb)),&
              score, nbij, rmsd, nbij_dtri, rmsd2, gdt_tm, gdt_ha, time1
    else
      write(*,"(t1,a,tr1,a,tr1,f12.6,2(tr1,i5,tr1,f12.6),2(tr1,f8.3),tr1,f12.6)")&
              protea(ic(protea):length(protea)),&
              proteb(ic(proteb):length(proteb)),&
              score, nbij, rmsd, nbij_dtri, rmsd2, gdt_tm, gdt_ha, time1
    end if
  end if
 
  ! Printing final data
 
  if(iprint.eq.1) then
    write(*,dash_line)
    write(*, *)' Time used in this alignment:', time1 
    write(*,dash_line)
    write(*,*) ' END OF ALIGNMENT '
    write(*,dash_line)
  endif

end subroutine protall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                                   !
! Subroutine initial: using pseudoproteins compute the bijection    !
! that maximizes the score using dynamic programming and return     !
! protein A ready for alignment                                     !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial(pseudoa,pseudob,na,nb,prota,protb,seqfix)

  use sizes
  implicit none
  integer :: na, nb, nbij, ngaps, bije(maxatom,2), i
  double precision :: pseudoa(maxatom,3), pseudob(maxatom,3), gap,&
                      dzero2, bijscore(maxatom), score,&
                      prota(maxatom,3), protb(maxatom,3)
  logical :: seqfix

  if ( seqfix ) return

  if(min(na,nb).le.5) then
    write(*,*) ' Too few atoms. Ignoring pseudoprot initial point.'
    return
  end if

  ! Parameters for scoring internal distances

      dzero2 = 100.
      gap = 1.

  ! Initialization based on internal distances and dynamic programming

  call structal(pseudoa,pseudob,na-3,nb-3,dzero2,gap,bije,nbij,&
                bijscore,ngaps,score,seqfix)

  do i = 1, nbij
    bije(i,1) = bije(i,1) + 1
    bije(i,2) = bije(i,2) + 1
  end do

  call procrustes(nbij,na,bije,prota,protb)

end subroutine initial

!
! Subroutine procrustes: given two sets of vectors, finds the best
!                   rotation that to align the two sets. 
!
!                 Method: S. K. Kearsley, 
!                         "On the orthogonal transformation used for
!                          structural comparisons"
!                         Acta Crystallog. A 45 (1989), 208-210 
!
!  Input: nbij: number of points of the bijection
!         bije: array containing the bijection indices
!         yref: reference vector set
!         xvar: variable vector set
!
!  Ouput: xvar vector aligned to y
!         If the original vector is xvar, then the transformed 
!         vector is xvar = u * (var - cmx) + cmy
!

subroutine procrustes(nbij,na,bije,xvar,yref)

  use sizes
  implicit none
  integer :: i, j, k, nbij, iq, na
  integer :: bije(maxatom,2)
  double precision :: x(maxatom,3), y(maxatom,3),&
                      cmx(3), cmy(3), q(4,4), xp(maxatom),&
                      yp(maxatom), zp(maxatom), xm(maxatom),&
                      ym(maxatom), zm(maxatom), qmin,&
                      a(4,4), u(3,3), xvar(maxatom, 3), vecaux(3),&
                      yref(maxatom,3)

  ! Safeguard for the case in which the bijection contains only one atom

  if(nbij.eq.1) then
    vecaux(1) = xvar(bije(1,2),1) - yref(bije(1,1),1)
    vecaux(2) = xvar(bije(1,2),2) - yref(bije(1,1),2)
    vecaux(3) = xvar(bije(1,2),3) - yref(bije(1,1),3)
    do i = 1, na
      xvar(i,1) = xvar(i,1) - vecaux(1)
      xvar(i,2) = xvar(i,2) - vecaux(2)
      xvar(i,3) = xvar(i,3) - vecaux(3)
    end do
    return
  end if

  ! Copying the arrays to have only the atoms of the bijections

  do i = 1, nbij
    x(i,1) = xvar(bije(i,1),1)
    x(i,2) = xvar(bije(i,1),2)
    x(i,3) = xvar(bije(i,1),3)
    y(i,1) = yref(bije(i,2),1)
    y(i,2) = yref(bije(i,2),2)
    y(i,3) = yref(bije(i,2),3)
  end do
         
  ! Computing the centroid of the structures

  do i = 1, 3
    cmx(i) = 0.d0
    cmy(i) = 0.d0
  end do

  do j = 1, 3
    do i = 1, nbij
      cmx(j) = cmx(j) + x(i,j)
      cmy(j) = cmy(j) + y(i,j)  
    end do
  end do

  do i = 1, 3
    cmx(i) = cmx(i) / dfloat(nbij) 
    cmy(i) = cmy(i) / dfloat(nbij)
  end do

  ! Moving structures to their baricenters

  do j = 1, 3
    do i = 1, nbij 
      x(i,j) = x(i,j) - cmx(j)    
      y(i,j) = y(i,j) - cmy(j) 
    end do
  end do

  ! Computing the quaternion matrix

  do i = 1, nbij
    xm(i) = y(i,1) - x(i,1) 
    ym(i) = y(i,2) - x(i,2) 
    zm(i) = y(i,3) - x(i,3) 
    xp(i) = y(i,1) + x(i,1) 
    yp(i) = y(i,2) + x(i,2) 
    zp(i) = y(i,3) + x(i,3)
  end do 
 
  do j = 1, 4
    do i = 1, 4
      q(i,j) = 0.d0
    end do
  end do
 
  do i = 1, nbij  
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)
     
  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix

  call jacobi(q,a,4)

  ! Choosing the quaternion that corresponds to the minimum

  iq = 1
  qmin = q(1,1)
  do i = 2, 4
    if(q(i,i).lt.qmin) then
      iq = i
      qmin = q(i,i)
    end if
  end do
  
  ! Computing the rotation matrix

  u(1,1) = a(1,iq)**2 + a(2,iq)**2 - a(3,iq)**2 - a(4,iq)**2
  u(1,2) = 2. * ( a(2,iq)*a(3,iq) + a(1,iq)*a(4,iq) )
  u(1,3) = 2. * ( a(2,iq)*a(4,iq) - a(1,iq)*a(3,iq) )  
  u(2,1) = 2. * ( a(2,iq)*a(3,iq) - a(1,iq)*a(4,iq) )  
  u(2,2) = a(1,iq)**2 + a(3,iq)**2 - a(2,iq)**2 - a(4,iq)**2 
  u(2,3) = 2. * ( a(3,iq)*a(4,iq) + a(1,iq)*a(2,iq) )  
  u(3,1) = 2. * ( a(2,iq)*a(4,iq) + a(1,iq)*a(3,iq) )  
  u(3,2) = 2. * ( a(3,iq)*a(4,iq) - a(1,iq)*a(2,iq) )  
  u(3,3) = a(1,iq)**2 + a(4,iq)**2 - a(2,iq)**2 - a(3,iq)**2 
 
  ! Rotating-translating the whole protein

  do i = 1, na
    do j = 1, 3
      vecaux(j) = cmy(j)
      do k = 1, 3
        vecaux(j) = vecaux(j) + u(j,k) * (xvar(i,k)-cmx(k))
      end do
    end do
    do j = 1, 3
      xvar(i,j) = vecaux(j)
    end do
  end do

end subroutine procrustes

!
! subroutine tocm: Put the center of mass of protein A
!                  together with the center of mass of protein B,
!                  used as a simple initial point if the 
!                  pseudoprotein initial point is not wanted
!

subroutine tocm(prota,protb,na,nb)

  use sizes
  implicit none
  integer :: na, nb, i
  double precision :: prota(maxatom,3), protb(maxatom,3), cma(3), cmb(3)

  cma(1) = 0.d0
  cma(2) = 0.d0
  cma(3) = 0.d0
  cmb(1) = 0.d0
  cmb(2) = 0.d0
  cmb(3) = 0.d0

  do i = 1, na
    cma(1) = cma(1) + prota(i,1)
    cma(2) = cma(2) + prota(i,2)
    cma(3) = cma(3) + prota(i,3)
  end do
  cma(1) = cma(1) / dfloat(na)
  cma(2) = cma(2) / dfloat(na)
  cma(3) = cma(3) / dfloat(na)

  do i = 1, nb
    cmb(1) = cmb(1) + protb(i,1)
    cmb(2) = cmb(2) + protb(i,2)
    cmb(3) = cmb(3) + protb(i,3)
  end do
  cmb(1) = cmb(1) / dfloat(nb)
  cmb(2) = cmb(2) / dfloat(nb)
  cmb(3) = cmb(3) / dfloat(nb)

  do i = 1, na 
    prota(i,1) = prota(i,1) + ( cmb(1) - cma(1) )
    prota(i,2) = prota(i,2) + ( cmb(2) - cma(2) )
    prota(i,3) = prota(i,3) + ( cmb(3) - cma(3) )
  end do

end subroutine tocm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
! Subroutines required for functional evaluation in different  !
! contexts                                                     !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine moveprot: Moves a protein relative to its initial
!                      orientation given the displacement and
!                      the euler angles of rotation
!
! On input: 
!          x(1)...x(3): Displacement of the center of mass of
!                       the protein to be applied
!          x(4)...x(6): Euler angles of the rotation to be
!                       be applied
!          na: Number of atoms of the protein
!          prota: Original coordinates of the protein
!
! On return:
!          prota: The coordinates rotated and translated
!

subroutine moveprot(x,na,prota)

  use sizes
  implicit none
  integer :: i, na 
  double precision :: x(6), prota(maxatom,3), ca, sa, cb, sb, cg, sg,&
                      xt, yt, zt

  ca = dcos(x(4))
  sa = dsin(x(4))
  cb = dcos(x(5))
  sb = dsin(x(5))
  cg = dcos(x(6))
  sg = dsin(x(6))

  do i = 1, na
    xt = prota(i,1) 
    yt = prota(i,2) 
    zt = prota(i,3) 
    prota(i,1) = x(1) + cb*ca*xt - cb*sa*yt - sb*zt
    prota(i,2) = x(2) + ( cg*sa - sg*sb*ca )*xt &
                      + ( cg*ca + sg*sb*sa )*yt &
                      - sg*cb*zt
    prota(i,3) = x(3) + ( sg*sa + cg*sb*ca )*xt &
                      + ( sg*ca - cg*sb*sa )*yt &
                      + cg*cb*zt
  end do

end subroutine moveprot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                !
! Subroutines for computing the STRUCTAL score   !
!                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine structal: Computes the structal score if the
!                      bijection needs to be determined 
!                      using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dzero2: square of the structal normalization factor
!           gap: penalization for gaps
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           bijscore: individual scores for each pair of the bijection
!           ngaps: number of gaps of the bijection
!           score: the structal score (maximum for current orientation)
!

subroutine structal(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                    bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2)
  double precision :: prota(maxatom,3), protb(maxatom,3), dzero2,&
                      dist, score, scorin(maxatom,maxatom), gap,&
                      bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    score = 0.d0
    do i = 1, na
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      bijscore(i) = 20.d0 / ( 1.d0 + dist / dzero2 )
      score = score + bijscore(i)
      bije(i,1) = i
      bije(i,2) = i
    end do
    nbij = na
    ngaps = 0
    return
  end if

  ! Computes individual scores for all pairs

  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = 20.d0 / ( 1.d0 + dist / dzero2 )
    end do
  end do

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,score)

end subroutine structal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                          !
! Subroutines for computing the TM-SCORE   !
!                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine tmscore: Computes the TM-SCORE if the
!                     bijection needs to be determined 
!                     using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dzero2: square of the structal normalization factor
!           gap: penalization for gaps (usually 0 for TM-SCORE)
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           ngaps: number of gaps of the bijection
!           score: the TM-SCORE (maximum for current orientation)
!

subroutine tmscore(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                   bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2)
  double precision prota(maxatom,3), protb(maxatom,3), dzero2,&
                   dist, score, scorin(maxatom,maxatom), gap,&
                   bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    score = 0.d0
    do i = 1, na
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      score = score + 1.d0 / ( 1.d0 + dist / dzero2 )
      bije(i,1) = i
      bije(i,2) = i
    end do
    score = score / dfloat(na)
    nbij = na
    ngaps = 0
    return
  end if

  ! Computes individual scores for all pairs

  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = 1.d0 / ( 1.d0 + dist / dzero2 )
    end do
  end do

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,score)
  score = score / dfloat(nb)

end subroutine tmscore

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                  !
! Subroutines for computing the TRIANGULAR score   !
!                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine triang: Computes the structal score if the
!                    bijection needs to be determined 
!                    using dynamic programming
!
!  On input: 
!           prota: coordinates of the first protein
!           protb: coordinates of the second protein
!           na: number of atoms of the first protein
!           nb: number of atoms of the second protein
!           dtri2: square of the triangular distance 
!           gap: penalization for gaps
!  On return:
!           bije: the bijection that maximizes the score for the
!                 the current orientation between the proteins
!           nbij: number of corresponding atoms in the bijection
!           bijscore: individual scores for each pair of the bijection
!           ngaps: number of gaps of the bijection
!           score: the triangular score (maximum for current orientation)
!

subroutine triang(prota,protb,na,nb,dtri2,gap,bije,nbij,&
                  bijscore,ngaps,score,seqfix)

  use sizes
  implicit none
  integer :: na, nb, i, j, nbij, ngaps, bije(maxatom,2), npos,&
            nbij_temp, bije_temp(maxatom,2)
  double precision :: prota(maxatom,3), protb(maxatom,3), dtri2,&
                      dist, score, scorin(maxatom,maxatom), gap,&
                      bijscore(maxatom)
  logical :: seqfix

  ! If using a fixed bijection, just compute score and return
  
  if ( seqfix ) then
    score = 0.d0
    do i = 1, na
      dist = (prota(i,1) - protb(i,1))**2 &
           + (prota(i,2) - protb(i,2))**2 &
           + (prota(i,3) - protb(i,3))**2
      score = score + dmax1(0.d0, 1.d0 - dist / dtri2)
      bije(i,1) = i
      bije(i,2) = i
    end do
    nbij = na
    ngaps = 0
    return
  end if

  ! Computes individual scores for all pairs

  npos = 0
  do j = 1, nb
    do i = 1, na
      dist = (prota(i,1) - protb(j,1))**2 &
           + (prota(i,2) - protb(j,2))**2 &
           + (prota(i,3) - protb(j,3))**2
      scorin(i,j) = dmax1(0.d0, 1.d0 - dist / dtri2)
      if(scorin(i,j).gt.0.d0) npos = npos + 1
    end do
  end do

  ! Test if the number of atoms in the bijection is greater than zero

  if(npos.eq.0) then
    nbij_temp = 1
    bije_temp(1,1) = 1
    bije_temp(1,2) = 1
    bijscore(1) = 1.d0
    score = 1.d0
    ngaps = 0
    return
  end if

  ! Perform dynamic programming to obtain the best bijection

  call prodin(na,nb,scorin,gap,ngaps,bije_temp,nbij_temp,&
              bijscore,score)

  ! Include ontly atoms for which the score is positive in the bijection

  nbij = 0
  do i = 1, nbij_temp
    if(bijscore(i).gt.0.d0) then
      nbij = nbij + 1
      bije(nbij,1) = bije_temp(i,1)
      bije(nbij,2) = bije_temp(i,2)
    end if
  end do

end subroutine triang

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
! Subroutines for computing the NON-BIJECTIVE TRIANGULAR score   !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine nonbscore: Subroutine that computes the non-bijective
!                       triangular score.
! 
! On input:
!   na: Number of atoms of protein A
!   nb: Number of atoms of protein B
!   prota: Coordinates of protein A
!   protb: Coordinates of protein B
!   dtri2: Squared threshold distance (above which atoms will
!          not be considered in the correspondence)
!   gap: penalization for gaps
!   disord: Matrix containing internal distances of protein B
!   indisord: Matrix containing the indexes of internal distances
!             of protein B, from the closer to the farther atom,
!             for each atom
!   it: Current iteration of the algorithms
!   pair: The atom o B associated with each atom of A
!
! On return:
!   score: The non-bijective triangular score.
!   bije: The correspondence between atoms obtained.
!   nbij: The number of pairs of the correspondence.
!   pair: The new atom to atom pair association (to be used
!         in the next iteration of the algorithm).
!   ngaps: Number of gaps of the correspondence
!

subroutine nonbscore(na, nb, prota, protb, dtri2, gap,&
                     disord, indisord, it, pair,& 
                     score, bije, nbij, ngaps)

  use sizes
  implicit none
  double precision :: prota(maxatom,3), protb(maxatom,3),&
                      disord(maxatom-1, maxatom),&
                      distpair(maxatom), score, dtri2, dmin, gap
  integer :: na, nb, nbij, bije(maxatom,2), pair(maxatom),&
             indisord(maxatom-1, maxatom),&
             i, nearest, it, ngaps

  ! If this is the first iteration, the first guess is not very good 

  if(it.eq.0) then

    ! First guess for first atom

    nearest = 1
    call get_nearest(protb,nb,indisord,disord,&
                     prota(1,1), prota(1,2), prota(1,3),&
                     nearest,dmin)
    pair(1) = nearest
    distpair(1) = dmin

    ! Guesses for other atoms (use the result for first atom to
    ! improve guess

    do i = 2, na
      nearest = pair(i-1)
      call get_nearest(protb,nb,indisord,disord,&
                       prota(i,1), prota(i,2), prota(i,3),&
                       nearest,dmin)
      pair(i) = nearest
      distpair(i) = dmin
    end do

    ! For the second iteration on, use the previous assignment as the 
    ! guess for the next one

  else
    do i = 1, na
      nearest = pair(i)
      call get_nearest(protb,nb,indisord,disord,&
                       prota(i,1), prota(i,2), prota(i,3),&
                       nearest,dmin)
      pair(i) = nearest
      distpair(i) = dmin
    end do
  end if

  ! Include in the correspondence only pairs for which the distance
  ! is smaller than dtri, and compute TRIANGULAR SCORE

  nbij = 0
  score = 0.d0
  do i = 1, na 
    if( distpair(i) .lt. dtri2 ) then
      nbij = nbij + 1
      bije(nbij,1) = i
      bije(nbij,2) = pair(i)
      score = score + 20.d0*( 1.d0 - distpair(i) / dtri2 )
    end if
  end do

  ! Penalization for gaps (not reasonable to use for this type of alignment)

  if( gap .gt. 1.d-10 ) then
    ngaps = 0
    do i = 1, nbij - 1
      if ( bije(i+1,1) .ne. bije(i,1) + 1 .or.&
           bije(i+1,2) .ne. bije(i,2) + 1 ) then
        ngaps = ngaps + 1
      end if
    end do
    score = score - gap*ngaps
  end if

end subroutine nonbscore

!
! Subroutine get_nearest: gets the nearest atom of protein B to point 
! X using the fast algorithm that uses the ordered internal distances
! of protein B
!
! On input:
!   protb: The coordinates of proteinB
!   nb: The number of atoms of proteinB
!   indisord: Matrix containing the indexes of of the atoms of 
!             proteinB such that, for each atom, the other 
!             atoms are ordered from closer to farther. 
!   disord: Matrix containing the internal distances of proteinB
!   x, y, z: coordinates of current atom (of proteinA)
!   nearest: First guess of what atom is the closest one
!
! On output:
!   nearest: The index of the closest atom.
!   dmin: The (squared) distance of the nearest atom to x.
!
      
subroutine get_nearest(protb,nb,indisord,disord,&
                       x, y, z,&
                       nearest, dmin)

  use sizes
  implicit none
  double precision :: protb(maxatom,3), disord(maxatom-1,maxatom),&
                      x, y, z, dmin, d, dbase
  integer :: nearest, i, nb, indisord(maxatom-1,maxatom), ibase

  d = ( x - protb(nearest,1) )**2 + &
      ( y - protb(nearest,2) )**2 + &
      ( z - protb(nearest,3) )**2
  dmin = d
  ibase = nearest
  dbase = dmin

  do i = 1, nb - 1

    if( disord(i,ibase) .ge. 4.d0*dbase ) return

    d = ( x - protb(indisord(i,ibase),1) )**2 + &
        ( y - protb(indisord(i,ibase),2) )**2 + &
        ( z - protb(indisord(i,ibase),3) )**2

    if( d .lt. dmin ) then
      dmin = d
      nearest = indisord(i,ibase)
    end if

  end do

end subroutine get_nearest

!
! Subroutine getrmsd: Computes the rmsd given the bijection
!

subroutine getrmsd(prota,protb,bije,nbij,rmsd)

  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2)
  double precision :: rmsd, dist, prota(maxatom,3), protb(maxatom,3)

  rmsd = 0.d0
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    rmsd = rmsd + dist
  end do
  rmsd = dsqrt(rmsd/dfloat(nbij))

end subroutine getrmsd

!
! Subroutine that computes the GDT scores
!

subroutine computegdt(na,nb,prota,protb,bije,nbij,gdt_threshold,gdt_tm,gdt_ha)
 
  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2), na, nb
  double precision :: gdt_tm, gdt_ha, gdt_threshold, gdt_ha_threshold
  double precision :: prota(maxatom,3), protb(maxatom,3), dist

  gdt_ha_threshold = gdt_threshold / 2.d0
  gdt_tm = 0.d0
  gdt_ha = 0.d0
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    dist = dsqrt(dist)

    ! GDT score (with threshold = 4.d0 Angs by default)
    if ( dist < (gdt_threshold/4.d0) ) gdt_tm = gdt_tm + 1.d0 
    if ( dist < (gdt_threshold/2.d0) ) gdt_tm = gdt_tm + 1.d0
    if ( dist < gdt_threshold ) gdt_tm = gdt_tm + 1.d0
    if ( dist < (2.d0*gdt_threshold) ) gdt_tm = gdt_tm + 1.d0

    ! GDT-"High Accuracy" score (with threshold = 2.d0 Angs by default)
    if ( dist < (gdt_ha_threshold/4.d0) ) gdt_ha = gdt_ha + 1.d0 
    if ( dist < (gdt_ha_threshold/2.d0) ) gdt_ha = gdt_ha + 1.d0
    if ( dist < gdt_ha_threshold ) gdt_ha = gdt_ha + 1.d0
    if ( dist < (2.d0*gdt_ha_threshold) ) gdt_ha = gdt_ha + 1.d0
  end do
  gdt_tm = 100.d0 * gdt_tm / (4.d0*min(na,nb))
  gdt_ha = 100.d0 * gdt_ha / (4.d0*min(na,nb))

end subroutine computegdt

!
! Subroutine getrmsd2: Computes the rmsd given the bijection,
!                      only for atoms which are closer than
!                      some tolerance deffined by the dtri parameter
!

subroutine getrmsd2(prota,protb,bije,nbij,rmsd,nbij_dtri,dtri)

  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2), nbij_dtri
  double precision :: rmsd, dist, prota(maxatom,3), protb(maxatom,3),&
                      dtri, dtri2

  dtri2 = dtri*dtri
  rmsd = 0.d0
  nbij_dtri = nbij
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    if(dist.le.dtri2) then
      rmsd = rmsd + dist
    else
      nbij_dtri = nbij_dtri - 1
    end if         
  end do
  if(nbij_dtri.gt.0) then
    rmsd = dsqrt(rmsd/dfloat(nbij_dtri))
  else
    rmsd = 0.d0
  end if

end subroutine getrmsd2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
! Subroutine prodin: Performs Dynamic Programming in a matrix of    !
!                    individual scores and returns the bijection    !
!                    that globally maximizes the score and the      !
!                    corresponding score                            !
!                                                                   !
! On input:                                                         !
!          na: number of rows of the individual score matrix        !
!          nb: number of columns of the individual score matrix     !
!          scorin: the individual score matrix                      !
!          gap: Penalization for gaps                               !
!                                                                   !
! On return:                                                        !
!          bije: the bijection that maximizes the score             !
!          nbij: number of correspondences in the bijection         !
!          ngaps: number of gaps of the bijection                   !
!                                                                   !
! This subroutine takes most of the time of the calculations.       !
! If you have a better way to perform this dynamic programming      !
! it would have an important impact on the total alignment time.    !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prodin(na,nb,scorin,gap,ngaps,bije,nbij,bijscore,scomax)

  use sizes
  implicit none
  integer :: i, j, na, nb, bije(maxatom,2), nbij, ngaps, nl, nc, maxj,&
             maxi(maxatom), seg(maxatom,maxatom,2), imax, jmax 
  double precision scorin(maxatom,maxatom), aux(maxatom,maxatom),&
                   gap, gaf(maxatom,maxatom), bijscore(maxatom),&
                   scomax, gamax

  ! As currently implemented, the cost of this procedure depends on
  ! which protein is the largest one. Therefore, we always set the
  ! the largest protein to be the first

  if(na.ge.nb) then
    do j = 1, nb
      do i = 1, na
        aux(i,j) = scorin(i,j)
      end do
    end do
    nc = na
    nl = nb
  else 
    do j = 1, nb
      do i = 1, na
        aux(j,i) = scorin(i,j)
      end do
    end do
    nc = nb
    nl = na
  end if

  ! Here begin the dynamic programming procedure

  do i = 1, nc
    gaf(i,nl) = aux(i,nl)
    maxi(i) = nl
  end do

  ! Dynamic programming in O(n^2)

  do j = nl - 1, 1, -1
    gaf(nc,j) = aux(nc,j)
    maxj = nc
    do i = nc - 1, 1, -1
      gamax = gaf(i+1,j+1)
      seg(i,j,1) = i + 1
      seg(i,j,2) = j + 1
      if(gaf(i+1,maxi(i+1))-gap.gt.gamax) then
        gamax = gaf(i+1, maxi(i+1)) - gap
        seg(i,j,1) = i+1
        seg(i,j,2) = maxi(i+1)
      end if
      if(gaf(maxj, j+1)-gap.gt.gamax) then
        gamax = gaf(maxj,j+1) - gap
        seg(i,j,1) = maxj
        seg(i,j,2) = j + 1
      end if
      if(gaf(i,j+1).ge.gaf(maxj,j+1)) maxj = i
      gaf(i,j) = aux(i,j) + gamax
    end do
    do i = 1, nc
      if(gaf(i,j).ge.gaf(i, maxi(i))) maxi(i)=j
    end do
  end do
      
  ! Now will compute the score to the end

  scomax = 0.d0
  do j = 1, nl
    if(gaf(1,j).ge.scomax) then
      imax = 1
      jmax = j
      scomax = gaf(1,j)
    end if
  end do
  do i = 1, nc
    if(gaf(i,1).ge.scomax) then
      imax = i
      jmax = 1
      scomax = gaf(i,1)
    endif
  end do

  ! Save maximum score obtained and the corredponding bijection      

  bije(1,1) = imax
  bije(1,2) = jmax
  nbij = 1
  i = 2
  do i = 2, nc
    if(bije(i-1,1).eq.nc.or.bije(i-1,2).eq.nl) exit
    bije(i,1) = seg(bije(i-1,1),bije(i-1,2),1)
    bije(i,2) = seg(bije(i-1,1),bije(i-1,2),2)
    nbij = i
  end do

  ngaps = 0
  do i = 1, nbij-1
    if(bije(i+1,1).ne.bije(i,1)+1.or.&
       bije(i+1,2).ne.bije(i,2)+1) then
      ngaps = ngaps+1
    endif
  end do

  do i = 1, nbij
    bijscore(i) = aux(bije(i,1),bije(i,2))
  end do

  ! If the largest protein was protein B, restore

  if(na.lt.nb) then
    do i = 1, nbij
      j = bije(i,1)
      bije(i,1) = bije(i,2)
      bije(i,2) = j
    end do
  end if

end subroutine prodin


!
! Subroutines used for heuristic initial approximations
!
 
!
! Subroutine pseudoprot: Computes pseudoprotein A for initial point
!                        generation
!

subroutine pseudoprot(prota,pseudoa,na)

  use sizes
  implicit none
  integer :: i, na
  double precision :: prota(maxatom,3), pseudoa(maxatom,3), dist

  if(na.le.5) then
    return
  end if

  do i = 1, na - 4
    pseudoa(i,1) = dist(prota,i,i+2)
    pseudoa(i,2) = dist(prota,i,i+3)
    pseudoa(i,3) = dist(prota,i,i+4)
  end do
  pseudoa(na-3,1) = dist(prota,na-3,na-1) 
  pseudoa(na-3,2) = dist(prota,na-3,na)
  pseudoa(na-3,3) = dist(prota,na-2,na)

end subroutine pseudoprot

! Function that computes the distance between two atoms

double precision function dist(prota,i,j)

  use sizes
  implicit none
  integer :: i, j
  double precision :: prota(maxatom,3)

  dist = dsqrt((prota(i,1)-prota(j,1))**2 + &
               (prota(i,2)-prota(j,2))**2 + &
               (prota(i,3)-prota(j,3))**2)
  
end function dist
 
!
! Next, the subroutine for file input and output
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
! SUBROUTINES FOR READING THE PARAMETERS FROM THE COMMAND LINE  !
!                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Subroutine that reads the parameters from the command line
! 

subroutine getpars(method,gap,maxit,dtri,gdt_threshold,&
                   protea,proteb,chaina,chainb,iprint,&
                   all,seqoff,seqfix,&
                   pdblist,output,pdbout,useini,beta1,beta2,&
                   ocup1,ocup2)
 
  use ioformat
  use sizes
  implicit none

  integer :: method, maxit, iprint, narg, length, i, ival, iargc, ioerr
  double precision :: gap, dval, dtri, gdt_threshold
  character(len=1) :: chaina, chainb
  character(len=200) :: protea, proteb, pdblist, keyword, value, pdbout
  logical :: output, all, seqoff, useini, beta1, beta2, ocup1, ocup2
  logical :: seqfix

  ! Reading the command line specifications

  narg = iargc() 
  i = 1
  gap = 1.d30
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
    else if(keyword(1:length(keyword)).eq.'-pdblist') then
      call getarg(i+1,value)
      pdblist = value
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
    else if(keyword(1:length(keyword)).eq.'-all') then
      all = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-seqoff') then
      seqoff = .true.
      i = i - 1
    else if(keyword(1:length(keyword)).eq.'-seqfix') then
      seqfix = .true.
      i = i - 1
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

!
! Subroutine that reads protein coordinates both in the pdb format or in
! simple clean coordinates file used for tests or other purposes
!

subroutine readfile(protea,prota,chaina,beta1,ocup1,na,resa,numa,all,error)

  use sizes
  implicit none
  integer :: na, length, ic, numa(maxatom),ioerr
  real :: beta, occupancy
  double precision :: prota(maxatom,3)
  logical :: all, error, beta1, ocup1
  character(len=1) :: chaina, resa(maxatom), letter
  character(len=3) :: resid
  character(len=200) :: protea, record

  error = .false. 

  !
  ! Trying to read coordinates in the simple xyz (only coordinates) format
  !

  open(10,file=protea(1:length(protea)),status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file:',protea(ic(protea):length(protea))
    error=.true.
    return
  end if

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

          numa(na) = na
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

!
! Function that gets the integer value from a variable string
!

integer function ival(string)

  implicit none
  integer :: ioerr
  character(len=200) :: string

  read(string,*,iostat=ioerr) ival
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read integer value from some'
    write(*,*) '        of the parameters. Some parameter with'
    write(*,*) '        a expected integer value was not set '
    write(*,*) '        using -keyword [integer]'
    stop
  end if

end function ival

!
! Function that gets the double precision value from a variable string
!

double precision function dval(string)

  implicit none
  integer :: ioerr
  character(len=200) :: string

  read(string,*,iostat=ioerr) dval
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read real value from some'
    write(*,*) '        of the parameters. Some parameter with'
    write(*,*) '        a expected real value was not set '
    write(*,*) '        using -keyword [real]'
    stop
  end if

end function dval

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

!
! Subroutine readlist: Reads names from the pdblist for running
! modes 1 and 2
!

subroutine readlist(pdblist,pdbfiles,nfiles)

  use sizes
  implicit none
  integer :: nfiles, length, ioerr
  character(len=200) :: pdblist, pdbfiles(maxfiles), record

  open(10,file=pdblist,status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', pdblist(1:length(pdblist))
    stop
  end if

  nfiles = 0
  do while(.true.)
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(record.gt.' ') then
      nfiles = nfiles + 1
      pdbfiles(nfiles) = record(1:length(record))
    end if
  end do
  close(10)

  if(nfiles.gt.maxfiles) then
    write(*,*) ' Number of files in list greater than MAXFILES. '
    write(*,*) ' Increase the maxfiles parameter. '
    stop
  end if

end subroutine readlist

!
! Subroutine writepdb: Reads all atoms of protein (not only CAs),
!                      translates and rotates them according to
!                      the best alignment obtained and prints the
!                      output pdb file with the aligned atoms. 
!

subroutine writepdb(pdbout,protea,prota,chaina,beta1,ocup1,&
                    na,resa,numa,proteb,all)

  use sizes
  implicit none
      
  integer :: i, j, na, nca, iq, length, ic, numa(maxatom), ioerr
  double precision :: prota(maxatom,3), aux(maxatom,3), cmo(3),&
                      cma(3), xm(maxatom), ym(maxatom), zm(maxatom),&
                      xp(maxatom), yp(maxatom), zp(maxatom), q(4,4),&
                      a(4,4), qmin, u(3,3), xn(3)
  character(len=200) pdbout, protea, record, proteb
  character(len=1) chaina, resa(maxatom)
  logical :: all, mark, error, beta1, ocup1

  ! Reading the original CA coordinates

  call readfile(protea,aux,chaina,beta1,ocup1,nca,resa,numa,all,error)
  if(nca.ne.na) then
    write(*,*) ' ERROR reading the file of the first'
    write(*,*) ' protein when going to write the output file.'
    write(*,*) ' Wrong number of atoms.'
    write(*,*) ' Output file not written.'
    return
  end if 

! If the coordinates were given in xyz format we can have an error here
 
  if(nca.lt.3) then
    write(*,*) ' ERROR: The output option is only available for files in PDB format. '
    return
  end if
 
! Translating the original and aligned proteins to the origin

  do i = 1, 3
    cmo(i) = 0.d0
    cma(i) = 0.d0
  end do
  do i = 1, na
    cmo(1) = cmo(1) + aux(i,1)
    cmo(2) = cmo(2) + aux(i,2)
    cmo(3) = cmo(3) + aux(i,3)
    cma(1) = cma(1) + prota(i,1)
    cma(2) = cma(2) + prota(i,2)
    cma(3) = cma(3) + prota(i,3)
  end do
  do i = 1, 3
    cmo(i) = cmo(i) / dfloat(na)
    cma(i) = cma(i) / dfloat(na)
  end do
  do i = 1, na
    aux(i,1) = aux(i,1) - cmo(1)
    aux(i,2) = aux(i,2) - cmo(2)
    aux(i,3) = aux(i,3) - cmo(3)
    prota(i,1) = prota(i,1) - cma(1)
    prota(i,2) = prota(i,2) - cma(2)
    prota(i,3) = prota(i,3) - cma(3)
  end do

  ! Obtaining the best rotation that alignes them (this is almost a whole
  ! procrustes iteration of subroutine procrustes.

  do i = 1, na
    xm(i) = prota(i,1) - aux(i,1) 
    ym(i) = prota(i,2) - aux(i,2) 
    zm(i) = prota(i,3) - aux(i,3) 
    xp(i) = prota(i,1) + aux(i,1) 
    yp(i) = prota(i,2) + aux(i,2) 
    zp(i) = prota(i,3) + aux(i,3)
  end do 

  do j = 1, 4
    do i = 1, 4
      q(i,j) = 0.
    end do
  end do

  do i = 1, na
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)
     
  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix

  call jacobi(q,a,4)

  ! Choosing the quaternion that corresponds to the minimum

  qmin = 1.e30
  do i = 1, 4
    if(q(i,i).lt.qmin) then
      iq = i
      qmin = q(i,i)
    end if
  end do
  
! Computing the rotation matrix

  u(1,1) = a(1,iq)**2 + a(2,iq)**2 - a(3,iq)**2 - a(4,iq)**2
  u(1,2) = 2. * ( a(2,iq)*a(3,iq) + a(1,iq)*a(4,iq) )
  u(1,3) = 2. * ( a(2,iq)*a(4,iq) - a(1,iq)*a(3,iq) )  
  u(2,1) = 2. * ( a(2,iq)*a(3,iq) - a(1,iq)*a(4,iq) )  
  u(2,2) = a(1,iq)**2 + a(3,iq)**2 - a(2,iq)**2 - a(4,iq)**2 
  u(2,3) = 2. * ( a(3,iq)*a(4,iq) + a(1,iq)*a(2,iq) )  
  u(3,1) = 2. * ( a(2,iq)*a(4,iq) + a(1,iq)*a(3,iq) )  
  u(3,2) = 2. * ( a(3,iq)*a(4,iq) - a(1,iq)*a(2,iq) )  
  u(3,3) = a(1,iq)**2 + a(4,iq)**2 - a(2,iq)**2 - a(3,iq)**2 
 
  ! Now we need to read the pdb file again and apply the transformations
  ! resulted in the best alignment for all atoms of the protein

  open(10,file=protea(1:length(protea)),status='old')
  open(20,file=pdbout(1:length(pdbout)))
  write(20,"('REMARK Protein ',a,' aligned to',/,&
            &'REMARK protein ',a,/,&
            &'REMARK with LovoAlign: ',&
            &'http://www.ime.unicamp.br/~martinez/lovoalign')")&
        protea(ic(protea):length(protea)), proteb(ic(proteb):length(proteb))

  mark = .false.
  do while(.true.)
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(length(record).lt.82) then
      do i = length(record)+1, 82
        record(i:i) = ' '
      end do
    end if
    if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
      read(record(31:38),*,iostat=ioerr) aux(1,1)
      if ( ioerr /= 0 ) error = .true.
      read(record(39:46),*,iostat=ioerr) aux(1,2)
      if ( ioerr /= 0 ) error = .true.
      read(record(47:54),*,iostat=ioerr) aux(1,3)
      if ( ioerr /= 0 ) error = .true.
      if ( error ) then
        write(*,*) ' ERROR: Failed reading coordinates in file ',&
                   protea(ic(protea):length(protea)),&
                   ': output file may be not complete.' 
        exit
      end if
      aux(1,1) = aux(1,1) - cmo(1)
      aux(1,2) = aux(1,2) - cmo(2)
      aux(1,3) = aux(1,3) - cmo(3)
      xn(1) = u(1,1)*aux(1,1)+u(1,2)*aux(1,2)+u(1,3)*aux(1,3)
      xn(2) = u(2,1)*aux(1,1)+u(2,2)*aux(1,2)+u(2,3)*aux(1,3)
      xn(3) = u(3,1)*aux(1,1)+u(3,2)*aux(1,2)+u(3,3)*aux(1,3)
      xn(1) = xn(1) + cma(1)
      xn(2) = xn(2) + cma(2)
      xn(3) = xn(3) + cma(3)
      write(20,"(a30,3(f8.3),a28)") record(1:30), xn(1), xn(2), xn(3), record(55:82)
    else
      write(20,"(t1,a82)") record(1:82)
    end if
  end do
  close(10)
  close(20)

end subroutine writepdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                                      !
! Subroutine orprot: Computes the ordered square distance matrices for !
! the fast evaluation of correspondences                               !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine orprot(nprot,proteina,disord,indisord)

  use sizes
  implicit none
  
  ! disprot(nb,nb) contains the internal square distances of protein B
  ! disord(nb,nb-1) contains, for each atom of B, the distances to the other
  ! atoms, in increasing order indisord(nb,nb-1) contains the indices of the
  ! atoms corresponding to disord.
 
  double precision :: proteina(maxatom, 3), disprot(maxatom, maxatom),&
                      disord(maxatom-1, maxatom), z
  real :: aux(maxatom)
  integer :: i, j, k, nprot
  integer :: indisord(maxatom-1, maxatom), lflash(maxatom), mflash, indflash(maxatom) 

  !  Compute disprot, squared internal distances

  mflash = 1 + int(float(nprot)/10.)
  do i = 1, nprot
    disprot(i, i) = 0.d0
  end do

  do j = 1, nprot-1
    do i = j+1, nprot
      z = 0.d0
      do k = 1, 3
        z = z + (proteina(i, k) - proteina(j, k))**2
      end do
      disprot(i, j) = z
      disprot(j, i) = disprot(i, j)
    end do
  end do

  ! Sort the distances to atom i of the protein

  do j = 1, nprot
    do i = 1, nprot
      aux(i) = real(disprot(i,j))
    end do
    call flash1(aux, nprot, lflash, mflash, indflash)
    do i = 1, nprot-1
      disord(i,j) = aux(i+1)
      indisord(i,j) = indflash(i+1)
    end do
  end do

end subroutine orprot

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

!        
! Subroutine writebije: Writes the bijection obtained at the
!                       solution

subroutine writebije(na,nb,resa,resb,numa,numb,bije,nbij)

  use sizes
  implicit none
      
  integer :: na, nb, bije(maxatom,2), nbij, nlines,&
             numa(maxatom), numb(maxatom), i, j, iamin, ibmin,&
             ia, ib, ina, inb, ic
  character(len=1) :: resa(maxatom), resb(maxatom)
  character(len=50) :: line1, line2

  write(*,"(/,'  ',25('-'),' SEQUENCE ALIGNMENT ',26('-'))")

  nlines = (na + nb - nbij - 1) / 50 + 1

  iamin = 1
  ibmin = 1
  ia = 1
  ib = 1
  ina = 1
  inb = 1
  do i = 1, nlines
    do j = 1, 50
      line1(j:j) = ' '
      line2(j:j) = ' '
    end do
    ic = 1
    do while(ic.le.50)
      if(ia.lt.bije(ina,1)) then
        line1(ic:ic) = resa(ia)
        line2(ic:ic) = '-'
        ia = ia + 1
      else if(ib.lt.bije(inb,2)) then
        line1(ic:ic) = '-'
        line2(ic:ic) = resb(ib)
        ib = ib + 1
      else
        if(ia.le.na) line1(ic:ic) = resa(ia)
        if(ib.le.nb) line2(ic:ic) = resb(ib)
        if(ia.gt.na.and.ib.le.nb) line1(ic:ic) = '-'
        if(ib.gt.nb.and.ia.le.na) line2(ic:ic) = '-'
        ia = ia + 1
        ib = ib + 1
        if(ina.lt.nbij) ina = ina + 1
        if(inb.lt.nbij) inb = inb + 1
      end if
      ic = ic + 1
    end do 
    write(*,"('            .         .         .         .         .')")
    write(*,"(tr6,i5,tr1,a50,tr1,i5)") numa(min0(iamin,na)),&
                                       line1,&
                                       numa(max0(1,min0(ia-1,na)))
    write(*,"(tr6,i5,tr1,a50,tr1,i5)") numb(min0(ibmin,nb)),&
                                       line2,&
                                       numb(max0(1,min0(ib-1,nb)))
    if(ia.gt.iamin) iamin = ia 
    if(ib.gt.ibmin) ibmin = ib 
  end do
  write(*,*)

end subroutine writebije

!
! One letter amino acid codes for printing sequence alignment
!

character function letter(resid)

  implicit none
  character*3 resid

  if(resid.eq.'ALA') then
    letter = 'A'
  else if(resid.eq.'ARG') then
    letter = 'R'
  else if(resid.eq.'ASN') then
    letter = 'N'
  else if(resid.eq.'ASP') then
    letter = 'D'
  else if(resid.eq.'ASX') then
    letter = 'B'
  else if(resid.eq.'CYS') then
    letter = 'C'
  else if(resid.eq.'GLU') then
    letter = 'E'
  else if(resid.eq.'GLN') then
    letter = 'Q'
  else if(resid.eq.'GLX') then
    letter = 'Z'
  else if(resid.eq.'GLY') then
    letter = 'G'
  else if(resid.eq.'HIS') then
    letter = 'H'
  else if(resid.eq.'HSD') then
    letter = 'H'
  else if(resid.eq.'HSP') then
    letter = 'H'
  else if(resid.eq.'HSE') then
    letter = 'H'
  else if(resid.eq.'ILE') then
    letter = 'I'
  else if(resid.eq.'LEU') then
    letter = 'L'
  else if(resid.eq.'LYS') then
    letter = 'K'
  else if(resid.eq.'MET') then
    letter = 'M'
  else if(resid.eq.'PHE') then
    letter = 'F'
  else if(resid.eq.'PRO') then
    letter = 'P'
  else if(resid.eq.'SER') then
    letter = 'S'
  else if(resid.eq.'THR') then
    letter = 'T'
  else if(resid.eq.'TRP') then
    letter = 'W'
  else if(resid.eq.'TYR') then
    letter = 'Y'
  else if(resid.eq.'VAL') then
    letter = 'V'
  else if(resid.eq.'XXX') then
    letter = 'X'
  else
    letter = '?'
  end if

end function letter

!
! Print problem specifications
!

subroutine printdata(protea,proteb,na,nb,chaina,&
                     chainb,method,gap,maxit,dtri,gdt_threshold,&
                     useini)

  use ioformat
  implicit none
  integer :: na, nb, maxit, method, length, ic
  double precision :: gap, dtri, gdt_threshold
  character(len=1) :: chaina, chainb
  character(len=200) :: protea, proteb
  logical :: useini
       
  write(*,*) ' Problem specifications: '
  write(*,dash_line) 
  write(*,*) ' Protein A: ', protea(ic(protea):length(protea))
  write(*,*) ' Protein B: ', proteb(ic(proteb):length(proteb))
  write(*,*) ' Number of atoms: A:', na, ' B:', nb
  if(chaina.ne.'#') write(*,*) ' Protein A chain: ', chaina
  if(chainb.ne.'#') write(*,*) ' Protein B chain: ', chainb
  if(method.eq.1) write(*,*) ' Will maximize the STRUCTAL score'
  if(method.eq.2) write(*,*) ' Will maximize the TM-SCORE '
  if(method.eq.3) write(*,*) ' Will maximize the TRIANGULAR score '
  if(method.eq.4) write(*,*) ' Will maximize the NON-BIJECTIVE TRIANGULAR score '
  write(*,*) ' Penalization for gaps: ', gap
  write(*,*) ' Maximum number of iterations: ', maxit
  if(useini) write(*,*) ' Using internal-distance initial point.'
  if(method.eq.3) then
    write(*,*) ' Triangular score with cutoff: ', dtri
  end if
  write(*,*) ' GDT Threshold: ', gdt_threshold
  write(*,dash_line) 

end subroutine printdata
 
!
! Subroutine help: Prints simple usage instructions
!

subroutine help

  use ioformat
  implicit none
  call title()
  write(*,*) 
  write(*,*) ' How to align two proteins: '
  write(*,*) ' ./lovoalign -p1 prot1.pdb -p2 prot2.pdb -o p1aligned.pdb'
  write(*,*)
  write(*,*) ' How to align a single protein to a database: '
  write(*,*) ' ./lovoalign -p1 prot1.pdb -pdblist files.dat '
  write(*,*) 
  write(*,*) ' Performing an all-on-all database comparison:'
  write(*,*) ' ./lovoalign -pdblist files.dat '
  write(*,*)
  write(*,*) ' Some options: '
  write(*,*) ' -m 1   Maximize STRUCTAL score '
  write(*,*) ' -m 2   Maximize TM-SCORE '
  write(*,*) ' -m 3   Maximize Triangular score '
  write(*,*) ' -c1 A  Consider only chain A of protein 1 '
  write(*,*) ' -c2 A  Consider only chian A of protein 2 '
  write(*,*) ' -beta1 Consider atoms with beta > 0 in protein 1'
  write(*,*) ' -beta2 Consider atoms with beta > 0 in protein 2'
  write(*,*) ' -ocup1 Consider atoms with occupancy > 0 in protein 1'
  write(*,*) ' -ocup2 Consider atoms with occupancy > 0 in protein 2'
  write(*,*) ' -g [real] Penalization for gaps '
  write(*,*) ' -dtri [real] Atoms farther than this will not be'
  write(*,*) '        considered. Distance for Triangular score.'
  write(*,*)
  write(*,*) ' Other options and instructions can be found at  '
  write(*,*) ' the initial source code comments and at: '
  write(*,*) ' http://www.ime.unicamp.br/~martinez/lovoalign'
  write(*,*) 
  write(*,*) ' Authors: L. Martinez, R. Andreani, J. M. Martinez. '
  write(*,*) ' University of Campinas (UNICAMP) - Brazil '
  write(*,dash_line)
  stop

end subroutine help
