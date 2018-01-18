!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
! SUBROUTINE PROTALL: This is the main subroutine of the program. The !
! actual subroutine that performs the alignment between two proteins. !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
subroutine protall(prota,protb,na,nb,disord,indisord,resa,resb,numa,numb) 

  use sizes
  use inputpars
  use ioformat
  use bijetype
  implicit none

  double precision :: dzero, prota(maxatom,3), protb(maxatom,3),&
                      bijscore(maxatom), gnor,&
                      score, dzero2, tol, scale,&
                      prevscore, rmsd, rmsd2, dtri2,&
                      disord(maxatom-1,maxatom), &
                      gdt_ts, gdt_ha, scorebest, prota_best(maxatom,3)
  real :: etime, tarray(2), time1
  integer :: na, nb, i, j, ii, jj, jlast, &
             bije(maxatom,2), bijebest(maxatom,2),&
             ngaps, nbij, nbijbest, nef, nbij_dtri,&
             length, ic, numa(maxatom),&
             numb(maxatom), it, iglobal, itrial, &
             indisord(maxatom-1,maxatom), pair(maxatom)
  character(len=1) :: resa(maxatom), resb(maxatom)
  external :: structal, tmscore

  ! Time used in this alignment is computed from here

  time1 = etime(tarray)

  ! This is the relative precision for score convergence

  tol = 1.d-2

  ! Define the fixed sequence alignment, if it is the case

  ! Fixed alignment: 1:1-2:2-etc.
  if ( seqtype > 0 ) then
    if ( seqtype == 1 ) then
      fixnbij = min(na,nb)
      do i = 1, fixnbij
        fixbije(i,1) = i 
        fixbije(i,2) = i 
      end do
    end if
    ! Fixed alignment based on residue numbers
    if ( seqtype == 2 ) then
      fixnbij = 0 
      jlast = 1
      aatoms : do i = 1, na
        do j = jlast, nb
          if ( numa(i) == numb(j) ) then
            fixnbij = fixnbij + 1
            fixbije(fixnbij,1) = i
            fixbije(fixnbij,2) = j
            jlast = j
            cycle aatoms
          end if
        end do
      end do aatoms
    end if
    ! Fasta alignment 
    if ( seqtype == 3 ) then
      fixnbij = 0
      i = 1
      j = 1
      ii = 1
      jj = 1
      do while( ii <= na .and. jj <= nb ) 
        if ( fasta(1)%seq(i:i) /= "-" .and. &
             fasta(2)%seq(j:j) /= "-") then
          fixnbij = fixnbij + 1
          fixbije(fixnbij,1) = ii
          fixbije(fixnbij,2) = jj
        end if
        if ( fasta(1)%seq(i:i) /= "-" ) ii = ii + 1
        if ( fasta(2)%seq(j:j) /= "-" ) jj = jj + 1
        i = i + 1
        j = j + 1
      end do
    end if
    if ( fixnbij < 4 ) then
      write(*,*) ' ERROR: Less than 4 atoms are associated in the sequence alignment '
      stop
    end if
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                !
  ! Method 1: Maximize the STRUCTAL score of Gerstein and Levitt   !
  !                                                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  if(method.eq.1) then

    ! Normalization of STRUCTAL score
  
    dzero = 2.24d0
    dzero2 = dzero*dzero
    scale = 20.d0

    scorebest = -1.d0
    iglobal = 0
    itrial = 0
    do while( iglobal < nglobal .and. itrial <= maxtrial  ) 
      itrial = itrial + 1

      ! Create random initial point of itrial > 1 

      if ( itrial > 1 ) call randomini(na,nb,prota,protb,nbij,bije)
 
      ! Writes the titles for regular iteration printing

      if(iprint.eq.2) write(*,title_format)

      ! Number iterations, functional evaluations and DP calls

      nef = 0 
      it = 0
      prevscore = -1.d0

      ! Compute the DP bijection and the score at initial point          

      call structal(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                    bijscore,ngaps,score)

      nef = nef + 1
      if(iprint.eq.2) write(*,data_format) it, score, 0.d0, nbij, ngaps,nef
            
      ! Here begin the iteration loop

      do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

        it = it + 1
        prevscore = score

        ! Perform a newton step to maximize the score given the bijection

        call newton(structal,na,nb,prota,protb,score,bije,bijscore,&
                    dzero2,scale,nbij,gap,ngaps,nef,gnor)

        ! Output regular iteration data
            
        if(iprint.eq.2) write(*,data_format) it, score, gnor, nbij, ngaps,nef

      end do
      if ( score > scorebest ) then
        if ( abs(score-scorebest) > tol*score ) then
          iglobal = 1
        end if
        nbijbest = nbij
        scorebest = score
        do i = 1, nbij
          bijebest(i,1) = bije(i,1)
          bijebest(i,2) = bije(i,2)
        end do
        do i = 1, na
          prota_best(i,1) = prota(i,1)
          prota_best(i,2) = prota(i,2)
          prota_best(i,3) = prota(i,3)
        end do
        if(iprint.eq.1) write(*,trial_format) itrial, score, nbij, ngaps, iglobal
      else
        if ( scorebest-score <= tol*score ) then
          iglobal = iglobal + 1
          if ( iglobal == nglobal .and. iprint == 1 ) then
            write(*,"(a,i5,a,i5)") '  Repeated best solution found ', nglobal,' times at trial ', itrial
          end if
        end if
      end if
      if(iprint.eq.2) write(*,trial_format) itrial, score, nbij, ngaps, iglobal

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

    ! Normalization of the TM-SCORE score

    dzero = 1.24d0 * (nb-15.d0)**(1.d0/3.d0) - 1.8d0
    dzero2 = dzero*dzero
    scale = 1.d0 / dfloat(nb)

    scorebest = -1.d0
    iglobal = 0
    itrial = 0
    do while( iglobal < nglobal .and. itrial <= maxtrial ) 
      itrial = itrial + 1

      ! Writes the titles for regular iteration printing

      if(iprint.eq.2) write(*,title_format)

      ! Create random initial point of itrial > 1 

      if ( itrial > 1 ) call randomini(na,nb,prota,protb,nbij,bije)

      ! Number iterations, functional evaluations and DP calls

      nef = 0 
      it = 0
      prevscore = -1.d0
 
      ! Compute the DP bijection and the score at initial point          

      call tmscore(prota,protb,na,nb,dzero2,gap,bije,nbij,&
                   bijscore,ngaps,score)
      nef = nef + 1
      if(iprint.eq.2) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
            
      ! Here begin the iteration loop

      do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

        it = it + 1
        prevscore = score

        ! Perform a newton step to maximize the score given the bijection

        call newton(tmscore,na,nb,prota,protb,score,bije,bijscore,&
                    dzero2,scale,nbij,gap,ngaps,nef,gnor)

        ! Output regular iteration data
            
        if(iprint.eq.2) write(*,data_format) it, score, gnor, nbij, ngaps, nef

      end do
      if ( score > scorebest ) then
        iglobal = 1
        nbijbest = nbij
        scorebest = score
        do i = 1, nbij
          bijebest(i,1) = bije(i,1)
          bijebest(i,2) = bije(i,2)
        end do
        do i = 1, na
          prota_best(i,1) = prota(i,1)
          prota_best(i,2) = prota(i,2)
          prota_best(i,3) = prota(i,3)
        end do
        if(iprint.eq.1) write(*,trial_format) itrial, score, nbij, ngaps, iglobal
      else
        if ( scorebest-score <= tol*score ) then
          iglobal = iglobal + 1
          if ( iglobal == nglobal .and. iprint == 1 ) then
            write(*,"(a,i5,a,i5)") '  Repeated best solution found ', nglobal,' times at trial ', itrial
          end if
        end if
      end if
      if(iprint.eq.2) write(*,trial_format) itrial, score, nbij, ngaps, iglobal

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
 
    scorebest = -1.d0
    iglobal = 0
    itrial = 0
    do while( iglobal < nglobal .and. itrial <= maxtrial ) 
      itrial = itrial + 1

      ! Writes the titles for regular iteration printing

      if(iprint.eq.2) write(*,title_format)

      ! Create random initial point of itrial > 1 

      if ( itrial > 1 ) call randomini(na,nb,prota,protb,nbij,bije)

      ! Number iterations, functional evaluations and DP calls

      nef = 0 
      it = 0
      prevscore = -1.d0
 
      ! Compute the correspondence and the score at initial point          

      call triang(prota,protb,na,nb,dtri2,gap,bije,nbij,&
                  bijscore,ngaps,score)
      nef = nef + 1
      if(iprint.eq.2) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
            
      ! Here begin the iteration loop

      do while(it.le.maxit.and.(score-prevscore).gt.abs(tol*score))

        it = it + 1
        prevscore = score
 
        ! Given the correspondence, perform a Procrustes RMSD alignment

        call procrustes(nbij,na,bije,prota,protb)

        ! Compute the DP bijection and the score at the new orientation

        call triang(prota,protb,na,nb,dtri2,gap,bije,nbij,bijscore,ngaps,&
                    score)

        ! Output regular iteration data
            
        if(iprint.eq.2) write(*,data_format) it, score, score-prevscore, nbij, ngaps, nef

      end do
      if ( score > scorebest ) then
        if ( abs(score-scorebest) > tol*score ) then
          iglobal = 1
        end if
        nbijbest = nbij
        scorebest = score
        do i = 1, nbij
          bijebest(i,1) = bije(i,1)
          bijebest(i,2) = bije(i,2)
        end do
        do i = 1, na
          prota_best(i,1) = prota(i,1)
          prota_best(i,2) = prota(i,2)
          prota_best(i,3) = prota(i,3)
        end do
        if(iprint.eq.1) write(*,trial_format) itrial, score, nbij, ngaps, iglobal
      else
        if ( scorebest-score <= tol*score ) then
          iglobal = iglobal + 1
          if ( iglobal == nglobal .and. iprint == 1 ) then
            write(*,"(a,i5,a,i5)") '  Repeated best solution found ', nglobal,' times at trial ', itrial
          end if
        end if
      end if
      if(iprint.eq.2) write(*,trial_format) itrial, score, nbij, ngaps, iglobal

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

    scorebest = -1.d0
    iglobal = 0
    itrial = 0
    do while( iglobal < nglobal .and. itrial <= maxtrial ) 
      itrial = itrial + 1

      ! Writes the titles for regular iteration printing

      if(iprint.eq.2) write(*,title_format)

      ! Create random initial point of itrial > 1 

      if ( itrial > 1 ) call randomini(na,nb,prota,protb,nbij,bije)

      ! Number iterations and functional evaluations

      nef = 0 
      it = 0
      prevscore = -1.d0
 
      ! Compute the correspondence and the score at initial point          

      call nonbscore(na, nb, prota, protb, dtri2, gap,&
                     disord, indisord, it, pair,&
                     score, bije, nbij, ngaps)
      nef = nef + 1
      if(iprint.eq.2) write(*,data_format) it, score, 0.d0, nbij, ngaps, nef
            
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
            
        if(iprint.eq.2) write(*,data_format) it, score, score-prevscore, nbij, ngaps, nef

      end do
      if ( score > scorebest ) then
        if ( abs(score-scorebest) > tol*score ) then
          iglobal = 1
        end if
        nbijbest = nbij
        scorebest = score
        do i = 1, nbij
          bijebest(i,1) = bije(i,1)
          bijebest(i,2) = bije(i,2)
        end do
        do i = 1, na
          prota_best(i,1) = prota(i,1)
          prota_best(i,2) = prota(i,2)
          prota_best(i,3) = prota(i,3)
        end do
        if(iprint.eq.1) write(*,trial_format) itrial, score, nbij, ngaps, iglobal
      else
        if ( scorebest-score <= tol*score ) then
          iglobal = iglobal + 1
          if ( iglobal == nglobal .and. iprint == 1 ) then
            write(*,"(a,i5,a,i5)") '  Repeated best solution found ', nglobal,' times at trial ', itrial
          end if
        end if
      end if
      if(iprint.eq.2) write(*,trial_format) itrial, score, nbij, ngaps, iglobal
    end do

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                               !
  ! POST ANALYSIS AND REPORT                      !
  !                                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Saving the prota_best file in prota to return best model

  do i = 1, na
    prota(i,1) = prota_best(i,1)
    prota(i,2) = prota_best(i,2)
    prota(i,3) = prota_best(i,3)
  end do
 
  ! Writting the final bijection obtained

  if(iprint.eq.1 .and. .not. seqoff) then
    if ( method /= 4 ) then
      call writebije(na,nb,resa,resb,numa,numb,bijebest,nbijbest)
    else
      call writenonbije(bijebest,nbijbest)
    end if
  end if
        
  ! Computing the RMSD of aligned residues at the solution

  call getrmsd(prota,protb,bijebest,nbijbest,rmsd)

  ! Computing the GDT scores at the solution

  call computegdt(na,nb,prota,protb,bijebest,nbijbest,gdt_threshold,gdt_ts,gdt_ha)
  call writermsf(na,nb,prota,protb,bijebest,nbijbest,&
                 numa,rmsf,rmsfout,rmsftrend,rmsftrendout)
 
  ! Printing the final score

  if(iprint.eq.1) then
    write(*,dash_line)
    write(*,"(a14,tr1,f12.6,tr1,a10,tr1,i5,tr1,a6,tr1,f10.6,tr1,a6,i4)")&
            '  FINAL SCORE:', scorebest,' COVERAGE:', nbijbest,' RMSD:', rmsd,' GAPS:', ngaps
    if ( method == 2 ) then
      write(*,"(a,f12.6)") "  Final score normalized by smallest protein: ",& 
                           scorebest*max(na,nb)/min(na,nb)
    end if
  endif
 
  ! Compute rmsd for atoms which are closer than some tolerance

  call getrmsd2(prota,protb,bijebest,nbijbest,rmsd2,nbij_dtri,dtri)
 
  ! Printing the final score

  if(iprint.eq.1) then
    write(*,dash_line)
    write(*,"(a,f8.4,a,f10.6,a,i6)")&
          '  ATOMS CLOSER THAN ',dtri,' Ang: RMSD: ',rmsd2,' COVERAGE: ', nbij_dtri
    write(*,"(a,f8.3,t34,a,f8.3)")&
          '  GDT_TS SCORE: ', gdt_ts, ' GDT_HA SCORE: ', gdt_ha

  endif

  ! Alignment time

  time1 = etime(tarray) - time1
  
  ! Printing concise output for database comparisons

  if(iprint.eq.0) then
    write(*,concise_format)&
            protea(ic(protea):length(protea)),&
            proteb(ic(proteb):length(proteb)),&
            scorebest, nbijbest, rmsd, nbij_dtri, rmsd2, gdt_ts, gdt_ha, time1
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
