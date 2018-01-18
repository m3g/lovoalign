!
! Print problem specifications
!

subroutine printdata(protea,proteb,na,nb,chaina,&
                     chainb,method,gap,maxit,dtri,gdt_threshold,&
                     useini)

  use ioformat
  use bijetype
  use warnings
  implicit none
  integer :: na, nb, maxit, method, length, ic, i
  double precision :: gap, dtri, gdt_threshold
  character(len=1) :: chaina, chainb
  character(len=200) :: protea, proteb
  logical :: useini
       
  write(*,*) ' Problem specifications: '
  write(*,dash_line) 
  write(*,*) ' Protein A: ', protea(ic(protea):length(protea))
  write(*,*) ' Protein B: ', proteb(ic(proteb):length(proteb))
  write(*,*) ' Number of atoms: A:', na, ' B:', nb
  if ( nwarn > 0 ) then
    do i = 1, nwarn
      write(*,"(a,a)") "  Warning: ", trim(warn(i))
    end do
  end if
  if(chaina.ne.'#') write(*,*) ' Protein A chain: ', chaina
  if(chainb.ne.'#') write(*,*) ' Protein B chain: ', chainb
  if(method.eq.1) write(*,*) ' Will maximize the STRUCTAL score'
  if(method.eq.2) write(*,*) ' Will maximize the TM-SCORE '
  if(method.eq.3) write(*,*) ' Will maximize the TRIANGULAR score '
  if(method.eq.4) write(*,*) ' Will maximize the NON-BIJECTIVE TRIANGULAR score '
  write(*,*) ' Penalization for gaps: ', gap
  if(seqtype == 3) write(*,*) ' Sequence alignment given in fasta file: ', trim(adjustl(fastafile))
  write(*,*) ' Maximum number of iterations: ', maxit
  if(useini) write(*,*) ' Using internal-distance initial point.'
  if(method.eq.3) then
    write(*,*) ' Triangular score with cutoff: ', dtri
  end if
  write(*,*) ' GDT Threshold: ', gdt_threshold
  write(*,dash_line) 

end subroutine printdata
