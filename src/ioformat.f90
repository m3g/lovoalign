!
! Module containing ahestetic formats
!

module ioformat

  implicit none
  integer, parameter :: max_string_length = 200
  character(len=*), parameter :: dash_line = "('  ',71('-'))"
  character(len=max_string_length) :: hash_line
  character(len=1000) :: header_list
  character(len=max_string_length) :: title_format, data_format, trial_format
  character(len=max_string_length) :: concise_format
  integer :: max_filename_size

  contains

    subroutine outputformats()

      implicit none
      integer i, i1, i2, imax
      character(len=3) :: ichar(3)

      ! Formats for single runs

      title_format = "(t3,'ITER',t20,'SCORE',t30,'GRADIENT NORM',&
                      &t45,'COVERAGE',t56,'GAPS',t64,'NEF')"
      data_format = "(i6,tr1,e17.10,tr1,e17.10,tr4,i6,tr1,i6,tr1,i6)"
      trial_format = "(t3,'TRIAL: ',i7,' SCORE: ',f12.5,' COVERAGE: ',i6,' GAPS: ',i6,' GLOB:',i2)"

      ! Header for concise output printing (note that dtri is printed here)

      i1 = max_filename_size + 2 
      i2 = 2*max_filename_size + 4
      imax = i2 + 80
      write(ichar(1),"(i3)") i1
      write(ichar(2),"(i3)") i2
      write(ichar(3),"(i3)") imax
      ichar(1) = adjustl(ichar(1))
      ichar(2) = adjustl(ichar(2))
      do i = 1, min(imax,max_string_length)
        hash_line(i:i) = "#"
      end do
      do i = imax + 1, max_string_length
        hash_line(i:i) = " "
      end do

      write(header_list,*) "(",ichar(3),"('#'),/,&
               & '# LOVOALIGN ',/&
               &,'# http://www.ime.unicamp.br/~martinez/lovoalign',/,"&
               &,ichar(3),"('#'),/,&
               &'# Prot A: Variable protein: ',a,/,&
               &'# Prot B: Target (fixed) protein: ',a,/,&
               &'# PDB file list: ',a,/,&
               &'# Number of files in list: ',i8,/,&
               &'# SCORE: Best ',a,' score obtained. ',/,&
               &'# COV: Coverage (number of corresponding atoms).',/,&
               &'# RMSD: Root mean square deviation of COV atoms.',/,&
               &'# COV2: Number of atoms closer than ',f8.3,' Angstroms.',/,&
               &'# RMSD2: Root mean square deviation of COV2 atoms.',/,&
               &'# GDT_TS: Global Distance Test (GDT) total score.',/,&
               &'# GDT_HA: High-accuracy GDT score.',/,&
               &'# TIME: Time used in this alignment.',/,"&
               &,ichar(3),"('#'),/,&
               &'# Prot A',t",ichar(1),",'Prot B',t",ichar(2),",tr7,'SCORE',tr3,'COV',tr9,'RMSD',&
               &tr2,'COV2',tr8,'RMSD2',tr3,'GDT_TS',tr3,'GDT_HA',tr9,'TIME')"

      write(concise_format,"(a,a3,a,a3,a)") &
                           "(t1,a,t",adjustl(ichar(1)),",a,t",adjustl(ichar(2)),&
                           "f12.6,2(tr1,i5,tr1,f12.6),2(tr1,f8.3),tr1,f12.6)"

    end subroutine outputformats

end module ioformat
