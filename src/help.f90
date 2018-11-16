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
  write(*,*) ' -c2 A  Consider only chain A of protein 2 '
  write(*,*) ' -beta1 Consider atoms with beta > 0 in protein 1'
  write(*,*) ' -beta2 Consider atoms with beta > 0 in protein 2'
  write(*,*) ' -ocup1 Consider atoms with occupancy > 0 in protein 1'
  write(*,*) ' -ocup2 Consider atoms with occupancy > 0 in protein 2'
  write(*,*) ' -g [real] Penalization for gaps '
  write(*,*) ' -dtri [real] Atoms farther than this will not be'
  write(*,*) '        considered. Distance for Triangular score.'
  write(*,*) ' -rmsf [file] Print RMSF plot to file.'
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
