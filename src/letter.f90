!
! One letter amino acid codes for printing sequence alignment
!

character function letter(resid_in)

  implicit none
  character*5 resid_in
  character*3 resid

  resid = trim(adjustl(resid_in))
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
