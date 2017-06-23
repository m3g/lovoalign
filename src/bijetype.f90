!
! Module to define the sequence alignment type
!
module bijetype

  use sizes
  implicit none
  integer :: seqtype, fixnbij
  integer :: fixbije(maxatom,2)
  character(len=200) :: fastafile
  type fasta_sequence
    character(len=maxatom) :: seq
  end type fasta_sequence
  type(fasta_sequence) :: fasta(2)

end module bijetype
