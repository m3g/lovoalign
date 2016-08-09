! Module containing the input parameters

module inputpars

  integer :: method, maxit, iprint, nglobal
  double precision :: gap, dtri, gdt_threshold 
  character(len=1) :: chaina, chainb
  character(len=200) :: protea, proteb, pdblist, pdbout
  logical :: output, useini
  logical :: all, seqoff, seqfix
  logical :: beta1, beta2, ocup1, ocup2
  logical :: rmsf, rmsftrend
  character(len=200) :: rmsfout, rmsftrendout

end module inputpars
