! Module containing the input parameters

module inputpars

  integer :: method, maxit, iprint, nglobal, maxtrial
  integer :: rmin1, rmin2, rmax1, rmax2
  double precision :: gap, dtri, gdt_threshold 
  character(len=1) :: chaina, chainb
  character(len=200) :: protea, proteb, pdblist, pdbout
  logical :: output, useini
  logical :: all, seqoff, seqfix
  logical :: beta1, beta2, ocup1, ocup2
  logical :: rmsf, rmsftrend
  logical :: skip
  character(len=200) :: rmsfout, rmsftrendout

end module inputpars
