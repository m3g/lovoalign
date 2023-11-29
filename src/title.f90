!
! Print title with version
!
subroutine title()

  implicit none
  write(*,"(/,&
     &'  ',71('#'),/,&
     &/,&
     &'                                LOVOALIGN',/,&
     &'             Convergent algorithms for protein structural alignment',/,&
     &'     Martinez, Andreani, and Martinez, BMC Bioinformatics, 8:306, 2007',/,&
     &'                              Version 22.0.2 ',/&
     &/,&
     &'  ',71('#'))")

end subroutine title



