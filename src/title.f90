!
! Print title with version
!
subroutine title()

  implicit none
  write(*,"(/,&
     &'  ',71('#'),/,&
     &/,&
     &'                           PROGRAM: LOVOALIGN',/,&
     &'                   Low Order Value Optimization Methods',/,&
     &'                          for Protein Alignment',/,&
     &'                              Version 16.327.2 ',/&
     &/,&
     &'  ',71('#'))")

end subroutine title

