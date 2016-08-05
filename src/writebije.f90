!        
! Subroutine writebije: Writes the bijection obtained at the
!                       solution

subroutine writebije(na,nb,resa,resb,numa,numb,bije,nbij)

  use sizes
  implicit none
      
  integer :: na, nb, bije(maxatom,2), nbij, nlines,&
             numa(maxatom), numb(maxatom), i, j, iamin, ibmin,&
             ia, ib, ina, inb, ic
  character(len=1) :: resa(maxatom), resb(maxatom)
  character(len=50) :: line1, line2

  write(*,"(/,'  ',25('-'),' SEQUENCE ALIGNMENT ',26('-'))")

  nlines = (na + nb - nbij - 1) / 50 + 1

  iamin = 1
  ibmin = 1
  ia = 1
  ib = 1
  ina = 1
  inb = 1
  do i = 1, nlines
    do j = 1, 50
      line1(j:j) = ' '
      line2(j:j) = ' '
    end do
    ic = 1
    do while(ic.le.50)
      if(ia.lt.bije(ina,1)) then
        line1(ic:ic) = resa(ia)
        line2(ic:ic) = '-'
        ia = ia + 1
      else if(ib.lt.bije(inb,2)) then
        line1(ic:ic) = '-'
        line2(ic:ic) = resb(ib)
        ib = ib + 1
      else
        if(ia.le.na) line1(ic:ic) = resa(ia)
        if(ib.le.nb) line2(ic:ic) = resb(ib)
        if(ia.gt.na.and.ib.le.nb) line1(ic:ic) = '-'
        if(ib.gt.nb.and.ia.le.na) line2(ic:ic) = '-'
        ia = ia + 1
        ib = ib + 1
        if(ina.lt.nbij) ina = ina + 1
        if(inb.lt.nbij) inb = inb + 1
      end if
      ic = ic + 1
    end do 
    write(*,"('            .         .         .         .         .')")
    write(*,"(tr6,i5,tr1,a50,tr1,i5)") numa(min0(iamin,na)),&
                                       line1,&
                                       numa(max0(1,min0(ia-1,na)))
    write(*,"(tr6,i5,tr1,a50,tr1,i5)") numb(min0(ibmin,nb)),&
                                       line2,&
                                       numb(max0(1,min0(ib-1,nb)))
    if(ia.gt.iamin) iamin = ia 
    if(ib.gt.ibmin) ibmin = ib 
  end do
  write(*,*)

end subroutine writebije
