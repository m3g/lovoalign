!
! Subroutine that writes RMSF data 
!

subroutine writermsf(na,nb,prota,protb,bije,nbij,&
                     numa,rmsf,rmsfout,rmsftrend,rmsftrendout)

  use sizes
  use ioformat
  implicit none
  integer i, j, nbij, bije(maxatom,2), na, nb, numa(maxatom), trendlist(100)
  double precision ::prota(maxatom,3), protb(maxatom,3), dist
  logical :: rmsf, rmsftrend
  character(len=200) :: rmsfout, rmsftrendout

  if ( rmsf ) then
    do i = 1, 100
      trendlist(i) = 0
    end do
    open(10,file=rmsfout)
    write(10,"('# RESIDUE_NUMBER             RMSF')")
    do i = 1, nbij
      dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
           + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
           + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
      dist = dsqrt(dist)

      j = 100
      do while( dble(j) / 2.d0 > dist )  
        trendlist(j) = trendlist(j) + 1
        j = j - 1
      end do
      write(10,"(tr8,i8,tr5,f12.5)") numa(bije(i,1)), dist 
    end do
    close(10)
    write(*,dash_line)
    write(*,*) ' Wrote RMSF data file: ', trim(adjustl(rmsfout))
  end if

  if ( rmsftrend ) then 
    open(10,file=rmsftrendout)
    write(10,"('#            RMSF         FRACTION')")
    j = 1
    do while(trendlist(j) < nbij)
      write(10,"(tr2,f15.5,tr2,f15.5)") dble(j)/2.d0, dble(trendlist(j)) / min(na,nb)
      j = j + 1
    end do
    write(*,dash_line)
    write(*,*) ' Wrote RMSF-trend data file: ', trim(adjustl(rmsftrendout))
    close(10)
  end if

end subroutine writermsf
