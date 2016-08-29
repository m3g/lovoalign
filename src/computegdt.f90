!
! Subroutine that computes the GDT scores
!

subroutine computegdt(na,nb,prota,protb,bije,nbij,gdt_threshold,gdt_ts,gdt_ha)
 
  use sizes
  implicit none
  integer :: i, nbij, bije(maxatom,2), na, nb
  double precision :: gdt_ts, gdt_ha, gdt_threshold, gdt_ha_threshold
  double precision :: prota(maxatom,3), protb(maxatom,3), dist

  gdt_ha_threshold = gdt_threshold / 2.d0
  gdt_ts = 0.d0
  gdt_ha = 0.d0
  do i = 1, nbij
    dist = (prota(bije(i,1),1) - protb(bije(i,2),1))**2 &
         + (prota(bije(i,1),2) - protb(bije(i,2),2))**2 &
         + (prota(bije(i,1),3) - protb(bije(i,2),3))**2
    dist = dsqrt(dist)

    ! GDT_TS score (with threshold = 4.d0 Angs by default)
    if ( dist < (gdt_threshold/4.d0) ) gdt_ts = gdt_ts + 1.d0 
    if ( dist < (gdt_threshold/2.d0) ) gdt_ts = gdt_ts + 1.d0
    if ( dist < gdt_threshold ) gdt_ts = gdt_ts + 1.d0
    if ( dist < (2.d0*gdt_threshold) ) gdt_ts = gdt_ts + 1.d0

    ! GDT-"High Accuracy" score (with threshold = 2.d0 Angs by default)
    if ( dist < (gdt_ha_threshold/4.d0) ) gdt_ha = gdt_ha + 1.d0 
    if ( dist < (gdt_ha_threshold/2.d0) ) gdt_ha = gdt_ha + 1.d0
    if ( dist < gdt_ha_threshold ) gdt_ha = gdt_ha + 1.d0
    if ( dist < (2.d0*gdt_ha_threshold) ) gdt_ha = gdt_ha + 1.d0
  end do
  gdt_ts = 100.d0 * gdt_ts / (4.d0*min(na,nb))
  gdt_ha = 100.d0 * gdt_ha / (4.d0*min(na,nb))

end subroutine computegdt
