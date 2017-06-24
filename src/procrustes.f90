!
! Subroutine procrustes: given two sets of vectors, finds the best
!                   rotation that to align the two sets. 
!
!                 Method: S. K. Kearsley, 
!                         "On the orthogonal transformation used for
!                          structural comparisons"
!                         Acta Crystallog. A 45 (1989), 208-210 
!
!  Input: nbij: number of points of the bijection
!         bije: array containing the bijection indices
!         yref: reference vector set
!         xvar: variable vector set
!
!  Ouput: xvar vector aligned to y
!         If the original vector is xvar, then the transformed 
!         vector is xvar = u * (var - cmx) + cmy
!

subroutine procrustes(nbij,na,bije,xvar,yref)

  use sizes
  implicit none
  integer :: i, j, k, nbij, iq, na
  integer :: bije(maxatom,2)
  double precision :: x(maxatom,3), y(maxatom,3),&
                      cmx(3), cmy(3), q(4,4), xp(maxatom),&
                      yp(maxatom), zp(maxatom), xm(maxatom),&
                      ym(maxatom), zm(maxatom), &
                      a(4,4), u(3,3), xvar(maxatom, 3), vecaux(3),&
                      yref(maxatom,3)
  ! For dsyev
  double precision :: work(12)
  integer :: info

  ! Safeguard for the case in which the bijection contains only one atom

  if(nbij.eq.1) then
    vecaux(1) = xvar(bije(1,2),1) - yref(bije(1,1),1)
    vecaux(2) = xvar(bije(1,2),2) - yref(bije(1,1),2)
    vecaux(3) = xvar(bije(1,2),3) - yref(bije(1,1),3)
    do i = 1, na
      xvar(i,1) = xvar(i,1) - vecaux(1)
      xvar(i,2) = xvar(i,2) - vecaux(2)
      xvar(i,3) = xvar(i,3) - vecaux(3)
    end do
    return
  end if

  ! Copying the arrays to have only the atoms of the bijections

  do i = 1, nbij
    x(i,1) = xvar(bije(i,1),1)
    x(i,2) = xvar(bije(i,1),2)
    x(i,3) = xvar(bije(i,1),3)
    y(i,1) = yref(bije(i,2),1)
    y(i,2) = yref(bije(i,2),2)
    y(i,3) = yref(bije(i,2),3)
  end do
         
  ! Computing the centroid of the structures

  do i = 1, 3
    cmx(i) = 0.d0
    cmy(i) = 0.d0
  end do

  do j = 1, 3
    do i = 1, nbij
      cmx(j) = cmx(j) + x(i,j)
      cmy(j) = cmy(j) + y(i,j)  
    end do
  end do

  do i = 1, 3
    cmx(i) = cmx(i) / dfloat(nbij) 
    cmy(i) = cmy(i) / dfloat(nbij)
  end do

  ! Moving structures to their baricenters

  do j = 1, 3
    do i = 1, nbij 
      x(i,j) = x(i,j) - cmx(j)    
      y(i,j) = y(i,j) - cmy(j) 
    end do
  end do

  ! Computing the quaternion matrix

  do i = 1, nbij
    xm(i) = y(i,1) - x(i,1) 
    ym(i) = y(i,2) - x(i,2) 
    zm(i) = y(i,3) - x(i,3) 
    xp(i) = y(i,1) + x(i,1) 
    yp(i) = y(i,2) + x(i,2) 
    zp(i) = y(i,3) + x(i,3)
  end do 
 
  do j = 1, 4
    do i = 1, 4
      q(i,j) = 0.d0
    end do
  end do
 
  do i = 1, nbij  
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)
     
  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix

  call dsyev('V','U',4,q,4,a,work,12,info)

  ! Computing the rotation matrix

  iq = 1
  u(1,1) = q(1,iq)**2 + q(2,iq)**2 - q(3,iq)**2 - q(4,iq)**2
  u(1,2) = 2. * ( q(2,iq)*q(3,iq) + q(1,iq)*q(4,iq) )
  u(1,3) = 2. * ( q(2,iq)*q(4,iq) - q(1,iq)*q(3,iq) )
  u(2,1) = 2. * ( q(2,iq)*q(3,iq) - q(1,iq)*q(4,iq) )
  u(2,2) = q(1,iq)**2 + q(3,iq)**2 - q(2,iq)**2 - q(4,iq)**2
  u(2,3) = 2. * ( q(3,iq)*q(4,iq) + q(1,iq)*q(2,iq) )
  u(3,1) = 2. * ( q(2,iq)*q(4,iq) + q(1,iq)*q(3,iq) )
  u(3,2) = 2. * ( q(3,iq)*q(4,iq) - q(1,iq)*q(2,iq) )
  u(3,3) = q(1,iq)**2 + q(4,iq)**2 - q(2,iq)**2 - q(3,iq)**2

  ! Rotating-translating the whole protein

  do i = 1, na
    do j = 1, 3
      vecaux(j) = cmy(j)
      do k = 1, 3
        vecaux(j) = vecaux(j) + u(j,k) * (xvar(i,k)-cmx(k))
      end do
    end do
    do j = 1, 3
      xvar(i,j) = vecaux(j)
    end do
  end do

end subroutine procrustes
