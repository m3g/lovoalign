!
! Subroutine writepdb: Reads all atoms of protein (not only CAs),
!                      translates and rotates them according to
!                      the best alignment obtained and prints the
!                      output pdb file with the aligned atoms. 
!

subroutine writepdb(pdbout,protea,prota,chaina,beta1,ocup1,&
                    rmin1,rmax1,na,resa,numa,proteb,all)

  use sizes
  implicit none
      
  integer :: i, j, na, nca, iq, length, ic, numa(maxatom), ioerr
  integer :: rmin1, rmax1
  double precision :: prota(maxatom,3), aux(maxatom,3), cmo(3),&
                      cma(3), xm(maxatom), ym(maxatom), zm(maxatom),&
                      xp(maxatom), yp(maxatom), zp(maxatom), q(4,4),&
                      a(4,4), qmin, u(3,3), xn(3)
  character(len=200) pdbout, protea, record, proteb
  character(len=1) chaina, resa(maxatom)
  logical :: all, mark, error, beta1, ocup1

  ! Reading the original CA coordinates

  call readfile(protea,aux,chaina,beta1,ocup1,rmin1,rmax1,nca,resa,numa,all,error)
  if(nca.ne.na) then
    write(*,*) ' ERROR reading the file of the first'
    write(*,*) ' protein when going to write the output file.'
    write(*,*) ' Wrong number of atoms.'
    write(*,*) ' Output file not written.'
    return
  end if 

! If the coordinates were given in xyz format we can have an error here
 
  if(nca.lt.3) then
    write(*,*) ' ERROR: The output option is only available for files in PDB format. '
    return
  end if
 
! Translating the original and aligned proteins to the origin

  do i = 1, 3
    cmo(i) = 0.d0
    cma(i) = 0.d0
  end do
  do i = 1, na
    cmo(1) = cmo(1) + aux(i,1)
    cmo(2) = cmo(2) + aux(i,2)
    cmo(3) = cmo(3) + aux(i,3)
    cma(1) = cma(1) + prota(i,1)
    cma(2) = cma(2) + prota(i,2)
    cma(3) = cma(3) + prota(i,3)
  end do
  do i = 1, 3
    cmo(i) = cmo(i) / dfloat(na)
    cma(i) = cma(i) / dfloat(na)
  end do
  do i = 1, na
    aux(i,1) = aux(i,1) - cmo(1)
    aux(i,2) = aux(i,2) - cmo(2)
    aux(i,3) = aux(i,3) - cmo(3)
    prota(i,1) = prota(i,1) - cma(1)
    prota(i,2) = prota(i,2) - cma(2)
    prota(i,3) = prota(i,3) - cma(3)
  end do

  ! Obtaining the best rotation that alignes them (this is almost a whole
  ! procrustes iteration of subroutine procrustes.

  do i = 1, na
    xm(i) = prota(i,1) - aux(i,1) 
    ym(i) = prota(i,2) - aux(i,2) 
    zm(i) = prota(i,3) - aux(i,3) 
    xp(i) = prota(i,1) + aux(i,1) 
    yp(i) = prota(i,2) + aux(i,2) 
    zp(i) = prota(i,3) + aux(i,3)
  end do 

  do j = 1, 4
    do i = 1, 4
      q(i,j) = 0.
    end do
  end do

  do i = 1, na
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

  call jacobi(q,a,4)

  ! Choosing the quaternion that corresponds to the minimum

  qmin = 1.e30
  do i = 1, 4
    if(q(i,i).lt.qmin) then
      iq = i
      qmin = q(i,i)
    end if
  end do
  
! Computing the rotation matrix

  u(1,1) = a(1,iq)**2 + a(2,iq)**2 - a(3,iq)**2 - a(4,iq)**2
  u(1,2) = 2. * ( a(2,iq)*a(3,iq) + a(1,iq)*a(4,iq) )
  u(1,3) = 2. * ( a(2,iq)*a(4,iq) - a(1,iq)*a(3,iq) )  
  u(2,1) = 2. * ( a(2,iq)*a(3,iq) - a(1,iq)*a(4,iq) )  
  u(2,2) = a(1,iq)**2 + a(3,iq)**2 - a(2,iq)**2 - a(4,iq)**2 
  u(2,3) = 2. * ( a(3,iq)*a(4,iq) + a(1,iq)*a(2,iq) )  
  u(3,1) = 2. * ( a(2,iq)*a(4,iq) + a(1,iq)*a(3,iq) )  
  u(3,2) = 2. * ( a(3,iq)*a(4,iq) - a(1,iq)*a(2,iq) )  
  u(3,3) = a(1,iq)**2 + a(4,iq)**2 - a(2,iq)**2 - a(3,iq)**2 
 
  ! Now we need to read the pdb file again and apply the transformations
  ! resulted in the best alignment for all atoms of the protein

  open(10,file=protea(1:length(protea)),status='old')
  open(20,file=pdbout(1:length(pdbout)))
  write(20,"('REMARK Protein ',a,' aligned to',/,&
            &'REMARK protein ',a,/,&
            &'REMARK with LovoAlign: ',&
            &'http://www.ime.unicamp.br/~martinez/lovoalign')")&
        protea(ic(protea):length(protea)), proteb(ic(proteb):length(proteb))

  mark = .false.
  do while(.true.)
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(length(record).lt.82) then
      do i = length(record)+1, 82
        record(i:i) = ' '
      end do
    end if
    if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
      read(record(31:38),*,iostat=ioerr) aux(1,1)
      if ( ioerr /= 0 ) error = .true.
      read(record(39:46),*,iostat=ioerr) aux(1,2)
      if ( ioerr /= 0 ) error = .true.
      read(record(47:54),*,iostat=ioerr) aux(1,3)
      if ( ioerr /= 0 ) error = .true.
      if ( error ) then
        write(*,*) ' ERROR: Failed reading coordinates in file ',&
                   protea(ic(protea):length(protea)),&
                   ': output file may be not complete.' 
        exit
      end if
      aux(1,1) = aux(1,1) - cmo(1)
      aux(1,2) = aux(1,2) - cmo(2)
      aux(1,3) = aux(1,3) - cmo(3)
      xn(1) = u(1,1)*aux(1,1)+u(1,2)*aux(1,2)+u(1,3)*aux(1,3)
      xn(2) = u(2,1)*aux(1,1)+u(2,2)*aux(1,2)+u(2,3)*aux(1,3)
      xn(3) = u(3,1)*aux(1,1)+u(3,2)*aux(1,2)+u(3,3)*aux(1,3)
      xn(1) = xn(1) + cma(1)
      xn(2) = xn(2) + cma(2)
      xn(3) = xn(3) + cma(3)
      write(20,"(a30,3(f8.3),a28)") record(1:30), xn(1), xn(2), xn(3), record(55:82)
    else
      write(20,"(t1,a82)") record(1:82)
    end if
  end do
  close(10)
  close(20)

end subroutine writepdb
