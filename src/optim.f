       
cccccccccccccccccccccccccccccccccc 
c                                c
c The Newton optimization method c
c                                c
cccccccccccccccccccccccccccccccccc 

      subroutine newton(evalf,na,nb,prota,protb,score,bije,bijscore,
     +                  dzero2,scale,nbij,gap,ngaps,nef,gnor)
        
      use sizes
      implicit none
      integer i, j, na, nb, bije(maxatom,2), nbij, ngaps, nef,
     +        isub
      double precision delta, grad(6),
     +                 hess(6,6), prota(maxatom,3), protb(maxatom,3), 
     +                 gnor,score, scortrial, ared, pred, factor, 
     +                 factor2, protrial(maxatom,3), x(6), gap, 
     +                 dzero2, tol, scale, bijscore(maxatom)
      external evalf

c Maximum size of step x and convergence tolerance

      delta = 12.d0
      tol = 1.d-9

c Factors used for gradient and hessian computations 

      factor = 1.d0/dzero2
      factor2 = factor*factor

c Main Newton loop, iterates until a minimum is found

      pred = 1.d20

c Compute the gradient and the hessian

      call gradhes(grad,hess,nbij,bije,prota,protb,factor,
     +             factor2,scale)

c Subproblem iteration, while sufficient decrease is not obtained

      ared = -1.d0
      isub = 1
      do while(ared.le.0.1d0*pred.and.pred.gt.abs(tol*score))
    
        call subprob(grad,hess,delta,x,pred,isub,ared,scale)
        isub = 2

c  Compute trial protein A, called protrial here

        do i = 1, na
          protrial(i,1) = prota(i,1)
          protrial(i,2) = prota(i,2)
          protrial(i,3) = prota(i,3)
        end do
        call moveprot(x,na,protrial)
      
c  Compute score at trial point
                                              
        call evalf(protrial,protb,na,nb,dzero2,gap,bije,nbij,
     +             bijscore,ngaps,scortrial)
        nef = nef + 1
        ared = scortrial - score 
      
      end do

c Update protein A as the accepted trial point

      do i = 1, na
        prota(i,1) = protrial(i,1)
        prota(i,2) = protrial(i,2)
        prota(i,3) = protrial(i,3)
      end do
      score = scortrial

c Computing the norm of the gradient for report

      gnor = 0.d0
      do j = 1, 6
        gnor = dmax1(gnor,dabs(grad(j)))
      end do
      
      return
      end

c
c Subroutine that computes the gradient and the hessian
c

      subroutine gradhes(grad,hess,nbij,bije,prota,protb,
     +                   factor,factor2,scale)

      use sizes
      implicit none

      double precision grad(6), hess(6, 6), d2,
     +  prota(maxatom,3), protb(maxatom,3),    
     +  factor, factor2, atoma(3), atomb(3), gd2(6),
     +  hd2(6, 6), scale
      
      integer bije(maxatom,2), nbij, i, j, ii    

c  Compute gradient and Hessian of  -(structal score) at x = 0.

      do j = 1, 6
      grad(j) = 0.d0
      end do
      do j = 1, 6
      do i = 1, 6
      hess(i, j) = 0.d0 
      end do
      end do

      do ii = 1, nbij 
      do j = 1, 3
      atoma(j) = prota(bije(ii,1), j)
      atomb(j) = protb(bije(ii,2), j)
      end do

      d2 = 0.d0
      do i = 1, 3
      d2 = d2 + (atoma(i)-atomb(i))**2
      end do

      do j = 1, 3
      gd2(j) = atoma(j)-atomb(j)
      end do
      gd2(4) = atoma(2)*atomb(1) - atoma(1)*atomb(2)
      gd2(5) = atoma(3)*atomb(1) - atoma(1)*atomb(3)
      gd2(6) = atoma(3)*atomb(2) - atoma(2)*atomb(3)

      do i = 1, 6
      gd2(i) = 2.d0 * gd2(i)
      end do

      do j = 1, 6
      do i = 1, 6
      hd2(i, j) = 0.d0
      end do
      end do

      do i = 1, 3
      hd2(i, i) = 1.d0
      end do

      hd2(1, 4) = -atoma(2)
      hd2(1, 5) = -atoma(3)
      hd2(2, 4) = atoma(1)
      hd2(2, 6) = -atoma(3)
      hd2(3, 5) = atoma(1)    
      hd2(3, 6) = atoma(2)

      hd2(4, 4) = atoma(1)*atomb(1) + atoma(2)*atomb(2)
      hd2(4, 5) = atoma(2)*atomb(3)
      hd2(4, 6) = -atoma(1)*atomb(3)
      hd2(5, 5) = atoma(1)*atomb(1) + atoma(3)*atomb(3)
      hd2(5, 6) = atoma(1)*atomb(2)
      hd2(6, 6) = atoma(2)*atomb(2) + atoma(3)*atomb(3)

      do j = 1, 5
      do i = j+1, 6
      hd2(i, j) = hd2(j, i)
      end do
      end do

      do j = 1, 6
      do i = 1, 6
      hd2(i, j) = 2.d0 * hd2(i, j)
      end do
      end do

      do i = 1, 6
      grad(i) = grad(i) + (-factor/(1.d0 + factor*d2)**2)*gd2(i)
      end do 

      do j = 1, 6
      do i = 1, 6
      hess(i, j) = hess(i, j) + 
     *   (2.d0*factor2/(1.d0+factor*d2)**3)*gd2(i)*gd2(j)
     *    - factor * hd2(i, j)/(1.d0+ factor * d2)**2
      end do
      end do

      end do

      do j = 1, 6
        do i = 1, 6
          hess(i,j)=-scale*hess(i,j)
        end do
        grad(j) =-scale*grad(j)
      end do
    
      return
      end

      subroutine subprob(grad,hess,delta,x,pred,isub,ared,scale)

      implicit none
      double precision grad(6), hess(6, 6), x(6), xnor,
     *     delta, pred, savehess(6, 6), auxdiag(6), alpha, atmp, ared,
     +     scale, scale40
      integer isub, i, j, ichol
 
c   Solves the Unconstrained quadratic subproblem
      if(isub.eq.1) then


      do j = 1, 6
      do i = 1, 6
      savehess(i, j) = hess(i, j)
      end do
      end do      

      call chole (6, hess, ichol, auxdiag, 6)

      scale40 = scale/40.d0 
      do while(ichol.ne.0)
      do j = 1, 6
      do i = 1, 6
      if(i.eq.j) then
      hess(i, j) = (1.d0 - scale40) * savehess(i, j) + scale40
      else
      hess(i, j) = (1.d0 - scale40) * savehess(i, j)
      endif
      savehess(i, j) = hess(i, j)
      end do
      end do
      call chole (6, hess, ichol, auxdiag, 6)
      end do

      call sicho (6, hess, x, grad, auxdiag, 6)
      pred = 0.d0
      do i = 1, 6
      pred = pred + x(i) * grad(i)
      x(i) = -x(i)
      end do
      
      call norma (x, xnor)
      if(xnor.gt.delta) then
      do i = 1, 6
      x(i) = x(i)*delta/xnor
      end do 
      pred = pred*delta/xnor
      endif

      else

      atmp =  pred  /( 2.0d0 * ( -ared +  pred ) )

      if ( atmp .lt. 1.d-2 .or. atmp .gt. 0.9d0 ) then
          alpha = 0.5d0 

      else
          alpha = atmp
      end if

      do i = 1, 6
      x(i) = x(i) * alpha
      end do

      pred = pred * alpha
c      write(*, *)' isub = ', isub,
c     *   ' atmp =', atmp, ' alpha:', alpha,' pred:', pred
      endif
      return
      end      

      subroutine norma(x, xnor)
      implicit none
      double precision x(6), xnor, xi
      integer i 
      xnor = 0.d0
      do i = 1, 6
      xi = x(i)
      xnor = xnor + xi * xi
      end do
      xnor = dsqrt(xnor)
      return
      end

c
c Subroutine cholesky
c

      subroutine chole (n, a, ier, diag, nlin)
      implicit none
      integer ier, i, nlin, j, k, ii, n
      double precision a(nlin, n), diag(n), z, temp
      ier = 0
c Save the diagonal of the matrix
      do i = 1, n
      diag(i) = a(i, i)
      end do
c Test nonnegativity of diagonal
      do i = 1, n
      if(diag(i).le.0) then
      ier = 1
      return
      endif
      end do

      a(1,1) = dsqrt(a(1, 1))
      if(n.eq.1)return

      do 1 i = 2, n
      do 2 j = 1 ,i-1
      z = 0.d0
      if(j.gt.1)then
      do 3 k=1,j-1
3      z = z + a(i,k) * a(j,k)
      endif
      a(i,j) = (a(i,j) - z)/a(j,j)
2      continue
      z = 0.d0
      do 4 j=1,i-1
4      z = z + a(i,j)**2
      temp = a(i, i) - z

c   Test positive definiteness
      if(temp.le.0.d0) then
      ier = i
c   Restore the diagonal
      do ii = 1, n
      a(ii, ii) = diag(ii)
      end do
c   Restore lower triangular part
      do ii = 2, n
      do j = 1, ii-1
      a(ii, j) = a(j, ii)
      end do
      end do
      return
      endif

      a(i,i) = dsqrt(temp)
1      continue
      return
      end

      subroutine sicho (n, a, x, b, aux, nlin)
      implicit none
      integer n, i, nlin, j          
      double precision a(nlin, n),x(n),b(n),aux(n), z

      aux(1) = b(1)/a(1,1)
      if(n.gt.1)then
      do 1 i=2,n
      z = 0.d0
      do 2 j=1,i-1
2      z = z + a(i,j)*aux(j)
      aux(i) = (b(i) - z) / a(i,i)
1      continue
      endif
      x(n) = aux(n)/a(n,n)
      if(n.eq.1)return
      do 3 i=n-1,1,-1
      z = 0.d0
      do 4 j=i+1,n
4      z = z + a(j,i)*x(j)
      x(i) = (aux(i) - z)/a(i,i)
3      continue
      return
      end



