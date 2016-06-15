       
      SUBROUTINE JACOBI(A,V,N) 
C*********************************************************************** 
C 
C     DIAGONALISATION OF REAL SYMMETRIC MATICES BY JACOBI METHOD 
C 
C*********************************************************************** 
      
      implicit none  
      integer i, j, k, n, m 
      double precision a(4, 4), v(4, 4), rho, tes, scl, omg, u, v1, v2
     *  ,v3, s, c, tem
      RHO=1.0d-12 
      TES=0.d0 
      SCL=0.d0 
      DO 10 I=1,N 
   10 SCL=SCL+A(I,I)**2 
      SCL=dsqrt(SCL)/DFLOAT(N) 
      DO 20 I=1,N 
      DO 20 J=1,N 
   20 A(I,J)=A(I,J)/SCL 
      DO 30 I=1,N 
      DO 30 J=1,N 
      V(I,J)=0.d0 
      IF(I.EQ.J)V(I,J)=1. 
   30 CONTINUE 
      DO 100 I=2,N 
      DO 100 J=1,I-1 
  100 TES=TES+2.*A(I,J)*A(I,J) 
      TES=dsqrt(TES) 
      M=0 
  105 TES=TES/DFLOAT(N) 
      IF(TES.LT.RHO)TES=RHO 
  110 DO 165 I=2,N 
      DO 165 J=1,I-1 
      IF(dABS(A(I,J))-TES)165,115,115 
  115 M=1 
      V1=A(J,J) 
      V2=A(I,J) 
      V3=A(I,I) 
      U=0.5*(V1-V3) 
      IF(dABS(U)-RHO)120,125,125 
  120 OMG=-1. 
      GO TO 130 
  125 OMG=-V2/dSQRT(V2*V2+U*U) 
      IF(U.LT.0.)OMG=-OMG 
  130 S=OMG/dSQRT(2.*(1.+dSQRT(1.-OMG*OMG))) 
      C=dSQRT(1.-S*S) 
      DO 160 K=1,N 
      IF(K-I)140,135,135 
  135 TEM=A(K,J)*C-A(K,I)*S 
      A(K,I)=A(K,J)*S+A(K,I)*C 
      A(K,J)=TEM 
      GO TO 155 
  140 IF(K-J)145,150,150 
  145 TEM=A(J,K)*C-A(I,K)*S 
      A(I,K)=A(J,K)*S+A(I,K)*C 
      A(J,K)=TEM 
      GO TO 155 
  150 TEM=A(K,J)*C-A(I,K)*S 
      A(I,K)=A(K,J)*S+A(I,K)*C 
      A(K,J)=TEM 
  155 TEM=V(K,J)*C-V(K,I)*S 
      V(K,I)=V(K,J)*S+V(K,I)*C 
      V(K,J)=TEM 
  160 CONTINUE 
      A(J,J)=V1*C*C+V3*S*S-2.*V2*S*C 
      A(I,I)=V1*S*S+V3*C*C+2.*V2*S*C 
      A(I,J)=(V1-V3)*S*C+V2*(C*C-S*S) 
  165 CONTINUE 
      IF(M-1)175,170,170 
  170 M=0 
      GO TO 110 
  175 IF(TES-RHO)180,180,105 
  180 DO 190 I=1,N 
      DO 190 J=1,N 
  190 A(I,J)=SCL*A(I,J) 
      RETURN 
      END 
 
