      subroutine diag(a,n,nn,w,z)
* "nn" represents the size of matrix stored in "a" whose size is "n"
* "n" should match the dimension of "a" as defined in the calling
* routine
*
* Compilation syntax:  ifort source.f diag.f -mkl
*
      DIMENSION ISUPPZ(N), IWORK(1000)
      DOUBLE PRECISION A(N,N), W(N), Z(N,N),
     $                 WORK(1000)
      ABSTOL = -1.0
      IL = 1
      IU = Nn
      LDA = N
      LDZ = N 
      LWMAX = 1000

c      write(*,*)'The following matrix was considered'
c      do i=1,nn
c        write (*,'(2x,100f9.4)')(a(i,j),j=1,nn)
c      end do
*
*     Query the optimal workspace.
*
      LWORK = -1
      LIWORK = -1
      CALL DSYEVR( 'Vectors', 'All', 'Upper', Nn, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL DSYEVR( 'Vectors', 'All', 'Upper', Nn, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        STOP
c      ELSE
c        WRITE(*,'(/A,I2)')' The total number of eigenvalues found:', M
c        WRITE(*,'(1000F8.4)')(W(I),I=1,M)
c        WRITE(*,'(/A)')' The eigenvectors are:'
c        DO I = 1, M
c          WRITE(*,'(1000F8.4)')(Z(I,J),J=1,M)
c        END DO
      END IF
      END
