PROGRAM shell_model
IMPLICIT NONE
  REAL,ALLOCATABLE,DIMENSION(:,:)::E
  INTEGER::i,j,N,l,k 
  REAL::C,D
  C=0.05
  D=0.0225
!To input the value of Quantum Number N
  PRINT*,"Enter the value of N:"
  READ*,N
  ALLOCATE(E(0:N,0:N+1))
  DO i=0,N
    E(i,0)=i+1.5  !contains energy values of SHO
    j=1
    DO l=i,0,-2
      IF(l.NE.0) THEN
!For Energy level splitting with j=l+1/2 and j=l-1/2 respectively
        E(i,j)=(i+1.5)-(C*l)-(D*l*(l+1))
        j=j+1
        E(i,j)=(i+1.5)+(C*(l+1))-(D*l*(l+1))
        j=j+1
      END IF
      IF(l.EQ.0) THEN
!For no splitting in energy levels
        E(i,j)=i+1.5
      END IF
    END DO
    PRINT*,(E(i,j),j=0,i+1)
  END DO
  PRINT*,"HEY"

  OPEN(10,FILE="single_particle_Shell_model.dat")
  DO i=0,N
    DO j=1,i+1
      DO k=0,3
        WRITE(10,*) k,E(i,0)
      END DO
      DO k=4,7
        WRITE(10,*) k,E(i,j)
      END DO
      WRITE(10,*)
    END DO
  END DO
  CLOSE(10)
  PRINT*,"HEY"
END PROGRAM shell_model
