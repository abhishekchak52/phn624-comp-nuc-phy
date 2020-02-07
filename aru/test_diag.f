      implicit real*8 (a-h,o-z)
      parameter(nmax=100)
      dimension a(nmax,nmax),w(nmax),z(nmax,nmax)
      read(*,*)nn
      do i=1,nn
        read(*,*)(a(i,j),j=1,nn)
      end do
      call diag(a,100,nn,w,z)
      write(*,'(/a,i2)')' The total number of eigenvalues found:', nn
      write(*,'(100f9.5)')(w(i),i=1,nn)
      write(*,'(/a)')' The eigenvectors are:'
      do i=1,nn
        write(*,'(2x,100f9.5)')(z(i,j),j=1,nn)
      end do
      end 
