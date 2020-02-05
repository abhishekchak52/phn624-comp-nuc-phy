program simpson_test_wf
      implicit none
      integer, parameter :: n = 1001, dp = selected_real_kind(14)
      integer :: i
      real(kind=dp), dimension(n) :: x, y
      real(kind=dp) :: sum, high, low, dx, x_sq

      interface
        function simpson(x,y,n)
                integer, parameter :: dp = selected_real_kind(14)
                integer, intent(in) :: n
                real(kind=dp) :: simpson
                real(kind=dp), dimension(n), intent(in) :: x, y
        end function simpson
      end interface


      
      sum=0.0d0
      high = 4.0d0; low = -high;
      dx = (high - low)/(n-1)
      do i=1, n
        x(i) = low + dfloat(i)*dx
        x_sq = x(i)*x(i)
        y(i) = (x_sq - 1.0d0)*exp(-0.5d0*x_sq)
        write(*, '(i, 2f10.3)') i, x(i), y(i)
      end do

      sum = simpson(x,y,n)
      write(*,*) sum, dsqrt(3.14d0)/2.0d0

end program simpson_test_wf
