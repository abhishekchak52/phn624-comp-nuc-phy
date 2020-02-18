program test_rk4
      use numerical_integration
      implicit none

      integer :: i
      integer, parameter :: n=100
      real(kind=dp), parameter :: pi = 4.0d0*atan(1.0d0)
      real(kind=dp), parameter :: dx = pi/(n-1)
      external :: deriv_sine, deriv_sho
      real(kind=dp) :: x(n) , y(n, 2)

      x(1)  = 0.0d0; y(1, :) = [ 1.0d0, 0.0d0] 
      call rk4(x, y, dx, n, 2, deriv_sho)
      do i=1,n
        write(*, *) x(i), y(i, 1), y(i, 2)
      end do
end program test_rk4
      
subroutine deriv_sine(x,y,dydx)
      use numerical_integration
      implicit none
      integer, parameter :: n_eq = 1
      real(kind=dp), intent(in) :: x, y(n_eq)
      real(kind=dp), dimension(n_eq), intent(out) :: dydx

      dydx = cos(x)

end subroutine deriv_sine

subroutine deriv_sho(x,y,dydx)
      use numerical_integration
      implicit none
      integer, parameter :: n_eq = 2
      real(kind=dp), intent(in) :: x, y(n_eq)
      real(kind=dp), dimension(n_eq), intent(out) :: dydx

      dydx(1) = y(2)
      dydx(2) = -y(1)

end subroutine deriv_sho
