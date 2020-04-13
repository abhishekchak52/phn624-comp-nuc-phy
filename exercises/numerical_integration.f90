module numerical_integration
        use numerical_constants
        implicit none
        ! Contains procedures for numerical integration
        ! From Numerical Recipes in FORTRAN
contains
        function simpson(x,y,n)
                implicit none
                integer, parameter :: dp = selected_real_kind(14)
                integer, intent(in) :: n
                integer :: i
                real(kind=dp) :: simpson
                real(kind=dp), dimension(n), intent(in) :: x, y
                real(kind=dp) :: sum, h

                sum = 0.0d0
                h = dabs(x(n) - x(1))/dfloat(n-1)
                do i=1, n/2 ! Since n is odd, n/2 is the same as (n-1)/2
                        sum = sum + y(2*i-1) + 4.0d0*y(2*i) + y(2*i+1)
                end do
                simpson = sum*h/3.0d0
        end function simpson

        subroutine rk4(x,y, dx,n_steps,  n_eq, derivs)
                implicit none
                integer, intent(in) :: n_eq, n_steps
                integer :: i
                real(kind=dp), intent(inout) :: x(n_steps), y(n_steps, n_eq)
                real(kind=dp), intent(in) ::  dx
                external :: derivs
                real(kind=dp) :: xt 
                real(kind=dp), dimension(n_eq) :: dydx, k1, k2, k3, k4, yt

                do i=2, n_steps
                        xt = x(i-1); yt = y(i-1,:)
                        call  derivs(xt,yt,dydx)
                        k1 = dx*dydx
                        call  derivs(xt+0.5d0*dx, yt+0.5d0*k1,dydx)
                        k2 = dx*dydx
                        call  derivs(xt+0.5d0*dx, yt+0.5d0*k2,dydx)
                        k3 = dx*dydx
                        call  derivs(xt+dx,yt+k3,dydx)
                        k4 = dx*dydx

                        x(i) = xt + dx
                        y(i,:) = yt +(0.5d0*(k1+k4) + k2 + k3)/3.0d0
                end do
        end subroutine rk4
end module numerical_integration
      
      


