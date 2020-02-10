module numerical_integration
        implicit none
        ! Contains procedures for numerical integration

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
end module numerical_integration
      
      


