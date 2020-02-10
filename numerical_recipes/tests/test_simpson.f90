program simpson_test
        use numerical_integration
        implicit none
        integer, parameter :: n = 101, dp = selected_real_kind(14)
        integer :: i
        real(kind=dp), dimension(0:n-1) :: x, y
        real(kind=dp) :: sum

        sum=0.0d0
        do i=0,n-1
                x(i) = dfloat(i)/(n-1)
                y(i) = x(i)
        end do

        sum = simpson(x,y,n)
        write(*,*) sum

end program simpson_test
