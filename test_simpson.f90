program simpson_test
      implicit none
      integer, parameter :: n = 100, dp = selected_real_kind(14)
      integer :: i
      real(kind=dp), dimension(0:n-1) :: x, y
      real(kind=dp) :: sum

      interface
        function simpson(x,y,n)
                integer, parameter :: dp = selected_real_kind(14)
                integer, intent(in) :: n
                real(kind=dp) :: simpson
                real(kind=dp), dimension(n), intent(in) :: x, y
        end function simpson
      end interface


      
      sum=0.0d0
      do i=1,100
        x(i) = dfloat(i)/100
        y(i) = x(i)
!        write(*, '(i, 2f10.3)') i, x(i), y(i)
      end do

      sum = simpson(x,y,100)
      write(*,*) sum

end program simpson_test
