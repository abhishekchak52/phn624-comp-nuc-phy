program test
      implicit none
      integer, parameter :: dp = selected_real_kind(14)
      integer :: i
      real(8), dimension(10) :: a = [ (1 + i/10.0d0, i = 1,10) ] 

      do i=0,10
        write(*,*) "hello"
      end do
      write(*,*) a
end program test
