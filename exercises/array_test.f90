program array_test
  implicit none
  real :: x(21)
  integer :: i, j, k

  x = [(i*0.1, i=0,20)]

  do i=2,4
     print*, conjg(cmplx(i))
  end do

end program array_test
