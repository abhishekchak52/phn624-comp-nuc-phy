program complex_test
  use numerical_constants
  implicit none

  complex(kind=dp) :: trig, polar
  real :: a, b(10)

  b(:) = 2.0d0

  print *, b

  !trig = (cos(pi), sin(pi))
  ! write(*,*) a, b, pi, trig, polar

end program complex_test
