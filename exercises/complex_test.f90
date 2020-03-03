program complex_test
  implicit none

  complex :: trig, polar
  real :: a, pi, b(10)

  pi = 4.0d0*atan(1.0d0)
  a = 2.0d0

  b = [ (a, a=1.0d0,4.0d0,0.5d0) ]

  !trig = (cos(pi), sin(pi))
  ! write(*,*) a, b, pi, trig, polar

end program complex_test
