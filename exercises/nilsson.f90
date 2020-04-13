program nilsson
  use numerical_constants
  use nuclear_models
  integer, parameter :: max_levels = 100
  real(kind=dp) :: spec(max_levels), delta = 0.0d0


  open(unit=25, file="nilsson_data.txt")
  do delta = -0.1d0,0.101d0,0.01d0
     call nilsson_model(7, spec, delta, max_levels)

     write(25,*) delta, spec
  end do
  close(25)

end program nilsson
