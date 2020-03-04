program assoc_leg_pol_test
  use special_functions
  integer, parameter :: num_pts = 100, lmax=4
  real(kind=dp) :: alp(num_pts, 0:lmax,-lmax:lmax), x_pts(num_pts)
  integer :: i, j, k
  character :: m

  x_pts = [ (-1+2*float(i)/(num_pts-1), i=0, 99) ]

  call gen_assoc_legendre_pol(x_pts, num_pts, lmax, alp)

  do i=1,num_pts
    print '(6f10.4)', x_pts(i), alp(i, :, 1)
  end do
end program assoc_leg_pol_test
