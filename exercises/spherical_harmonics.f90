program spherical_harmonics
  use special_functions
  integer, parameter :: num_theta_pts = 100, num_phi_pts = 100, lmax=4
  real(kind=dp) :: alp(num_theta_pts, 0:lmax,-lmax:lmax), theta_pts(num_theta_pts), dtheta, dphi
  real(kind=dp) :: fac, phi
  integer :: i, j, k, l, m, fact(2*lmax)
  complex(kind=dp) :: spha_ylm(num_theta_pts, num_phi_pts)


  dtheta = pi/(num_theta_pts-1)
  dphi = 2*pi/(num_phi_pts-1)
  theta_pts = [ (-pi/2+i*dtheta, i=0, num_theta_pts-1) ]
  l=2
  m=0


  call calc_factorials(2*lmax, fact)
  call gen_assoc_legendre_pol(cos(theta_pts), num_theta_pts, lmax, alp)

  do i=0,num_phi_pts-1
    phi = i*dphi
    fac=dsqrt((2*l+1)*fact(l-m)/(4.0d0*pi*fact(l+m)))
    spha_ylm(:,i)=fac*alp(:,l,m)*cmplx(cos(m*phi),sin(m*phi))
  end do



  do i=1,num_theta_pts
     do j=1,num_phi_pts
    print '(5f10.4)', theta_pts(i), (j-1)*dphi,  spha_ylm(i,j)%re, spha_ylm(i,j)%im, abs(spha_ylm(i, j))
     end do
  end do
end program spherical_harmonics
