
program multipole_moment
  use special_functions
  use numerical_integration
  integer, parameter :: num_theta_pts = 101, num_phi_pts = 101, num_r_pts = 101, lmax=4
  real(kind=dp) :: alp(num_theta_pts, 0:lmax,-lmax:lmax), theta_pts(num_theta_pts), dtheta, dphi
  real(kind=dp) :: fac, phi, r_pts(num_r_pts), rho(num_r_pts), r_cutoff, dr, phi_pts(num_phi_pts)
  real(kind=dp) :: re_integral, im_integ, re_phi_integ(num_r_pts, num_theta_pts), re_theta_integ(num_r_pts)
  integer :: i, j, k, l, m, fact(2*lmax)
  complex(kind=dp) :: spha_ylm(num_theta_pts, num_phi_pts), integrand(num_r_pts, num_theta_pts, num_phi_pts)
  complex(kind=dp) :: integral


  dtheta = pi/(num_theta_pts-1)
  dphi = 2*pi/(num_phi_pts-1)
  theta_pts = [ (-pi/2 + i*dtheta, i=0, num_theta_pts-1) ]
  r_cutoff = 10.0d0
  dr = r_cutoff/(num_r_pts-1)
  r_pts = [ (i*dr, i=0, num_r_pts-1) ]
  phi_pts = [(i*dphi, i=0, num_phi_pts-1)]
  rho_pts= 1.0d0
  l=2
  m=0


  call calc_factorials(2*lmax, fact)
  call gen_assoc_legendre_pol(cos(theta_pts), num_theta_pts, lmax, alp)

  do i=0,num_phi_pts-1
    phi = i*dphi
    fac=dsqrt((2*l+1)*fact(l-m)/(4.0d0*pi*fact(l+m)))
    spha_ylm(:,i)=fac*alp(:,l,m)*cmplx(cos(m*phi),sin(m*phi))
  end do

  do i=1,num_r_pts
     do j=1,num_theta_pts
        do k=1,num_phi_pts
            integrand(i, j, k) = r_pts(i)**(l+2)*conjg(spha_ylm(j,k))*rho(i)*sin(theta_pts(j))*dr*dtheta*dphi
        end do
     end do
  end do


  do i=1,num_r_pts
    do j=1,num_theta_pts
        re_phi_integ(i,j) = simpson(phi_pts, real(integrand(i,j,:)), num_phi_pts)
    end do
    re_theta_integ(i) = simpson(theta_pts, re_phi_integ(i, :), num_theta_pts)
  end do
  re_integral = simpson(r_pts, re_theta_integ(:), num_r_pts)

  print*, re_integral

  ! do i=1,num_theta_pts
  !    do j=1,num_phi_pts
  !   print '(5f10.4)', theta_pts(i), (j-1)*dphi,  spha_ylm(i,j)%re, spha_ylm(i,j)%im, abs(spha_ylm(i, j))
  !    end do
  ! end do
end program multipole_moment
