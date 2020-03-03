program rfsc_pol
      use numerical_integration, only: dp, pi
      implicit none

      integer, parameter :: n_pts = 100
      real(kind=dp), parameter :: Mn = 939.5654133, Mz = 938.2720813
      integer :: Zp, Zt, Np, Nt, Ap, At, i, j, k
      real(kind=dp) :: Mp, Mt, Rt, imp_param, max_imp_param, KEp, dtheta
      real(kind=dp) :: theta(n_pts), u(n_pts), kappa, t0, u0, x0, theta_i, theta_f

      read(*,*) Zp, Np, Zt, Nt, imp_param
      x0=-100.0d0

      Ap = Zp+Np; At = Zt+Nt;
      Mp = Zp*Mz + Np*Mn; Mt = Zt*Mz + Np*Mn;
      Rt = 1.2*At**(1.0d0/3.0d0)
      KEp = 5.0d0

      kappa = 1.44d0*Zp*Zt/(2*KEp*imp_param**2)
      t0 = pi/2+atan(imp_param*kappa)
      u0=-kappa/cos(t0)
      theta_i = atan2(imp_param,x0)
      theta_f = 2*t0-pi
      dtheta = (theta_f - theta_i)/(n_pts-1)
      theta(:) = [ (theta_i + i*dtheta, i=0,n_pts-1)]
      u(:) = u0*cos(theta-t0) - kappa
      do i=1,n_pts
         print '(2f10.4)', theta(i), 1.0d0/u(i)
      end do

end program rfsc_pol
