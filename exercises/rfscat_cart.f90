program rutherford_scattering
      use numerical_integration
      implicit none

      integer, parameter :: n_pts = 100
      real(kind=dp), parameter :: Mn = 939.5654133, Mz = 938.2720813
      integer :: Zp, Zt, Np, Nt, Ap, At, i, j, k
      real(kind=dp) :: Mp, Mt, Rt, imp_param, KEp, Vp, dt
      real(kind=dp) :: t_list(n_pts), z(n_pts, 4) 

      read(*,*) Zp, Np, Zt, Nt, imp_param

      dt = 100
      Ap = Zp+Np; At = Zt+Nt;
      Mp = Zp*Mz + Np*Mn; Mt = Zt*Mz + Np*Mn;
      Rt = 1.2*At**(1.0d0/3.0d0)
      KEp = 5.0d0
      Vp = sqrt(2*KEp/Mp)

      ! Set inital conditions
      z(1,:) = [-100.0d0, vp, imp_param, 0.0d0]

      call rk4(t_list, z, dt, n_pts, 4, rf_scattering_deriv)

      do i=1, n_pts
        write(*, "(3f14.5)") t_list(i), z(i,1), z(i,3)
      end do


contains
      subroutine rf_scattering_deriv(t, z, dzdt)
              implicit none
              real(kind=dp), intent(in) :: t, z(4)
              real(kind=dp), intent(out) :: dzdt(4)
              real(kind=dp) :: r
              real(kind=dp), parameter :: e2 = 1.44
              integer ::  qQ

              qQ = Zp*Zt
              r = max(sqrt(z(1)*z(1) + z(3)*z(3)), Rt)

              dzdt(1) = z(2); dzdt(3) = z(4);
              dzdt(2) = qQ*e2*z(1)/(Mp*r*r*r)
              dzdt(4) = qQ*e2*z(3)/(Mp*r*r*r)

      end subroutine rf_scattering_deriv
      
end program rutherford_scattering
