program rutherford_scattering
        implicit none

        integer, parameter :: dp = selected_real_kind(14)
        real(kind=dp), parameter :: pi = real(4, kind=dp)*atan(real(1, kind=dp))
        integer, parameter :: n_pts = 100
        real(kind=dp), parameter :: Mn = 939.5654133, Mz = 938.2720813
        integer :: Zp, Zt, Np, Nt, Ap, At, i, j, k
        real(kind=dp) :: Mp, Mt, Rt, imp_param, KEp, Vp, dt
        real(kind=dp) :: t_list(n_pts), z(n_pts, 4) 
        ! Bisection variables
        real(kind=dp) :: R_hi, R_lo, R_mid, Rt_th

        read(*,*) KEp
        imp_param = 0
        Zp=2;Np=2;Zt=79;Nt=118

        dt = 100
        Ap = Zp+Np; At = Zt+Nt;
        Mp = Zp*Mz + Np*Mn; Mt = Zt*Mz + Np*Mn;
        Vp = sqrt(2*KEp/Mp)

        Rt_th = 1.2*At**(1.0d0/3.0d0) !Theoretical value
        ! Set inital conditions
        z(1,:) = [-100.0d0, vp, imp_param, 0.0d0]


        ! Setting initial bounds for bisection
        R_hi = 100.0d0*Rt_th; R_lo=0.01d0*Rt_th;

        ! Bisection search for max radius
        do i=1,100
                R_mid = 0.5d0*(R_hi + R_lo)
        !        print *, i, R_mid
                Rt = R_mid 
                if( abs(R_hi - R_lo) < 0.001d0) exit
                call rk4(t_list, z, dt, n_pts, 4, rf_scattering_deriv)
                if (maxval(z(:,1)) > -1.0d0*Rt) then
                        R_hi = R_mid
                else 
                        R_lo = R_mid
                end if
        end do

        write(*, '(3f7.1,i3)') Rt, maxval(z(:,1)), Rt_th, i
!        do i=1, n_pts
!        write(*, "(3f14.5)") t_list(i), z(i,1), z(i,3)
!        end do


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
      


        subroutine rk4(x,y, dx,n_steps,  n_eq, derivs)
                implicit none
                integer, intent(in) :: n_eq, n_steps
                integer :: i
                real(kind=dp), intent(inout) :: x(n_steps), y(n_steps, n_eq)
                real(kind=dp), intent(in) ::  dx
                external :: derivs
                real(kind=dp) :: xt 
                real(kind=dp), dimension(n_eq) :: dydx, k1, k2, k3, k4, yt

                do i=2, n_steps
                        xt = x(i-1); yt = y(i-1,:)
                        call  derivs(xt,yt,dydx)
                        k1 = dx*dydx
                        call  derivs(xt+0.5d0*dx, yt+0.5d0*k1,dydx)
                        k2 = dx*dydx
                        call  derivs(xt+0.5d0*dx, yt+0.5d0*k2,dydx)
                        k3 = dx*dydx
                        call  derivs(xt+dx,yt+k3,dydx)
                        k4 = dx*dydx

                        x(i) = xt + dx
                        y(i,:) = yt +(0.5d0*(k1+k4) + k2 + k3)/3.0d0
                end do
        end subroutine rk4

end program rutherford_scattering
