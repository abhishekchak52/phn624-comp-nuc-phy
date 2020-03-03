program ho_wf
        use numerical_integration
        implicit none
        integer, parameter :: nmax = 10, num_pts = 101 ! Number of H.O. basis functions
        integer :: m, n, i, j
        real(kind=dp), parameter :: pi_4_root = sqrt(sqrt(pi))
        real(kind=dp) :: x_min, x_max, dx, x_pts(num_pts)
        real(kind=dp), dimension(0:nmax-1, 0:nmax-1) :: inner_prod, pot_e, kin_e, ham
        real(kind=dp), dimension(num_pts, 0:nmax-1) :: wf, wf_deriv2

        x_min = -10.0d0; x_max = 10.0d0
        dx = (x_max-x_min)/(num_pts-1)
        x_pts = [ (x_min+i*dx, i=0, num_pts-1) ]


        call gen_wf(x_pts, num_pts, nmax, wf, wf_deriv2)

        ! Calculate second derivative

        ! Print Wavefunctions upto order N
        open(unit=50, file='howf_data.txt')
        do i=1,num_pts
           write(50, *) x_pts(i), wf(i, :)
        end do
        close(50)

        ! Calculate orthogonality
        do i=0,nmax-1
           do j=0,nmax-1
                inner_prod(i,j) = simpson(x_pts, wf(:,i)*wf(:,j), num_pts)
                pot_e(i, j) = simpson(x_pts, wf(:,i)*wf(:,j)*pot_func(x_pts,num_pts), num_pts)
                kin_e(i,j) = simpson(x_pts, -wf(:,i)*wf_deriv2(:,j), num_pts)
           end do
        end do

        ham = kin_e + pot_e

        open(unit=50, file='ortho_data.txt')
        open(unit=51, file='pot_e_data.txt')
        open(unit=52, file='kin_e_data.txt')
        open(unit=53, file='ham_data.txt')
        do i=0,nmax-1
                write(50,*) inner_prod(i,:)
                write(51,*) pot_e(i,:)
                write(52,*) kin_e(i,:)
                write(53,*) ham(i,:)
        end do
        close(50)
        close(51)
        close(52)
        close(53)
contains
        function pot_func(x, n)
                integer, intent(in) :: n
                real(kind=dp) :: pot_func(n)
                real(kind=dp), intent(in) :: x(n)

                pot_func = x**2
        end function pot_func

        subroutine gen_wf(x_pts, num_pts, max_order, wf, wf_deriv2)
          ! Generates harmonic oscillator wavefunctions upto order max_order
          ! Wavefunction order starts from 0
                implicit none
                integer, intent(in) :: max_order, num_pts
                real(kind=dp), intent(in) :: x_pts(num_pts)
                real(kind=dp), intent(inout), dimension(num_pts, 0:max_order-1) :: wf, wf_deriv2
                integer :: i, j, k

                wf(:,0) = exp(-0.5d0*x_pts**2)/pi_4_root
                wf_deriv2(:,0) = (x_pts**2 - 1)*wf(:,0)
                if (nmax > 1) then
                        wf(:,1) = sqrt(2.0d0)*x_pts*wf(:,0)
                        wf_deriv2(:,1) = (x_pts**2 - 1)*wf(:,1) - sqrt(8.0d0)*x_pts*wf(:,0)
                        do i=2, max_order-1
                           wf(:,i) = sqrt(2.0d0/i)*x_pts*wf(:,i-1) - sqrt(float(i-1)/i)*wf(:,i-2)
                           wf_deriv2(:,i) = 2*sqrt(float(i*(i-1)))*wf(:,i-2) - sqrt(8.0d0*i)*x_pts*wf(:,i-1) &
                                + (x_pts**2 - 1)*wf(:,i)
                        end do
                end if

        end subroutine gen_wf        

end program ho_wf
