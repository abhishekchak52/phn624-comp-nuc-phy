program swe_inf_pot_well
      use numerical_integration
      implicit none

      integer, parameter :: n_pts = 1000, n_energy = 40, max_iter_bisec = 50
      real(kind=dp) :: x(n_pts), psi(n_pts,2) ! Contains psi and psi-dot
      real(kind=dp) :: E, E_high, E_low, xmax, xmin, dx, E_list(n_energy)
      real(kind=dp), allocatable :: E_levels(:), wf_list(:,:)
      integer :: i, j, n_zero_low, n_zero_high, n_zero_mid, n_evals
      
      xmax = 10
      xmin = -10
      dx = (xmax-xmin)/(n_pts-1)
      x(1) = xmin
      psi(1,:) = [0.0d0, 0.01d0]
      E_list = [(i*0.025d0, i=1,n_energy)]

      ! Count energy levels
      E = E_list(n_energy)
      call rk4(x, psi, dx, n_pts, 2, swe_deriv)
      n_evals = n_zeros(psi, n_pts)

      allocate(E_levels(n_evals))
      allocate(wf_list(n_evals, n_pts))
!      write(*,*) E_list
      do i=1,n_energy
              E_low = E_list(i); E_high = E_list(i+1);

              E = E_low
              call rk4(x, psi, dx, n_pts, 2, swe_deriv)
              n_zero_low = n_zeros(psi, n_pts)

!              write(*,*) E_low, n_zero_low
              E = E_high
              call rk4(x, psi, dx, n_pts, 2, swe_deriv)
              n_zero_high = n_zeros(psi, n_pts)

              if (n_zero_high == n_zero_low + 1) then
                      do j=1,max_iter_bisec
                        E = 0.5d0*(E_low + E_high)
                        call rk4(x, psi, dx, n_pts, 2, swe_deriv)
                        n_zero_mid = n_zeros(psi, n_pts)

                        if (n_zero_mid == n_zero_high) then
                                E_high = E
                        else
                                E_low = E
                        end if 
                      end do
                      E_levels(i) = E
                      wf_list(i,:) = psi(:,1)
!                      write(*, *) E*938.272
              end if
              !do i=1,n_pts
              !        write(*,*) x(i), psi(i,1), psi(i,2)
              !end do  
      end do
      write (*,*) E_levels
!      do i=1, n_pts
!        write(*,*) x(i), wf_list(:, i)
!      end do
!      write(*,*) n_zeros(psi, n_pts)
contains

        function n_zeros(wf, n_steps)
                implicit none
                integer, intent(in) :: n_steps
                real(kind=dp), dimension(n_steps), intent(in) :: wf
                integer :: i,n_zeros

                n_zeros = 0


                do i=2,n_steps
                        if (wf(i-1)*wf(i) < 0) n_zeros = n_zeros + 1
                end do

        end function n_zeros


        subroutine swe_deriv(x, psi, dpsidx)
              implicit none
              real(kind=dp) :: pot
              real(kind=dp), intent(in)  :: x, psi(2)
              real(kind=dp), dimension(2), intent(out) :: dpsidx

              pot  = 0
              dpsidx(1) = psi(2)
              dpsidx(2) = -(E-pot)*psi(1)

        end subroutine swe_deriv

end program swe_inf_pot_well
