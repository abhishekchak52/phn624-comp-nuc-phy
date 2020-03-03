program swe_inf_pot_well
      use numerical_integration
      implicit none

      integer, parameter :: n_pts = 10000, n_energy = 40, max_iter_bisec = 100
      real(kind=dp), parameter :: l0 = 0.14871085d0, E0 = 938.27208816d0
      real(kind=dp) :: x(n_pts), psi(n_pts,2) ! Contains psi and psi-dot
      real(kind=dp) :: E, E_high, E_low, xmax, xmin, dx, E_list(n_energy)
      real(kind=dp) :: E_levels(5), wf_list(5,n_pts), wf_norm
      integer :: i, j, n_zero_low, n_zero_high, n_zero_mid, n_evals
      ! length and energy scales (fm and MeV)
      xmax = 12.0d0/l0
      xmin = -12.0d0/l0
      dx = (xmax-xmin)/(n_pts-1)
      x(1) = xmin
      psi(1,:) = [0.0d0, 0.01d0]
      E_list = [(i*0.001d0, i=0,n_energy-1)]

      ! Count energy levels
!      E = E_list(n_energy)
!      call rk4(x, psi, dx, n_pts, 2, swe_deriv)
!      n_evals = n_zeros(psi, n_pts)

!      write(*,*) E_list
      n_evals = 1
      do i=1,n_energy
              if (n_evals > 5) exit;
              E_low = E_list(i); E_high = E_list(i+1);

              E = E_low
              call rk4(x, psi, dx, n_pts, 2, swe_deriv)
              n_zero_low = n_zeros(psi, n_pts)

!              writ(*,*) E_low, n_zero_low
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
!                      E_levels(i) = E
                      wf_list(n_evals, :) = psi(:,1)
!                      write(*, *) E
                      E_levels(n_evals) = E
                      n_evals = n_evals + 1
              end if
              !do i=1,n_pts
              !        write(*,*) x(i), psi(i,1), psi(i,2)
              !end do  
      end do
      write (*,*) E_levels*E0
      ! Normalize the wavefunctions
      do i=1, 5
         wf_norm = 0.0d0
         do j=1,n_pts
            wf_norm = wf_norm + wf_list(i, j)*wf_list(i, j)*dx
         end do
         wf_list(i, :) = wf_list(i, :)/sqrt(wf_norm)
      end do

      open(unit=25, file='wf_data.txt')
      write(25, "(A,T12,5f11.4)") 'x(t)', E_levels*E0
      do i=1, n_pts
        write(25,"(6f11.4)") x(i)*l0, wf_list(:, i)
      end do
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

              ! Infinite Well
              ! pot = 0.0d0

              ! Finite Well
              ! if (abs(x) > 10.0d0/l0) then
              !         pot = 25.0d0/E0
              ! else
              !         pot = 0.0d0
              ! end if

              ! Woods-Saxon(Correct factors not yet calculated)
              pot = 25.0d0*(1.0d0-1.0d0/(1.0d0+exp(abs(2*l0*x)-10)))/E0

              dpsidx(1) = psi(2)
              dpsidx(2) = -1*(E-pot)*psi(1)

        end subroutine swe_deriv

end program swe_inf_pot_well
