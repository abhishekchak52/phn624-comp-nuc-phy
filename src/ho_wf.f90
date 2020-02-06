module ho_wf
        implicit none
        save ! To retain wavefunction data across subroutines and programs
        !        integer, parameter :: nmax = 100 ! Number of H.O. basis functions
!        integer, parameter :: num_pts = 1000 
        integer, parameter :: dp = selected_real_kind(14)
        integer :: m, n, i, nmax, num_pts
        real(kind=dp), parameter :: pi = real(4, kind=dp)*atan(real(1, kind=dp)), &
                pi_4_root = sqrt(sqrt(pi))
                
        real(kind=dp) :: x_min, x_max, dx
        real(kind=dp), allocatable,  dimension(:,:) :: wf, wf_deriv2

contains
        subroutine gen_wf(nmax, x_min, x_max, num_pts, wf)
                implicit none
                integer, intent(in) :: nmax, num_pts
                integer, parameter :: dp = selected_real_kind(14)
                real(kind=dp), intent(in) :: x_min, x_max
                real(kind=dp) ::  dx, x_pts(num_pts)
                real(kind=dp), dimension(num_pts, nmax), intent(inout) :: wf
                
                dx = (x_max-x_min)/(num_pts-1)
                x_pts = [ (x_min+i*dx, i=0, num_pts-1) ]

                write(*, *) x_pts

        end subroutine gen_wf        

end module ho_wf
