program ortho_test
! Checks for orthogonality of Harmonic Oscillator 
! position-space wavefunctions.       
        use ho_wf
        implicit none

        nmax = 10; num_pts=1000
        x_min = -5.0d0; x_max = 5.0d0
       ! allocate(wf(nmax, num_pts))
        
!        call gen_wf(nmax, x_min, x_max, num_pts, wf)

end program ortho_test
