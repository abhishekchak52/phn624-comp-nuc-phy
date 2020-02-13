program test_simpson_3d
        use numerical_integration
        implicit none
        integer, parameter :: dp = selected_real_kind(14)
        integer, parameter :: num_pts = 2 
        integer :: i, j, k

        real(kind=dp), dimension(num_pts) :: x, y, z       
        real(kind=dp) :: fvals(num_pts,num_pts,num_pts), z1(num_pts,num_pts), &
z2(num_pts), z3

        x =  [ (float(i)/num_pts, i=1,num_pts) ] 
        y =  [ (float(i)/num_pts, i=1,num_pts) ] 
        z =  [ (float(i)/num_pts, i=1,num_pts) ] 

        do i=1, num_pts
                do j=1,num_pts
                        z1(i, j) =  f(x(i), y(j), z)
                end do 
        end do

!        write(*,*) f(x, y, z)

contains

       elemental function f(x,y,z)
               implicit none 
               real(kind=dp), intent(in) ::  x,y,z
               real(kind=dp) :: r, f
               

               r = sqrt(x*x + y*y + z*z)
               f = r

       end function f


end program test_simpson_3d
