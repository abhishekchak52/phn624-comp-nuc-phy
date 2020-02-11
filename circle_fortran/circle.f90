program circle
implicit none
integer, parameter :: num_pts = 100 
integer :: i
real, parameter :: pi = 4.0*atan(1.0)
real :: x,y, theta

do i=1, num_pts
        x = cos(i*2*pi/num_pts)
        y = sin(i*2*pi/num_pts)
        write(*,*) x, y, atan2(y,x)
        
end do
end program circle




