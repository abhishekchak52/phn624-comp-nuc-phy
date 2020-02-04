program simpson_test
      implicit double precision (a-h,o-z)
      dimension x(100), y(100)
!      parameter n=100
      
      sum=0.0d0
      do i=1,100
        x(i) = dfloat(i)/100
        y(i) = x(i)
        !write(*, *) i, x(i), y(i)
      end do

      sum = simpson(x,y,100)
      write(*,*) sum

end program simpson_test
