       implicit double precision (a-h,o-z)
       common /factl/ fact(0:100)
       common /herm/ h(0:100)
       dimension wf(0:100)
       pi=4.0d0*datan(1.0d0)
       read(*,*)nmax
       call factorial(nmax)
       nx=100
       xmin=-5.d0
       xmax=+5.d0
       xinc=(xmax-xmin)/nx
       x=xmin
       do i=1,nx+1
          call hermite(nmax,x)
          do n=0,nmax
            wf(n)=dsqrt(1.0d0/(2.0d0**n*dsqrt(pi)*fact(n)))*
     &         dexp(-x*x/2.0d0)*h(n) 
          end do
          write(*,'(1x,f5.2,101f12.6)')x,(wf(n),n=0,nmax)
          x=x+xinc
       end do

       end

       subroutine factorial(n)
       implicit double precision (a-h,o-z)
       common /factl/ fact(0:100)
       fact(0)=1
       do i=1,n
         fact(i)=fact(i-1)*i
       end do
       end

      subroutine hermite(n,x)
      implicit real*8 (a-h,o-z)
      common /herm/ h(0:100)
      h(0)=1.0d0
      h(1)=2.0d0*x
      do i=2,n
        h(i)=2.0d0*x*h(i-1)-2.0d0*(i-1)*h(i-2)
      end do
      end

   
