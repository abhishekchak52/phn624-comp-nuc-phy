      subroutine splderiv(x,y,nd,n,u)
c Calculates the numerical derivative using spline interpolation
c Calls the Numerical Recipes routine TRIDIAG (double precision)
c REMARK
c End point behaviour is not so good
c INPUTS
c x(n),y(n) - arrays comprising x & y data (equispaced x)
c n - array size
c OUTPUT
c u(n) - first order derivative (dy/dx)
c
c P. Arumugam, IIT Roorkee - 07-June-2010
c
      implicit double precision (a-h,o-z)
      dimension x(nd),y(nd)
      dimension a(n),b(n),c(n),r(n),u(n)
      h=(x(n)-x(1))/dfloat(n-1)
      a(n)=1.d0
      b(1)=2.d0
      b(n)=2.d0
      c(n-1)=1.d0
      r(1)=3.d0*(y(2)-y(1))/h
      r(n)=3.d0*(y(n)-y(n-1))/h
      do i=2,n-1
        a(i)=1.d0
        b(i)=4.d0
        c(i-1)=1.d0
        r(i)=3.d0*(y(i+1)-y(i-1))/h
      end do
      call tridag(a,b,c,r,u,n)
      end
