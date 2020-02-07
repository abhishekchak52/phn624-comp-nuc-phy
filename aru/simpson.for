      real*8 function simp(x,y,n)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      h=dabs(x(n)-x(1))/dfloat(n-1)
      sum=0.0d0
      do i=1,n/2
         sum=sum+y(2*i-1)+4.0d0*y(2*i)+y(2*i+1)
c         write(*,*)i,x(i),y(i)
      end do
      simp=sum*h/3.0d0
      end

