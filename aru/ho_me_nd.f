       implicit double precision (a-h,o-z)
       common /factl/ fact(0:100)
       common /herm/ h(0:100)
       dimension wf(0:100,1000),ddwf(0:100,1000)
       dimension xx(1000),yy(1000)
       dimension dy(1000),ddy(1000)
       dimension ham(200,200)
       pi=4.0d0*datan(1.0d0)
       read(*,*)nmax
       call factorial(nmax)
       nx=301
       xmin=-15.d0
       xmax=+15.d0
       xinc=(xmax-xmin)/dfloat(nx-1)
       x=xmin
       do i=1,nx
         xx(i)=x
         call hermite(nmax,x)
         do n=0,nmax
           wf(n,i)=dsqrt(1.0d0/(2.0d0**n*dsqrt(pi)*fact(n)))*
     &       dexp(-x*x/2.0d0)*h(n) 
         end do
         ddwf(0,i)=(x**2-1.0d0)*wf(0,i) 
         ddwf(1,i)=-dsqrt(8.0d0)*x*wf(0,i)+(x**2-1.0d0)*wf(1,i) 
         do n=2,nmax
           ddwf(n,i)=2*dsqrt(dfloat(n*(n-1)))*wf(n-2,i)
     &             -dsqrt(dfloat(8*n))*x*wf(n-1,i)
     &             +(x**2-1.0d0)*wf(n,i)
         end do
         x=x+xinc
       end do
       do n=0,nmax
         do i=1,nx
           yy(i)=wf(n,i)
         end do
         call splderiv(xx,yy,1000,nx,dy)
         call splderiv(xx,dy,1000,nx,ddy)
         do i=1,nx
c           print *, i,ddwf(n,i),ddy(i)
c           ddwf(n,i)=ddy(i)
         end do
       end do
       do n=0,nmax
         do n1=0,nmax
           do i=1,nx
             yy(i)=wf(n1,i)*xx(i)**2/2.0d0*wf(n,i)
           end do
           ppe=simp(xx,yy,nx)
           do i=1,nx
             yy(i)=-wf(n1,i)*ddwf(n,i)/2.0d0
c             write(*,*)i,xx(i),yy(i)
           end do
           rke=simp(xx,yy,nx)
           write(*,'(1x,2i2,2f8.4)')n,n1,ppe,rke
           ham(n+1,n1+1)=rke+ppe
         end do
       end do
       do n=0,nmax
         write(*,'(1x,200f8.4)')(ham(n+1,n1+1),n1=0,nmax)
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

   
      real*8 function simp(x,y,n)
      implicit real*8 (a-h,o-z)
      dimension x(1000),y(1000)
      h=dabs(x(n)-x(1))/dfloat(n-1)
      sum=0.0d0
      do i=1,n/2
         sum=sum+y(2*i-1)+4.0d0*y(2*i)+y(2*i+1)
      end do
      simp=sum*h/3.0d0
c      write(*,'(1x,i5,6f12.6)')n,x(1),y(1),x(n),y(n),h,sum
      end
