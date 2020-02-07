      implicit real*8 (a-h,o-z)
      parameter(n=15)
      dimension xi(n),yi(n),dy(n),dy1(n)
      pi=4.d0*datan(1.d0)
      xmin=-pi
      xmax= pi
      dx=(xmax-xmin)/dfloat(n-1)
      do i=1,n
        xi(i)=xmin+float(i-1)*dx
        yi(i)=sin(xi(i))
      end do
      read *, ndiv
      call nderiv(xi,yi,n,n,ndiv,dy)
      call splderiv(xi,yi,n,n,dy1)
      do i=1,n
        print '(1x,i5,8f12.6)', i,xi(i),yi(i),dy(i)
     &          ,dabs(cos(xi(i))-dy(i))
     &          ,dabs(cos(xi(i))-dy1(i))
      end do
      end

      include 'nderiv.f'
