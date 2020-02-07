      subroutine nderiv(x,y,nd,n,ndiv,dy)
* nd -- dimension of x,y,dy.  n--size,
* ndiv -- no. of divisions between succesive x 
*         used for interpolation (ndiv~20)
      implicit real*8 (a-h,o-z)
      dimension x(nd),y(nd),dy(nd)
      dimension x2(n*ndiv),y2(n*ndiv)
      dx=(x(n)-x(1))/dfloat(n-1)
      dx1=dx/dfloat(ndiv)
      nres=n*ndiv
      x2(1)=x(1)-dx1
      y2(1)=polyint(x2(1),x,y,nd,n)
      do i=2,nres
        x2(i)=x2(i-1)+dx1
        y2(i)=polyint(x2(i),x,y,nd,n)
c        print '(1x,i5,5f12.6)', i,x2(i),y2(i),sin(x2(i))
      end do
      j=1
      do i=1,nres
        if(x2(i).ge.x(j))then
          dy(j)=(y2(i+1)-y2(i-1))/(x2(i+1)-x2(i-1))
          j=j+1
          if(j.eq.n+1)exit
        end if
      end do
      end

      real*8 function polyint(x,xi,yi,ndi,ni)
      implicit real*8 (a-h,o-z)
      dimension xi(ndi),yi(ndi)
      polyint=0.d0
      do i=1,ni
        prod=1.d0
        do j=1,ni
          if (i.eq.j) cycle
          prod=prod*(x-xi(j))/(xi(i)-xi(j))
        end do
        polyint=polyint+yi(i)*prod
      end do
      end



