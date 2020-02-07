      implicit double precision (a-h,o-z)
      complex*8 ylm
      character*1 aal,aam
      common /alegp/ alp(0:100,-100:100)
      common /factl/ fact(0:100)
      dimension ylm(0:100,0:100)
      pi=4.0d0*datan(1.0d0)
      read(*,*)lmax
      call factorial(2*lmax)
      ifno=0
      do l=0,lmax
        write(aal,'(i1)')l 
        do m=0,l
          write(aam,'(i1)')m 
          ifno=ifno+1
          open(ifno,file='y'//aal//aam//'r.out',status='unknown')
          ifno=ifno+1
          open(ifno,file='y'//aal//aam//'i.out',status='unknown')
        end do
      end do
      theta=0.d0
      ainc=pi/10.d0
      do while (theta.le.pi)
        x=dcos(theta)
        call alegendrep(lmax,x)
        phi=0.d0
        do while (phi.le.2.0*pi)
          xx=dsin(theta)*dcos(phi)
          yy=dsin(theta)*dsin(phi)
          ifno=0
          do l=0,lmax
            do m=0,l
              fac=dsqrt((2*l+1)*fact(l-m)/(4.0d0*pi*fact(l+m)))
              ylm(l,m)=fac*alp(l,m)*dcmplx(dcos(m*phi),dsin(m*phi))
              ifno=ifno+1
              rr=dreal(ylm(l,m))
              write(ifno,'(1x,3f12.6)') rr*xx,rr*yy,rr*dcos(theta)
              ifno=ifno+1
              rr=imag(ylm(l,m))
              write(ifno,'(1x,3f12.6)') rr*xx,rr*yy,rr*dcos(theta)
            end do
          end do
          phi=phi+ainc
        end do
        theta=theta+ainc
      end do

      end

c Calculates assosicated Legendre polynomials upto lmax for a given x
c -1 < x < 1; Condon-Shortley phase is included
c i/p lmax, x
c o/p alp(0:lmax,-lmax:lmax)
      subroutine alegendrep(lmax,x)
      implicit real*8 (a-h,o-z)
      common /alegp/ alp(0:100,-100:100)
      common /factl/ fact(0:100)
      alp(0,0)=1.0d0
      do l=1,lmax
        alp(l,l-1)=x*(2.0d0*(l-1.0d0)+1.0d0)*alp(l-1,l-1)
        alp(l,-l+1)=(-1)**(l-1)*alp(l,l-1)/fact(l+l-1)  !m=l-1
        alp(l,l)=-(2.0d0*(l-1.0d0)+1.0d0)*dsqrt(1.d0-x*x)*alp(l-1,l-1)
        alp(l,-l)=(-1)**l*alp(l,l)/fact(l+l)
        do m=l-2,1,-1
          alp(l,m)=(x*(2*l-1)*alp(l-1,m)-(l+m-1)*alp(l-2,m))/dfloat(l-m)
          alp(l,-m)=(-1)**m*fact(l-m)*alp(l,m)/fact(l+m)
        end do
        alp(l,0)=(x*(2*l-1)*alp(l-1,0)-(l-1)*alp(l-2,0))/dfloat(l)
      end do
      end


c Calculates and stores 0!,1!,2!,...,n!
c i/p n
c o/p fact(0:n)
      subroutine factorial(n)
      implicit double precision (a-h,o-z)
      common /factl/ fact(0:100)
      fact(0)=1
      do i=1,n
        fact(i)=fact(i-1)*i
      end do
      end
