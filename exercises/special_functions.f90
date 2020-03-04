module special_functions
  use numerical_constants
  implicit none

contains
  subroutine gen_assoc_legendre_pol(x_pts, n_pts, lmax, alp)
    implicit none
    integer, intent(in) :: n_pts, lmax
    real(kind=dp), intent(in) :: x_pts(n_pts)
    real(kind=dp), intent(out) :: alp(n_pts, 0:lmax, -lmax:lmax)
    integer :: l, m, fact(lmax)

    call calc_factorials(lmax, fact)

    alp(:, 0, 0) = 1.0d0
    do l=1,lmax
        alp(:,l,l-1)=x_pts*(2.0d0*(l-1.0d0)+1.0d0)*alp(:,l-1,l-1)
        alp(:,l,-l+1)=(-1)**(l-1)*alp(:,l,l-1)/fact(l+l-1)  !m=l-1
        alp(:,l,l)=-(2.0d0*(l-1.0d0)+1.0d0)*sqrt(1.d0-x_pts*x_pts)*alp(:,l-1,l-1)
        alp(:,l,-l)=(-1)**l*alp(:,l,l)/fact(l+l)
        do m=l-2,1,-1
            alp(:,l,m)=(x_pts*(2*l-1)*alp(:,l-1,m)-(l+m-1)*alp(:,l-2,m))/float(l-m)
            alp(:,l,-m)=(-1)**m*fact(l-m)*alp(:,l,m)/fact(l+m)
        end do
        alp(:,l,0)=(x_pts*(2*l-1)*alp(:,l-1,0)-(l-1)*alp(:,l-2,0))/float(l)
    end do
  end subroutine gen_assoc_legendre_pol

  subroutine calc_factorials(nmax, fact)
    implicit none
    integer, intent(in) :: nmax
    integer, intent(out) :: fact(nmax)
    integer :: i

    fact(0) = 1
    do i=1,nmax
       fact(i) = fact(i-1)*i
    end do

  end subroutine calc_factorials
end module special_functions
