module nuclear_models
  use numerical_constants
  implicit none

contains
  subroutine nilsson_model(num_levels, E_levels, delta, max_num_levels)
    implicit none
    integer, intent(in) :: num_levels, max_num_levels
    real(kind=dp), intent(in) :: delta
    real(kind=dp), intent(out) :: E_levels(max_num_levels)
    real(kind=dp) :: ham(max_num_levels, max_num_levels), evals(max_num_levels), evecs(max_num_levels)
    real(kind=dp) :: rkappa = 0.05d0, rmu(0:7) = [0.0d0,0.0d0, 0.0d0, 0.35d0, 0.625d0, 0.63d0, 0.448d0, 0.434d0]
    real(kind=dp) :: fdel
    real(kind=dp) :: hw0, hw00, c, d, h00, hl2, hls, hr2, hy20, hdel
    ! Quantum numbers and basis states calculations
    integer :: n, l, omega, lambda, sigma, l_prime, lambda_prime, sigma_prime, nbas=0, i, j, k, icol
    integer, dimension(num_levels) :: l_list, lambda_list, sigma_list
    real(kind=dp), external :: cgc, factlog
    external :: diag

    ! All quantum numbers multiplied by 2
    fdel=((1.0d0 + (2.0d0/3.0d0)*delta)**2)*(1.0d0-(4.0d0/3.0d0)*delta)**(-1.0d0/6.0d0)
    icol = 0
    do n=0,num_levels*2,2
       do omega=1,n+1, 2
          nbas=0
          do l=n,0,-4
             do lambda=-l,l,2
                sigma=omega-lambda
                if(abs(sigma) .ne. 1) cycle ! Only spin half
                nbas=nbas+1
                l_list(nbas) = l
                lambda_list(nbas) = lambda
                sigma_list(nbas) = sigma
             end do
          end do
          hw0=1.0d0; hw00=hw0*fdel;
          c=-2.0d0*hw00*rkappa
          d=rmu(n/2)*c/2.0d0
          do i=1, nbas
             l=l_list(i)
             lambda = lambda_list(i)
             sigma = sigma_list(i)
             do j = i, nbas
                l_prime = l_list(j)
                lambda_prime = lambda_list(j)
                sigma_prime = sigma_list(j)

                h00 = 0.0d0; hl2=0.0d0; hls=0.0d0; hr2=0.0d0; hy20=0.0d0;

                if(i .eq. j) then
                   h00 = float(n+3)/2.0d0*hw0
                   hl2 = float(l*(l+2))/4.0d0
                end if

                if (l .eq. l_prime) then
                   if ((lambda_prime .eq. lambda+2 ) .and. (sigma_prime .eq. sigma-2)) then
                      hls = sqrt(float((l-lambda)*(l+lambda+2)))/4.0d0
                   else if ((lambda_prime .eq. lambda) .and. (sigma_prime .eq. sigma)) then
                      hls = lambda*sigma/4.0d0
                   else if ((lambda_prime .eq. lambda-2) .and. (sigma_prime .eq. sigma+2)) then
                      hls = sqrt(float((l+lambda)*(l-lambda+2)))/4.0d0
                   end if
                   hr2 = float(n+3)/2.0d0
                else if (l_prime .eq. l-4) then
                   hr2 = sqrt(float((n-l+4)*(n+l+2)))/2.0d0
                else if (l_prime .eq. l+4) then
                   hr2 = sqrt(float((n-l)*(n+l+6)))/2.0d0
                end if

                if ((abs(hr2) .gt. 1e-5) .and. (lambda_prime .eq. lambda)) then
                   hy20 = sqrt(float(l+1)/float(l_prime+1))*&
                   cgc(l,4,l_prime, lambda, 0)*cgc(l, 4, l_prime, 0,0)
                end if
                hdel = -delta*hw0*sqrt(4.0d0)/3.0d0*hr2*hy20

                ham(i,j) = h00+hdel+c*hls + d*hl2
                ham(j,i) = ham(i,j)
             end do
          end do

          call diag(ham, max_num_levels, nbas, evals, evecs)

          do i=1,nbas
             E_levels(icol+i) = evals(i)
          end do
          icol = icol+nbas
       end do
    end do




  end subroutine nilsson_model
end module nuclear_models
