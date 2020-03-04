module numerical_constants
  implicit none

        integer, parameter :: dp = selected_real_kind(14)
        real(kind=dp), parameter :: pi = real(4, kind=dp)*atan(real(1, kind=dp))

end module numerical_constants
