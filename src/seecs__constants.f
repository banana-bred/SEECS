! ================================================================================================================================ !
module seecs__constants

  use seecs__kinds, only: dp

  implicit none (type, external)

  private

  real(dp),    parameter, public :: zero  = 0.0_dp
  real(dp),    parameter, public :: one   = 1.0_dp
  real(dp),    parameter, public :: two   = 2.0_dp
  real(dp),    parameter, public :: three = 3.0_dp
  real(dp),    parameter, public :: five  = 5.0_dp
  real(dp),    parameter, public :: six   = 6.0_dp
  real(dp),    parameter, public :: pi = acos(-1.0_dp)
  complex(dp), parameter, public :: im = (zero, one)
    !! $\sqrt{-1}$
  character(13), parameter, public :: numeric = "0123456789.+-"
    !! characters considered "numeric"

  real(dp), parameter, public :: macheps = epsilon(1._dp)
    !! Machine epsilon

  ! -- ASCII constants
  integer, parameter, public :: uppercase_a = ichar('A')
  integer, parameter, public :: uppercase_z = ichar('Z')
  integer, parameter, public :: lowercase_a = ichar('a')
  integer, parameter, public :: lowercase_z = ichar('z')

! ================================================================================================================================ !
end module seecs__constants
! ================================================================================================================================ !
