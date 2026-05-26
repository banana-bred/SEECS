! ================================================================================================================================ !
module seecs__utils

  implicit none (type, external)

  private

  public :: operator(.aprx.)

  interface operator(.aprx.)
    module procedure :: is_approx_real
  end interface operator(.aprx.)

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

 ! ------------------------------------------------------------------------------------------------------------------------------- !
 pure module function is_approx_real(a, b) result(res)
   !! Determine if a is approximately equal to b, i.e., |a-b| <= 10*ε*max(1, |a|, |b|)
   use seecs__kinds, only: dp
   use seecs__constants, only: macheps
   real(dp), intent(in) :: a, b
   logical :: res
   res = .false.
   if(abs(b - a) .gt. 10*macheps*max(1.0_dp, abs(a), abs(b))) return
   res = .true.
 end function is_approx_real

! ================================================================================================================================ !
end module seecs__utils
! ================================================================================================================================ !
