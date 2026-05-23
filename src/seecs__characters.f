! ================================================================================================================================ !
module seecs__characters

  implicit none (type, external)

  private

  public :: int2char
  public :: ndigits
  public :: to_lower
  public :: lower

  character(13), parameter, public :: numeric = "0123456789.+-"

  interface int2char
    module procedure :: scalar_int2char
    module procedure :: vector_int2char
  end interface int2char

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! -------------------------------------------------------------------------------------------------------------------------------!
  pure elemental function ndigits(n) result(num)
    !! Returns number of characters an integer will occupy
    use seecs__kinds, only: dp
    implicit none (type, external)
    integer, intent(in) :: n
    integer :: num
    real(dp), parameter :: one = 1._dp
    num = 1
    if(n .eq. 0) return
    num = floor(log10(abs(n) * one)) + 1
    ! -- account for minus sign
    if(n.lt.1) num = num + 1
  end function ndigits

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function scalar_int2char(i) result(res)
    !! Writes the value i to a character as I0
    implicit none (type, external)
    integer, intent(in) :: i
    character(:), allocatable :: res
    allocate(character(ndigits(i)) :: res)
    write(res, '(I0)') i
  end function scalar_int2char

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function vector_int2char(i) result(res)
    !! Writes the value i to a character as I0
    implicit none (type, external)
    integer, intent(in) :: i(:)
    character(:), allocatable :: res
    integer :: j, n
    if(size(i,1) .eq. 0) then
      res = "()"
      return
    endif
    ! -- n-1 commas, and all the digits
    n = (size(i,1) - 1) + sum(ndigits(i))
    allocate(character(n) :: res)
    write(res, '(*(I0,:,","))') i
    res = "(" // res // ")"
  end function vector_int2char

  ! -------------------------------------------------------------------------------------------------------------------------------!
  pure elemental subroutine to_lower(chr)
    !! converts a character to lower case
    use seecs__constants, only: uppercase_a, uppercase_z
    implicit none (type, external)
    character(*), intent(inout) :: chr
    integer :: i
    integer :: n
    integer :: ic
    n = len(chr)
    do i = 1, n
      ic = ichar(chr(i:i))
      if(ic .lt. uppercase_a) cycle
      if(ic .gt. uppercase_z) cycle
      chr(i:i) = char(ic + 32)
    enddo
  end subroutine to_lower

  ! -------------------------------------------------------------------------------------------------------------------------------!
  pure function lower(chr) result(res)
    !! returns a lower case character
    implicit none (type, external)
    character(*), intent(in) :: chr
    character(:), allocatable :: res
    integer, parameter :: shift = ichar('a') - ichar("A")
    integer, parameter :: uppercase_a = ichar('A')
    integer, parameter :: uppercase_z = ichar('Z')
    integer :: i, n, ic
    n = len(chr)
    res = chr
    do i = 1, n
      ic = ichar(res(i:i))
      ! -- cycle if the character isn't in [A,Z]
      if(ic .lt. uppercase_a) cycle
      if(ic .gt. uppercase_z) cycle
      res(i:i) = char(ic + shift)
    enddo
  end function lower

! ================================================================================================================================ !
end module seecs__characters
! ================================================================================================================================ !
