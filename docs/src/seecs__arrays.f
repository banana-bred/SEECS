! ================================================================================================================================ !
module seecs__arrays
  use seecs__kinds, only: dp

  implicit none (type, external)

  private

  public :: realloc
  public :: size_check
  public :: norm_frob
  public :: is_symmetric
  public :: is_unitary
  public :: eye
  public :: adjoint
  public :: sort_index

  interface realloc
    module procedure :: realloc_1d_real
    module procedure :: realloc_1d_cmplx
  end interface realloc

  interface size_check
    module procedure :: size_check_1d
    module procedure :: size_check_2d
  end interface size_check

  interface norm_frob
    module procedure :: norm_frob_i
    module procedure :: norm_frob_r
    module procedure :: norm_frob_c
  end interface norm_frob

  interface adjoint
    module procedure :: adjoint_i
    module procedure :: adjoint_r
    module procedure :: adjoint_c
  end interface adjoint

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_i(A) result(res)
    !! Returns the adjoint of an integer-valued matrix
    implicit none (type, external)
    integer, intent(in) :: A(:,:)
    integer :: res(size(A, 2), size(A, 1))
    res = transpose(A)
  end function adjoint_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_r(A) result(res)
    !! Returns the adjoint of a real-valued matrix
    implicit none (type, external)
    real(dp), intent(in) :: A(:,:)
    real(dp) :: res(size(A, 2), size(A, 1))
    res = transpose(A)
  end function adjoint_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_c(A) result(res)
    !! Returns the adjoint of a complex-valued matrix
    implicit none (type, external)
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: res(size(A, 2), size(A, 1))
    res = conjg(transpose(A))
  end function adjoint_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_real(arr, n)
    implicit none (type, external)
    real(dp), intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_real
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_cmplx(arr, n)
    implicit none (type, external)
    complex(dp), intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_cmplx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine size_check_1d(arr, larr, name)
    !! Check that the size of the array arr is of length larr
    use seecs__system,     only: die
    use seecs__characters, only: i2c => int2char
    implicit none (type, external)
    class(*), intent(in) :: arr(:)
    integer, intent(in) :: larr
    character(*), intent(in) :: name
    if(size(arr, 1) .ne. larr) &
      call die("Array " // name // "(:) " // i2c(shape(arr)) // " must have the shape " // i2c([larr]))
  end subroutine size_check_1d
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine size_check_2d(arr, larr, name)
    !! Check that the size of the array arr is of length larr
    use seecs__system,     only: die
    use seecs__characters, only: i2c => int2char
    implicit none (type, external)
    class(*), intent(in) :: arr(:,:)
    integer, intent(in) :: larr(:)
    character(*), intent(in) :: name
    if(any(shape(arr) .ne. larr)) &
      call die("Array " // name // "(:,:) " // i2c(shape(arr)) // " must have the shape " // i2c(larr))
  end subroutine size_check_2d

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function unitary_defect(A) result(rF)
    !! Return the unitary defect with respect to the Frobenius norm
    !! \(rF = ||A^{\dagger}A-I||_F / sqrt{n}\)
    use seecs__system, only: die
    implicit none (type, external)
    class(*), intent(in) :: A(:,:)
    real(dp) :: rF
    integer :: n, m
    integer :: i, j
    real(dp) :: err2

    n = size(A, 1)
    m = size(A, 2)

    if(n .ne. m) call die("Cannot determine the unitarity of a nonsquare matrix")

    err2 = 0

    select type(A)

    type is (real(dp))

    realmat: block
      real(dp), allocatable :: G(:,:)
      real(dp) :: d
      allocate(G(n,m))
      ! -- Gram matrix
      G = matmul(adjoint(A), A)
      do j=1,n
        do i=1,n
          ! -- subtract identity if i==j
          d = G(i,j) - merge(1, 0, i .eq. j)
          ! -- accumulate norm squared
          err2 = err2 + d*d
        enddo
      enddo
    end block realmat

    type is (complex(dp))

    cmplxmat: block
      complex(dp), allocatable :: G(:,:)
      complex(dp) :: d
      allocate(G(n,m))
      ! -- Gram matrix
      G = matmul(adjoint(A), A)
      do j=1,n
        do i=1,n
          ! -- subtract identity if i==j
          d = G(i,j) - merge(1, 0, i .eq. j)
          ! -- accumulate norm squared
          err2 = err2 + d*d
        enddo
      enddo
    end block cmplxmat

    class default

      call die("Trying to determine the unitary of a matrix that is neither real nor complex")

    end select

    ! -- relative size-aware defect
    rF = sqrt(err2) / sqrt(real(n, kind = dp))

  end function unitary_defect

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_unitary(A, rtol) result(res)
    !! Determine whether a matrix is unitary w.r.t the Frobenius norm
    use seecs__constants, only: macheps
    class(*), intent(in) :: A(:,:)
      !! The matrix
    real(dp), intent(in), optional :: rtol
      !! The optional relative tolerance, default 1e-10
    logical :: res
    integer :: n, m
    real(dp) :: rF, tol
    tol = 1e-10_dp ; if(present(rtol)) tol = rtol
    res = .false.
    n = size(A, 1)
    m = size(A, 2)
    if(n .ne. m) return
    rF = unitary_defect(A)
    tol = max(tol, 100*macheps*sqrt(real(n, kind=dp)))
    res = rF .lt. tol
  end function is_unitary

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_symmetric(A, rtol) result(res)
    !! Determine whether a matrix is symmetric
    use seecs__constants, only: macheps
    class(*), intent(in) :: A(:,:)
      !! The matrix
    real(dp), intent(in), optional :: rtol
      !! The optional relative tolerance, default 1e-10
    logical :: res
    integer :: n, m
    real(dp) :: rF, tol, diff, denom
    tol = 1e-10_dp ; if(present(rtol)) tol = rtol
    res = .false.
    n = size(A, 1)
    m = size(A, 2)
    if(n .ne. m) return
    select type(AA => A)
    type is (integer)
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    type is (real(dp))
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    type is (complex(dp))
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    end select
    res = diff .le. tol*denom + sqrt(macheps)
  end function is_symmetric

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function eye(n) result(res)
    !! Return an n x n identity matrix
    use seecs__system, only: die
    integer, intent(in) :: n
    integer :: res(n,n)
    integer :: i
    if(n .lt. 1) call die("Attempt to make an identity matrix with dims < 1")
    res = 0
    do concurrent(i=1:n)
      res(i,i) = 1
    enddo
  end function eye

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_i(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none (type, external)
    integer, intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(real(sum(A*A), kind=dp))
  end function norm_frob_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_r(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none (type, external)
    real(dp), intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(sum(A*A))
  end function norm_frob_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_c(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none (type, external)
    complex(dp), intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(sum(A*adjoint(A)))
  end function norm_frob_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine sort_index(vals, idx)
    !! Sort the array vals and return the permutation indices
    implicit none (type, external)
    complex(dp), intent(in)  :: vals(:)
    integer,     intent(out) :: idx(:)
    integer :: i, j, n
    integer :: key
    real(dp) :: targ
    n = size(vals, 1)
    call size_check(idx, n, "IDX")
    do concurrent(i=1:n) ; idx(i) = i ; enddo
    if(n .le. 1) return
    do i = 2, n
      key = idx(i)
      j = i - 1
      targ = vals(key) % re
      do
        if(j .lt. 1) exit
        if(vals(idx(j))%re .le. targ) exit
        idx(j+1) = idx(j)
        j = j - 1
      enddo
      idx(j+1) = key
    enddo
  end subroutine sort_index

! ================================================================================================================================ !
end module seecs__arrays
! ================================================================================================================================ !
