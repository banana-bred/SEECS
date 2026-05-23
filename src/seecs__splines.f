! ================================================================================================================================ !
module seecs__splines
  !! General splines module

  use seecs__kinds,       only: dp
  use bspline_sub_module, only: db1ink, db1val, dintrv, dbspvd

  implicit none (type, external)

  private

  public :: spline_data_type
  public :: db1ink, db1val, dintrv, dbspvd
  public :: build_knots
  public :: gauss_legendre
  public :: wf_eval

  type spline_data_type
    integer :: nx = 0
      !! The number of interpolation points in \(x\)
    integer :: kx
      !! The B-spline order (k = p+1)
    integer :: inbvx
      !! Initialization parameter that should not be changed by the user
    real(dp), allocatable :: tx(:)
      !! Sequence of knots defining the piecewise polynomial representing the potential
    real(dp), allocatable :: bcoef(:)
      !! The B-spline coefficients computed by [[db1ink]]
  end type spline_data_type

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine build_knots(xmin, xmax, kx, nelem, tx, nbasis)
    !! Build clamped B-spline knots on [xmin,xmax]
    use seecs__system, only: die
    implicit none (type, external)
    real(dp), intent(in) :: xmin, xmax
      !! Interval endpoints
    integer, intent(in) :: kx
      !! The spline order
    integer, intent(in) :: nelem
      !! The number of knot spans in [xmin,xmax]
    real(dp), intent(out), allocatable :: tx(:)
      !! The knot vector
    integer, intent(out) :: nbasis
      !! The number of basis elements, nelem + (kx - 1)
    integer :: i
    integer :: nknots
    real(dp) :: hx
    if(nelem .lt. 1)    call die("Must have at least one knot span")
    if(xmax  .le. xmin) call die("Must have xmin < xmax")
    if(kx .lt. 1) call die("The order kx must be > 0")
    nknots = nelem + 2*kx - 1
    nbasis = nelem + kx - 1
    allocate(tx(nknots))
    ! -- xmin
    tx(1:kx)= xmin
    ! -- xmax
    tx(nknots-kx+1:nknots) = xmax
    if(nelem .eq. 1) return
    ! -- (xmin,xmax)
    hx = (xmax - xmin) / nelem
    do i=1,nelem-1 ; tx(kx + i) = xmin + hx*i ; enddo
  end subroutine build_knots

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine gauss_legendre(n, x, w)
    !! n-point Golub-Welsch Gauss-Legendre quadrature on [-1,1]
    use seecs__system, only: die
    use seecs__linalg, only: dstev
    implicit none (type, external)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n), w(n)
      ! -- dstev variables
      integer               :: info
      real(dp), allocatable :: e(:)
      real(dp), allocatable :: J(:,:)
      real(dp), allocatable :: work(:)
    integer :: i
    if(n.lt.1) call die("Need n >= 1")
    allocate(e(n-1), J(n,n), work(2*n-2))
    ! -- build diagonal (x) and off-diagonal (e)
    x = 0 ; do i=1,n-1 ; e(i) = i / sqrt(real(4*i*i - 1, kind=dp)) ; enddo
    call dstev('V', n, x, e, J, n, work, info)
    ! -- weights
    do i=1,n ; w(i) = 2*J(1,i)**2 ; enddo
    deallocate(e, J, work)
  end subroutine gauss_legendre

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine wf_eval(kx, tx, idx, wfcoeffs, Rgrid, wfs, R0, theta)
    !! Evaluate the wavefunctions ψ(R) on a real/complex grid using the B-spline coefficients

    use seecs__arrays, only: size_check
    use seecs__system, only: die

    implicit none (type, external)

    integer, intent(in) :: kx
      !! B-spline order
    real(dp), intent(in) :: tx(:)
      !! knot vector
    integer, intent(in) :: idx(:)
      !! List of nactive active B-spline indices. Enforcing BCs like ψ(0) = 0 or ψ'(Rmax)=0
      !! reduces the full basis to a smaller effective basis  that is used to solve
      !! the reduced problem
    class(*), intent(in) :: wfcoeffs(:,:)
      !! The nactive x nm bound-state wavefunction coefficients in the B-spline basis.
      !! There are nm ≤ nactive bound states to be returned by the user
    real(dp), intent(in) :: Rgrid(:)
      !! The nR real-valued grid on which to evaluate the wavefunctions
    class(*), intent(inout) :: wfs(:,:)
      !! The nR x nm wavefunctions to be returned
    real(dp), intent(in), optional :: R0
      !! The scaling radius
    real(dp), intent(in), optional :: theta
      !! The scaling angle

    logical, parameter :: extrap = .false.
    logical :: use_ecs

    integer  :: nt, nbasis, nactive, nm, nR
    integer  :: i, ilo, iR, m, icoeff, gi
    integer  :: lwork, iflag, mflag, ileft, ldv
    real(dp) :: R
    real(dp) :: tleft, tright, tshift
    integer,  allocatable :: wfmap(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: vnikx(:,:)

    nt      = size(tx, 1)
    nbasis  = nt - kx
    nactive = size(idx, 1)
    nm      = size(wfcoeffs, 2)
    nR      = size(Rgrid, 1)
    lwork   = kx + 1
    ldv     = kx

    ! -- size checks
    call size_check(wfcoeffs, [nactive, nm], "WFCOEFFS")
    call size_check(wfs,      [nR,      nm], "WFS")

    ! -- basis index → rows of wfcoeffs
    allocate(wfmap(nbasis)) ; wfmap = 0
    do concurrent (i=1:nactive)
      wfmap(idx(i)) = i
    enddo

    ! -- knot bounds, with tiny shift to ensure we don't get ileft out of bounds for evaluating on the endpoints
    tleft  = tx(kx)
    tright = tx(nbasis + 1)
    tshift = 10 * epsilon(1._dp) * max(1.0_dp, abs(tleft), abs(tright))

    allocate(vnikx(ldv,1))
    allocate(work((lwork*(lwork+1))/2))
    ilo = 1

    use_ecs = .false.
    if(present(R0) .AND. present(theta)) then
      use_ecs = .true.
    elseif(present(R0) .OR. present(theta)) then
      call die("R0 XOR THETA was passed to the wavefunction evaluator. Either both or neither must be present")
    endif

    ! -- wfs = 0
    call initialize_wfs(wfs)

    ! -- loop over gemetries R
    do iR = 1, nR

      ! -- the grid point R, shfited by a tiny amount to avoid evaluation on xmin or xmax resulting in out of bounds ILEFT
      R = Rgrid(iR)

      ! -- R > R0, intepret the value as ψ(z(R)), basis evaluation remains at x0 = R
      R = min(max(R, tleft + tshift), tright - tshift)

      call dintrv(tx, nt, R, ilo, ileft, mflag, extrap)

      ! -- check bounds on ileft from dintrv
      if(ileft .lt. kx)      call die("ILEFT < kx, possibly trying to evaluate the wavefunctions at x < xmin")
      if(ileft .gt. nt - kx) call die("ILEFT > nt - kx, possibly trying to evaluate the wavefunctions at x > xmax")

      call dbspvd(tx, kx, 1, R, ileft, ldv, vnikx, work, iflag)
      if(iflag .ne. 0) call die("DBSPVD failed in wf evaluation")

        do i=1,ldv

          ! -- active basis index at R
          gi = ileft - kx + i

          ! -- the corresponding row in wfcoeffs
          icoeff = wfmap(gi)

          if(icoeff .eq. 0) cycle

          ! -- wfs(iR, m=1:nm) += vnikx(i, 1) * wfcoeffs(icoeff, m=1:nm)
          call wfcoeff2wf(wfs(iR, 1:nm), vnikx(i, 1), wfcoeffs(icoeff, 1:nm))

        enddo

    enddo

  contains

    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure subroutine initialize_wfs(wfs)
      !! Initialize the given wavefunction to zero
      implicit none (type, external)
      class(*), intent(inout) :: wfs(:,:)
      select type(wfs)
      type is(real(dp))
        wfs = 0
      type is (complex(dp))
        wfs = 0
      class default
        call die("Trying to initalize WFS, but it's neither real nor complex")
      end select
    end subroutine initialize_wfs

    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure elemental subroutine wfcoeff2wf(wf, vnikx, wfcoeff)
      !! wf += vnikx * wfcoeff
      implicit none (type, external)
      class(*), intent(inout) :: wf
      real(dp), intent(in)  :: vnikx
      class(*), intent(in)  :: wfcoeff
      select type(wfcoeff)
      type is (real(dp))
        select type(wf)
        type is (real(dp))
          wf = wf + vnikx * wfcoeff
        class default
          call die("Trying to evaluate a nonreal WFS with a real WFCOEFFS")
        end select
      type is (complex(dp))
        select type(wf)
        type is (complex(dp))
          wf = wf + vnikx * wfcoeff
        class default
          call die("Trying to evaluate a non-complex WFS with a complex WFCOEFFS")
        end select
      class default
        call die("WFCOEFFS is neither real nor complex")
      end select
    end subroutine wfcoeff2wf

  end subroutine wf_eval


! ================================================================================================================================ !
end module seecs__splines
! ================================================================================================================================ !
