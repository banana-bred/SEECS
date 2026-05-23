! ================================================================================================================================ !
module seecs__hamilton

  use seecs__kinds,  only: dp
  use seecs__system, only: die

  implicit none (type, external)

  private

  save

  public :: build_hamiltonian
  public :: solve_schrodinger
  public :: Vinterp
  public :: Vspline

  interface build_hamiltonian
    module procedure :: build_hamiltonian_real
    module procedure :: build_hamiltonian_ecs
  end interface build_hamiltonian

  interface Vspline
    module procedure :: Vspline_real
    module procedure :: Vspline_cmplx
  end interface Vspline

  interface solve_schrodinger
    module procedure :: solve_schrodinger_real
    module procedure :: solve_schrodinger_cmplx
  end interface solve_schrodinger

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine Vinterp(xarr, Varr, k, Vspline_data)
    !! Interpolate the potential V(x) with a degree-k interpolant

    use seecs__arrays,  only: realloc
    use seecs__splines, only: spline_data_type, db1ink

    implicit none (type, external)

    real(dp), intent(in) :: xarr(:), Varr(:)
    integer,  intent(in) :: k
      !! The order k = p + 1, where p is the degree
    type(spline_data_type), intent(inout) :: Vspline_data
    integer, parameter :: iknot = 0 ! <-- let db1ink choose knots
    integer :: iflag

    associate(                        &
        nx    => Vspline_data % nx    &
      , kx    => Vspline_data % kx    &
      , inbvx => Vspline_data % inbvx &
      )

      nx = size(xarr, 1)
      if(nx .ne. size(Varr, 1)) call die("x-grid and V array need to be of the same size")
      kx = k

      ! -- potential knots and knot coefficients to be stored as module variables
      call realloc(Vspline_data % tx, nx+kx)
      call realloc(Vspline_data % bcoef, nx)

      ! -- interpolate
      call db1ink(xarr, Vspline_data % nx, Varr, Vspline_data % kx, iknot, Vspline_data%tx, Vspline_data%bcoef, iflag)

      if(iflag .ne. 0) call die("db1ink failed in potential interpolation")

      ! -- must be 1 on first call of db1val and not changed after
      inbvx = 1

    end associate

  end subroutine Vinterp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function Vspline_real(R, Vspline_data, w0) result(VR)
    !! Return the potential at point R after it has been interpolated
    use seecs__splines, only: spline_data_type, db1val
    implicit none (type, external)
    real(dp),          intent(in)    :: R
    type(spline_data_type), intent(inout) :: Vspline_data
    real(dp),          intent(inout) :: w0(3*Vspline_data%kx)
    real(dp) :: VR
    logical, parameter :: extrap = .false. ! -- do not extrapolate
    integer, parameter :: idx = 0 ! -- evaluate the interpolant (0th order derivative)
    integer :: iflag
    associate(                        &
        nx    => Vspline_data % nx    &
      , kx    => Vspline_data % kx    &
      , tx    => Vspline_data % tx    &
      , inbvx => Vspline_data % inbvx &
      , bcoef => Vspline_data % bcoef &
      )
      if(nx .eq. 0) call die("The potential must fist be interpolated before it can be evaluated")
      call db1val(R, idx, tx, nx, kx, bcoef, VR, iflag, inbvx, w0, extrap)
    end associate
    if(iflag .ne. 0) call die("Error detected in db1val")
  end function Vspline_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function Vspline_cmplx(z, Vspline_data, w0) result(Vz)
    !! Return the potential V(z) for complex z after the potential has been
    !! interpolated on a real-valued grid
    use seecs__system,  only: die
    use seecs__splines, only: spline_data_type, dintrv
    implicit none (type, external)
    complex(dp),          intent(in)    :: z
      !! The complex-valued evaluation point of the potential,
      !! \(z  = R_0 + (x_0 - R_0)e^{i\theta}\),
      !! where \(x_0\) is the physical quadrature point inside the current knot element
      !! and \(R_0\) is the scaling radius for the ECS grid
    type(spline_data_type), intent(inout) :: Vspline_data
    real(dp),          intent(inout) :: w0(3*Vspline_data%kx)
    complex(dp) :: Vz
    integer :: kx, nt
    integer :: r, j, lidx, ridx, ileft, ilo, mflag
    real(dp) :: tl, tr
    complex(dp) :: a
    complex(dp), allocatable :: d(:)
    kx = Vspline_data % kx

    ! -- allocation checks
    if(.not. allocated(Vspline_data % tx)) call die("Vspline knot array TX has not been initialized,&
      & but we're trying to evaluate the potential at a complex point z. Make sure that the potential is&
      & interpolated before it is evaluated.")
    nt = size(Vspline_data % tx, 1)
    if(.not. allocated(Vspline_data % bcoef)) call die("Vspline coefficient array BCOEF has not been initialized,&
      & but somehow its know array TX has been. Something is very strange")

    nt = size(Vspline_data % tx, 1)

    ! -- make sure we're evaluating the potential at Re(z) and not x0
    call dintrv(Vspline_data%tx, nt, z%re, ilo, ileft, mflag, extrap = .false.)

    ! -- ILEFT bounds check to make sure we're not evaluating where we don't have a potential
    if(ileft .lt. kx)      call die("ILEFT < kx, possibly trying to evaluate the wavefunctions at Re(z) < xmin")
    if(ileft .gt. nt - kx) call die("ILEFT > nt - kx, possibly trying to evaluate the wavefunctions at Re(z) > xmax")

    d = cmplx(Vspline_data % bcoef(ileft-kx+1:ileft), kind=dp)

    ! -- de Boor recursion
    do r = 1, kx-1
      do j = kx, r+1, -1
        lidx = ileft - (kx - 1) + j ! left index:  i-k+1+j
        ridx = ileft + j   + 1  - r  ! right index: i+j+1-r
        tl = Vspline_data % tx(lidx)
        tr = Vspline_data % tx(ridx)
        if(tr .eq. tl) then
          ! -- multiple knots, carry over the left value
          d(j) = d(j-1)
        else
          a = (z-tl) / (tr-tl)
          d(j) = (1._dp-a)*d(j-1) + a*d(j)
        endif
      enddo
    enddo

    Vz = d(kx)

  end function Vspline_cmplx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine  build_hamiltonian_ecs(kx, tx, mu, j, nq, xq, wq, R0, theta, Vspline_data, w0, H, S)
    !! Build the complex-valued matrices H and S, i.e. the general case of ECS
    use seecs__system,    only: die
    use seecs__arrays,    only: size_check
    use seecs__splines,   only: spline_data_type, dintrv, dbspvd
    use seecs__constants, only: im
    implicit none (type, external)
    integer, intent(in) :: kx
      !! B-spline order kx = p+1
    real(dp), intent(in) :: tx(:)
      !! the knot vector
    real(dp), intent(in) :: mu
      !! the reduced mass in atomic units (not atomic mass units)
    integer, intent(in) :: j
      !! the rotational quantum number
    integer, intent(in) :: nq
      !! number of Gauss-Legendre points per element
    real(dp), intent(in) :: xq(nq)
      !! Gauss-Legendre nodes on [-1,1]
    real(dp), intent(in) :: wq(nq)
      !! Gauss-Legendre weights on [-1,1]
    real(dp), intent(in) :: R0
      !! The scaling radius \(R_0\)
    real(dp), intent(in) :: theta
      !! The scaling angle \(\theta\)
    type(spline_data_type), intent(inout) :: Vspline_data
      !! the inteprolated internuclear potential
    real(dp), intent(inout) :: w0(3*kx)
      !! workspace array for spline evaluation
    complex(dp), intent(out) :: H(:,:)
      !! The Hamiltonian
    complex(dp), intent(out) :: S(:,:)
      !! The overlap integrals

    integer :: e, i, k, q
    integer :: gi, gk, g0
    integer :: iflag, ileft, ilo, mflag
    integer :: p, nt, nbasis, nelem
    integer :: ldvnik, lwork
    real(dp) :: x0
    real(dp) :: Bi, Bk, dBi, dBk
    real(dp) :: tr, tl
    real(dp) :: jacobi
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: vnikx(:,:)
    complex(dp) :: z, Vz
    complex(dp) :: expo, expoinv

    p      = kx-1
    nt     = size(tx, 1)
    nbasis = nt - kx
    nelem  = nt - 2*kx + 1
    ldvnik = kx
    lwork  = kx + 1

    ! -- must be initialized to 1 the first time [[dintrv]] is called, then is controlled by
    !    subsequent [[dintrv]] calls
    ilo = 1

    call size_check(H, [nbasis, nbasis], "H")
    call size_check(S, [nbasis, nbasis], "S")
    H = 0; S = 0

    allocate(vnikx(ldvnik, 2))
    allocate(work((lwork*(lwork+1))/2))

    elements: do e = 1, nelem

      ! -- pick left and right knots of the element e
      tl = tx(p+e)
      tr = tx(p+e+1)

      ! -- the Jacobian J of the map [-1,1] → [tl, tr]
      jacobi = (tr - tl)/2

      nodes: do q = 1, nq

        ! -- xq(q) → x0
        x0 = (tr+tl)/2 + jacobi*xq(q)

        if(x0 .le. R0) then
          ! -- below the scaling radius R0
          expo    = 1
          z       = x0
          expoinv = 1
        else
          ! -- above the scaling radius R0
          expo    = exp(im*theta)
          expoinv = exp(-im*theta)
          z       = R0 + (x0-R0)*expo
        endif

        ! -- find the knot span for the real value x0
        call dintrv(tx, nt, x0, ilo, ileft, mflag, extrap = .false.)

        if(ileft .lt. kx) call die("ileft < kx ⇒  S and H will not be properly indexed")

        ! -- evaluate kx nonzero B-splines and first derivatives at x0
        call dbspvd(tx, kx, 2, x0, ileft, ldvnik, vnikx, work, iflag)

        if(iflag .ne. 0) call die("Failed in spline evaluation while building the complex Hamiltonian")

        ! -- V(r) + j(j+1)/r² at z = r + iy
        Vz = Vspline(z, Vspline_data, w0) + j*(j+1)/(2*mu*z*z)

        g0 = ileft - kx
        do k=1, ldvnik
          gk = g0 + k
          Bk = vnikx(k, 1)
          dBk = vnikx(k, 2)

          do i=1, ldvnik
            gi = g0 + i
            Bi = vnikx(i, 1)
            dBi = vnikx(i, 2)

            S(gi, gk) = S(gi, gk) + wq(q)*jacobi * expo*Bi*Bk
            H(gi, gk) = H(gi, gk) + wq(q)*jacobi * ( expoinv*(dBi*dBk / (2*mu)) + expo*(Bi*Vz*Bk) )
          enddo
        enddo

      enddo nodes
    enddo elements

    deallocate(vnikx, work)

  end subroutine build_hamiltonian_ecs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine build_hamiltonian_real(kx, tx, mu, j, nq, xq, wq, Vspline_data, w0, H, S)
    !! Build the real and symmetric matrices H and S, i.e. no ECS for real-valued bound/box states only
    use seecs__system,  only: die
    use seecs__arrays,  only: size_check
    use seecs__splines, only: spline_data_type, dintrv, dbspvd

    implicit none (type, external)

    integer, intent(in) :: kx
      !! B-spline order kx = p+1
    real(dp), intent(in) :: tx(:)
      !! the knot vector
    real(dp), intent(in) :: mu
      !! the reduced mass in atomic units (not atomic mass units)
    integer, intent(in) :: j
      !! the rotational quantum number
    integer, intent(in) :: nq
      !! number of Gauss-Legendre points per element
    real(dp), intent(in) :: xq(nq)
      !! Gauss-Legendre nodes on [-1,1]
    real(dp), intent(in) :: wq(nq)
      !! Gauss-Legendre weights on [-1,1]
    type(spline_data_type), intent(inout) :: Vspline_data
      !! the inteprolated internuclear potential
    real(dp), intent(inout) :: w0(3*kx)
      !! workspace array for spline evaluation
    real(dp), intent(out) :: H(:,:)
      !! The Hamiltonian
    real(dp), intent(out) :: S(:,:)
      !! The overlap integrals

    logical, parameter :: extrap = .false.
    integer :: nbasis, nelem
    integer :: p, nt
    integer :: ldvnik, lwork
    integer :: e, q, i, k
    integer :: ilo
    integer :: iflag, ileft, mflag
    integer :: gi0, gi, gk
    real(dp) :: tr, tl, jacobi, R, V
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: vnikx(:,:)

    p = kx - 1
    ldvnik = kx
    lwork  = kx + 1

    nt = size(tx, 1)
    nbasis = nt - kx
    nelem = nt - 2*kx + 1

    ! -- must be initialized to 1 the first time [[dintrv]] is called, then is controlled by
    !    subsequent [[dintrv]] calls
    ilo = 1

    call size_check(H, [nbasis, nbasis], "H")
    call size_check(S, [nbasis, nbasis], "S")
    ! if( any(shape(H) .ne. [nbasis, nbasis]) ) call die("Hamiltonian array H has the wrong size")
    ! if( any(shape(S) .ne. [nbasis, nbasis]) ) call die("Overlap array S has the wrong size")

    H = 0 ; S = 0
    allocate(vnikx(ldvnik, 2))
    allocate(work((lwork*(lwork+1))/2))

    ! -- loop over the knot intervals (elements). H and S are evaluated by taking the integrals
    !    element by element, then summing over elements
    elements: do e = 1, nelem

      ! -- pick left and right knots of the element e
      tl = tx(p+e)
      tr = tx(p+e+1)

      ! -- the Jacobian J of the map [-1,1] → [tl, tr]
      jacobi = (tr - tl)/2

      nodes: do q = 1, nq

        ! -- xq(q) → R
        R = (tr+tl)/2 + jacobi*xq(q)

        ! -- V(R) + j(j+1)/2μR²
        V = Vspline(R, Vspline_data, w0) + j*(j+1)/(2*mu*R*R)

        ! -- find the knot  span containing R. Returns ileft such that
        !        tx(ileft) ≤ R < tx(lieft+1)
        !    good for knowing which kx splines are nonzero at x
        call dintrv(tx, nt, R, ilo, ileft, mflag, extrap)

        if(ileft .lt. kx) call die("ileft < kx ⇒  S and H will not be properly indexed")

        ! -- evaluate the kx nonzero splines vnikx(:,1) and their first derivatives vnikx(:,2) at R
        call dbspvd(tx, kx, 2, R, ileft, ldvnik, vnikx, work, iflag)

        ! -- build S and H
        ! S += w Bi(R) Bk(R) J
        ! H += T + V
        !   T = w (Bi'(R) Bk'(R) / 2μ) J
        !   V = w Bi(R) V Bk(R) J
        gi0 = ileft - kx
        do k = 1, ldvnik
          gk  = gi0 + k
          do i = 1, ldvnik
            gi  = gi0 + i
            S(gi, gk) = S(gi, gk) + wq(q) * vnikx(i,1) * vnikx(k,1) * jacobi
            H(gi, gk) = H(gi, gk) + wq(q) * ( vnikx(i,2)*vnikx(k,2)/(2*mu) + vnikx(i,1)*V*vnikx(k,1) ) * jacobi
          enddo
        enddo

      enddo nodes

    enddo elements

    deallocate(vnikx, work)

  end subroutine build_hamiltonian_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine solve_schrodinger_real(H, S, idx, energies, wfcoeffs)
    !! Solve the real-valued time-independent rovibrational schrödinger equation
    !! for the bound-state wavefunctions and energies given the Hamiltonian matrix H and the
    !! B-spline overlap matrix S.

    use seecs__system,     only: die
    use seecs__arrays,     only: size_check
    use seecs__linalg,     only: dsygv
    use seecs__characters, only: i2c => int2char

    implicit none (type, external)

    real(dp), intent(in)  :: H(:,:)
      !! The full Hamiltonian matrix
    real(dp), intent(in)  :: S(:,:)
      !! The full B-spline overlap matrix
    integer,  intent(in)  :: idx(:)
      !! List of n active B-spline indices. Enforcing BCs like ψ(0) = 0 or ψ'(Rmax)=0
      !! reduces the full basis to a smaller effective basis  that is used to solve
      !! the reduced problem
    real(dp), intent(out) :: energies(:)
      !! The n bound-state energies
    real(dp), intent(out) :: wfcoeffs(:,:)
      !! The n x n bound-state wavefunction coefficients in the B-spline basis.
      !! To be evaluated later for ψ(R)

    integer :: n
    integer :: info, lwork
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: Heff(:,:), Seff(:,:)

    n = size(idx, 1)
    lwork = 3*n-1

    ! -- size checks
    call size_check(energies, n, "ENERGIES")
    call size_check(wfcoeffs, [n,n], "WFCOEFFS")

    allocate(work(lwork))

    ! -- build the effective matrices
    Heff = H(idx(:), idx(:))
    Seff = S(idx(:), idx(:))

    ! -- solve the generalized eigenvalue problem Hc = ESc
    call dsygv(1, "V", "U", n, Heff, n, Seff, n, energies, work, lwork, info)

    if(info .ne. 0) call die("DSYGV returned a nonzero value for INFO")

    wfcoeffs = Heff

    deallocate(work, Heff, Seff)

  end subroutine solve_schrodinger_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine solve_schrodinger_cmplx(H, S, idx, energies, reigvecs)
    !! Solve the complex-valued time-independent rovibrational schrödinger equation
    !! for the bound-state wavefunctions and energies given the Hamiltonian matrix H and the
    !! B-spline overlap matrix S.

    use seecs__system,     only: die, stderr, stdout
    use seecs__arrays,     only: size_check, is_symmetric, realloc, sort_index
    use seecs__linalg,     only: zggev
    use seecs__characters, only: i2c => int2char

    implicit none (type, external)

    complex(dp), intent(in)  :: H(:,:)
      !! The full Hamiltonian matrix
    complex(dp), intent(in)  :: S(:,:)
      !! The full B-spline overlap matrix
    integer,  intent(in)  :: idx(:)
      !! List of n active B-spline indices. Enforcing BCs like ψ(0) = 0 or ψ'(Rmax)=0
      !! reduces the full basis to a smaller effective basis  that is used to solve
      !! the reduced problem
    complex(dp), intent(out) :: energies(:)
      !! The n bound-state energies
    complex(dp), intent(out) :: reigvecs(:,:)
      !! The n x n right eigenvectors; the coefficients of the wavefunctions in the B-spline basis.
      !! The Hamiltonian is complex-symmetric, so these double as the left eigenvectors too

    logical :: S_is_symmetric
    integer :: n
    integer :: info, lwork
    integer, allocatable :: perm(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: alpha(:), beta(:)
    complex(dp), allocatable :: work(:)
    complex(dp), allocatable :: Heff(:,:), Seff(:,:)
    complex(dp), allocatable :: A(:,:), B(:,:)
    complex(dp), allocatable :: leigvecs(:,:)
    character(1), parameter :: jobvr = "V"
    character(1), parameter :: jobvl = "N"

    n = size(idx, 1)
    lwork = 3*n-1

    ! -- size checks
    call size_check(energies, n, "ENERGIES")
    call size_check(reigvecs, [n,n], "REIGVECS")
    allocate(leigvecs(n,n))

    lwork = -1
    allocate(work(abs(lwork)))
    allocate(rwork(8*n))
    allocate(alpha(n), beta(n))

    ! -- build the effective matrices
    Heff = H(idx(:), idx(:)) ; A = Heff
    Seff = S(idx(:), idx(:)) ; B = Seff

    if(is_symmetric(Heff) .eqv. .false.) then
      write(stderr, '(A, "F9.6")') "maxval(abs(Heff - transpose(Heff))): ", maxval(abs(Heff - transpose(Heff)))
      call die("The Hamiltonian is not symmetric")
    endif
    if(is_symmetric(Seff) .eqv. .false.) then
      write(stderr, '(A, "F9.6")') "maxval(abs(Seff - transpose(Seff))): ", maxval(abs(Seff - transpose(Seff)))
      call die("The overlap matrix is not symmetric")
    endif

    ! -- determine the optimal LWORK
    print*, n
    call zggev(jobvl, jobvr, n, A, n, B, n, alpha, beta, leigvecs, n, reigvecs, n, work, lwork, rwork, info)
    lwork = nint(work(1) % re)

    call realloc(work, lwork)

    ! -- solve the generalized eigenvalue problem Hc = ESc
    call zggev(jobvl, jobvr, n, A, n, B, n, alpha, beta, leigvecs, n, reigvecs, n, work, lwork, rwork, info)

    if(info .ne. 0) call die("ZGGEV returned a nonzero value for INFO")

    energies = alpha/beta

    call biorthonormalize_R(Seff, reigvecs)
    ! if(is_symmetric(Seff)) then
    !   write(stdout, '(A)') "Detected complex-symmetric B-spline overlap matrix.&
    !     & Normalizing right eigenvectors because they are the left eigenvectors."
    !   call biorthonormalize_R(Seff, reigvecs)
    ! else
    !   write(stdout, '(A)') "Detected complex-nonsymmetric B-spline overlap matrix.&
    !     & Normalizing left eigenvectors."
    !   call biorthonormalize_LR(Seff, leigvecs, reigvecs)
    ! endif

    ! -- sort the ECS energies/states because they are unordered from ZGGEV
    allocate(perm(n)) ; perm = 0
    call sort_index(energies, perm)
    energies = energies(perm)
    reigvecs = reigvecs(:, perm)

    deallocate(work, Heff, Seff, leigvecs)

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !

    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure subroutine biorthonormalize_R(M, rvecs)
      !! Given the right eigenvectors of a complex-symmetric eigenvalue problem,
      !! scale the right eigenvectors such that
      !!   \(R^T_i M R_j = \delta_{ij}\)
      use seecs__system,     only: die
      use seecs__characters, only: i2c => int2char
      implicit none (type, external)
      complex(dp), intent(in) :: M(:,:)
      complex(dp), intent(inout) :: rvecs(:,:)
      integer :: n, k
      complex(dp) :: prod
      complex(dp), allocatable :: MR(:,:)
      n = size(M, 1)
      MR = matmul(M,  rvecs)
      do concurrent (k=1:n)
        prod = sum(rvecs(:,k) * MR(:,k))
        if(abs(prod) .eq. 0) call die("Identically zero biorthogonal norm for mode " // i2c(k))
        rvecs(:,k) = rvecs(:,k) / sqrt(prod)
      enddo
    end subroutine biorthonormalize_R
    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure subroutine biorthonormalize_LR(M, lvecs, rvecs)
      !! Given the left and right eigenvectors of a generalized eigenvalue problem,
      !! scale the left eigenvectors such that
      !!   \(L^\dagger_i M R_j = \delta_{ij}\)
      use seecs__system,     only: die
      use seecs__characters, only: i2c => int2char
      implicit none (type, external)
      complex(dp), intent(in) :: M(:,:)
      complex(dp), intent(inout) :: lvecs(:,:)
      complex(dp), intent(in) :: rvecs(:,:)
      integer :: n, k
      complex(dp) :: prod
      complex(dp), allocatable :: MR(:,:)
      n = size(M, 1)
      MR = matmul(M,  rvecs)
      do concurrent (k=1:n)
        ! -- dot_product of complex arrays is acutally the inner product. confusing name !
        prod = dot_product(lvecs(:,k), MR(:,k))
        if(abs(prod) .eq. 0) call die("Identically zero biorthogonal norm for mode " // i2c(k))
        lvecs(:,k) = lvecs(:,k) / sqrt(prod)
      enddo
    end subroutine biorthonormalize_LR

  end subroutine solve_schrodinger_cmplx

! ================================================================================================================================ !
end module seecs__hamilton
! ================================================================================================================================ !
