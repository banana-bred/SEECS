! ================================================================================================================================ !
module seecs__drivers

  use seecs__kinds, only: dp

  implicit none (type, external)

  private

  public :: rvsolve

  interface rvsolve
    module procedure :: rvsolve_noecs
    module procedure :: rvsolve_ecs
  end interface rvsolve

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine rvsolve_noecs( rvals, vvals, jrot, kx, nelem, redmass, ndropl, ndropr, nwf, wf &
                                 , energies, brv, nR_wf, R_wf, rmin, rmax)
    !! Real-valued driver (no ECS) to return NWF (ro)vibrational wavefunctions (and their energies)
    !! solved on the grid RVALS, on which the potential VVALS was evaluated for a given rotational
    !! level JROT. No ECS means that this only calculates bound/box states only.
    !!
    !! Wavefunctions are evaluated on a uniform grid R_WF of NR_WF points on [RMIN, RMAX].
    !! When RMIN/RMAX are not supplied, they are taken as the endpoints of the input grid.

    use seecs__arrays,    only: size_check
    use seecs__splines,   only: spline_data_type, wf_eval
    use seecs__hamilton,  only: build_hamiltonian, solve_schrodinger

    implicit none (type, external)

    real(dp), intent(in)  :: rvals(:)       !! Grid of R-values
    real(dp), intent(in)  :: vvals(:)       !! Grid of V(R)-values
    integer,  intent(in)  :: jrot           !! Rotational quantum number j/N
    integer,  intent(in)  :: kx             !! B-spline order (kx = p + 1)
    integer,  intent(in)  :: nelem          !! Number of knot intervals (should be >> NWF)
    real(dp), intent(in)  :: redmass        !! Reduced mass μ in atomic (electron-mass) units
    integer,  intent(in)  :: ndropl         !! Number of B-spline functions dropped on the left (small R) edge
    integer,  intent(in)  :: ndropr         !! Number of B-spline functions dropped on the right (large R) edge
    integer,  intent(in)  :: nwf            !! Number of wavefunctions/energies to return (lowest first)
    real(dp), intent(out) :: wf(:,:)        !! The NR_WF x NWF wavefunctions ψ(R), evaluated on R_WF
    real(dp), intent(out) :: energies(:)    !! The NWF energies
    real(dp), intent(out) :: brv(:)         !! Rotational constant B per state
    integer,  intent(in)  :: nR_wf          !! Number of wavefunction evaluation points
    real(dp), intent(out) :: R_wf(:)        !! The NR_WF uniform evaluation grid
    real(dp), intent(in), optional :: rmin  !! Lower bound of eval grid [default RVALS(1)]
    real(dp), intent(in), optional :: rmax  !! Upper bound of eval grid [default RVALS(N)]

    integer  :: nrvals, nbasis, nbasis_actual, nq, iv
    real(dp) :: w0(3*kx)
    integer,  allocatable :: idx(:)
    real(dp), allocatable :: tx(:), xq(:), wq(:)
    real(dp), allocatable :: H(:,:), S(:,:), energies_full(:), wfcoeffs(:,:)
    real(dp), allocatable :: W(:,:)
    complex(dp), allocatable :: Wfull(:,:)

    type(spline_data_type) :: vspline_data

    nrvals = size(rvals, 1)
    call size_check(vvals,    nrvals,        "VVALS")
    call size_check(energies, nwf,           "ENERGIES")
    call size_check(brv,      nwf,           "BRV")
    call size_check(wf,       [nR_wf, nwf],  "WF")

    call build_basis(rvals, vvals, kx, nelem, vspline_data, tx, nbasis, nq, xq, wq)
    call build_eval_grid(rvals, nR_wf, R_wf, rmin, rmax)
    call make_idx(nbasis, ndropl, ndropr, nwf, idx, nbasis_actual)

    ! -- real, symmetric H and S
    allocate(H(nbasis, nbasis), S(nbasis, nbasis))
    call build_hamiltonian(kx, tx, redmass, jrot, nq, xq, wq, vspline_data, w0, H, S)

    ! -- solve Hc = ESc in the reduced subspace
    allocate(energies_full(nbasis_actual), wfcoeffs(nbasis_actual, nbasis_actual))
    call solve_schrodinger(H, S, idx, energies_full, wfcoeffs)

    ! -- return the NWF lowest states, evaluated on R_WF
    energies(1:nwf) = energies_full(1:nwf)
    call wf_eval(kx, tx, idx, wfcoeffs(:, 1:nwf), R_wf, wf)

    ! -- rotational constant Bv= cv⁺ W cv
    allocate(Wfull(nbasis, nbasis))

    ! -- too lazy to split this up into real/cmplx
    call centrifugal_matrix(kx, tx, redmass, nq, xq, wq, Wfull)
    W = real(Wfull(idx, idx), kind=dp)

    do iv=1, nwf
      brv(iv) = sum(wfcoeffs(:, iv) * matmul(W, wfcoeffs(:, iv)))
    enddo

  end subroutine rvsolve_noecs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine rvsolve_ecs( rvals, vvals, jrot, kx, nelem, redmass, ndropl, ndropr, nwf, R0, theta &
                               , wf, energies, brv, nR_wf, R_wf, rmin, rmax)
    !! Complex-valued driver (allows ECS) to return NWF (ro)vibrational wavefunctions (and their energies)
    !! solved on the grid RVALS, on which the potential VVALS was evaluated for a given rotational
    !! level JROT. ECS means that this can calculate bound states in the continuum that are not just
    !! box states.
    !!
    !! Wavefunctions are evaluated on a uniform grid R_WF of NR_WF points on [RMIN, RMAX].
    !! When RMIN/RMAX are not supplied, they are taken as the endpoints of the input grid.

    use seecs__arrays,    only: size_check
    use seecs__splines,   only: spline_data_type, wf_eval
    use seecs__hamilton,  only: build_hamiltonian, solve_schrodinger

    implicit none (type, external)

    real(dp),    intent(in)  :: rvals(:)       !! Grid of R-values
    real(dp),    intent(in)  :: vvals(:)       !! Grid of V(R)-values
    integer,     intent(in)  :: jrot           !! Rotational quantum number j/N
    integer,     intent(in)  :: kx             !! B-spline order (kx = p + 1)
    integer,     intent(in)  :: nelem          !! Number of knot intervals (should be >> NWF)
    real(dp),    intent(in)  :: redmass        !! Reduced mass μ in atomic (electron-mass) units
    integer,     intent(in)  :: ndropl         !! Number of B-spline functions dropped on the left (small R) edge
    integer,     intent(in)  :: ndropr         !! Number of B-spline functions dropped on the right (large R) edge
    integer,     intent(in)  :: nwf            !! Number of wavefunctions/energies to return (lowest first)
    real(dp),    intent(in)  :: R0             !! ECS scaling radius (snapped to nearest interior knot)
    real(dp),    intent(in)  :: theta          !! ECS scaling angle, in radians
    complex(dp), intent(out) :: wf(:,:)        !! The NR_WF x NWF wavefunctions ψ(R), evaluated on R_WF
    complex(dp), intent(out) :: energies(:)    !! The NWF energies
    complex(dp), intent(out) :: brv(:)         !! Rotational constant B per state.
    integer,     intent(in)  :: nR_wf          !! Number of wavefunction evaluation points
    real(dp),    intent(out) :: R_wf(:)        !! The NR_WF uniform evaluation grid
    real(dp),    intent(in), optional :: rmin  !! Lower bound of eval grid [default RVALS(1)]
    real(dp),    intent(in), optional :: rmax  !! Upper bound of eval grid [default RVALS(N)]

    integer  :: nrvals, nbasis, nbasis_actual, nq, iv
    real(dp) :: R0_use
    real(dp) :: w0(3*kx)
    integer,     allocatable :: idx(:)
    real(dp),    allocatable :: tx(:), xq(:), wq(:)
    complex(dp), allocatable :: H(:,:), S(:,:), energies_full(:), wfcoeffs(:,:)
    complex(dp), allocatable :: W(:,:), Wfull(:,:)
    type(spline_data_type) :: vspline_data

    nrvals = size(rvals, 1)
    call size_check(vvals,    nrvals,        "VVALS")
    call size_check(energies, nwf,           "ENERGIES")
    call size_check(brv,      nwf,           "BRV")
    call size_check(wf,       [nR_wf, nwf],  "WF")

    call build_basis(rvals, vvals, kx, nelem, vspline_data, tx, nbasis, nq, xq, wq)
    call build_eval_grid(rvals, nR_wf, R_wf, rmin, rmax)
    call make_idx(nbasis, ndropl, ndropr, nwf, idx, nbasis_actual)

    ! -- snap R0 to the nearest interior knot
    R0_use = snap2knot(tx, kx, nbasis, R0)

    ! -- complex H and S
    allocate(H(nbasis, nbasis), S(nbasis, nbasis))
    call build_hamiltonian(kx, tx, redmass, jrot, nq, xq, wq, R0_use, theta, vspline_data, w0, H, S)

    ! -- solve the complex generalized problem Hc = ESc in the reduced subspace
    allocate(energies_full(nbasis_actual), wfcoeffs(nbasis_actual, nbasis_actual))
    call solve_schrodinger(H, S, idx, energies_full, wfcoeffs)

    ! -- return the NWF lowest states, evaluated on R_WF
    energies(1:nwf) = energies_full(1:nwf)
    call wf_eval(kx, tx, idx, wfcoeffs(:, 1:nwf), R_wf, wf, R0_use, theta)

    ! -- rotational constant Bv= cv⁺ W cv
    allocate(Wfull(nbasis, nbasis))

    call centrifugal_matrix(kx, tx, redmass, nq, xq, wq, Wfull, R0_use, theta)
    W = Wfull(idx, idx)

    do iv=1, nwf
      ! -- NOT the inner product
      brv(iv) = sum(wfcoeffs(:, iv) * matmul(W, wfcoeffs(:, iv)))
    enddo

  end subroutine rvsolve_ecs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine build_basis(rvals, vvals, kx, nelem, vspline_data, tx, nbasis, nq, xq, wq)
    !! Interpolate the potential, build the clamped knot vector
    !! on the full input grid, and generate the Gauss-Legendre rule.

    use seecs__splines,   only: spline_data_type, build_knots, gauss_legendre
    use seecs__hamilton,  only: Vinterp

    implicit none (type, external)

    real(dp), intent(in)  :: rvals(:), vvals(:)
    integer,  intent(in)  :: kx, nelem
    type(spline_data_type), intent(out) :: vspline_data
      !! The interpolated potential
    real(dp), allocatable, intent(out) :: tx(:)
      !! Knot vector
    integer,  intent(out) :: nbasis
      !! Number of B-spline basis functions
    integer,  intent(out) :: nq
      !! Number of Gauss-Legendre nodes
    real(dp), allocatable, intent(out) :: xq(:), wq(:)
      !! Gauss-Legendre nodes and weights on [-1,1]

    integer :: px, nrvals

    nrvals = size(rvals, 1)
    px     = kx - 1

    ! -- interpolate the potential V(R)
    call Vinterp(rvals, vvals, kx, vspline_data)

    ! -- clamped knots on the full input grid [rvals(1), rvals(nrvals)]
    call build_knots(rvals(1), rvals(nrvals), kx, nelem, tx, nbasis)

    ! -- Gauss-Legendre rule sufficient to integrate the kx-1 degree products exactly
    nq = ceiling(real(3*px + 1, kind=dp)/2._dp) + 1
    allocate(xq(nq), wq(nq))
    call gauss_legendre(nq, xq, wq)

  end subroutine build_basis

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine build_eval_grid(rvals, nR_wf, R_wf, rmin, rmax)
    !! Build the NR_WF-point wavefunction evaluation grid on [RMIN, RMAX]. RMIN/RMAX are optional and
    !! default to the endpoints of the input grid RVALS.

    use seecs__arrays, only: size_check
    use seecs__system, only: die

    implicit none (type, external)

    real(dp), intent(in)  :: rvals(:)
    integer,  intent(in)  :: nR_wf
    real(dp), intent(out) :: R_wf(:)
    real(dp), intent(in), optional :: rmin, rmax

    integer  :: i, nrvals
    real(dp) :: rmin_, rmax_, dR

    nrvals = size(rvals, 1)

    if(nR_wf .lt. 1) call die("NR_WF must be >= 1")
    call size_check(R_wf, nR_wf, "R_WF")

    ! -- resolve the range, defaulting to the input-grid endpoints
    rmin_ = rvals(1)      ; if(present(rmin)) rmin_ = rmin
    rmax_ = rvals(nrvals) ; if(present(rmax)) rmax_ = rmax

    associate(rmin => rmin_, rmax=> rmax_)

      if(rmin .lt. rvals(1))      call die("RMIN lies below the potential grid RVALS(1)")
      if(rmax .gt. rvals(nrvals)) call die("RMAX lies above the potential grid RVALS(NRVALS)")
      if(rmin .gt. rmax) call die("RMIN > RMAX detected, not allowed")

      if(nR_wf .eq. 1) then
        R_wf(1) = rmin
      else
        dR = (rmax - rmin) / real(nR_wf - 1, kind=dp)
        R_wf = [( rmin + real(i - 1, kind=dp)*dR, i = 1, nR_wf )]
      endif

    end associate

  end subroutine build_eval_grid

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine make_idx(nbasis, ndropl, ndropr, nwf, idx, nbasis_actual)
    !! Build the index map enforcing the boundary conditions by dropping NDROPL functions on the
    !! left and NDROPR on the right, and verify that the reduced basis can actually supply the
    !! requested NWF wavefunctions.

    use seecs__system, only: die

    implicit none (type, external)

    integer, intent(in)  :: nbasis, ndropl, ndropr, nwf
    integer, allocatable, intent(out) :: idx(:)
    integer, intent(out) :: nbasis_actual
    integer :: i

    if(nwf .lt. 1) call die("NWF must be ≥ 1")

    idx           = [(i, i = 1 + ndropl, nbasis - ndropr)]
    nbasis_actual = size(idx, 1)

    if(nwf .gt. nbasis_actual) call die("B-spline basis is smaller than the desired number of wavefunctions !")

  end subroutine make_idx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function snap2knot(tx, kx, nbasis, R0) result(res)
    !! Snap the ECS scaling radius R0 to the nearest interior knot so that the e^{iθ} kink lands on an element boundary.

    implicit none (type, external)

    real(dp), intent(in) :: tx(:)
    integer,  intent(in) :: kx, nbasis
    real(dp), intent(in) :: R0
    real(dp) :: res
    integer  :: m

    m = max(kx, min(nbasis + 1, count(tx .le. R0)))
    if(abs(R0 - tx(m)) .le. abs(R0 - tx(m+1))) then
      res = tx(m)
    else
      res = tx(m+1)
    endif

  end function snap2knot

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine centrifugal_matrix(kx, tx, redmass, nq, xq, wq, W, R0, theta)
    !! Build the matrix  W_{ik} = ∫ B_i(R) [1/(2 μ R²)] B_k(R) dR  in the B-spline basis.
    !! This is the rotational term without the j(j+1) factor. The rotational constant of vibrational state v is
    !! then given by B_v = c_v⁺ W c_v. Distortion constants could be obtained with ⟨w|1/(2μR²)|v⟩ = c_wᵀ W c_v,
    !! but we can burn that bridge when we get there.
    !!
    !! When R0 and THETA are both supplied, the ECS contour is used. Because 1/(2μR²) is just a multiplicative operator,
    !! it carries the same e^{iθ} weighting as the potential, and W is complex.
    !!
    !! With neither R0 nor THETA are supplied, the integral is still complex-typed, but real-valued, which is exactly
    !! what one gets when the THETA = 0.

    use seecs__system,    only: die
    use seecs__arrays,    only: size_check
    use seecs__splines,   only: dintrv, dbspvd
    use seecs__constants, only: im

    implicit none (type, external)

    integer,     intent(in)  :: kx             !! B-spline order kx = p+1
    real(dp),    intent(in)  :: tx(:)          !! the knot vector
    real(dp),    intent(in)  :: redmass        !! reduced mass in au
    integer,     intent(in)  :: nq             !! number of Gauss-Legendre nodes per element
    real(dp),    intent(in)  :: xq(nq)         !! Gauss-Legendre nodes on [-1,1]
    real(dp),    intent(in)  :: wq(nq)         !! Gauss-Legendre weights on [-1,1]
    complex(dp), intent(out) :: W(:,:)         !! the centrifugal matrix
    real(dp),    intent(in), optional :: R0    !! ECS scaling radius (must be passed or omitted with THETA)
    real(dp),    intent(in), optional :: theta !! ECS scaling angle in radians (must be passed or omitted with R0)

    logical  :: use_ecs
    integer  :: e, i, k, q, g0, gi, gk
    integer  :: p, nt, nbasis, nelem, ldvnik, lwork, ileft, ilo, mflag, iflag
    real(dp) :: tl, tr, jacobi, x0
    real(dp),    allocatable :: work(:), vnikx(:,:)
    complex(dp) :: z, crot, expo

    p      = kx - 1
    nt     = size(tx, 1)
    nbasis = nt - kx
    nelem  = nt - 2*kx + 1
    ldvnik = kx
    lwork  = kx + 1
    ilo    = 1

    call size_check(W, [nbasis, nbasis], "W")
    W = 0

    use_ecs = present(R0) .AND. present(theta)
    if(present(R0) .neqv. present(theta)) call die("either pass both R0 and THETA, or neither")

    allocate(vnikx(ldvnik, 1))   ! <-- only need B-spline values, no derivatives
    allocate(work((lwork*(lwork+1))/2))

    elements: do e = 1, nelem

      tl = tx(p+e) ; tr = tx(p+e+1)
      jacobi = (tr - tl)/2.0_dp

      nodes: do q = 1, nq

        x0 = (tr+tl)/2.0_dp + jacobi*xq(q)

        if(use_ecs) then
          if(x0 .gt. R0) then
            expo = exp(im*theta)
            z    = R0 + (x0 - R0)*expo
          else
            expo = 1
            z    = x0
          endif
        else
          expo = 1
          z    = x0
        endif

        call dintrv(tx, nt, x0, ilo, ileft, mflag, extrap = .false.)
        if(ileft .lt. kx) call die("ileft < kx")

        call dbspvd(tx, kx, 1, x0, ileft, ldvnik, vnikx, work, iflag)
        if(iflag .ne. 0) call die("dbspvd failed")

        crot = 1._dp / (2*redmass*z*z)

        g0 = ileft - kx
        do k = 1, ldvnik
          gk = g0 + k
          do i = 1, ldvnik
            gi = g0 + i
            W(gi, gk) = W(gi, gk) + wq(q)*jacobi * expo * vnikx(i,1)*crot*vnikx(k,1)
          enddo
        enddo

      enddo nodes
    enddo elements

    deallocate(vnikx, work)

  end subroutine centrifugal_matrix

! ================================================================================================================================ !
end module seecs__drivers
! ================================================================================================================================ !
