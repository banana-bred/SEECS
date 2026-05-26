! ================================================================================================================================ !
program seecs
  !! solve the (rovibrational) Schroedinger Equation with optional Exterior Complex Scaling

  use seecs__kinds,      only: dp
  use seecs__system,     only: die, stdin, stderr
  use seecs__drivers,    only: rvsolve
  use seecs__splines,    only: gauss_legendre, build_knots
  use seecs__units,      only: convert_energy, convert_length, convert_mass, convert_brot

  implicit none (type, external)

  type config_type
    !! Global configuration type
    integer  :: k, nelem, ndropl, ndropr, jrot, nwf, nR_wf
    logical  :: do_ecs
    real(dp) :: redmass, theta, R0, rmin, rmax
    character(:), allocatable :: potential_input_file
    character(:), allocatable :: energies_output_file, wfs_output_file
    character(:), allocatable :: munits_in, eunits_in, runits_in
    character(:), allocatable :: bunits_out, eunits_out, runits_out
  end type config_type

  integer :: nrvals
  real(dp),    allocatable :: rvals(:), vvals(:), R_wf(:), bvrot(:)
  real(dp),    allocatable :: energies_real(:), brv_real(:)
  real(dp),    allocatable :: wfs_real(:,:)
  complex(dp), allocatable :: energies_cmplx(:), brv_cmplx(:)
  complex(dp), allocatable :: wfs_cmplx(:,:)

  type(config_type) :: cfg

  call read_namelist(stdin, cfg)

  call read_potential(cfg%potential_input_file, rvals, vvals)

  nrvals = size(rvals, 1)

  ! -- convert inputs to au
  call convert_mass(cfg%redmass, cfg%munits_in, "au")
  call convert_length(rvals,     cfg%runits_in, "au")
  call convert_energy(vvals,     cfg%eunits_in, "au")
  if(cfg%do_ecs) call convert_length(cfg%R0, cfg%runits_in, "au")
  ! -- resolve rmin/rmax, convert if  necessary
  if(cfg%rmin .lt. 0.0_dp) then ; cfg%rmin = rvals(1)      ; else ; call convert_length(cfg%rmin, cfg%runits_in, "au") ; endif
  if(cfg%rmax .lt. 0.0_dp) then ; cfg%rmax = rvals(nrvals) ; else ; call convert_length(cfg%rmax, cfg%runits_in, "au") ; endif

  allocate(R_wf(cfg%nR_wf), source=0.0_dp)

  if(cfg%do_ecs) then

    allocate(wfs_cmplx(cfg%nR_wf, cfg%nwf), energies_cmplx(cfg%nwf), brv_cmplx(cfg%nwf), source=(0.0_dp, 0.0_dp))

    call rvsolve(      &
        rvals          &
      , vvals          &
      , cfg%jrot       &
      , cfg%k          &
      , cfg%nelem      &
      , cfg%redmass    &
      , cfg%ndropl     &
      , cfg%ndropr     &
      , cfg%nwf        &
      , cfg%R0         &
      , cfg%theta      &
      , wfs_cmplx      &
      , energies_cmplx &
      , brv_cmplx      &
      , cfg%nR_wf      &
      , R_wf           &
      , cfg%rmin       &
      , cfg%rmax       &
    )

  else

    allocate(wfs_real(cfg%nR_wf, cfg%nwf), energies_real(cfg%nwf), brv_real(cfg%nwf), source = 0.0_dp)

    call rvsolve(     &
        rvals         &
      , vvals         &
      , cfg%jrot      &
      , cfg%k         &
      , cfg%nelem     &
      , cfg%redmass   &
      , cfg%ndropl    &
      , cfg%ndropr    &
      , cfg%nwf       &
      , wfs_real      &
      , energies_real &
      , brv_real      &
      , cfg%nR_wf     &
      , R_wf          &
      , cfg%rmin      &
      , cfg%rmax      &
    )

  endif

  call convert_length(R_wf, "au", cfg%runits_out)

  if(cfg%do_ecs) then
    call convert_energy(energies_cmplx, "au", cfg%eunits_out)
    call convert_brot  (brv_cmplx,      "au", cfg%bunits_out)
    call write_energies_cmplx(cfg%energies_output_file, energies_cmplx, brv_cmplx, cfg%jrot, cfg%eunits_out, cfg%bunits_out)
    call write_wfs_cmplx(cfg%wfs_output_file, R_wf, wfs_cmplx, cfg%runits_out)
  else
    call convert_energy(energies_real, "au", cfg%eunits_out)
    call convert_brot  (brv_real,      "au", cfg%bunits_out)
    call write_energies_real(cfg%energies_output_file, energies_real, brv_real, cfg%jrot, cfg%eunits_out, cfg%bunits_out)
    call write_wfs_real(cfg%wfs_output_file, R_wf, wfs_real, cfg%runits_out)
  endif

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  subroutine read_namelist(funit, cfg)
    !! Reads the CONTROL namelist (from FUNIT, e.g. stdin) that controls program execution. All physical-unit
    !! selectors must be specified explicitly; there are no implicit unit defaults. THETA is read in degrees
    !! and converted to radians here. RMIN/RMAX are optional and left at a negative sentinel when omitted, to
    !! be resolved against the input grid by MAIN. They are in the same units as the input R grid

    use seecs__constants,  only: pi
    use seecs__system,     only: stderr, die
    use seecs__characters, only: to_lower

    implicit none (type, external)

    integer,           intent(in)  :: funit !! File unit from which to read the namelist
    type(config_type), intent(out) :: cfg   !! Global configuration type to return

    integer,  parameter :: DEFAULT_INT  = -1
    real(dp), parameter :: DEFAULT_REAL = -huge(1._dp)
    integer,  parameter :: CHARLEN = 1000 ! -- default character length
    integer,  parameter :: ULEN = 10     ! -- spechial unit-selector buffer length

    logical :: flag
    integer :: io

    ! -- namelist variables (with defaults; unit selectors deliberately have NO default)
    integer  :: k       = DEFAULT_INT
    integer  :: nelem   = DEFAULT_INT
    integer  :: ndropl  = DEFAULT_INT
    integer  :: ndropr  = DEFAULT_INT
    integer  :: jrot    = DEFAULT_INT
    integer  :: nwf     = DEFAULT_INT
    integer  :: nR_wf   = DEFAULT_INT
    logical  :: do_ecs  = .false.
    real(dp) :: redmass = DEFAULT_REAL
    real(dp) :: theta   = DEFAULT_REAL
    real(dp) :: R0      = DEFAULT_REAL
    real(dp) :: rmin    = DEFAULT_REAL   ! -- optional; default (<0) is to use grid lower endpoint
    real(dp) :: rmax    = DEFAULT_REAL   ! -- optional; default (<0) is to use grid upper endpoint

    character(CHARLEN) :: potential_input_file = ""
    character(CHARLEN) :: energies_output_file = ""
    character(CHARLEN) :: wfs_output_file      = ""

    character(ULEN) :: runits_in  = ""    ! -- distance units of the potential file R column (and R0/RMIN/RMAX)
    character(ULEN) :: eunits_in  = ""    ! -- energy units of the potential file V column
    character(ULEN) :: munits_in  = ""    ! -- units of REDMASS
    character(ULEN) :: runits_out = ""    ! -- distance units for the R_WF grid in the wavefunction file
    character(ULEN) :: eunits_out = ""    ! -- energy units for the energies file
    character(ULEN) :: bunits_out = ""    ! -- units for the rotational constants

    namelist / control /                                            &
        k, nelem, ndropl, ndropr, jrot, nwf, nR_wf, redmass         &
      , do_ecs, R0, theta                                           &
      , rmin, rmax                                                  &
      , potential_input_file, energies_output_file, wfs_output_file &
      , runits_in,  eunits_in,  munits_in                           &
      , runits_out, eunits_out, bunits_out

    read(funit, nml=control, iostat=io)
    if(io .ne. 0) call die("Error reading the CONTROL namelist")

    flag = .false.

    if(k       .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the B-spline order K" ; endif
    if(nelem   .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the number of elements NELEM" ; endif
    if(ndropl  .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the left BC drop NDROPL" ; endif
    if(ndropr  .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the right BC drop NDROPR" ; endif
    if(jrot    .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the rotational number JROT" ; endif
    if(nwf     .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the number of wavefunctions NWF" ; endif
    if(nR_wf   .eq. DEFAULT_INT)  then ; flag=.true. ; write(stderr,'(A)') "Please define the number of eval points NR_WF" ; endif
    if(redmass .eq. DEFAULT_REAL) then ; flag=.true. ; write(stderr,'(A)') "Please define the reduced mass REDMASS" ; endif

    ! -- ECS only
    if(do_ecs .AND. theta .eq. DEFAULT_REAL) then ; flag=.true. ; write(stderr,'(A)') "Please define the scaling angle THETA (deg)" ; endif
    if(do_ecs .AND. R0    .eq. DEFAULT_REAL) then ; flag=.true. ; write(stderr,'(A)') "Please define the scaling radius R0" ; endif

    ! -- files
    if(len_trim(potential_input_file) .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define POTENTIAL_INPUT_FILE" ; endif
    if(len_trim(energies_output_file) .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define ENERGIES_OUTPUT_FILE" ; endif
    if(len_trim(wfs_output_file)      .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define WFS_OUTPUT_FILE" ; endif

    ! -- units
    if(len_trim(runits_in)  .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define input distance units RUNITS_IN" ; endif
    if(len_trim(eunits_in)  .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define input energy units EUNITS_IN" ; endif
    if(len_trim(munits_in)  .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define input mass units MUNITS_IN" ; endif
    if(len_trim(runits_out) .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define output distance units RUNITS_OUT" ; endif
    if(len_trim(eunits_out) .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define output energy units EUNITS_OUT" ; endif
    if(len_trim(bunits_out) .eq. 0) then ; flag=.true. ; write(stderr,'(A)') "Please define output B units BUNITS_OUT" ; endif

    if(flag) then
      write(stderr, *)
      call die("Critically undefined namelist variables detected. See above messages")
    endif

    call to_lower(runits_in)  ; call to_lower(eunits_in)  ; call to_lower(munits_in)
    call to_lower(runits_out) ; call to_lower(eunits_out) ; call to_lower(bunits_out)

    call check_runit(runits_in,  "RUNITS_IN") ; call check_runit(runits_out, "RUNITS_OUT")
    call check_eunit(eunits_in,  "EUNITS_IN") ; call check_eunit(eunits_out, "EUNITS_OUT")
    call check_munit(munits_in,  "MUNITS_IN") ; call check_bunit(bunits_out, "BUNITS_OUT")

    ! -- R-grid verification
    if(rmin .ne. DEFAULT_REAL .AND. rmax .ne. DEFAULT_REAL) then
      if(rmin .ge. rmax) call die("RMIN must be < RMAX")
    endif

    ! -- θ: ° -> rads
    if(do_ecs) theta = theta * pi / 180._dp

    ! -- populate the config
    cfg % k       = k
    cfg % nelem   = nelem
    cfg % ndropl  = ndropl
    cfg % ndropr  = ndropr
    cfg % jrot    = jrot
    cfg % nwf     = nwf
    cfg % nR_wf   = nR_wf
    cfg % do_ecs  = do_ecs
    cfg % redmass = redmass
    cfg % theta   = theta
    cfg % R0      = R0
    cfg % rmin    = rmin
    cfg % rmax    = rmax

    cfg % potential_input_file = trim(potential_input_file)
    cfg % energies_output_file = trim(energies_output_file)
    cfg % wfs_output_file      = trim(wfs_output_file)

    cfg % runits_in  = trim(runits_in)
    cfg % eunits_in  = trim(eunits_in)
    cfg % munits_in  = trim(munits_in)
    cfg % runits_out = trim(runits_out)
    cfg % eunits_out = trim(eunits_out)
    cfg % bunits_out = trim(bunits_out)

  end subroutine read_namelist

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine check_runit(u, name)
    character(*), intent(in) :: u, name
    select case(trim(u))
    case("au", "bohr", "angstrom")
      ! ok
    case default
      call die(name//" must be 'au', 'bohr', or 'angstrom' (got '"//trim(u)//"')")
    end select
  end subroutine check_runit

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine check_eunit(u, name)
    character(*), intent(in) :: u, name
    select case(trim(u))
    case("hartree", "au", "ev", "invcm", "rydberg")
      ! ok
    case default
      call die(name//" must be 'hartree', 'au', 'ev', 'invcm', or 'rydberg' (got '"//trim(u)//"')")
    end select
  end subroutine check_eunit

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine check_munit(u, name)
    character(*), intent(in) :: u, name
    select case(trim(u))
    case("amu", "au")
      ! ok
    case default
      call die(name//" must be 'amu', or 'au' (got '"//trim(u)//"')")
    end select
  end subroutine check_munit

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine check_bunit(u, name)
    character(*), intent(in) :: u, name
    select case(trim(u))
    case("invcm", "mhz", "hartree", "au")
      ! ok
    case default
      call die(name//" must be 'invcm', 'mhz', 'au', or 'hartree' (got '"//trim(u)//"')")
    end select
  end subroutine check_bunit

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_potential(filename, rvals, vvals)
    !! Read the potential at FILENAME
    use iso_fortran_env, only: iostat_end

    implicit none (type, external)

    character(*), intent(in) :: filename
      !! The potential lives here as R V columns (space-separated)
    real(dp), intent(out), allocatable :: rvals(:)
      !! R-grid
    real(dp), intent(out), allocatable :: vvals(:)
      !! V(R)

    integer, parameter :: iostat_ok = 0

    integer :: funit, io

    inquire(file=filename, iostat=io)
    if(io .ne. iostat_ok) call die("Problem reading potential file " // filename)

    call read_commented_file(filename, rvals, vvals, '#')

  end subroutine read_potential

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_commented_file(fname, a, b, comment_char_in)
    !! Read a file's contents into arrays a and b, with the capability to skip lines that are properly commented.
    !! Refer to the inerface definition for an explanation
    !! The comment character defaults to "!", but can be set to anything not in the character 'numeric' from the 'characters' module and
    !! also not whitespace. This version reads real numbers.

    use iso_fortran_env,   only: iostat_end
    use seecs__characters, only: numeric

    implicit none (type, external)

    real(dp), intent(inout), allocatable :: a(:)
    real(dp), intent(inout), allocatable :: b(:)

    character(*), intent(in) :: fname

    character(1), intent(in), optional :: comment_char_in

    real(dp) :: aElement
    real(dp) :: bElement

    character(1) :: comment_char

    integer :: count
    integer :: inunit
    integer :: io
    integer :: commentStart
    integer :: numericStart
    integer :: numericEnd

    character(1000) :: line

    ! -- default comment character
    comment_char = "!"

    ! -- don't know a priori how much we will read, build arrays as we go. Make sure' they're deallocated
    if(allocated(a)) deallocate(a)
    if(allocated(b)) deallocate(b)

    ! -- set different comment character maybe
    if(present(comment_char_in)) comment_char = comment_char_in

    open(newunit = inunit, file = fname, action="read")

    count = 0
    do
      read(inunit, "(A)", iostat = io) line
      if(io .eq. iostat_end) exit
      if(io .ne. 0) call die("Problem reading data from file " // "'" // fname // "'")
      commentStart = scan(line, comment_char)    ! -- position of first comment character
      numericStart = scan(line, numeric)         ! -- position of first numeric character
      numericEnd   = scan(line, numeric, .true.) ! -- position of last  numeric character
      ! -- remove comments that appear after all numeric chars
      if(commentStart .gt. numericEnd) line = line(numericStart:numericEnd)
      ! -- cycle reading if the line appears commented out
      if(commentStart .gt. 0 .AND. commentStart .lt. numericStart) cycle
      read(line, *) aElement, bElement
      count = count + 1
      ! call append(a, aElement)
      ! call append(b, bElement)
    enddo

    rewind(inunit)
    allocate(a(count), b(count))
    count = 0

    do
      read(inunit, "(A)", iostat = io) line
      if(io .eq. iostat_end) exit
      if(io .ne. 0) then
        write(stderr, '("IOSTAT = ", I0)') io
        call die("Problem reading data from file " // "'" // fname // "'")
      endif
      commentStart = scan(line, comment_char)    ! -- position of first comment character
      numericStart = scan(line, numeric)         ! -- position of first numeric character
      numericEnd   = scan(line, numeric, .true.) ! -- position of last  numeric character
      ! -- remove comments that appear after all numeric chars
      if(commentStart .gt. numericEnd) line = line(numericStart:numericEnd)
      ! -- cycle reading if the line appears commented out
      if(commentStart .gt. 0 .AND. commentStart .lt. numericStart) cycle
      read(line, *) aElement, bElement
      count = count + 1
      a(count) = aElement
      b(count) = bElement
    enddo

  end subroutine read_commented_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_energies_real(filename, energies, brv, jrot, eunits, bunits)
    use seecs__arrays, only: size_check
    implicit none (type, external)
    character(*), intent(in) :: filename
    real(dp),     intent(in) :: energies(:), brv(:)
    integer,      intent(in) :: jrot
    character(*), intent(in) :: eunits, bunits
    integer :: funit, iv, v, n, io
    n = size(energies, 1)
    call size_check(brv, n, "BRV")
    open(newunit=funit, file=filename, action="write", status="replace", iostat=io)
    if(io .ne. 0) then
      write(stderr, '("IOSTAT: ", I0)') io
      call die("Trouble opening energies output file at "//filename//". Please ensure its parent directory exists.")
    endif
    write(funit, '(A)')     "# SEECS (ro)vibrational energies (no ECS)"
    write(funit, '(A)')     "# energy units: "//trim(eunits)//"  |  B units: "//trim(bunits)
    write(funit, '(A, I0)') "# rotational quantum number: j = ", jrot
    write(funit, '("# ", A6, 2(1X, A26) )') "v", "E", "B"
    do iv = 1, n
      v = iv - 1
      write(funit, '(2X, I6, 2(1X, ES26.16E3) )') v, energies(iv), brv(iv)
    enddo
    close(funit)
  end subroutine write_energies_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_energies_cmplx(filename, energies, brv, jrot, eunits, bunits)
    use seecs__arrays, only: size_check
    implicit none (type, external)
    character(*), intent(in) :: filename
    complex(dp),  intent(in) :: energies(:), brv(:)
    integer,      intent(in) :: jrot
    character(*), intent(in) :: eunits, bunits
    integer :: funit, iv, v, n, io
    n = size(energies, 1)
    call size_check(brv, n, "BRV")
    open(newunit=funit, file=filename, action="write", status="replace", iostat=io)
    if(io .ne. 0) call die("Had trouble opening energies output file at "//filename)
    write(funit, '(A)')     "# SEECS (ro)vibrational energies (no ECS)"
    write(funit, '(A)')     "# energy units: "//trim(eunits)//"  |  B units: "//trim(bunits)
    write(funit, '(A, I0)') "# rotational quantum number: j = ", jrot
    write(funit, '("# ", A6, 4(1X, A26) )') "v", "Re(E)", "Im(E)", "Re(B)", "Im(B)"
    do iv = 1, n
      v = iv - 1
      write(funit, '(2X, I6, 4(1X, ES26.16E3) )') v, energies(iv), brv(iv)
    enddo
    close(funit)
  end subroutine write_energies_cmplx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_wfs_real(filename, R_wf, wfs, runits)
    use seecs__arrays, only: size_check
    implicit none (type, external)
    character(*), intent(in) :: filename
    real(dp),     intent(in) :: R_wf(:)
    real(dp),  intent(in) :: wfs(:,:)
    character(*), intent(in) :: runits
    integer :: funit, iR, iwf, nR, nwf, io
    nR  = size(R_wf, 1)
    nwf = size(wfs, 2)
    call size_check(wfs, [nR, nwf], "WFS")
    open(newunit=funit, file=filename, action="write", status="replace", iostat=io)
    if(io .ne. 0) call die("Had trouble opening wavefunction output file at "//filename)
    write(funit, '(A)') "# SEECS wavefunctions (ECS)"
    write(funit, '(A)') "# distace units: "//trim(runits)//"  |  Ψ(R) in au"
    write(funit, '("# ", 3A26, " ...")') "R", "ψ₀(R)", "ψ₁(R)"
    do iR = 1, nR
      write(funit, '(2X, *(ES26.16E3, 1X))') R_wf(iR), ( wfs(iR, iwf), iwf = 1, nwf )
    enddo
    close(funit)
  end subroutine write_wfs_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_wfs_cmplx(filename, R_wf, wfs, runits)
    use seecs__arrays, only: size_check
    implicit none (type, external)
    character(*), intent(in) :: filename
    real(dp),     intent(in) :: R_wf(:)
    complex(dp),  intent(in) :: wfs(:,:)
    character(*), intent(in) :: runits
    integer :: funit, iR, iwf, nR, nwf, io
    nR  = size(R_wf, 1)
    nwf = size(wfs, 2)
    call size_check(wfs, [nR, nwf], "WFS")
    open(newunit=funit, file=filename, action="write", status="replace", iostat=io)
    if(io .ne. 0) call die("Had trouble opening wavefunction output file at "//filename)
    write(funit, '(A)') "# SEECS wavefunctions (ECS)"
    write(funit, '(A)') "# distace units: "//trim(runits)//"  |  Ψ(R) in au"
    write(funit, '("# ", 5A26, " ...")') "R", "Re[ψ₀(R)]", "Im[ψ₀(R)]", "Re[ψ₁(R)]", "Im[ψ₁(R)]"
    do iR = 1, nR
      write(funit, '(2X, *(ES26.16E3, 1X))') R_wf(iR), ( wfs(iR, iwf), iwf = 1, nwf )
    enddo
    close(funit)
  end subroutine write_wfs_cmplx

! ================================================================================================================================ !
end program seecs
! ================================================================================================================================ !
