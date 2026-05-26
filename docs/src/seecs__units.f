! ================================================================================================================================ !
module seecs__units

  use seecs__kinds,      only: dp
  use seecs__system,     only: die
  use seecs__characters, only: lower

  implicit none (type, external)

  private

  public :: convert_length
  public :: convert_mass
  public :: convert_energy
  public :: convert_brot

  interface convert_energy
    module procedure convert_energy_real
    module procedure convert_energy_cmplx
  end interface convert_energy

  interface convert_brot
    module procedure convert_brot_real
    module procedure convert_brot_cmplx
  end interface convert_brot


  real(dp), parameter, public :: au2invcm = 219474.6313710_dp
    !! Atomic units of energy -> cm¯¹
  real(dp), parameter, public :: au2ev = 27.2113834e0_dp
    !! Atomic units of energy -> electron volts
  real(dp), parameter, public :: au2amu = 5.4858010860603975e-4_dp
    !! Atomic units of (electron) mass -> atomic mass units (Daltons)
  real(dp), parameter, public :: au2ang = 0.5291772083_dp
    !! Atomic units of legth (Bohr) -> Angstrom
  real(dp), parameter, public :: au2mhz = 6579683920.5018206_dp
    !! Atomic units of energy -> Mhz

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_length(length, units_in, units_out)
    implicit none (type, external)
    real(dp),     intent(inout) :: length
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au", "bohr")
      continue
    case("angstrom")
      length = length / au2ang
    case default
      call die("Unacceptable input length unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au", "bohr")
      continue
    case("angstrom")
      length = length * au2ang
    case default
      call die("Unacceptable output length unit: "//units_out)
    end select
  end subroutine convert_length

  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_energy_real(energy, units_in, units_out)
    implicit none (type, external)
    real(dp),     intent(inout) :: energy
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au", "hartree")
      continue
    case("invcm")
      energy = energy / au2invcm
    case("ev")
      energy = energy / au2ev
    case("rydberg")
      energy = energy / 2.0_dp
    case default
      call die("Unacceptable input energy unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au", "hartree")
      continue
    case("invcm")
      energy = energy * au2invcm
    case("ev")
      energy = energy * au2ev
    case("rydberg")
      energy = energy * 2.0_dp
    case default
      call die("Unacceptable output energy unit: "//units_out)
    end select
  end subroutine convert_energy_real
  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_energy_cmplx(energy, units_in, units_out)
    implicit none (type, external)
    complex(dp),     intent(inout) :: energy
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au", "hartree")
      continue
    case("invcm")
      energy = energy / au2invcm
    case("ev")
      energy = energy / au2ev
    case("rydberg")
      energy = energy / 2.0_dp
    case default
      call die("Unacceptable input energy unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au", "hartree")
      continue
    case("invcm")
      energy = energy * au2invcm
    case("ev")
      energy = energy * au2ev
    case("rydberg")
      energy = energy * 2.0_dp
    case default
      call die("Unacceptable output energy unit: "//units_out)
    end select
  end subroutine convert_energy_cmplx

  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_brot_real(brot, units_in, units_out)
    implicit none (type, external)
    real(dp),     intent(inout) :: brot
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au", "hartree")
      continue
    case("invcm")
      brot = brot / au2invcm
    case("mhz")
      brot = brot / au2mhz
    case default
      call die("Unacceptable input brot unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au", "hartree")
      continue
    case("invcm")
      brot = brot * au2invcm
    case("mhz")
      brot = brot * au2mhz
    case default
      call die("Unacceptable output brot unit: "//units_out)
    end select
  end subroutine convert_brot_real
  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_brot_cmplx(brot, units_in, units_out)
    implicit none (type, external)
    complex(dp),     intent(inout) :: brot
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au", "hartree")
      continue
    case("invcm")
      brot = brot / au2invcm
    case("mhz")
      brot = brot / au2mhz
    case default
      call die("Unacceptable input brot unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au", "hartree")
      continue
    case("invcm")
      brot = brot * au2invcm
    case("mhz")
      brot = brot * au2mhz
    case default
      call die("Unacceptable output brot unit: "//units_out)
    end select
  end subroutine convert_brot_cmplx

  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental subroutine convert_mass(mass, units_in, units_out)
    implicit none (type, external)
    real(dp),     intent(inout) :: mass
    character(*), intent(in)    :: units_in, units_out
    ! -- in
    select case(lower(units_in))
    case("au")
      continue
    case("amu")
      mass = mass / au2amu
    case default
      call die("Unacceptable input mass unit: "//units_in)
    end select
    ! -- out
    select case(lower(units_out))
    case("au")
      continue
    case("amu")
      mass = mass * au2amu
    case default
      call die("Unacceptable output mass unit: "//units_out)
    end select
  end subroutine convert_mass

! ================================================================================================================================ !
end module seecs__units
! ================================================================================================================================ !
