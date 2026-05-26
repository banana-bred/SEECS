! ================================================================================================================================ !
module seecs__linalg

  use seecs__kinds, only: dp
#ifndef USE_EXTERNAL_LAPACK
  use stdlib_linalg_lapack, only: dsygv => stdlib_dsygv, dstev => stdlib_dstev, zggev => stdlib_zggev
#endif

  implicit none (type, external)

  private

  public :: dsygv
  public :: dstev
  public :: zggev

#ifdef USE_EXTERNAL_LAPACK
  interface
    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      import dp
      implicit none (type, external)
      integer   :: itype
      character :: jobz
      character :: uplo
      integer   :: n
      real(dp)  :: a(lda,*)
      integer   :: lda
      real(dp)  :: b(ldb,*)
      integer   :: ldb
      real(dp)  :: w(*)
      real(dp)  :: work(*)
      integer   :: lwork
      integer   :: info
    end subroutine dsygv
  end interface
  interface
    subroutine dstev(jobz, n, d, e, z, ldz, work, info)
      import dp
      implicit none (type, external)
      character :: jobz
      integer   :: n
      real(dp)  :: d(*)
      real(dp)  :: e(*)
      real(dp)  :: z(ldz,*)
      integer   :: ldz
      real(dp)  :: work(*)
      integer   :: info
    end subroutine dstev
  end interface
  interface
    subroutine zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      import dp
      implicit none (type, external)
      character   :: jobvl
      character   :: jobvr
      integer     :: n
      complex(dp) :: a(lda, *)
      integer     :: lda
      complex(dp) :: b(ldb, *)
      integer     :: ldb
      complex(dp) :: alpha(*)
      complex(dp) :: beta(*)
      complex(dp) :: vl(ldvl, *)
      integer     :: ldvl
      complex(dp) :: vr(ldvr, *)
      integer     :: ldvr
      complex(dp) :: work(*)
      integer     :: lwork
      real(dp)    :: rwork(*)
      integer     :: info
    end subroutine zggev
  end interface
#endif

! ================================================================================================================================ !
end module seecs__linalg
! ================================================================================================================================ !
