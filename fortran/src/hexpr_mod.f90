! hexpr_mod.f90 -- Fortran 95 + iso_c_binding interface to libhexpr
!
! Naming convention:
!   c_hexpr_*  -- bind(C) interface declarations (private, exact C ABI)
!   hexpr_*    -- public Fortran wrappers (thin layer over c_hexpr_*)
!
! Version constant naming note:
!   Fortran is case-insensitive.  The C macros HEXPR_VERSION_MAJOR etc.
!   would clash with the runtime functions hexpr_version_major() etc.
!   The parameter constants are therefore spelled HEXPR_VMAJOR /
!   HEXPR_VMINOR / HEXPR_VPATCH (short form, no clash).
module hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! ----------------------------------------------------------------
  ! Status codes  (hexpr_status_t)
  ! ----------------------------------------------------------------
  integer(c_int), parameter, public :: HEXPR_OK              = 0_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_INVALID_ARG = 1_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_OOM         = 2_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_IO          = 3_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_RANGE       = 4_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_NOT_BUILT   = 5_c_int
  integer(c_int), parameter, public :: HEXPR_ERR_INTERNAL    = 99_c_int

  ! ----------------------------------------------------------------
  ! DRT kind  (hexpr_drt_kind_t)
  ! ----------------------------------------------------------------
  integer(c_int), parameter, public :: HEXPR_DRT_SOCI   = 0_c_int
  integer(c_int), parameter, public :: HEXPR_DRT_FOCI   = 1_c_int
  integer(c_int), parameter, public :: HEXPR_DRT_CSFSET = 2_c_int

  ! ----------------------------------------------------------------
  ! Expression selector  (hexpr_expr_which_t)
  ! ----------------------------------------------------------------
  integer(c_int), parameter, public :: HEXPR_EXPR_FULL = 0_c_int
  integer(c_int), parameter, public :: HEXPR_EXPR_ONE  = 1_c_int

  ! ----------------------------------------------------------------
  ! CSF step codes  (HEXPR_STEP_*)
  ! ----------------------------------------------------------------
  integer(c_int), parameter, public :: HEXPR_STEP_E = 0_c_int
  integer(c_int), parameter, public :: HEXPR_STEP_U = 1_c_int
  integer(c_int), parameter, public :: HEXPR_STEP_D = 2_c_int
  integer(c_int), parameter, public :: HEXPR_STEP_F = 3_c_int

  ! ----------------------------------------------------------------
  ! Compile-time version constants
  ! Mirrors HEXPR_VERSION_MAJOR/MINOR/PATCH macros in hexpr.h.
  ! Short-form names (VMAJOR/VMINOR/VPATCH) avoid the case-insensitive
  ! clash with the runtime functions hexpr_version_major/minor/patch.
  ! ----------------------------------------------------------------
  integer(c_int), parameter, public :: HEXPR_VMAJOR = 1_c_int
  integer(c_int), parameter, public :: HEXPR_VMINOR = 0_c_int
  integer(c_int), parameter, public :: HEXPR_VPATCH = 0_c_int

  ! ----------------------------------------------------------------
  ! Module save: version string cache
  !
  ! hexpr_version_string() (C) returns a static const char * that
  ! lives for the library lifetime.  We copy it once into a fixed-
  ! length Fortran character so callers never hold a C pointer.
  ! Fixed-length character (no deferred) is required for F95.
  ! ----------------------------------------------------------------
  character(len=64), save :: cached_version = ''
  logical,           save :: version_cached  = .false.

  ! ----------------------------------------------------------------
  ! C stdout FILE* for the *_show* wrappers.
  !
  ! The hexpr_*_show* C functions take a FILE * stream.  Aliasing the
  ! libc "stdout" data symbol from Fortran is unreliable (copy
  ! relocation produces a zero-valued local definition), so we obtain
  ! the real FILE* through the C accessor hexpr_fortran_stdout().  The
  ! public *_show* wrappers always emit to stdout and take no stream
  ! argument.
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! C interface block
  !
  ! All names below are private (module default).
  ! bind(C, name="...") uses the exact symbol from hexpr.h.
  ! Opaque handles are type(c_ptr); VALUE matches C pointer-by-value.
  ! ----------------------------------------------------------------
  interface

    ! ---- Version -------------------------------------------------

    function c_hexpr_version_string() &
        bind(C, name="hexpr_version_string") result(ptr)
      import :: c_ptr
      implicit none
      type(c_ptr) :: ptr          ! const char * (static C string)
    end function c_hexpr_version_string

    function c_hexpr_version_major() &
        bind(C, name="hexpr_version_major") result(v)
      import :: c_int
      implicit none
      integer(c_int) :: v
    end function c_hexpr_version_major

    function c_hexpr_version_minor() &
        bind(C, name="hexpr_version_minor") result(v)
      import :: c_int
      implicit none
      integer(c_int) :: v
    end function c_hexpr_version_minor

    function c_hexpr_version_patch() &
        bind(C, name="hexpr_version_patch") result(v)
      import :: c_int
      implicit none
      integer(c_int) :: v
    end function c_hexpr_version_patch

    ! ---- Error ---------------------------------------------------

    function c_hexpr_last_error() &
        bind(C, name="hexpr_last_error") result(ptr)
      import :: c_ptr
      implicit none
      ! Thread-local; valid only until the next hexpr call on this
      ! thread.  The public wrapper copies it immediately -- do NOT
      ! cache this pointer.
      type(c_ptr) :: ptr          ! const char *
    end function c_hexpr_last_error

    ! ---- Lifecycle -----------------------------------------------

    function c_hexpr_init() &
        bind(C, name="hexpr_init") result(status)
      import :: c_int
      implicit none
      integer(c_int) :: status    ! hexpr_status_t
    end function c_hexpr_init

    subroutine c_hexpr_shutdown() &
        bind(C, name="hexpr_shutdown")
      implicit none
    end subroutine c_hexpr_shutdown

    ! ---- DRT lifecycle -------------------------------------------

    function c_hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next) &
        bind(C, name="hexpr_drt_create") result(handle)
      import :: c_ptr, c_int
      implicit none
      integer(c_int), value :: kind     ! hexpr_drt_kind_t
      integer(c_int), value :: nve
      integer(c_int), value :: spin_x2
      integer(c_int), value :: ncore
      integer(c_int), value :: nval
      integer(c_int), value :: next
      type(c_ptr) :: handle             ! hexpr_drt_t *, NULL on failure
    end function c_hexpr_drt_create

    ! C: void hexpr_drt_destroy(hexpr_drt_t *drt)
    ! Pointer passed by value; the public wrapper nullifies the
    ! Fortran handle with c_null_ptr after this call.
    subroutine c_hexpr_drt_destroy(drt) &
        bind(C, name="hexpr_drt_destroy")
      import :: c_ptr
      implicit none
      type(c_ptr), value :: drt
    end subroutine c_hexpr_drt_destroy

    ! C: hexpr_drt_t *hexpr_drt_from_csfset(int norb,
    !                                        const unsigned char *csfs,
    !                                        int n_csfs)
    ! csfs is a packed buffer of n_csfs * norb bytes (step codes 0..3).
    ! The public wrapper passes c_loc(csfs) for the byte-buffer argument,
    ! matching the array-pointer convention used by c_hexpr_expr_read.
    function c_hexpr_drt_from_csfset(norb, csfs, n_csfs) &
        bind(C, name="hexpr_drt_from_csfset") result(p)
      import :: c_int, c_ptr
      implicit none
      integer(c_int), value :: norb
      type(c_ptr),    value :: csfs     ! const unsigned char *
      integer(c_int), value :: n_csfs
      type(c_ptr)           :: p        ! hexpr_drt_t *, NULL on failure
    end function c_hexpr_drt_from_csfset

    ! ---- DRT accessors -------------------------------------------

    function c_hexpr_drt_norb(drt) &
        bind(C, name="hexpr_drt_norb") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_norb

    function c_hexpr_drt_nnodes(drt) &
        bind(C, name="hexpr_drt_nnodes") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_nnodes

    function c_hexpr_drt_ncsfs(drt) &
        bind(C, name="hexpr_drt_ncsfs") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_ncsfs

    function c_hexpr_drt_walk_ncsfs(drt) &
        bind(C, name="hexpr_drt_walk_ncsfs") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_walk_ncsfs

    function c_hexpr_drt_kind(drt) &
        bind(C, name="hexpr_drt_kind") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v         ! hexpr_drt_kind_t cast to int
    end function c_hexpr_drt_kind

    function c_hexpr_drt_ncore(drt) &
        bind(C, name="hexpr_drt_ncore") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_ncore

    function c_hexpr_drt_nval(drt) &
        bind(C, name="hexpr_drt_nval") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_nval

    function c_hexpr_drt_next(drt) &
        bind(C, name="hexpr_drt_next") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_next

    function c_hexpr_drt_nve_state(drt) &
        bind(C, name="hexpr_drt_nve_state") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_nve_state

    function c_hexpr_drt_spin_x2(drt) &
        bind(C, name="hexpr_drt_spin_x2") result(v)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      integer(c_int) :: v
    end function c_hexpr_drt_spin_x2

    ! ---- DRT node info -------------------------------------------

    function c_hexpr_drt_node_info(drt, node_index, &
                                   out_orb, out_nve, out_is) &
        bind(C, name="hexpr_drt_node_info") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: node_index
      integer(c_int), intent(out) :: out_orb, out_nve, out_is
      integer(c_int) :: status    ! hexpr_status_t
    end function c_hexpr_drt_node_info

    ! ---- CSF lookup ----------------------------------------------

    ! C: hexpr_status_t hexpr_csf_index(const hexpr_drt_t *drt,
    !                                   const unsigned char *step,
    !                                   int *out_index)
    ! step is a packed buffer of norb step-code bytes; the public
    ! wrapper passes c_loc(step) for the byte-buffer argument.
    function c_hexpr_csf_index(drt, step, out_index) &
        bind(C, name="hexpr_csf_index") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      type(c_ptr),    value       :: step       ! const unsigned char *
      integer(c_int), intent(out) :: out_index
      integer(c_int)              :: status
    end function c_hexpr_csf_index

    ! C: hexpr_status_t hexpr_csf_steps(const hexpr_drt_t *drt,
    !                                   int index, unsigned char *out_step)
    ! Writes norb step-code bytes into the caller-supplied buffer; the
    ! public wrapper passes c_loc(out_step).
    function c_hexpr_csf_steps(drt, index, out_step) &
        bind(C, name="hexpr_csf_steps") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value :: drt
      integer(c_int), value :: index
      type(c_ptr),    value :: out_step    ! unsigned char *
      integer(c_int)        :: status
    end function c_hexpr_csf_steps

    ! C: hexpr_status_t hexpr_step_parse(const char *str, int norb,
    !                                    unsigned char *out_step)
    ! Parses a NUL-terminated "eudf"-style string of length norb into
    ! norb step-code bytes.
    function c_hexpr_step_parse(str, norb, out_step) &
        bind(C, name="hexpr_step_parse") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value :: str         ! const char *
      integer(c_int), value :: norb
      type(c_ptr),    value :: out_step    ! unsigned char *
      integer(c_int)        :: status
    end function c_hexpr_step_parse

    ! ---- CSF -> determinant expansion ----------------------------

    ! C: hexpr_status_t hexpr_csf_ndets(const hexpr_drt_t *drt,
    !                                   int index, int ms_x2,
    !                                   int *out_ndets)
    ! Exact determinant count, not an upper bound.
    function c_hexpr_csf_ndets(drt, index, ms_x2, out_ndets) &
        bind(C, name="hexpr_csf_ndets") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: index
      integer(c_int), value       :: ms_x2      ! 2*Ms
      integer(c_int), intent(out) :: out_ndets
      integer(c_int)              :: status
    end function c_hexpr_csf_ndets

    ! C: hexpr_status_t hexpr_csf_to_dets(const hexpr_drt_t *drt,
    !                                     int index, int ms_x2,
    !                                     int *out_count,
    !                                     uint64_t *out_alpha,
    !                                     uint64_t *out_beta,
    !                                     double   *out_coef)
    ! The three array arguments are type(c_ptr), value for the same
    ! reason as c_hexpr_expr_read above; the public wrapper passes
    ! c_loc(arr).  Masks are uint64_t on the C side and arrive as
    ! signed integer(c_int64_t) in Fortran -- see hexpr_csf_to_dets.
    function c_hexpr_csf_to_dets(drt, index, ms_x2, out_count, &
                                 out_alpha, out_beta, out_coef) &
        bind(C, name="hexpr_csf_to_dets") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: index
      integer(c_int), value       :: ms_x2      ! 2*Ms
      integer(c_int), intent(out) :: out_count
      type(c_ptr),    value       :: out_alpha  ! uint64_t *
      type(c_ptr),    value       :: out_beta   ! uint64_t *
      type(c_ptr),    value       :: out_coef   ! double *
      integer(c_int)              :: status
    end function c_hexpr_csf_to_dets

    ! ---- DRT bulk ------------------------------------------------

    ! C: const unsigned char *hexpr_drt_csfs(const hexpr_drt_t *)
    ! Returns pointer to ncsfs*norb bytes owned by the DRT (do not free).
    ! The public wrapper re-associates it as a 2-D Fortran pointer.
    function c_hexpr_drt_csfs(drt) &
        bind(C, name="hexpr_drt_csfs") result(ptr)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: drt
      type(c_ptr) :: ptr          ! const unsigned char * (DRT-owned)
    end function c_hexpr_drt_csfs

    ! C: const unsigned char *hexpr_drt_walk_csfs(const hexpr_drt_t *)
    ! Canonical GUGA-walk CSF buffer (walk_ncsfs*norb bytes, DRT-owned).
    function c_hexpr_drt_walk_csfs(drt) &
        bind(C, name="hexpr_drt_walk_csfs") result(ptr)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: drt
      type(c_ptr) :: ptr          ! const unsigned char * (DRT-owned)
    end function c_hexpr_drt_walk_csfs

    ! ---- DRT show ------------------------------------------------

    ! C: void *hexpr_fortran_stdout(void)
    ! Returns the C stdout FILE* for the *_show* wrappers.
    function c_hexpr_fortran_stdout() &
        bind(C, name="hexpr_fortran_stdout") result(p)
      import :: c_ptr
      implicit none
      type(c_ptr) :: p          ! FILE *
    end function c_hexpr_fortran_stdout

    ! C: hexpr_status_t hexpr_drt_show(const hexpr_drt_t *, FILE *out)
    ! The public wrapper passes the C stdout FILE* (c_stdout).
    function c_hexpr_drt_show(drt, out) &
        bind(C, name="hexpr_drt_show") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: drt
      type(c_ptr), value :: out         ! FILE *
      integer(c_int)     :: status
    end function c_hexpr_drt_show

    ! ---- Expr lifecycle ------------------------------------------

    function c_hexpr_expr_build(drt) &
        bind(C, name="hexpr_expr_build") result(handle)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: drt
      type(c_ptr) :: handle       ! hexpr_expr_t *, NULL on failure
    end function c_hexpr_expr_build

    ! C: void hexpr_expr_destroy(hexpr_expr_t *)
    ! The public wrapper nullifies the Fortran handle after this call.
    subroutine c_hexpr_expr_destroy(expr) &
        bind(C, name="hexpr_expr_destroy")
      import :: c_ptr
      implicit none
      type(c_ptr), value :: expr
    end subroutine c_hexpr_expr_destroy

    ! ---- Expr accessors ------------------------------------------

    ! Both return int64_t -- use c_int64_t, not c_int.
    function c_hexpr_expr_nterms_full(expr) &
        bind(C, name="hexpr_expr_nterms_full") result(n)
      import :: c_ptr, c_int64_t
      implicit none
      type(c_ptr), value :: expr
      integer(c_int64_t) :: n
    end function c_hexpr_expr_nterms_full

    function c_hexpr_expr_nterms_one(expr) &
        bind(C, name="hexpr_expr_nterms_one") result(n)
      import :: c_ptr, c_int64_t
      implicit none
      type(c_ptr), value :: expr
      integer(c_int64_t) :: n
    end function c_hexpr_expr_nterms_one

    ! ---- Expr bulk -----------------------------------------------

    ! C: hexpr_status_t hexpr_expr_read(expr, which,
    !                                   int32_t *ij, int32_t *pqrs,
    !                                   double  *coef)
    ! The three array arguments are declared as type(c_ptr), value so
    ! that the Fortran public wrapper can pass c_loc(arr) for assumed-
    ! shape arrays.  On the C side this is identical to passing a
    ! pointer by value (both are 8-byte addresses on LP64 targets).
    ! Assumed-size dummies in a bind(C) interface would not accept
    ! elements of assumed-shape arrays as actual arguments (F95 section12.4.1.7).
    function c_hexpr_expr_read(expr, which, out_ij, out_pqrs, out_coef) &
        bind(C, name="hexpr_expr_read") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value :: expr
      integer(c_int), value :: which     ! HEXPR_EXPR_FULL / _ONE
      type(c_ptr),    value :: out_ij    ! int32_t *
      type(c_ptr),    value :: out_pqrs  ! int32_t *
      type(c_ptr),    value :: out_coef  ! double *
      integer(c_int) :: status           ! hexpr_status_t
    end function c_hexpr_expr_read

    ! ---- Expr show -----------------------------------------------

    ! C: hexpr_status_t hexpr_expr_show_full(const hexpr_expr_t *, FILE *out)
    !    hexpr_status_t hexpr_expr_show_one (const hexpr_expr_t *, FILE *out)
    ! The public wrappers pass the C stdout FILE* (c_stdout).
    function c_hexpr_expr_show_full(expr, out) &
        bind(C, name="hexpr_expr_show_full") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: expr
      type(c_ptr), value :: out         ! FILE *
      integer(c_int)     :: status
    end function c_hexpr_expr_show_full

    function c_hexpr_expr_show_one(expr, out) &
        bind(C, name="hexpr_expr_show_one") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: expr
      type(c_ptr), value :: out         ! FILE *
      integer(c_int)     :: status
    end function c_hexpr_expr_show_one

    ! ---- Per-pair / block / subspace ----------------------------

    function c_hexpr_pair_count(drt, bra, ket, which, out_count) &
        bind(C, name="hexpr_pair_count") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: bra, ket, which
      integer(c_int), intent(out) :: out_count
      integer(c_int)              :: status
    end function c_hexpr_pair_count

    function c_hexpr_pair_max_terms(drt, which, out_max) &
        bind(C, name="hexpr_pair_max_terms") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: which
      integer(c_int), intent(out) :: out_max
      integer(c_int)              :: status
    end function c_hexpr_pair_max_terms

    function c_hexpr_pair_eval(drt, bra, ket, which, out_count, &
                                out_ij, out_pqrs, out_coef) &
        bind(C, name="hexpr_pair_eval") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: bra, ket, which
      integer(c_int), intent(out) :: out_count
      type(c_ptr),    value       :: out_ij, out_pqrs, out_coef
      integer(c_int)              :: status
    end function c_hexpr_pair_eval

    function c_hexpr_block_count(drt, bra_indices, n_bra, &
                                  ket_indices, n_ket, which, out_count) &
        bind(C, name="hexpr_block_count") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      type(c_ptr),    value       :: bra_indices
      integer(c_int), value       :: n_bra
      type(c_ptr),    value       :: ket_indices
      integer(c_int), value       :: n_ket, which
      integer(c_int), intent(out) :: out_count
      integer(c_int)              :: status
    end function c_hexpr_block_count

    function c_hexpr_block_max_terms(drt, n_bra, n_ket, which, out_max) &
        bind(C, name="hexpr_block_max_terms") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: n_bra, n_ket, which
      integer(c_int), intent(out) :: out_max
      integer(c_int)              :: status
    end function c_hexpr_block_max_terms

    function c_hexpr_block_eval(drt, bra_indices, n_bra, ket_indices, n_ket, &
                                 which, out_count, &
                                 out_bra_pos, out_ket_pos, &
                                 out_ij, out_pqrs, out_coef) &
        bind(C, name="hexpr_block_eval") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      type(c_ptr),    value       :: bra_indices
      integer(c_int), value       :: n_bra
      type(c_ptr),    value       :: ket_indices
      integer(c_int), value       :: n_ket, which
      integer(c_int), intent(out) :: out_count
      type(c_ptr),    value       :: out_bra_pos, out_ket_pos
      type(c_ptr),    value       :: out_ij, out_pqrs, out_coef
      integer(c_int)              :: status
    end function c_hexpr_block_eval

    function c_hexpr_subspace_count(drt, indices, n, which, out_count) &
        bind(C, name="hexpr_subspace_count") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      type(c_ptr),    value       :: indices
      integer(c_int), value       :: n, which
      integer(c_int), intent(out) :: out_count
      integer(c_int)              :: status
    end function c_hexpr_subspace_count

    function c_hexpr_subspace_max_terms(drt, n, which, out_max) &
        bind(C, name="hexpr_subspace_max_terms") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: n, which
      integer(c_int), intent(out) :: out_max
      integer(c_int)              :: status
    end function c_hexpr_subspace_max_terms

    function c_hexpr_subspace_eval(drt, indices, n, which, out_count, &
                                    out_bra_pos, out_ket_pos, &
                                    out_ij, out_pqrs, out_coef) &
        bind(C, name="hexpr_subspace_eval") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      type(c_ptr),    value       :: indices
      integer(c_int), value       :: n, which
      integer(c_int), intent(out) :: out_count
      type(c_ptr),    value       :: out_bra_pos, out_ket_pos
      type(c_ptr),    value       :: out_ij, out_pqrs, out_coef
      integer(c_int)              :: status
    end function c_hexpr_subspace_eval

  end interface

  ! ----------------------------------------------------------------
  ! Public procedure declarations (forward refs to contains section)
  ! These names are defined below in contains; the access attribute
  ! is established here in the module specification part.
  ! ----------------------------------------------------------------
  public :: hexpr_init, hexpr_shutdown
  public :: hexpr_version_string, hexpr_version_major, &
            hexpr_version_minor,  hexpr_version_patch
  public :: hexpr_last_error
  public :: hexpr_drt_create,    hexpr_drt_destroy, hexpr_drt_from_csfset
  public :: hexpr_drt_norb,      hexpr_drt_nnodes,    hexpr_drt_ncsfs,  &
            hexpr_drt_kind,      hexpr_drt_ncore,     hexpr_drt_nval,   &
            hexpr_drt_next,      hexpr_drt_nve_state, hexpr_drt_spin_x2
  public :: hexpr_drt_walk_ncsfs
  public :: hexpr_drt_node_info
  public :: hexpr_csf_index, hexpr_csf_steps, hexpr_step_parse
  public :: hexpr_csf_ndets, hexpr_csf_to_dets
  public :: hexpr_drt_csfs
  public :: hexpr_drt_walk_csfs
  public :: hexpr_drt_show
  public :: hexpr_expr_build,       hexpr_expr_destroy
  public :: hexpr_expr_nterms_full, hexpr_expr_nterms_one
  public :: hexpr_expr_read, hexpr_expr_read_full, hexpr_expr_read_one
  public :: hexpr_expr_show_full, hexpr_expr_show_one
  public :: hexpr_pair_count, hexpr_pair_max_terms, hexpr_pair_eval
  public :: hexpr_block_count, hexpr_block_max_terms, hexpr_block_eval
  public :: hexpr_subspace_count, hexpr_subspace_max_terms, hexpr_subspace_eval

contains

  ! ================================================================
  ! Private helper: copy a NUL-terminated C string into a
  ! fixed-length Fortran character variable.
  !
  ! F95-safe: uses c_f_pointer over a rank-1 character(c_char) array,
  ! then copies element by element until c_null_char is found.
  ! No allocatable deferred-length strings; no Fortran 2003 features.
  ! ================================================================
  subroutine c_string_to_fortran(cptr, fstr)
    type(c_ptr),      value,    intent(in)  :: cptr
    character(len=*),           intent(out) :: fstr
    character(kind=c_char), pointer :: buf(:)
    integer :: i, n
    fstr = ''
    if (.not. c_associated(cptr)) return
    n = len(fstr)
    call c_f_pointer(cptr, buf, [n])
    do i = 1, n
      if (buf(i) == c_null_char) exit
      fstr(i:i) = buf(i)
    end do
  end subroutine c_string_to_fortran

  ! ================================================================
  ! (a) Lifecycle
  ! ================================================================

  function hexpr_init() result(status)
    integer(c_int) :: status
    status = c_hexpr_init()
  end function hexpr_init

  subroutine hexpr_shutdown()
    call c_hexpr_shutdown()
  end subroutine hexpr_shutdown

  ! ================================================================
  ! (b) Version / error
  ! ================================================================

  function hexpr_version_string() result(s)
    character(len=64) :: s
    type(c_ptr) :: cptr
    if (.not. version_cached) then
      cptr = c_hexpr_version_string()
      call c_string_to_fortran(cptr, cached_version)
      version_cached = .true.
    end if
    s = cached_version
  end function hexpr_version_string

  function hexpr_version_major() result(v)
    integer(c_int) :: v
    v = c_hexpr_version_major()
  end function hexpr_version_major

  function hexpr_version_minor() result(v)
    integer(c_int) :: v
    v = c_hexpr_version_minor()
  end function hexpr_version_minor

  function hexpr_version_patch() result(v)
    integer(c_int) :: v
    v = c_hexpr_version_patch()
  end function hexpr_version_patch

  ! last_error: NOT cached -- thread-local C pointer is valid only
  ! until the next hexpr call; copy immediately on every call.
  function hexpr_last_error() result(s)
    character(len=256) :: s
    type(c_ptr) :: cptr
    cptr = c_hexpr_last_error()
    call c_string_to_fortran(cptr, s)
  end function hexpr_last_error

  ! ================================================================
  ! (c) DRT lifecycle
  ! ================================================================

  function hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next) &
      result(drt)
    integer(c_int), intent(in) :: kind, nve, spin_x2, ncore, nval, next
    type(c_ptr) :: drt
    drt = c_hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next)
    ! drt is c_null_ptr when C returned NULL (allocation failure)
  end function hexpr_drt_create

  subroutine hexpr_drt_destroy(drt)
    type(c_ptr), intent(inout) :: drt
    if (.not. c_associated(drt)) return
    call c_hexpr_drt_destroy(drt)
    drt = c_null_ptr
  end subroutine hexpr_drt_destroy

  ! csfs must be a contiguous integer(c_int8_t) array of length n_csfs * norb.
  function hexpr_drt_from_csfset(norb, csfs, n_csfs) result(drt)
    integer(c_int),              intent(in) :: norb, n_csfs
    integer(c_int8_t), target,  intent(in) :: csfs(:)
    type(c_ptr)                            :: drt
    drt = c_hexpr_drt_from_csfset(norb, c_loc(csfs), n_csfs)
  end function hexpr_drt_from_csfset

  ! ================================================================
  ! (d) DRT accessors
  ! ================================================================

  function hexpr_drt_norb(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_norb(drt)
  end function hexpr_drt_norb

  function hexpr_drt_nnodes(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_nnodes(drt)
  end function hexpr_drt_nnodes

  function hexpr_drt_ncsfs(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_ncsfs(drt)
  end function hexpr_drt_ncsfs

  ! Canonical GUGA-walk CSF count. CSFSET: may differ from hexpr_drt_ncsfs
  ! (superset M>=N, or dedup); SOCI/FOCI: equals it. See hexpr_drt_walk_csfs.
  function hexpr_drt_walk_ncsfs(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_walk_ncsfs(drt)
  end function hexpr_drt_walk_ncsfs

  function hexpr_drt_kind(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_kind(drt)
  end function hexpr_drt_kind

  function hexpr_drt_ncore(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_ncore(drt)
  end function hexpr_drt_ncore

  function hexpr_drt_nval(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_nval(drt)
  end function hexpr_drt_nval

  function hexpr_drt_next(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_next(drt)
  end function hexpr_drt_next

  function hexpr_drt_nve_state(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_nve_state(drt)
  end function hexpr_drt_nve_state

  function hexpr_drt_spin_x2(drt) result(v)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: v
    v = c_hexpr_drt_spin_x2(drt)
  end function hexpr_drt_spin_x2

  ! ================================================================
  ! (e) DRT node info
  ! ================================================================

  function hexpr_drt_node_info(drt, node_index, out_orb, out_nve, out_is) &
      result(status)
    type(c_ptr),    value,    intent(in)  :: drt
    integer(c_int),           intent(in)  :: node_index
    integer(c_int),           intent(out) :: out_orb, out_nve, out_is
    integer(c_int) :: status
    status = c_hexpr_drt_node_info(drt, node_index, out_orb, out_nve, out_is)
  end function hexpr_drt_node_info

  ! ================================================================
  ! (e2) CSF lookup
  ! ================================================================

  ! step must be a contiguous integer(c_int8_t) array of length norb.
  ! out_index receives the 0-based CSF row index.
  function hexpr_csf_index(drt, step, out_index) result(status)
    type(c_ptr),       value,    intent(in)  :: drt
    integer(c_int8_t), target,   intent(in)  :: step(:)
    integer(c_int),              intent(out) :: out_index
    integer(c_int) :: status
    status = c_hexpr_csf_index(drt, c_loc(step), out_index)
  end function hexpr_csf_index

  ! out_step must be a contiguous integer(c_int8_t) buffer of length
  ! at least hexpr_drt_norb(drt); the C side writes norb step codes.
  function hexpr_csf_steps(drt, index, out_step) result(status)
    type(c_ptr),       value,    intent(in)  :: drt
    integer(c_int),              intent(in)  :: index
    integer(c_int8_t), target,   intent(out) :: out_step(:)
    integer(c_int) :: status
    status = c_hexpr_csf_steps(drt, index, c_loc(out_step))
  end function hexpr_csf_steps

  ! Parse an "eudf"-style string into out_step (length >= norb).
  ! The string is copied into a NUL-terminated, contiguous c_char buffer
  ! (c_loc requires the target to be contiguous) before the C call.
  function hexpr_step_parse(str, norb, out_step) result(status)
    character(len=*),            intent(in)  :: str
    integer(c_int),              intent(in)  :: norb
    integer(c_int8_t), target,   intent(out) :: out_step(:)
    integer(c_int) :: status
    character(kind=c_char), target :: cbuf(len(str) + 1)
    integer :: i, n
    n = len(str)
    do i = 1, n
      cbuf(i) = str(i:i)
    end do
    cbuf(n + 1) = c_null_char
    status = c_hexpr_step_parse(c_loc(cbuf), norb, c_loc(out_step))
  end function hexpr_step_parse

  ! ================================================================
  ! (e3) CSF -> determinant expansion
  !
  ! Array shape: the C buffers are out_alpha[k*nwords + w] (row-major
  ! in C, word index fastest).  Fortran's first subscript is fastest,
  ! so the same bytes are alpha(w, k) with shape (nwords, ndets) --
  ! the identical reinterpretation already used by hexpr_drt_csfs.
  ! No copy, no transpose.
  !
  ! Unsigned masks: iso_c_binding has no c_uint64_t.  The C side's
  ! uint64_t arrives as SIGNED integer(c_int64_t), so a determinant
  ! occupying spatial orbital 63 sets the sign bit and the value reads
  ! as negative.  Test occupation with BTEST, which works on the bit
  ! pattern; never with a magnitude comparison, a decimal write, or a
  ! use as an array subscript.  Spatial orbital p (0-based) is bit
  ! IAND(p, 63) of word ISHFT(p, -6) + 1.
  !
  ! Conventions of the expansion itself (coupling order, CG phase,
  ! operator order, determinant order, normalization, and the fact
  ! that the OVERALL sign of a CSF is a convention rather than a
  ! measured quantity) are stated in hexpr.h and are not restated per
  ! wrapper; see c/src/csf/csf2det.c for the derivation.
  ! ================================================================

  ! index follows the same convention as hexpr_csf_steps: the input
  ! position, in [0, hexpr_drt_ncsfs).  So do hexpr_pair_eval and
  ! friends -- there is one index space across the whole API.
  ! ms_x2 is 2*Ms; any |Ms| <= S of the same parity as spin_x2 is
  ! accepted, negative values included.  out_ndets is exact.
  function hexpr_csf_ndets(drt, index, ms_x2, out_ndets) result(status)
    type(c_ptr),       value,    intent(in)  :: drt
    integer(c_int),              intent(in)  :: index
    integer(c_int),              intent(in)  :: ms_x2
    integer(c_int),              intent(out) :: out_ndets
    integer(c_int) :: status
    status = c_hexpr_csf_ndets(drt, index, ms_x2, out_ndets)
  end function hexpr_csf_ndets

  ! Size contract: alpha and beta must each be at least
  ! (nwords, ndets) and coef at least (ndets), where ndets comes from
  ! hexpr_csf_ndets and nwords = (hexpr_drt_norb(drt) + 63) / 64.
  ! A short buffer is a caller error and is not diagnosed here, as
  ! with hexpr_expr_read and hexpr_csf_steps.
  ! Callers must pass contiguous arrays (allocatables always satisfy
  ! this; pointer arrays may need CONTIGUOUS).
  function hexpr_csf_to_dets(drt, index, ms_x2, out_count, &
                             alpha, beta, coef) result(status)
    type(c_ptr),        value,       intent(in)  :: drt
    integer(c_int),                  intent(in)  :: index
    integer(c_int),                  intent(in)  :: ms_x2
    integer(c_int),                  intent(out) :: out_count
    integer(c_int64_t), target,      intent(out) :: alpha(:,:)
    integer(c_int64_t), target,      intent(out) :: beta(:,:)
    real(c_double),     target,      intent(out) :: coef(:)
    integer(c_int) :: status
    status = c_hexpr_csf_to_dets(drt, index, ms_x2, out_count, &
                                 c_loc(alpha), c_loc(beta), c_loc(coef))
  end function hexpr_csf_to_dets

  ! ================================================================
  ! (f) DRT bulk: csfs as 2-D pointer
  !
  ! C buffer layout: csfs[csf_i * norb + orb_j]  (row-major in C).
  ! orb_j (0..norb-1) is the fastest-changing index in C.
  !
  ! Fortran column-major convention: the FIRST subscript is fastest.
  ! We therefore declare shape (norb, ncsfs): buf(orb, csf).
  ! buf(orb, csf) corresponds to C csfs[(csf-1)*norb + (orb-1)].
  ! In memory the two layouts are byte-for-byte identical -- no
  ! data copy, no transpose, only index interpretation differs.
  ! ================================================================
  function hexpr_drt_csfs(drt, buf) result(status)
    type(c_ptr),                     value,    intent(in)  :: drt
    integer(c_int8_t), pointer,                intent(out) :: buf(:,:)
    integer(c_int) :: status
    type(c_ptr)    :: cptr
    integer(c_int) :: norb, ncsfs_count

    norb        = c_hexpr_drt_norb(drt)
    ncsfs_count = c_hexpr_drt_ncsfs(drt)
    cptr        = c_hexpr_drt_csfs(drt)

    if (.not. c_associated(cptr)) then
      buf    => null()
      status = HEXPR_ERR_NOT_BUILT
      return
    end if

    call c_f_pointer(cptr, buf, [norb, ncsfs_count])
    status = HEXPR_OK
  end function hexpr_drt_csfs

  ! ================================================================
  ! Canonical GUGA-walk CSF list (CSFSET: superset/canonical view).
  ! Same column-major (norb, walk_ncsfs) convention as hexpr_drt_csfs;
  ! re-associates the DRT-owned buffer as a 2-D Fortran pointer.
  ! ================================================================
  function hexpr_drt_walk_csfs(drt, buf) result(status)
    type(c_ptr),                     value,    intent(in)  :: drt
    integer(c_int8_t), pointer,                intent(out) :: buf(:,:)
    integer(c_int) :: status
    type(c_ptr)    :: cptr
    integer(c_int) :: norb, ncsfs_count

    norb        = c_hexpr_drt_norb(drt)
    ncsfs_count = c_hexpr_drt_walk_ncsfs(drt)
    cptr        = c_hexpr_drt_walk_csfs(drt)

    if (.not. c_associated(cptr)) then
      buf    => null()
      status = HEXPR_ERR_NOT_BUILT
      return
    end if

    call c_f_pointer(cptr, buf, [norb, ncsfs_count])
    status = HEXPR_OK
  end function hexpr_drt_walk_csfs

  ! ================================================================
  ! (f2) DRT show: dump to stdout (no stream argument)
  ! ================================================================

  function hexpr_drt_show(drt) result(status)
    type(c_ptr), value, intent(in) :: drt
    integer(c_int) :: status
    status = c_hexpr_drt_show(drt, c_hexpr_fortran_stdout())
  end function hexpr_drt_show

  ! ================================================================
  ! (g) Expr lifecycle
  ! ================================================================

  function hexpr_expr_build(drt) result(expr)
    type(c_ptr), value, intent(in) :: drt
    type(c_ptr) :: expr
    expr = c_hexpr_expr_build(drt)
    ! expr is c_null_ptr when C returned NULL (build failure)
  end function hexpr_expr_build

  subroutine hexpr_expr_destroy(expr)
    type(c_ptr), intent(inout) :: expr
    if (.not. c_associated(expr)) return
    call c_hexpr_expr_destroy(expr)
    expr = c_null_ptr
  end subroutine hexpr_expr_destroy

  ! ================================================================
  ! (h) Expr accessors
  ! ================================================================

  function hexpr_expr_nterms_full(expr) result(n)
    type(c_ptr), value, intent(in) :: expr
    integer(c_int64_t) :: n
    n = c_hexpr_expr_nterms_full(expr)
  end function hexpr_expr_nterms_full

  function hexpr_expr_nterms_one(expr) result(n)
    type(c_ptr), value, intent(in) :: expr
    integer(c_int64_t) :: n
    n = c_hexpr_expr_nterms_one(expr)
  end function hexpr_expr_nterms_one

  ! ================================================================
  ! (i) Expr bulk
  !
  ! Raw: hexpr_expr_read takes the which selector explicitly.
  ! Sugar: hexpr_expr_read_full / hexpr_expr_read_one fix which.
  !
  ! Array passing: the c_hexpr_expr_read bind(C) interface takes the
  ! array arguments as type(c_ptr) by value (i.e. raw C pointers).
  ! We obtain those addresses with c_loc() on the assumed-shape
  ! dummies, which requires the dummy arrays to have the TARGET
  ! attribute.  Callers must pass contiguous arrays (allocatables
  ! always satisfy this; pointer arrays may need CONTIGUOUS).
  ! ================================================================

  function hexpr_expr_read(expr, which, ij, pqrs, coef) result(status)
    type(c_ptr),        value,          intent(in)  :: expr
    integer(c_int),                     intent(in)  :: which
    integer(c_int32_t), target,         intent(out) :: ij(:)
    integer(c_int32_t), target,         intent(out) :: pqrs(:)
    real(c_double),     target,         intent(out) :: coef(:)
    integer(c_int) :: status
    status = c_hexpr_expr_read(expr, which, &
                               c_loc(ij), c_loc(pqrs), c_loc(coef))
  end function hexpr_expr_read

  function hexpr_expr_read_full(expr, ij, pqrs, coef) result(status)
    type(c_ptr),        value,          intent(in)  :: expr
    integer(c_int32_t), target,         intent(out) :: ij(:)
    integer(c_int32_t), target,         intent(out) :: pqrs(:)
    real(c_double),     target,         intent(out) :: coef(:)
    integer(c_int) :: status
    status = hexpr_expr_read(expr, HEXPR_EXPR_FULL, ij, pqrs, coef)
  end function hexpr_expr_read_full

  function hexpr_expr_read_one(expr, ij, pqrs, coef) result(status)
    type(c_ptr),        value,          intent(in)  :: expr
    integer(c_int32_t), target,         intent(out) :: ij(:)
    integer(c_int32_t), target,         intent(out) :: pqrs(:)
    real(c_double),     target,         intent(out) :: coef(:)
    integer(c_int) :: status
    status = hexpr_expr_read(expr, HEXPR_EXPR_ONE, ij, pqrs, coef)
  end function hexpr_expr_read_one

  ! ================================================================
  ! (i2) Expr show: dump to stdout (no stream argument)
  ! ================================================================

  function hexpr_expr_show_full(expr) result(status)
    type(c_ptr), value, intent(in) :: expr
    integer(c_int) :: status
    status = c_hexpr_expr_show_full(expr, c_hexpr_fortran_stdout())
  end function hexpr_expr_show_full

  function hexpr_expr_show_one(expr) result(status)
    type(c_ptr), value, intent(in) :: expr
    integer(c_int) :: status
    status = c_hexpr_expr_show_one(expr, c_hexpr_fortran_stdout())
  end function hexpr_expr_show_one

  ! ================================================================
  ! (j) Per-pair / block / subspace evaluation
  !
  ! Every CSF index below is an INPUT-ORDER index in
  ! [0, hexpr_drt_ncsfs) -- the same space as hexpr_drt_csfs,
  ! hexpr_csf_index and hexpr_csf_steps.  There is one index space
  ! across the whole API.  For SOCI and FOCI it is also the canonical
  ! walk order; for a CSFSET sub-DRT the canonical walk may reorder,
  ! extend or dedup the input, and the library translates internally.
  !
  ! out_ij is the one exception and is not an identifier a caller
  ! needs here.  Per docs/PQRS_INDEX_SPEC.md (normative for v1.0.0)
  ! ij is triangle() of the two Shavitt PATH WEIGHTS, so it stays in
  ! the canonical space; an input position is not a path weight.  To
  ! attribute a term to a pair use out_bra_pos / out_ket_pos from the
  ! block and subspace forms: they are positions within the arrays
  ! you passed.  Do NOT decode out_ij and feed the result back to
  ! these calls or to hexpr_csf_steps.
  ! ================================================================

  function hexpr_pair_count(drt, bra, ket, which, out_count) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int),         intent(in)  :: bra, ket, which
    integer(c_int),         intent(out) :: out_count
    integer(c_int) :: status
    status = c_hexpr_pair_count(drt, bra, ket, which, out_count)
  end function hexpr_pair_count

  function hexpr_pair_max_terms(drt, which, out_max) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int),         intent(in)  :: which
    integer(c_int),         intent(out) :: out_max
    integer(c_int) :: status
    status = c_hexpr_pair_max_terms(drt, which, out_max)
  end function hexpr_pair_max_terms

  ! out_ij, out_pqrs, out_coef must be allocated to at least pair_max_terms
  ! elements by the caller. out_count is set to the number of terms written.
  function hexpr_pair_eval(drt, bra, ket, which, out_count, &
                            out_ij, out_pqrs, out_coef) result(status)
    type(c_ptr),        value,  intent(in)  :: drt
    integer(c_int),             intent(in)  :: bra, ket, which
    integer(c_int),             intent(out) :: out_count
    integer(c_int32_t), target, intent(out) :: out_ij(:), out_pqrs(:)
    real(c_double),     target, intent(out) :: out_coef(:)
    integer(c_int) :: status
    status = c_hexpr_pair_eval(drt, bra, ket, which, out_count, &
                                c_loc(out_ij), c_loc(out_pqrs), c_loc(out_coef))
  end function hexpr_pair_eval

  function hexpr_block_count(drt, bra_indices, n_bra, ket_indices, n_ket, &
                              which, out_count) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int), target, intent(in)  :: bra_indices(:)
    integer(c_int),         intent(in)  :: n_bra
    integer(c_int), target, intent(in)  :: ket_indices(:)
    integer(c_int),         intent(in)  :: n_ket, which
    integer(c_int),         intent(out) :: out_count
    integer(c_int) :: status
    status = c_hexpr_block_count(drt, c_loc(bra_indices), n_bra, &
                                  c_loc(ket_indices), n_ket, which, out_count)
  end function hexpr_block_count

  function hexpr_block_max_terms(drt, n_bra, n_ket, which, out_max) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int),         intent(in)  :: n_bra, n_ket, which
    integer(c_int),         intent(out) :: out_max
    integer(c_int) :: status
    status = c_hexpr_block_max_terms(drt, n_bra, n_ket, which, out_max)
  end function hexpr_block_max_terms

  ! Output arrays must be allocated to at least block_max_terms elements.
  function hexpr_block_eval(drt, bra_indices, n_bra, ket_indices, n_ket, &
                             which, out_count, &
                             out_bra_pos, out_ket_pos, &
                             out_ij, out_pqrs, out_coef) result(status)
    type(c_ptr),        value,  intent(in)  :: drt
    integer(c_int),     target, intent(in)  :: bra_indices(:)
    integer(c_int),             intent(in)  :: n_bra
    integer(c_int),     target, intent(in)  :: ket_indices(:)
    integer(c_int),             intent(in)  :: n_ket, which
    integer(c_int),             intent(out) :: out_count
    integer(c_int32_t), target, intent(out) :: out_bra_pos(:), out_ket_pos(:)
    integer(c_int32_t), target, intent(out) :: out_ij(:), out_pqrs(:)
    real(c_double),     target, intent(out) :: out_coef(:)
    integer(c_int) :: status
    status = c_hexpr_block_eval(drt, c_loc(bra_indices), n_bra, &
                                 c_loc(ket_indices), n_ket, which, out_count, &
                                 c_loc(out_bra_pos), c_loc(out_ket_pos), &
                                 c_loc(out_ij), c_loc(out_pqrs), c_loc(out_coef))
  end function hexpr_block_eval

  function hexpr_subspace_count(drt, indices, n, which, out_count) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int), target, intent(in)  :: indices(:)
    integer(c_int),         intent(in)  :: n, which
    integer(c_int),         intent(out) :: out_count
    integer(c_int) :: status
    status = c_hexpr_subspace_count(drt, c_loc(indices), n, which, out_count)
  end function hexpr_subspace_count

  function hexpr_subspace_max_terms(drt, n, which, out_max) result(status)
    type(c_ptr),    value,  intent(in)  :: drt
    integer(c_int),         intent(in)  :: n, which
    integer(c_int),         intent(out) :: out_max
    integer(c_int) :: status
    status = c_hexpr_subspace_max_terms(drt, n, which, out_max)
  end function hexpr_subspace_max_terms

  ! Output arrays must be allocated to at least subspace_max_terms elements.
  function hexpr_subspace_eval(drt, indices, n, which, out_count, &
                                out_bra_pos, out_ket_pos, &
                                out_ij, out_pqrs, out_coef) result(status)
    type(c_ptr),        value,  intent(in)  :: drt
    integer(c_int),     target, intent(in)  :: indices(:)
    integer(c_int),             intent(in)  :: n, which
    integer(c_int),             intent(out) :: out_count
    integer(c_int32_t), target, intent(out) :: out_bra_pos(:), out_ket_pos(:)
    integer(c_int32_t), target, intent(out) :: out_ij(:), out_pqrs(:)
    real(c_double),     target, intent(out) :: out_coef(:)
    integer(c_int) :: status
    status = c_hexpr_subspace_eval(drt, c_loc(indices), n, which, out_count, &
                                    c_loc(out_bra_pos), c_loc(out_ket_pos), &
                                    c_loc(out_ij), c_loc(out_pqrs), c_loc(out_coef))
  end function hexpr_subspace_eval

end module hexpr_mod
