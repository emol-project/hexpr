! Test expression build, term counts, and bulk read via the hexpr F90 module.
! Reference case: s_full_4o4e (SOCI, nve=4, spin=0, ncore=0, nval=4, next=0)
!   nterms_full=638, nterms_one=124
program test_expr
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr)                          :: drt, expr
  integer(c_int)                       :: st
  integer(c_int64_t)                   :: nf, no
  integer(c_int32_t), allocatable, target :: ij_arr(:), pqrs_arr(:)
  real(c_double),     allocatable, target :: coef_arr(:)

  st = hexpr_init()
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_init returned', st
    stop 1
  end if

  drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4_c_int, 0_c_int, &
                         0_c_int, 4_c_int, 0_c_int)
  if (.not. c_associated(drt)) then
    write(*,*) 'FAIL: hexpr_drt_create returned NULL'
    stop 1
  end if

  expr = hexpr_expr_build(drt)
  if (.not. c_associated(expr)) then
    write(*,*) 'FAIL: hexpr_expr_build returned NULL'
    stop 1
  end if

  nf = hexpr_expr_nterms_full(expr)
  no = hexpr_expr_nterms_one(expr)

  if (nf /= 638_c_int64_t) then
    write(*,'(a,i0,a)') 'FAIL: nterms_full = ', nf, ' (expected 638)'
    stop 1
  end if
  if (no /= 124_c_int64_t) then
    write(*,'(a,i0,a)') 'FAIL: nterms_one = ', no, ' (expected 124)'
    stop 1
  end if

  ! Read full expression into allocatable arrays
  allocate(ij_arr(nf), pqrs_arr(nf), coef_arr(nf))
  st = hexpr_expr_read_full(expr, ij_arr, pqrs_arr, coef_arr)
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_expr_read_full returned', st
    stop 1
  end if
  deallocate(ij_arr, pqrs_arr, coef_arr)

  ! Read one-electron expression
  allocate(ij_arr(no), pqrs_arr(no), coef_arr(no))
  st = hexpr_expr_read_one(expr, ij_arr, pqrs_arr, coef_arr)
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_expr_read_one returned', st
    stop 1
  end if
  deallocate(ij_arr, pqrs_arr, coef_arr)

  call hexpr_expr_destroy(expr)
  call hexpr_drt_destroy(drt)
  call hexpr_shutdown()

  write(*,*) 'test_expr                                          PASS'
  stop 0
end program test_expr
