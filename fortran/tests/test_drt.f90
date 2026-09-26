! Test DRT construction and accessors via the hexpr F90 module.
! Reference case: s_full_4o4e (SOCI, nve=4, spin=0, ncore=0, nval=4, next=0)
!   norb=4, nnodes=14, ncsfs=20
!   first CSF: [E,E,F,F] = step codes [0,0,3,3]
program test_drt
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr)                     :: drt
  integer(c_int)                  :: st
  integer(c_int8_t), pointer      :: csfs(:,:)

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

  if (hexpr_drt_norb(drt) /= 4) then
    write(*,'(a,i0,a)') 'FAIL: norb = ', hexpr_drt_norb(drt), ' (expected 4)'
    stop 1
  end if
  if (hexpr_drt_nnodes(drt) /= 14) then
    write(*,'(a,i0,a)') 'FAIL: nnodes = ', hexpr_drt_nnodes(drt), ' (expected 14)'
    stop 1
  end if
  if (hexpr_drt_ncsfs(drt) /= 20) then
    write(*,'(a,i0,a)') 'FAIL: ncsfs = ', hexpr_drt_ncsfs(drt), ' (expected 20)'
    stop 1
  end if
  if (hexpr_drt_kind(drt) /= HEXPR_DRT_SOCI) then
    write(*,*) 'FAIL: kind mismatch'
    stop 1
  end if
  if (hexpr_drt_ncore(drt) /= 0) then
    write(*,'(a,i0)') 'FAIL: ncore = ', hexpr_drt_ncore(drt)
    stop 1
  end if
  if (hexpr_drt_nval(drt) /= 4) then
    write(*,'(a,i0)') 'FAIL: nval = ', hexpr_drt_nval(drt)
    stop 1
  end if
  if (hexpr_drt_next(drt) /= 0) then
    write(*,'(a,i0)') 'FAIL: next = ', hexpr_drt_next(drt)
    stop 1
  end if
  if (hexpr_drt_nve_state(drt) /= 4) then
    write(*,'(a,i0)') 'FAIL: nve_state = ', hexpr_drt_nve_state(drt)
    stop 1
  end if
  if (hexpr_drt_spin_x2(drt) /= 0) then
    write(*,'(a,i0)') 'FAIL: spin_x2 = ', hexpr_drt_spin_x2(drt)
    stop 1
  end if

  ! CSF buffer: shape (norb, ncsfs), Fortran order (orb fastest)
  ! First CSF is [E,E,F,F] = [0,0,3,3] in orb order 1..4
  st = hexpr_drt_csfs(drt, csfs)
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_drt_csfs returned', st
    stop 1
  end if
  if (csfs(1,1) /= int(HEXPR_STEP_E, c_int8_t) .or. &
      csfs(2,1) /= int(HEXPR_STEP_E, c_int8_t) .or. &
      csfs(3,1) /= int(HEXPR_STEP_F, c_int8_t) .or. &
      csfs(4,1) /= int(HEXPR_STEP_F, c_int8_t)) then
    write(*,'(a,4i3)') 'FAIL: first CSF step codes = ', csfs(1:4,1)
    stop 1
  end if

  call hexpr_drt_destroy(drt)
  if (c_associated(drt)) then
    write(*,*) 'FAIL: handle not nullified after destroy'
    stop 1
  end if

  call hexpr_shutdown()
  write(*,*) 'test_drt                                           PASS'
  stop 0
end program test_drt
