program test_smoke
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int) :: st, maj, min, pat
  character(len=64) :: vs

  st = hexpr_init()
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_init returned', st
    stop 1
  end if

  vs  = hexpr_version_string()
  maj = hexpr_version_major()
  min = hexpr_version_minor()
  pat = hexpr_version_patch()

  if (maj /= HEXPR_VMAJOR) then
    write(*,'(a,i0,a,i0)') 'FAIL: version_major = ', maj, ', expected ', HEXPR_VMAJOR
    stop 1
  end if
  if (min /= HEXPR_VMINOR) then
    write(*,'(a,i0,a,i0)') 'FAIL: version_minor = ', min, ', expected ', HEXPR_VMINOR
    stop 1
  end if
  if (pat /= HEXPR_VPATCH) then
    write(*,'(a,i0,a,i0)') 'FAIL: version_patch = ', pat, ', expected ', HEXPR_VPATCH
    stop 1
  end if

  write(*,'(2a)')    '  version_string : ', trim(vs)
  write(*,'(a,i0,a,i0,a,i0)') '  version M.m.p  : ', maj, '.', min, '.', pat

  call hexpr_shutdown()

  write(*,*) 'test_smoke                                         PASS'
  stop 0
end program test_smoke
