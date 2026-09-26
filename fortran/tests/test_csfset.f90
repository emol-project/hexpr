! Test hexpr_drt_from_csfset via the hexpr F90 module.
! Reference case: s_full_4o4e (SOCI, nve=4, spin=0, ncore=0, nval=4, next=0)
!   norb=4, nnodes=14, ncsfs=20
!
! Checks:
!   T1: round-trip -- sub-DRT from full CSF set has same norb/ncsfs/nnodes
!   T2: 2-CSF sub-DRT has fewer nnodes and ncsfs==2
!   T3: hexpr_pair_eval(0,0) on 2-CSF sub-DRT returns HEXPR_OK with
!       a sensible term count (0 <= n_terms <= max_terms upper bound)
program test_csfset
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  ! Local bind(C) interfaces for pair evaluation (not yet in hexpr_mod)
  interface
    function c_pair_max_terms(drt, which, out_max) &
        bind(C, name="hexpr_pair_max_terms") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: which
      integer(c_int), intent(out) :: out_max
      integer(c_int)              :: status
    end function c_pair_max_terms

    function c_pair_eval(drt, bra, ket, which, &
                         out_count, out_ij, out_pqrs, out_coef) &
        bind(C, name="hexpr_pair_eval") result(status)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr),    value       :: drt
      integer(c_int), value       :: bra, ket, which
      integer(c_int), intent(out) :: out_count
      type(c_ptr),    value       :: out_ij, out_pqrs, out_coef
      integer(c_int)              :: status
    end function c_pair_eval
  end interface

  type(c_ptr)                              :: soci_drt, sub_full, sub_2, sub_super
  integer(c_int)                           :: st
  integer(c_int8_t), pointer               :: csf_ptr(:,:)
  integer(c_int8_t), pointer               :: walk_ptr(:,:), csfs_chk(:,:)
  integer(c_int8_t), allocatable, target   :: csf_copy(:), csf_2(:), csf_super(:)
  integer(c_int)                           :: norb, ncsfs, nnodes, walk_m, n_super
  integer(c_int)                           :: max_terms, n_terms
  integer(c_int32_t), allocatable, target  :: out_ij(:), out_pqrs(:)
  real(c_double),     allocatable, target  :: out_coef(:)
  integer :: nfail

  nfail = 0

  st = hexpr_init()
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') 'FAIL: hexpr_init returned ', st
    stop 1
  end if

  ! Build the reference s_full_4o4e SOCI DRT
  soci_drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4_c_int, 0_c_int, &
                               0_c_int, 4_c_int, 0_c_int)
  if (.not. c_associated(soci_drt)) then
    write(*,*) 'FAIL: hexpr_drt_create returned NULL'
    stop 1
  end if

  norb   = hexpr_drt_norb(soci_drt)
  ncsfs  = hexpr_drt_ncsfs(soci_drt)
  nnodes = hexpr_drt_nnodes(soci_drt)

  ! Get CSF buffer and copy to a local contiguous array
  st = hexpr_drt_csfs(soci_drt, csf_ptr)
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') 'FAIL: hexpr_drt_csfs returned ', st
    stop 1
  end if
  allocate(csf_copy(norb * ncsfs))
  csf_copy = reshape(csf_ptr, [norb * ncsfs])

  ! ----------------------------------------------------------------
  ! T1: round-trip on full CSF set
  ! ----------------------------------------------------------------
  sub_full = hexpr_drt_from_csfset(norb, csf_copy, ncsfs)
  if (.not. c_associated(sub_full)) then
    write(*,*) '  T1 FAIL: sub_full NULL:', trim(hexpr_last_error())
    nfail = nfail + 1
  else
    if (hexpr_drt_norb(sub_full) /= norb) then
      write(*,'(a,i0,a,i0)') '  T1 FAIL: norb=', hexpr_drt_norb(sub_full), &
                              ' expected', norb
      nfail = nfail + 1
    else if (hexpr_drt_ncsfs(sub_full) /= ncsfs) then
      write(*,'(a,i0,a,i0)') '  T1 FAIL: ncsfs=', hexpr_drt_ncsfs(sub_full), &
                              ' expected', ncsfs
      nfail = nfail + 1
    else if (hexpr_drt_nnodes(sub_full) /= nnodes) then
      write(*,'(a,i0,a,i0)') '  T1 FAIL: nnodes=', hexpr_drt_nnodes(sub_full), &
                              ' expected', nnodes
      nfail = nfail + 1
    else
      write(*,*) '  T1 round-trip (full CSF set): PASS'
    end if
    call hexpr_drt_destroy(sub_full)
  end if

  ! ----------------------------------------------------------------
  ! T2: 2-CSF sub-DRT
  ! ----------------------------------------------------------------
  allocate(csf_2(2 * norb))
  csf_2 = csf_copy(1 : 2 * norb)
  sub_2 = hexpr_drt_from_csfset(norb, csf_2, 2_c_int)

  if (.not. c_associated(sub_2)) then
    write(*,*) '  T2 FAIL: sub_2 NULL:', trim(hexpr_last_error())
    nfail = nfail + 1
    ! sub_2 unavailable -- skip T3
  else
    if (hexpr_drt_nnodes(sub_2) >= nnodes) then
      write(*,'(a,i0,a,i0)') '  T2 FAIL: sub_2 nnodes=', &
                              hexpr_drt_nnodes(sub_2), &
                              ' not < soci nnodes=', nnodes
      nfail = nfail + 1
    end if
    if (hexpr_drt_ncsfs(sub_2) /= 2) then
      write(*,'(a,i0)') '  T2 FAIL: sub_2 ncsfs=', hexpr_drt_ncsfs(sub_2)
      nfail = nfail + 1
    end if
    if (nfail == 0) write(*,*) '  T2 2-CSF sub-DRT shape: PASS'

    ! ----------------------------------------------------------------
    ! T3: pair_eval(0,0) FULL on the 2-CSF sub-DRT
    ! ----------------------------------------------------------------
    st = c_pair_max_terms(sub_2, HEXPR_EXPR_FULL, max_terms)
    if (st /= HEXPR_OK) then
      write(*,'(a,i0)') '  T3 FAIL: pair_max_terms returned ', st
      nfail = nfail + 1
    else
      allocate(out_ij(max_terms), out_pqrs(max_terms), out_coef(max_terms))
      st = c_pair_eval(sub_2, 0_c_int, 0_c_int, HEXPR_EXPR_FULL, &
                       n_terms, c_loc(out_ij), c_loc(out_pqrs), c_loc(out_coef))
      if (st /= HEXPR_OK) then
        write(*,'(a,i0)') '  T3 FAIL: pair_eval returned ', st
        nfail = nfail + 1
      else if (n_terms < 0 .or. n_terms > max_terms) then
        write(*,'(a,i0,a,i0)') '  T3 FAIL: n_terms=', n_terms, &
                                ' max=', max_terms
        nfail = nfail + 1
      else
        write(*,'(a,i0,a,i0,a)') '  T3 pair_eval(0,0) FULL: n_terms=', &
                                  n_terms, ' (max=', max_terms, '): PASS'
      end if
      deallocate(out_ij, out_pqrs, out_coef)
    end if

    call hexpr_drt_destroy(sub_2)
  end if

  ! ----------------------------------------------------------------
  ! T4: input-order csfs() + walk_* superset on a non-canonical input.
  !     Input [uefd, eudf] (N=2) builds a canonical superset (M=4).
  ! ----------------------------------------------------------------
  allocate(csf_super(2 * 4))
  csf_super = int([1,0,3,2,  0,1,2,3], c_int8_t)   ! row0=uefd, row1=eudf
  sub_super = hexpr_drt_from_csfset(4_c_int, csf_super, 2_c_int)
  if (.not. c_associated(sub_super)) then
    write(*,*) '  T4 FAIL: sub_super NULL:', trim(hexpr_last_error())
    nfail = nfail + 1
  else
    ! csfs()/ncsfs == input set in input order
    if (hexpr_drt_ncsfs(sub_super) /= 2) then
      write(*,'(a,i0)') '  T4 FAIL: ncsfs=', hexpr_drt_ncsfs(sub_super)
      nfail = nfail + 1
    end if
    st = hexpr_drt_csfs(sub_super, csfs_chk)   ! shape (norb, ncsfs)
    if (st /= HEXPR_OK) then
      write(*,'(a,i0)') '  T4 FAIL: csfs returned ', st
      nfail = nfail + 1
    else if (any(csfs_chk(:,1) /= int([1,0,3,2], c_int8_t)) .or. &
             any(csfs_chk(:,2) /= int([0,1,2,3], c_int8_t))) then
      write(*,*) '  T4 FAIL: csfs() not in input order'
      nfail = nfail + 1
    end if
    ! walk_* == canonical superset M >= N
    walk_m = hexpr_drt_walk_ncsfs(sub_super)
    n_super = hexpr_drt_ncsfs(sub_super)
    if ((walk_m /= 4) .or. (walk_m < n_super)) then
      write(*,'(a,i0)') '  T4 FAIL: walk_ncsfs=', walk_m
      nfail = nfail + 1
    end if
    st = hexpr_drt_walk_csfs(sub_super, walk_ptr)  ! shape (norb, walk_ncsfs)
    if (st /= HEXPR_OK) then
      write(*,'(a,i0)') '  T4 FAIL: walk_csfs returned ', st
      nfail = nfail + 1
    else if (size(walk_ptr, 2) /= 4) then
      write(*,'(a,i0)') '  T4 FAIL: walk_csfs cols=', size(walk_ptr, 2)
      nfail = nfail + 1
    end if
    ! SOCI sanity: walk_* equals csfs/ncsfs (no input set)
    if (hexpr_drt_walk_ncsfs(soci_drt) /= hexpr_drt_ncsfs(soci_drt)) then
      write(*,*) '  T4 FAIL: SOCI walk_ncsfs /= ncsfs'
      nfail = nfail + 1
    end if
    if (nfail == 0) write(*,*) '  T4 input-order csfs + walk_* superset: PASS'
    call hexpr_drt_destroy(sub_super)
  end if

  call hexpr_drt_destroy(soci_drt)
  deallocate(csf_copy, csf_2, csf_super)
  call hexpr_shutdown()

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_csfset: ', nfail, ' check(s) FAILED'
    stop 1
  end if
  write(*,*) 'test_csfset                                        PASS'
  stop 0
end program test_csfset
