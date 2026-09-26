! F90 functional test for the binding-completeness procedures added to
! hexpr_mod: hexpr_csf_index, hexpr_csf_steps, hexpr_step_parse,
! hexpr_drt_show, hexpr_expr_show_full, hexpr_expr_show_one.
!
! Reference case: s_full_4o4e (SOCI, nve=4, spin=0, ncore=0, nval=4,
! next=0): norb=4, ncsfs=20.  *_show* tests only assert HEXPR_OK return;
! stdout content is not part of the contract.
!
! The csf2det section runs twice: once on that norb=4 case, where the
! values are small enough to write out by hand, and once on a 70-orbital
! FOCI DRT, where nwords == 2 and the mask word index is exercised.  The
! wide case cannot be reached by any of the natural reference cases,
! which all sit at norb <= 6.
program test_completeness
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr)                            :: drt, expr
  integer(c_int)                         :: st, i, idx, norb, ncsfs
  integer(c_int)                         :: nd, nwords, occ, want
  integer(c_int8_t), allocatable, target :: steps(:)
  integer(c_int8_t), target              :: parsed(4)
  integer(c_int64_t), allocatable, target :: alpha(:,:), beta(:,:)
  real(c_double),     allocatable, target :: coef(:)
  real(c_double)                         :: rhalf

  ! Wide (norb > 64) case; see the block near the end.
  type(c_ptr)                             :: wdrt
  integer(c_int)                          :: wnorb, wncsfs, wnwords
  integer(c_int)                          :: k, p, wnd, wbad
  integer(c_int8_t), allocatable, target  :: wsteps(:)
  integer(c_int64_t), allocatable, target :: walpha(:,:), wbeta(:,:)
  real(c_double),     allocatable, target :: wcoef(:)
  real(c_double)                          :: wnrm
  logical                                 :: whigh

  integer :: nfail

  nfail = 0

  st = hexpr_init()
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') 'FAIL: hexpr_init returned ', st
    stop 1
  end if

  drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4_c_int, 0_c_int, &
                          0_c_int, 4_c_int, 0_c_int)
  if (.not. c_associated(drt)) then
    write(*,'(a)') 'FAIL: hexpr_drt_create returned NULL'
    stop 1
  end if

  norb  = hexpr_drt_norb(drt)
  ncsfs = hexpr_drt_ncsfs(drt)
  if (norb /= 4)   nfail = nfail + 1
  if (ncsfs /= 20) nfail = nfail + 1

  ! step_parse("fude", 4) -> [3,1,2,0]
  st = hexpr_step_parse("fude", norb, parsed)
  if (st /= HEXPR_OK) nfail = nfail + 1
  if (.not. (parsed(1) == 3 .and. parsed(2) == 1 .and. &
             parsed(3) == 2 .and. parsed(4) == 0)) then
    write(*,'(a,4i3)') 'FAIL: step_parse got ', parsed
    nfail = nfail + 1
  end if

  ! csf_steps / csf_index round-trip for every CSF row
  allocate(steps(norb))
  do i = 0, ncsfs - 1
    st = hexpr_csf_steps(drt, i, steps)
    if (st /= HEXPR_OK) nfail = nfail + 1
    st = hexpr_csf_index(drt, steps, idx)
    if (st /= HEXPR_OK) nfail = nfail + 1
    if (idx /= i) then
      write(*,'(a,i0,a,i0)') 'FAIL: round-trip csf ', i, ' got ', idx
      nfail = nfail + 1
    end if
  end do
  deallocate(steps)

  ! csf2det: closed-shell 'ffee' expands into ONE determinant with
  ! orbitals 0 and 1 doubly occupied.  No open shells means no
  ! coupling and no sign convention, so the masks are forced.
  nwords = (norb + 63) / 64
  if (nwords /= 1) nfail = nfail + 1
  allocate(alpha(nwords, 2), beta(nwords, 2), coef(2))

  st = hexpr_step_parse("ffee", norb, parsed)
  if (st /= HEXPR_OK) nfail = nfail + 1
  st = hexpr_csf_index(drt, parsed, idx)
  if (st /= HEXPR_OK) nfail = nfail + 1
  st = hexpr_csf_ndets(drt, idx, 0_c_int, nd)
  if (st /= HEXPR_OK) nfail = nfail + 1
  if (nd /= 1) then
    write(*,'(a,i0)') 'FAIL: ndets(ffee) = ', nd
    nfail = nfail + 1
  end if
  st = hexpr_csf_to_dets(drt, idx, 0_c_int, nd, alpha, beta, coef)
  if (st /= HEXPR_OK) nfail = nfail + 1
  if (nd /= 1) nfail = nfail + 1
  ! Read the occupation back with BTEST, the only safe test: the mask
  ! is signed here, so a magnitude comparison would be wrong at bit 63.
  do i = 0, norb - 1
    occ = 0
    if (btest(alpha(ishft(i, -6) + 1, 1), iand(i, 63))) occ = occ + 1
    if (btest(beta (ishft(i, -6) + 1, 1), iand(i, 63))) occ = occ + 1
    want = 0
    if (i < 2) want = 2
    if (occ /= want) then
      write(*,'(a,i0,a,i0,a,i0)') 'FAIL: ffee orbital ', i, &
                                  ' occ ', occ, ' want ', want
      nfail = nfail + 1
    end if
  end do
  if (abs(abs(coef(1)) - 1.0_c_double) > 1.0e-12_c_double) then
    write(*,'(a,es23.15)') 'FAIL: ffee coef ', coef(1)
    nfail = nfail + 1
  end if

  ! csf2det: open-shell singlet 'udfe' expands into TWO determinants.
  ! Asserted: the count, the magnitudes, the SAME relative sign
  ! (singlet), and that the two determinants are each other's spin
  ! flip.  The sign of an individual coefficient is a convention and
  ! is deliberately not asserted.
  st = hexpr_step_parse("udfe", norb, parsed)
  if (st /= HEXPR_OK) nfail = nfail + 1
  st = hexpr_csf_index(drt, parsed, idx)
  if (st /= HEXPR_OK) nfail = nfail + 1
  st = hexpr_csf_ndets(drt, idx, 0_c_int, nd)
  if (st /= HEXPR_OK) nfail = nfail + 1
  if (nd /= 2) then
    write(*,'(a,i0)') 'FAIL: ndets(udfe) = ', nd
    nfail = nfail + 1
  else
    st = hexpr_csf_to_dets(drt, idx, 0_c_int, nd, alpha, beta, coef)
    if (st /= HEXPR_OK) nfail = nfail + 1
    rhalf = 1.0_c_double / sqrt(2.0_c_double)
    if (abs(abs(coef(1)) - rhalf) > 1.0e-12_c_double) nfail = nfail + 1
    if (abs(abs(coef(2)) - rhalf) > 1.0e-12_c_double) nfail = nfail + 1
    if (coef(1) * coef(2) <= 0.0_c_double) then
      write(*,'(a,2es23.15)') 'FAIL: udfe relative sign ', coef
      nfail = nfail + 1
    end if
    if (alpha(1,1) /= beta(1,2) .or. beta(1,1) /= alpha(1,2)) then
      write(*,'(a)') 'FAIL: udfe determinants are not a spin flip'
      nfail = nfail + 1
    end if
    if (alpha(1,1) == alpha(1,2)) then
      write(*,'(a)') 'FAIL: udfe determinants coincide'
      nfail = nfail + 1
    end if
  end if
  deallocate(alpha, beta, coef)

  ! ------------------------------------------------------------------
  ! csf2det, wide case: norb = 70, so nwords = 2 and the word index of
  ! the masks is finally exercised.
  !
  ! Everything above runs at norb = 4, where nwords == 1 and the first
  ! mask subscript is always 1; a wrong word index or a wrong stride
  ! cannot show up there.  This walks every CSF of a 70-orbital FOCI
  ! DRT and reconstructs the occupation of every spatial orbital from
  ! the masks, so a bit that lands in the wrong word is a mismatch
  ! rather than a silently different answer.  Same DRT as the Python,
  ! Julia and Ruby wide tests.
  !
  ! BTEST again: the masks are signed here, and orbitals >= 64 live in
  ! word 2 where bit 63 is the sign bit of that word.
  ! ------------------------------------------------------------------
  wdrt = hexpr_drt_create(HEXPR_DRT_FOCI, 2_c_int, 0_c_int, &
                          1_c_int, 1_c_int, 68_c_int)
  if (.not. c_associated(wdrt)) then
    write(*,'(a)') 'FAIL: wide hexpr_drt_create returned NULL'
    nfail = nfail + 1
  else
    wnorb   = hexpr_drt_norb(wdrt)
    wncsfs  = hexpr_drt_ncsfs(wdrt)
    wnwords = (wnorb + 63) / 64
    if (wnorb /= 70) then
      write(*,'(a,i0)') 'FAIL: wide norb ', wnorb
      nfail = nfail + 1
    end if
    if (wnwords /= 2) then
      write(*,'(a,i0)') 'FAIL: wide nwords ', wnwords
      nfail = nfail + 1
    end if

    allocate(wsteps(wnorb))
    whigh = .false.
    wbad  = 0

    do i = 0, wncsfs - 1
      st = hexpr_csf_steps(wdrt, i, wsteps)
      if (st /= HEXPR_OK) nfail = nfail + 1
      st = hexpr_csf_ndets(wdrt, i, 0_c_int, wnd)
      if (st /= HEXPR_OK) nfail = nfail + 1

      allocate(walpha(wnwords, wnd), wbeta(wnwords, wnd), wcoef(wnd))
      st = hexpr_csf_to_dets(wdrt, i, 0_c_int, nd, walpha, wbeta, wcoef)
      if (st /= HEXPR_OK) nfail = nfail + 1
      if (nd /= wnd) then
        write(*,'(a,i0,a,i0,a,i0)') 'FAIL: wide csf ', i, &
              ': ndets said ', wnd, ' but to_dets wrote ', nd
        nfail = nfail + 1
      end if

      wnrm = 0.0_c_double
      do k = 1, nd
        wnrm = wnrm + wcoef(k) * wcoef(k)
      end do
      if (abs(wnrm - 1.0_c_double) > 1.0e-12_c_double) then
        write(*,'(a,i0,a,es23.15)') 'FAIL: wide csf ', i, &
              ' norm ', wnrm
        nfail = nfail + 1
      end if

      do k = 1, nd
        do p = 0, wnorb - 1
          occ = 0
          if (btest(walpha(ishft(p, -6) + 1, k), iand(p, 63))) &
              occ = occ + 1
          if (btest(wbeta (ishft(p, -6) + 1, k), iand(p, 63))) &
              occ = occ + 1
          want = 0
          if (wsteps(p + 1) == HEXPR_STEP_U .or. &
              wsteps(p + 1) == HEXPR_STEP_D) want = 1
          if (wsteps(p + 1) == HEXPR_STEP_F) want = 2
          if (occ /= want) then
            if (wbad == 0) then
              write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') &
                'FAIL: wide csf ', i, ' det ', k, ' orbital ', p, &
                ' occ ', occ, ' want ', want
            end if
            wbad = wbad + 1
          end if
        end do
        if (wnwords >= 2) then
          if (walpha(wnwords, k) /= 0_c_int64_t .or. &
              wbeta (wnwords, k) /= 0_c_int64_t) whigh = .true.
        end if
      end do

      deallocate(walpha, wbeta, wcoef)
    end do

    deallocate(wsteps)

    if (wbad /= 0) then
      write(*,'(a,i0,a)') 'FAIL: wide occupation mismatches: ', wbad, &
            ' (first one shown above)'
      nfail = nfail + 1
    end if
    if (.not. whigh) then
      write(*,'(a)') 'FAIL: no determinant occupied an orbital >= 64;'&
                  // ' the wide path was not taken'
      nfail = nfail + 1
    end if

    call hexpr_drt_destroy(wdrt)
  end if

  ! *_show* : callable, return HEXPR_OK (no assertion on stdout content)
  expr = hexpr_expr_build(drt)
  if (.not. c_associated(expr)) then
    write(*,'(a)') 'FAIL: hexpr_expr_build returned NULL'
    nfail = nfail + 1
  else
    if (hexpr_drt_show(drt)        /= HEXPR_OK) nfail = nfail + 1
    if (hexpr_expr_show_full(expr) /= HEXPR_OK) nfail = nfail + 1
    if (hexpr_expr_show_one(expr)  /= HEXPR_OK) nfail = nfail + 1
    call hexpr_expr_destroy(expr)
  end if

  call hexpr_drt_destroy(drt)
  call hexpr_shutdown()

  if (nfail == 0) then
    write(*,'(a)') 'test_completeness                          PASS'
  else
    write(*,'(a,i0,a)') 'test_completeness: FAIL (', nfail, ' failures)'
    stop 1
  end if
end program test_completeness
