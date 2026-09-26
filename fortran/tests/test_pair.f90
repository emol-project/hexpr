! Test per-pair / block / subspace evaluation via the hexpr F90 module.
! Reference case: s_full_4o4e (SOCI, nve=4, spin=0, ncore=0, nval=4, next=0)
!   nc=20, nterms_full=638
!
! Checks:
!   T1: pair_count(0,1) returns HEXPR_OK and a non-negative count
!   T2: pair_max_terms returns HEXPR_OK and a positive bound
!   T3: pair_eval(0,1) returns n_terms == pair_count; all ij == triangle(0,1)==1
!   T4: aggregate test -- sweep all (bra,ket) pairs, sort+dedup, compare vs Expression
!   T5: block_eval bra_pos/ket_pos are within [0,n_bra) and [0,n_ket)
!   T6: block_eval per-(ip,jp) count matches pair_count
!   T7: subspace_eval matches block_eval(indices, indices)
program test_pair
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr)    :: drt, expr
  integer(c_int) :: st
  integer(c_int) :: nc, max_terms, n_count, n_terms
  integer(c_int) :: bra_i, ket_i
  integer        :: nfail, i, j, k, m

  ! Single-pair buffers (sized to max_terms after T2)
  integer(c_int32_t), allocatable, target :: ij_buf(:), pq_buf(:)
  real(c_double),     allocatable, target :: co_buf(:)

  ! Aggregate (T4)
  integer(c_int64_t) :: nterms_full
  integer(c_int32_t), allocatable, target :: ref_ij(:), ref_pq(:)
  real(c_double),     allocatable, target :: ref_co(:)
  integer(c_int32_t), allocatable :: all_ij(:), all_pq(:), ded_ij(:), ded_pq(:)
  real(c_double),     allocatable :: all_co(:), ded_co(:)
  integer :: n_all, n_ded

  ! Block (T5, T6)
  integer(c_int), target :: bra_idx(3), ket_idx(3)
  integer(c_int32_t), allocatable, target :: bp(:), kp(:), b_ij(:), b_pq(:)
  real(c_double),     allocatable, target :: b_co(:)
  integer(c_int) :: n_bra_t, n_ket_t, n_sub, pc

  ! Subspace (T7)
  integer(c_int), target :: sub_idx(3)
  integer(c_int32_t), allocatable, target :: s_bp(:), s_kp(:), s_ij(:), s_pq(:)
  real(c_double),     allocatable, target :: s_co(:)
  integer(c_int32_t), allocatable, target :: b2_bp(:), b2_kp(:), b2_ij(:), b2_pq(:)
  real(c_double),     allocatable, target :: b2_co(:)
  integer(c_int) :: ns, nb

  nfail = 0

  st = hexpr_init()
  if (st /= HEXPR_OK) then
    write(*,*) 'FAIL: hexpr_init'
    stop 1
  end if

  drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4_c_int, 0_c_int, &
                         0_c_int, 4_c_int, 0_c_int)
  if (.not. c_associated(drt)) then
    write(*,*) 'FAIL: hexpr_drt_create NULL'
    stop 1
  end if
  nc = hexpr_drt_ncsfs(drt)   ! 20

  ! ---- T1: pair_count basics ----------------------------------
  st = hexpr_pair_count(drt, 0_c_int, 1_c_int, HEXPR_EXPR_FULL, n_count)
  if (st /= HEXPR_OK .or. n_count < 0) then
    write(*,'(a,2i4)') '  T1 FAIL: pair_count st=', st, n_count
    nfail = nfail + 1
  else
    write(*,'(a,i0,a)') '  T1 pair_count(0,1)=', n_count, ': PASS'
  end if

  ! ---- T2: pair_max_terms -------------------------------------
  st = hexpr_pair_max_terms(drt, HEXPR_EXPR_FULL, max_terms)
  if (st /= HEXPR_OK .or. max_terms <= 0) then
    write(*,'(a,2i4)') '  T2 FAIL: pair_max_terms st=', st, max_terms
    nfail = nfail + 1
  else
    write(*,'(a,i0,a)') '  T2 pair_max_terms=', max_terms, ': PASS'
  end if

  ! ---- T3: pair_eval(0,1) -------------------------------------
  allocate(ij_buf(max_terms), pq_buf(max_terms), co_buf(max_terms))
  st = hexpr_pair_eval(drt, 0_c_int, 1_c_int, HEXPR_EXPR_FULL, &
                       n_terms, ij_buf, pq_buf, co_buf)
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') '  T3 FAIL: pair_eval st=', st
    nfail = nfail + 1
  else if (n_terms /= n_count) then
    write(*,'(a,2i4)') '  T3 FAIL: n_terms /= n_count: ', n_terms, n_count
    nfail = nfail + 1
  else if (n_terms < 0 .or. n_terms > max_terms) then
    write(*,'(a,i0)') '  T3 FAIL: n_terms out of range: ', n_terms
    nfail = nfail + 1
  else
    ! All ij values must equal triangle(0,1) = 1
    j = 0
    do k = 1, n_terms
      if (ij_buf(k) /= 1_c_int32_t) j = j + 1
    end do
    if (j > 0) then
      write(*,'(a,i0,a)') '  T3 FAIL: ', j, ' ij != 1'
      nfail = nfail + 1
    else
      write(*,'(a,i0,a)') '  T3 pair_eval(0,1): n_terms=', n_terms, ': PASS'
    end if
  end if
  deallocate(ij_buf, pq_buf, co_buf)

  ! ---- T4: aggregate test -------------------------------------
  expr = hexpr_expr_build(drt)
  if (.not. c_associated(expr)) then
    write(*,*) '  T4 FAIL: hexpr_expr_build NULL'
    nfail = nfail + 1
  else
    nterms_full = hexpr_expr_nterms_full(expr)
    allocate(ref_ij(nterms_full), ref_pq(nterms_full), ref_co(nterms_full))
    st = hexpr_expr_read_full(expr, ref_ij, ref_pq, ref_co)
    call hexpr_expr_destroy(expr)

    ! Sweep all pairs; accumulate raw terms
    n_all = 0
    allocate(all_ij(nc*nc*max_terms), all_pq(nc*nc*max_terms), all_co(nc*nc*max_terms))
    allocate(ij_buf(max_terms), pq_buf(max_terms), co_buf(max_terms))
    do bra_i = 0, nc-1
      do ket_i = 0, nc-1
        st = hexpr_pair_eval(drt, bra_i, ket_i, HEXPR_EXPR_FULL, &
                             n_terms, ij_buf, pq_buf, co_buf)
        if (st /= HEXPR_OK) then
          write(*,'(a,2i3)') '  T4 FAIL: pair_eval bra/ket=', bra_i, ket_i
          nfail = nfail + 1
        end if
        do k = 1, n_terms
          all_ij(n_all+k) = ij_buf(k)
          all_pq(n_all+k) = pq_buf(k)
          all_co(n_all+k) = co_buf(k)
        end do
        n_all = n_all + n_terms
      end do
    end do
    deallocate(ij_buf, pq_buf, co_buf)

    ! Sort all terms by (ij, pq)
    call sort_terms(n_all, all_ij(1:n_all), all_pq(1:n_all), all_co(1:n_all))

    ! Dedup: keep first occurrence of each (ij, pq) key
    allocate(ded_ij(n_all), ded_pq(n_all), ded_co(n_all))
    n_ded = 0
    do k = 1, n_all
      if (k == 1 .or. all_ij(k) /= all_ij(k-1) .or. all_pq(k) /= all_pq(k-1)) then
        n_ded = n_ded + 1
        ded_ij(n_ded) = all_ij(k)
        ded_pq(n_ded) = all_pq(k)
        ded_co(n_ded) = all_co(k)
      end if
    end do
    deallocate(all_ij, all_pq, all_co)

    ! Sort reference by (ij, pq) too
    call sort_terms(int(nterms_full), ref_ij, ref_pq, ref_co)

    if (n_ded /= int(nterms_full)) then
      write(*,'(a,i0,a,i0)') '  T4 FAIL: n_dedup=', n_ded, &
                              ' expected=', int(nterms_full)
      nfail = nfail + 1
    else
      j = 0
      do i = 1, n_ded
        if (ded_ij(i) /= ref_ij(i) .or. ded_pq(i) /= ref_pq(i)) then
          j = j + 1
        else if (abs(ded_co(i) - ref_co(i)) > 1.0e-10_c_double) then
          j = j + 1
        end if
      end do
      if (j > 0) then
        write(*,'(a,i0,a)') '  T4 FAIL: ', j, ' term mismatch(es)'
        nfail = nfail + 1
      else
        write(*,'(a,i0,a)') '  T4 aggregate (', n_ded, ' terms): PASS'
      end if
    end if
    deallocate(ref_ij, ref_pq, ref_co, ded_ij, ded_pq, ded_co)
  end if

  ! ---- T5: block_eval bra_pos/ket_pos within-array positions --
  bra_idx = [5_c_int, 9_c_int, 12_c_int]
  ket_idx  = [2_c_int, 7_c_int, 3_c_int]
  n_bra_t = 3_c_int;  n_ket_t = 3_c_int
  st = hexpr_block_max_terms(drt, n_bra_t, n_ket_t, HEXPR_EXPR_FULL, max_terms)
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') '  T5 FAIL: block_max_terms st=', st
    nfail = nfail + 1
  else
    allocate(bp(max(max_terms,1)), kp(max(max_terms,1)), &
             b_ij(max(max_terms,1)), b_pq(max(max_terms,1)), b_co(max(max_terms,1)))
    st = hexpr_block_eval(drt, bra_idx, n_bra_t, ket_idx, n_ket_t, &
                          HEXPR_EXPR_FULL, n_terms, bp, kp, b_ij, b_pq, b_co)
    if (st /= HEXPR_OK) then
      write(*,'(a,i0)') '  T5 FAIL: block_eval st=', st
      nfail = nfail + 1
    else
      j = 0
      do k = 1, n_terms
        if (bp(k) < 0_c_int32_t .or. bp(k) >= n_bra_t) j = j + 1
        if (kp(k) < 0_c_int32_t .or. kp(k) >= n_ket_t) j = j + 1
      end do
      if (j > 0) then
        write(*,'(a,i0,a)') '  T5 FAIL: ', j, ' pos value(s) out of range'
        nfail = nfail + 1
      else
        write(*,'(a,i0,a)') '  T5 block bra/ket_pos range (', n_terms, ' terms): PASS'
      end if
    end if
    deallocate(bp, kp, b_ij, b_pq, b_co)
  end if

  ! ---- T6: block_eval per-(ip,jp) count matches pair_count ----
  bra_idx(1:2) = [0_c_int, 1_c_int]
  ket_idx(1:2) = [2_c_int, 3_c_int]
  n_bra_t = 2_c_int;  n_ket_t = 2_c_int
  st = hexpr_block_max_terms(drt, n_bra_t, n_ket_t, HEXPR_EXPR_FULL, max_terms)
  allocate(bp(max(max_terms,1)), kp(max(max_terms,1)), &
           b_ij(max(max_terms,1)), b_pq(max(max_terms,1)), b_co(max(max_terms,1)))
  st = hexpr_block_eval(drt, bra_idx(1:2), n_bra_t, ket_idx(1:2), n_ket_t, &
                        HEXPR_EXPR_FULL, n_terms, bp, kp, b_ij, b_pq, b_co)
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') '  T6 FAIL: block_eval st=', st
    nfail = nfail + 1
  else
    j = 0
    do i = 0, n_bra_t-1
      do m = 0, n_ket_t-1
        ! Count block terms for position (i, m)
        n_sub = 0
        do k = 1, n_terms
          if (bp(k) == int(i, c_int32_t) .and. kp(k) == int(m, c_int32_t)) &
              n_sub = n_sub + 1
        end do
        ! Compare with pair_count
        st = hexpr_pair_count(drt, bra_idx(i+1), ket_idx(m+1), &
                              HEXPR_EXPR_FULL, pc)
        if (n_sub /= pc) then
          write(*,'(a,2i2,a,2i4)') '  T6 FAIL: (ip,jp)=', i, m, &
                                    ' sub/pair=', n_sub, pc
          j = j + 1
        end if
      end do
    end do
    if (j > 0) then
      nfail = nfail + 1
    else
      write(*,'(a,i0,a)') '  T6 block matches pair_count (', n_terms, ' total): PASS'
    end if
  end if
  deallocate(bp, kp, b_ij, b_pq, b_co)

  ! ---- T7: subspace_eval == block_eval(idx, idx) --------------
  sub_idx = [0_c_int, 3_c_int, 7_c_int]
  st = hexpr_subspace_max_terms(drt, 3_c_int, HEXPR_EXPR_FULL, max_terms)
  allocate(s_bp(max(max_terms,1)), s_kp(max(max_terms,1)), &
           s_ij(max(max_terms,1)), s_pq(max(max_terms,1)), s_co(max(max_terms,1)))
  st = hexpr_subspace_eval(drt, sub_idx, 3_c_int, HEXPR_EXPR_FULL, &
                            ns, s_bp, s_kp, s_ij, s_pq, s_co)
  if (st /= HEXPR_OK) then
    write(*,'(a,i0)') '  T7 FAIL: subspace_eval st=', st
    nfail = nfail + 1
  else
    st = hexpr_block_max_terms(drt, 3_c_int, 3_c_int, HEXPR_EXPR_FULL, max_terms)
    allocate(b2_bp(max(max_terms,1)), b2_kp(max(max_terms,1)), &
             b2_ij(max(max_terms,1)), b2_pq(max(max_terms,1)), b2_co(max(max_terms,1)))
    st = hexpr_block_eval(drt, sub_idx, 3_c_int, sub_idx, 3_c_int, &
                          HEXPR_EXPR_FULL, nb, b2_bp, b2_kp, b2_ij, b2_pq, b2_co)
    if (st /= HEXPR_OK) then
      write(*,'(a,i0)') '  T7 FAIL: block_eval st=', st
      nfail = nfail + 1
    else if (ns /= nb) then
      write(*,'(a,2i4)') '  T7 FAIL: subspace/block count=', ns, nb
      nfail = nfail + 1
    else
      j = 0
      do k = 1, ns
        if (s_ij(k)  /= b2_ij(k)  .or. s_pq(k) /= b2_pq(k) .or. &
            s_bp(k)  /= b2_bp(k)  .or. s_kp(k) /= b2_kp(k) .or. &
            abs(s_co(k) - b2_co(k)) > 1.0e-14_c_double) j = j + 1
      end do
      if (j > 0) then
        write(*,'(a,i0,a)') '  T7 FAIL: ', j, ' mismatch(es) vs block_eval'
        nfail = nfail + 1
      else
        write(*,'(a,i0,a)') '  T7 subspace == block_eval (', ns, ' terms): PASS'
      end if
    end if
    deallocate(b2_bp, b2_kp, b2_ij, b2_pq, b2_co)
  end if
  deallocate(s_bp, s_kp, s_ij, s_pq, s_co)

  call hexpr_drt_destroy(drt)
  call hexpr_shutdown()

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_pair: ', nfail, ' check(s) FAILED'
    stop 1
  end if
  write(*,*) 'test_pair                                          PASS'
  stop 0

contains

  ! Shell sort on three parallel arrays by (ij, pq) key.
  subroutine sort_terms(n, ij_a, pq_a, co_a)
    integer,            intent(in)    :: n
    integer(c_int32_t), intent(inout) :: ij_a(n), pq_a(n)
    real(c_double),     intent(inout) :: co_a(n)
    integer            :: gap, ii, jj
    integer(c_int32_t) :: tij, tpq
    real(c_double)     :: tco
    gap = n / 2
    do while (gap > 0)
      do ii = gap + 1, n
        tij = ij_a(ii);  tpq = pq_a(ii);  tco = co_a(ii)
        jj = ii
        do while (jj > gap)
          if (ij_a(jj-gap) > tij .or. &
              (ij_a(jj-gap) == tij .and. pq_a(jj-gap) > tpq)) then
            ij_a(jj) = ij_a(jj-gap)
            pq_a(jj) = pq_a(jj-gap)
            co_a(jj) = co_a(jj-gap)
            jj = jj - gap
          else
            exit
          end if
        end do
        ij_a(jj) = tij;  pq_a(jj) = tpq;  co_a(jj) = tco
      end do
      gap = gap / 2
    end do
  end subroutine sort_terms

end program test_pair
