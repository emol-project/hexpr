! Test hexpr_rank1_pair (reduced matrix elements of (a+_p a~_q)^1, since 1.2)
! via the hexpr F90 module: the worked examples of the spin formulation, the
! total-spin identity sum_p R_pp = -sqrt(2) sqrt(S(S+1)(2S+1)) and the adjoint
! relation R_qp(ket,bra) = (-1)^(S'-S) R_pq(bra,ket), as in c/tests/test_rank1.c.
program test_rank1
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none

  integer :: nfail
  integer(c_int) :: st

  nfail = 0
  st = hexpr_init()
  call check(st == HEXPR_OK, "init")

  ! Worked examples of the per-orbital factors.
  call expect("udf", "fuu", 2, 0, sqrt(1.5d0), "singlet <- triplet, R_31")
  call expect("uuu", "udu", 0, 0, sqrt(2.0d0), "quartet <- doublet, R_11")
  call expect("ud",  "uu",  0, 0, -sqrt(1.5d0), "singlet <- triplet, R_11")
  call expect("ud",  "uu",  1, 1, sqrt(1.5d0), "singlet <- triplet, R_22")

  ! Total spin: sum_p R_pp(I, I) = -sqrt(2) sqrt(S(S+1)(2S+1)).
  call total_spin("uuu", 1.5d0)
  call total_spin("udu", 0.5d0)
  call total_spin("fuu", 1.0d0)
  call total_spin("uudu", 1.0d0)

  ! Adjoint relation on pairs of different S.
  call adjoint("udf", "fuu", 0, 2)
  call adjoint("uuu", "udu", 3, 1)
  call adjoint("uduu", "fuue", 2, 2)
  call adjoint("ufud", "uuuf", 1, 3)

  ! Bad input: an invalid walk (d before any u) is rejected.
  call check(rank1_status("du", "uu") /= HEXPR_OK, "invalid walk rejected")

  call hexpr_shutdown()
  if (nfail /= 0) then
    print '(a,i0,a)', "test_rank1: ", nfail, " failure(s)"
    stop 1
  end if
  print '(a)', "test_rank1: all checks passed"

contains

  subroutine check(ok, what)
    logical, intent(in) :: ok
    character(*), intent(in) :: what
    if (.not. ok) then
      print '(2a)', "FAIL: ", what
      nfail = nfail + 1
    end if
  end subroutine check

  function steps(s) result(v)
    character(*), intent(in) :: s
    integer(c_int8_t) :: v(len(s))
    integer :: i
    do i = 1, len(s)
      v(i) = int(index("eudf", s(i:i)) - 1, c_int8_t)
    end do
  end function steps

  ! Terms of <bra||T^1||ket>; n < 0 on a library error.
  subroutine terms(bra, ket, n, p, q, r)
    character(*), intent(in) :: bra, ket
    integer(c_int), intent(out) :: n
    integer(c_int32_t), allocatable, intent(out) :: p(:), q(:)
    real(c_double), allocatable, intent(out) :: r(:)
    integer(c_int8_t), allocatable, target :: b(:), k(:)
    integer(c_int) :: nmax, st
    nmax = max(hexpr_rank1_max_terms(int(len(bra), c_int)), 1)
    allocate(p(nmax), q(nmax), r(nmax))
    b = steps(bra)
    k = steps(ket)
    st = hexpr_rank1_pair(int(len(bra), c_int), b, k, n, p, q, r)
    if (st /= HEXPR_OK) n = -1
  end subroutine terms

  function rank1_status(bra, ket) result(st)
    character(*), intent(in) :: bra, ket
    integer(c_int) :: st, n
    integer(c_int8_t), allocatable, target :: b(:), k(:)
    integer(c_int32_t), allocatable, target :: p(:), q(:)
    real(c_double), allocatable, target :: r(:)
    integer(c_int) :: nmax
    nmax = max(hexpr_rank1_max_terms(int(len(bra), c_int)), 1)
    allocate(p(nmax), q(nmax), r(nmax))
    b = steps(bra)
    k = steps(ket)
    st = hexpr_rank1_pair(int(len(bra), c_int), b, k, n, p, q, r)
  end function rank1_status

  ! R_pq of the pair, 0 if absent.
  function value_of(bra, ket, ip, iq) result(v)
    character(*), intent(in) :: bra, ket
    integer, intent(in) :: ip, iq
    real(c_double) :: v
    integer(c_int) :: n
    integer(c_int32_t), allocatable :: p(:), q(:)
    real(c_double), allocatable :: r(:)
    integer :: i
    call terms(bra, ket, n, p, q, r)
    v = 0.0d0
    do i = 1, n
      if (p(i) == ip .and. q(i) == iq) v = r(i)
    end do
  end function value_of

  subroutine expect(bra, ket, ip, iq, val, what)
    character(*), intent(in) :: bra, ket, what
    integer, intent(in) :: ip, iq
    real(c_double), intent(in) :: val
    call check(abs(value_of(bra, ket, ip, iq) - val) < 1d-12, what)
  end subroutine expect

  subroutine total_spin(s, spin)
    character(*), intent(in) :: s
    real(c_double), intent(in) :: spin
    integer(c_int) :: n
    integer(c_int32_t), allocatable :: p(:), q(:)
    real(c_double), allocatable :: r(:)
    real(c_double) :: tot
    integer :: i
    call terms(s, s, n, p, q, r)
    tot = 0.0d0
    do i = 1, n
      if (p(i) == q(i)) tot = tot + r(i)
    end do
    call check(abs(tot + sqrt(2.0d0) * sqrt(spin * (spin + 1) * (2 * spin + 1))) &
               < 1d-12, "total spin of "//s)
  end subroutine total_spin

  ! sb, sk: 2S of bra and ket.
  subroutine adjoint(bra, ket, sb, sk)
    character(*), intent(in) :: bra, ket
    integer, intent(in) :: sb, sk
    integer(c_int) :: n, m
    integer(c_int32_t), allocatable :: p(:), q(:)
    real(c_double), allocatable :: r(:)
    real(c_double) :: sign
    integer :: i
    sign = 1.0d0
    if (mod(abs(sb - sk) / 2, 2) == 1) sign = -1.0d0
    call terms(bra, ket, n, p, q, r)
    m = 0
    do i = 1, n
      call check(abs(value_of(ket, bra, int(q(i)), int(p(i))) - sign * r(i)) &
                 < 1d-12, "adjoint "//bra//" "//ket)
      m = m + 1
    end do
    call check(n > 0 .and. m == n, "adjoint pair has terms: "//bra//" "//ket)
  end subroutine adjoint

end program test_rank1
