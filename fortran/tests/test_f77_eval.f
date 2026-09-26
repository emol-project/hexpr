C     F77 fixed-form end-to-end test: drive the full CSFSET -> evaluation
C     path through the trailing-underscore wrapper ABI only (no module USE).
C
C     It builds a 2-CSF sub-DRT with hexpr_drt_from_csfset_ and calls
C     hexpr_pair_count_ / hexpr_pair_max_terms_ / hexpr_pair_eval_ on three
C     CSF pairs, asserting the exact (ij,pqrs,coef) terms.
C
C     INDEX SPACE.  The pair indices passed below are INPUT positions, 0 and
C     1, the same space as hexpr_csf_index_ and hexpr_csf_steps_.  This case
C     is chosen so the two spaces differ: input row 0 is canonical 3 and
C     input row 1 is canonical 0.  The test asserts that, via the
C     hexpr_drt_walk_csfs_ buffer, so a build where the library stopped
C     translating would not slip through as a coincidence.  The expected ij
C     values stay CANONICAL (ij is triangle() of Shavitt path weights, per
C     docs/PQRS_INDEX_SPEC.md), which is why they are 9, 6 and 0 rather than
C     the 0, 1 and 2 the input positions would give.
C
C     Reference case O2 (norb=4, ncsfs=2):
C       input row 0 steps = 1 0 3 2   (u e f d)
C       input row 1 steps = 0 1 2 3   (e u d f)
C     Expected terms are the frozen oracle reference_csfset/O2.txt
C     (which=full), which c/tests/test_csfset_pair verifies bit-exact against
C     the current library ("O2: pair_eval full == file blocks  PASS").
C     Coefficients are all +-1 or +-2 to within ~1e-15; compared here to a
C     1e-9 absolute tolerance.
C
C     gfortran appends '_' to call names, so CALL hexpr_init maps to the
C     library symbol hexpr_init_.  Do NOT add a trailing underscore in the
C     source.  All arguments are passed by reference; the DRT handle is an
C     INTEGER*8, ij/pqrs are INTEGER*4, coef is DOUBLE PRECISION, and CSF
C     indices are 0-based.
      program test_f77_eval
      implicit none
      integer status, norb, ncsfs, nin, m, mt, nt
      integer which, c0, c1, w0, w1, nfail
      integer*8 drt
      integer*8 hexpr_drt_from_csfset
      external hexpr_init, hexpr_shutdown, hexpr_drt_destroy
      external hexpr_drt_from_csfset
      external hexpr_drt_ncsfs, hexpr_drt_walk_ncsfs
      external hexpr_drt_walk_csfs
      external hexpr_pair_count, hexpr_pair_max_terms, hexpr_pair_eval

C     Input CSF set (case O2): row0 in csfs(1:4), row1 in csfs(5:8).
      integer*1 csfs(8), wbuf(64)
C     Emitted-term buffers (sized well above any 4-orbital max_terms).
      integer ij(256), pq(256)
      double precision co(256)

C     Expected terms from reference_csfset/O2.txt (which=full).
C     ij is canonical: row0 is canonical 3, row1 is canonical 0, so
C     triangle(3,3)=9, triangle(3,0)=6, triangle(0,0)=0.
C     Pair (row0,row0): ij=9, 10 terms.
      integer e00ij, e00pq(10)
      double precision e00co(10)
C     Pair (row0,row1): ij=6, 2 terms.
      integer e01ij, e01pq(2)
      double precision e01co(2)
C     Pair (row1,row1): ij=0, 10 terms.
      integer e11ij, e11pq(10)
      double precision e11co(10)

      data csfs /1, 0, 3, 2,  0, 1, 2, 3/

      data e00ij /9/
      data e00pq /0, 5, 9, 19, 25, 30, 37, 54, 55, 60/
      data e00co /1.0d0, 2.0d0, 1.0d0, -1.0d0, 2.0d0,
     &            1.0d0, 1.0d0, -1.0d0, 1.0d0, 2.0d0/

      data e01ij /6/
      data e01pq /35, 47/
      data e01co /2.0d0, -1.0d0/

      data e11ij /0/
      data e11pq /2, 5, 9, 24, 27, 45, 54, 57, 60, 64/
      data e11co /1.0d0, 1.0d0, 2.0d0, 1.0d0, 1.0d0,
     &            -1.0d0, -1.0d0, 2.0d0, 2.0d0, 1.0d0/

      nfail = 0
      norb  = 4
      nin   = 2
      which = 0

      call hexpr_init(status)
      if (status .ne. 0) then
        write(*,*) 'FAIL: hexpr_init status=', status
        stop 1
      end if

C     ---- Build the CSFSET sub-DRT -----------------------------------
      drt = hexpr_drt_from_csfset(norb, csfs, nin)
      if (drt .eq. 0) then
        write(*,*) 'FAIL: hexpr_drt_from_csfset returned NULL'
        call hexpr_shutdown()
        stop 1
      end if

C     Input-order count must equal the supplied set size.
      call hexpr_drt_ncsfs(drt, ncsfs)
      if (ncsfs .ne. 2) then
        write(*,*) 'FAIL: ncsfs=', ncsfs, ' expected 2'
        nfail = nfail + 1
      end if

C     ---- Pair indices are INPUT positions --------------------------
      c0 = 0
      c1 = 1

C     ---- Confirm the two index spaces really differ here -----------
C     Not used for indexing: this pins the fixture.  If input order ever
C     coincided with canonical order the pair calls would still pass while
C     proving nothing about the library's translation.
      call hexpr_drt_walk_ncsfs(drt, m)
      if (m .lt. 1 .or. m * norb .gt. 64) then
        write(*,*) 'FAIL: walk_ncsfs=', m
        nfail = nfail + 1
      else
        call hexpr_drt_walk_csfs(drt, wbuf)
        call findcx(wbuf, m, norb, csfs(1), w0)
        call findcx(wbuf, m, norb, csfs(5), w1)
        if (w0 .ne. 3 .or. w1 .ne. 0) then
          write(*,*) 'FAIL: canonical indices w0/w1=', w0, w1,
     &               ' expected 3 and 0'
          nfail = nfail + 1
        end if
      end if

C     ---- pair_max_terms upper bound --------------------------------
      call hexpr_pair_max_terms(drt, which, mt, status)
      if (status .ne. 0 .or. mt .le. 0 .or. mt .gt. 256) then
        write(*,*) 'FAIL: pair_max_terms status/mt=', status, mt
        nfail = nfail + 1
      end if

      if (nfail .eq. 0) then
C       ---- Pair (input 0, input 0): 10 terms, canonical ij=9 -------
        call hexpr_pair_count(drt, c0, c0, which, nt, status)
        if (status .ne. 0 .or. nt .ne. 10) then
          write(*,*) 'FAIL: pair_count(0,0)=', nt, ' status=', status
          nfail = nfail + 1
        end if
        call hexpr_pair_eval(drt, c0, c0, which, nt,
     &                       ij, pq, co, status)
        if (status .ne. 0) then
          write(*,*) 'FAIL: pair_eval(0,0) status=', status
          nfail = nfail + 1
        end if
        call chkpr(1, nt, ij, pq, co, 10, e00ij, e00pq, e00co, nfail)

C       ---- Pair (input 0, input 1): 2 terms, canonical ij=6 --------
        call hexpr_pair_count(drt, c0, c1, which, nt, status)
        if (status .ne. 0 .or. nt .ne. 2) then
          write(*,*) 'FAIL: pair_count(0,1)=', nt, ' status=', status
          nfail = nfail + 1
        end if
        call hexpr_pair_eval(drt, c0, c1, which, nt,
     &                       ij, pq, co, status)
        if (status .ne. 0) then
          write(*,*) 'FAIL: pair_eval(0,1) status=', status
          nfail = nfail + 1
        end if
        call chkpr(2, nt, ij, pq, co, 2, e01ij, e01pq, e01co, nfail)

C       ---- Pair (input 1, input 1): 10 terms, canonical ij=0 -------
        call hexpr_pair_count(drt, c1, c1, which, nt, status)
        if (status .ne. 0 .or. nt .ne. 10) then
          write(*,*) 'FAIL: pair_count(1,1)=', nt, ' status=', status
          nfail = nfail + 1
        end if
        call hexpr_pair_eval(drt, c1, c1, which, nt,
     &                       ij, pq, co, status)
        if (status .ne. 0) then
          write(*,*) 'FAIL: pair_eval(1,1) status=', status
          nfail = nfail + 1
        end if
        call chkpr(3, nt, ij, pq, co, 10, e11ij, e11pq, e11co, nfail)
      end if

      call hexpr_drt_destroy(drt)
      call hexpr_shutdown()

      if (nfail .eq. 0) then
        write(*,*) 'test_f77_eval                              PASS'
        stop 0
      else
        write(*,*) 'test_f77_eval FAIL nfail=', nfail
        stop 1
      end if
      end

C     ----------------------------------------------------------------
C     Canonical (walk-order) index of an input row: locate the length-
C     norb row in the row-major hexpr_drt_walk_csfs_ buffer.  Used ONLY to
C     assert that input order differs from canonical order for this case;
C     the pair calls address CSFs by input position and the library does
C     the translation.  Returns -1 if not found.
      subroutine findcx(wbuf, m, norb, row, cidx)
      implicit none
      integer m, norb, cidx
      integer*1 wbuf(*), row(*)
      integer i, b, ok
      cidx = -1
      do i = 0, m - 1
        ok = 1
        do b = 1, norb
          if (wbuf(i * norb + b) .ne. row(b)) ok = 0
        end do
        if (ok .eq. 1 .and. cidx .lt. 0) cidx = i
      end do
      return
      end

C     ----------------------------------------------------------------
C     Assert that the emitted (ij,pqrs,coef) terms of one pair match the
C     expected set (order-independent): exactly nexp terms, every emitted
C     ij == eij, and each expected pqrs appears once with |coef| matching
C     to 1e-9.  Increments nfail on any discrepancy.
      subroutine chkpr(tag, nt, ij, pq, co, nexp, eij, epq, eco,
     &                 nfail)
      implicit none
      integer tag, nt, nexp, eij, nfail
      integer ij(*), pq(*), epq(*)
      double precision co(*), eco(*)
      integer k, mm, found, bad
      bad = 0
      if (nt .ne. nexp) then
        write(*,*) 'FAIL pair', tag, ' nterms=', nt, ' expected', nexp
        bad = bad + 1
      end if
      do k = 1, nt
        if (ij(k) .ne. eij) then
          write(*,*) 'FAIL pair', tag, ' ij=', ij(k), ' expected', eij
          bad = bad + 1
        end if
      end do
      do mm = 1, nexp
        found = 0
        do k = 1, nt
          if (pq(k) .eq. epq(mm) .and.
     &        abs(co(k) - eco(mm)) .lt. 1.0d-9) found = found + 1
        end do
        if (found .ne. 1) then
          write(*,*) 'FAIL pair', tag, ' pqrs=', epq(mm),
     &               ' matches=', found
          bad = bad + 1
        end if
      end do
      if (bad .eq. 0) then
        write(*,*) '  pair', tag, ' (', nt, ' terms): PASS'
      end if
      nfail = nfail + bad
      return
      end
