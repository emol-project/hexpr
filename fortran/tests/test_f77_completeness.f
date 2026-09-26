C     F77 functional test for the binding-completeness aliases added in
C     hexpr_fortran.c: hexpr_drt_from_csfset_, hexpr_csf_index_,
C     hexpr_csf_steps_, hexpr_step_parse_, and the scalar DRT getters
C     hexpr_drt_ncore_/nval_/next_/nve_state_/spin_x2_.
C
C     gfortran appends '_' to call names, so 'hexpr_init' maps to the
C     library symbol 'hexpr_init_'.  Reference case: s_full_4o4e
C     (SOCI, nve=4, spin=0, ncore=0, nval=4, next=0): norb=4, ncsfs=20.
C
C     The csf2det section runs twice: on that norb=4 case, and on a
C     70-orbital FOCI DRT where nwords = 2 and the mask word index is
C     exercised.  No natural reference case reaches norb > 64.
      program test_f77_completeness
      implicit none
      integer status, i, j, idx, norb, ncsfs, nfail
      integer ncore, nval, nextv, spinx2, nvest
      integer kind, nve, spin, nc, nv, nx, n4, nsub, subnc
      integer ms, nd
      integer*8 drt, sub
      integer*8 alpha(2), beta(2)
      double precision coef(2), rhalf
      integer*1 step(4), parsed(4), buf(8)
C     Wide (norb > 64) case; see the block near the end.  F77 cannot
C     ALLOCATE, so the determinant buffers are sized from the case:
C     nve=2 with spin=0 gives at most two open shells, hence ndets <= 2.
C     The arrays hold 4 and the count is checked before they are used.
      integer*8 wdrt
      integer*8 walpha(2,4), wbeta(2,4)
      double precision wcoef(4), wnrm
      integer*1 wstep(70)
      integer wnorb, wncsf, wnw, wnd, wcnt, wocc, wwant, wbad
      integer k, p
      logical whigh
      integer*8 hexpr_drt_create, hexpr_drt_from_csfset
      external hexpr_init, hexpr_shutdown
      external hexpr_drt_create, hexpr_drt_from_csfset
      external hexpr_drt_norb, hexpr_drt_ncsfs
      external hexpr_drt_ncore, hexpr_drt_nval, hexpr_drt_next
      external hexpr_drt_spin_x2, hexpr_drt_nve_state
      external hexpr_csf_steps, hexpr_csf_index, hexpr_step_parse
      external hexpr_csf_ndets, hexpr_csf_to_dets
      external hexpr_drt_destroy

      nfail = 0

      call hexpr_init(status)
      if (status .ne. 0) then
        write(*,*) 'FAIL: hexpr_init status=', status
        stop 1
      end if

      kind = 0
      nve  = 4
      spin = 0
      nc   = 0
      nv   = 4
      nx   = 0
      drt = hexpr_drt_create(kind, nve, spin, nc, nv, nx)
      if (drt .eq. 0) then
        write(*,*) 'FAIL: hexpr_drt_create returned NULL'
        stop 1
      end if

      call hexpr_drt_norb(drt, norb)
      call hexpr_drt_ncsfs(drt, ncsfs)
      if (norb .ne. 4)   nfail = nfail + 1
      if (ncsfs .ne. 20) nfail = nfail + 1

C     Scalar DRT getters
      call hexpr_drt_ncore(drt, ncore)
      call hexpr_drt_nval(drt, nval)
      call hexpr_drt_next(drt, nextv)
      call hexpr_drt_spin_x2(drt, spinx2)
      call hexpr_drt_nve_state(drt, nvest)
      if (ncore  .ne. 0) nfail = nfail + 1
      if (nval   .ne. 4) nfail = nfail + 1
      if (nextv  .ne. 0) nfail = nfail + 1
      if (spinx2 .ne. 0) nfail = nfail + 1
      if (nvest  .ne. 4) nfail = nfail + 1

C     step_parse("fude", 4) -> [3,1,2,0]
      n4 = 4
      call hexpr_step_parse('fude', n4, parsed, status)
      if (status .ne. 0) nfail = nfail + 1
      if (parsed(1).ne.3 .or. parsed(2).ne.1 .or.
     &    parsed(3).ne.2 .or. parsed(4).ne.0) then
        write(*,*) 'FAIL: step_parse got', parsed
        nfail = nfail + 1
      end if

C     csf_steps / csf_index round-trip for every CSF row
      do i = 0, ncsfs - 1
        call hexpr_csf_steps(drt, i, step, status)
        if (status .ne. 0) nfail = nfail + 1
        call hexpr_csf_index(drt, step, idx, status)
        if (status .ne. 0) nfail = nfail + 1
        if (idx .ne. i) then
          write(*,*) 'FAIL: round-trip csf', i, 'got', idx
          nfail = nfail + 1
        end if
      end do

C     csf2det: closed-shell CSF 'ffee' expands into ONE determinant,
C     with orbitals 1 and 2 doubly occupied.  No open shells means no
C     coupling and no sign convention -- the masks are forced.
      n4 = 4
      ms = 0
      call hexpr_step_parse('ffee', n4, parsed, status)
      if (status .ne. 0) nfail = nfail + 1
      call hexpr_csf_index(drt, parsed, idx, status)
      if (status .ne. 0) nfail = nfail + 1
      call hexpr_csf_ndets(drt, idx, ms, nd, status)
      if (status .ne. 0) nfail = nfail + 1
      if (nd .ne. 1) then
        write(*,*) 'FAIL: ndets(ffee) =', nd
        nfail = nfail + 1
      end if
      call hexpr_csf_to_dets(drt, idx, ms, nd, alpha, beta, coef,
     &                       status)
      if (status .ne. 0) nfail = nfail + 1
      if (nd .ne. 1) nfail = nfail + 1
      if (alpha(1) .ne. 3 .or. beta(1) .ne. 3) then
        write(*,*) 'FAIL: ffee masks', alpha(1), beta(1)
        nfail = nfail + 1
      end if
      if (abs(abs(coef(1)) - 1.0d0) .gt. 1.0d-12) then
        write(*,*) 'FAIL: ffee coef', coef(1)
        nfail = nfail + 1
      end if

C     csf2det: open-shell singlet 'udfe' expands into TWO determinants.
C     Asserted here: the count, the magnitudes, the SAME relative sign
C     (singlet), and that the two determinants are each other's spin
C     flip.  The sign of an individual coefficient is a convention and
C     is deliberately not asserted.
      rhalf = 1.0d0 / dsqrt(2.0d0)
      call hexpr_step_parse('udfe', n4, parsed, status)
      if (status .ne. 0) nfail = nfail + 1
      call hexpr_csf_index(drt, parsed, idx, status)
      if (status .ne. 0) nfail = nfail + 1
      call hexpr_csf_ndets(drt, idx, ms, nd, status)
      if (status .ne. 0) nfail = nfail + 1
      if (nd .ne. 2) then
        write(*,*) 'FAIL: ndets(udfe) =', nd
        nfail = nfail + 1
      else
        call hexpr_csf_to_dets(drt, idx, ms, nd, alpha, beta, coef,
     &                         status)
        if (status .ne. 0) nfail = nfail + 1
        if (abs(abs(coef(1)) - rhalf) .gt. 1.0d-12) nfail = nfail + 1
        if (abs(abs(coef(2)) - rhalf) .gt. 1.0d-12) nfail = nfail + 1
        if (coef(1) * coef(2) .le. 0.0d0) then
          write(*,*) 'FAIL: udfe relative sign', coef(1), coef(2)
          nfail = nfail + 1
        end if
        if (alpha(1) .ne. beta(2) .or. beta(1) .ne. alpha(2)) then
          write(*,*) 'FAIL: udfe not a spin flip', alpha(1), beta(1),
     &               alpha(2), beta(2)
          nfail = nfail + 1
        end if
        if (alpha(1) .eq. alpha(2)) then
          write(*,*) 'FAIL: udfe determinants coincide'
          nfail = nfail + 1
        end if
      end if

C     csf2det, wide case: norb = 70, so nwords = 2 and the mask word
C     index is finally exercised.  Everything above runs at norb = 4,
C     where nwords = 1 and the first mask subscript is always 1, so a
C     wrong word index or stride cannot show up.  This walks every CSF
C     and reconstructs each spatial orbital's occupation from the
C     masks.  Same DRT as the Python, Julia and Ruby wide tests.
      kind = 1
      nve  = 2
      spin = 0
      nc   = 1
      nv   = 1
      nx   = 68
      ms   = 0
      wdrt = hexpr_drt_create(kind, nve, spin, nc, nv, nx)
      if (wdrt .eq. 0) then
        write(*,*) 'FAIL: wide hexpr_drt_create returned NULL'
        nfail = nfail + 1
      else
        call hexpr_drt_norb(wdrt, wnorb)
        call hexpr_drt_ncsfs(wdrt, wncsf)
        wnw = (wnorb + 63) / 64
        if (wnorb .ne. 70) then
          write(*,*) 'FAIL: wide norb', wnorb
          nfail = nfail + 1
        end if
        if (wnw .ne. 2) then
          write(*,*) 'FAIL: wide nwords', wnw
          nfail = nfail + 1
        end if
        whigh = .false.
        wbad  = 0
        do i = 0, wncsf - 1
          call hexpr_csf_steps(wdrt, i, wstep, status)
          if (status .ne. 0) nfail = nfail + 1
          call hexpr_csf_ndets(wdrt, i, ms, wnd, status)
          if (status .ne. 0) nfail = nfail + 1
          if (wnd .gt. 4) then
            if (wbad .eq. 0) write(*,*) 'FAIL: wide ndets', i, wnd
            wbad = wbad + 1
          else
            call hexpr_csf_to_dets(wdrt, i, ms, wcnt, walpha, wbeta,
     &                             wcoef, status)
            if (status .ne. 0) nfail = nfail + 1
            if (wcnt .ne. wnd) then
              if (wbad .eq. 0) write(*,*) 'FAIL: wide count', i,
     &                                    wnd, wcnt
              wbad = wbad + 1
            end if
            wnrm = 0.0d0
            do k = 1, wcnt
              wnrm = wnrm + wcoef(k) * wcoef(k)
            end do
            if (abs(wnrm - 1.0d0) .gt. 1.0d-12) then
              if (wbad .eq. 0) write(*,*) 'FAIL: wide norm', i, wnrm
              wbad = wbad + 1
            end if
            do k = 1, wcnt
              do p = 0, wnorb - 1
                wocc = 0
                if (btest(walpha(ishft(p,-6)+1,k), iand(p,63)))
     &            wocc = wocc + 1
                if (btest(wbeta(ishft(p,-6)+1,k), iand(p,63)))
     &            wocc = wocc + 1
                wwant = 0
                if (wstep(p+1) .eq. 1 .or. wstep(p+1) .eq. 2)
     &            wwant = 1
                if (wstep(p+1) .eq. 3) wwant = 2
                if (wocc .ne. wwant) then
                  if (wbad .eq. 0) write(*,*) 'FAIL: wide occ csf',
     &              i, ' det', k, ' orb', p, wocc, wwant
                  wbad = wbad + 1
                end if
              end do
              if (walpha(2,k) .ne. 0 .or. wbeta(2,k) .ne. 0)
     &          whigh = .true.
            end do
          end if
        end do
        if (wbad .ne. 0) then
          write(*,*) 'FAIL: wide mismatches', wbad
          nfail = nfail + 1
        end if
        if (.not. whigh) then
          write(*,*) 'FAIL: no determinant reached orbital 64 or above'
          nfail = nfail + 1
        end if
        call hexpr_drt_destroy(wdrt)
      end if

C     from_csfset: build a 2-CSF sub-DRT from the first two CSF rows
      call hexpr_csf_steps(drt, 0, step, status)
      do j = 1, 4
        buf(j) = step(j)
      end do
      call hexpr_csf_steps(drt, 1, step, status)
      do j = 1, 4
        buf(4 + j) = step(j)
      end do
      nsub = 2
      sub = hexpr_drt_from_csfset(norb, buf, nsub)
      if (sub .eq. 0) then
        write(*,*) 'FAIL: hexpr_drt_from_csfset returned NULL'
        nfail = nfail + 1
      else
        call hexpr_drt_ncsfs(sub, subnc)
        if (subnc .lt. 2) nfail = nfail + 1
        call hexpr_drt_destroy(sub)
      end if

      call hexpr_drt_destroy(drt)
      call hexpr_shutdown()

      if (nfail .eq. 0) then
        write(*,*) 'test_f77_completeness                      PASS'
        stop 0
      else
        write(*,*) 'test_f77_completeness FAIL nfail=', nfail
        stop 1
      end if
      end
