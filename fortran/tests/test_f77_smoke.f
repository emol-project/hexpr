C     F77 fixed-form smoke test: verify that the F77 aliases in
C     hexpr_fortran.c (hexpr_init_, hexpr_version_major_, hexpr_shutdown_)
C     link and run correctly from Fortran.
C
C     Fortran symbol mangling: gfortran appends '_' to the call name, so
C     a call to 'hexpr_init' maps to the symbol 'hexpr_init_' in the
C     shared library.  Do NOT add a trailing underscore in the Fortran
C     source -- the compiler adds it.
      program test_f77_smoke
      implicit none
      integer status, maj
      external hexpr_init, hexpr_version_major, hexpr_shutdown

      call hexpr_init(status)
      if (status .ne. 0) then
        write(*,*) 'FAIL: hexpr_init status=', status
        stop 1
      end if

      call hexpr_version_major(maj)
      if (maj .ne. 1) then
        write(*,*) 'FAIL: version_major=', maj, '(expected 1)'
        stop 1
      end if

      call hexpr_shutdown()

      write(*,*) 'test_f77_smoke                             PASS'
      stop 0
      end
