# hexpr

C library and multi-language bindings for generating expressions of the
CSF-based Hamiltonian matrix used in quantum chemistry CI calculations.
Input CSFs can be supplied directly, or alternatively the full CSF set of a
first-order CI (FOCI) or second-order CI (SOCI) space can be specified by
its orbital classification and the library will generate the corresponding
CSFs internally.

Internally, hexpr uses the Distinct Row Table (DRT) representation and
Brooks-style loop traversal. The matrix element formulation follows
F. Sasaki's tensor recoupling and reduction technique, originally
developed for atomic systems in LS coupling and here applied to the
molecular case without assuming spatial symmetry.

hexpr is part of an ongoing research project that also produced
[EMOL](https://github.com/hh076/EMOL). Both projects are supported by the
same JSPS KAKENHI grant. The expression module developed here was sufficiently
mature and well-formulated to warrant separate release as a standalone library,
accompanied by a methodology paper.

## Repository layout

```text
hexpr/
+-- LICENSE
+-- README.md
+-- c/                  C library and test suite
|   +-- include/        Public header: hexpr.h
|   +-- src/            Library sources (common/, csf/, loopgen/, hexpr_*.c)
|   +-- wig/            Wigner 3j/6j/9j sub-library (F. Sasaki)
|   +-- tests/          C test suites
|   +-- Makefile
|   `-- README.md
+-- python/             Python ctypes binding
+-- fortran/            Fortran binding (iso_c_binding + F77 aliases)
+-- julia/              Julia ccall binding
+-- ruby/               Ruby ffi gem binding
+-- reference/          regression test files (C test harness; 38 of 42 in-tree)
+-- validation/         PySCF + OpenMolcas cross-validation scripts
+-- example_ci/         worked end-to-end CI examples (PySCF / Psi4 drivers)
`-- docs/               Design notes and index specifications
```

- **c/** -- Primary C library. Builds `libhexpr.so` and a standalone
  `run_reference` driver. All C test suites live in `c/tests/`.
  See [c/README.md](c/README.md).
- **python/** -- ctypes binding; requires NumPy and nothing else beyond
  the standard library. See [python/README.md](python/README.md).
- **fortran/** -- F90 module (`iso_c_binding`) and F77 underscore-suffixed
  aliases. See [fortran/README.md](fortran/README.md).
- **julia/** -- ccall binding via the `HExpr` module. See
  [julia/README.md](julia/README.md).
- **ruby/** -- ffi gem binding. See [ruby/README.md](ruby/README.md).
- **reference/** -- 42 regression test files total (14 cases x 3 file
  types: `.drt.txt`, `.expr.txt`, `.expr_one.txt`). The 4 largest
  (`*_med.expr.txt`) are distributed separately to keep clone size down,
  so 38 are committed in-tree. Used by the C test harness; see
  [c/README.md](c/README.md) for how to fetch the remaining 4.
- **validation/** -- PySCF + OpenMolcas cross-validation scripts. See
  [validation/README.md](validation/README.md).
- **example_ci/** -- Worked end-to-end examples that drive hexpr to CI
  energies, with PySCF- and Psi4-based integral drivers (`run_pyscf.py`
  / `run_psi4.py`) alongside Julia and Ruby versions. See
  [example_ci/README.md](example_ci/README.md).
- **docs/** -- Specifications and design notes, including the index
  encoding ([docs/PQRS_INDEX_SPEC.md](docs/PQRS_INDEX_SPEC.md)), the
  per-pair API contract
  ([docs/PAIR_API_SEMANTICS.md](docs/PAIR_API_SEMANTICS.md)), and how to
  compare hexpr with other CSF codes
  ([docs/COMPARING_WITH_OTHER_CSF_CODES.md](docs/COMPARING_WITH_OTHER_CSF_CODES.md)).

## Building and testing

For a complete install across all language bindings into a chosen
prefix (so that Python, Julia, Ruby, Fortran, and C all work from
one directory tree), see [INSTALL.md](INSTALL.md). The instructions
in this section are the quick reference for an in-tree build and
test cycle.

> macOS: substitute `libhexpr.dylib` for `libhexpr.so` and
> `DYLD_LIBRARY_PATH` for `LD_LIBRARY_PATH` throughout.

### C library

Requires: GCC (C11), GNU Make.

```sh
make -C c                  # build libhexpr.so and run_reference
make -C c/tests check      # run the core test suite
```

Expected output ends with `ALL TESTS PASSED`. See [c/README.md](c/README.md)
for full details.

### Bindings

All bindings load `libhexpr.so` built above. Each binding's README gives full
prerequisites and test instructions; quick reference:

- **Python** -- `cd python && PYTHONPATH=$PWD python3 -m pytest tests/ -v`
  -- see [python/README.md](python/README.md)
- **Fortran** -- `make -C fortran/tests check`
  -- see [fortran/README.md](fortran/README.md)
- **Julia** -- `julia --project=julia -e 'using Pkg; Pkg.test()'`
  -- see [julia/README.md](julia/README.md)
- **Ruby** -- `cd ruby && bundle exec ruby -Ilib -Itest -e 'Dir["test/test_*.rb"].each { |f| require_relative f }'`
  -- see [ruby/README.md](ruby/README.md)

## License

hexpr is released under the MIT License -- see [LICENSE](LICENSE).

The Wigner 3j/6j/9j coefficient library in `c/wig/` is separately copyrighted
by F. Sasaki under the same MIT license terms -- see
[c/wig/LICENSE](c/wig/LICENSE).

## Authors and acknowledgments

hexpr was developed by H. Honda and T. Noro.

This work originated from research supported by JSPS KAKENHI
Grant Number JP16K05660 ("Research and development of a
prototyping environment for high-performance electronic
structure calculations", PI: H. Honda, Kyushu University,
FY2016-2017), Grant-in-Aid for Scientific Research (C). The
Ruby-based reference implementation that hexpr was ported
from was developed under that project.

The Wigner 3j/6j/9j coefficient library bundled in `c/wig/` is an
original implementation written in Fortran by F. Sasaki, refactored
into modern C for use in hexpr. It is included with permission; its
license terms are given under "License" above.

## Development notes

Parts of this codebase were developed with AI assistance (Claude,
Anthropic): code review, refactoring, test design, source comments and
documentation. The numerical method, the original Ruby implementation
from which the C library was ported, and the reference data derived
from it are the authors' own. Every AI-assisted change was reviewed by
the authors before being committed, is covered by the test suite, and
is described in its commit message.

## Citation

If you use hexpr in published work, please cite both the methodology
paper and the software:

> H. Honda and T. Noro, "hexpr: A Configuration State Function Based
> Hamiltonian Matrix Element Expression Library", preprint (2026),
> Zenodo. doi:[10.5281/zenodo.22977329](https://doi.org/10.5281/zenodo.22977329)

> H. Honda and T. Noro, hexpr, software, Zenodo.
> doi:[10.5281/zenodo.22974832](https://doi.org/10.5281/zenodo.22974832)

The software DOI above covers all releases and resolves to the latest.
To cite the exact release you used, take its version DOI from the
Zenodo record; v1.0.0 is
[10.5281/zenodo.22974833](https://doi.org/10.5281/zenodo.22974833).
GitHub's "Cite this repository" button gives the paper reference in APA
and BibTeX form (from `CITATION.cff`).
