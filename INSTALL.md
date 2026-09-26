# Installing hexpr

This document describes a single end-to-end install path that produces a
working setup for all language bindings of hexpr from one chosen prefix,
and is intended for users who want the full install. If you only need
one binding, read that binding's README directly:
[C](c/README.md), [Python](python/README.md), [Julia](julia/README.md),
[Ruby](ruby/README.md), or [Fortran](fortran/README.md). The install
model is a plain copy to a prefix plus a handful of environment
variables; no package manager is required.

## 1. Prerequisites

**Primary supported environment**

Ubuntu 22.04 LTS and Ubuntu 24.04 LTS. macOS (Intel, Apple-clang) is also
supported: the bindings load `libhexpr.dylib`, and the Fortran tests run
without any environment variable via a baked rpath. Other Linux
distributions and WSL are expected to work but are untested. Please report
issues if anything breaks on your platform.

> **macOS note.** Throughout this document, wherever a Linux example shows
> `libhexpr.so` or `LD_LIBRARY_PATH`, the macOS equivalents are
> `libhexpr.dylib` and `DYLD_LIBRARY_PATH`. The install script and bindings
> select the right names automatically.

**Already assumed to be present**

These are taken as given on any reasonably set-up Linux development
machine, and this document does not provide install instructions for
them:

- a C11 compiler (GCC 9+ or Clang 11+)
- GNU Make
- `git`
- the standard build essentials (`build-essential` on Debian/Ubuntu)

If you are on a fresh Ubuntu install and lack these, the following
one-liner is enough:

```sh
sudo apt install build-essential git
```

**Additional dependencies for the components you plan to use**

Only what hexpr specifically requires on top of the baseline above:

- Python binding: `sudo apt install python3 python3-numpy` (Python 3.9
  or later).
- Julia binding: install the official Julia >= 1.12 binary from
  https://julialang.org. Distro packages are usually older.
- Ruby binding: `sudo apt install ruby ruby-dev`, then
  `gem install ffi`. Ruby 3.1 or later.
- Fortran binding: `sudo apt install gfortran` (9 or later). LLVM
  Flang 20+ and Intel ifx are also supported; see
  [fortran/README.md](fortran/README.md) for details.
- `example_ci` Python examples (optional): `pip install pyscf h5py`.
- `validation/` end-to-end run (optional): OpenMolcas. Most users do
  not need this.

## 2. Build the C library

1. Clone the repository:

   ```sh
   git clone https://github.com/emol-project/hexpr.git hexpr
   cd hexpr
   ```

2. Build the shared library:

   ```sh
   make -C c
   ```

   This produces `c/libhexpr.so` and the `c/run_reference` CLI driver.

3. (Optional sanity check) Run the C test suite:

   ```sh
   make -C c/tests check
   ```

   See [c/README.md](c/README.md) for details on the test suite.

## 3. Install to a prefix

Pick any writable directory as the prefix. This document uses
`$HOME/local` as the example:

```sh
bash scripts/install.sh $HOME/local
```

The script does no compilation; it is a plain copy plus directory
creation. It copies `libhexpr.so` to `$PREFIX/lib/`, the public
header(s) to `$PREFIX/include/`, and each language binding to
`$PREFIX/share/hexpr/<lang>/` (Python, Julia, Ruby, Fortran). After the
copy completes, the install no longer references any file inside the
source repository. The script prints the environment variable block at
the end; the next section explains those.

## 4. Set environment variables

Add the following to `~/.bashrc` (or your shell's rc file) for
persistence:

```sh
export HEXPR_PREFIX=$HOME/local
export HEXPR_LIBRARY=$HEXPR_PREFIX/lib/libhexpr.so
export LD_LIBRARY_PATH=$HEXPR_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HEXPR_PREFIX/share/hexpr/python:$PYTHONPATH
export JULIA_LOAD_PATH=$HEXPR_PREFIX/share/hexpr/julia:$JULIA_LOAD_PATH
export RUBYLIB=$HEXPR_PREFIX/share/hexpr/ruby/lib:$RUBYLIB
```

| Variable | Used by | Purpose |
|---|---|---|
| `HEXPR_PREFIX` | the user, as a convenience | shorthand for the install root |
| `HEXPR_LIBRARY` | Python, Julia, Ruby bindings | full path to `libhexpr.so` |
| `LD_LIBRARY_PATH` | the dynamic linker | finds `libhexpr.so` for C programs |
| `PYTHONPATH` | Python | finds the `hexpr` package |
| `JULIA_LOAD_PATH` | Julia | finds the `HExpr` package |
| `RUBYLIB` | Ruby | finds `hexpr.rb` |

`HEXPR_LIBRARY` is the authoritative source for the bindings: each
binding prefers it over any other discovery path. `LD_LIBRARY_PATH` is
for C and Fortran programs that link against `libhexpr.so` directly.

The Fortran binding is distributed as source
(`$HEXPR_PREFIX/share/hexpr/fortran/hexpr_mod.f90`); compile it together
with your own Fortran project. See [fortran/README.md](fortran/README.md).

## 5. Verify the installation

**Quick check.** Run the included verification script from the repo:

```sh
bash scripts/install_test.sh
```

This creates a scratch `install_test_<timestamp>/` directory, installs
into it, runs a smoke test for C, Python, Julia, and Ruby, and reports
pass/fail. Expected output ends with `PASS: 4/4`. Delete the scratch
directory when done. Use this to verify the installer works on your
system.

**Live check against your real prefix.** With the environment variables
from section 4 set, run these four commands:

```sh
# C: build and run a one-line smoke test
echo '#include <stdio.h>
#include "hexpr.h"
int main(void){hexpr_init();printf("%s\n",hexpr_version_string());return 0;}' > /tmp/smoke.c
gcc -I$HEXPR_PREFIX/include /tmp/smoke.c -L$HEXPR_PREFIX/lib -lhexpr -o /tmp/smoke
/tmp/smoke

# Python
python3 -c 'import hexpr; print(hexpr.version())'

# Julia
julia --startup-file=no -e 'using HExpr; HExpr.init(); println(HExpr.version_string())'

# Ruby
ruby -e 'require "hexpr"; HExpr.init; puts HExpr.version_string'
```

Each command should print the hexpr version string.

Note: four large reference files (`reference/*_med.expr.txt`) are
not bundled with the repository. The test suite still reports
`ALL TESTS PASSED` without them; the affected comparisons print
`SKIP`. To restore full regression coverage, download
[reference.tar.gz](https://github.com/emol-project/hexpr/releases/download/v1.0.0/reference.tar.gz)
and extract it on top of the repository. See [c/README.md](c/README.md)
for details.

## 6. Optional components

### 6.1 example_ci

The [example_ci/](example_ci/README.md) directory ships worked
end-to-end examples (H2 and H2O SOCI) in Python, Julia, and Ruby.
The examples are installed alongside the bindings at
`$HEXPR_PREFIX/share/hexpr/example_ci/`, so they can be run from
either the installed tree or the source repository.

To run them, install the prerequisites listed in section 1 (`pyscf`,
`h5py`), then see [example_ci/README.md](example_ci/README.md). The
Julia examples use their own project file
(`example_ci/Project.toml`) so that the `HExpr` binding itself does
not depend on HDF5.

### 6.2 validation

The [validation/](validation/README.md) directory contains an
OpenMolcas-anchored end-to-end validation of the PySCF + hexpr pipeline
for H2O / 6-31G(d) SOCI. It is not required to use hexpr. See
[validation/README.md](validation/README.md) to reproduce the reference
numbers.

## 7. Daily use

### Starting a REPL

With the environment variables set:

```sh
python3                # then `import hexpr; hexpr.init()`
julia                  # then `using HExpr; HExpr.init()`
irb                    # then `require "hexpr"; HExpr.init`
```

For C and Fortran, see the `Quick start` sections in their respective
READMEs.

### Updating the install after a rebuild

After editing source and rerunning `make -C c`:

```sh
bash scripts/install.sh $HOME/local
```

The script overwrites whatever is already at the prefix. No uninstall
step is required first.

### Note on Julia precompile caching

Julia evaluates `const _LIBHEXPR = get(ENV, "HEXPR_LIBRARY", ...)` at
module precompile time, so the value of `HEXPR_LIBRARY` is baked into
the precompiled `HExpr` module image. If you later change
`HEXPR_LIBRARY` to point to a different `libhexpr.so`, Julia may
still use the cached value. Force recompilation with:

```sh
julia --startup-file=no -e 'using Pkg; Pkg.precompile()'
```

Or delete `~/.julia/compiled/v*/HExpr/` to be sure.

## 8. Troubleshooting

1. **`cannot open shared object file: libhexpr.so`**
   `LD_LIBRARY_PATH` is missing `$HEXPR_PREFIX/lib`. Re-export the
   variables from section 4, or restart the shell after editing
   `~/.bashrc`.

2. **Python: `ImportError: cannot import name ...` or empty namespace package**
   `PYTHONPATH` likely contains `.` or `..`, which causes Python's
   namespace-package lookup to shadow the real hexpr install. See the
   `PYTHONPATH caution` section in [python/README.md](python/README.md)
   for the fix.

3. **Julia: changes to `HEXPR_LIBRARY` are ignored**
   Julia's precompile cache holds the value of `HEXPR_LIBRARY` that
   was in effect when the `HExpr` module was first compiled. Run
   `julia -e 'using Pkg; Pkg.precompile()'` or delete the
   `~/.julia/compiled/v*/HExpr/` cache to pick up a new value.

4. **Ruby: `cannot load such file -- ffi (LoadError)`**
   The `ffi` gem is not installed. Run `gem install ffi`. If you use
   `rbenv` or similar, ensure you install into the active Ruby.

5. **Fortran: `.mod` file mismatch when linking**
   `hexpr_mod.mod` is compiler-specific. Rebuild it from
   `$HEXPR_PREFIX/share/hexpr/fortran/hexpr_mod.f90` with the same
   Fortran compiler used for the rest of your project. See
   [fortran/README.md](fortran/README.md).

## 9. Uninstall

Run the uninstall script:

```sh
bash scripts/uninstall.sh $HEXPR_PREFIX
```

This reads `$HEXPR_PREFIX/share/hexpr/install_manifest.txt` and
removes every file listed there. Empty directories under
`$HEXPR_PREFIX/share/hexpr/` are also cleaned up; `$HEXPR_PREFIX/lib/`
and `$HEXPR_PREFIX/include/` are left in place in case other
packages share them.

Then remove the export lines from `~/.bashrc` (or wherever you put them).

## See also

- Project overview: [README.md](README.md)
- C library: [c/README.md](c/README.md)
- Public API header: [c/include/hexpr.h](c/include/hexpr.h)
- Python binding: [python/README.md](python/README.md)
- Julia binding: [julia/README.md](julia/README.md)
- Ruby binding: [ruby/README.md](ruby/README.md)
- Fortran binding: [fortran/README.md](fortran/README.md)
- End-to-end examples: [example_ci/README.md](example_ci/README.md)
- Validation: [validation/README.md](validation/README.md)
