module HExpr

# In-tree shared-library extension: "so" on Linux, "dylib" on macOS.
# Base.Libc.Libdl.dlext avoids needing Libdl as a declared package dependency.
const _LIBHEXPR = get(ENV, "HEXPR_LIBRARY",
                      joinpath(@__DIR__, "..", "..", "c",
                               "libhexpr." * Base.Libc.Libdl.dlext))

# ---- Status codes ----
const HEXPR_OK              = Cint(0)
const HEXPR_ERR_INVALID_ARG = Cint(1)
const HEXPR_ERR_OOM         = Cint(2)
const HEXPR_ERR_IO          = Cint(3)
const HEXPR_ERR_RANGE       = Cint(4)
const HEXPR_ERR_NOT_BUILT   = Cint(5)
const HEXPR_ERR_INTERNAL    = Cint(99)

# ---- DRT kind ----
const HEXPR_DRT_SOCI   = Cint(0)
const HEXPR_DRT_FOCI   = Cint(1)
const HEXPR_DRT_CSFSET = Cint(2)

# ---- Expression selector ----
const HEXPR_EXPR_FULL = Cint(0)
const HEXPR_EXPR_ONE  = Cint(1)

# ---- Step codes ----
const HEXPR_STEP_E = Cint(0)
const HEXPR_STEP_U = Cint(1)
const HEXPR_STEP_D = Cint(2)
const HEXPR_STEP_F = Cint(3)

# ---- Lifecycle ----

function init()
    status = ccall((:hexpr_init, _LIBHEXPR), Cint, ())
    status == HEXPR_OK || error("hexpr_init failed: $status")
    nothing
end

function shutdown()
    ccall((:hexpr_shutdown, _LIBHEXPR), Cvoid, ())
    nothing
end

# ---- Version / error ----

function version_string()
    ptr = ccall((:hexpr_version_string, _LIBHEXPR), Cstring, ())
    unsafe_string(ptr)
end

version_major() = ccall((:hexpr_version_major, _LIBHEXPR), Cint, ())
version_minor() = ccall((:hexpr_version_minor, _LIBHEXPR), Cint, ())
version_patch() = ccall((:hexpr_version_patch, _LIBHEXPR), Cint, ())

function last_error()
    ptr = ccall((:hexpr_last_error, _LIBHEXPR), Cstring, ())
    ptr == C_NULL ? "" : unsafe_string(ptr)
end

# ---- C stdout FILE* for the *_show* wrappers ----
# The hexpr_*_show* C functions take a FILE * stream.  libc exposes the
# standard streams as real globals; cglobal(:stdout) resolves the symbol
# from the running process and we load the FILE* it holds.  Cached.
# libc names this symbol `stdout` on Linux but `__stdoutp` on macOS.
const _STDOUT_SYM = Sys.isapple() ? :__stdoutp : :stdout
const _C_STDOUT = Ref{Ptr{Cvoid}}(C_NULL)

function _c_stdout()
    if _C_STDOUT[] == C_NULL
        _C_STDOUT[] = unsafe_load(cglobal(_STDOUT_SYM, Ptr{Cvoid}))
    end
    _C_STDOUT[]
end

# ---- Drt struct ----

mutable struct Drt
    handle::Ptr{Cvoid}
    function Drt(kind::Integer, nve::Integer, spin_x2::Integer,
                 ncore::Integer, nval::Integer, next::Integer)
        h = ccall((:hexpr_drt_create, _LIBHEXPR), Ptr{Cvoid},
                  (Cint, Cint, Cint, Cint, Cint, Cint),
                  kind, nve, spin_x2, ncore, nval, next)
        h == C_NULL && error("hexpr_drt_create failed: $(last_error())")
        d = new(h)
        finalizer(d) do drt
            if drt.handle != C_NULL
                ccall((:hexpr_drt_destroy, _LIBHEXPR), Cvoid,
                      (Ptr{Cvoid},), drt.handle)
                drt.handle = C_NULL
            end
        end
        d
    end
    function Drt(csfs::AbstractMatrix{UInt8}, norb::Integer)
        ncsfs_count, ncols = size(csfs)
        ncols == norb || error("csfs second dim ($ncols) != norb ($norb)")
        buf = permutedims(csfs, (2, 1))  # (norb, ncsfs) contiguous copy for C row-major
        h = ccall((:hexpr_drt_from_csfset, _LIBHEXPR), Ptr{Cvoid},
                  (Cint, Ptr{UInt8}, Cint),
                  Cint(norb), buf, Cint(ncsfs_count))
        h == C_NULL && error("hexpr_drt_from_csfset failed: $(last_error())")
        d = new(h)
        finalizer(d) do drt
            if drt.handle != C_NULL
                ccall((:hexpr_drt_destroy, _LIBHEXPR), Cvoid,
                      (Ptr{Cvoid},), drt.handle)
                drt.handle = C_NULL
            end
        end
        d
    end
end

function Drt(kind::AbstractString; nve, spin_x2, ncore, nval, next)
    k = kind == "soci"   ? HEXPR_DRT_SOCI   :
        kind == "foci"   ? HEXPR_DRT_FOCI   :
        kind == "csfset" ? HEXPR_DRT_CSFSET :
        error("unknown DRT kind: $kind")
    Drt(k, nve, spin_x2, ncore, nval, next)
end

# ---- Drt accessors ----

norb(d::Drt)      = ccall((:hexpr_drt_norb,      _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
nnodes(d::Drt)    = ccall((:hexpr_drt_nnodes,    _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
ncsfs(d::Drt)     = ccall((:hexpr_drt_ncsfs,     _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
# Canonical GUGA-walk CSF count. CSFSET: may differ from ncsfs (superset
# M>=N, or dedup); SOCI/FOCI: equals ncsfs. See walk_csfs.
walk_ncsfs(d::Drt) = ccall((:hexpr_drt_walk_ncsfs, _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
kind(d::Drt)      = ccall((:hexpr_drt_kind,      _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
ncore(d::Drt)     = ccall((:hexpr_drt_ncore,     _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
nval(d::Drt)      = ccall((:hexpr_drt_nval,      _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
next_orb(d::Drt)  = ccall((:hexpr_drt_next,      _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
nve_state(d::Drt) = ccall((:hexpr_drt_nve_state, _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)
spin_x2(d::Drt)   = ccall((:hexpr_drt_spin_x2,   _LIBHEXPR), Cint, (Ptr{Cvoid},), d.handle)

function csfs(d::Drt)
    nc  = Int(ncsfs(d))
    no  = Int(norb(d))
    ptr = ccall((:hexpr_drt_csfs, _LIBHEXPR), Ptr{UInt8}, (Ptr{Cvoid},), d.handle)
    mat = unsafe_wrap(Array, ptr, (no, nc))  # (norb, ncsfs): column-major == C row-major
    permutedims(mat, (2, 1))                 # copy as (ncsfs, norb)
end

# Canonical GUGA-walk CSF list, shape (walk_ncsfs, norb). CSFSET: the internal
# canonical view, which may be a superset of csfs() (walk_ncsfs >= ncsfs) or a
# dedup; SOCI/FOCI: equals csfs().
function walk_csfs(d::Drt)
    nc  = Int(walk_ncsfs(d))
    no  = Int(norb(d))
    ptr = ccall((:hexpr_drt_walk_csfs, _LIBHEXPR), Ptr{UInt8}, (Ptr{Cvoid},), d.handle)
    mat = unsafe_wrap(Array, ptr, (no, nc))  # (norb, walk_ncsfs): column-major == C row-major
    permutedims(mat, (2, 1))                 # copy as (walk_ncsfs, norb)
end

# Node details by node index; returns a named tuple (orb, nve, is)
# matching the F90 hexpr_drt_node_info argument meaning.
function drt_node_info(d::Drt, node_index::Integer)
    out_orb = Ref{Cint}(0)
    out_nve = Ref{Cint}(0)
    out_is  = Ref{Cint}(0)
    status = ccall((:hexpr_drt_node_info, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}),
                   d.handle, Cint(node_index), out_orb, out_nve, out_is)
    status == HEXPR_OK || error("hexpr_drt_node_info failed: $(last_error())")
    (orb = Int(out_orb[]), nve = Int(out_nve[]), is = Int(out_is[]))
end

# Dump a human-readable DRT summary to stdout.
function drt_show(d::Drt)
    status = ccall((:hexpr_drt_show, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Cvoid}), d.handle, _c_stdout())
    status == HEXPR_OK || error("hexpr_drt_show failed: $(last_error())")
    nothing
end

# ---- Expression struct ----
# Named Expression (not Expr) to avoid collision with Base.Expr.

mutable struct Expression
    handle::Ptr{Cvoid}
    parent::Drt   # keep Drt alive while Expression exists
    function Expression(d::Drt)
        h = ccall((:hexpr_expr_build, _LIBHEXPR), Ptr{Cvoid},
                  (Ptr{Cvoid},), d.handle)
        h == C_NULL && error("hexpr_expr_build failed: $(last_error())")
        e = new(h, d)
        finalizer(e) do expr
            if expr.handle != C_NULL
                ccall((:hexpr_expr_destroy, _LIBHEXPR), Cvoid,
                      (Ptr{Cvoid},), expr.handle)
                expr.handle = C_NULL
            end
        end
        e
    end
end

# ---- Expression accessors ----
# nterms_full/nterms_one return int64_t -- use Int64, not Cint.

nterms_full(e::Expression) =
    ccall((:hexpr_expr_nterms_full, _LIBHEXPR), Int64, (Ptr{Cvoid},), e.handle)

nterms_one(e::Expression) =
    ccall((:hexpr_expr_nterms_one, _LIBHEXPR), Int64, (Ptr{Cvoid},), e.handle)

# ---- Bulk read ----

function read_full(e::Expression)
    n = nterms_full(e)
    ij   = Vector{Int32}(undef, n)
    pqrs = Vector{Int32}(undef, n)
    coef = Vector{Float64}(undef, n)
    status = ccall((:hexpr_expr_read, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
                   e.handle, HEXPR_EXPR_FULL, ij, pqrs, coef)
    status == HEXPR_OK || error("hexpr_expr_read failed: $(last_error())")
    ij, pqrs, coef
end

function read_one(e::Expression)
    n = nterms_one(e)
    ij   = Vector{Int32}(undef, n)
    pqrs = Vector{Int32}(undef, n)
    coef = Vector{Float64}(undef, n)
    status = ccall((:hexpr_expr_read, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
                   e.handle, HEXPR_EXPR_ONE, ij, pqrs, coef)
    status == HEXPR_OK || error("hexpr_expr_read failed: $(last_error())")
    ij, pqrs, coef
end

# ---- Expression show ----
# Dump the full / 1-electron expression to stdout.

function expr_show_full(e::Expression)
    status = ccall((:hexpr_expr_show_full, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Cvoid}), e.handle, _c_stdout())
    status == HEXPR_OK || error("hexpr_expr_show_full failed: $(last_error())")
    nothing
end

function expr_show_one(e::Expression)
    status = ccall((:hexpr_expr_show_one, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Cvoid}), e.handle, _c_stdout())
    status == HEXPR_OK || error("hexpr_expr_show_one failed: $(last_error())")
    nothing
end

# ---- Per-pair / block / subspace evaluation ----

function _which_to_cint(which::Symbol)
    which === :full && return HEXPR_EXPR_FULL
    which === :one  && return HEXPR_EXPR_ONE
    throw(ArgumentError("which must be :full or :one, got $(repr(which))"))
end

function csf_index(drt::Drt, steps::AbstractVector{UInt8})
    length(steps) == Int(norb(drt)) ||
        error("steps length $(length(steps)) != norb $(norb(drt))")
    buf = Vector{UInt8}(steps)
    out_idx = Ref{Cint}(0)
    status = ccall((:hexpr_csf_index, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{UInt8}, Ref{Cint}),
                   drt.handle, buf, out_idx)
    status == HEXPR_OK || error("hexpr_csf_index failed: $(last_error())")
    Int(out_idx[])
end

# Step codes for CSF row `index`; returns a Vector{UInt8} of length norb(drt).
function csf_steps(drt::Drt, index::Integer)
    no  = Int(norb(drt))
    out = Vector{UInt8}(undef, no)
    status = ccall((:hexpr_csf_steps, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Ptr{UInt8}),
                   drt.handle, Cint(index), out)
    status == HEXPR_OK || error("hexpr_csf_steps failed: $(last_error())")
    out
end

# Parse an "eudf"-style string of length norb into step codes.
function step_parse(s::AbstractString, norb::Integer)
    out = Vector{UInt8}(undef, norb)
    status = ccall((:hexpr_step_parse, _LIBHEXPR), Cint,
                   (Cstring, Cint, Ptr{UInt8}),
                   s, Cint(norb), out)
    status == HEXPR_OK || error("hexpr_step_parse failed: $(last_error())")
    out
end

# ---- CSF -> determinant expansion ----

# Number of determinants in the expansion of CSF `index` at `ms_x2`.
# Exact, not an upper bound: it is C(n_open, n_alpha_open).
#
# `index` follows the same convention as csf_steps -- the input position,
# in 0:(ncsfs-1), 0-based as every index in this binding is.  So do
# pair_eval/block_eval/subspace_eval: there is one index space across the
# whole API.
#
# `ms_x2` is 2*Ms.  Any |Ms| <= S with the same parity as spin_x2 is
# accepted, negative values included; Ms is not assumed to equal S.
function csf_ndets(drt::Drt, index::Integer, ms_x2::Integer)
    out = Ref{Cint}(0)
    status = ccall((:hexpr_csf_ndets, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Ref{Cint}),
                   drt.handle, Cint(index), Cint(ms_x2), out)
    status == HEXPR_OK || error("hexpr_csf_ndets failed: $(last_error())")
    Int(out[])
end

# Expand CSF `index` into Slater determinants at `ms_x2`.
#
# Returns (alpha, beta, coef).  alpha and beta are Matrix{UInt64} of shape
# (ndets, nwords) with nwords = (norb + 63) >> 6; coef is a Vector{Float64}
# of length ndets.  Bit p of word w is spatial orbital 64*(w-1) + p, so
# orbital p (0-based) is occupied by an alpha electron in determinant k iff
#
#     (alpha[k, (p >> 6) + 1] >> (p & 63)) & 1 == 1
#
# The second axis is kept even when nwords == 1, because collapsing it
# would change the caller's indexing at norb > 64.
#
# Conventions, all fixed by the C API and derived in c/src/csf/csf2det.c:
#
#   * Spin coupling: genealogical, in the Eq. (5) order -- the new
#     one-orbital state couples on the LEFT of the accumulated partial
#     CSF.  A CSF in another coupling scheme will not reproduce these
#     coefficients.
#   * Clebsch-Gordan phase: Condon-Shortley.
#   * Operator order: all alpha in ascending spatial orbital, then all
#     beta in ascending spatial orbital.  An interleaved-by-orbital
#     convention (some qubit mappings) differs by a determinant-dependent
#     permutation sign.
#   * Determinant order: lexicographically ascending in the set of
#     open-shell positions carrying alpha.
#   * Normalization: sum(coef .^ 2) == 1 for every CSF at every Ms.
#   * Zero coefficients are emitted, not dropped, so the count is exact
#     and no threshold enters this API.
#
# The OVERALL sign of a CSF is a convention, not a measured quantity: it
# is whatever the coupling order and the operator order produce.  Relative
# signs -- within one CSF, and between CSFs of one DRT -- are meaningful
# and are verified against PySCF by validation/csf2det_vs_pyscf.py.  The
# overall sign is not something to assert against.
function csf_to_dets(drt::Drt, index::Integer, ms_x2::Integer)
    ndets  = csf_ndets(drt, index, ms_x2)
    nwords = (Int(norb(drt)) + 63) >> 6
    # C writes out_alpha[k*nwords + w].  A Julia (nwords, ndets) matrix is
    # column-major, so element (w+1, k+1) sits at exactly that offset --
    # the same reinterpretation csfs() relies on.  Transposed on the way
    # out so that one row is one determinant, matching csfs() returning
    # one row per CSF.
    a    = Matrix{UInt64}(undef, nwords, ndets)
    b    = Matrix{UInt64}(undef, nwords, ndets)
    coef = Vector{Float64}(undef, ndets)
    out_count = Ref{Cint}(0)
    status = ccall((:hexpr_csf_to_dets, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Ref{Cint},
                    Ptr{UInt64}, Ptr{UInt64}, Ptr{Float64}),
                   drt.handle, Cint(index), Cint(ms_x2), out_count,
                   a, b, coef)
    status == HEXPR_OK || error("hexpr_csf_to_dets failed: $(last_error())")
    n = Int(out_count[])
    n == ndets || error("hexpr_csf_to_dets wrote $n determinants but " *
                        "csf_ndets promised $ndets")
    permutedims(a, (2, 1)), permutedims(b, (2, 1)), coef
end

function pair_count(drt::Drt, bra::Integer, ket::Integer; which::Symbol=:full)
    out = Ref{Cint}(0)
    status = ccall((:hexpr_pair_count, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Cint, Ref{Cint}),
                   drt.handle, Cint(bra), Cint(ket), _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_pair_count failed: $(last_error())")
    Int(out[])
end

function pair_max_terms(drt::Drt; which::Symbol=:full)
    out = Ref{Cint}(0)
    status = ccall((:hexpr_pair_max_terms, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Ref{Cint}),
                   drt.handle, _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_pair_max_terms failed: $(last_error())")
    Int(out[])
end

# `bra` and `ket` are input-order CSF indices in 0:(ncsfs-1), the same space
# as csf_steps and csf_index.  The returned `ij` is a packed CANONICAL pair
# address, not a pair of indices into this DRT (docs/PQRS_INDEX_SPEC.md);
# do not decode it and feed the result back.
function pair_eval(drt::Drt, bra::Integer, ket::Integer; which::Symbol=:full)
    which_int = _which_to_cint(which)
    max_t = pair_max_terms(drt; which=which)
    ij   = Vector{Int32}(undef, max_t)
    pqrs = Vector{Int32}(undef, max_t)
    coef = Vector{Float64}(undef, max_t)
    out_count = Ref{Cint}(0)
    status = ccall((:hexpr_pair_eval, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Cint, Ref{Cint},
                    Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
                   drt.handle, Cint(bra), Cint(ket), which_int,
                   out_count, ij, pqrs, coef)
    status == HEXPR_OK || error("hexpr_pair_eval failed: $(last_error())")
    n = Int(out_count[])
    ij[1:n], pqrs[1:n], coef[1:n]
end

function block_count(drt::Drt, bra_indices, ket_indices; which::Symbol=:full)
    bra = Int32.(collect(bra_indices))
    ket = Int32.(collect(ket_indices))
    out = Ref{Cint}(0)
    status = ccall((:hexpr_block_count, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Int32}, Cint, Ptr{Int32}, Cint, Cint, Ref{Cint}),
                   drt.handle, bra, Cint(length(bra)), ket, Cint(length(ket)),
                   _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_block_count failed: $(last_error())")
    Int(out[])
end

function block_max_terms(drt::Drt, n_bra::Integer, n_ket::Integer; which::Symbol=:full)
    out = Ref{Cint}(0)
    status = ccall((:hexpr_block_max_terms, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Cint, Ref{Cint}),
                   drt.handle, Cint(n_bra), Cint(n_ket), _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_block_max_terms failed: $(last_error())")
    Int(out[])
end

# `bra_indices`/`ket_indices` hold input-order CSF indices in 0:(ncsfs-1).
# The returned bra_pos/ket_pos are positions within those arrays, and they,
# not `ij`, are how a term is attributed to a pair.
function block_eval(drt::Drt, bra_indices, ket_indices; which::Symbol=:full)
    which_int = _which_to_cint(which)
    bra = Int32.(collect(bra_indices))
    ket = Int32.(collect(ket_indices))
    n_bra, n_ket = length(bra), length(ket)
    max_t = block_max_terms(drt, n_bra, n_ket; which=which)
    bra_pos = Vector{Int32}(undef, max_t)
    ket_pos = Vector{Int32}(undef, max_t)
    ij      = Vector{Int32}(undef, max_t)
    pqrs    = Vector{Int32}(undef, max_t)
    coef    = Vector{Float64}(undef, max_t)
    out_count = Ref{Cint}(0)
    status = ccall((:hexpr_block_eval, _LIBHEXPR), Cint,
                   (Ptr{Cvoid},
                    Ptr{Int32}, Cint, Ptr{Int32}, Cint,
                    Cint, Ref{Cint},
                    Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
                   drt.handle,
                   bra, Cint(n_bra), ket, Cint(n_ket),
                   which_int, out_count,
                   bra_pos, ket_pos, ij, pqrs, coef)
    status == HEXPR_OK || error("hexpr_block_eval failed: $(last_error())")
    n = Int(out_count[])
    ij[1:n], pqrs[1:n], coef[1:n], bra_pos[1:n], ket_pos[1:n]
end

function subspace_count(drt::Drt, indices; which::Symbol=:full)
    idx = Int32.(collect(indices))
    out = Ref{Cint}(0)
    status = ccall((:hexpr_subspace_count, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Int32}, Cint, Cint, Ref{Cint}),
                   drt.handle, idx, Cint(length(idx)), _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_subspace_count failed: $(last_error())")
    Int(out[])
end

function subspace_max_terms(drt::Drt, n::Integer; which::Symbol=:full)
    out = Ref{Cint}(0)
    status = ccall((:hexpr_subspace_max_terms, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Cint, Cint, Ref{Cint}),
                   drt.handle, Cint(n), _which_to_cint(which), out)
    status == HEXPR_OK || error("hexpr_subspace_max_terms failed: $(last_error())")
    Int(out[])
end

# `indices` holds input-order CSF indices in 0:(ncsfs-1); bra_pos/ket_pos come
# back as positions within that array.
function subspace_eval(drt::Drt, indices; which::Symbol=:full)
    which_int = _which_to_cint(which)
    idx = Int32.(collect(indices))
    n = length(idx)
    max_t = subspace_max_terms(drt, n; which=which)
    bra_pos = Vector{Int32}(undef, max_t)
    ket_pos = Vector{Int32}(undef, max_t)
    ij      = Vector{Int32}(undef, max_t)
    pqrs    = Vector{Int32}(undef, max_t)
    coef    = Vector{Float64}(undef, max_t)
    out_count = Ref{Cint}(0)
    status = ccall((:hexpr_subspace_eval, _LIBHEXPR), Cint,
                   (Ptr{Cvoid}, Ptr{Int32}, Cint,
                    Cint, Ref{Cint},
                    Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
                   drt.handle, idx, Cint(n),
                   which_int, out_count,
                   bra_pos, ket_pos, ij, pqrs, coef)
    status == HEXPR_OK || error("hexpr_subspace_eval failed: $(last_error())")
    n_out = Int(out_count[])
    ij[1:n_out], pqrs[1:n_out], coef[1:n_out], bra_pos[1:n_out], ket_pos[1:n_out]
end

# ---- High-level CSF step-code functions ----

# The three functions below build a sub-DRT from the supplied rows and then
# address those rows by their INPUT positions, 0..n-1. Earlier versions kept a
# {row -> canonical index} lookup here because the eval API took canonical
# (walk) indices, which differ from input order for a from_csfset sub-DRT.
# The library does that translation now, so the lookup is gone.

function csf_pair_eval(bra_steps::AbstractVector{UInt8},
                       ket_steps::AbstractVector{UInt8};
                       which::Symbol=:full)
    norb_val = length(bra_steps)
    length(ket_steps) == norb_val ||
        error("bra_steps length $norb_val != ket_steps length $(length(ket_steps))")
    stacked = Matrix{UInt8}(undef, 2, norb_val)
    stacked[1, :] = bra_steps
    stacked[2, :] = ket_steps
    sub = Drt(stacked, norb_val)
    # Row 0 is the bra, row 1 is the ket, whatever the canonical order is.
    pair_eval(sub, Int32(0), Int32(1); which=which)
end

function csf_block_eval(bra_csfs::AbstractMatrix{UInt8},
                        ket_csfs::AbstractMatrix{UInt8};
                        which::Symbol=:full)
    n_bra, norb_b = size(bra_csfs)
    n_ket, norb_k = size(ket_csfs)
    norb_b == norb_k ||
        error("bra_csfs and ket_csfs have different norb: $norb_b vs $norb_k")
    norb_val = norb_b
    all_csfs = vcat(bra_csfs, ket_csfs)
    sub = Drt(all_csfs, norb_val)
    # The sub-DRT's input rows are bra_csfs followed by ket_csfs, so the two
    # index sets are the two halves of 0:(n_bra + n_ket - 1).
    bra_indices = Int32[i for i in 0:(n_bra - 1)]
    ket_indices = Int32[j for j in n_bra:(n_bra + n_ket - 1)]
    block_eval(sub, bra_indices, ket_indices; which=which)
end

function csf_subspace_eval(csfs::AbstractMatrix{UInt8}; which::Symbol=:full)
    n, norb_val = size(csfs)
    sub = Drt(csfs, norb_val)
    indices = Int32[i for i in 0:(n - 1)]
    subspace_eval(sub, indices; which=which)
end

# ---- Exports ----

export Drt, Expression
export init, shutdown, version_string
export version_major, version_minor, version_patch, last_error
export norb, nnodes, ncsfs, kind, ncore, nval, next_orb, nve_state, spin_x2, csfs
export walk_csfs, walk_ncsfs
export drt_node_info, drt_show
export nterms_full, nterms_one, read_full, read_one, expr_show_full, expr_show_one
export csf_index, csf_steps, step_parse
export csf_ndets, csf_to_dets
export pair_count, pair_max_terms, pair_eval
export block_count, block_max_terms, block_eval
export subspace_count, subspace_max_terms, subspace_eval
export csf_pair_eval, csf_block_eval, csf_subspace_eval
export HEXPR_OK, HEXPR_DRT_SOCI, HEXPR_DRT_FOCI, HEXPR_DRT_CSFSET
export HEXPR_EXPR_FULL, HEXPR_EXPR_ONE
export HEXPR_STEP_E, HEXPR_STEP_U, HEXPR_STEP_D, HEXPR_STEP_F
export HEXPR_ERR_INVALID_ARG, HEXPR_ERR_OOM, HEXPR_ERR_IO,
       HEXPR_ERR_RANGE, HEXPR_ERR_NOT_BUILT, HEXPR_ERR_INTERNAL

end # module HExpr
