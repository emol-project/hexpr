require "ffi"
require "pathname"

module HExpr
  # In-tree shared-library name: libhexpr.so on Linux, libhexpr.dylib on macOS.
  LIB_PATH = ENV["HEXPR_LIBRARY"] ||
             Pathname.new(__dir__).join("..", "..", "c",
               "libhexpr.#{FFI::Platform::LIBSUFFIX}").expand_path.to_s

  extend FFI::Library
  ffi_lib LIB_PATH

  # Status codes
  HEXPR_OK              = 0
  HEXPR_ERR_INVALID_ARG = 1
  HEXPR_ERR_OOM         = 2
  HEXPR_ERR_IO          = 3
  HEXPR_ERR_RANGE       = 4
  HEXPR_ERR_NOT_BUILT   = 5
  HEXPR_ERR_INTERNAL    = 99

  # DRT kind
  HEXPR_DRT_SOCI   = 0
  HEXPR_DRT_FOCI   = 1
  HEXPR_DRT_CSFSET = 2

  # Expression selector
  HEXPR_EXPR_FULL = 0
  HEXPR_EXPR_ONE  = 1

  # Step codes
  HEXPR_STEP_E = 0
  HEXPR_STEP_U = 1
  HEXPR_STEP_D = 2
  HEXPR_STEP_F = 3

  # Lifecycle and version
  attach_function :hexpr_init,           [], :int
  attach_function :hexpr_shutdown,       [], :void
  attach_function :hexpr_version_string, [], :string
  attach_function :hexpr_version_major,  [], :int
  attach_function :hexpr_version_minor,  [], :int
  attach_function :hexpr_version_patch,  [], :int
  attach_function :hexpr_last_error,     [], :string

  # DRT
  attach_function :hexpr_drt_create,
                  [:int, :int, :int, :int, :int, :int], :pointer
  attach_function :hexpr_drt_from_csfset,
                  [:int, :pointer, :int], :pointer
  attach_function :hexpr_drt_destroy,    [:pointer], :void
  attach_function :hexpr_drt_norb,       [:pointer], :int
  attach_function :hexpr_drt_nnodes,     [:pointer], :int
  attach_function :hexpr_drt_ncsfs,      [:pointer], :int
  attach_function :hexpr_drt_walk_ncsfs, [:pointer], :int
  attach_function :hexpr_drt_kind,       [:pointer], :int
  attach_function :hexpr_drt_ncore,      [:pointer], :int
  attach_function :hexpr_drt_nval,       [:pointer], :int
  attach_function :hexpr_drt_next,       [:pointer], :int
  attach_function :hexpr_drt_nve_state,  [:pointer], :int
  attach_function :hexpr_drt_spin_x2,    [:pointer], :int
  attach_function :hexpr_drt_csfs,       [:pointer], :pointer
  attach_function :hexpr_drt_walk_csfs,  [:pointer], :pointer

  # Expression
  attach_function :hexpr_expr_build,       [:pointer], :pointer
  attach_function :hexpr_expr_destroy,     [:pointer], :void
  attach_function :hexpr_expr_nterms_full, [:pointer], :int64
  attach_function :hexpr_expr_nterms_one,  [:pointer], :int64
  attach_function :hexpr_expr_read,
                  [:pointer, :int, :pointer, :pointer, :pointer], :int

  # CSF index / steps / parse
  attach_function :hexpr_csf_index,
                  [:pointer, :pointer, :pointer], :int
  attach_function :hexpr_csf_steps,
                  [:pointer, :int, :pointer], :int
  attach_function :hexpr_step_parse,
                  [:string, :int, :pointer], :int

  # CSF -> determinant expansion
  attach_function :hexpr_csf_ndets,
                  [:pointer, :int, :int, :pointer], :int
  attach_function :hexpr_csf_to_dets,
                  [:pointer, :int, :int, :pointer,
                   :pointer, :pointer, :pointer], :int

  # DRT node info
  attach_function :hexpr_drt_node_info,
                  [:pointer, :int, :pointer, :pointer, :pointer], :int

  # Show / dump (to a C FILE* stream)
  attach_function :hexpr_drt_show,        [:pointer, :pointer], :int
  attach_function :hexpr_expr_show_full,  [:pointer, :pointer], :int
  attach_function :hexpr_expr_show_one,   [:pointer, :pointer], :int

  # C stdout FILE* for the *_show* wrappers (libc exposes `stdout` as a
  # real global holding a FILE*; the show wrappers always emit there).
  # That symbol is named `stdout` on Linux but `__stdoutp` on macOS.
  _stdout_sym = FFI::Platform.mac? ? "__stdoutp" : "stdout"
  attach_variable :c_stdout, _stdout_sym, :pointer

  attach_function :hexpr_pair_count,
                  [:pointer, :int, :int, :int, :pointer], :int
  attach_function :hexpr_pair_eval,
                  [:pointer, :int, :int, :int, :pointer,
                   :pointer, :pointer, :pointer], :int
  attach_function :hexpr_pair_max_terms,
                  [:pointer, :int, :pointer], :int

  attach_function :hexpr_block_count,
                  [:pointer, :pointer, :int, :pointer, :int, :int, :pointer], :int
  attach_function :hexpr_block_eval,
                  [:pointer,
                   :pointer, :int, :pointer, :int,
                   :int, :pointer,
                   :pointer, :pointer, :pointer, :pointer, :pointer], :int
  attach_function :hexpr_block_max_terms,
                  [:pointer, :int, :int, :int, :pointer], :int

  attach_function :hexpr_subspace_count,
                  [:pointer, :pointer, :int, :int, :pointer], :int
  attach_function :hexpr_subspace_eval,
                  [:pointer, :pointer, :int,
                   :int, :pointer,
                   :pointer, :pointer, :pointer, :pointer, :pointer], :int
  attach_function :hexpr_subspace_max_terms,
                  [:pointer, :int, :int, :pointer], :int

  module_function

  def init
    status = hexpr_init
    raise "hexpr_init failed: #{status}" unless status == HEXPR_OK
    nil
  end

  def shutdown
    hexpr_shutdown
    nil
  end

  def version_string = hexpr_version_string
  def version_major  = hexpr_version_major
  def version_minor  = hexpr_version_minor
  def version_patch  = hexpr_version_patch
  def last_error     = hexpr_last_error || ""

  # Parse an "eudf"-style step string of length norb into an
  # Array<Integer> of step codes (0..3, E/U/D/F). For example,
  # step_parse("fude", 4) => [3, 1, 2, 0].
  def step_parse(str, norb)
    buf = FFI::MemoryPointer.new(:uint8, norb)
    status = hexpr_step_parse(str, norb, buf)
    raise "hexpr_step_parse failed: #{last_error}" unless status == HEXPR_OK
    buf.read_array_of_uint8(norb)
  end

  # The three methods below build a sub-DRT from the supplied rows and then
  # address those rows by their INPUT positions, 0..n-1. Earlier versions kept
  # a {row => canonical index} lookup here because the eval API took canonical
  # (walk) indices, which differ from input order for a from_csfset sub-DRT.
  # The library does that translation now, so the lookup is gone.

  def csf_pair_eval(bra_steps, ket_steps, which: :full)
    norb = bra_steps.size
    unless ket_steps.size == norb
      raise ArgumentError, "ket_steps length #{ket_steps.size} != norb #{norb}"
    end
    sub = Drt.from_csfset([bra_steps.to_a, ket_steps.to_a], norb)
    # Row 0 is the bra, row 1 is the ket, whatever the canonical order is.
    sub.pair_eval(0, 1, which: which)
  end

  def csf_block_eval(bra_csfs, ket_csfs, which: :full)
    raise ArgumentError, "bra_csfs must not be empty" if bra_csfs.empty?
    raise ArgumentError, "ket_csfs must not be empty" if ket_csfs.empty?
    norb = bra_csfs[0].size
    unless bra_csfs.all? { |r| r.size == norb } && ket_csfs.all? { |r| r.size == norb }
      raise ArgumentError, "all CSF rows must have the same length (norb=#{norb})"
    end
    sub = Drt.from_csfset(bra_csfs.map(&:to_a) + ket_csfs.map(&:to_a), norb)
    # The sub-DRT's input rows are bra_csfs followed by ket_csfs, so the two
    # index sets are the two halves of 0...(bra_csfs.size + ket_csfs.size).
    bra_indices = (0...bra_csfs.size).to_a
    ket_indices = (bra_csfs.size...(bra_csfs.size + ket_csfs.size)).to_a
    sub.block_eval(bra_indices, ket_indices, which: which)
  end

  def csf_subspace_eval(csfs, which: :full)
    raise ArgumentError, "csfs must not be empty" if csfs.empty?
    norb = csfs[0].size
    unless csfs.all? { |r| r.size == norb }
      raise ArgumentError, "all CSF rows must have the same length (norb=#{norb})"
    end
    sub = Drt.from_csfset(csfs.map(&:to_a), norb)
    sub.subspace_eval((0...csfs.size).to_a, which: which)
  end

  class Drt
    KINDS = { "soci" => HEXPR_DRT_SOCI, "foci" => HEXPR_DRT_FOCI }.freeze

    def initialize(kind:, nve:, spin_x2:, ncore:, nval:, next:)
      kind_int =
        case kind
        when Integer then kind
        when String  then KINDS.fetch(kind) {
          raise ArgumentError, "kind must be \"soci\" or \"foci\", got #{kind.inspect}"
        }
        else
          raise ArgumentError, "kind must be Integer or String"
        end

      ptr = HExpr.hexpr_drt_create(kind_int, nve, spin_x2, ncore, nval,
                                   binding.local_variable_get(:next))
      raise "hexpr_drt_create failed: #{HExpr.last_error}" if ptr.null?

      @ptr = FFI::AutoPointer.new(ptr, Drt.method(:release))
    end

    def self.release(ptr)
      HExpr.hexpr_drt_destroy(ptr) unless ptr.null?
    end

    def self.from_csfset(csfs, norb)
      unless csfs.is_a?(Array) && csfs.all? { |row| row.is_a?(Array) }
        raise ArgumentError, "csfs must be Array<Array<Integer>>"
      end
      unless csfs.all? { |row| row.size == norb }
        raise ArgumentError, "each csfs row must have length norb=#{norb}"
      end

      n_csfs = csfs.size
      buf = FFI::MemoryPointer.new(:uint8, n_csfs * norb)
      buf.put_array_of_uint8(0, csfs.flatten)

      ptr = HExpr.hexpr_drt_from_csfset(norb, buf, n_csfs)
      raise "hexpr_drt_from_csfset failed: #{HExpr.last_error}" if ptr.null?

      instance = allocate
      instance.instance_variable_set(:@ptr,
        FFI::AutoPointer.new(ptr, Drt.method(:release)))
      instance
    end

    def to_ptr    = @ptr
    def norb      = HExpr.hexpr_drt_norb(@ptr)
    def nnodes    = HExpr.hexpr_drt_nnodes(@ptr)
    def ncsfs     = HExpr.hexpr_drt_ncsfs(@ptr)
    # Canonical GUGA-walk CSF count. CSFSET: may differ from #ncsfs
    # (superset M >= N, or dedup); SOCI/FOCI: equals #ncsfs. See #walk_csfs.
    def walk_ncsfs = HExpr.hexpr_drt_walk_ncsfs(@ptr)
    def kind      = HExpr.hexpr_drt_kind(@ptr)
    def ncore     = HExpr.hexpr_drt_ncore(@ptr)
    def nval      = HExpr.hexpr_drt_nval(@ptr)
    def next_orb  = HExpr.hexpr_drt_next(@ptr)
    def nve_state = HExpr.hexpr_drt_nve_state(@ptr)
    def spin_x2   = HExpr.hexpr_drt_spin_x2(@ptr)

    # --- CSF step-code matrix ---

    # Returns the CSF step codes as an Array<Array<Integer>>, ncsfs rows x
    # norb columns, each element in 0..3 (E/U/D/F). Mirrors the input format
    # accepted by Drt.from_csfset, so the output round-trips back into it.
    # For a from_csfset sub-DRT the rows ARE the input set in input order:
    # this returns exactly the CSFs you supplied, in the order you supplied
    # them, and #ncsfs is that input count. The bra_pos/ket_pos from
    # block_eval/subspace_eval index this same input array. The internal DRT
    # may expand to a canonical superset (or dedup); that canonical view, in
    # GUGA order, is available via #walk_csfs / #walk_ncsfs.
    def csfs
      n = ncsfs
      return [] if n.zero?
      no  = norb
      ptr = HExpr.hexpr_drt_csfs(@ptr)
      raise "hexpr_drt_csfs failed: #{HExpr.last_error}" if ptr.null?
      flat = ptr.read_array_of_uint8(n * no)
      flat.each_slice(no).to_a
    end

    # Returns the canonical GUGA-walk CSF list as an Array<Array<Integer>>,
    # walk_ncsfs rows x norb columns, each element in 0..3. For a from_csfset
    # sub-DRT this is the internal canonical view, which may be a superset of
    # #csfs (walk_ncsfs >= ncsfs, in GUGA order e<u<d<f) or a dedup; for
    # SOCI/FOCI it equals #csfs.
    def walk_csfs
      n = walk_ncsfs
      return [] if n.zero?
      no  = norb
      ptr = HExpr.hexpr_drt_walk_csfs(@ptr)
      raise "hexpr_drt_walk_csfs failed: #{HExpr.last_error}" if ptr.null?
      flat = ptr.read_array_of_uint8(n * no)
      flat.each_slice(no).to_a
    end

    # --- CSF index lookup ---

    def csf_index(steps)
      buf = FFI::MemoryPointer.new(:uint8, steps.size)
      buf.put_array_of_uint8(0, steps.map(&:to_i))
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_csf_index(@ptr, buf, out)
      raise "hexpr_csf_index failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    # Returns the step codes for CSF row index as an Array<Integer> of
    # length norb (each 0..3). Round-trips back into csf_index.
    def csf_steps(index)
      no  = norb
      buf = FFI::MemoryPointer.new(:uint8, no)
      status = HExpr.hexpr_csf_steps(@ptr, index, buf)
      raise "hexpr_csf_steps failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      buf.read_array_of_uint8(no)
    end

    # --- CSF -> determinant expansion ---

    # Number of determinants in the expansion of CSF index at ms_x2.
    # Exact, not an upper bound: it is C(n_open, n_alpha_open).
    #
    # index follows the same convention as #csf_steps -- the input
    # position, in 0...ncsfs, 0-based as every index in this binding is.
    # So do #pair_eval / #block_eval / #subspace_eval: there is one index
    # space across the whole API.
    #
    # ms_x2 is 2*Ms. Any |Ms| <= S with the same parity as #spin_x2 is
    # accepted, negative values included; Ms is not assumed to equal S.
    def csf_ndets(index, ms_x2)
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_csf_ndets(@ptr, index, ms_x2, out)
      raise "hexpr_csf_ndets failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    # Expand CSF index into Slater determinants at ms_x2.
    #
    # Returns [alpha, beta, coef]. alpha and beta are Array<Array<Integer>>,
    # ndets rows x nwords columns, where nwords = (norb + 63) / 64; coef is
    # an Array<Float> of length ndets. Bit p of word w is spatial orbital
    # 64*w + p, so orbital p (0-based) carries an alpha electron in
    # determinant k iff
    #
    #     alpha[k][p / 64][p % 64] == 1
    #
    # Word splitting is kept even at nwords == 1, and even though Ruby's
    # Integer is arbitrary-precision and could hold the whole mask in one
    # value. That is deliberate, not an oversight: this is a wrapper, and
    # its job is to present the C contract as the C API states it. The
    # other four bindings expose the same (ndets, nwords) shape, and a
    # Ruby-only packed form would have to be read back against a different
    # contract from every other language. A packed accessor can be added
    # later without breaking this one; collapsing this one later could not
    # be undone.
    #
    # Conventions, all fixed by the C API and derived in
    # c/src/csf/csf2det.c:
    #
    #   * Spin coupling: genealogical, in the Eq. (5) order -- the new
    #     one-orbital state couples on the LEFT of the accumulated partial
    #     CSF. A CSF in another coupling scheme will not reproduce these
    #     coefficients.
    #   * Clebsch-Gordan phase: Condon-Shortley.
    #   * Operator order: all alpha in ascending spatial orbital, then all
    #     beta in ascending spatial orbital. An interleaved-by-orbital
    #     convention (some qubit mappings) differs by a
    #     determinant-dependent permutation sign.
    #   * Determinant order: lexicographically ascending in the set of
    #     open-shell positions carrying alpha.
    #   * Normalization: the coefficients of any one CSF square-sum to 1
    #     at every Ms.
    #   * Zero coefficients are emitted, not dropped, so the count is
    #     exact and no threshold enters this API.
    #
    # The OVERALL sign of a CSF is a convention, not a measured quantity:
    # it is whatever the coupling order and the operator order produce.
    # Relative signs -- within one CSF, and between CSFs of one DRT -- are
    # meaningful and are verified against PySCF by
    # validation/csf2det_vs_pyscf.py. The overall sign is not something to
    # assert against.
    def csf_to_dets(index, ms_x2)
      ndets = csf_ndets(index, ms_x2)
      return [[], [], []] if ndets.zero?
      nwords = (norb + 63) / 64
      n = ndets * nwords
      alpha_buf = FFI::MemoryPointer.new(:uint64, n)
      beta_buf  = FFI::MemoryPointer.new(:uint64, n)
      coef_buf  = FFI::MemoryPointer.new(:double, ndets)
      out_count = FFI::MemoryPointer.new(:int32,  1)
      status = HExpr.hexpr_csf_to_dets(@ptr, index, ms_x2, out_count,
                                       alpha_buf, beta_buf, coef_buf)
      raise "hexpr_csf_to_dets failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      got = out_count.read_int32
      unless got == ndets
        raise "hexpr_csf_to_dets wrote #{got} determinants but " \
              "csf_ndets promised #{ndets}"
      end
      [alpha_buf.read_array_of_uint64(n).each_slice(nwords).to_a,
       beta_buf.read_array_of_uint64(n).each_slice(nwords).to_a,
       coef_buf.read_array_of_double(ndets)]
    end

    # --- Node info ---

    # Returns [orb, nve, is] for the given node index.
    def node_info(index)
      orb = FFI::MemoryPointer.new(:int32, 1)
      nve = FFI::MemoryPointer.new(:int32, 1)
      is  = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_drt_node_info(@ptr, index, orb, nve, is)
      raise "hexpr_drt_node_info failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      [orb.read_int32, nve.read_int32, is.read_int32]
    end

    # --- Show ---

    # Dump a human-readable DRT summary to stdout.
    def show
      status = HExpr.hexpr_drt_show(@ptr, HExpr.c_stdout)
      raise "hexpr_drt_show failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      nil
    end

    # --- Per-pair evaluation ---

    def pair_count(bra, ket, which: :full)
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_pair_count(@ptr, bra, ket, which_int(which), out)
      raise "hexpr_pair_count failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    def pair_max_terms(which: :full)
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_pair_max_terms(@ptr, which_int(which), out)
      raise "hexpr_pair_max_terms failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    # bra and ket are input-order CSF indices in 0...ncsfs, the same space
    # as #csf_steps and #csf_index. The returned ij is a packed CANONICAL
    # pair address, not a pair of indices into this DRT
    # (docs/PQRS_INDEX_SPEC.md); do not decode it and feed the result back.
    def pair_eval(bra, ket, which: :full)
      wi    = which_int(which)
      max_t = pair_max_terms(which: which)
      ij_buf    = FFI::MemoryPointer.new(:int32,  max_t)
      pqrs_buf  = FFI::MemoryPointer.new(:int32,  max_t)
      coef_buf  = FFI::MemoryPointer.new(:double, max_t)
      out_count = FFI::MemoryPointer.new(:int32,  1)
      status = HExpr.hexpr_pair_eval(@ptr, bra, ket, wi, out_count,
                                      ij_buf, pqrs_buf, coef_buf)
      raise "hexpr_pair_eval failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      n = out_count.read_int32
      [ij_buf.read_array_of_int32(n),
       pqrs_buf.read_array_of_int32(n),
       coef_buf.read_array_of_double(n)]
    end

    # --- Block evaluation ---

    def block_count(bra_indices, ket_indices, which: :full)
      n_bra = bra_indices.size
      n_ket = ket_indices.size
      out   = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_block_count(@ptr,
                                        int32_buf(bra_indices), n_bra,
                                        int32_buf(ket_indices), n_ket,
                                        which_int(which), out)
      raise "hexpr_block_count failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    def block_max_terms(n_bra, n_ket, which: :full)
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_block_max_terms(@ptr, n_bra, n_ket, which_int(which), out)
      raise "hexpr_block_max_terms failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    # bra_indices/ket_indices hold input-order CSF indices in 0...ncsfs.
    # The returned bra_pos/ket_pos are positions within those arrays, and
    # they, not ij, are how a term is attributed to a pair.
    def block_eval(bra_indices, ket_indices, which: :full)
      wi    = which_int(which)
      n_bra = bra_indices.size
      n_ket = ket_indices.size
      alloc = [block_max_terms(n_bra, n_ket, which: which), 1].max
      bra_pos_buf = FFI::MemoryPointer.new(:int32,  alloc)
      ket_pos_buf = FFI::MemoryPointer.new(:int32,  alloc)
      ij_buf      = FFI::MemoryPointer.new(:int32,  alloc)
      pqrs_buf    = FFI::MemoryPointer.new(:int32,  alloc)
      coef_buf    = FFI::MemoryPointer.new(:double, alloc)
      out_count   = FFI::MemoryPointer.new(:int32,  1)
      status = HExpr.hexpr_block_eval(@ptr,
                                       int32_buf(bra_indices), n_bra,
                                       int32_buf(ket_indices), n_ket,
                                       wi, out_count,
                                       bra_pos_buf, ket_pos_buf,
                                       ij_buf, pqrs_buf, coef_buf)
      raise "hexpr_block_eval failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      n = out_count.read_int32
      [ij_buf.read_array_of_int32(n),
       pqrs_buf.read_array_of_int32(n),
       coef_buf.read_array_of_double(n),
       bra_pos_buf.read_array_of_int32(n),
       ket_pos_buf.read_array_of_int32(n)]
    end

    # --- Subspace evaluation ---

    def subspace_count(indices, which: :full)
      n   = indices.size
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_subspace_count(@ptr, int32_buf(indices), n,
                                           which_int(which), out)
      raise "hexpr_subspace_count failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    def subspace_max_terms(n, which: :full)
      out = FFI::MemoryPointer.new(:int32, 1)
      status = HExpr.hexpr_subspace_max_terms(@ptr, n, which_int(which), out)
      raise "hexpr_subspace_max_terms failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      out.read_int32
    end

    # indices holds input-order CSF indices in 0...ncsfs; bra_pos/ket_pos
    # come back as positions within that array.
    def subspace_eval(indices, which: :full)
      wi    = which_int(which)
      n     = indices.size
      alloc = [subspace_max_terms(n, which: which), 1].max
      bra_pos_buf = FFI::MemoryPointer.new(:int32,  alloc)
      ket_pos_buf = FFI::MemoryPointer.new(:int32,  alloc)
      ij_buf      = FFI::MemoryPointer.new(:int32,  alloc)
      pqrs_buf    = FFI::MemoryPointer.new(:int32,  alloc)
      coef_buf    = FFI::MemoryPointer.new(:double, alloc)
      out_count   = FFI::MemoryPointer.new(:int32,  1)
      status = HExpr.hexpr_subspace_eval(@ptr, int32_buf(indices), n,
                                          wi, out_count,
                                          bra_pos_buf, ket_pos_buf,
                                          ij_buf, pqrs_buf, coef_buf)
      raise "hexpr_subspace_eval failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      n_out = out_count.read_int32
      [ij_buf.read_array_of_int32(n_out),
       pqrs_buf.read_array_of_int32(n_out),
       coef_buf.read_array_of_double(n_out),
       bra_pos_buf.read_array_of_int32(n_out),
       ket_pos_buf.read_array_of_int32(n_out)]
    end

    private

    def which_int(which)
      case which
      when :full then HExpr::HEXPR_EXPR_FULL
      when :one  then HExpr::HEXPR_EXPR_ONE
      else raise ArgumentError, "which must be :full or :one, got #{which.inspect}"
      end
    end

    def int32_buf(arr)
      a = arr.map(&:to_i)
      buf = FFI::MemoryPointer.new(:int32, [a.size, 1].max)
      buf.put_array_of_int32(0, a) unless a.empty?
      buf
    end
  end

  class Expression
    def initialize(drt)
      @parent = drt
      ptr = HExpr.hexpr_expr_build(drt.to_ptr)
      raise "hexpr_expr_build failed: #{HExpr.last_error}" if ptr.null?
      @ptr = FFI::AutoPointer.new(ptr, Expression.method(:release))
    end

    def self.release(ptr)
      HExpr.hexpr_expr_destroy(ptr) unless ptr.null?
    end

    def nterms_full = HExpr.hexpr_expr_nterms_full(@ptr)
    def nterms_one  = HExpr.hexpr_expr_nterms_one(@ptr)

    def read_full = read_arrays(HEXPR_EXPR_FULL, nterms_full)
    def read_one  = read_arrays(HEXPR_EXPR_ONE,  nterms_one)

    # Dump the full (1e+2e) expression to stdout.
    def show_full
      status = HExpr.hexpr_expr_show_full(@ptr, HExpr.c_stdout)
      raise "hexpr_expr_show_full failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      nil
    end

    # Dump the 1-electron expression subset to stdout.
    def show_one
      status = HExpr.hexpr_expr_show_one(@ptr, HExpr.c_stdout)
      raise "hexpr_expr_show_one failed: #{HExpr.last_error}" unless status == HExpr::HEXPR_OK
      nil
    end

    private

    def read_arrays(selector, n)
      ij_buf   = FFI::MemoryPointer.new(:int32,  n)
      pqrs_buf = FFI::MemoryPointer.new(:int32,  n)
      coef_buf = FFI::MemoryPointer.new(:double, n)
      status = HExpr.hexpr_expr_read(@ptr, selector, ij_buf, pqrs_buf, coef_buf)
      raise "hexpr_expr_read failed: #{HExpr.last_error}" unless status == HEXPR_OK
      [
        ij_buf.read_array_of_int32(n),
        pqrs_buf.read_array_of_int32(n),
        coef_buf.read_array_of_double(n),
      ]
    end
  end
end
