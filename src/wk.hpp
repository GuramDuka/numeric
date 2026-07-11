/*-
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Guram Duka
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
//------------------------------------------------------------------------------
/// @file wk.hpp
/// @brief Cryptographic sponge constructions (cdc256, cdc512) and a
///        condition-variable-based thread pool (ThreadPoolT)
///
/// @details This header provides two unrelated components used by the Numeric
///   arbitrary-precision arithmetic library:
///
///   1. **cdc256 / cdc512** — Custom cryptographic sponge constructions by
///      Guram Duka.  These provide a 256-bit (resp. 512-bit) internal state
///      with an ARX-based permutation.  Used for deterministic pseudo-random
///      number generation within the library (currently commented-out usage
///      in ii.hpp).  Not a NIST-standardised construction (cf. SHA-3).
///
///   2. **ThreadPoolT\<FnType\>** — A simple condition-variable-based thread
///      pool.  Each worker thread waits on its own condition variable.  Tasks
///      are enqueued from a single producer and distributed round-robin across
///      workers for parallel big-integer multiplication.
///
/// Complexity and thread-safety notes are provided per-class.
///
/// @see Guram Duka, "CDC Cryptographic Sponge" (custom design)
/// @see ThreadPoolT implementation based on standard condition-variable
///      signalling pattern
//------------------------------------------------------------------------------
#ifndef WK_HPP_INCLUDED
#define WK_HPP_INCLUDED
//------------------------------------------------------------------------------
#include <functional>
#include <thread>
#include <queue>
#include <mutex>
#include <shared_mutex>
#include <memory>
#include <cstring>
#include <condition_variable>
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------

//##############################################################################
// CDC256 — 256-bit Cryptographic Sponge
//##############################################################################

//------------------------------------------------------------------------------
/// @brief 256-bit cryptographic sponge construction
///
/// @details A custom ARX (Add-Rotate-XOR) based sponge with 256 bits of
///   internal state (eight 64-bit words).  The sponge follows the standard
///   absorb-squeeze paradigm:
///     - init()   : initialise state to hardcoded constants
///     - update() : absorb arbitrary-length input (padded with zeros)
///     - (final)  : squeeze (commented out; digest is in data_ directly)
///
///   The core permutation is `shuffle()` — 8 rounds of ARX operations
///   interleaving subtraction, XOR with rotated values, and addition.
///   Each 64-byte input block is first XOR-merged into the state via
///   `update(const data &)` (the "absorb" step), then the state is
///   permuted via `shuffle()`.
///
///   The input length (in bytes) is also absorbed after data to provide
///   domain separation (a simple "length padding").
///
/// Complexity:
///   - Per 64-byte block: O(1) — 8 rounds × 8 word operations
///   - Total: O(n) for n bytes of input
///
/// Thread safety: NOT thread-safe.  The context holds mutable state
///   (`data_` and `position_`).  Do not call update() concurrently.
///
/// @warning This is a custom sponge, NOT a standardised hash like SHA-3.
///   It is used for deterministic pseudo-randomness within the library,
///   not for security-critical hashing.
///
/// @see Guram Duka — CDC cryptographic sponge design
/// @see Bertoni et al., "Keccak Sponge Construction" (SHA-3) for the
///      standardised sponge paradigm
//------------------------------------------------------------------------------
struct cdc256 {
    //--------------------------------------------------------------------------
    /// @brief Internal state word bundle (256 bits = 8 × uint64_t)
    ///
    /// @details The eight words a–h form the 256-bit sponge state.
    ///   The `shuffle()` method applies one round of the ARX permutation.
    ///   The `update(const data &)` method absorbs an input block into
    ///   the state by subtracting the input words and XOR-rotating.
    ///
    /// Complexity:
    ///   - shuffle() : O(1), 8 rounds of 3 word-operations each
    ///   - update()  : O(1), same structure with input words
    ///
    /// Thread safety: Not thread-safe; operates on `this` directly.
    //--------------------------------------------------------------------------
    struct data {
        uint64_t a, b, c, d, e, f, g, h;

        //----------------------------------------------------------------------
        /// @brief Apply one round of the ARX permutation to the state
        ///
        /// @details Each round performs 8 parallelised ARX triples:
        ///   ```
        ///   word[i] -= word[i+4];           // subtraction (additive)
        ///   word[i+4] ^= word[i-1] >>/<< k; // XOR with rotated neighbour
        ///   word[i-1] += word[i];           // addition
        ///   ```
        ///   where indices wrap modulo 8.  Rotation amounts vary per step
        ///   (9, 9, 23, 15, 14, 20, 17, 14).  This is NOT the Keccak-f
        ///   permutation; it is a custom design.
        ///
        /// Precondition: none (works on any uint64_t values)
        /// Postcondition: state words are mixed according to the ARX schedule
        ///
        /// Complexity: O(1), 24 ALU operations
        //----------------------------------------------------------------------
        constexpr void shuffle() {
            a -= e; f ^= h >>  9; h += a;
            b -= f; g ^= a <<  9; a += b;
            c -= g; h ^= b >> 23; b += c;
            d -= h; a ^= c << 15; c += d;
            e -= a; b ^= d >> 14; d += e;
            f -= b; c ^= e << 20; e += f;
            g -= c; d ^= f >> 17; f += g;
            h -= d; e ^= g << 14; g += h;
        }

        //----------------------------------------------------------------------
        /// @brief Absorb one 256-bit input block into the state
        ///
        /// @details The input block @p v is XOR-merged into the state using
        ///   the same pattern as shuffle() but with input words in place of
        ///   the shift-rotation.  This is the "absorb" phase: the input
        ///   word is subtracted, the shifted/XORed state word is taken from
        ///   the corresponding input word, and the addition uses the input
        ///   word.
        ///
        /// Precondition: none
        /// Postcondition: state reflects the absorbed input block
        ///
        /// Complexity: O(1), 24 ALU operations
        ///
        /// @param v The 256-bit input block to absorb
        //----------------------------------------------------------------------
        constexpr void update(const data & v) {
            a -= v.e; f ^= v.h >>  9; h += v.a;
            b -= v.f; g ^= v.a <<  9; a += v.b;
            c -= v.g; h ^= v.b >> 23; b += v.c;
            d -= v.h; a ^= v.c << 15; c += v.d;
            e -= v.a; b ^= v.d >> 14; d += v.e;
            f -= v.b; c ^= v.e << 20; e += v.f;
            g -= v.c; d ^= v.f >> 17; f += v.g;
            h -= v.d; e ^= v.g << 14; g += v.h;
        }
    };

    //--------------------------------------------------------------------------
    /// @brief Construct a default-initialised cdc256 context
    ///
    /// @details Initialises state with the hardcoded constant vector and
    ///   resets the byte-position counter to zero.  Equivalent to calling
    ///   init().
    //--------------------------------------------------------------------------
    cdc256() { init(); }

    //--------------------------------------------------------------------------
    /// @brief Initialise/reset the sponge to its initial state
    ///
    /// @details Sets the eight 64-bit state words to hardcoded constants
    ///   (nothing-up-my-sleeve numbers from fractional parts of square
    ///   roots or similar — though the provenance is undocumented) and
    ///   zeroes the byte-position counter.
    ///
    /// Complexity: O(1) — 8 assignments + 1 assignment
    ///
    /// Postcondition: data_ == init_value, position_ == 0
    //--------------------------------------------------------------------------
    void init() {
        static constexpr data init_value = {
            UINT64_C(0xA640524A5B44F1FC),
            UINT64_C(0xC535059705F0BB7E),
            UINT64_C(0xC8ED76CF6B6EA626),
            UINT64_C(0x531D1E8E254EA59E),
            UINT64_C(0x8C0FE7F3E46E2A80),
            UINT64_C(0x1C53F41FD1E3A7F8),
            UINT64_C(0x08D4DEAAA1C33335),
            UINT64_C(0x4C592980FBE9B011)
        };
        data_     = init_value;
        position_ = 0;
    }

    //--------------------------------------------------------------------------
    /// @brief Absorb input bytes into the sponge
    ///
    /// @details Processes the input buffer of @p size bytes in 64-byte
    ///   ( = sizeof(data) ) blocks.  Each full block is absorbed via
    ///   data::update() followed by data::shuffle() (absorb-permute).
    ///   Any trailing partial block is zero-padded before absorbing.
    ///   Finally, the total absorbed byte count (@p position_) is itself
    ///   absorbed as a domain-separation tag.
    ///
    /// Precondition: @p p must point to @p size readable bytes, or
    ///   size == 0 (in which case no read occurs).
    /// Postcondition: The sponge state reflects the entire input, with
    ///   length padding appended.
    ///
    /// Complexity: O(n) where n = size bytes processed
    ///
    /// Thread safety: NOT thread-safe; modifies internal mutable state.
    ///
    /// @param p    Pointer to the input byte buffer
    /// @param size Number of bytes to absorb
    //--------------------------------------------------------------------------
    void update(const void * p, uintptr_t size) {
        data pad;

        auto dig = [](data & d, const data & p) {
            d.update(p);
            d.shuffle();
        };

        position_ += size;

        while (size >= sizeof(data)) {
            dig(data_, *static_cast<const data *>(p));
            p = static_cast<const uint8_t *>(p) + sizeof(data);
            size -= sizeof(data);
        }

        if (size > 0) {
            std::memcpy(&pad, p, size);
            std::memset(reinterpret_cast<uint8_t *>(&pad) + size, 0,
                        sizeof(data) - size);
            dig(data_, pad);
        }

        pad.a = pad.b = pad.c = pad.d =
        pad.e = pad.f = pad.g = pad.h = position_;
        dig(data_, pad);
    }

    //--------------------------------------------------------------------------
    /// @brief Absorb the contents of an STL-like container
    ///
    /// @details Convenience overload: calls update(container.data(),
    ///   container.size()).  Works with std::vector, std::array,
    ///   std::string, or any contiguous container providing .data()
    ///   and .size().
    ///
    /// @tparam T Contiguous container type with data()/size()
    /// @param container The container whose bytes to absorb
    //--------------------------------------------------------------------------
    template <typename T>
    void update(T & container) {
        update(container.data(), container.size());
    }

    // -- Commented-out finalisation -----------------------------------------
    // The original design included:
    //   std::array<uint8_t, 64> final() {
    //       return *reinterpret_cast<std::array<uint8_t, 64> *>(&data_);
    //   }
    // The digest is still accessible via the data_ member.  This method was
    // left unused and is preserved only for reference.
    // -----------------------------------------------------------------------

    data     data_;      ///< Sponge internal state (256 bits)
    uint64_t position_;  ///< Total bytes absorbed so far
};

//------------------------------------------------------------------------------
//##############################################################################
// CDC512 — 512-bit Cryptographic Sponge
//##############################################################################

//------------------------------------------------------------------------------
/// @brief 512-bit cryptographic sponge construction
///
/// @details An enlarged variant of cdc256 with 512 bits of internal state
///   (sixteen 64-bit words a–p).  Follows the same absorb-permute sponge
///   paradigm.  The `shuffle()` permutation uses two interleaved 256-bit
///   halves with cross-lane diffusion.
///
///   The input is processed in 128-byte blocks.  The length padding scheme
///   is identical to cdc256.
///
/// Complexity:
///   - Per 128-byte block: O(1) — 16 rounds of ARX operations
///   - Total: O(n) for n bytes of input
///
/// Thread safety: NOT thread-safe.
///
/// @see cdc256 for the general sponge design notes
//------------------------------------------------------------------------------
struct cdc512 {
    //--------------------------------------------------------------------------
    /// @brief Internal state word bundle (512 bits = 16 × uint64_t)
    ///
    /// @details The sixteen words a–p form the 512-bit sponge state.
    ///   The permutation applies two 256-bit halves with cross-coupling
    ///   between the lower half (a–h) and upper half (i–p): the first
    ///   8 steps use upper-half inputs to modify the lower half, then the
    ///   next 8 steps use lower-half inputs to modify the upper half.
    ///
    /// Complexity: O(1), 48 ALU operations per call
    ///
    /// Thread safety: Not thread-safe.
    //--------------------------------------------------------------------------
    struct data {
        uint64_t a, b, c, d, e, f, g, h;
        uint64_t i, j, k, l, m, n, o, p;

        //------------------------------------------------------------------
        /// @brief Apply one round of the 512-bit ARX permutation
        ///
        /// @details The permutation operates in two halves of 8 steps each:
        ///   - Steps 1–8  (a–h ← f(a–h, i–p)): lower half absorbs upper
        ///   - Steps 9–16 (i–p ← f(i–p, a–h)): upper half absorbs lower
        ///
        ///   Each step follows the same ARX triple pattern as cdc256 but
        ///   with cross-half sources for the XOR-rotate operand.
        ///
        /// Precondition: none
        /// Postcondition: all 16 state words are mixed
        ///
        /// Complexity: O(1), 48 ALU operations
        //------------------------------------------------------------------
        constexpr void shuffle() {
            a -= i; i ^= h >> 10; p += a;
            b -= j; j ^= g << 10; o += b;
            c -= k; k ^= f >> 23; n += c;
            d -= l; l ^= e << 15; m += d;
            e -= m; m ^= d >> 14; l += e;
            f -= n; n ^= c << 20; k += f;
            g -= o; o ^= b >> 17; j += g;
            h -= p; p ^= a << 14; i += h;
            i -= a; a ^= p >> 10; h += i;
            j -= b; b ^= o << 10; g += j;
            k -= c; c ^= n >> 23; f += k;
            l -= d; d ^= m << 15; e += l;
            m -= e; e ^= l >> 14; d += m;
            n -= f; f ^= k << 20; c += n;
            o -= g; g ^= j >> 17; b += o;
            p -= h; h ^= i << 14; a += p;
        }

        //------------------------------------------------------------------
        /// @brief Absorb one 512-bit input block into the state
        ///
        /// @details XOR-merges the input block @p v into the 512-bit state.
        ///   The pattern mirrors shuffle() but substitutes input words for
        ///   rotated state words.
        ///
        /// Precondition: none
        /// Postcondition: state reflects the absorbed input block
        ///
        /// Complexity: O(1), 48 ALU operations
        ///
        /// @param v The 512-bit input block to absorb
        //------------------------------------------------------------------
        constexpr void update(const data & v) {
            a -= v.i; i ^= v.h >> 10; p += v.a;
            b -= v.j; j ^= v.g << 10; o += v.b;
            c -= v.k; k ^= v.f >> 23; n += v.c;
            d -= v.l; l ^= v.e << 15; m += v.d;
            e -= v.m; m ^= v.d >> 14; l += v.e;
            f -= v.n; n ^= v.c << 20; k += v.f;
            g -= v.o; o ^= v.b >> 17; j += v.g;
            h -= v.p; p ^= v.a << 14; i += v.h;
            i -= v.a; a ^= v.p >> 10; h += v.i;
            j -= v.b; b ^= v.o << 10; g += v.j;
            k -= v.c; c ^= v.n >> 23; f += v.k;
            l -= v.d; d ^= v.m << 15; e += v.l;
            m -= v.e; e ^= v.l >> 14; d += v.m;
            n -= v.f; f ^= v.k << 20; c += v.n;
            o -= v.g; g ^= v.j >> 17; b += v.o;
            p -= v.h; h ^= v.i << 14; a += v.p;
        }
    };

    //--------------------------------------------------------------------------
    /// @brief Construct a default-initialised cdc512 context
    ///
    /// @details Initialises state with hardcoded constants and resets the
    ///   byte-position counter to zero.  Equivalent to calling init().
    //--------------------------------------------------------------------------
    cdc512() { init(); }

    //--------------------------------------------------------------------------
    /// @brief Initialise/reset the 512-bit sponge to its initial state
    ///
    /// @details Sets the sixteen 64-bit state words to hardcoded constants
    ///   (the first eight are the same as cdc256; the remaining eight are
    ///   additional arbitrary-looking constants).  Resets position_ to 0.
    ///
    /// Complexity: O(1) — 16 assignments + 1 assignment
    ///
    /// Postcondition: data_ == init_value, position_ == 0
    //--------------------------------------------------------------------------
    void init() {
        static constexpr data init_value = {
            UINT64_C(0xA640524A5B44F1FC),
            UINT64_C(0xC535059705F0BB7E),
            UINT64_C(0xC8ED76CF6B6EA626),
            UINT64_C(0x531D1E8E254EA59E),
            UINT64_C(0x8C0FE7F3E46E2A80),
            UINT64_C(0x1C53F41FD1E3A7F8),
            UINT64_C(0x08D4DEAAA1C33335),
            UINT64_C(0x4C592980FBE9B011),
            UINT64_C(0x992E367BE6F0EA1E),
            UINT64_C(0x71DCF41FFACC283F),
            UINT64_C(0xC9581F48D85ABD75),
            UINT64_C(0xE4B93335FF1CE990),
            UINT64_C(0xE51D6424EFEC1E01),
            UINT64_C(0x353867A0E66C2A39),
            UINT64_C(0xA8DBF7B782226B67),
            UINT64_C(0x9F8B7F0DC254488E)
        };
        data_     = init_value;
        position_ = 0;
    }

    //--------------------------------------------------------------------------
    /// @brief Absorb input bytes into the 512-bit sponge
    ///
    /// @details Processes the input in 128-byte blocks.  Identical logic
    ///   to cdc256::update() but with the larger block size.
    ///
    /// Precondition: @p p must point to @p size readable bytes, or
    ///   size == 0.
    /// Postcondition: The sponge state reflects the entire input with
    ///   length padding.
    ///
    /// Complexity: O(n) where n = size bytes processed
    ///
    /// Thread safety: NOT thread-safe.
    ///
    /// @param p    Pointer to the input byte buffer
    /// @param size Number of bytes to absorb
    //--------------------------------------------------------------------------
    void update(const void * p, uintptr_t size) {
        data pad;

        auto dig = [](data & d, const data & p) {
            d.update(p);
            d.shuffle();
        };

        position_ += size;

        while (size >= sizeof(data)) {
            dig(data_, *static_cast<const data *>(p));
            p = static_cast<const uint8_t *>(p) + sizeof(data);
            size -= sizeof(data);
        }

        if (size > 0) {
            std::memcpy(&pad, p, size);
            std::memset(reinterpret_cast<uint8_t *>(&pad) + size, 0,
                        sizeof(data) - size);
            dig(data_, pad);
        }

        pad.a = pad.b = pad.c = pad.d =
        pad.e = pad.f = pad.g = pad.h =
        pad.i = pad.j = pad.k = pad.l =
        pad.m = pad.n = pad.o = pad.p = position_;
        dig(data_, pad);
    }

    //--------------------------------------------------------------------------
    /// @brief Absorb the contents of an STL-like container
    ///
    /// @tparam T Contiguous container type with data()/size()
    /// @param container The container whose bytes to absorb
    //--------------------------------------------------------------------------
    template <typename T>
    void update(T & container) {
        update(container.data(), container.size());
    }

    // -- Commented-out finalisation -----------------------------------------
    // Original design:
    //   std::array<uint8_t, 128> final() {
    //       return *reinterpret_cast<std::array<uint8_t, 128> *>(&data_);
    //   }
    // -----------------------------------------------------------------------

    data     data_;      ///< Sponge internal state (512 bits)
    uint64_t position_;  ///< Total bytes absorbed so far
};

//------------------------------------------------------------------------------
//##############################################################################
// ThreadPoolT — Condition-Variable-Based Thread Pool
//##############################################################################

//------------------------------------------------------------------------------
/// @brief A thread pool for parallel task execution with round-robin
///        distribution
///
/// @details Workers are created at construction time, each waiting on its
///   own condition variable.  Tasks are enqueued by a single producer,
///   then a barrier (`wait()`) signals the first N workers to begin.
///   Each worker processes tasks at indices congruent to its worker index
///   modulo the number of threads (round-robin striding).
///
///   This pool is designed specifically for parallel big-integer
///   multiplication where independent sub-products can be computed
///   concurrently.  The single-producer, barrier-synchronised model
///   suits this workload well.
///
/// @tparam FnType Function signature type (default usage: void())
///   Internally wrapped in std::function<FnType>.
///
/// Thread safety:
///   - `run_task()`: NOT thread-safe; call from a single thread before
///     `wait()`.
///   - `wait()`: Can be called from any thread (serialises on the mutex).
///   - `run()` (private): Each worker runs on its own thread.
///
/// Known issues:
///   - `counter_` is set in `wait()` BEFORE the mutex is acquired
///     (pre-existing race-condition note; the producer is assumed to
///     be single-threaded so this is benign in practice).
///
/// Complexity:
///   - Task dispatch: O(1) per task
///   - Work distribution: O(n / threads) per worker
///   - Destruction: O(threads)
///
/// @see std::thread, std::condition_variable for the underlying primitives
/// @see Michael, "Scalable Lock-Free Dynamic Memory Allocation" for the
///      work-stealing inspiration (though this pool does NOT steal tasks)
//------------------------------------------------------------------------------
template <typename FnType>
class ThreadPoolT {
public:
    /// Type alias for the stored task wrapper
    using TaskType = std::function<FnType>;

    //--------------------------------------------------------------------------
    /// @brief Destructor — signals shutdown and joins all worker threads
    ///
    /// @details Sets the shutdown flag, notifies ALL workers (they will
    ///   wake, see shutdown_ == true, and exit), then joins every thread.
    ///
    /// Precondition: No threads may be processing tasks at destruction
    ///   (call wait() first if tasks were enqueued).
    /// Postcondition: All worker threads have exited and been joined.
    ///
    /// Complexity: O(threads) — notify + join each thread
    //--------------------------------------------------------------------------
    ~ThreadPoolT() {
        shutdown_ = true;

        for (auto & cv : cvs_)
            cv->notify_all();

        for (auto & t : threads_)
            t->join();
    }

    //--------------------------------------------------------------------------
    /// @brief Construct the thread pool and spawn worker threads
    ///
    /// @details Creates @p threads worker threads (auto-detected via
    ///   std::thread::hardware_concurrency() when 0 is passed).  Each
    ///   worker is given its own condition_variable and a run-flag.
    ///   Task storage is pre-reserved to avoid reallocations.
    ///
    /// Precondition: none
    /// Postcondition: `threads_` workers are alive and waiting.
    ///
    /// Complexity: O(threads) — creation and setup
    ///
    /// @param threads Number of worker threads.  Default 0 means
    ///   std::thread::hardware_concurrency().
    //--------------------------------------------------------------------------
    ThreadPoolT(uintptr_t threads = 0)
        : threads_()
        , cvs_()
        , runs_()
        , wv_()
        , mutex_()
        , tasks_()
        , counter_(0)
        , done_(false)
        , shutdown_(false)
    {
        if (threads == 0)
            threads = std::thread::hardware_concurrency();

        tasks_.reserve(threads * 2);

        for (uintptr_t i = 0; i < threads; ++i) {
            cvs_.emplace_back(std::make_unique<std::condition_variable>());
            runs_.emplace_back(false);
        }

        while (threads-- > 0) {
            threads_.emplace_back(
                std::make_unique<std::thread>(
                    &ThreadPoolT<FnType>::run, this, threads_.size()));
        }
    }

    //--------------------------------------------------------------------------
    /// @brief Enqueue a task for parallel execution
    ///
    /// @details Appends @p task to the shared task vector.  Tasks are
    ///   distributed to workers in round-robin fashion during `wait()`.
    ///   This method does NOT notify workers; call `wait()` to begin
    ///   processing.
    ///
    /// Precondition: May only be called from one thread at a time (NOT
    ///   thread-safe).  Must be called before `wait()`.
    /// Postcondition: @p task is appended to the internal queue.
    ///
    /// Complexity: O(1) amortised
    ///
    /// @param task The callable to enqueue
    //--------------------------------------------------------------------------
    void run_task(const TaskType & task) {
        // NOTE: The original code had a commented-out mutex lock here.
        // run_task() is deliberately NOT thread-safe — the single-producer
        // model avoids the locking overhead on the hot path.
        tasks_.emplace_back(task);
    }

    //--------------------------------------------------------------------------
    /// @brief Signal worker threads and wait for all tasks to complete
    ///
    /// @details Calculates the number of workers to engage (at most
    ///   min(tasks.size(), threads.size())).  Sets the corresponding
    ///   run-flags and notifies the first `counter_` workers.  Then
    ///   blocks on the done condition variable until all workers finish.
    ///   Finally clears the task vector and resets the done flag.
    ///
    /// Precondition:
    ///   - All tasks have been enqueued via run_task()
    ///   - No concurrent calls to wait()
    /// Postcondition:
    ///   - All enqueued tasks have been executed
    ///   - Task vector is cleared
    ///   - Internal flags are reset for the next batch
    ///
    /// Known issue: `counter_` is assigned (line ~279) BEFORE the mutex
    ///   is acquired.  The original author's intent was to let workers
    ///   check `counter_` without serialising; however, this is a data
    ///   race under the C++ memory model.  In practice the single-
    ///   producer design prevents concurrent access.
    ///
    /// Complexity: O(threads) for notifications + O(p) for producer block
    ///
    /// @see https://en.cppreference.com/w/cpp/thread/condition_variable
    //--------------------------------------------------------------------------
    void wait() {
        counter_ = tasks_.size() >= threads_.size()
                     ? threads_.size()
                     : tasks_.size();

        if (counter_ != 0) {
            std::unique_lock<decltype(mutex_)> locker(mutex_);
            for (uintptr_t i = 0; i < counter_; ++i) {
                runs_[i] = true;
                cvs_[i]->notify_one();
            }
            wv_.wait(locker, [&] { return done_; });
            tasks_.clear();
        }

        done_ = false;
    }

private:
    // -- Members ---------------------------------------------------------

    std::vector<std::unique_ptr<std::thread>>              threads_;  ///< Worker threads
    std::vector<std::unique_ptr<std::condition_variable>>  cvs_;      ///< Per-worker condition variables
    std::vector<bool>  runs_;     ///< Per-worker run-flag (signalled by producer)
    std::condition_variable wv_;  ///< Condition variable for wait() completion
    std::mutex             mutex_;  ///< Mutex protecting shared state
    std::vector<TaskType>  tasks_;  ///< Queued tasks
    uintptr_t  counter_;   ///< Number of workers engaged in current batch
    bool       done_;      ///< Completion flag (set by last finishing worker)
    bool       shutdown_;  ///< Shutdown flag (set by destructor)

    //--------------------------------------------------------------------------
    /// @brief Worker thread entry point
    ///
    /// @details Each worker loops until shutdown:
    ///   1. Wait on its condition variable until its run-flag is set or
    ///      shutdown is requested.
    ///   2. If the run-flag is set:
    ///      a. Release the lock (so other workers can proceed).
    ///      b. Execute tasks at indices congruent to `index` modulo
    ///         `threads_.size()` (round-robin stride).
    ///      c. Re-acquire the lock.
    ///      d. Decrement `counter_`.  If counter_ reaches 0, the last
    ///         worker sets `done_ = true` and notifies the producer.
    ///      e. Clear its run-flag.
    ///   3. If shutdown is requested, exit the loop.
    ///
    /// Precondition: `index < threads_.size()`
    /// Postcondition: On shutdown, thread exits.  On task completion,
    ///   worker returns to waiting state.
    ///
    /// @param index The worker's index in the thread vector (0-based)
    //--------------------------------------------------------------------------
    void run(uintptr_t index) {
        std::unique_lock<decltype(mutex_)> locker(mutex_);

        while (!shutdown_) {
            cvs_[index]->wait(locker,
                              [&] { return runs_[index] || shutdown_; });

            if (runs_[index]) {
                locker.unlock();

                bool did_work = false;

                for (auto i = index; i < tasks_.size(); i += threads_.size()) {
                    tasks_[i]();
                    did_work = true;
                }

                locker.lock();

                if (did_work && --counter_ == 0) {
                    done_ = true;
                    wv_.notify_all();
                }

                runs_[index] = false;
            }
        }
    }
};

//------------------------------------------------------------------------------
/// @brief Concrete thread pool type for void-returning tasks
///
/// @details This is the global pool type used by `operator*` in ii.hpp
///   for parallel big-integer multiplication.
///
/// @see ii.hpp line ~1135:
///   static ThreadPool pool(thread_hardware_concurrency);
//------------------------------------------------------------------------------
using ThreadPool = ThreadPoolT<void()>;

//------------------------------------------------------------------------------
} // namespace nn (NaturalNumbers)
//------------------------------------------------------------------------------
#endif // WK_HPP_INCLUDED
//------------------------------------------------------------------------------
