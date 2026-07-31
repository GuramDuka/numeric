/*-
 * The MIT License (MIT)
 *
 * Copyright (c) 2014, 2015, 2016 Guram Duka
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
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
/// @file id.hpp
/// @namespace nn (NaturalNumbers)
/// @brief Arbitrary-precision integer data structures and core arithmetic
///
/// This file defines the foundational data types and algorithms for the
/// Numeric library's arbitrary-precision integer support. The design follows
/// principles described in Knuth's "The Art of Computer Programming",
/// Volume 2: Seminumerical Algorithms (3rd ed.), Section 4.3.
///
/// Key components:
///   1. Word type selection via SIZEOF_WORD — compile-time dispatch to
///      the appropriate unsigned/signed integer widths based on the target
///      architecture's native word size (64-bit on x86_64, 32-bit on x86, etc.).
///   2. nn_integer_data — ref-counted, sign-magnitude representation of
///      arbitrary-precision integers with a variable-length data array
///      allocated through a custom TLSF (Two-Level Segregated Fit) allocator.
///   3. Schoolbook addition/subtraction algorithms (O(n) time) with
///      length-dispatching for optimal performance across operand sizes.
///   4. Arithmetic shift operations (isal/isar) using overlapping dword
///      accesses to propagate bits across word boundaries in a single pass.
///   5. Sign-extending bitwise operations (iand/ior/ixor/inot) that
///      handle mismatched operand lengths by extending the shorter operand's
///      sign word.
///   6. count_trailing_zeros using __builtin_ctz[ll] (GCC/Clang) with an
///      efficient binary-search fallback for platforms without hardware
///      count-trailing-zeros.
///   7. constexpr hex string conversion utilities for compile-time
///      initialization of numeric constants from hexadecimal literals.
///
/// Data flow:
///   tlsf.hpp  →  id.hpp  →  wk.hpp  →  ii.hpp  →  nn.hpp
///
/// References:
///   - Knuth, D.E. "The Art of Computer Programming", Vol. 2 (3rd ed.),
///     Section 4.3.1: "The Classical Algorithms" (addition, subtraction).
///   - Warren, H.S. "Hacker's Delight" (2nd ed.), Chapter 2: "Basics"
///     (sign-magnitude representation, shifts, bitwise operations).
///   - Massey, J.L. "Shift-register synthesis and BCH decoding" (1969)
///     for de Bruijn sequence concepts in bit-scanning.
//------------------------------------------------------------------------------
#ifndef ID_HPP_INCLUDED
#define ID_HPP_INCLUDED
//------------------------------------------------------------------------------
#include "tlsf.hpp"
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Returns the larger of two values
///
/// @details A type-safe, constexpr-compatible maximum function. Unlike
///   std::max (which requires both arguments to have the same type), this
///   variant accepts mixed types A and B, returning by const reference to A.
///   The comparison promotes B to A's type via the built-in > operator.
///
/// @tparam A The type of the first operand (return type)
/// @tparam B The type of the second operand (must be comparable with A)
/// @param a First operand
/// @param b Second operand
/// @return `a` if `a > b`, otherwise a reference to `a` (effectively `a`)
/// @note The return is always a const reference to `a`, regardless of which
///   value is logically larger. This is safe because both parameters are
///   const references with lifetimes exceeding the call.
///
/// @complexity O(1)
// ---------------------------------------------------------------------------
template <typename A, typename B>
static inline constexpr const A & imax(const A & a, const B & b)
{
    return a > b ? a : b;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Returns the smaller of two values
///
/// @details A type-safe, constexpr-compatible minimum function. Analogous to
///   imax but returns the smaller operand. Accepts mixed types.
///
/// @tparam A The type of the first operand (return type)
/// @tparam B The type of the second operand
/// @param a First operand
/// @param b Second operand
/// @return `a` if `a < b`, otherwise a reference to `a`
///
/// @complexity O(1)
// ---------------------------------------------------------------------------
template <typename A, typename B>
static inline constexpr const A & imin(const A & a, const B & b)
{
    return a < b ? a : b;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Compile-time word size selection
///
/// @details The SIZEOF_WORD macro controls the fundamental integer width used
///   for arbitrary-precision arithmetic. The word size is chosen to match the
///   target architecture's native register width:
///
///   | SIZEOF_WORD | word type | dword type    | Platform              |
///   |-------------|-----------|---------------|-----------------------|
///   | 1           | uint8_t   | uint16_t      | 8-bit microcontrollers|
///   | 2           | uint16_t  | uint32_t      | 16-bit embedded       |
///   | 4           | uint32_t  | uint64_t      | x86 (32-bit), ARM     |
///   | 8           | uint64_t  | __uint128_t   | x86_64, aarch64       |
///
///   On 64-bit architectures with GCC/Clang (__x86_64__ macro), the double-word
///   type is __uint128_t, a compiler-provided 128-bit extension that the
///   back-end maps to the `mul`/`div` instruction pair for efficient
///   carry/borrow propagation. This is critical for the schoolbook algorithms
///   which rely on dword-wide intermediate accumulators.
///
///   The fallback path (lines 92-100) defaults to 32-bit words for
///   unrecognized architectures, ensuring portability at the cost of
///   performance (more word iterations per operation).
///
///   Reference: Knuth, TAOCP Vol. 2, Section 4.3.1 discusses radix choice,
///   where the radix b = 2^(8*SIZEOF_WORD) is chosen as a power of 2 matching
///   the machine word. This maximizes the use of hardware ALU carry flags.
// ---------------------------------------------------------------------------
#if SIZEOF_WORD == 1
using word    = uint8_t;
using sword   = int8_t;
using dword   = uint16_t;
using sdword  = int16_t;
#elif SIZEOF_WORD == 2
using word    = uint16_t;
using sword   = int16_t;
using dword   = uint32_t;
using sdword  = int32_t;
#elif SIZEOF_WORD == 4
using word    = uint32_t;
using sword   = int32_t;
using dword   = uint64_t;
using sdword  = int64_t;
#elif SIZEOF_WORD == 8 && __GNUC__ && __x86_64__
using word    = uint64_t;
using sword   = int64_t;
using dword   = __uint128_t;
using sdword  = __int128_t;
#elif __GNUC__ && __x86_64__
using word    = uint64_t;
using sword   = int64_t;
using dword   = __uint128_t;
using sdword  = __int128_t;
#ifdef SIZEOF_WORD
#undef SIZEOF_WORD
#endif
#define SIZEOF_WORD 8
#elif _MSC_VER && _M_IX86
using word    = uint32_t;
using sword   = int32_t;
using dword   = uint64_t;
using sdword  = int64_t;
#ifdef SIZEOF_WORD
#undef SIZEOF_WORD
#endif
#define SIZEOF_WORD 4
#elif _MSC_VER && _M_X64
using word    = uint32_t;
using sword   = int32_t;
using dword   = uint64_t;
using sdword  = int64_t;
#ifdef SIZEOF_WORD
#undef SIZEOF_WORD
#endif
#define SIZEOF_WORD 4
#else
using word    = uint32_t;
using sword   = int32_t;
using dword   = uint64_t;
using sdword  = int64_t;
#ifdef SIZEOF_WORD
#undef SIZEOF_WORD
#endif
#define SIZEOF_WORD 4
#endif
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Maximum-width unsigned integer type for accumulator operations
///
/// @details When SIZEOF_WORD < 8 (i.e., dword fits in uint64_t), the
///   maximum accumulator width is uintmax_t (typically uint64_t). When
///   SIZEOF_WORD == 8, the dword is __uint128_t, and umaxword_t mirrors
///   this for operations that need the widest available accumulator.
// ---------------------------------------------------------------------------
#if SIZEOF_WORD < 8
using umaxword_t = uintmax_t;
#else
using umaxword_t = __uint128_t;
#endif
//------------------------------------------------------------------------------
class integer;
class numeric;
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Core arbitrary-precision integer data structure
///
/// @details nn_integer_data provides the fundamental storage and arithmetic
///   primitives for the Numeric library. Key design decisions:
///
///   **Sign-magnitude representation**:
///   The integer value is stored as a sequence of `word` units in `data_[]`,
///   with a dedicated sign word at `data_[length_]`. The sign word is 0 for
///   positive values and all-ones (~0) for negative values (two's complement
///   sign extension at the word level). This is equivalent to the "signed
///   magnitude with separate sign" approach discussed in Knuth TAOCP Vol. 2,
///   Section 4.3.1, but using the high word for sign rather than a dedicated
///   sign field.
///
///   **Memory layout**:
///   The struct's declared `data_[4]` (or variant depending on SIZEOF_WORD)
///   is merely a placeholder. Actual allocations use `get_size()` which
///   computes `offsetof(nn_integer_data, data_) + (length + 2) * sizeof(word)`,
///   allocating the struct via the TLSF allocator with additional trailing
///   words beyond the declared array size. This is a well-known flexible
///   array member pattern in C++ (akin to C's FAM) using allocator-level
///   over-allocation.
///
///   **Reference counting**:
///   `ref_count_` enables shared ownership semantics without deep copies.
///   The `add_ref()` / `release()` pair follows the COM-style intrusive
///   reference counting pattern. This is crucial for the integer class:
///   assignment and parameter passing are O(1) pointer copies rather than
///   O(n) data copies.
///
///   **dummy_ field**:
///   Ensures `data_` has at least some padding for shift operations that
///   read one word past the declared array boundary. On 64-bit builds
///   (SIZEOF_WORD == 8), `dummy_` is a single word; on smaller word sizes
///   it is a uintptr_t sized to provide equivalent overflow protection.
///
///   References:
///     - Knuth, TAOCP Vol. 2, Section 4.3.1 (signed-magnitude arithmetic)
///     - ISO C++ Standard [expr.reinterpret.cast] for pointer aliasing
///     - TLSF: Two-Level Segregated Fit memory allocator (M. Masmano et al.)
// ---------------------------------------------------------------------------
struct nn_integer_data {
    protected:
        // -------------------------------------------------------------------
        /// @brief Constructs an empty nn_integer_data with given ref count and length
        ///
        /// @param ref_count Initial reference count (typically 1)
        /// @param length Number of word-sized digits (excluding sign word)
        ///
        /// @note This constructor is protected; only nn_new() and the
        ///   umaxword_t constructor may create instances.
        // -------------------------------------------------------------------
        nn_integer_data(intptr_t ref_count, uintptr_t length)
            : ref_count_(ref_count), length_(length), dummy_(0) {}

    public:
        ~nn_integer_data() {}

        // -------------------------------------------------------------------
        /// @brief Constructs from a umaxword_t value (small integer constant)
        ///
        /// @details Extracts individual words from the umaxword_t value by
        ///   shifting and masking. The sign word is implicitly 0 (positive).
        ///
        /// @param data Source value to initialize from
        // -------------------------------------------------------------------
        nn_integer_data(umaxword_t data)
            : ref_count_(1), length_(sizeof(data) / sizeof(data_[0])), dummy_(0) {
            for( size_t i = 0; i < sizeof(data) / sizeof(data_[0]); i++ ) {
                data_[i] = static_cast<word>(data);
                data >>= sizeof(data_[0]) * CHAR_BIT;
            }
        }

        mutable uintptr_t ref_count_;
        mutable uintptr_t length_;
#if SIZEOF_WORD == 1
        mutable uintptr_t dummy_;       // for shift down overflow bits
        mutable word data_[16];
#elif SIZEOF_WORD == 2
        mutable uintptr_t dummy_;       // for shift down overflow bits
        mutable word data_[8];
#elif SIZEOF_WORD == 4
        mutable uintptr_t dummy_;
        mutable word data_[4];
#elif SIZEOF_WORD == 8
        mutable word dummy_;
        mutable word data_[4];
#else
#error Invalid macro SIZEOF_WORD
#endif

        /// @brief Pointer type for nn_integer_data (convention: nn_integer is a pointer)
        using nn_integer = nn_integer_data *;

        // -------------------------------------------------------------------
        /// @brief Increments reference count and returns a non-const this pointer
        ///
        /// @details Used when creating an additional reference to the same
        ///   underlying data, avoiding a deep copy.
        ///
        /// @return Non-const nn_integer pointing to this object
        // -------------------------------------------------------------------
        nn_integer add_ref() const {
            ref_count_++;
            return const_cast<nn_integer>(this);
        }

        // -------------------------------------------------------------------
        /// @brief Computes the allocated size for an integer with the given length
        ///
        /// @details The formula accounts for:
        ///   - The base struct size up to `data_` (via offsetof)
        ///   - `length` word-sized digits
        ///   - 1 sign word at data_[length]
        ///   - 1 extra word for shift-overflow safety
        ///   Total: offsetof + (length + 2) * sizeof(word)
        ///
        /// @param length Requested number of digits (excluding sign and guard)
        /// @return Size in bytes suitable for TLSF::malloc
        ///
        /// @complexity O(1)
        // -------------------------------------------------------------------
        static constexpr size_t get_size(uintptr_t length) {
            return offsetof(nn_integer_data, data_) + (length + 2) * sizeof(data_[0]);
        }

        // -------------------------------------------------------------------
        /// @brief Returns the singleton TLSF allocator for integer data
        ///
        /// @details Uses the Meyer singleton pattern (function-local static)
        ///   to provide thread-safe, lazy initialization of the TLSF allocator.
        ///   All nn_integer_data allocations share this single allocator.
        ///
        /// @return Reference to the global TLSF allocator instance
        // -------------------------------------------------------------------
        static tlsf::TLSF_Impl & static_allocator() {
            static tlsf::TLSF_Impl allocator; // singleton
            return allocator;
        }

        // -------------------------------------------------------------------
        /// @brief Allocates and constructs an nn_integer_data with the given digit length
        ///
        /// @details Uses placement new on TLSF-allocated memory. The returned
        ///   object has ref_count_ = 1 and length_ = length.
        ///
        /// @param length Number of digits required
        /// @return Pointer to a newly constructed nn_integer_data, or nullptr
        ///   if allocation fails
        ///
        /// @complexity O(1) allocation, O(1) construction
        // -------------------------------------------------------------------
        static nn_integer nn_new(uintptr_t length) {
            nn_integer p = static_cast<nn_integer>(static_allocator().malloc(get_size(length)));

            if( p != nullptr )
                new (p) nn_integer_data(1, length);

            return p;
        }

        // -------------------------------------------------------------------
        /// @brief Decrements the reference count and frees memory if it reaches zero
        ///
        /// @details Following COM-style intrusive reference counting:
        ///   - If ref_count_ drops to 0, the destructor is called explicitly
        ///     and the memory is returned to the TLSF allocator.
        ///   - Otherwise, the object remains alive for other references.
        ///
        /// @complexity O(1) amortized
        // -------------------------------------------------------------------
        void release() {
            if( --ref_count_ == 0 ) {
                this->~nn_integer_data();
                static_allocator().free(this);
            }
        }

        // -------------------------------------------------------------------
        /// @brief Extracts the sign of the integer
        ///
        /// @details The sign word at data_[length_] is 0 for non-negative
        ///   and all-ones (~0) for negative. Right-shifting by
        ///   (sizeof(word)*CHAR_BIT - 1) produces -1 for negative, 0 for
        ///   non-negative. This is the arithmetic (sign-extending) right
        ///   shift applied to the sign word interpreted as a signed integer.
        ///
        /// @return 0 if non-negative, -1 (all-ones) if negative
        ///
        /// @complexity O(1)
        // -------------------------------------------------------------------
        sword isign() const {
            return static_cast<sword>(data_[length_]) >> (sizeof(word) * CHAR_BIT - 1);
        }

        // -------------------------------------------------------------------
        /// @brief Returns the most significant digit (highest word in magnitude)
        ///
        /// @details Assumes canonical form: the sign word is the only extra
        ///   word beyond the magnitude. The most significant magnitude word
        ///   is at index length_ - 1.
        ///
        /// @return The highest word of the magnitude
        ///
        /// @complexity O(1)
        // -------------------------------------------------------------------
        word high_word() const {
            return data_[length_ - 1];
        }

        // -------------------------------------------------------------------
        /// @brief Normalizes the representation by removing redundant sign words
        ///
        /// @details After arithmetic operations, the result may have a trailing
        ///   sign word that duplicates the previous word (e.g., 0x000...00 or
        ///   0xFFF...FF). This function shrinks length_ until the top two words
        ///   differ, ensuring canonical form. Normalization is O(n) worst-case
        ///   but O(1) amortized since most operations add at most one redundant
        ///   word.
        ///
        ///   Reference: Knuth TAOCP Vol. 2, Section 4.3.1 discusses the need
        ///   for canonicalization after addition/subtraction to maintain the
        ///   invariant that the leading digit is non-redundant.
        ///
        /// @complexity O(n) worst-case, O(1) typical
        // -------------------------------------------------------------------
        void normalize() {
            while( length_ > 1 && data_[length_] == data_[length_ - 1] )
                length_--;
        }

        // -------------------------------------------------------------------
        /// @brief Adds two big integers where `length_ <= p1->length_`
        ///   (this is the shorter operand, p1 is the longer)
        ///
        /// @details Implements schoolbook addition for the case where the
        ///   right operand (p1) is at least as long as the left operand (this).
        ///   The algorithm proceeds in three phases:
        ///
        ///   Phase 1 (common overlap): Iterate through `length_` words,
        ///     adding this->data_[i] + p1->data_[i] + carry.
        ///   Phase 2 (remainder of p1): Iterate through the remaining words of
        ///     p1, adding s0 (this's sign word) + p1->data_[i] + carry.
        ///   Phase 3 (final carry): Add s0 + s1 with final carry to produce
        ///     the result's two highest words.
        ///
        ///   The result length is initially p1->length_, then incremented
        ///   and normalized in case of overflow into the sign position.
        ///
        ///   Reference: Knuth, TAOCP Vol. 2, Section 4.3.1, Algorithm A
        ///   (addition of non-negative integers), extended here for signed
        ///   two's complement via sign-word inclusion in all phases.
        ///
        /// @param p1 The addend (must satisfy p1->length_ >= length_)
        /// @return A new nn_integer containing the sum
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(p1->length_) time, O(p1->length_) space
        // -------------------------------------------------------------------
        nn_integer addm(const nn_integer p1) const {
            nn_integer result = nn_new(p1->length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vsum = _mm256_add_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vsum);

                for( int i = 0; i < 4; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vsum = _mm_add_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vsum);

                for( int i = 0; i < 2; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) + *d1++ + cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            e = result->data_ + p1->length_;

            while( r < e ) {
                q = static_cast<dword>(s0) + *d1++ + cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) + s1 + cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) + s1 + cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Adds two big integers of equal length
        ///
        /// @details Specialized case where `length_ == p1->length_`. Unlike
        ///   addm/addp, there is no remainder phase — both operands are
        ///   exhausted simultaneously. The algorithm runs a single loop over
        ///   all words, then processes the two sign words with carry.
        ///
        /// @param p1 The addend (must satisfy p1->length_ == length_)
        /// @return A new nn_integer containing the sum
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(length_) time, O(length_) space
        // -------------------------------------------------------------------
        nn_integer addz(const nn_integer p1) const {
            nn_integer result = nn_new(length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + result->length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vsum = _mm256_add_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vsum);

                for( int i = 0; i < 4; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vsum = _mm_add_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vsum);

                for( int i = 0; i < 2; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) + *d1++ + cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) + s1 + cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) + s1 + cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Adds two big integers where `length_ >= p1->length_`
        ///   (this is the longer operand, p1 is the shorter)
        ///
        /// @details The left operand (this) is longer or equal. Phase 1 adds
        ///   the overlapping p1->length_ words. Phase 2 propagates carries
        ///   through this's remaining words using s1 (p1's sign word).
        ///   Phase 3 computes the final carry with s0 + s1.
        ///
        /// @param p1 The addend (must satisfy p1->length_ <= length_)
        /// @return A new nn_integer containing the sum
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(length_) time, O(length_) space
        // -------------------------------------------------------------------
        nn_integer addp(const nn_integer p1) const {
            nn_integer result = nn_new(length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + p1->length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vsum = _mm256_add_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vsum);

                for( int i = 0; i < 4; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vsum = _mm_add_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vsum);

                for( int i = 0; i < 2; ++i ) {
                    word carry0 = r[i] < d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) + cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | carry0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) + *d1++ + cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            e = result->data_ + length_;

            while( r < e ) {
                q = static_cast<dword>(*d0++) + s1 + cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) + s1 + cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) + s1 + cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Subtracts p1 from this, where `length_ <= p1->length_`
        ///   (the subtrahend is longer than the minuend)
        ///
        /// @details Implements schoolbook subtraction with borrow propagation.
        ///   The three-phase structure mirrors addm but uses subtraction:
        ///
        ///   Phase 1: this->data_[i] - p1->data_[i] - borrow (overlap)
        ///   Phase 2: s0 - p1->data_[i] - borrow (this exhausted)
        ///   Phase 3: s0 - s1 - final borrow
        ///
        ///   The borrow is extracted from bit sizeof(word)*CHAR_BIT of the
        ///   dword result (1 if the subtraction underflowed).
        ///
        ///   Reference: Knuth, TAOCP Vol. 2, Section 4.3.1, Algorithm S
        ///   (subtraction), extended for signed two's complement.
        ///
        /// @param p1 The subtrahend (must satisfy p1->length_ >= length_)
        /// @return A new nn_integer containing (this - p1)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(p1->length_) time, O(p1->length_) space
        // -------------------------------------------------------------------
        nn_integer subm(const nn_integer p1) const {
            nn_integer result = nn_new(p1->length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vdiff = _mm256_sub_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vdiff);

                for( int i = 0; i < 4; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vdiff = _mm_sub_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vdiff);

                for( int i = 0; i < 2; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) - *d1++ - cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            e = result->data_ + p1->length_;

            while( r < e ) {
                q = static_cast<dword>(s0) - *d1++ - cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) - s1 - cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) - s1 - cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Subtracts p1 from this, where both have equal length
        ///
        /// @details Specialized case for equal-length operands. A single loop
        ///   processes all digits, followed by the two-word sign computation.
        ///
        /// @param p1 The subtrahend (must satisfy p1->length_ == length_)
        /// @return A new nn_integer containing (this - p1)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(length_) time, O(length_) space
        // -------------------------------------------------------------------
        nn_integer subz(const nn_integer p1) const {
            nn_integer result = nn_new(length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + result->length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vdiff = _mm256_sub_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vdiff);

                for( int i = 0; i < 4; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vdiff = _mm_sub_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vdiff);

                for( int i = 0; i < 2; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) - *d1++ - cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) - s1 - cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) - s1 - cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Subtracts p1 from this, where `length_ >= p1->length_`
        ///   (the minuend is longer than the subtrahend)
        ///
        /// @details The left operand (this) is longer. Phase 1 processes the
        ///   overlap with p1. Phase 2 borrow-propagates through this's
        ///   remaining words using s1 (p1's sign word). Phase 3 computes
        ///   s0 - s1 with final borrow.
        ///
        /// @param p1 The subtrahend (must satisfy p1->length_ <= length_)
        /// @return A new nn_integer containing (this - p1)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(length_) time, O(length_) space
        // -------------------------------------------------------------------
        nn_integer subp(const nn_integer p1) const {
            nn_integer result = nn_new(length_);

            dword q = 0;
            word cf = 0;
            word * r = result->data_;
            const word * e = r + p1->length_;
            const word * d0 = data_;
            const word * d1 = p1->data_;
            word s0 = d0[length_];
            word s1 = d1[p1->length_];

#ifdef NN_HAVE_AVX2
            while( r + 4 <= e ) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d0));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(d1));
                __m256i vdiff = _mm256_sub_epi64(va, vb);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(r), vdiff);

                for( int i = 0; i < 4; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 4; d1 += 4; r += 4;
            }
#elif defined(NN_HAVE_AVX)
            while( r + 2 <= e ) {
                __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d0));
                __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(d1));
                __m128i vdiff = _mm_sub_epi64(va, vb);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(r), vdiff);

                for( int i = 0; i < 2; ++i ) {
                    word borrow0 = r[i] > d0[i] ? 1 : 0;
                    dword s = static_cast<dword>(r[i]) - cf;
                    cf = ((s >> (sizeof(word) * CHAR_BIT)) != 0) | borrow0;
                    r[i] = static_cast<word>(s);
                }

                d0 += 2; d1 += 2; r += 2;
            }
#endif

            while( r < e ) {
                q = static_cast<dword>(*d0++) - *d1++ - cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            e = result->data_ + length_;

            while( r < e ) {
                q = static_cast<dword>(*d0++) - s1 - cf;
                cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
                *r++ = static_cast<word>(q);
            }

            q = static_cast<dword>(s0) - s1 - cf;
            *r++ = static_cast<word>(q);
            cf = static_cast<word>((q >> sizeof(word) * CHAR_BIT) & 1);
            *r = static_cast<dword>(s0) - s1 - cf;

            // If overflow then grow size
            result->length_++;
            result->normalize();

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Dispatches to the correct addition variant based on length
        ///   comparison between this and v
        ///
        /// @details Uses the three-way comparison `length_ - v->length_` to
        ///   select the optimal addition strategy:
        ///     - length_ > v->length_: addp (this is longer)
        ///     - length_ < v->length_: addm (v is longer)
        ///     - length_ == v->length_: addz (equal lengths)
        ///
        ///   The difference is stored in a signed intptr_t to handle the
        ///   case where v is longer (producing a negative value from unsigned
        ///   subtraction). This is well-defined: uintptr_t subtraction wraps
        ///   modulo 2^N, and the resulting bit pattern reinterpreted as
        ///   intptr_t (two's complement) gives the correct sign.
        ///
        /// @param v The addend
        /// @return A new nn_integer containing (this + v)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(max(length_, v->length_)) time, same for space
        // -------------------------------------------------------------------
        nn_integer iadd(const nn_integer v) const {
            intptr_t c = static_cast<intptr_t>(length_ - v->length_);

            if( c > 0 )
                return addp(v);

            if( c < 0 )
                return addm(v);

            return addz(v);
        }

        // -------------------------------------------------------------------
        /// @brief Dispatches to the correct subtraction variant based on
        ///   length comparison between this and v
        ///
        /// @details Same dispatch pattern as iadd, using the three length
        ///   cases to select subm/subz/subp.
        ///
        /// @param v The subtrahend
        /// @return A new nn_integer containing (this - v)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(max(length_, v->length_)) time, same for space
        // -------------------------------------------------------------------
        nn_integer isub(const nn_integer v) const {
            intptr_t c = static_cast<intptr_t>(length_ - v->length_);

            if( c > 0 )
                return subp(v);

            if( c < 0 )
                return subm(v);

            return subz(v);
        }

        // -------------------------------------------------------------------
        /// @brief Arithmetic left shift (declaration)
        ///
        /// @details Defined out-of-line due to #pragma GCC optimize("O3")
        ///   guard. See the implementation below for details.
        // -------------------------------------------------------------------
        nn_integer isal(uintptr_t bit_count) const;
        // -------------------------------------------------------------------
        /// @brief Arithmetic right shift (declaration)
        ///
        /// @details Defined out-of-line below. See the implementation for details.
        // -------------------------------------------------------------------
        nn_integer isar(uintptr_t bit_count) const;

        // -------------------------------------------------------------------
        /// @brief Tests whether the integer is zero
        ///
        /// @details Checks all data words from 0 to length_ (inclusive of
        ///   sign word) for non-zero values. Uses an early exit: if data_[0]
        ///   is non-zero, returns false immediately without scanning the rest.
        ///
        /// @return true if the integer is zero, false otherwise
        ///
        /// @complexity O(length_) worst-case, O(1) average
        // -------------------------------------------------------------------
        bool z() const {
            if( data_[0] != 0 )
                return false;

            for( uintptr_t i = length_; i > 0; i-- )
                if( data_[i] != 0 )
                    return false;

            return true;
        }

        // -------------------------------------------------------------------
        /// @brief Tests whether the integer is non-zero
        ///
        /// @details Linear scan through all words including sign word.
        ///   Returns true at the first non-zero word found.
        ///
        /// @return true if any word is non-zero, false if the integer is zero
        ///
        /// @complexity O(length_) worst-case, O(1) average (early exit)
        // -------------------------------------------------------------------
        bool nz() const {
            for( uintptr_t i = length_; i != static_cast<uintptr_t>(-1); i-- )
                if( data_[i] != 0 )
                    return true;

            return false;
        }

        // -------------------------------------------------------------------
        /// @brief Compares this integer with v1
        ///
        /// @details Computes the sign of (this - v1) by performing a full
        ///   subtraction and examining the result. This is O(n) in the
        ///   worst case but avoids branching on the sign of the operands.
        ///
        ///   Note: The commented-out code suggests an earlier direct sign
        ///   comparison was considered but replaced with the subtract-and-
        ///   check approach for correctness with mixed-sign operands.
        ///
        /// @param v1 The integer to compare against
        /// @return -1 if this < v1, 0 if this == v1, 1 if this > v1
        ///
        /// @complexity O(max(length_, v1->length_)) time
        /// @note The temporary result from isub must be released
        // -------------------------------------------------------------------
        intptr_t icompare(const nn_integer v1) const {
            intptr_t c = 0;

            if( c == 0 ) {
                nn_integer a = isub(v1);

                c = a->isign() < 0 ? -1 : a->z() ? 0 : 1;

                a->release();
            }

            return c;
        }

        // -------------------------------------------------------------------
        /// @brief Adds a word value to a specific digit position with carry
        ///   propagation
        ///
        /// @details Starting at the given index, adds `value` to data_[index].
        ///   If overflow occurs (value wraps past 2^wordsize), the carry is
        ///   propagated to successive positions. This is the inner loop of
        ///   the multiplication algorithm.
        ///
        /// @param index Starting digit position within data_[]
        /// @param value The word value to add (may be a partial product)
        ///
        /// @complexity O(length_ - index) in pathological case, O(1) typical
        // -------------------------------------------------------------------
        void accumulate_with_carry(uintptr_t index, word value) {
            dword q;

            while( value != 0 && index < length_ ) {
                q = static_cast<dword>(data_[index]) + value;
                data_[index] = static_cast<word>(q);
                value = static_cast<word>(q >> (sizeof(word) * CHAR_BIT));
                index++;
            }
        }

        // -------------------------------------------------------------------
        /// @brief Multiplies a single word m by all words of b, accumulating
        ///   into this at position i
        ///
        /// @details Classical double-and-add multiplication step. For each
        ///   digit j in b where b->data_[j] != 0, computes the dword product
        ///   m * b->data_[j] and accumulates the low and high halves into
        ///   positions i+j and i+j+1.
        ///
        ///   This is the inner loop of the O(n*m) schoolbook multiplication.
        ///
        /// @param b The multiplicand (multi-word integer)
        /// @param i Offset into this's data_ for accumulation
        /// @param m The multiplier word
        ///
        /// @complexity O(b->length_) time
        // -------------------------------------------------------------------
        void imul(const nn_integer b, uintptr_t i, word m) {
            for( uintptr_t j = 0; j < b->length_; j++ ) {
                if( b->data_[j] == 0 )
                    continue;

                dword q = static_cast<dword>(m) * b->data_[j];

                accumulate_with_carry(i + j, static_cast<word>(q));
                accumulate_with_carry(i + j + 1, static_cast<word>(q >> (sizeof(word) * CHAR_BIT)));
            }
        }

        // -------------------------------------------------------------------
        /// @brief Multiplies two big integers using schoolbook multiplication
        ///
        /// @details For each word a->data_[i], calls the single-word
        ///   multiplier imul(b, i, a->data_[i]) to add a->data_[i] * b
        ///   shifted by i word positions.
        ///
        ///   The commented-out code shows an equivalent direct nested loop
        ///   that was replaced by the factored single-word multiply for
        ///   readability and reuse.
        ///
        ///   Reference: Knuth, TAOCP Vol. 2, Section 4.3.1, Algorithm M
        ///   (multiplication of non-negative integers).
        ///
        /// @param a First operand
        /// @param b Second operand
        ///
        /// @complexity O(a->length_ * b->length_) time
        /// @pre this->data_ must be zero-initialized and have sufficient
        ///   length to hold the product (size >= a->length_ + b->length_)
        // -------------------------------------------------------------------
        void imul(const nn_integer a, const nn_integer b) {
            for( uintptr_t i = 0; i < a->length_; i++ ) {
                if( a->data_[i] == 0 )
                    continue;

                imul(b, i, a->data_[i]);
            }
        }

        // -------------------------------------------------------------------
        /// @brief Bitwise NOT (one's complement) with sign extension
        ///
        /// @details Produces ~this by complementing every data word including
        ///   the sign word. Since two's complement negation is (~x + 1),
        ///   this operation is one step short of arithmetic negation.
        ///
        /// @return A new nn_integer containing the bitwise NOT
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(length_) time, O(length_) space
        // -------------------------------------------------------------------
        nn_integer inot() const {
            nn_integer result = nn_new(length_);

            for( uintptr_t i = length_ + 1; i != static_cast<uintptr_t>(-1); i-- )
                result->data_[i] = ~data_[i];

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Bitwise AND with sign extension
        ///
        /// @details Computes this & v. When operand lengths differ, the
        ///   shorter operand's sign word is extended to match the longer
        ///   operand's length. This ensures correct two's complement
        ///   semantics: for example, ANDing a positive value (sign = 0)
        ///   with any value produces the positive value limited to the
        ///   shorter length's range.
        ///
        /// @param v The second operand
        /// @return A new nn_integer containing (this & v)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(max(length_, v->length_)) time and space
        // -------------------------------------------------------------------
        nn_integer iand(const nn_integer v) const {
            nn_integer result;
            uintptr_t i;

            if( length_ > v->length_ ) {
                result = nn_new(length_);

                for( i = 0; i < v->length_; i++ )
                    result->data_[i] = data_[i] & v->data_[i];

                while( i < length_ + 2 ) {
                    result->data_[i] = data_[i] & v->data_[v->length_];
                    i++;
                }
            }
            else if( length_ < v->length_ ) {
                result = nn_new(v->length_);

                for( i = 0; i < length_; i++ )
                    result->data_[i] = data_[i] & v->data_[i];

                while( i < v->length_ + 2 ) {
                    result->data_[i] = data_[length_] & v->data_[i];
                    i++;
                }
            }
            else {
                result = nn_new(length_);

                for( i = 0; i < length_ + 2; i++ )
                    result->data_[i] = data_[i] & v->data_[i];
            }

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Bitwise OR with sign extension
        ///
        /// @details Computes this | v with the same sign-extension strategy
        ///   as iand. The sign word of the shorter operand is extended to
        ///   fill the difference in length.
        ///
        /// @param v The second operand
        /// @return A new nn_integer containing (this | v)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(max(length_, v->length_)) time and space
        // -------------------------------------------------------------------
        nn_integer ior(const nn_integer v) const {
            nn_integer result;
            uintptr_t i;

            if( length_ > v->length_ ) {
                result = nn_new(length_);

                for( i = 0; i < v->length_; i++ )
                    result->data_[i] = data_[i] | v->data_[i];

                while( i < length_ + 2 ) {
                    result->data_[i] = data_[i] | v->data_[v->length_];
                    i++;
                }
            }
            else if( length_ < v->length_ ) {
                result = nn_new(v->length_);

                for( i = 0; i < length_; i++ )
                    result->data_[i] = data_[i] | v->data_[i];

                while( i < v->length_ + 2 ) {
                    result->data_[i] = data_[length_] | v->data_[i];
                    i++;
                }
            }
            else {
                result = nn_new(length_);

                for( i = 0; i < length_ + 2; i++ )
                    result->data_[i] = data_[i] | v->data_[i];
            }

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Bitwise XOR with sign extension
        ///
        /// @details Computes this ^ v with the same sign-extension strategy
        ///   as iand/ior. XOR with sign-extension produces the correct
        ///   two's complement result for mixed-length operands.
        ///
        /// @param v The second operand
        /// @return A new nn_integer containing (this ^ v)
        /// @note Caller must release() the returned object when done
        ///
        /// @complexity O(max(length_, v->length_)) time and space
        // -------------------------------------------------------------------
        nn_integer ixor(const nn_integer v) const {
            nn_integer result;
            uintptr_t i;

            if( length_ > v->length_ ) {
                result = nn_new(length_);

                for( i = 0; i < v->length_; i++ )
                    result->data_[i] = data_[i] ^ v->data_[i];

                while( i < length_ + 2 ) {
                    result->data_[i] = data_[i] ^ v->data_[v->length_];
                    i++;
                }
            }
            else if( length_ < v->length_ ) {
                result = nn_new(v->length_);

                for( i = 0; i < length_; i++ )
                    result->data_[i] = data_[i] ^ v->data_[i];

                while( i < v->length_ + 2 ) {
                    result->data_[i] = data_[length_] ^ v->data_[i];
                    i++;
                }
            }
            else {
                result = nn_new(length_);

                for( i = 0; i < length_ + 2; i++ )
                    result->data_[i] = data_[i] ^ v->data_[i];
            }

            return result;
        }

        // -------------------------------------------------------------------
        /// @brief Reads the i-th bit of the integer
        ///
        /// @details Computes data_[i / wordsize] >> (i % wordsize) & 1.
        ///   The assertion checks that i is within the allocated range.
        ///
        /// @param i Bit index (0 = least significant bit)
        /// @return 0 or 1
        ///
        /// @complexity O(1)
        // -------------------------------------------------------------------
        uint8_t ibit(uintptr_t i) const {
            assert( i < length_ * sizeof(word) * CHAR_BIT );
            return static_cast<uint8_t>(data_[i / (sizeof(word) * CHAR_BIT)] >> (i & (sizeof(word) * CHAR_BIT - 1)) & 1);
        }

        // -------------------------------------------------------------------
        /// @brief Sets the i-th bit of the integer (via OR)
        ///
        /// @details ORs the value v (by default 1) into bit position i.
        ///   The assertion checks bounds.
        ///
        /// @tparam T Type of the bit value (default: word)
        /// @param i Bit index
        /// @param v Value to OR into the bit position (default: 1)
        ///
        /// @complexity O(1)
        // -------------------------------------------------------------------
        template <typename T = word>
        void sbit(uintptr_t i, const T v = 1) const {
            assert( i < length_ * sizeof(word) * CHAR_BIT );
            data_[i / (sizeof(word) * CHAR_BIT)] |= static_cast<word>(v) << (i & (sizeof(word) * CHAR_BIT - 1));
        }

        // -------------------------------------------------------------------
        /// @brief Counts trailing zero bits across all words
        ///
        /// @details Iterates through words from LSB upward. For each word,
        ///   uses platform-specific intrinsics if available:
        ///
        ///     - GCC/Clang (SIZEOF_WORD <= 4): __builtin_ctz(x)
        ///     - GCC/Clang (SIZEOF_WORD == 8): __builtin_ctzll(x)
        ///     - Fallback: binary-search divide-and-conquer using 32/16/8/4/2
        ///       bit masks. This is equivalent to the de Bruijn sequence
        ///       approach described in Warren's "Hacker's Delight".
        ///
        ///   The search continues until the first non-zero word is found,
        ///   adding word-width accumulators for preceding zero words.
        ///
        ///   Used by binary GCD (Stein's algorithm) to extract the
        ///   factor 2^k common to both operands.
        ///
        /// @return Number of consecutive zero bits starting from LSB
        ///
        /// @complexity O(n) worst-case (all words zero), O(n/wordsize + 1)
        ///   average, O(1) with builtins
        ///
        ///   Reference: Stein, J. "Computational problems associated with
        ///   Racah algebra" (1967) — binary GCD algorithm.
        ///   Warren, H.S. "Hacker's Delight" (2nd ed.), Section 5-4.
        // -------------------------------------------------------------------
        uintptr_t count_trailing_zeros() const {
            uintptr_t n;
            word * p = data_;

            for( n = 0; n < length_ * sizeof(word) * CHAR_BIT; n += sizeof(word) * CHAR_BIT, p++ ) {
                word x = *p;

                if( x != 0 ) {
#if __GNUC__ && SIZEOF_WORD <= 4
                    n += __builtin_ctz(x);
#elif __GNUC__ && SIZEOF_WORD == 8
                    n += __builtin_ctzll(x);
#else
                    uintptr_t a = 1;
#if SIZEOF_WORD >= 8
                    if( (x & 0xFFFFFFFF) == 0 ) { a = a + 32; x = x >> 32; }
#endif
#if SIZEOF_WORD >= 4
                    if( (x & 0x0000FFFF) == 0 ) { a = a + 16; x = x >> 16; }
#endif
#if SIZEOF_WORD >= 2
                    if( (x & 0x000000FF) == 0 ) { a = a +  8; x = x >>  8; }
#endif
                    if( (x & 0x0000000F) == 0 ) { a = a +  4; x = x >>  4; }
                    if( (x & 0x00000003) == 0 ) { a = a +  2; x = x >>  2; }

                    n += a - (x & 1);
#endif
                    break;
                }
            }

            return n;
        }
};
//------------------------------------------------------------------------------
/// @brief Public type alias for nn_integer_data pointer (C++17 `using` alias)
using nn_integer = nn_integer_data *;
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Converts a single hex character to its 4-bit integer value
///
/// @details Processes '0'-'9', 'A'-'F', 'a'-'f' using subtract-based
///   conversion. Invalid characters cause an exception via the ternary's
///   throw arm. The function is constexpr, enabling compile-time evaluation
///   of hex string literals.
///
/// @param input Hex character (ASCII)
/// @return 4-bit unsigned integer value (0-15)
/// @throws std::invalid_argument if input is not a valid hex digit
///
/// @complexity O(1)
// ---------------------------------------------------------------------------
constexpr uint8_t hex_char2int(char input)
{
    return
    ((input >= 'a') && (input <= 'f'))
    ? static_cast<uint8_t>(input - 'a')
    : ((input >= 'A') && (input <= 'F'))
    ? static_cast<uint8_t>(input - 'A')
    : ((input >= '0') && (input <= '9'))
    ? static_cast<uint8_t>(input - '0')
    : throw std::invalid_argument{"hex_char2int: invalid hex character"};
}
// ---------------------------------------------------------------------------
/// @brief Converts two hex characters into a single byte
///
/// @details Combines `hex_char2int(high) << 4 | hex_char2int(low)` to
///   produce a byte from a pair of hex nibbles. The high nibble character
///   is the first character in reading order.
///
/// @param high Upper nibble hex character
/// @param low Lower nibble hex character
/// @return Byte value (0-255) formed from the two hex digits
///
/// @complexity O(1)
// ---------------------------------------------------------------------------
constexpr uint8_t hex_char(char high, char low)
{
    return static_cast<uint8_t>((hex_char2int(high) << 4) | (hex_char2int(low)));
}
// ---------------------------------------------------------------------------
/// @brief Adapter that converts a hex string literal to an array of bytes
///   at compile time
///
/// @details Uses std::index_sequence (C++14) to unpack pairs of characters
///   from the input string, calling hex_char on each pair and collecting
///   the results in a uniform initialization list for T's constructor.
///
/// @tparam T The destination type (typically std::array<uint8_t, N>)
/// @tparam length Length of the input character array (including null)
/// @tparam index std::index_sequence of indices for compile-time unpacking
/// @param input Hex string literal
/// @return T initialized with the decoded bytes
///
/// @complexity O(length), evaluated at compile time
// ---------------------------------------------------------------------------
template <typename T, size_t length, size_t ... index>
constexpr T hex_string(const char (&input)[length], const std::index_sequence<index...>&)
{
    return T{hex_char(input[index * 2], input[index * 2 + 1])...};
}
// ---------------------------------------------------------------------------
/// @brief Entry function for compile-time hex string to byte array conversion
///
/// @details Creates a std::index_sequence of half the input length and
///   forwards to the implementation that unpacks character pairs.
///
/// @tparam T The destination type
/// @tparam length Length of input (including null terminator)
/// @param input Hex string literal
/// @return T initialized with decoded bytes
///
/// @complexity O(length), evaluated at compile time
// ---------------------------------------------------------------------------
template <typename T, std::size_t length>
constexpr T hex_string(const char (&input)[length])
{
    return hex_string<T>(input, std::make_index_sequence<(length / 2)>{});
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Compile-time computation of the maximum value for a given unsigned
///   integer type (MSVC variant)
///
/// @details Uses repeated multiplication by 10 until overflow occurs.
///   The multiplication by 10 is decomposed as (m << 3) + (m << 1) which
///   prevents the compiler from issuing a warning about overflow in a
///   constexpr context. When the result wraps past the type's maximum,
///   the condition n > m fails, and the previous m is returned.
///
/// @tparam T The unsigned integer type
/// @param m Seed value (default: 1)
/// @return Maximum representable value for type T, rounded down to the
///   nearest multiple of 10
///
/// @complexity O(log_{10}(2^{sizeof(T)*8})) iterations
///
/// @note This is functionally equivalent to
///   std::numeric_limits<T>::max() - std::numeric_limits<T>::max() % 10.
///   The magic division by 10 is avoided in favor of portable multiplication.
// ---------------------------------------------------------------------------
#if _MSC_VER
template <typename T> constexpr const T uint_max(const T m = static_cast<T>(1))
{
    return (m << 3) + (m << 1) > m ? uint_max<T>((m << 3) + (m << 1)) : m;
}
#else
// ---------------------------------------------------------------------------
/// @brief Compile-time computation of the maximum value for a given unsigned
///   integer type (GCC/Clang variant)
///
/// @details Iteratively multiplies by 10 until overflow is detected
///   (when n <= m after multiplication). Returns the largest value that
///   does not overflow. Uses an infinite for(;;) loop with overflow
///   break condition, which is constexpr-valid in C++14 and later.
///
/// @tparam T The unsigned integer type
/// @return Maximum representable value for type T
///
/// @complexity O(log_{10}(2^{sizeof(T)*8})) iterations
// ---------------------------------------------------------------------------
template <typename T> constexpr const T uint_max()
{
    T m = 1u;

    for( ;; ) {
        T n = static_cast<T>(m * 10u);
        if( n <= m ) // overflow
            break;
        m = n;
    }

    return m;
}
#endif
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @name Global integer constants (pre-allocated singletons)
///
/// @details These extern const objects are immutable instances of
///   nn_integer_data for commonly used small integers (0, 1, 2, 4, 5, 6,
///   8, 10). They are allocated once at program startup and shared via
///   reference counting (each has ref_count_ initialized to 1, but as
///   const objects their ref_count_ is never modified).
///
///   Using pre-allocated singletons avoids repeated heap allocations for
///   small integers and enables pointer-equality checks (e.g., checking
///   if a value is zero by comparing its data pointer to &nn_izero).
///
///   nn_maxull is initialized at compile time via uint_max<umaxword_t>()
///   and represents the maximum value of umaxword_t rounded down to a
///   multiple of 10. It is used as an overflow sentinel for string-to-
///   integer conversion.
///
/// @{
// ---------------------------------------------------------------------------
extern const nn_integer_data nn_izero(0);
extern const nn_integer_data nn_ione(1);
extern const nn_integer_data nn_itwo(2);
extern const nn_integer_data nn_ifour(4);
extern const nn_integer_data nn_ifive(5);
extern const nn_integer_data nn_isix(6);
extern const nn_integer_data nn_ieight(8);
extern const nn_integer_data nn_iten(10);
// 1000000000u                         == 0x3B9ACA00
// 10000000000000000000u               == 0x8AC7230489E80000
// 100000000000000000000000000000000000000u == 0x4B3B4CA85A86C47A098A224000000000
#if SIZEOF_WORD == 1
extern const nn_integer_data nn_maxull(uint_max<umaxword_t>());
#elif SIZEOF_WORD == 2
extern const nn_integer_data nn_maxull(uint_max<umaxword_t>());
#elif SIZEOF_WORD == 4
extern const nn_integer_data nn_maxull(uint_max<umaxword_t>());
#elif SIZEOF_WORD == 8
extern const nn_integer_data nn_maxull(uint_max<umaxword_t>());
#endif
/// @}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Allocates a new nn_integer with the specified digit length
///
/// @details Wraps nn_integer_data::nn_new with mandatory null check.
///   If TLSF allocation fails, throws std::bad_alloc instead of returning
///   nullptr. This ensures all callers can rely on a valid pointer.
///
/// @param length Number of word-sized digits
/// @return A newly allocated nn_integer (ref_count_ == 1)
/// @throws std::bad_alloc if TLSF memory allocation fails
///
/// @complexity O(1) allocation
// ---------------------------------------------------------------------------
static inline nn_integer nn_new(uintptr_t length)
{
    nn_integer p = nn_integer_data::nn_new(length);
    if( p == nullptr )
        throw std::bad_alloc();
    return p;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Creates an nn_integer from a signed long long, with small-integer
///   fast path
///
/// @details Common small values (0, 1, 2, 4, 5, 6, 8, 10) return pre-allocated
///   singletons via add_ref(), avoiding heap allocation. Larger values are
///   allocated via nn_new, populated with memcpy, and their sign word is
///   set via arithmetic right shift of the source value.
///
/// @param v The source value
/// @return nn_integer representing v
///
/// @complexity O(1) for small values, O(sizeof(v)/sizeof(word)) otherwise
// ---------------------------------------------------------------------------
static inline nn_integer nn_init_illong(long long v)
{
    switch( v ) {
        case   0 : return nn_izero.add_ref();
        case   1 : return nn_ione.add_ref();
        case   2 : return nn_itwo.add_ref();
        case   4 : return nn_ifour.add_ref();
        case   5 : return nn_ifive.add_ref();
        case   6 : return nn_isix.add_ref();
        case   8 : return nn_ieight.add_ref();
        case  10 : return nn_iten.add_ref();
        default  : ;
    }

    nn_integer p = nn_new(sizeof(v) / sizeof(word));

    memcpy(p->data_, &v, sizeof(v));
    p->data_[p->length_] = static_cast<word>(v >> (sizeof(v) * CHAR_BIT - 1));
    p->normalize();

    return p;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Creates an nn_integer from an unsigned long long, with small-integer
///   fast path
///
/// @details Like nn_init_illong but for unsigned values. The sign word is
///   explicitly set to 0 since unsigned values are always non-negative.
///   Also includes overflow-case fast paths for nn_maxull at the appropriate
///   SIZEOF_WORD thresholds.
///
/// @param v The source value
/// @return nn_integer representing v
///
/// @complexity O(1) for small values, O(sizeof(v)/sizeof(word)) otherwise
// ---------------------------------------------------------------------------
static inline nn_integer nn_init_iullong(unsigned long long v)
{
    switch( v ) {
        case   0 : return nn_izero.add_ref();
        case   1 : return nn_ione.add_ref();
        case   2 : return nn_itwo.add_ref();
        case   4 : return nn_ifour.add_ref();
        case   5 : return nn_ifive.add_ref();
        case   6 : return nn_isix.add_ref();
        case   8 : return nn_ieight.add_ref();
        case  10 : return nn_iten.add_ref();
#if SIZEOF_WORD < 2
        case 1000000000u :
            return nn_maxull.add_ref();
#elif SIZEOF_WORD < 8
        case 10000000000000000000ull :
            return nn_maxull.add_ref();
#endif
        default  : ;
    }

    nn_integer p = nn_new(sizeof(v) / sizeof(word));

    memcpy(p->data_, &v, sizeof(v));
    p->data_[p->length_] = 0;
    p->normalize();

    return p;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Creates an nn_integer from a signed long (delegates to illong)
///
/// @param v The source value
/// @return nn_integer representing v
// ---------------------------------------------------------------------------
static inline nn_integer nn_init_ulong(long v)
{
    return nn_init_illong(v);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Creates an nn_integer from an unsigned long (delegates to iullong)
///
/// @param v The source value
/// @return nn_integer representing v
// ---------------------------------------------------------------------------
static inline nn_integer nn_init_iulong(unsigned long v)
{
    return nn_init_iullong(v);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @name Shift operation optimization pragma
///
/// @brief Apply -O3 optimization to isal for raw pointer performance
///
/// @details The arithmetic left shift (isal) uses overlapping dword-access
///   loops that benefit significantly from aggressive optimization.
///   The -O3 flag enables auto-vectorization and improved register
///   allocation for the loop body. This pragma is only applied in release
///   builds (NDEBUG defined), as -O3 can degrade debugger experience.
///
///   The right shift (isar) does not use this optimization guard in the
///   original code, as its loop direction and access pattern may interact
///   differently with the optimizer.
///
/// @{
// ---------------------------------------------------------------------------
#if __GNUC__ && NDEBUG
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Arithmetic left shift (x << k)
///
/// @details Shifts the integer left by `bit_count` bits, preserving the
///   sign bit and growing the representation as necessary. The algorithm
///   uses overlapping double-word accesses to propagate overflow bits
///   across word boundaries in a single pass right-to-left.
///
///   Procedure:
///   1. Handle zero-bit shift as identity (returns add_ref()).
///   2. Compute new bit size aligned to word boundary, allocate.
///   3. Zero-fill the lower `offset` words (bits shifted in from right).
///   4. If `shift == 0` (word-aligned), memcpy the data directly.
///   5. Otherwise, iterate from MSW to LSW, reading dwords at each word
///      position, shifting left, and writing the result through a dword
///      pointer. The overlapping dword accesses naturally propagate the
///      shifted-out bits into the next higher word.
///   6. Set the sign words from isign() and normalize.
///
///   The overlapping dword technique is a well-known optimization for
///   multi-precision shifts: each dword read covers two adjacent words
///   at addresses spaced by sizeof(word), so the high bits of position i
///   become the low bits of position i+1 after the shift.
///
///   Reference: Warren, H.S. "Hacker's Delight" (2nd ed.), Section 2-6
///   (shifts of multi-word integers).
///
/// @param bit_count Number of bit positions to shift left
/// @return A new nn_integer with value (this << bit_count)
/// @note Caller must release() the returned object when done
///
/// @complexity O(length_) time, O(length_ + bit_count/wordsize) space
// ---------------------------------------------------------------------------
inline nn_integer nn_integer_data::isal(uintptr_t bit_count) const
{
    if( bit_count == 0 )
        return add_ref();

    uintptr_t new_bit_size = length_ * sizeof(word) * CHAR_BIT + bit_count;
    new_bit_size += -static_cast<intptr_t>(new_bit_size) & (sizeof(word) * CHAR_BIT - 1);

    nn_integer result = nn_new(new_bit_size / (sizeof(word) * CHAR_BIT));

    uintptr_t offset = bit_count / (sizeof(word) * CHAR_BIT);
    uintptr_t shift = bit_count & (sizeof(word) * CHAR_BIT - 1);

    // Establish word pointers for the overlapping dword-access loop.
    // src_w points to the most significant word of the source (data_[length_-1]).
    // dst_w points to the word at result->length_ - 2, which is aligned to
    //   receive the most significant shifted pair.
    const word * src_w = data_ + length_ - 1;
    word * dst_w = result->data_ + result->length_ - 2;

    memset(result->data_, 0, offset * sizeof(word));

    if( shift == 0 ) {
        // Word-aligned shift: direct copy, no bit-level manipulation needed.
        memcpy(result->data_ + offset, data_, length_ * sizeof(word));
    }
    else {
        // Sub-word shift: iterate right-to-left, reading/writing overlapping
        // dwords via memcpy to propagate carry bits across word boundaries.
        // This is the safe alternative to union type-punning (strict aliasing).
        for( uintptr_t i = length_; i-- > 0; dst_w--, src_w-- ) {
            dword tmp;
            memcpy(&tmp, src_w, sizeof(dword));
            tmp <<= shift;
            memcpy(dst_w, &tmp, sizeof(dword));
        }
    }

    result->data_[result->length_ + 1] = result->data_[result->length_] = isign();
    result->normalize();

    return result;
}
//------------------------------------------------------------------------------
#if __GNUC__ && NDEBUG
#pragma GCC pop_options
#endif
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Arithmetic right shift (x >> k)
///
/// @details Shifts the integer right by `bit_count` bits, performing sign
///   extension. If bit_count >= total bit size, returns zero.
///
///   Procedure:
///   1. Handle zero-bit shift as identity.
///   2. If bit_count >= bit_size, return zero (all bits shifted out).
///   3. Compute new bit size aligned to word boundary, allocate.
///   4. Starting from the offset word in the source, iterate left-to-right,
///      reading overlapping dwords, shifting right, and writing to the
///      result. The overlapping access pattern ensures that bits shifted
///      out of the low end of word i are captured from the high end of
///      word i+1.
///   5. Set sign words from isign() and normalize.
///
///   The right-to-left direction means the MSW is processed first, and
///   subsequent iterations handle progressively lower words. The
///   overlapping dword read at position i reads words i and i+1, so the
///   right shift naturally pulls in bits from the next higher word.
///
///   Reference: Warren, H.S. "Hacker's Delight" (2nd ed.), Section 2-6.
///
/// @param bit_count Number of bit positions to shift right
/// @return A new nn_integer with value (this >> bit_count)
/// @note Caller must release() the returned object when done
///
/// @complexity O(result->length_) time, O(result->length_) space
// ---------------------------------------------------------------------------
inline nn_integer nn_integer_data::isar(uintptr_t bit_count) const
{
    if( bit_count == 0 )
        return add_ref();

    uintptr_t bit_size = length_ * sizeof(word) * CHAR_BIT;

    if( bit_count >= bit_size )
        return nn_izero.add_ref();

    intptr_t new_bit_size = static_cast<intptr_t>(bit_size - bit_count);
    new_bit_size += -static_cast<intptr_t>(new_bit_size) & (sizeof(word) * CHAR_BIT - 1);

    nn_integer result = nn_new(static_cast<uintptr_t>(new_bit_size) / (sizeof(word) * CHAR_BIT));

    uintptr_t offset = bit_count / (sizeof(word) * CHAR_BIT);
    uintptr_t shift = bit_count & (sizeof(word) * CHAR_BIT - 1);

    // src_w points to the first word in the source that contributes to
    // the result (after skipping `offset` whole words). dst_w points to
    // the beginning of the result buffer.
    const word * src_w = data_ + offset;
    word * dst_w = result->data_;

    if( shift == 0 ) {
        // Word-aligned shift: direct copy of the remaining words.
        memcpy(result->data_, data_ + offset, result->length_ * sizeof(word));
    }
    else {
        // Sub-word shift: iterate left-to-right, reading overlapping dwords
        // via memcpy to propagate borrowed bits across word boundaries.
        // This is the safe alternative to union type-punning (strict aliasing).
        for( uintptr_t i = result->length_; i-- > 0; dst_w++, src_w++ ) {
            dword tmp;
            memcpy(&tmp, src_w, sizeof(dword));
            tmp >>= shift;
            memcpy(dst_w, &tmp, sizeof(dword));
        }
    }

    result->dummy_ = 0;
    result->data_[result->length_ + 1] = result->data_[result->length_] = isign();
    result->normalize();

    return result;
}
//------------------------------------------------------------------------------
} // namespace NaturalNumbers
//------------------------------------------------------------------------------
#endif // ID_HPP_INCLUDED
//------------------------------------------------------------------------------
