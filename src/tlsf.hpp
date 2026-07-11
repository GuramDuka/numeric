/*-
 * The MIT License (MIT)
 *
 * Copyright (c) 2014, 2015, 2016, 2017 Guram Duka
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
#ifndef TLSF_HPP_INCLUDED
#define TLSF_HPP_INCLUDED
//------------------------------------------------------------------------------
#if __GNUC__
#ifndef EMSCRIPTEN
#include <bits/c++config.h>
#endif
//#include <ansidecl.h>
//#ifndef HAVE_LONG_DOUBLE
//#define HAVE_LONG_DOUBLE 1
//#endif
#endif
#if __GNUC__
#include <sys/param.h>
#endif
//------------------------------------------------------------------------------
#ifdef EMSCRIPTEN
#include <unistd.h>
#define _isatty isatty
#define _fileno fileno
#elifdef _WIN32
#include <windows.h>
#include <io.h>
#endif
#include <immintrin.h>
#include <emmintrin.h>
#include <cstddef>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <limits>
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <map>
#include <new>
#include <memory.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <locale>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <atomic>
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @namespace tlsf
/// @brief  Two-Level Segregated Fit (TLSF) memory allocator
///
/// @details
///   TLSF is a dynamic memory allocator that guarantees O(1) allocation
///   and deallocation time. It achieves this through a two-level segregated
///   free-list structure:
///
///   - **First level (fl)**: Divides block sizes into ranges using a
///     logarithmic (base-2) classification. Each first-level index
///     corresponds to a size range [2^i, 2^{i+1}).
///
///   - **Second level (sl)**: Subdivides each first-level range into
///     2^SL_INDEX_COUNT_LOG2 equally-spaced size classes, providing
///     finer granularity within each logarithmic bin.
///
///   Free blocks are stored in a 2D array \c blocks[fl][sl], and
///   occupancy is tracked via bitmaps (\c fl_bitmap and \c sl_bitmap[])
///   for O(1) lookup of the smallest suitable free block.
///
///   **Algorithm reference**:
///     - Masmano et al., "TLSF: A New Dynamic Memory Allocator for
///       Real-Time Systems", ECRTS 2004.
///     - Matthew Conte, TLSF implementation version 3.0:
///       http://tlsf.baisoku.org
///
///   **Complexity**:
///     - malloc: O(1) worst-case
///     - free:   O(1) worst-case (with coalescing)
///     - realloc: O(1) best-case (in-place growth), O(n) worst-case (copy)
///
///   **Key data structures**:
///     - block_header_t: in-band metadata for each block (size, free flags,
///       doubly-linked free-list pointers)
///     - control_t: allocator state (bitmaps, free-list array, sentinel)
///
///   This implementation is a singleton via function-local static
///   (Meyer's singleton pattern).
// ---------------------------------------------------------------------------
namespace tlsf {  // Two Level Segregated Fit memory allocator
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief TLSF algorithm implementation notes
///
/// Based on the original documentation by Miguel Masmano:
///   http://rtportal.upv.es/rtmalloc/allocators/tlsf/index.shtml
///
/// This implementation was written to the specification of the document,
/// therefore no GPL restrictions apply.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief Architecture-specific bit manipulation routines
///
/// TLSF achieves O(1) cost for malloc and free operations by limiting
/// the search for a free block to a free list of guaranteed size
/// adequate to fulfill the request, combined with efficient free list
/// queries using bitmasks and architecture-specific bit-manipulation
/// routines.
///
/// Most modern processors provide instructions to count leading zeroes
/// in a word, find the lowest and highest set bit, etc. These
/// specific implementations will be used when available, falling back
/// to a reasonably efficient generic implementation.
///
/// NOTE: TLSF spec relies on ffs/fls returning value 0..31.
/// ffs/fls return 1-32 by default, returning 0 for error.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief Detect 32- vs 64-bit architecture
///
/// There is no reliable portable method at compile-time.
// ---------------------------------------------------------------------------
#if defined (__alpha__) || defined (__ia64__) || defined (__x86_64__) \
	|| defined (_WIN64) || defined (__LP64__) || defined (__LLP64__)
#define TLSF_64BIT
#endif

// ---------------------------------------------------------------------------
/// @name tlsf_ffs / tlsf_fls — Find First Set and Find Last Set
///
/// @details
///   These functions find the position of the lowest (ffs) and highest (fls)
///   set bit in a word. They are the critical primitive for converting a
///   size to its (fl, sl) bin indices.
///
///   - tlsf_ffs(word): returns 0-based index of least significant set bit
///   - tlsf_fls(word): returns 0-based index of most significant set bit
///
///   Returns -1 if no bit is set (word == 0).
///
///   Platform implementations:
///     - GCC/Clang: __builtin_ffs / 32 - __builtin_clz
///     - MSVC x86/x64: _BitScanForward / _BitScanReverse
///     - MSVC PowerPC: _CountLeadingZeros
///     - ARMCC: __clz
///     - Green Hills: __CLZ32
///     - Generic fallback: binary search on byte/halfword/word boundaries
///
/// @see Masmano et al., Section 3.1 for mapping functions
// ---------------------------------------------------------------------------
//@{

// ---------------------------------------------------------------------------
/// @brief gcc 3.4+ have builtin support, specialized for architecture.
/// Some compilers masquerade as gcc; patchlevel test filters them out.
// ---------------------------------------------------------------------------
#if defined (__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)) \
	&& defined (__GNUC_PATCHLEVEL__)

/// @brief Find First Set (lowest set bit) using GCC builtin
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	return __builtin_ffs(word) - 1;
}

/// @brief Find Last Set (highest set bit) using GCC builtin
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __builtin_clz(word) : 0;
	return bit - 1;
}

#elif defined (_MSC_VER) && (_MSC_VER >= 1400) && (defined (_M_IX86) || defined (_M_X64))
// ---------------------------------------------------------------------------
/// @brief Microsoft Visual C++ support on x86/X64 architectures.
// ---------------------------------------------------------------------------

#include <intrin.h>

#pragma intrinsic(_BitScanReverse)
#pragma intrinsic(_BitScanForward)

/// @brief Find Last Set using MSVC _BitScanReverse
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	unsigned long index;
	return _BitScanReverse(&index, word) ? static_cast<int>(index) : -1;
}

/// @brief Find First Set using MSVC _BitScanForward
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	unsigned long index;
	return _BitScanForward(&index, word) ? static_cast<int>(index) : -1;
}

#elif defined (_MSC_VER) && defined (_M_PPC)
// ---------------------------------------------------------------------------
/// @brief Microsoft Visual C++ support on PowerPC architectures.
// ---------------------------------------------------------------------------

#include <ppcintrinsics.h>

/// @brief Find Last Set using MSVC PowerPC _CountLeadingZeros
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	const int bit = 32 - _CountLeadingZeros(word);
	return bit - 1;
}

/// @brief Find First Set via reverse-then-count leading zeros
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - _CountLeadingZeros(reverse);
	return bit - 1;
}

#elif defined (__ARMCC_VERSION)
// ---------------------------------------------------------------------------
/// @brief RealView Compilation Tools for ARM
// ---------------------------------------------------------------------------

/// @brief Find First Set using ARM __clz on reversed word
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - __clz(reverse);
	return bit - 1;
}

/// @brief Find Last Set using ARM __clz
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __clz(word) : 0;
	return bit - 1;
}

#elif defined (__ghs__)
// ---------------------------------------------------------------------------
/// @brief Green Hills support for PowerPC
// ---------------------------------------------------------------------------

#include <ppc_ghs.h>

/// @brief Find First Set using Green Hills __CLZ32 on reversed word
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - __CLZ32(reverse);
	return bit - 1;
}

/// @brief Find Last Set using Green Hills __CLZ32
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __CLZ32(word) : 0;
	return bit - 1;
}

#else
// ---------------------------------------------------------------------------
/// @brief Fall back to generic implementation using binary decomposition.
/// @details Uses a series of conditional shifts to narrow the position
///   of the highest set bit, checking 16-bit, 8-bit, 4-bit, 2-bit, and
///   1-bit windows sequentially.
// ---------------------------------------------------------------------------

/// @brief Generic Find Last Set implementation (returns 0..32)
/// @param word Input word
/// @return Bit position (1-indexed), or 0 if word == 0
inline int tlsf_fls_generic(unsigned int word)
{
	int bit = 32;

	if (!word) bit -= 1;
	if (!(word & 0xffff0000)) { word <<= 16; bit -= 16; }
	if (!(word & 0xff000000)) { word <<= 8; bit -= 8; }
	if (!(word & 0xf0000000)) { word <<= 4; bit -= 4; }
	if (!(word & 0xc0000000)) { word <<= 2; bit -= 2; }
	if (!(word & 0x80000000)) { word <<= 1; bit -= 1; }

	return bit;
}

/// @brief tlsf_ffs implemented in terms of tlsf_fls_generic
/// @param word Input word
/// @return 0-based index of lowest set bit, or -1 if word == 0
inline int tlsf_ffs(unsigned int word)
{
	return tlsf_fls_generic(word & (~word + 1)) - 1;
}

/// @brief tlsf_fls implemented in terms of tlsf_fls_generic
/// @param word Input word
/// @return 0-based index of highest set bit, or -1 if word == 0
inline int tlsf_fls(unsigned int word)
{
	return tlsf_fls_generic(word) - 1;
}

#endif

// ---------------------------------------------------------------------------
/// @brief 64-bit version of tlsf_fls for size_t arguments.
///
/// On 64-bit platforms, size_t is 64 bits wide; this function extends the
/// 32-bit fls to handle the upper 32 bits.
///
/// @param size Input size value
/// @return 0-based index of highest set bit, or -1 if size == 0
// ---------------------------------------------------------------------------
#if defined (TLSF_64BIT)
inline int tlsf_fls_sizet(size_t size)
{
	unsigned int high = static_cast<unsigned int>(size >> 32);
	int bits = 0;

	if( high ) {
		bits = 32 + tlsf_fls(high);
	}
	else {
		bits = tlsf_fls(static_cast<unsigned int>(size & 0xffffffff));
	}

	return bits;
}
#else
inline int tlsf_fls_sizet(size_t size){ return tlsf_fls(size); }
#endif
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @name TLSF Tunable Constants
///
/// @details These enums define the granularity and range of the two-level
///   free-list structure:
///
///   - ALIGN_SIZE_LOG2 / ALIGN_SIZE: Minimum alignment (4 bytes on 32-bit,
///     8 bytes on 64-bit). All blocks are rounded up to this alignment.
///   - SL_INDEX_COUNT_LOG2: log2 of the number of second-level bins per
///     first-level bin (default 5 → 32 second-level bins).
///   - SL_INDEX_COUNT = 1 << SL_INDEX_COUNT_LOG2
///   - FL_INDEX_SHIFT = SL_INDEX_COUNT_LOG2 + ALIGN_SIZE_LOG2: shift to
///     compute first-level index from size.
///   - FL_INDEX_MAX: Maximum first-level index (capped by max representable
///     block size).
///   - FL_INDEX_COUNT = FL_INDEX_MAX - FL_INDEX_SHIFT + 1
///   - SMALL_BLOCK_SIZE = 1 << FL_INDEX_SHIFT: Blocks below this size go
///     to fl=0 with linear second-level mapping.
// ---------------------------------------------------------------------------
//@{

/// @brief Public TLSF tunable constants
enum tlsf_public {
	SL_INDEX_COUNT_LOG2 = 5,  ///< log2 of second-level index count
};

/// @brief Private (internal) TLSF tunable constants
enum tlsf_private {
#if defined (TLSF_64BIT)
	ALIGN_SIZE_LOG2 = 3,  ///< log2 of alignment (8 bytes on 64-bit)
#else
	ALIGN_SIZE_LOG2 = 2,  ///< log2 of alignment (4 bytes on 32-bit)
#endif
	ALIGN_SIZE = (1 << ALIGN_SIZE_LOG2),  ///< Minimum alignment in bytes

#if defined (TLSF_64BIT)
	FL_INDEX_MAX = 32,    ///< Max first-level index (64-bit: 2^32 max block)
#else
	FL_INDEX_MAX = 30,    ///< Max first-level index (32-bit: 2^30 max block)
#endif
	SL_INDEX_COUNT = (1 << SL_INDEX_COUNT_LOG2),  ///< Number of second-level bins
	FL_INDEX_SHIFT = (SL_INDEX_COUNT_LOG2 + ALIGN_SIZE_LOG2),  ///< Shift for first-level mapping
	FL_INDEX_COUNT = (FL_INDEX_MAX - FL_INDEX_SHIFT + 1),       ///< Number of first-level bins

	SMALL_BLOCK_SIZE = (1 << FL_INDEX_SHIFT),  ///< Threshold for small-block path
};
//@}
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief Block header structure — in-band metadata for each memory block
///
/// @details Each memory block managed by TLSF has a header immediately
///   preceding the user-accessible data. The layout is:
///
///   | Offset | Field             | Description                                |
///   |--------|-------------------|--------------------------------------------|
///   | 0      | prev_phys_block   | Pointer to previous physical block         |
///   | 8/4    | size              | Block size (low 2 bits used as flags)      |
///   | 16/8   | next_free         | Next free block in same bin (free list)    |
///   | 24/12  | prev_free         | Previous free block in same bin (free list)|
///
///   Size field layout:
///   - bit 0: free bit (1 = free, 0 = used)
///   - bit 1: prev_free bit (1 = previous block is free)
///   - bits 2+: actual block size (always aligned)
///
///   The next_free and prev_free fields are only valid when the block is
///   free. When allocated, these bytes are part of the user payload.
///
///   Complexity: O(1) to read size or flags.
// ---------------------------------------------------------------------------
struct block_header_t {
	block_header_t * prev_phys_block;  ///< Pointer to previous block in physical memory order
	size_t size;                        ///< Block size with embedded flags (free, prev_free)
	block_header_t * next_free;         ///< Next free block in same bin (free-list linkage)
	block_header_t * prev_free;         ///< Previous free block in same bin (free-list linkage)
};

// ---------------------------------------------------------------------------
/// @brief Allocator control structure — allocator state
///
/// @details The control structure holds:
///   - block_null: sentinel block that terminates free lists (points to
///     itself when empty)
///   - fl_bitmap: one bit per first-level index (1 = list non-empty)
///   - sl_bitmap[FL_INDEX_COUNT]: one bit per second-level index per
///     first-level bin
///   - blocks[FL_INDEX_COUNT][SL_INDEX_COUNT]: 2D array of free-list heads
///
///   Memory layout: control_t is placed at the start of each TLSF pool.
// ---------------------------------------------------------------------------
struct control_t {
	block_header_t block_null;                    ///< Sentinel block (self-referential when empty)
	unsigned int fl_bitmap;                        ///< First-level bitmap (non-empty bins)
	unsigned int sl_bitmap[FL_INDEX_COUNT];        ///< Second-level bitmaps per first-level bin
	block_header_t * blocks[FL_INDEX_COUNT][SL_INDEX_COUNT];  ///< Free-list heads
};

// ---------------------------------------------------------------------------
/// @name Basic type aliases
// ---------------------------------------------------------------------------
//@{

using tlsfptr_t = ptrdiff_t;  ///< Pointer-width integer type for address arithmetic
using tlsf_t = void *;         ///< Opaque handle to a TLSF allocator instance
using pool_t = void *;         ///< Opaque handle to a memory pool
//@}

// ---------------------------------------------------------------------------
/// @class TLSF_Impl
/// @brief Singleton TLSF memory allocator implementation
///
/// @details Provides malloc, realloc, free with O(1) cost. Internally
///   maintains a set of memory pools and a single control structure.
///   Thread-safety is not provided by this class; callers must provide
///   external synchronization.
///
///   The singleton pattern uses a function-local static instance, ensuring
///   thread-safe initialization (C++11 guarantee).
///
///   @see http://tlsf.baisoku.org
///   @see Masmano et al., ECRTS 2004
// ---------------------------------------------------------------------------
class TLSF_Impl { // must be singleton
    private:
        TLSF_Impl(const TLSF_Impl &) = delete;
        void operator = (const TLSF_Impl &) = delete;
	public:
		~TLSF_Impl() {
			assert( ref_count_ == 0 );
			assert( pools_ == nullptr );
			assert( tlsf_ == nullptr );
		}

		TLSF_Impl() : ref_count_(0) {
			assert( ref_count_ == 0 );
			assert( pools_ == nullptr );
			assert( tlsf_ == nullptr );
		}

		/// @brief Allocate a block of memory
		/// @param sz Requested size in bytes
		/// @return Pointer to allocated memory, or nullptr on failure
		void * malloc(size_t sz);

		/// @brief Reallocate a memory block
		/// @param p Previously allocated pointer (may be nullptr)
		/// @param sz New size in bytes
		/// @return Pointer to reallocated memory, or nullptr on failure
		void * realloc(void * p, size_t sz);

		/// @brief Free a memory block
		/// @param p Pointer to previously allocated memory
		void free(void * p);

	protected:
		uintptr_t ref_count_;   ///< Reference count (number of outstanding allocations)
		pool_t * pools_ = nullptr;  ///< Array of pool pointers (slot 0 = count)
		tlsf_t tlsf_ = nullptr;     ///< Opaque pointer to TLSF control structure

		// -----------------------------------------------------------------------
		/// @brief Inline quicksort implementation
		///
		/// @details Uses an explicit stack (instead of recursion) to avoid
		///   stack overflow. The stack capacity is sizeof(void*) * CHAR_BIT
		///   entries, sufficient for any array representable by a pointer.
		///   Always sorts the smaller partition first to minimize stack depth.
		///
		/// @tparam T Element type
		/// @tparam F Comparator functor type (returning -1, 0, or 1)
		/// @param p Pointer to the array (offset by lb)
		/// @param lb Lower bound index (inclusive)
		/// @param ub Upper bound index (inclusive)
		/// @param f Comparator functor
		///
		/// @complexity O(n log n) average, O(n^2) worst-case
		/// @space O(log n) stack space
		// -----------------------------------------------------------------------
		template <typename T,typename F>
		static void qsort(T * p, intptr_t lb, intptr_t ub, F f) {
			constexpr size_t MAX_LEVELS = sizeof(void *) * CHAR_BIT;
			intptr_t b[MAX_LEVELS], e[MAX_LEVELS], i = 0, L, R, sw;

			p += lb;

			b[0] = 0;
			e[0] = ub - lb + 1;

			while( i >= 0 ){
    			L = b[i];
				R = e[i] - 1;

				if( L < R ){
					T piv(p[L]);

					while ( L < R ){
						while( f(p[R],piv) >= 0 && L < R ) R--;
						if( L < R ) p[L++] = p[R];
						while( f(p[L],piv) <= 0 && L < R) L++;
						if( L < R ) p[R--] = p[L];
					}

					p[L] = piv;

					b[i + 1] = L + 1;
					e[i + 1] = e[i];
					e[i++] = L;

					if( e[i] - b[i] > e[i - 1] - b[i - 1] ){
						sw = b[i];
						b[i] = b[i - 1];
						b[i-1] = sw;

						sw = e[i];
						e[i] = e[i - 1];
						e[i - 1] = sw;
					}
				}
				else {
					i--;
				}
			}
		}

		// -----------------------------------------------------------------------
		/// @brief Inline quicksort with default comparator (ascending order)
		/// @tparam T Element type (must support operator< and operator>)
		/// @param p Pointer to the array
		/// @param lb Lower bound index
		/// @param ub Upper bound index
		// -----------------------------------------------------------------------
		template <typename T>
		static void qsort(T * p, intptr_t lb, intptr_t ub)	{
			qsort<T>(p,lb,ub,[] (const T & a,const T & b) {
				return a > b ? 1 : a < b ? -1 : 0;
			});
		}

		// -----------------------------------------------------------------------
		/// @brief Inline binary search with custom comparator
		///
		/// @tparam T Key type
		/// @tparam F Comparator functor type
		/// @param keys Sorted array of keys
		/// @param low Lower bound index (inclusive)
		/// @param high Upper bound index (inclusive)
		/// @param key Key to search for
		/// @param f Comparator functor
		/// @return Index of matching element, or -1 if not found
		///
		/// @complexity O(log n)
		// -----------------------------------------------------------------------
		template <typename T,typename F>
		static intptr_t bsearch(T * keys,intptr_t low, intptr_t high, const T & key, F f) {
			for(;;){
				uintptr_t p = (low + high) / 2;

				if( low > high ) break;

				intptr_t c = f(key,keys[p]);

				if( c > 0 ){
					low = static_cast<intptr_t>(p) + 1;
				}
				else if( c < 0 ){
					high = static_cast<intptr_t>(p) - 1;
				}
				else
					return static_cast<intptr_t>(p);
			}

			return -1;
		}

		// -----------------------------------------------------------------------
		/// @brief Inline binary search with default comparator
		/// @tparam T Key type
		/// @param keys Sorted array of keys
		/// @param low Lower bound index
		/// @param high Upper bound index
		/// @param key Key to search for
		/// @return Index of matching element, or -1 if not found
		// -----------------------------------------------------------------------
		template <typename T>
		static intptr_t bsearch(T * keys,intptr_t low, intptr_t high, const T & key) {
			return bsearch<T>(keys,low,high,key,[] (const T & key,const T & b) {
				return key > b ? 1 : key < b ? -1 : 0;
			});
		}

	private:
		// -----------------------------------------------------------------------
		/// @brief Constexpr minimum of two values
		// -----------------------------------------------------------------------
		template <typename T> static constexpr const T & tlsf_min(const T & a,const T & b) { return a < b ? a : b; }

		// -----------------------------------------------------------------------
		/// @brief Constexpr maximum of two values
		// -----------------------------------------------------------------------
		template <typename T> static constexpr const T & tlsf_max(const T & a,const T & b) { return a > b ? a : b; }

		// -----------------------------------------------------------------------
		/// @name Block size constants (all constexpr)
		// -----------------------------------------------------------------------
		//@{

		static constexpr size_t pool_size()						{ return size_t(64) * 1024u * 1024u; }
		static constexpr size_t block_header_free_bit()			{ return size_t(1) << 0; }
		static constexpr size_t block_header_prev_free_bit()	{ return size_t(1) << 1; }
		static constexpr size_t block_header_overhead()			{ return sizeof(size_t); }
		static constexpr size_t block_start_offset()			{ return offsetof(block_header_t, size) + sizeof(size_t); }
		static constexpr size_t block_size_min()				{ return sizeof(block_header_t) - sizeof(block_header_t *); }
		static constexpr size_t block_size_max()				{ return size_t(1) << FL_INDEX_MAX; }
		//@}

		// -----------------------------------------------------------------------
		/// @brief Extract the usable size from a block header (mask out flag bits)
		///
		/// The low 2 bits of the size field are used as free/prev_free flags.
		/// This function returns size with those bits cleared.
		///
		/// @param block Pointer to block header
		/// @return Block size (always aligned)
		// -----------------------------------------------------------------------
		static size_t block_size(const block_header_t * block) {
			return block->size & ~(block_header_free_bit() | block_header_prev_free_bit());
		}

		// -----------------------------------------------------------------------
		/// @brief Set the size field of a block header, preserving flag bits
		///
		/// @param block Pointer to block header
		/// @param size New size value (flag bits will be preserved)
		// -----------------------------------------------------------------------
		static void block_set_size(block_header_t* block, size_t size) {
			const size_t oldsize = block->size;
			block->size = size | (oldsize & (block_header_free_bit() | block_header_prev_free_bit()));
		}

		// -----------------------------------------------------------------------
		/// @brief Check if a block is the last block in the pool (size == 0)
		/// @param block Pointer to block header
		/// @return 1 if last block, 0 otherwise
		// -----------------------------------------------------------------------
		static int block_is_last(const block_header_t * block) {
			return 0 == block_size(block);
		}

		// -----------------------------------------------------------------------
		/// @brief Check if a block is free
		/// @param block Pointer to block header
		/// @return 1 if free, 0 if used
		// -----------------------------------------------------------------------
		static int block_is_free(const block_header_t * block) {
			return static_cast<int>(block->size & block_header_free_bit());
		}

		// -----------------------------------------------------------------------
		/// @brief Mark a block as free (set the free bit)
		/// @param block Pointer to block header
		// -----------------------------------------------------------------------
		static void block_set_free(block_header_t * block) {
			block->size |= block_header_free_bit();
		}

		// -----------------------------------------------------------------------
		/// @brief Mark a block as used (clear the free bit)
		/// @param block Pointer to block header
		// -----------------------------------------------------------------------
		static void block_set_used(block_header_t* block) {
			block->size &= ~block_header_free_bit();
		}

		// -----------------------------------------------------------------------
		/// @brief Check if the previous block is free
		/// @param block Pointer to block header
		/// @return 1 if previous block is free, 0 otherwise
		// -----------------------------------------------------------------------
		static int block_is_prev_free(const block_header_t * block)	{
			return static_cast<int>(block->size & block_header_prev_free_bit());
		}

		// -----------------------------------------------------------------------
		/// @brief Mark the previous block as free
		/// @param block Pointer to block header
		// -----------------------------------------------------------------------
		static void block_set_prev_free(block_header_t* block) {
			block->size |= block_header_prev_free_bit();
		}

		// -----------------------------------------------------------------------
		/// @brief Mark the previous block as used
		/// @param block Pointer to block header
		// -----------------------------------------------------------------------
		static void block_set_prev_used(block_header_t* block) {
			block->size &= ~block_header_prev_free_bit();
		}

		// -----------------------------------------------------------------------
		/// @brief Given a user pointer, recover the block header
		///
		/// The block header immediately precedes the user data at a fixed
		/// negative offset (block_start_offset).
		///
		/// @param ptr Pointer to user data
		/// @return Pointer to the block_header_t for this allocation
		///
		/// @pre ptr was returned by a previous TLSF allocation
		// -----------------------------------------------------------------------
		static block_header_t * block_from_ptr(const void * ptr) {
			return reinterpret_cast<block_header_t *>(
				reinterpret_cast<uintptr_t>(ptr) - block_start_offset());
		}

		// -----------------------------------------------------------------------
		/// @brief Given a block header, return the user data pointer
		/// @param block Pointer to block header
		/// @return Pointer to the user data area
		// -----------------------------------------------------------------------
		static void * block_to_ptr(const block_header_t * block) {
			return reinterpret_cast<void *>(
				reinterpret_cast<uintptr_t>(block) + block_start_offset());
		}

		// -----------------------------------------------------------------------
		/// @brief Given a base pointer and offset, compute pointer to a block header
		/// @param ptr Base pointer
		/// @param size Offset in bytes
		/// @return Pointer to block_header_t at the computed address
		// -----------------------------------------------------------------------
		static block_header_t * offset_to_block(const void * ptr, size_t size) {
			return reinterpret_cast<block_header_t *>(
				reinterpret_cast<tlsfptr_t>(ptr) + static_cast<tlsfptr_t>(size));
		}

		// -----------------------------------------------------------------------
		/// @brief Get the previous physical block
		/// @param block Pointer to current block
		/// @return Pointer to previous block header
		// -----------------------------------------------------------------------
		static block_header_t * block_prev(const block_header_t * block) {
			return block->prev_phys_block;
		}

		// -----------------------------------------------------------------------
		/// @brief Get the next physical block
		///
		/// The next block is located at: block_data + block_size - overhead.
		/// The overhead (sizeof(size_t)) accounts for the size field of the
		/// current block that is "shared" with the next block's prev_phys_block.
		///
		/// @param block Pointer to current block
		/// @return Pointer to next block header
		// -----------------------------------------------------------------------
		static block_header_t * block_next(const block_header_t * block) {
			block_header_t * next = offset_to_block(block_to_ptr(block),
				block_size(block) - block_header_overhead());
			return next;
		}

		// -----------------------------------------------------------------------
		/// @brief Link the next block's prev_phys_block to the current block
		/// @param block Pointer to current block
		/// @return Pointer to the next block
		// -----------------------------------------------------------------------
		static block_header_t * block_link_next(block_header_t * block) {
			block_header_t * next = block_next(block);
			next->prev_phys_block = block;
			return next;
		}

		// -----------------------------------------------------------------------
		/// @brief Mark a block as free and update the next block's prev_free flag
		/// @param block Pointer to block to mark as free
		// -----------------------------------------------------------------------
		static void block_mark_as_free(block_header_t* block) {
			block_header_t * next = block_link_next(block);
			block_set_prev_free(next);
			block_set_free(block);
		}

		// -----------------------------------------------------------------------
		/// @brief Mark a block as used and update the next block's prev_free flag
		/// @param block Pointer to block to mark as used
		// -----------------------------------------------------------------------
		static void block_mark_as_used(block_header_t * block) {
			block_header_t * next = block_next(block);
			block_set_prev_used(next);
			block_set_used(block);
		}

		// -----------------------------------------------------------------------
		/// @brief Align a size upward to the given alignment
		/// @param x Input size
		/// @param align Alignment (must be power of 2)
		/// @return Aligned size
		// -----------------------------------------------------------------------
		static size_t align_up(size_t x, size_t align) {
			return (x + (align - 1)) & ~(align - 1);
		}

		// -----------------------------------------------------------------------
		/// @brief Align a size downward to the given alignment
		/// @param x Input size
		/// @param align Alignment (must be power of 2)
		/// @return Aligned size
		// -----------------------------------------------------------------------
		static size_t align_down(size_t x, size_t align) {
			return x - (x & (align - 1));
		}

		// -----------------------------------------------------------------------
		/// @brief Align a pointer upward to the given alignment
		/// @param ptr Input pointer
		/// @param align Alignment (must be power of 2)
		/// @return Aligned pointer
		// -----------------------------------------------------------------------
		static void * align_ptr(const void* ptr, size_t align) {
			return reinterpret_cast<void *>(
				(reinterpret_cast<tlsfptr_t>(ptr) + static_cast<tlsfptr_t>(align - 1))
				& ~static_cast<tlsfptr_t>(align - 1));
		}

		// -----------------------------------------------------------------------
		/// @brief Adjust a request size to the minimum block size and alignment
		///
		/// If the requested size is valid (non-zero and below max), rounds up
		/// to alignment and ensures at least block_size_min().
		///
		/// @param size Requested size
		/// @param align Alignment
		/// @return Adjusted size, or 0 if request is invalid
		// -----------------------------------------------------------------------
		static size_t adjust_request_size(size_t size, size_t align) {
			size_t adjust = 0;

			if( size && size < block_size_max() )
				adjust = tlsf_max(align_up(size, align), block_size_min());

			return adjust;
		}

		// -----------------------------------------------------------------------
		/// @brief Mapping function: size → (first-level, second-level) bin indices
		///
		/// @details This is the core mapping that translates a block size into
		///   its (fl, sl) bin coordinates.
		///
		///   For small blocks (size < SMALL_BLOCK_SIZE), fl = 0 and sl is a
		///   linear division of the small-block range.
		///
		///   For larger blocks, fl is derived from tlsf_fls (the position of
		///   the highest set bit in the size), and sl is derived from the next
		///   SL_INDEX_COUNT_LOG2 bits below the MSB.
		///
		///   Formally:
		///     fl = floor(log2(size))        for size >= SMALL_BLOCK_SIZE
		///     sl = (size >> (fl - SL_INDEX_COUNT_LOG2)) ^ (1 << SL_INDEX_COUNT_LOG2)
		///
		/// @param[in]  size Block size (adjusted, aligned)
		/// @param[out] fli  First-level index
		/// @param[out] sli  Second-level index
		///
		/// @complexity O(1)
		/// @see Masmano et al., Section 3.2: "Mapping Function"
		// -----------------------------------------------------------------------
		static void mapping_insert(size_t size, int * fli, int * sli) {
			int fl, sl;

			if( size < SMALL_BLOCK_SIZE ) {
				fl = 0;
				sl = static_cast<int>(size) / (SMALL_BLOCK_SIZE / SL_INDEX_COUNT);
			}
			else {
				fl = tlsf_fls_sizet(size);
				sl = static_cast<int>(size >> (fl - SL_INDEX_COUNT_LOG2)) ^ (1 << SL_INDEX_COUNT_LOG2);
				fl -= (FL_INDEX_SHIFT - 1);
			}

			*fli = fl;
			*sli = sl;
		}

		// -----------------------------------------------------------------------
		/// @brief Search mapping: round up the size before mapping
		///
		/// @details When searching for a free block, we need to ensure we
		///   round up the requested size so that the mapping gives us a bin
		///   that can definitely satisfy the request. Adds a rounding term
		///   before calling mapping_insert.
		///
		/// @param[in]  size Requested size
		/// @param[out] fli  First-level index
		/// @param[out] sli  Second-level index
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void mapping_search(size_t size, int * fli, int * sli) {
			if( size >= (static_cast<size_t>(1) << SL_INDEX_COUNT_LOG2) ) {
				const size_t round = (static_cast<size_t>(1) << (tlsf_fls_sizet(size) - SL_INDEX_COUNT_LOG2)) - 1;
				size += round;
			}
			mapping_insert(size, fli, sli);
		}

		// -----------------------------------------------------------------------
		/// @brief Search for a suitable free block given (fl, sl) hints
		///
		/// @details Starting from the given (fl, sl) bin, checks if a free
		///   block exists. If not, scans upward through second-level bitmaps
		///   within the same first-level bin, then to higher first-level bins,
		///   using bitmap operations for O(1) search.
		///
		///   - Uses sl_bitmap[fl] with mask to find the smallest sl >= requested
		///   - If no block in current fl, uses fl_bitmap to find next non-empty
		///     first-level bin
		///
		/// @param[in]  control Pointer to control structure
		/// @param[in,out] fli  First-level index (may be updated)
		/// @param[in,out] sli  Second-level index (may be updated)
		/// @return Pointer to a suitable free block, or nullptr if none found
		///
		/// @complexity O(1) worst-case (bounded by bitmap scan)
		// -----------------------------------------------------------------------
		static block_header_t * search_suitable_block(control_t * control, int * fli, int * sli) {
			int fl = *fli;
			int sl = *sli;

			unsigned int sl_map = control->sl_bitmap[fl] & (~0u << sl);

			if( !sl_map ) {
				/* No block exists. Search in the next largest first-level list. */
				const unsigned int fl_map = control->fl_bitmap & (~0u << (fl + 1));

				if( !fl_map )
					return nullptr;

				fl = tlsf_ffs(fl_map);
				*fli = fl;
				sl_map = control->sl_bitmap[fl];
			}

			sl = tlsf_ffs(sl_map);
			*sli = sl;

			return control->blocks[fl][sl];
		}

		// -----------------------------------------------------------------------
		/// @brief Remove a free block from its free list
		///
		/// @details Unlinks the block from the doubly-linked free list at
		///   (fl, sl). If the block was the head of the list, updates the
		///   bitmap (clears sl_bitmap entry, and if the bin becomes empty,
		///   clears fl_bitmap entry).
		///
		/// @param control Pointer to control structure
		/// @param block Block to remove
		/// @param fl First-level index
		/// @param sl Second-level index
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void remove_free_block(control_t * control, block_header_t * block, int fl, int sl) {
			block_header_t* prev = block->prev_free;
			block_header_t* next = block->next_free;
			next->prev_free = prev;
			prev->next_free = next;

			if( control->blocks[fl][sl] == block ) {
				control->blocks[fl][sl] = next;

				if( next == &control->block_null ) {
					control->sl_bitmap[fl] &= ~(static_cast<unsigned int>(1) << sl);
					if( !control->sl_bitmap[fl])
						control->fl_bitmap &= ~(static_cast<unsigned int>(1) << fl);
				}
			}
		}

		// -----------------------------------------------------------------------
		/// @brief Insert a free block into the appropriate free list
		///
		/// @details Inserts the block at the head of the free list for bin
		///   (fl, sl) and updates the bitmaps.
		///
		/// @param control Pointer to control structure
		/// @param block Block to insert
		/// @param fl First-level index
		/// @param sl Second-level index
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void insert_free_block(control_t * control, block_header_t * block, int fl, int sl) {
			block_header_t * current = control->blocks[fl][sl];
			block->next_free = current;
			block->prev_free = &control->block_null;
			current->prev_free = block;
			control->blocks[fl][sl] = block;
			control->fl_bitmap |= (static_cast<unsigned int>(1) << fl);
			control->sl_bitmap[fl] |= (static_cast<unsigned int>(1) << sl);
		}

		// -----------------------------------------------------------------------
		/// @brief Remove a block from the free lists (mapping + remove)
		/// @param control Pointer to control structure
		/// @param block Block to remove
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void block_remove(control_t * control, block_header_t * block) {
			int fl, sl;
			mapping_insert(block_size(block), &fl, &sl);
			remove_free_block(control, block, fl, sl);
		}

		// -----------------------------------------------------------------------
		/// @brief Insert a block into the free lists (mapping + insert)
		/// @param control Pointer to control structure
		/// @param block Block to insert
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void block_insert(control_t * control, block_header_t * block) {
			int fl, sl;
			mapping_insert(block_size(block), &fl, &sl);
			insert_free_block(control, block, fl, sl);
		}

		// -----------------------------------------------------------------------
		/// @brief Check if a block can be split to satisfy a request
		///
		/// A block can be split if the remaining portion (after removing
		/// the requested size) is at least as large as a block header.
		///
		/// @param block Block to check
		/// @param size Requested size (adjusted, aligned)
		/// @return 1 if splittable, 0 otherwise
		// -----------------------------------------------------------------------
		static int block_can_split(block_header_t * block, size_t size) {
			return block_size(block) >= sizeof(block_header_t) + size;
		}

		// -----------------------------------------------------------------------
		/// @brief Split a block into allocated and free portions
		///
		/// @details Given a block larger than the requested size, splits off
		///   the remainder as a new free block. The remainder is placed
		///   immediately after the allocated portion in physical memory.
		///
		/// @param block Block to split (becomes the allocated portion)
		/// @param size Requested size
		/// @return Pointer to the new free remainder block
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_split(block_header_t* block, size_t size) {
			block_header_t * remaining =
				offset_to_block(block_to_ptr(block), size - block_header_overhead());

			const size_t remain_size = block_size(block) - (size + block_header_overhead());

			block_set_size(remaining, remain_size);

			block_set_size(block, size);
			block_mark_as_free(remaining);

			return remaining;
		}

		// -----------------------------------------------------------------------
		/// @brief Absorb the next block into the current block (coalescing)
		///
		/// @details Adds the next block's size (plus header overhead) to the
		///   current block. Used for both forward and backward coalescing.
		///
		/// @param prev Block that will absorb (the survivor)
		/// @param block Block to absorb
		/// @return Pointer to prev (the coalesced block)
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_absorb(block_header_t * prev, block_header_t * block) {
			prev->size += block_size(block) + block_header_overhead();
			block_link_next(prev);
			return prev;
		}

		// -----------------------------------------------------------------------
		/// @brief Coalesce a block with the previous block if it is free
		///
		/// @details Checks the prev_free flag. If set, removes the previous
		///   block from the free lists, absorbs it into the current block.
		///
		/// @param control Pointer to control structure
		/// @param block Current block (may be merged with predecessor)
		/// @return Pointer to the (possibly merged) block
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_merge_prev(control_t * control, block_header_t * block) {
			if( block_is_prev_free(block) )	{
				block_header_t* prev = block_prev(block);
				block_remove(control, prev);
				block = block_absorb(prev, block);
			}
			return block;
		}

		// -----------------------------------------------------------------------
		/// @brief Coalesce a block with the next block if it is free
		///
		/// @details Checks the free bit of the next block. If set, removes
		///   the next block from the free lists, absorbs it into the current
		///   block.
		///
		/// @param control Pointer to control structure
		/// @param block Current block (may be merged with successor)
		/// @return Pointer to the (possibly merged) block
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_merge_next(control_t * control, block_header_t * block) {
			block_header_t * next = block_next(block);

			if( block_is_free(next) ) {
				block_remove(control, next);
				block = block_absorb(block, next);
			}

			return block;
		}

		// -----------------------------------------------------------------------
		/// @brief Trim a free block by splitting off the used portion
		///
		/// @details If the block is large enough to split, splits off the
		///   remainder as a new free block and inserts it into the free lists.
		///
		/// @param control Pointer to control structure
		/// @param block Block to trim (will be the allocated portion)
		/// @param size Requested size
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void block_trim_free(control_t * control, block_header_t * block, size_t size) {
			if( block_can_split(block, size) ) {
				block_header_t* remaining_block = block_split(block, size);
				block_link_next(block);
				block_set_prev_free(remaining_block);
				block_insert(control, remaining_block);
			}
		}

		// -----------------------------------------------------------------------
		/// @brief Trim a used block by splitting off excess and freeing it
		///
		/// @details If the block is larger than needed, splits the excess as
		///   a new free block, then merges it with any following free block
		///   and inserts into free lists.
		///
		/// @param control Pointer to control structure
		/// @param block Block to trim
		/// @param size New size for the block
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void block_trim_used(control_t * control, block_header_t * block, size_t size) {
			if( block_can_split(block, size) ) {
				block_header_t* remaining_block = block_split(block, size);
				block_set_prev_used(remaining_block);
				remaining_block = block_merge_next(control, remaining_block);
				block_insert(control, remaining_block);
			}
		}

		// -----------------------------------------------------------------------
		/// @brief Trim leading space from a block to create an aligned allocation
		///
		/// @details Used by memalign to discard the gap between the block
		///   start and the aligned user pointer. The discarded leading portion
		///   is inserted into free lists.
		///
		/// @param control Pointer to control structure
		/// @param block Block to trim
		/// @param size Size of the gap/trim
		/// @return Pointer to the remainder (the aligned portion)
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_trim_free_leading(control_t * control, block_header_t * block, size_t size) {
			block_header_t * remaining_block = block;
			if( block_can_split(block, size) ) {
				remaining_block = block_split(block, size - block_header_overhead());
				block_set_prev_free(remaining_block);
				block_link_next(block);
				block_insert(control, block);
			}

			return remaining_block;
		}

		// -----------------------------------------------------------------------
		/// @brief Locate a free block of at least the given size
		///
		/// @details This is the primary search entry point for malloc. Maps
		///   the requested size to (fl, sl) via mapping_search, then finds a
		///   suitable block via search_suitable_block. Removes the found block
		///   from free lists.
		///
		/// @param control Pointer to control structure
		/// @param size Adjusted request size
		/// @return Pointer to a free block, or nullptr if none found
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static block_header_t * block_locate_free(control_t * control, size_t size) {
			int fl = 0, sl = 0;
			block_header_t * block = nullptr;

			if( size ) {
				mapping_search(size, &fl, &sl);
				block = search_suitable_block(control, &fl, &sl);
			}

			if( block != nullptr )
				remove_free_block(control, block, fl, sl);

			return block;
		}

		// -----------------------------------------------------------------------
		/// @brief Prepare a block for use: trim, mark used, compute user pointer
		///
		/// @param control Pointer to control structure
		/// @param block Block to prepare
		/// @param size Requested size
		/// @return User data pointer, or nullptr if block is null
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void * block_prepare_used(control_t * control, block_header_t * block, size_t size) {
			void * p = nullptr;

			if( block != nullptr ) {
				block_trim_free(control, block, size);
				block_mark_as_used(block);
				p = block_to_ptr(block);
			}

			return p;
		}

		// -----------------------------------------------------------------------
		/// @brief Initialize (construct) the control structure
		///
		/// @details Sets up the sentinel block (block_null) to be
		///   self-referential, clears bitmaps, and initializes all free-list
		///   heads to point to block_null.
		///
		/// @param control Pointer to uninitialized control structure memory
		///
		/// @complexity O(FL_INDEX_COUNT * SL_INDEX_COUNT)
		// -----------------------------------------------------------------------
		static void control_construct(control_t * control) {
			control->block_null.next_free = &control->block_null;
			control->block_null.prev_free = &control->block_null;

			control->fl_bitmap = 0;

			for( intptr_t i = 0; i < static_cast<intptr_t>(FL_INDEX_COUNT); ++i ) {
				control->sl_bitmap[i] = 0;

				for( intptr_t j = 0; j < static_cast<intptr_t>(SL_INDEX_COUNT); ++j )
					control->blocks[i][j] = &control->block_null;
			}
		}

		// -----------------------------------------------------------------------
		/// @name TLSF size/layout queries
		// -----------------------------------------------------------------------
		//@{

		static constexpr size_t tlsf_size() {
			return sizeof(control_t);
		}

		static constexpr size_t tlsf_align_size() {
			return ALIGN_SIZE;
		}

		static constexpr size_t tlsf_block_size_min() {
			return block_size_min();
		}

		static constexpr size_t tlsf_block_size_max() {
			return block_size_max();
		}

		static constexpr size_t tlsf_pool_overhead() {
			return 2 * block_header_overhead();
		}

		static constexpr size_t tlsf_alloc_overhead() {
			return block_header_overhead();
		}

		static constexpr size_t tlsf_pool_max_block() {
			return pool_size() - tlsf_pool_overhead() - block_header_overhead() - tlsf_size();
		}
		//@}

		// -----------------------------------------------------------------------
		/// @brief Add a memory pool to the TLSF allocator
		///
		/// @details Carves a pool from the given memory region, creating an
		///   initial free block and a terminating sentinel block. The pool
		///   memory must be aligned to ALIGN_SIZE.
		///
		/// @param tlsf Opaque TLSF handle
		/// @param mem Pointer to pool memory (must be ALIGN_SIZE-aligned)
		/// @param bytes Size of the pool in bytes
		/// @return Pool handle (same as mem) on success, nullptr on failure
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static pool_t tlsf_add_pool(tlsf_t tlsf, void * mem, size_t bytes) {
			block_header_t * block;
			block_header_t * next;

			const size_t pool_overhead = tlsf_pool_overhead();
			const size_t pool_bytes = align_down(bytes - pool_overhead, ALIGN_SIZE);

			if( reinterpret_cast<ptrdiff_t>(mem) % static_cast<ptrdiff_t>(ALIGN_SIZE) != 0)
				return nullptr;

			if( pool_bytes < block_size_min() || pool_bytes > block_size_max() )
				return nullptr;

			block = offset_to_block(mem, -static_cast<tlsfptr_t>(block_header_overhead()));
			block_set_size(block, pool_bytes);
			block_set_free(block);
			block_set_prev_used(block);
			block_insert(reinterpret_cast<control_t *>(tlsf), block);

			next = block_link_next(block);
			block_set_size(next, 0);
			block_set_used(next);
			block_set_prev_free(next);

			return mem;
		}

		// -----------------------------------------------------------------------
		/// @brief Remove a memory pool from the TLSF allocator
		///
		/// @param tlsf Opaque TLSF handle
		/// @param pool Pool handle to remove
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void tlsf_remove_pool(tlsf_t tlsf, pool_t pool) {
			control_t * control = reinterpret_cast<control_t *>(tlsf);
			block_header_t * block = offset_to_block(pool, -static_cast<int>(block_header_overhead()));

			int fl = 0, sl = 0;

			mapping_insert(block_size(block), &fl, &sl);
			remove_free_block(control, block, fl, sl);
		}

		// -----------------------------------------------------------------------
		/// @brief Create a TLSF control structure in the given memory region
		///
		/// @param mem Pointer to memory for control structure (must be aligned)
		/// @return TLSF handle, or nullptr if alignment is wrong
		///
		/// @complexity O(FL_INDEX_COUNT * SL_INDEX_COUNT) (due to control_construct)
		// -----------------------------------------------------------------------
		static tlsf_t tlsf_create(void * mem) {
			if( reinterpret_cast<tlsfptr_t>(mem) % static_cast<tlsfptr_t>(ALIGN_SIZE) != 0 )
				return nullptr;

			control_construct(reinterpret_cast<control_t *>(mem));

			return static_cast<tlsf_t>(mem);
		}

		// -----------------------------------------------------------------------
		/// @brief Create a TLSF allocator with an initial pool
		///
		/// @param mem Memory for both control structure and pool
		/// @param bytes Total size of memory region
		/// @return TLSF handle, or nullptr on failure
		///
		/// @complexity O(FL_INDEX_COUNT * SL_INDEX_COUNT)
		// -----------------------------------------------------------------------
		static tlsf_t tlsf_create_with_pool(void * mem, size_t bytes) {
			tlsf_t tlsf = tlsf_create(mem);
			tlsf_add_pool(tlsf, reinterpret_cast<char*>(mem) + tlsf_size(), bytes - tlsf_size());
			return tlsf;
		}

		// -----------------------------------------------------------------------
		/// @brief Destroy a TLSF allocator (no-op in this implementation)
		/// @param tlsf TLSF handle (unused)
		// -----------------------------------------------------------------------
		static void tlsf_destroy(tlsf_t tlsf) {
			(void)tlsf;
		}

		// -----------------------------------------------------------------------
		/// @brief Get the pool handle from a TLSF allocator
		/// @param tlsf TLSF handle
		/// @return Pool handle
		// -----------------------------------------------------------------------
		pool_t tlsf_get_pool(tlsf_t tlsf) {
			return reinterpret_cast<pool_t>(reinterpret_cast<char *>(tlsf) + tlsf_size());
		}

		// -----------------------------------------------------------------------
		/// @brief TLSF malloc implementation
		///
		/// @details Adjusts the request size, locates a free block, and
		///   prepares it for use (trim, mark used).
		///
		/// @param tlsf TLSF handle
		/// @param size Requested size in bytes
		/// @return User data pointer, or nullptr on allocation failure
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void * tlsf_malloc(tlsf_t tlsf, size_t size) {
			control_t * control = reinterpret_cast<control_t *>(tlsf);
			const size_t adjust = adjust_request_size(size, ALIGN_SIZE);
			block_header_t * block = block_locate_free(control, adjust);
			return block_prepare_used(control, block, adjust);
		}

		// -----------------------------------------------------------------------
		/// @brief TLSF aligned memory allocation
		///
		/// @details Allocates memory with a specified alignment. Uses a
		///   worst-case allocation that accounts for alignment overhead,
		///   then trims the leading gap to achieve the desired alignment.
		///
		/// @param tlsf TLSF handle
		/// @param align Required alignment (must be power of 2)
		/// @param size Requested size in bytes
		/// @return Aligned user data pointer, or nullptr on failure
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void * tlsf_memalign(tlsf_t tlsf, size_t align, size_t size) {
			control_t * control = reinterpret_cast<control_t *>(tlsf);
			const size_t adjust = adjust_request_size(size, ALIGN_SIZE);
			const size_t gap_minimum = sizeof(block_header_t);
			const size_t size_with_gap = adjust_request_size(adjust + align + gap_minimum, align);
			const size_t aligned_size = (align <= ALIGN_SIZE) ? adjust : size_with_gap;
			block_header_t * block = block_locate_free(control, aligned_size);

			if( block != nullptr ) {
				void * ptr = block_to_ptr(block);
				void * aligned_ptr = align_ptr(ptr, align);
				size_t gap = static_cast<size_t>(reinterpret_cast<tlsfptr_t>(aligned_ptr) - reinterpret_cast<tlsfptr_t>(ptr));

				if( gap && gap < gap_minimum ) {
					const size_t gap_remain = gap_minimum - gap;
					const size_t offset = tlsf_max(gap_remain, align);
					const void * next_aligned = reinterpret_cast<void *>(reinterpret_cast<tlsfptr_t>(aligned_ptr) + static_cast<tlsfptr_t>(offset));

					aligned_ptr = align_ptr(next_aligned, align);
					gap = static_cast<size_t>(reinterpret_cast<tlsfptr_t>(aligned_ptr) - reinterpret_cast<tlsfptr_t>(ptr));
				}

				if( gap )
					block = block_trim_free_leading(control, block, gap);
			}

			return block_prepare_used(control, block, adjust);
		}

		// -----------------------------------------------------------------------
		/// @brief Get the usable size of a TLSF allocation
		///
		/// @param ptr Pointer to user data
		/// @return Size of the allocated block (excluding header)
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static size_t tlsf_block_size(void * ptr) {
			block_header_t * block = block_from_ptr(ptr);
			return block_size(block);
		}

		// -----------------------------------------------------------------------
		/// @brief TLSF free implementation
		///
		/// @details Marks the block as free, then coalesces with adjacent
		///   free blocks (both previous and next) before re-inserting into
		///   the free lists.
		///
		/// @param tlsf TLSF handle
		/// @param ptr Pointer to user data to free
		///
		/// @complexity O(1)
		// -----------------------------------------------------------------------
		static void tlsf_free(tlsf_t tlsf, void * ptr) {
			if( ptr != nullptr ) {
				control_t * control = reinterpret_cast<control_t *>(tlsf);
				block_header_t * block = block_from_ptr(ptr);
				block_mark_as_free(block);
				block = block_merge_prev(control, block);
				block = block_merge_next(control, block);
				block_insert(control, block);
			}
		}

		// -----------------------------------------------------------------------
		/// @brief TLSF realloc implementation
		///
		/// @details Attempts to grow/shrink the block in-place by merging
		///   with a free next block. If in-place growth is not possible,
		///   allocates a new block, copies data, and frees the old block.
		///
		/// @param tlsf TLSF handle
		/// @param ptr Pointer to previously allocated user data
		/// @param size New requested size
		/// @param pcursize Optional output pointer for current block size
		/// @return User data pointer (possibly moved), or nullptr on failure
		///
		/// @complexity O(1) best-case (in-place), O(n) worst-case (copy)
		// -----------------------------------------------------------------------
		static void * tlsf_realloc(tlsf_t tlsf, void* ptr, size_t size, size_t * pcursize) {
			control_t * control = reinterpret_cast<control_t *>(tlsf);
			void * p = nullptr;

			if( ptr != nullptr && size == 0 ) {
				tlsf_free(tlsf, ptr);
			}
			else if( ptr == nullptr ) {
				p = tlsf_malloc(tlsf, size);
			}
			else {
				block_header_t * block = block_from_ptr(ptr);
				block_header_t * next = block_next(block);

				const size_t cursize = block_size(block);
				if( pcursize != nullptr )
					*pcursize = cursize;
				const size_t combined = cursize + block_size(next) + block_header_overhead();
				const size_t adjust = adjust_request_size(size, ALIGN_SIZE);

				if( adjust > cursize && (!block_is_free(next) || adjust > combined) ) {
					p = tlsf_malloc(tlsf, size);
					if( p != nullptr ) {
						const size_t minsize = tlsf_min(cursize, size);
						memcpy(p, ptr, minsize);
						tlsf_free(tlsf, ptr);
					}
				}
				else {
					if( adjust > cursize ) {
						block_merge_next(control, block);
						block_mark_as_used(block);
					}

					block_trim_used(control, block, adjust);
					p = ptr;
				}
			}

			return p;
		}
};
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief TLSF_Impl::malloc — Public allocation entry point
///
/// @details
///   The allocation path operates as follows:
///   1. If no TLSF instance exists yet, create one with a default pool
///      (64 MiB via VirtualAlloc on Windows, ::malloc on other platforms).
///   2. Attempt TLSF allocation via tlsf_malloc.
///   3. If TLSF allocation fails (pool exhausted), add a new pool and retry.
///   4. If still fails, fall back to OS-level allocation (VirtualAlloc or
///      ::malloc with hidden size header).
///   5. Increment reference count and optionally fill with debug pattern.
///
/// @param sz Requested allocation size
/// @return Pointer to allocated memory, or nullptr on failure
///
/// @complexity O(1) typical, O(FL_INDEX_COUNT * SL_INDEX_COUNT) on first call
///
/// @post ref_count_ incremented
/// @post errno = ENOMEM on failure
// ---------------------------------------------------------------------------
inline void * TLSF_Impl::malloc(size_t sz)
{
	void * p = nullptr;
	void * mem = nullptr;
	const bool is_tlsf_size = sz <= tlsf_pool_max_block();

	if( tlsf_ == nullptr && is_tlsf_size ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
		mem = VirtualAlloc(nullptr, pool_size(), MEM_COMMIT, PAGE_READWRITE);
#else
		mem = ::malloc(pool_size());
#endif

		if( mem == nullptr ){
			errno = ENOMEM;
			return nullptr;
		}

		tlsf_ = tlsf_create_with_pool(mem, pool_size());

		if( tlsf_ == nullptr ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
			VirtualFree(mem, 0, MEM_RELEASE);
#else
			::free(mem);
#endif
			errno = ENOMEM;
			return nullptr;
		}
	}

	if( p == nullptr && is_tlsf_size )
		p = tlsf_malloc(tlsf_,sz);

	// try add new pool
	if( p == nullptr && is_tlsf_size ) {
#if defined(_INC_WINDOWS) && defined(_WIN32)
		//SYSTEM_INFO si;
		//GetNativeSystemInfo(&si);
		//const size_t page_size = si.dwPageSize;
		const size_t page_size = 65536;
#else
		const size_t page_size = 4096;
#endif

		uintptr_t pools = pools_ == nullptr ? 0 : uintptr_t(pools_[0]);
		pool_t * p_pool = nullptr;
		uintptr_t csz = ((pools + 1) * sizeof(pools_[0])), nsz = csz + sizeof(pools_[0]);
		uintptr_t cp = (csz % page_size != 0 ? 1 : 0), np = (nsz % page_size != 0 ? 1 : 0);
		uintptr_t pc = csz / page_size + cp, npc = nsz / page_size + np;

		if( pools_ == nullptr || pc != npc )
#if defined(_INC_WINDOWS) && defined(_WIN32)
			p_pool = reinterpret_cast<pool_t *>(VirtualAlloc(nullptr, npc * page_size, MEM_COMMIT, PAGE_READWRITE));
#else
			p_pool = reinterpret_cast<pool_t *>(::malloc(npc * page_size));
#endif
		else
			p_pool = pools_;

		if( p_pool == nullptr ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
			if( mem != nullptr )
				VirtualFree(mem, 0, MEM_RELEASE);
#else
			::free(mem);
#endif
			errno = ENOMEM;
        }
        else {
			if( p_pool != pools_ ){
				if( pools_ != nullptr )
					memcpy(p_pool,pools_,csz);
#if defined(_INC_WINDOWS) && defined(_WIN32)
				VirtualFree(pools_, 0, MEM_RELEASE);
#else
				::free(pools_);
#endif
				pools_ = p_pool;
			}

#if defined(_INC_WINDOWS) && defined(_WIN32)
			mem = VirtualAlloc(nullptr, pool_size(), MEM_COMMIT, PAGE_READWRITE);
#else
			mem = ::malloc(pool_size());
#endif
			if( mem == nullptr ){
				errno = ENOMEM;
				return nullptr;
			}

			pool_t pool = tlsf_add_pool(tlsf_,mem,pool_size());

			if( pool == nullptr ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
				VirtualFree(mem, 0, MEM_RELEASE);
#else
				::free(mem);
#endif
				errno = ENOMEM;
				return nullptr;
			}

			pools_[0] = reinterpret_cast<pool_t *>(++pools);
			pools_[pools] = pool;
			qsort(pools_,1,pools);
			p = tlsf_malloc(tlsf_,sz);
		}
	}

	if( p == nullptr ) {
#if defined(_INC_WINDOWS) && defined(_WIN32)
		p = VirtualAlloc(nullptr, sz, MEM_COMMIT, PAGE_READWRITE);
#else
		p = ::malloc(sz + sizeof(size_t));
		*reinterpret_cast<size_t *>(p) = sz;
		p = reinterpret_cast<uint8_t *>(p) + sizeof(size_t);
#endif
	}

	if( p == nullptr ) {
		errno = ENOMEM;
	}
	else {
		ref_count_++;
#if _DEBUG
		memset(p, 0xCC, sz);
#endif
	}

	return p;
}
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief TLSF_Impl::realloc — Public reallocation entry point
///
/// @details Reallocates memory, preserving contents up to the minimum of
///   old and new sizes. Attempts TLSF realloc for pool-managed pointers;
///   falls back to malloc+copy+free for OS-level allocations.
///
/// @param pp Previously allocated pointer (may be nullptr → malloc behavior)
/// @param sz New size (may be 0 → free behavior)
/// @return Pointer to reallocated memory, or nullptr on failure
///
/// @complexity O(1) best-case, O(n) worst-case
/// @post ref_count_ adjusted as needed
// ---------------------------------------------------------------------------
inline void * TLSF_Impl::realloc(void * pp, size_t sz)
{
	void * p = nullptr;

	if( pp == nullptr ){
		p = this->malloc(sz);
	}
	else if( sz == 0 ){
		this->free(pp);
		p = nullptr;
	}
	else {
		bool in_pool = reinterpret_cast<uintptr_t>(pp) >= reinterpret_cast<uintptr_t>(tlsf_) && reinterpret_cast<uintptr_t>(pp) < reinterpret_cast<uintptr_t>(tlsf_) + pool_size();
		uintptr_t pools = pools_ == nullptr ? 0 : uintptr_t(pools_[0]);

		in_pool = in_pool ||
			bsearch(pools_, 1, pools, pp, [] (const pool_t & key,const pool_t & b) {
				return reinterpret_cast<uintptr_t>(key) >= reinterpret_cast<uintptr_t>(b) + pool_size() ? 1 : key < b ? -1 : 0;
			}) >= 0;

		size_t cursize = 0;

		if( in_pool )
			p = tlsf_realloc(tlsf_, pp, sz, &cursize);

		if( p == nullptr ){
			p = this->malloc(sz);

			if( p != nullptr ){
				if( cursize == 0 ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
					MEMORY_BASIC_INFORMATION info;
					VirtualQuery(pp, &info, sizeof(info));
					cursize = info.RegionSize;
#else
					cursize = reinterpret_cast<size_t *>(p)[-1];
#endif
				}
				memcpy(p, pp, cursize > sz ? sz : cursize);
				this->free(pp);
			}
		}
	}

	return p;
}
//------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
/// @brief TLSF_Impl::free — Public deallocation entry point
///
/// @details
///   Determines whether the pointer belongs to a TLSF pool or was
///   allocated directly from the OS, then frees accordingly. When the
///   reference count drops to zero, all pools and the TLSF instance
///   are destroyed.
///
/// @param p Pointer to previously allocated memory
///
/// @complexity O(1) typical, O(n) during final cleanup
/// @post ref_count_ decremented (may reach 0)
/// @post errno = 0
// ---------------------------------------------------------------------------
inline void TLSF_Impl::free(void * p)
{
	if( p != nullptr ) {
		bool in_pool = reinterpret_cast<uintptr_t>(p) >= reinterpret_cast<uintptr_t>(tlsf_) && reinterpret_cast<uintptr_t>(p) < reinterpret_cast<uintptr_t>(tlsf_) + pool_size();
		uintptr_t pools = pools_ == nullptr ? 0 : uintptr_t(pools_[0]);

		in_pool = in_pool ||
			bsearch(pools_, 1, pools, p, [](const pool_t & key, const pool_t & b) {
				return reinterpret_cast<uintptr_t>(key) >= reinterpret_cast<uintptr_t>(b) + pool_size() ? 1 : key < b ? -1 : 0;
			}) >= 0;


		if( in_pool ){
#if _DEBUG
			memset(p, 0xCD, tlsf_block_size(p));
#endif
			tlsf_free(tlsf_, p);
		}
		else {
#if defined(_INC_WINDOWS) && defined(_WIN32)
#if _DEBUG
			MEMORY_BASIC_INFORMATION info;
			VirtualQuery(p, &info, sizeof(info));
			memset(p, 0xCD, info.RegionSize);
#endif
			VirtualFree(p, 0, MEM_RELEASE);
#else
			p = reinterpret_cast<uint8_t *>(p) - sizeof(size_t);
			::free(p);
#endif
		}

		if( --ref_count_ == 0 && tlsf_ != nullptr ){
			while( pools > 0 ){
				tlsf_remove_pool(tlsf_, pools_[pools]);
#if defined(_INC_WINDOWS) && defined(_WIN32)
				VirtualFree(pools_[pools], 0, MEM_RELEASE);
#else
				::free(pools_[pools]);
#endif
				pools--;
			}

#if defined(_INC_WINDOWS) && defined(_WIN32)
			if( pools_ != nullptr )
				VirtualFree(pools_, 0, MEM_RELEASE);
#else
			::free(pools_);
#endif
			pools_ = nullptr;
			tlsf_destroy(tlsf_);

#if defined(_INC_WINDOWS) && defined(_WIN32)
			VirtualFree(tlsf_, 0, MEM_RELEASE);
#else
			::free(tlsf_);
#endif
			tlsf_ = nullptr;
		}
	}

	errno = 0;
}
//------------------------------------------------------------------------------
} // namespace tlsf - Two Level Segregated Fit memory allocator
//------------------------------------------------------------------------------
#endif // TLSF_HPP_INCLUDED
//------------------------------------------------------------------------------
