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
#elif _WIN32
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
namespace tlsf {  // Two Level Segregated Fit memory allocator
//------------------------------------------------------------------------------
/*
** Two Level Segregated Fit memory allocator, version 3.0.
** Written by Matthew Conte, and placed in the Public Domain.
**	http://tlsf.baisoku.org
**
** Based on the original documentation by Miguel Masmano:
**	http://rtportal.upv.es/rtmalloc/allocators/tlsf/index.shtml
**
** Please see the accompanying Readme.txt for implementation
** notes and caveats.
**
** This implementation was written to the specification
** of the document, therefore no GPL restrictions apply.
*/
/*
** Architecture-specific bit manipulation routines.
**
** TLSF achieves O(1) cost for malloc and free operations by limiting
** the search for a free block to a free list of guaranteed size
** adequate to fulfill the request, combined with efficient free list
** queries using bitmasks and architecture-specific bit-manipulation
** routines.
**
** Most modern processors provide instructions to count leading zeroes
** in a word, find the lowest and highest set bit, etc. These
** specific implementations will be used when available, falling back
** to a reasonably efficient generic implementation.
**
** NOTE: TLSF spec relies on ffs/fls returning value 0..31.
** ffs/fls return 1-32 by default, returning 0 for error.
*/

/*
** Detect whether or not we are building for a 32- or 64-bit (LP/LLP)
** architecture. There is no reliable portable method at compile-time.
*/
#if defined (__alpha__) || defined (__ia64__) || defined (__x86_64__) \
	|| defined (_WIN64) || defined (__LP64__) || defined (__LLP64__)
#define TLSF_64BIT
#endif

/*
** gcc 3.4 and above have builtin support, specialized for architecture.
** Some compilers masquerade as gcc; patchlevel test filters them out.
*/
#if defined (__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)) \
	&& defined (__GNUC_PATCHLEVEL__)

inline int tlsf_ffs(unsigned int word)
{
	return __builtin_ffs(word) - 1;
}

inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __builtin_clz(word) : 0;
	return bit - 1;
}

#elif defined (_MSC_VER) && (_MSC_VER >= 1400) && (defined (_M_IX86) || defined (_M_X64))
/* Microsoft Visual C++ support on x86/X64 architectures. */

#include <intrin.h>

#pragma intrinsic(_BitScanReverse)
#pragma intrinsic(_BitScanForward)

inline int tlsf_fls(unsigned int word)
{
	unsigned long index;
	return _BitScanReverse(&index, word) ? index : -1;
}

inline int tlsf_ffs(unsigned int word)
{
	unsigned long index;
	return _BitScanForward(&index, word) ? index : -1;
}

#elif defined (_MSC_VER) && defined (_M_PPC)
/* Microsoft Visual C++ support on PowerPC architectures. */

#include <ppcintrinsics.h>

inline int tlsf_fls(unsigned int word)
{
	const int bit = 32 - _CountLeadingZeros(word);
	return bit - 1;
}

inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - _CountLeadingZeros(reverse);
	return bit - 1;
}

#elif defined (__ARMCC_VERSION)
/* RealView Compilation Tools for ARM */

inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - __clz(reverse);
	return bit - 1;
}

inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __clz(word) : 0;
	return bit - 1;
}

#elif defined (__ghs__)
/* Green Hills support for PowerPC */

#include <ppc_ghs.h>

inline int tlsf_ffs(unsigned int word)
{
	const unsigned int reverse = word & (~word + 1);
	const int bit = 32 - __CLZ32(reverse);
	return bit - 1;
}

inline int tlsf_fls(unsigned int word)
{
	const int bit = word ? 32 - __CLZ32(word) : 0;
	return bit - 1;
}

#else
/* Fall back to generic implementation. */

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

/* Implement ffs in terms of fls. */
inline int tlsf_ffs(unsigned int word)
{
	return tlsf_fls_generic(word & (~word + 1)) - 1;
}

inline int tlsf_fls(unsigned int word)
{
	return tlsf_fls_generic(word) - 1;
}

#endif

/* Possibly 64-bit version of tlsf_fls. */
#if defined (TLSF_64BIT)
inline int tlsf_fls_sizet(size_t size)
{
	unsigned int high = (unsigned int) (size >> 32);
	int bits = 0;

	if( high ) {
		bits = 32 + tlsf_fls(high);
	}
	else {
		bits = tlsf_fls((unsigned int) size & 0xffffffff);
	}

	return bits;
}
#else
inline int tlsf_fls_sizet(size_t size){ return tlsf_fls(size); }
#endif
//------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
enum tlsf_public {
	SL_INDEX_COUNT_LOG2 = 5,
};
//------------------------------------------------------------------------------
enum tlsf_private {
#if defined (TLSF_64BIT)
	ALIGN_SIZE_LOG2 = 3,
#else
	ALIGN_SIZE_LOG2 = 2,
#endif
	ALIGN_SIZE = (1 << ALIGN_SIZE_LOG2),

#if defined (TLSF_64BIT)
	FL_INDEX_MAX = 32,
#else
	FL_INDEX_MAX = 30,
#endif
	SL_INDEX_COUNT = (1 << SL_INDEX_COUNT_LOG2),
	FL_INDEX_SHIFT = (SL_INDEX_COUNT_LOG2 + ALIGN_SIZE_LOG2),
	FL_INDEX_COUNT = (FL_INDEX_MAX - FL_INDEX_SHIFT + 1),

	SMALL_BLOCK_SIZE = (1 << FL_INDEX_SHIFT),
};
//------------------------------------------------------------------------------
typedef struct block_header_t {
	struct block_header_t * prev_phys_block;
	size_t size;
	struct block_header_t * next_free;
	struct block_header_t * prev_free;
} block_header_t;
//------------------------------------------------------------------------------
typedef struct control_t {
	block_header_t block_null;
	unsigned int fl_bitmap;
	unsigned int sl_bitmap[FL_INDEX_COUNT];
	block_header_t * blocks[FL_INDEX_COUNT][SL_INDEX_COUNT];
} control_t;
//------------------------------------------------------------------------------
typedef ptrdiff_t tlsfptr_t;
typedef void * tlsf_t;
typedef void * pool_t;
//------------------------------------------------------------------------------
class TLSF_Impl { // must be singleton
    private:
        TLSF_Impl(const TLSF_Impl &);
        void operator = (const TLSF_Impl &);
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

		void * malloc(size_t sz);
		void * realloc(void * p, size_t sz);
		void free(void * p);
	protected:
		uintptr_t ref_count_;
		pool_t * pools_ = nullptr;
		tlsf_t tlsf_ = nullptr;

		template <typename T,typename F>
		static void qsort(T * p, intptr_t lb, intptr_t ub, F f) {
			const size_t MAX_LEVELS = sizeof(void *) * CHAR_BIT;
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

		template <typename T>
		static void qsort(T * p, intptr_t lb, intptr_t ub)	{
			qsort<T>(p,lb,ub,[] (const T & a,const T & b) {
				return a > b ? 1 : a < b ? -1 : 0;
			});
		}

		template <typename T,typename F>
		static intptr_t bsearch(T * keys,intptr_t low, intptr_t high, const T & key, F f) {
			for(;;){
				uintptr_t p = (low + high) / 2;

				if( low > high ) break;

				intptr_t c = f(key,keys[p]);

				if( c > 0 ){
					low = p + 1;
				}
				else if( c < 0 ){
					high = p - 1;
				}
				else
					return p;
			}

			return -1;
		}

		template <typename T>
		static intptr_t bsearch(T * keys,intptr_t low, intptr_t high, const T & key) {
			return bsearch<T>(keys,low,high,key,[] (const T & key,const T & b) {
				return key > b ? 1 : key < b ? -1 : 0;
			});
		}

	private:
		template <typename T> static const T & tlsf_min(const T & a,const T & b) { return a < b ? a : b; }
		template <typename T> static const T & tlsf_max(const T & a,const T & b) { return a > b ? a : b; }

		static constexpr size_t pool_size()						{ return size_t(64) * 1024u * 1024u; }
		static constexpr size_t block_header_free_bit()			{ return size_t(1) << 0; }
		static constexpr size_t block_header_prev_free_bit()	{ return size_t(1) << 1; }
		static constexpr size_t block_header_overhead()			{ return sizeof(size_t); }
		static constexpr size_t block_start_offset()			{ return offsetof(block_header_t, size) + sizeof(size_t); }
		static constexpr size_t block_size_min()				{ return sizeof(block_header_t) - sizeof(block_header_t *); }
		static constexpr size_t block_size_max()				{ return size_t(1) << FL_INDEX_MAX; }

		static size_t block_size(const block_header_t * block) {
			return block->size & ~(block_header_free_bit() | block_header_prev_free_bit());
		}

		static void block_set_size(block_header_t* block, size_t size) {
			const size_t oldsize = block->size;
			block->size = size | (oldsize & (block_header_free_bit() | block_header_prev_free_bit()));
		}

		static int block_is_last(const block_header_t * block) {
			return 0 == block_size(block);
		}

		static int block_is_free(const block_header_t * block) {
			return static_cast<int>(block->size & block_header_free_bit());
		}

		static void block_set_free(block_header_t * block) {
			block->size |= block_header_free_bit();
		}

		static void block_set_used(block_header_t* block) {
			block->size &= ~block_header_free_bit();
		}

		static int block_is_prev_free(const block_header_t * block)	{
			return static_cast<int>(block->size & block_header_prev_free_bit());
		}

		static void block_set_prev_free(block_header_t* block) {
			block->size |= block_header_prev_free_bit();
		}

		static void block_set_prev_used(block_header_t* block) {
			block->size &= ~block_header_prev_free_bit();
		}

		static block_header_t * block_from_ptr(const void * ptr) {
			return (block_header_t *) ((unsigned char *) ptr - block_start_offset());
		}

		static void * block_to_ptr(const block_header_t * block) {
			return (void *) ((unsigned char *) block + block_start_offset());
		}

		static block_header_t * offset_to_block(const void * ptr, size_t size) {
			return (block_header_t *) ((tlsfptr_t) ptr + size);
		}

		static block_header_t * block_prev(const block_header_t * block) {
			return block->prev_phys_block;
		}

		static block_header_t * block_next(const block_header_t * block) {
			block_header_t * next = offset_to_block(block_to_ptr(block),
				block_size(block) - block_header_overhead());
			return next;
		}

		static block_header_t * block_link_next(block_header_t * block) {
			block_header_t * next = block_next(block);
			next->prev_phys_block = block;
			return next;
		}

		static void block_mark_as_free(block_header_t* block) {
			block_header_t * next = block_link_next(block);
			block_set_prev_free(next);
			block_set_free(block);
		}

		static void block_mark_as_used(block_header_t * block) {
			block_header_t * next = block_next(block);
			block_set_prev_used(next);
			block_set_used(block);
		}

		static size_t align_up(size_t x, size_t align) {
			return (x + (align - 1)) & ~(align - 1);
		}

		static size_t align_down(size_t x, size_t align) {
			return x - (x & (align - 1));
		}

		static void * align_ptr(const void* ptr, size_t align) {
			return (void *) (((tlsfptr_t) ptr + (align - 1)) & ~(align - 1));
		}

		static size_t adjust_request_size(size_t size, size_t align) {
			size_t adjust = 0;

			if( size && size < block_size_max() )
				adjust = tlsf_max(align_up(size, align), block_size_min());

			return adjust;
		}

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

		static void mapping_search(size_t size, int * fli, int * sli) {
			if( size >= (1 << SL_INDEX_COUNT_LOG2) ) {
				const size_t round = (1 << (tlsf_fls_sizet(size) - SL_INDEX_COUNT_LOG2)) - 1;
				size += round;
			}
			mapping_insert(size, fli, sli);
		}

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

		static void remove_free_block(control_t * control, block_header_t * block, int fl, int sl) {
			block_header_t* prev = block->prev_free;
			block_header_t* next = block->next_free;
			next->prev_free = prev;
			prev->next_free = next;

			if( control->blocks[fl][sl] == block ) {
				control->blocks[fl][sl] = next;

				if( next == &control->block_null ) {
					control->sl_bitmap[fl] &= ~(1 << sl);
					if( !control->sl_bitmap[fl])
						control->fl_bitmap &= ~(1 << fl);
				}
			}
		}

		static void insert_free_block(control_t * control, block_header_t * block, int fl, int sl) {
			block_header_t * current = control->blocks[fl][sl];
			block->next_free = current;
			block->prev_free = &control->block_null;
			current->prev_free = block;
			control->blocks[fl][sl] = block;
			control->fl_bitmap |= (1 << fl);
			control->sl_bitmap[fl] |= (1 << sl);
		}

		static void block_remove(control_t * control, block_header_t * block) {
			int fl, sl;
			mapping_insert(block_size(block), &fl, &sl);
			remove_free_block(control, block, fl, sl);
		}

		static void block_insert(control_t * control, block_header_t * block) {
			int fl, sl;
			mapping_insert(block_size(block), &fl, &sl);
			insert_free_block(control, block, fl, sl);
		}

		static int block_can_split(block_header_t * block, size_t size) {
			return block_size(block) >= sizeof(block_header_t) + size;
		}

		static block_header_t * block_split(block_header_t* block, size_t size) {
			block_header_t * remaining =
				offset_to_block(block_to_ptr(block), size - block_header_overhead());

			const size_t remain_size = block_size(block) - (size + block_header_overhead());

			block_set_size(remaining, remain_size);

			block_set_size(block, size);
			block_mark_as_free(remaining);

			return remaining;
		}

		static block_header_t * block_absorb(block_header_t * prev, block_header_t * block) {
			prev->size += block_size(block) + block_header_overhead();
			block_link_next(prev);
			return prev;
		}

		static block_header_t * block_merge_prev(control_t * control, block_header_t * block) {
			if( block_is_prev_free(block) )	{
				block_header_t* prev = block_prev(block);
				block_remove(control, prev);
				block = block_absorb(prev, block);
			}
			return block;
		}

		static block_header_t * block_merge_next(control_t * control, block_header_t * block) {
			block_header_t * next = block_next(block);

			if( block_is_free(next) ) {
				block_remove(control, next);
				block = block_absorb(block, next);
			}

			return block;
		}

		static void block_trim_free(control_t * control, block_header_t * block, size_t size) {
			if( block_can_split(block, size) ) {
				block_header_t* remaining_block = block_split(block, size);
				block_link_next(block);
				block_set_prev_free(remaining_block);
				block_insert(control, remaining_block);
			}
		}

		static void block_trim_used(control_t * control, block_header_t * block, size_t size) {
			if( block_can_split(block, size) ) {
				block_header_t* remaining_block = block_split(block, size);
				block_set_prev_used(remaining_block);
				remaining_block = block_merge_next(control, remaining_block);
				block_insert(control, remaining_block);
			}
		}

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

		static void * block_prepare_used(control_t * control, block_header_t * block, size_t size) {
			void * p = nullptr;

			if( block != nullptr ) {
				block_trim_free(control, block, size);
				block_mark_as_used(block);
				p = block_to_ptr(block);
			}

			return p;
		}

		static void control_construct(control_t * control) {
			int i, j;

			control->block_null.next_free = &control->block_null;
			control->block_null.prev_free = &control->block_null;

			control->fl_bitmap = 0;

			for( i = 0; i < FL_INDEX_COUNT; ++i ) {
				control->sl_bitmap[i] = 0;

				for( j = 0; j < SL_INDEX_COUNT; ++j )
					control->blocks[i][j] = &control->block_null;
			}
		}

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

		static pool_t tlsf_add_pool(tlsf_t tlsf, void * mem, size_t bytes) {
			block_header_t * block;
			block_header_t * next;

			const size_t pool_overhead = tlsf_pool_overhead();
			const size_t pool_bytes = align_down(bytes - pool_overhead, ALIGN_SIZE);

			if( (ptrdiff_t) mem % ALIGN_SIZE != 0)
				return nullptr;

			if( pool_bytes < block_size_min() || pool_bytes > block_size_max() )
				return nullptr;

			block = offset_to_block(mem, -(tlsfptr_t) block_header_overhead());
			block_set_size(block, pool_bytes);
			block_set_free(block);
			block_set_prev_used(block);
			block_insert((control_t *) tlsf, block);

			next = block_link_next(block);
			block_set_size(next, 0);
			block_set_used(next);
			block_set_prev_free(next);

			return mem;
		}

		static void tlsf_remove_pool(tlsf_t tlsf, pool_t pool) {
			control_t * control = (control_t *) tlsf;
			block_header_t * block = offset_to_block(pool, -(int) block_header_overhead());

			int fl = 0, sl = 0;

			mapping_insert(block_size(block), &fl, &sl);
			remove_free_block(control, block, fl, sl);
		}

		static tlsf_t tlsf_create(void * mem) {
			if( (tlsfptr_t) mem % ALIGN_SIZE != 0 )
				return nullptr;

			control_construct((control_t *) mem);

			return (tlsf_t) mem;
		}

		static tlsf_t tlsf_create_with_pool(void * mem, size_t bytes) {
			tlsf_t tlsf = tlsf_create(mem);
			tlsf_add_pool(tlsf, (char*)mem + tlsf_size(), bytes - tlsf_size());
			return tlsf;
		}

		static void tlsf_destroy(tlsf_t tlsf) {
			(void)tlsf;
		}

		pool_t tlsf_get_pool(tlsf_t tlsf) {
			return (pool_t) ((char *) tlsf + tlsf_size());
		}

		static void * tlsf_malloc(tlsf_t tlsf, size_t size) {
			control_t * control = (control_t *) tlsf;
			const size_t adjust = adjust_request_size(size, ALIGN_SIZE);
			block_header_t * block = block_locate_free(control, adjust);
			return block_prepare_used(control, block, adjust);
		}

		static void * tlsf_memalign(tlsf_t tlsf, size_t align, size_t size) {
			control_t * control = (control_t *) tlsf;
			const size_t adjust = adjust_request_size(size, ALIGN_SIZE);
			const size_t gap_minimum = sizeof(block_header_t);
			const size_t size_with_gap = adjust_request_size(adjust + align + gap_minimum, align);
			const size_t aligned_size = (align <= ALIGN_SIZE) ? adjust : size_with_gap;
			block_header_t * block = block_locate_free(control, aligned_size);

			if( block != nullptr ) {
				void * ptr = block_to_ptr(block);
				void * aligned_ptr = align_ptr(ptr, align);
				size_t gap = (size_t) ((tlsfptr_t) aligned_ptr - (tlsfptr_t) ptr);

				if( gap && gap < gap_minimum ) {
					const size_t gap_remain = gap_minimum - gap;
					const size_t offset = tlsf_max(gap_remain, align);
					const void * next_aligned = (void *) ((tlsfptr_t) aligned_ptr + offset);

					aligned_ptr = align_ptr(next_aligned, align);
					gap = (size_t) ((tlsfptr_t) aligned_ptr - (tlsfptr_t) ptr);
				}

				if( gap )
					block = block_trim_free_leading(control, block, gap);
			}

			return block_prepare_used(control, block, adjust);
		}

		static size_t tlsf_block_size(void * ptr) {
			block_header_t * block = block_from_ptr(ptr);
			return block_size(block);
		}

		static void tlsf_free(tlsf_t tlsf, void * ptr) {
			if( ptr != nullptr ) {
				control_t * control = (control_t *) tlsf;
				block_header_t * block = block_from_ptr(ptr);
				block_mark_as_free(block);
				block = block_merge_prev(control, block);
				block = block_merge_next(control, block);
				block_insert(control, block);
			}
		}

		static void * tlsf_realloc(tlsf_t tlsf, void* ptr, size_t size, size_t * pcursize) {
			control_t * control = (control_t *) tlsf;
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
inline void * TLSF_Impl::malloc(size_t sz)
{
	void * p = nullptr;
	void * mem = nullptr;
	const bool is_tlsf_size = sz <= tlsf_pool_max_block();

	if( tlsf_ == nullptr && is_tlsf_size ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
		mem = VirtualAlloc(NULL, pool_size(), MEM_COMMIT, PAGE_READWRITE);
#else
		mem = ::malloc(pool_size());
#endif

		if( mem == NULL ){
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
		pool_t * p_pool = NULL;
		uintptr_t csz = ((pools + 1) * sizeof(pools_[0])), nsz = csz + sizeof(pools_[0]);
		uintptr_t cp = (csz % page_size != 0 ? 1 : 0), np = (nsz % page_size != 0 ? 1 : 0);
		uintptr_t pc = csz / page_size + cp, npc = nsz / page_size + np;

		if( pools_ == nullptr || pc != npc )
#if defined(_INC_WINDOWS) && defined(_WIN32)
			p_pool = (pool_t *) VirtualAlloc(NULL, npc * page_size, MEM_COMMIT, PAGE_READWRITE);
#else
			p_pool = (pool_t *) ::malloc(npc * page_size);
#endif
		else
			p_pool = pools_;

		if( p_pool == NULL ){
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
			mem = VirtualAlloc(NULL, pool_size(), MEM_COMMIT, PAGE_READWRITE);
#else
			mem = ::malloc(pool_size());
#endif
			if( mem == NULL ){
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

			pools_[0] = (pool_t *) (++pools);
			pools_[pools] = pool;
			qsort(pools_,1,pools);
			p = tlsf_malloc(tlsf_,sz);
		}
	}

	if( p == nullptr ) {
#if defined(_INC_WINDOWS) && defined(_WIN32)
		p = VirtualAlloc(NULL, sz, MEM_COMMIT, PAGE_READWRITE);
#else
		p = ::malloc(sz + sizeof(size_t));
		*(size_t *) p = sz;
		p = (uint8_t *) p + sizeof(size_t);
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
		bool in_pool = (uintptr_t) pp >= (uintptr_t) tlsf_ && (uintptr_t) pp < (uintptr_t) tlsf_ + pool_size();
		uintptr_t pools = pools_ == nullptr ? 0 : uintptr_t(pools_[0]);

		in_pool = in_pool ||
			bsearch(pools_, 1, pools, pp, [] (const pool_t & key,const pool_t & b) {
				return (uintptr_t) key >= (uintptr_t) b + pool_size() ? 1 : key < b ? -1 : 0;
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
					cursize = ((size_t *) p)[-1];
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
inline void TLSF_Impl::free(void * p)
{
	if( p != nullptr ) {
		bool in_pool = (uintptr_t) p >= (uintptr_t) tlsf_ && (uintptr_t) p < (uintptr_t) tlsf_ + pool_size();
		uintptr_t pools = pools_ == nullptr ? 0 : uintptr_t(pools_[0]);

		in_pool = in_pool ||
			bsearch(pools_, 1, pools, p, [](const pool_t & key, const pool_t & b) {
				return (uintptr_t) key >= (uintptr_t) b + pool_size() ? 1 : key < b ? -1 : 0;
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
			p = (uint8_t *) p - sizeof(size_t);
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
/*extern "C" void __cdecl free(void *_Memory);
extern "C" void * __cdecl calloc(size_t _NumOfElements,size_t _SizeOfElements);
extern "C" void * __cdecl malloc(size_t _Size);
extern "C" void * __cdecl realloc(void *_Memory,size_t _NewSize);
extern "C" void * __cdecl _recalloc(void *_Memory,size_t _Count,size_t _Size);
// Make sure that X86intrin.h doesn't produce here collisions.
#if (!defined (_XMMINTRIN_H_INCLUDED) && !defined (_MM_MALLOC_H_INCLUDED)) || defined(_aligned_malloc)
#pragma push_macro("_aligned_free")
#pragma push_macro("_aligned_malloc")
#undef _aligned_free
#undef _aligned_malloc
extern "C" void __cdecl _aligned_free(void *_Memory);
extern "C" void * __cdecl _aligned_malloc(size_t _Size,size_t _Alignment);
#pragma pop_macro("_aligned_free")
#pragma pop_macro("_aligned_malloc")
#endif
extern "C" void * __cdecl _aligned_offset_malloc(size_t _Size,size_t _Alignment,size_t _Offset);
extern "C" void * __cdecl _aligned_realloc(void *_Memory,size_t _Size,size_t _Alignment);
extern "C" void * __cdecl _aligned_recalloc(void *_Memory,size_t _Count,size_t _Size,size_t _Alignment);
extern "C" void * __cdecl _aligned_offset_realloc(void *_Memory,size_t _Size,size_t _Alignment,size_t _Offset);
extern "C" void * __cdecl _aligned_offset_recalloc(void *_Memory,size_t _Count,size_t _Size,size_t _Alignment,size_t _Offset);*/
//------------------------------------------------------------------------------
#endif // TLSF_HPP_INCLUDED
//------------------------------------------------------------------------------
