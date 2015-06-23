/*-
 * The MIT License (MIT)
 * 
 * Copyright (c) 2014, 2015 Guram Duka
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
#ifndef NN_HPP_INCLUDED
#define NN_HPP_INCLUDED
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
#if _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4244)
#pragma warning(disable : 4996)
#pragma warning(disable : 4290)
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
	public:
		~TLSF_Impl() {
			assert( ref_count_ == 0 );
			assert( pools_ == nullptr );
			assert( tlsf_ == nullptr );
		}

		TLSF_Impl() {
			assert( ref_count_ == 0 );
			assert( pools_ == nullptr );
			assert( tlsf_ == nullptr );
		}

		void * malloc(size_t sz);
		void * realloc(void * p, size_t sz);
		void free(void * p);
	protected:
		intptr_t ref_count_ = 0;
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

		static const size_t pool_size()						{ return size_t(256) * 1024u * 1024u; }
		static const size_t block_header_free_bit()			{ return size_t(1) << 0; }
		static const size_t block_header_prev_free_bit()	{ return size_t(1) << 1; }
		static const size_t block_header_overhead()			{ return sizeof(size_t); }
		static const size_t block_start_offset()			{ return offsetof(block_header_t, size) + sizeof(size_t); }
		static const size_t block_size_min()				{ return sizeof(block_header_t) - sizeof(block_header_t *); }
		static const size_t block_size_max()				{ return size_t(1) << FL_INDEX_MAX; }

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

			unsigned int sl_map = control->sl_bitmap[fl] & (~0 << sl);

			if( !sl_map ) {
				/* No block exists. Search in the next largest first-level list. */
				const unsigned int fl_map = control->fl_bitmap & (~0 << (fl + 1));

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

		static size_t tlsf_size() {
			return sizeof(control_t);
		}

		static size_t tlsf_align_size() {
			return ALIGN_SIZE;
		}

		static size_t tlsf_block_size_min() {
			return block_size_min();
		}

		static size_t tlsf_block_size_max() {
			return block_size_max();
		}

		static size_t tlsf_pool_overhead() {
			return 2 * block_header_overhead();
		}

		static size_t tlsf_alloc_overhead() {
			return block_header_overhead();
		}

		static size_t tlsf_pool_max_block() {
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
	errno = 0;
	void * p = nullptr;
	void * mem = nullptr;
	const bool is_tlsf_size = sz <= tlsf_pool_max_block();

	if( tlsf_ == nullptr && is_tlsf_size ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
		mem = VirtualAlloc(NULL,pool_size(),MEM_COMMIT,PAGE_READWRITE);
#else
		mem = ::malloc(pool_size());
#endif

		if( mem == NULL ){
			errno = ENOMEM;
			return nullptr;
		}

		tlsf_ = tlsf_create_with_pool(mem,pool_size());

		if( tlsf_ == nullptr ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
			VirtualFree(mem,0,MEM_RELEASE);
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
			p_pool = (pool_t *) VirtualAlloc(NULL,npc * page_size,MEM_COMMIT,PAGE_READWRITE);
#else
			p_pool = (pool_t *) ::malloc(npc * page_size);
#endif
		else
			p_pool = pools_;

		if( p_pool == NULL ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
			VirtualFree(mem,0,MEM_RELEASE);
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
				VirtualFree(pools_,0,MEM_RELEASE);
#else
				::free(pools_);
#endif
				pools_ = p_pool;
			}

#if defined(_INC_WINDOWS) && defined(_WIN32)
			mem = VirtualAlloc(NULL,pool_size(),MEM_COMMIT,PAGE_READWRITE);
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
				VirtualFree(mem,0,MEM_RELEASE);
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
		p = VirtualAlloc(NULL,sz,MEM_COMMIT,PAGE_READWRITE);
#else
		p = ::malloc(sz + sizeof(size_t));
		*(size_t *) p = sz;
		p = (uint8_t *) p + sizeof(size_t);
#endif
	}

	if( p == nullptr )
		errno = ENOMEM;
	else
		ref_count_++;

	return p;
}
//------------------------------------------------------------------------------
inline void * TLSF_Impl::realloc(void * pp,size_t sz)
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
			bsearch(pools_,1,pools,pp,[] (const pool_t & key,const pool_t & b) {
				return (uintptr_t) key >= (uintptr_t) b + pool_size() ? 1 : key < b ? -1 : 0;
			}) >= 0;

		size_t cursize = 0;

		if( in_pool )
			p = tlsf_realloc(tlsf_,pp,sz,&cursize);

		if( p == nullptr ){
			p = this->malloc(sz);

			if( p != nullptr ){
				if( cursize == 0 ){
#if defined(_INC_WINDOWS) && defined(_WIN32)
					MEMORY_BASIC_INFORMATION info;
					VirtualQuery(pp,&info,sizeof(info));
					cursize = info.RegionSize;
#else
					cursize = ((size_t *) p)[-1];
#endif
				}
				memcpy(p,pp,cursize > sz ? sz : cursize);
				this->free(pp);
			}
		}
		else {
			errno = 0;
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
			tlsf_free(tlsf_, p);
		}
		else {
#if defined(_INC_WINDOWS) && defined(_WIN32)
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
TLSF_Impl static_allocator;
//------------------------------------------------------------------------------
} // namespace tlsf - Two Level Segregated Fit memory allocator
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------
template <typename A,typename B>
static inline const A & imax(const A & a,const B & b)
{
	return a > b ? a : b;
}
//------------------------------------------------------------------------------
template <typename A,typename B>
static inline const A & imin(const A & a,const B & b)
{
	return a < b ? a : b;
}
//------------------------------------------------------------------------------
#if SIZEOF_WORD == 1
typedef uint8_t word;
typedef int8_t sword;
typedef uint16_t dword;
typedef int16_t sdword;
#elif SIZEOF_WORD == 2
typedef uint16_t word;
typedef int16_t sword;
typedef uint32_t dword;
typedef int32_t sdword;
#elif SIZEOF_WORD == 4
typedef uint32_t word;
typedef int32_t sword;
typedef uint64_t dword;
typedef int64_t sdword;
#elif SIZEOF_WORD == 8 && __GNUC__ && __x86_64__
typedef uint64_t word;
typedef int64_t sword;
typedef __uint128_t dword;
typedef __int128_t sdword;
#elif __GNUC__ && __x86_64__
typedef uint64_t word;
typedef int64_t sword;
typedef __uint128_t dword;
typedef __int128_t sdword;
#define SIZEOF_WORD 8
#elif _MSC_VER && _M_IX86
typedef uint32_t word;
typedef int32_t sword;
typedef uint64_t dword;
typedef int64_t sdword;
#define SIZEOF_WORD 4
#elif _MSC_VER && _M_X64
typedef uint32_t word;
typedef int32_t sword;
typedef uint64_t dword;
typedef int64_t sdword;
#define SIZEOF_WORD 4
#else
typedef uint32_t word;
typedef int32_t sword;
typedef uint64_t dword;
typedef int64_t sdword;
#define SIZEOF_WORD 4
#endif
//------------------------------------------------------------------------------
class integer;
class numeric;
//------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
typedef struct nn_integer_data {
	mutable intptr_t ref_count_;
	mutable uintptr_t length_;
#if SIZEOF_WORD == 1
	mutable uintptr_t dummy_;		// for shift down overflow bits
	mutable word data_[16];
#elif SIZEOF_WORD == 2
	mutable uintptr_t dummy_;		// for shift down overflow bits
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

	typedef nn_integer_data * nn_integer;

	nn_integer add_ref() const {
	  ref_count_++;
	  return const_cast<nn_integer>(this);
	}

	static uintptr_t get_size(uintptr_t length) {
		nn_integer_data pp;
		return (uintptr_t) pp.data_ - (uintptr_t) &pp + (length + 2) * sizeof(pp.data_[0]);
	}

	static nn_integer nn_new(uintptr_t length) {
		nn_integer p = (nn_integer) tlsf::static_allocator.malloc(get_size(length));
		p->ref_count_ = 1;
		p->length_ = length;
		p->dummy_ = 0;
		return p;
	}

	void release() {
		if( --ref_count_ == 0 )
			tlsf::static_allocator.free(this);
	}

	sword isign() const {
		return ((sword) data_[length_]) >> (sizeof(word) * CHAR_BIT - 1);
	}

	void normalize() {
		while( length_ > 1 && data_[length_] == data_[length_ - 1] )
			length_--;
	}

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

		while( r < e ){
			q = (dword) *d0++ + *d1++ + cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		e = result->data_ + p1->length_;

		while( r < e ){
			q = (dword) s0 + *d1++ + cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 + s1 + cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 + s1 + cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

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

		while( r < e ){
			q = (dword) *d0++ + *d1++ + cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 + s1 + cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 + s1 + cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

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

		while( r < e ){
			q = (dword) *d0++ + *d1++ + cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		e = result->data_ + length_;

		while( r < e ){
			q = (dword) *d0++ + s1 + cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 + s1 + cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 + s1 + cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

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

		while( r < e ){
			q = (dword) *d0++ - *d1++ - cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		e = result->data_ + p1->length_;

		while( r < e ){
			q = (dword) s0 - *d1++ - cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 - s1 - cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 - s1 - cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

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

		while( r < e ){
			q = (dword) *d0++ - *d1++ - cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 - s1 - cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 - s1 - cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

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

		while( r < e ){
			q = (dword) *d0++ - *d1++ - cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		e = result->data_ + length_;

		while( r < e ){
			q = (dword) *d0++ - s1 - cf;
			cf = (q >> sizeof(word) * CHAR_BIT) & 1;
			*r++ = q;
		}

		*r++ = q = (dword) s0 - s1 - cf;
		cf = (q >> sizeof(word) * CHAR_BIT) & 1;
		*r = (dword) s0 - s1 - cf;

		// if overflow then grow size
		result->length_++;
		result->normalize();

		return result;
	}

	nn_integer iadd(const nn_integer v) const {
		intptr_t c = length_ - v->length_;

		if( c > 0 )
			return addp(v);

		if( c < 0 )
			return addm(v);

		return addz(v);
	}

	nn_integer isub(const nn_integer v) const {
		intptr_t c = length_ - v->length_;

		if( c > 0 )
			return subp(v);

		if( c < 0 )
			return subm(v);

		return subz(v);
	}

	nn_integer isal(uintptr_t bit_count) const;
	nn_integer isar(uintptr_t bit_count) const;

	bool z() const { // is zero
		if( data_[0] != 0 )
			return false;

		for( intptr_t i = length_; i > 0; i-- )
			if( data_[i] != 0 )
				return false;

		return true;
	}

	bool nz() const { // iz not zero
		for( intptr_t i = length_; i >= 0; i-- )
			if( data_[i] != 0 )
				return true;

		return false;
	}

	intptr_t icompare(const nn_integer v1) const {
		intptr_t c = 0;
		nn_integer a = isub(v1);

		c = a->isign() < 0 ? -1 : a->z() ? 0 : 1;

		a->release();

		return c;
	}

	nn_integer inot() const {
		nn_integer result = nn_new(length_);

		for( intptr_t i = length_ + 1; i >= 0; i-- )
			result->data_[i] = ~data_[i];

		return result;
	}

	nn_integer iand(const nn_integer v) const {
		nn_integer result;
		uintptr_t i;

		if( length_ > v->length_ ){
			result = nn_new(length_);

			for( i = 0; i < v->length_; i++ )
				result->data_[i] = data_[i] & v->data_[i];

			while( i < length_ + 2 ){
				result->data_[i] = data_[i] & v->data_[v->length_];
				i++;
			}
		}
		else if( length_ < v->length_ ){
			result = nn_new(v->length_);

			for( i = 0; i < length_; i++ )
				result->data_[i] = data_[i] & v->data_[i];

			while( i < v->length_ + 2 ){
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

	nn_integer ior(const nn_integer v) const {
		nn_integer result;
		uintptr_t i;

		if( length_ > v->length_ ){
			result = nn_new(length_);

			for( i = 0; i < v->length_; i++ )
				result->data_[i] = data_[i] | v->data_[i];

			while( i < length_ + 2 ){
				result->data_[i] = data_[i] | v->data_[v->length_];
				i++;
			}
		}
		else if( length_ < v->length_ ){
			result = nn_new(v->length_);

			for( i = 0; i < length_; i++ )
				result->data_[i] = data_[i] | v->data_[i];

			while( i < v->length_ + 2 ){
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

	nn_integer ixor(const nn_integer v) const {
		nn_integer result;
		uintptr_t i;

		if( length_ > v->length_ ){
			result = nn_new(length_);

			for( i = 0; i < v->length_; i++ )
				result->data_[i] = data_[i] ^ v->data_[i];

			while( i < length_ + 2 ){
				result->data_[i] = data_[i] ^ v->data_[v->length_];
				i++;
			}
		}
		else if( length_ < v->length_ ){
			result = nn_new(v->length_);

			for( i = 0; i < length_; i++ )
				result->data_[i] = data_[i] ^ v->data_[i];

			while( i < v->length_ + 2 ){
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

	uint8_t ibit(uintptr_t i) const {
		assert( i < length_ * sizeof(word) * CHAR_BIT );
		return data_[i / (sizeof(word) * CHAR_BIT)] >> (i & (sizeof(word) * CHAR_BIT - 1)) & 1;
	}

	template <typename T = word>
	void sbit(uintptr_t i,const T v = 1) const {
		assert( i < length_ * sizeof(word) * CHAR_BIT );
		data_[i / (sizeof(word) * CHAR_BIT)] |= word(v) << (i & (sizeof(word) * CHAR_BIT - 1));
	}

	uintptr_t count_trailing_zeros() const {
		uintptr_t n;
		word * p = data_;

		for( n = 0; n < length_ * sizeof(word) * CHAR_BIT; n += sizeof(word) * CHAR_BIT, p++ ){
			word x = *p;

			if( x != 0 ){
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

} nn_integer_data;
//------------------------------------------------------------------------------
typedef nn_integer_data::nn_integer nn_integer;
//------------------------------------------------------------------------------
nn_integer_data nn_izero = { 1, 1, 0, { 0 } };
nn_integer_data nn_ione = { 1, 1, 0, { 1 } };
nn_integer_data nn_itwo = { 1, 1, 0, { 2 } };
nn_integer_data nn_ifour = { 1, 1, 0, { 4 } };
nn_integer_data nn_ifive = { 1, 1, 0, { 5 } };
nn_integer_data nn_isix = { 1, 1, 0, { 6 } };
nn_integer_data nn_ieight = { 1, 1, 0, { 8 } };
nn_integer_data nn_iten = { 1, 1, 0, { 10 } };
// 1000000000u								== 0x3B9ACA00
// 10000000000000000000u					== 0x8AC7230489E80000
// 100000000000000000000000000000000000000u	== 0x4B3B4CA85A86C47A098A224000000000
#if SIZEOF_WORD == 1
nn_integer_data nn_maxull = { 1, 4, 0, { 0x00, 0xCA, 0x9A, 0x3B } };
#elif SIZEOF_WORD == 2
nn_integer_data nn_maxull = { 1, 4, 0, { 0x0000, 0x89E8, 0x2304, 0x8AC7 } };
#elif SIZEOF_WORD == 4
nn_integer_data nn_maxull = { 1, 2, 0, { 0x89E80000, 0x8AC72304 } };
#elif SIZEOF_WORD == 8
nn_integer_data nn_maxull = { 1, 2, 0, { 0x098A224000000000ull, 0x4B3B4CA85A86C47Aull } };
#endif
//------------------------------------------------------------------------------
static inline nn_integer nn_new(uintptr_t length)
{
	nn_integer p = nn_integer_data::nn_new(length);
	if( p == nullptr )
		throw std::bad_alloc();
	return p;
}
//------------------------------------------------------------------------------
static inline nn_integer nn_init_illong(long long v)
{
	switch( v ){
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
	p->data_[p->length_] = v >> (sizeof(v) * CHAR_BIT - 1);
	p->normalize();

	return p;
}
//------------------------------------------------------------------------------
static inline nn_integer nn_init_iullong(unsigned long long v)
{
	switch( v ){
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
static inline nn_integer nn_init_ulong(long v)
{
	return nn_init_illong(v);
}
//------------------------------------------------------------------------------
static inline nn_integer nn_init_iulong(unsigned long v)
{
	return nn_init_iullong(v);
}
//------------------------------------------------------------------------------
#if __GNUC__ && NDEBUG
#pragma GCC push_options
//#pragma GCC optimize("O2")
#endif
//------------------------------------------------------------------------------
inline nn_integer nn_integer_data::isal(uintptr_t bit_count) const
{
	if( bit_count == 0 )
		return add_ref();

	uintptr_t new_bit_size = length_ * sizeof(word) * CHAR_BIT + bit_count;
	new_bit_size += -(intptr_t) new_bit_size & (sizeof(word) * CHAR_BIT - 1);

	nn_integer result = nn_new(new_bit_size / (sizeof(word) * CHAR_BIT));

	uintptr_t offset = bit_count / (sizeof(word) * CHAR_BIT);
	uintptr_t shift = bit_count & (sizeof(word) * CHAR_BIT - 1);

	union {
		word * w;
		dword * d;
	} dst, src;

	src.w = data_ + length_ - 1;
	dst.w = result->data_ + result->length_ - 2;

	memset(result->data_, 0, offset * sizeof(word));

	if( shift == 0 ){
		//dst.w++;
		memcpy(result->data_ + offset,data_,length_ * sizeof(word));
	}
	else {
		for( intptr_t i = length_ - 1; i >= 0; i--, dst.w--, src.w-- )
		*dst.d = *src.d << shift;
	}

	result->data_[result->length_ + 1] = result->data_[result->length_] = isign();
	result->normalize();

	//#if BYTE_ORDER == BIG_ENDIAN
	//#error Not implemented, ENOSYS
	//#elif BYTE_ORDER == LITTLE_ENDIAN
	//#error Not implemented, ENOSYS
	//#endif // BYTE_ORDER

	return result;
}
//------------------------------------------------------------------------------
#if __GNUC__ && NDEBUG
#pragma GCC pop_options
#endif
//------------------------------------------------------------------------------
inline nn_integer nn_integer_data::isar(uintptr_t bit_count) const
{
	if( bit_count == 0 )
		return add_ref();

	uintptr_t bit_size = length_ * sizeof(word) * CHAR_BIT;

	if( bit_count >= bit_size )
		return nn_izero.add_ref();

	intptr_t new_bit_size = bit_size - bit_count;
	new_bit_size += -intptr_t(new_bit_size) & (sizeof(word) * CHAR_BIT - 1);

	nn_integer result = nn_new(new_bit_size / (sizeof(word) * CHAR_BIT));

	uintptr_t offset = bit_count / (sizeof(word) * CHAR_BIT);
	uintptr_t shift = bit_count & (sizeof(word) * CHAR_BIT - 1);

	union {
		word * w;
		dword * d;
	} dst, src;

	src.w = data_ + offset;
	dst.w = result->data_;

	//if( offset > 0 )
	//  src_word++;

	if( shift == 0 ){
		memcpy(result->data_,data_ + offset,result->length_ * sizeof(word));
	}
	else {
		for( intptr_t i = result->length_ - 1; i >= 0; i--, dst.w++, src.w++ )
			*dst.d = *src.d >> shift;
	}

	result->dummy_ = 0;
	result->data_[result->length_ + 1] = result->data_[result->length_] = isign();
	result->normalize();

	return result;
}
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
class integer {
	friend class numeric;
	public:
		static uintptr_t stat_iadd_;
		static uintptr_t stat_isub_;
		static uintptr_t stat_imul_;
		static uintptr_t stat_idiv_;

		~integer(){
			proxy_->release();
		}

		// at user level, must not be used
		integer(nn_integer p) : proxy_(p) {}
		integer(nn_integer p, int) : proxy_(p->add_ref()) {}

		// at user level, must not be used
		nn_integer reset(nn_integer p){
			nn_integer r = proxy_;
			proxy_ = p->add_ref();
			return r;
		}

		integer() : integer((long long) 0) {}
		integer(int v) : integer((long long) v) {}
		integer(unsigned int v) : integer((unsigned long long) (v)) {}
		integer(long v) : integer((long long) v) {}
		integer(unsigned long v) : integer((unsigned long long) v) {}
		integer(long long v) : proxy_(nn_init_illong(v)) {}
		integer(unsigned long long v) : proxy_(nn_init_iullong(v)) {}
#if _MSC_VER || HAVE_LONG_DOUBLE
		integer(double v);
		integer(long double v);
#elif defined(LONG_DOUBLE)
		integer(double v);
		integer(LONG_DOUBLE v);
#else
		integer(double v);
#endif
		integer(const std::string & s);
		integer(const integer & v);
		integer(const numeric & v);

		integer & operator = (int v){
			return *this = integer(v);
		}

		integer & operator = (unsigned int v){
			return *this = integer(v);
		}

		integer & operator = (long v){
			return *this = integer(v);
		}

		integer & operator = (unsigned long v){
			return *this = integer(v);
		}

		integer & operator = (long long v){
			return *this = integer(v);
		}

		integer & operator = (unsigned long long v){
			return *this = integer(v);
		}

		integer & operator = (const integer & v);
		integer & operator = (const numeric & v);

		integer operator + (const integer & v) const {
			stat_iadd_++;
			return proxy_->iadd(v.proxy_);
		}

		integer & operator += (const integer & v){
			return *this = *this + v;
		}

		integer operator - (const integer & v) const {
			stat_isub_++;
			return proxy_->isub(v.proxy_);
		}

		integer & operator -= (const integer & v){
			return *this = *this - v;
		}

		integer operator << (uintptr_t bit_count) const {
			return proxy_->isal(bit_count);
		}

		integer & operator <<= (uintptr_t bit_count){
			return *this = operator << (bit_count);
		}

		integer operator >> (uintptr_t bit_count) const {
			return proxy_->isar(bit_count);
		}

		integer & operator >>= (uintptr_t bit_count){
			return *this = operator >> (bit_count);
		}

		integer operator * (const integer & v) const;
		integer & operator *= (const integer & v){
			return *this = *this * v;
		}

		integer operator / (const integer & v) const;
		integer & operator /= (const integer & v){
			return *this = div(v);
		}

		integer operator % (const integer & v) const {
			integer mod;
			div(v, &mod);
			return mod;
		}

		integer & operator %= (const integer & v){
			integer mod;
			div(v, &mod);
			return *this = mod;
		}

		integer operator & (const integer & v) const {
			return proxy_->iand(v.proxy_);
		}

		integer & operator &= (const integer & v) {
			return *this = *this & v;
		}

		integer operator | (const integer & v) const {
			return proxy_->ior(v.proxy_);
		}

		integer & operator |= (const integer & v) {
			return *this = *this | v;
		}

		integer operator ^ (const integer & v) const {
			return proxy_->ixor(v.proxy_);
		}

		integer & operator ^= (const integer & v) {
			return *this = *this ^ v;
		}

		bool operator >  (const integer & v) const {
			return compare(v) > 0;
		}

		bool operator >= (const integer & v) const {
			return compare(v) >= 0;
		}

		bool operator <  (const integer & v) const {
			return compare(v) < 0;
		}

		bool operator <= (const integer & v) const {
			return compare(v) <= 0;
		}

		bool operator == (const integer & v) const {
			return compare(v) == 0;
		}

		bool operator != (const integer & v) const {
			return compare(v) != 0;
		}

		integer operator ++ () {
			return *this += 1;
		}

		// postfix form
		integer operator ++ (int) {
			return *this += 1;
		}

		integer operator -- () {
			return *this -= 1;
		}

		// postfix form
		integer operator -- (int) {
			return *this -= 1;
		}

		integer operator + () const {
			return *this;
		}

		integer operator - () const {
			return integer(0) - *this;
		}

		integer operator ~ () const {
			return proxy_->inot();
		}

		bool operator ! () const {
			return is_zero();
		}

		friend std::ostream & operator << (std::ostream & out, const nn::integer & v);

		void base8digits(std::string & s,const bool keep_leading_zeros = false) const;
		void base8string(std::string & s,const bool keep_leading_zeros = false) const;
		void base10string(std::string & s,const bool keep_leading_zeros = false) const;
		void base16string(std::string & s,const bool uppercase = false) const;

		const std::string to_string(uintptr_t width = 0,uintptr_t base = 10) const;

		operator const std::string () const {
			return to_string();
		}

#if _GLIBCXX_USE_WCHAR_T
		operator std::wstring () const;
#endif

#if ((__cplusplus >= 201103L) && defined(_GLIBCXX_USE_C99_STDINT_TR1))
		operator std::u16string () const;
		operator std::u32string () const;
#endif

		integer abs() const;

		bool is_neg() const {
			return sign() < 0;
		}

		bool is_zero() const;
		bool is_one() const;
		bool is_ten() const;

		uintptr_t bit(uintptr_t i) const {
			assert(i < proxy_->length_ * sizeof(word) * CHAR_BIT);

			return proxy_->data_[i / (sizeof(word) * CHAR_BIT)]
				>> (i & (sizeof(word) * CHAR_BIT - 1)) & 1;
		}

		sword sign() const {
			return proxy_->isign();
		}

		intptr_t compare(const integer & v) const {
			return proxy_->icompare(v.proxy_);
		}

		integer div(const integer & divider, integer * p_mod = NULL) const;

		integer nod_nok(const integer & a, integer * p_nok = NULL) const;

		integer & swap(integer & v) {
			nn_integer t = proxy_;
			proxy_ = v.proxy_;
			v.proxy_ = t;
			return *this;
		}

		// greatest common divisor (наибольший общий делитель - НОД)
		integer gcd(const integer & v) const {
			// GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0
			if( is_zero() )
				return v;
			if( v.is_zero() )
				return *this;

			return gcd_impl(this->abs(),v.abs());
		}

		// least common multiple (наименьшее общее кратное - НОК)
		integer lcm(const integer & v) const {
			return *this * v / gcd(v);
		}

		integer pow(uintptr_t power) const;

		static integer pow(uintptr_t power, uintptr_t base);
		static integer factorial(uintptr_t degree);

	protected:
		mutable nn_integer proxy_;

		template <typename T = word>
		void sbit(uintptr_t i,const T v = 1) const {
			if( proxy_->ref_count_ > 1 || i >= proxy_->length_ * sizeof(word) * CHAR_BIT ){
				uintptr_t new_bit_size = (i + 1) + (-intptr_t(i + 1) & (sizeof(word) * CHAR_BIT - 1));
				uintptr_t new_size = new_bit_size / (sizeof(word) * CHAR_BIT);

				nn_integer result = nn_new(imax(new_size,proxy_->length_));

				memcpy(result->data_,proxy_->data_,sizeof(word) * proxy_->length_);

				for( uintptr_t a = proxy_->length_; a < result->length_ + 2; a++ )
					result->data_[a] = proxy_->data_[proxy_->length_];

				proxy_->release();
				proxy_ = result;
			}
			proxy_->data_[i / (sizeof(word) * CHAR_BIT)] |= word(v) << (i & (sizeof(word) * CHAR_BIT - 1));
		}

		static integer gcd_impl(integer u, integer v) {
			uintptr_t u_ctz = u.proxy_->count_trailing_zeros();
			uintptr_t v_ctz = v.proxy_->count_trailing_zeros();
			uintptr_t shift = imin(u_ctz,v_ctz);

			u >>= u_ctz;

			for(;;) {
				v >>= v_ctz;
				if( u > v )
					u.swap(v);
				v -= u;
				if( v.is_zero() )
					break;
				v_ctz = v.proxy_->count_trailing_zeros();
			}

			return u << shift;
		}

	private:
		integer unique() const;

		static void accumulate_with_carry(nn_integer r,uintptr_t index,word value);

		//static std::vector<integer> one_cache_;
		//static std::vector<integer> ten_cache_;
		//static std::vector<integer> factorial_cache_;
};
//------------------------------------------------------------------------------
uintptr_t integer::stat_iadd_ = 0;
uintptr_t integer::stat_isub_ = 0;
uintptr_t integer::stat_imul_ = 0;
uintptr_t integer::stat_idiv_ = 0;
//------------------------------------------------------------------------------
#if _MSC_VER || HAVE_LONG_DOUBLE
inline integer::integer(double v) : integer((long double) v) {}
inline integer::integer(long double v) : integer(0)
#elif defined(LONG_DOUBLE)
inline integer::integer(double v) : integer((LONG_DOUBLE) v) {}
inline integer::integer(LONG_DOUBLE v) : integer(0)
#else
inline integer::integer(double v) : integer(0)
#endif
{
#if _MSC_VER || HAVE_MODFL || _GLIBCXX_HAVE_MODFL
  long double ipart, m;
  modfl(v,&ipart);
#else
  double ipart, m;
  modf(v,&ipart);
#endif
  integer power(1);
#if _MSC_VER || HAVE_FABSL || _GLIBCXX_HAVE_FABSL
  for( ipart = fabsl(ipart); ipart >= 1; ipart /= 10 ){
#else
  for( ipart = fabs(ipart); ipart >= 1; ipart /= 10 ){
#endif
#if _MSC_VER || HAVE_MODFL || _GLIBCXX_HAVE_MODFL
    m = fmodl(ipart,10);
#else
    m = fmod(ipart,10);
#endif
    *this += power * sword(m);
    power *= integer(10);
  }

  if( v < 0 )
    *this = -*this;
}
//------------------------------------------------------------------------------
inline integer::integer(const integer & v) : proxy_(v.proxy_->add_ref())
{
}
//------------------------------------------------------------------------------
inline integer::integer(const std::string & s) : integer(0)
{
	std::string::const_iterator i(s.cend());
	std::string::const_iterator j(s.cbegin());

	if( (s.length() > 1 && j[0] == '0' && (j[1] == 'x' || j[1] == 'X'))
		|| (s.length() > 0 && (j[0] == 'x' || j[0] == 'X')) ){

		j += j[0] == '0' ? 2 : 1;
		uintptr_t l = (i - j) * 4, p = l;
		l += -intptr_t(l) & (sizeof(word) * CHAR_BIT - 1);

		nn_integer result = nn_new(l / (sizeof(word) * CHAR_BIT));

		memset(result->data_, 0, result->length_ * sizeof(word));

		while( j < i ){

			word v;
			std::string::value_type a = *j;

			if( a >= '0' && a <= '9' )
				v = a - '0';
			else if( a >= 'a' && a <= 'f' )
				v = a - 'a';
			else if( a >= 'A' && a <= 'F' )
				v = a - 'A';
			else
				throw std::invalid_argument(s);

			p -= 4;
			result->data_[p / (sizeof(word) * CHAR_BIT)] |= word(v & 0xF) << (p & (sizeof(word) * CHAR_BIT - 1));

			j++;

		}

		result->data_[result->length_] = sword(result->data_[result->length_ - 1]) >> (sizeof(word) * CHAR_BIT - 1);
		result->data_[result->length_ + 1] = result->data_[result->length_];

		*this = integer(result);
	}
	else {
		integer m(1);

		while( i > j ){

			i--;

			if( *i == '-' ){
				*this = -*this;
				break;
			}

			if( *i == '+' )
				break;

			if( isdigit(*i) ){
				*this += m * (*i - '0');
				m *= integer(10);
			}
			else {
				throw std::invalid_argument(s);
			}

		}
	}
}
//------------------------------------------------------------------------------
inline integer & integer::operator = (const integer & v)
{
	v.proxy_->add_ref();
	proxy_->release();
	proxy_ = v.proxy_;
	return *this;
}
//------------------------------------------------------------------------------
inline void integer::base8digits(std::string & s,const bool keep_leading_zeros) const
{
	uintptr_t bits = sizeof(word) * CHAR_BIT * proxy_->length_;
	if( bits % 3 )
		bits += 3 - bits % 3;
	s.reserve(bits / 3);
	uint8_t * p = (uint8_t *) proxy_->data_;

	for( intptr_t trit = bits - 3; trit >= 0; trit -= 3 ){
		uint16_t c = (*(uint16_t *) (p + (trit >> 3))) >> (trit & 3);
		s.push_back((char) c & 0x3);
	}

	if( !keep_leading_zeros ){
		// strip leading zeros
		for( uintptr_t index = 0; index < s.length(); index++ )
			if( s[index] != '\0' && index > 0 ){
				if( index >= s.length() )
					index = s.length() - 1;
				s.erase(0,index);
				break;
			}
	}
}
//------------------------------------------------------------------------------
inline void integer::base8string(std::string & s,const bool keep_leading_zeros) const
{
	base8digits(s,keep_leading_zeros);

	for( intptr_t i = s.length() - 1; i >= 0; i-- )
		s[i] += '0';
}
//------------------------------------------------------------------------------
inline void integer::base10string(std::string & s,const bool keep_leading_zeros) const
{
    base8digits(s,keep_leading_zeros);

    //------------------------------------------------------------------
    //  Convert the base 8 string to a base 10 string using the
    //  algorithm from "Semi-Numerical Methods", by Knuth.
    //  Double the K leading octal digits using decimal arithmetic
    //  and subtract them from the K + 1 leading digits using
    //  decimal arithmetic.
    //------------------------------------------------------------------

    char * digit_array_ptr = const_cast<char *>(s.data());
    uintptr_t digit_length = s.length();

    char * subtrahend_ptr = digit_length < 4096 ?
		(char *) alloca(digit_length) : new char [digit_length];

    for( uintptr_t k = 0; k < digit_length - 1; k++ ) {
        //--------------------------------------------------------------
        //  Double the K leading octal digits using base 10 arithmetic
        //  and copy these digits into the K + 1'st location in the
        //  subtrahend array.
        //--------------------------------------------------------------

        intptr_t j = 0;

        for( j = k + 1; j >= 0; j-- )
            subtrahend_ptr[j] = 0;

        for( j = k; j >= 0; j-- ) {
            char doubled_digit = digit_array_ptr[j] << 1;

            if( doubled_digit > 9 ) {
                subtrahend_ptr[j + 1] += doubled_digit - 10;
                subtrahend_ptr[j] += 1;
            }
            else {
                subtrahend_ptr[j + 1] += doubled_digit;
            }
        }

        //--------------------------------------------------------------
        //  Subtract the doubled digits from the original number
        //  using decimal arithmetic.
        //--------------------------------------------------------------

        for( intptr_t m = k + 1; m >= 0; m-- ) {
            char difference = digit_array_ptr[m] - subtrahend_ptr[m];

            if( difference < 0 ) {
                digit_array_ptr[m] =  difference + 10;

                if( m - 1 >= 0 )
                    digit_array_ptr[m - 1] -= 1;
            }
            else {
                digit_array_ptr[m] = difference;
            }
        }
    }

	if( digit_length >= 4096 )
		delete [] subtrahend_ptr;

    //------------------------------------------------------------------
    //  Convert the digits to characters. First skip all leading zeros.
    //------------------------------------------------------------------

    uintptr_t fnzp = 0;

	if( !keep_leading_zeros )
		for( fnzp = 0; fnzp < digit_length; fnzp++ )
			if( s[fnzp] != 0 )
	            break;

	if( fnzp > 0 )
		s.erase(0,fnzp - 1);

    for( intptr_t j = s.length() - 1; j >= 0; j-- )
        s[j] += '0';
}
//------------------------------------------------------------------------------
inline void integer::base16string(std::string & s,bool uppecase) const
{
	s.reserve(proxy_->length_ * sizeof(word) * 2);
	const char * syms = uppecase ? "0123456789ABCDEF" : "0123456789abcdef";

	for( intptr_t i = proxy_->length_; i >= 0; i-- ){
		word c = proxy_->data_[i];
#if SIZEOF_WORD >= 8
		s.push_back(syms[(c >> 60) & 0xF]);
		s.push_back(syms[(c >> 56) & 0xF]);
		s.push_back(syms[(c >> 52) & 0xF]);
		s.push_back(syms[(c >> 48) & 0xF]);
		s.push_back(syms[(c >> 44) & 0xF]);
		s.push_back(syms[(c >> 40) & 0xF]);
		s.push_back(syms[(c >> 36) & 0xF]);
		s.push_back(syms[(c >> 32) & 0xF]);
#endif
#if SIZEOF_WORD >= 4
		s.push_back(syms[(c >> 28) & 0xF]);
		s.push_back(syms[(c >> 24) & 0xF]);
		s.push_back(syms[(c >> 20) & 0xF]);
		s.push_back(syms[(c >> 16) & 0xF]);
#endif
#if SIZEOF_WORD >= 2
		s.push_back(syms[(c >> 12) & 0xF]);
		s.push_back(syms[(c >>  8) & 0xF]);
#endif
		s.push_back(syms[(c >>  4) & 0xF]);
		s.push_back(syms[(c >>  0) & 0xF]);
	}

	// strip leading zeros
	uintptr_t index;
	for( index = 0; index < s.length() && s[index] == '0'; index++ );
	if( index > 0 ){
		if( index >= s.length() )
			index = s.length() - 1;
		s.erase(0,index);
	}
}
//------------------------------------------------------------------------------
inline const std::string integer::to_string(uintptr_t width,uintptr_t base) const
{
	if( is_zero() )
		return "0";

	std::string t;

	if( base == 8 ){
		abs().base8string(t);
	}
	else if( base == 16 ){
		abs().base16string(t,true);
	}
	else {
#if 1
		t.reserve(imax(((proxy_->length_ * sizeof(word)) * 2375) / 1000,width));

		integer zero(0), m(abs()), q;
		integer pt(&nn_maxull, 0);

		do {
			m = m.div(pt, &q);
			void * p = q.proxy_->data_;
#if SIZEOF_WORD < 2
			// 9 digits
			unsigned long v = *(unsigned long *) p;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
#elif SIZEOF_WORD < 8
			// 19 digits
			unsigned long long v = *(unsigned long long *) p;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
#else
			// 39 digits
			// volatile because gcc optimizer generate runtime segmentation fault
			volatile __uint128_t v = *(__uint128_t *) p;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
			t.push_back((v % 10) + '0'); v /= 10;
#endif

		} while( !m.is_zero() );

		while( t.size() < width )
			t.push_back('0');

		std::reverse(t.begin(),t.end());
		//t.erase();
#else
		abs().base10string(t,true);
#endif
		intptr_t w = t.size() - width;

		if( w > 0 ){
			std::string::iterator i(t.begin()), e(t.end());

			while( i < e && w > 0 && *i == '0' ){
				w--;
				i++;
			}

			t.erase(t.begin(), i);
		}
	}

	if( is_neg() )
		t.insert(0, "-");

	return t;
}
//------------------------------------------------------------------------------
inline integer integer::abs() const
{
	return sign() < 0 ? -*this : *this;
}
//------------------------------------------------------------------------------
inline integer integer::unique() const
{
	nn_integer result = nn_new(proxy_->length_);

	memcpy(result->data_, proxy_->data_, sizeof(word) * (result->length_ + 2));

	return result;
}
//------------------------------------------------------------------------------
inline bool integer::is_zero() const
{
	if( proxy_ == &nn_izero )
		return true;

	if( proxy_->z() ){
		nn_izero.add_ref();
		proxy_->release();
		proxy_ = &nn_izero;
		return true;
	}

	return false;
}
//------------------------------------------------------------------------------
inline bool integer::is_one() const
{
	if( proxy_ == &nn_ione )
		return true;

	if( proxy_->data_[0] != 1 )
		return false;

	for( intptr_t i = proxy_->length_; i > 0; i-- )
		if( proxy_->data_[i] != 0 )
			return false;

	nn_ione.add_ref();
	proxy_->release();
	proxy_ = &nn_ione;

	return true;
}
//------------------------------------------------------------------------------
inline bool integer::is_ten() const
{
	if( proxy_ == &nn_iten )
		return true;

	if( proxy_->data_[0] != 10 )
		return false;

	for( intptr_t i = proxy_->length_; i > 0; i-- )
		if( proxy_->data_[i] != 0 )
			return false;

	nn_iten.add_ref();
	proxy_->release();
	proxy_ = &nn_iten;

	return true;
}
//------------------------------------------------------------------------------
inline integer integer::nod_nok(const integer & a, integer * p_nok) const
{
	if( a.is_zero() )
		::div(1, 0);
	//throw range_error("Integer divide by zero");

	integer x(abs()), y(a.abs()), u, v, q, g;

	if( is_zero() ){
		x = 0;
		if( p_nok != NULL ) *p_nok = 0;
	}
	else {
		intptr_t c;

		u = y;
		v = x;

		for( ;; ){
			c = x.compare(y);

			if( c > 0 ){
				if( p_nok != NULL ) g = u;
				for( q = y; x - (q << 1) >= y; q <<= 1 ) if( p_nok != NULL ) g <<= 1;
				x -= q;
				if( p_nok != NULL ) v += g;
			}
			else if( c < 0 ){
				if( p_nok != NULL ) g = v;
				for( q = x; y - (q << 1) >= x; q <<= 1 ) if( p_nok != NULL ) g <<= 1;
				y -= q;
				if( p_nok != NULL ) u += g;
			}
			else
				break;
		}

		if( p_nok != NULL ) *p_nok = (u + v) >> 1;
	}
	return x;
}
//------------------------------------------------------------------------------
inline integer integer::pow(uintptr_t power) const
{
	if( is_zero() )
		return 0;
	if( is_one() || power == 0 )
		return 1;

	if( power == 1 )
		return *this;

	integer t(1), x(*this);

	if( is_ten() ){
		for( uintptr_t i = 1; uintptr_t(i) <= power; i++ )
			t = (t << 3) + (t << 1);
	}
	else {
		if( power == 2 )
			return x * x;

		if( power == 3 )
			return x * x * x;

		if( power == 4 ){
			x = x * x;
			return x * x;
		}

		while( power != 0 ){

			if( power % 2 != 0 ){
				t *= x;
				power -= 1;
			}

			x *= x;
			power /= 2;
		}
	}

	return t;
}
//------------------------------------------------------------------------------
inline integer integer::pow(uintptr_t power, uintptr_t base)
{
	return integer(base).pow(power);
}
//------------------------------------------------------------------------------
inline integer integer::factorial(uintptr_t degree)
{
#if 0
	if( factorial_cache_.empty() ){
		factorial_cache_.push_back(1);
		factorial_cache_.push_back(1);
	}

	uintptr_t i = factorial_cache_.size() - 1;
	integer t(factorial_cache_[i]);

	for( ++i; i <= degree; i++ ){
		t *= i;
		factorial_cache_.push_back(t);
	}

	return factorial_cache_[degree];
#else
	if( degree == 0 )
		return 1;
	if( degree == 1 )
		return 1;

	integer t(1);

	for( uintptr_t i = 2; i <= degree; i++ )
		t *= i;

	return t;
#endif
}
//------------------------------------------------------------------------------
inline std::ostream & operator << (std::ostream & out, const nn::integer & v)
{
	if( out.flags() & std::ios_base::hex ){
		if( out.flags() & std::ios_base::showbase )
			out << "0x";
		std::string s;
		v.base16string(s,(out.flags() & std::ios_base::uppercase) != 0);
		out << s;
	}
	else if( out.flags() & std::ios_base::oct ){
		if( out.flags() & std::ios_base::showbase )
			out << '0';
		std::string s;
		v.base8string(s);
		out << s;
	}
	else if( out.flags() & std::ios_base::dec ){
		if( (out.flags() & std::ios_base::showpos) || v.is_neg() ){
			out << (v.is_neg() ? '-' : '+');
		}
		else if( v.is_neg() ){
			out << '-';
		}
		out << v.abs().to_string(out.width(),10);
	}
	return out;
}
//------------------------------------------------------------------------------
inline void integer::accumulate_with_carry(nn_integer r,uintptr_t index,word value)
{
	union {
		struct {
			word lo;
			word hi;
		};
		dword temp;
	} d;

    while( value != 0 && index < r->length_ ) {
		d.temp = r->data_[index];
        d.temp += value;
        r->data_[index] = d.lo;
        value = d.hi;
        index++;
    }
}
//------------------------------------------------------------------------------
inline integer integer::operator * (const integer & v) const
{
	stat_imul_++;
#if 0
	if( is_zero() || v.is_zero() )
		return 0;
	if( is_one() )
		return v;
	if( v.is_one() )
		return *this;

	if( proxy_->length_ < v.proxy_->length_ )
		return v * *this;

	// hard way, but simple
	integer sum(0);
	integer a(abs());
	integer n(v.abs());
	uintptr_t bit_size = n.proxy_->length_ * sizeof(word) * CHAR_BIT;

	if( (is_neg() ^ v.is_neg()) != 0 ){
		for( uintptr_t i = 0; i < bit_size; i++ )
			if( n.bit(i) != 0 )
				sum -= a << i;
	}
	else {
		for( uintptr_t i = 0; i < bit_size; i++ )
			if( n.bit(i) != 0 )
				sum += a << i;
	}

	return sum;
#endif
	integer b(v.abs());

	if( b.proxy_->length_ <= 2 ){
		dword * p = (dword *) b.proxy_->data_;

		switch( *p ){
			case  0 : return integer(0);
			case  1 : return *this;
			case  2 : return *this << 1;
			case  3 : return (*this << 1) + *this;
			case  4 : return *this << 2;
			case  5 : return (*this << 2) + *this;
			case  6 : return (*this << 2) + (*this << 1);
			case  7 : return (*this << 2) + (*this << 1) + *this;
			case  8 : return *this << 3;
			case  9 : return (*this << 3) + *this;
			case 10 : return (*this << 3) + (*this << 1);
			case 11 : return (*this << 3) + (*this << 1) + *this;
			case 12 : return (*this << 3) + (*this << 2);
			case 13 : return (*this << 3) + (*this << 2) + *this;
			case 14 : return (*this << 3) + (*this << 2) + (*this << 1);
			case 15 : return (*this << 3) + (*this << 2) + (*this << 1) + *this;
			case 16 : return *this << 4;
		}
	}

	integer a(abs());
	nn_integer r = nn_new(a.proxy_->length_ + b.proxy_->length_ + 1);

	memset(r->data_,0,(r->length_ + 2) * sizeof(word));

#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
	register __m128i aa;
	__m128i bb, cc;
	uint32_t * bbp = (uint32_t *) &bb, * ccp = (uint32_t *) &cc;
#endif

	word * ap = a.proxy_->data_, * bp = b.proxy_->data_;

	for( uintptr_t jj = b.proxy_->length_, ii = a.proxy_->length_, i = 0; i < ii; i++ ){
#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
		ccp[2] = ccp[0] = ap[i];
		aa = _mm_load_si128(&cc);
#else
		dword h = ap[i];
#endif

		for( uintptr_t j = 0; j < jj; ){
#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
			bbp[0] = bp[j];
			bbp[2] = bp[j + 1];
			cc = _mm_mul_epu32(aa,bb);

			uintptr_t index = i + j;
			accumulate_with_carry(r, index + 0, ccp[0]);
			accumulate_with_carry(r, index + 1, ccp[1]);
			accumulate_with_carry(r, index + 1, ccp[2]);
			accumulate_with_carry(r, index + 2, ccp[3]);
			j += 2;
#else
			union {
				struct {
					word p0;
					word p1;
				};
				dword q;
			} d;
			d.q = h * bp[j];
			uintptr_t index = i + j;
			// accumulate the the low bit sum
			accumulate_with_carry(r, index, d.p0);
			// accumulate the the high bit sum
			accumulate_with_carry(r, index + 1, d.p1);
			j++;
#endif
		}
	}

	r->normalize();

	if( (is_neg() ^ v.is_neg()) != 0 )
		return -integer(r);

	return r;
}
//------------------------------------------------------------------------------
inline integer integer::operator / (const integer & v) const
{
	stat_idiv_++;
	return div(v);
}
//------------------------------------------------------------------------------
inline integer integer::div(const integer & divider, integer * p_mod) const
{
	if( divider.is_zero() ) //throw range_error("Integer divide by zero");
		::div(1, 0);

	if( divider.is_one() )
		return *this;

	integer temp_r, & r = p_mod == NULL ? temp_r : *p_mod;

	r = 0;

	if( is_zero() )
		return 0;

	// http://en.wikipedia.org/wiki/Division_algorithm
	// divide N by D, placing the quotient in Q and the remainder in R
	// Q := 0                 -- initialize quotient and remainder to zero
	// R := 0
	// for i = n-1...0 do     -- where n is number of bits in N
	//	 R := R << 1          -- left-shift R by 1 bit
	//	 R(0) := N(i)         -- set the least-significant bit of R equal to bit i of the numerator
	//	 if R >= D then
	//		R := R - D
	//		Q(i) := 1
	//	 end
	// end
	integer n(abs()), d(divider.abs()), q(0);

	for( intptr_t i = n.proxy_->length_ * sizeof(word) * CHAR_BIT - 1; i >= 0; i-- ){
		r <<= 1;
		r.sbit(0,n.bit(i));

		if( r >= d ){
			r -= d;
			q.sbit(i);
		}
	}

	if( (is_neg() ^ divider.is_neg()) != 0 )
		q = -q;

	return q;
}
//------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
class numeric {
	friend class integer;
	public:
		~numeric();
		numeric();
		numeric(int v) : numeric(integer(v)) {}
		numeric(unsigned int v) : numeric(integer(v)) {}
		numeric(long v) : numeric(integer(v)) {}
		numeric(unsigned long v) : numeric(integer(v)) {}
		numeric(long long v) : numeric(integer(v)) {}
		numeric(unsigned long long v) : numeric(integer(v)) {}
		numeric(const integer & v);
		numeric(const integer & numerator, const integer & denominator);
		numeric(const numeric & v);
#if _MSC_VER || HAVE_LONG_DOUBLE
		numeric(double v);
		numeric(long double v);
#elif defined(LONG_DOUBLE)
		numeric(double v);
		numeric(LONG_DOUBLE v);
#else
		numeric(double v);
#endif
		numeric(const std::string & s);
		numeric(const std::string & numerator, const std::string & denominator) :
			numeric(integer(numerator), integer(denominator)) {}

		numeric & operator = (int v){
			return *this = integer(v);
		}

		numeric & operator = (unsigned int v){
			return *this = integer(v);
		}

		numeric & operator = (long v){
			return *this = integer(v);
		}

		numeric & operator = (unsigned long v){
			return *this = integer(v);
		}

		numeric & operator = (long long v){
			return *this = integer(v);
		}

		numeric & operator = (unsigned long long v){
			return *this = integer(v);
		}

		numeric & operator = (const integer & v);
		numeric & operator = (const numeric & v);

		numeric operator + (const integer & v) const;
		numeric & operator += (const integer & v){
			return *this = *this + v;
		}

		numeric operator + (const numeric & v) const;
		numeric & operator += (const numeric & v){
			return *this = *this + v;
		}

		numeric operator - (const integer & v) const;
		numeric & operator -= (const integer & v){
			return *this = *this - v;
		}

		numeric operator - (const numeric & v) const;
		numeric & operator -= (const numeric & v){
			return *this = *this - v;
		}

		numeric operator * (const integer & v) const;
		numeric & operator *= (const integer & v){
			return *this = *this * v;
		}

		numeric operator * (const numeric & v) const;
		numeric & operator *= (const numeric & v){
			return *this = *this * v;
		}

		numeric operator / (const integer & v) const;
		numeric & operator /= (const integer & v){
			return *this = *this / v;
		}

		numeric operator / (const numeric & v) const;
		numeric & operator /= (const numeric & v){
			return *this = *this / v;
		}

		numeric operator ++ () {
			return *this += integer(1);
		}

		// postfix form
		numeric operator ++ (int) {
			return *this += integer(1);
		}

		numeric operator -- () {
			return *this -= integer(1);
		}

		// postfix form
		numeric operator -- (int) {
			return *this -= integer(1);
		}

		numeric operator + () const {
			return *this;
		}

		bool operator ! () const {
			return numerator_.is_zero();
		}

		numeric operator - () const;

		bool operator >  (const numeric & v) const;
		bool operator >= (const numeric & v) const;
		bool operator <  (const numeric & v) const;
		bool operator <= (const numeric & v) const;
		bool operator == (const numeric & v) const;
		bool operator != (const numeric & v) const;

		bool is_neg() const {
			return numerator_.is_neg();
		}

		bool is_zero() const {
			return numerator_.is_zero();
		}

		bool is_one() const {
			return numerator_.is_one() && denominator_.is_one();
		}

		numeric abs() const {
			return numeric(numerator_.abs(), denominator_);
		}

		numeric mod(integer * p_ipart = NULL) const;
		numeric pow(intptr_t power) const;

		numeric sin(uintptr_t iter = 2) const;
		numeric cos(uintptr_t iter = 2) const;
		numeric atan(uintptr_t iter = 2) const;
		numeric root(uintptr_t power = 2, uintptr_t iter = 16) const;

		enum round_type { toLess, toGreater };

		numeric round(uintptr_t digits, round_type type) const;

		friend std::ostream & operator << (std::ostream & out, const nn::numeric & v);

		const std::string to_string(uintptr_t width = 0, uintptr_t precision = 0) const;

		operator const std::string () const {
			return to_string();
		}

#if _GLIBCXX_USE_WCHAR_T
		operator std::wstring () const;
#endif

#if ((__cplusplus >= 201103L) && defined(_GLIBCXX_USE_C99_STDINT_TR1))
		operator std::u16string () const;
		operator std::u32string () const;
#endif

		//numeric & normalize(uintptr_t threshold = 32, int how = 1);
		numeric & normalize(int how = 1);

		const std::string stat_length() const {
			std::stringstream s;
			s << numerator_.proxy_->length_ << ", " << denominator_.proxy_->length_;
			return s.str();
		}

	protected:
	private:
		integer numerator_;   // числитель
		integer denominator_; // знаменатель
};
//------------------------------------------------------------------------------
inline integer::integer(const numeric & v) : integer(v.numerator_.div(v.denominator_))
{
}
//------------------------------------------------------------------------------
inline integer & integer::operator = (const numeric & v)
{
	return *this = v.numerator_.div(v.denominator_);
}
//------------------------------------------------------------------------------
inline numeric::~numeric()
{
}
//------------------------------------------------------------------------------
inline numeric::numeric() : numerator_(integer(0)), denominator_(integer(1))
{
}
//------------------------------------------------------------------------------
inline numeric::numeric(const integer & v) : numerator_(v), denominator_(integer(1))
{
}
//------------------------------------------------------------------------------
inline numeric::numeric(const numeric & v) : numerator_(v.numerator_), denominator_(v.denominator_)
{
	normalize();
}
//------------------------------------------------------------------------------
inline numeric::numeric(const std::string & s) : numeric(0)
{
	std::string::const_iterator i(s.cend()), fp(i);
	integer m(1), ipart, fraction;
	integer * p_i = s.rfind('.') == std::string::npos ? &ipart : &fraction;

	while( i > s.cbegin() ){

		i--;

		if( *i == '-' ){
			ipart = -ipart;
			break;
		}

		if( *i == '+' )
			break;

		if( *i == '.' ){
			fp = i;
			m = 1;
			p_i = &ipart;
		}
		else if( isdigit(*i) ){
			*p_i += m * (*i - '0');
			m *= 10;
		}
		else {
			throw std::invalid_argument(s);
		}
	}

	numerator_ = ipart;

	if( s.end() - fp > 0 ){
		denominator_ = integer(10).pow(s.end() - fp - 1);
		numerator_ *= denominator_;
		numerator_ += fraction;
	}
}
//------------------------------------------------------------------------------
#if _MSC_VER || HAVE_LONG_DOUBLE
inline numeric::numeric(double v) : numeric((long double) v) {}
inline numeric::numeric(long double v) : numeric()
#elif defined(LONG_DOUBLE)
inline numeric::numeric(double v) : numeric((LONG_DOUBLE) v) {}
inline numeric::numeric(LONG_DOUBLE v) : numeric()
#else
inline numeric::numeric(double v) : numeric()
#endif
{
	uintptr_t power;

#if _MSC_VER || HAVE_MODFL || _GLIBCXX_HAVE_MODFL
	long double ipart, fraction, fipart, a;

	fraction = modfl(v, &ipart);

	for( power = 0;; fraction *= 10, power++ ){
		a = modfl(fraction, &fipart);
		if( a == 0 ) break;
	}
#else
	double ipart, fraction, fipart;

	fraction = modf(v, &ipart);

	for( power = 0; modf(fraction, &fipart) != 0; fraction *= 10, power++ );
#endif

#if __GNUC__
	numerator_ = integer((double) ipart);
#else
	numerator_ = integer(ipart);
#endif

	if( power > 0 ){
		denominator_ = integer(10).pow(power);
		numerator_ *= denominator_;
#if __GNUC__
		numerator_ += integer((double) fraction);
#else
		numerator_ += integer(fraction);
#endif
	}
}
//------------------------------------------------------------------------------
inline numeric::numeric(const integer & numerator, const integer & denominator) :
	numerator_(numerator), denominator_(denominator)
{
	// WARNING: don't call here normalize(), here because the comparison stops working properly
}
//------------------------------------------------------------------------------
inline numeric & numeric::operator = (const integer & v)
{
	numerator_ = v;
	denominator_ = 1;
	return *this;
}
//------------------------------------------------------------------------------
inline numeric & numeric::operator = (const numeric & v)
{
	numerator_ = v.numerator_;
	denominator_ = v.denominator_;
	return normalize();
}
//------------------------------------------------------------------------------
inline numeric numeric::operator + (const integer & v) const
{
	return numeric(numerator_ + v * denominator_, denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator + (const numeric & v) const
{
	return numeric(numerator_ * v.denominator_ + v.numerator_ * denominator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator - (const integer & v) const
{
	return numeric(numerator_ - v * denominator_, denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator - (const numeric & v) const
{
	return numeric(numerator_ * v.denominator_ - v.numerator_ * denominator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator * (const integer & v) const
{
	return numeric(numerator_ * v, denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator * (const numeric & v) const
{
	return numeric(numerator_ * v.numerator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
inline numeric numeric::operator / (const integer & v) const
{
	if( v.is_zero() )
		::div(1, 0);
	//throw range_error("Numeric divide by zero");

	numeric t(numerator_.abs(), denominator_ * v.abs());

	if( (numerator_.is_neg() ^ v.is_neg()) != 0 )
		t.numerator_ = -t.numerator_;

	return t;
}
//------------------------------------------------------------------------------
inline numeric numeric::operator / (const numeric & v) const
{
	if( v.numerator_.is_zero() )
		::div(1, 0);
	//throw range_error("Numeric divide by zero");

	numeric t(numerator_.abs() * v.denominator_, denominator_ * v.numerator_.abs());

	if( (numerator_.is_neg() ^ v.numerator_.is_neg()) != 0 )
		t.numerator_ = -t.numerator_;

	return t;
}
//------------------------------------------------------------------------------
inline numeric numeric::operator - () const
{
	return numeric(-numerator_, denominator_);
}
//------------------------------------------------------------------------------
inline bool numeric::operator >  (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ > n.numerator_;
}
//------------------------------------------------------------------------------
inline bool numeric::operator >= (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ >= n.numerator_;
}
//------------------------------------------------------------------------------
inline bool numeric::operator <  (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ < n.numerator_;
}
//------------------------------------------------------------------------------
inline bool numeric::operator <= (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ <= n.numerator_;
}
//------------------------------------------------------------------------------
inline bool numeric::operator == (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ == n.numerator_;
}
//------------------------------------------------------------------------------
inline bool numeric::operator != (const numeric & v) const
{
	// приводим к общему знаменателю
	numeric u(numerator_ * v.denominator_, denominator_ * v.denominator_);
	numeric n(v.numerator_ * denominator_, v.denominator_ * denominator_);

	return u.numerator_ != n.numerator_;
}
//------------------------------------------------------------------------------
#if 0
inline numeric & numeric::normalize(uintptr_t threshold, int how)
{
	uintptr_t nl = numerator_.proxy_->length_;
	uintptr_t dl = denominator_.proxy_->length_;

	if( (how & 1) != 0 && (threshold == 0 || nl >= threshold || dl >= threshold) ){
		// fast path
		uintptr_t i = 0;
		uintptr_t shift = 0;
		union {
			const word * pw;
			const uint8_t * pb;
		} p0, p1;
		p0.pw = numerator_.proxy_->data_;
		p1.pw = denominator_.proxy_->data_;

		while( i < nl && i < dl ){
			if( *p0.pw != 0 || *p1.pw != 0 ) break;
			shift += sizeof(word) * CHAR_BIT;
			p0.pw++;
			p1.pw++;
			i++;
		}

		i *= sizeof(word);
		nl *= sizeof(word);
		dl *= sizeof(word);

		while( i < nl && i < dl ){
			if( *p0.pb != 0 || *p1.pb != 0 ) break;
			shift += sizeof(uint8_t) * CHAR_BIT;
			p0.pb++;
			p1.pb++;
			i++;
		}

		i *= CHAR_BIT;
		nl *= CHAR_BIT;
		dl *= CHAR_BIT;

		while( i < nl && i < dl && numerator_.bit(i) == 0 && denominator_.bit(i) == 0 )
			i++;

		if( shift > 0 ){
			numerator_ >>= shift;
			denominator_ >>= shift;
		}

		nl = numerator_.proxy_->length_;
		dl = denominator_.proxy_->length_;
	}

	if( (how & 2) != 0 && (threshold == 0 || nl >= threshold || dl >= threshold) ){

		/*if( (nl >= dl && dl * 1.618 >= nl)
		  || nl < dl && nl * 1.618 >= dl )*/{

			// упрощаем дробь через наибольший общий делитель
			integer nod(numerator_.nod_nok(denominator_, NULL));

			if( !nod.is_zero() && !nod.is_one() ){
				numerator_ /= nod;
				denominator_ /= nod;
			}
		}
	}
#else
inline numeric & numeric::normalize(int how)
{
	if( how & 2 ){
		integer d(numerator_.gcd(denominator_));
		//integer nod(numerator_.nod_nok(denominator_, NULL));

		if( !d.is_zero() ){
			numerator_ /= d;
			denominator_ /= d;
		}
	}

	if( how & 1 ){
		uintptr_t i = 0;
		uintptr_t shift = 0;
		union {
			const word * pw;
			const uint8_t * pb;
		} p0, p1;
		p0.pw = numerator_.proxy_->data_;
		p1.pw = denominator_.proxy_->data_;

		uintptr_t nl = numerator_.proxy_->length_;
		uintptr_t dl = denominator_.proxy_->length_;

		while( i < nl && i < dl ){
			if( *p0.pw != 0 || *p1.pw != 0 ) break;
			shift += sizeof(word) * CHAR_BIT;
			p0.pw++;
			p1.pw++;
			i++;
		}

		i *= sizeof(word);
		nl *= sizeof(word);
		dl *= sizeof(word);

		while( i < nl && i < dl ){
			if( *p0.pb != 0 || *p1.pb != 0 ) break;
			shift += sizeof(uint8_t) * CHAR_BIT;
			p0.pb++;
			p1.pb++;
			i++;
		}

		i *= CHAR_BIT;
		nl *= CHAR_BIT;
		dl *= CHAR_BIT;

		while( i < nl && i < dl && numerator_.bit(i) == 0 && denominator_.bit(i) == 0 )
			i++;

		if( shift > 0 ){
			numerator_ >>= shift;
			denominator_ >>= shift;
		}
	}
#endif
	return *this;
}
//------------------------------------------------------------------------------
inline numeric numeric::mod(integer * p_ipart) const
{
	integer ipart(numerator_.div(denominator_));

	if( p_ipart != NULL )
		*p_ipart = ipart;

	return *this - ipart;
}
//------------------------------------------------------------------------------
inline numeric numeric::pow(intptr_t power) const
{
	if( power == 0 )
		return 1;

	if( power == 1 )
		return *this;

	numeric x(*this);

	if( power == 2 )
		return x * x;

	if( power == 3 )
		return x * x * x;

	if( power == 4 ){
		x = x * x;
		return x * x;
	}

	uintptr_t n = power >= 0 ? power : -power;
	numeric t(1);

	while( n != 0 ){

		if( n % 2 != 0 ){
			t *= x;
			n -= 1;
		}

		x *= x;
		n /= 2;
	}

	if( power < 0 )
		t = numeric(1) / t;

	return t;
}
//------------------------------------------------------------------------------
inline numeric numeric::round(uintptr_t digits, round_type type) const
{
	integer ipart, q, p(integer(10).pow(digits));
	numeric fraction(mod(&ipart)), fq;

	fraction *= p;
	fq = fraction.mod(&q) * integer(10);

	switch( type ){
		case toLess:
			if( q.is_neg() ){
				if( q < -5 ) fraction += integer(-10) - q; else fraction -= q;
			}
			else {
				if( q > 5 ) fraction += integer(10) - q; else fraction -= q;
			}
			break;
		case toGreater:
			if( q.is_neg() ){
				if( q <= -5 ) fraction += integer(-10) - q; else fraction -= q;
			}
			else {
				if( q >= 5 ) fraction += integer(10) - q; else fraction -= q;
			}
			break;
	}

	fraction /= p;

	return ipart + fraction;
}
//------------------------------------------------------------------------------
inline numeric numeric::sin(uintptr_t iter) const
{
	// https://ru.wikipedia.org/wiki/Ряд_Тейлора
	intptr_t s = 1;
	numeric result(0), powerv(*this);
	uintptr_t cur_pow = 2;
	integer factorial(1);
	uintptr_t cur_degree = 2;

	auto f = [&] (uintptr_t degree) -> integer {
		if( degree == 0 )
			return factorial;
		if( degree == 1 )
			return factorial;

		for( ; cur_degree <= degree; cur_degree++ )
			factorial *= cur_degree;

		return factorial;
	};

	auto p = [&] (uintptr_t power) -> numeric {
		if( power == 0 )
			return 1;

		if( power == 1 || is_zero() || is_one() )
			return powerv;

		for( ; cur_pow <= power; cur_pow++ )
			powerv *= *this;

		return powerv.normalize(3);
	};

	for( uintptr_t n = 0; n < iter; n++, s = -s ){
		if( s >= 0 )
			result += numeric(1) / f(2 * n + 1) * p(2 * n + 1);
		else
			result -= numeric(1) / f(2 * n + 1) * p(2 * n + 1);
		result.normalize(3);
	}

	return result;
}
//------------------------------------------------------------------------------
inline numeric numeric::cos(uintptr_t iter) const
{
	// https://ru.wikipedia.org/wiki/Ряд_Тейлора
	intptr_t s = 1;
	numeric result(0), powerv(*this);
	uintptr_t cur_pow = 2;
	integer factorial(1);
	uintptr_t cur_degree = 2;

	auto f = [&] (uintptr_t degree) -> integer {
		if( degree == 0 )
			return factorial;
		if( degree == 1 )
			return factorial;

		for( ; cur_degree <= degree; cur_degree++ )
			factorial *= cur_degree;

		return factorial;
	};

	auto p = [&] (uintptr_t power) -> numeric {
		if( power == 0 )
			return 1;

		if( power == 1 || is_zero() || is_one() )
			return powerv;

		for( ; cur_pow <= power; cur_pow++ )
			powerv *= *this;

		return powerv.normalize(3);
	};

	for( uintptr_t n = 0; n < iter; n++, s = -s ){
		if( s >= 0 )
			result += numeric(1) / f(2 * n) * p(2 * n);
		else
			result -= numeric(1) / f(2 * n) * p(2 * n);
		result.normalize(3);
	}

	return result;
}
//------------------------------------------------------------------------------
inline numeric numeric::atan(uintptr_t iter) const
{
	numeric result(0);
	numeric pow2np1(abs()), n2np1(1), pow2(*this * *this);

	for( uintptr_t n = 0; n < iter; n++ ){
		//    atan += one.Pow(n) * Pow(2 * n + 1) / (2 * n + 1);
		if( (n & 1) == 0 )
			result += pow2np1 / n2np1;
		else
			result -= pow2np1 / n2np1;
		pow2np1 *= pow2;
		n2np1 += integer(2);
	}

	return result;
}
//------------------------------------------------------------------------------
inline numeric numeric::root(uintptr_t power, uintptr_t iter) const
{
	if( *this == 1 )
		return *this;

	numeric r((power & 1) != 0 ? abs() : *this);
	r.normalize(3);

	if( r.is_zero() ){
	}
	else if( power == 2 ){
		// https://ru.wikipedia.org/wiki/Ряд_Тейлора
		numeric x(r), hone(1, 2);

		for( uintptr_t n = 0; n < iter; n++ ){
			x = hone * (x + (r / x));
			x.normalize(3);
		}

		r = x;
	}
	else { // методом деления отрезка пополам
		numeric a(0), b(r);
		numeric m, mm;

		for( uintptr_t n = 0; a < b && n < iter; n++ ){
			m = a + (b - a) / integer(2);
			mm = m.pow(power);
			//mm.normalize(3);

			if( mm < r ){
				a = m;
			}
			else if( mm > r ){
				b = m;
			}
			else
				break;
		}

		r = m;
	}

	if( (power & 1) != 0 && numerator_.sign() < 0 )
		r = -r;

	return r;
}
//------------------------------------------------------------------------------
inline const std::string numeric::to_string(uintptr_t width, uintptr_t precision) const
{
	integer ipart;
	numeric fraction(abs().mod(&ipart));
	std::string s(ipart.to_string(width)), f, ns, ds;

	if( numerator_.is_neg() )
		s.insert(0, "-");

	uintptr_t pre = precision == 0 ? 9 : precision;

	fraction *= integer::pow(pre, 10);
	fraction.mod(&ipart);
	f = ipart.to_string(pre);

	while( f.size() < pre )
		f.push_back('0');

	if( f.size() > pre )
		f = f.substr(0, pre);

	if( !f.empty() && precision == 0 ){
		std::string::iterator j(f.end()), i(j);

		while( i > f.begin() && *(i - 1) == '0' ) i--;

		f.erase(i, j);
	}

	return f.empty() ? s : s + "." + f;
}
//------------------------------------------------------------------------------
inline std::ostream & operator << (std::ostream & out, const nn::numeric & v)
{
	if( out.flags() & std::ios_base::hex ){
		out << (v.numerator_.is_neg() ? "-" : "")
			<< "x" << std::hex << v.numerator_.abs()
			<< "/" << "x" << std::hex << v.denominator_;
	}
	else if( out.flags() & std::ios_base::dec ){
		out << v.to_string(out.width(), out.precision());
	}

	return out;
}
//------------------------------------------------------------------------------
struct pi_cont {
	numeric e, four, five, six, one, two;
	integer sixteen, eight;
	uintptr_t i;

	pi_cont();

	pi_cont & to_iter(uintptr_t iter);
};
//------------------------------------------------------------------------------
// Continued fractions
static inline numeric cfpi(const integer & q, const integer & v, uintptr_t i, uintptr_t iter)
{
	if( i >= iter )
		return 1;
	return numeric(q) + numeric(v.pow(2)) / cfpi(q + 2, v + 1, i + 1, iter);
}
//------------------------------------------------------------------------------
inline pi_cont::pi_cont() :
	e(0), four(4), five(5), six(6), one(1), two(2), sixteen(16), eight(8), i(1)
{
	e = four - two / four - one / five - one / six;
}
//------------------------------------------------------------------------------
inline pi_cont & pi_cont::to_iter(uintptr_t iter)
{
	while( i < iter ){
		integer v8i(eight * i);
		numeric delta(
			(one / sixteen)
			* (four / (v8i + one)
			- two / (v8i + four)
			- one / (v8i + five)
			- one / (v8i + six))
			);

		e += delta;
		sixteen <<= 4;

		i++;
	}

	return *this;
}
//------------------------------------------------------------------------------
inline numeric pi(uintptr_t iter, pi_cont * p_cnt)
{
	/* формула Мэчина */
	//numeric n1_5(1,5), n1_239(1,239);
	//numeric result(n1_5.atan(iter) * integer(4) - n1_239.atan(iter));
	//return result * integer(4);

	// формула Рамануджана
	//numeric e(zero);

	//for( uintptr_t k = 0; k < iter; k++ ){
	//  numeric q(integer::factorial(4 * k) * (integer(1103) + integer(26390) * k));
	//  numeric m(integer::factorial(k).pow(4) * integer(396).pow(4 * k));

	//  e += q / m;
	//}

	//
	//return one / (numeric(two * two.sqrt(2,iter)) / numeric(9801));

	// формула Бэйли-Борвайна-Плаффа
	if( p_cnt == NULL ){
		pi_cont cnt;
		cnt.to_iter(iter);
		return cnt.e;
	}

	p_cnt->to_iter(iter);
	return p_cnt->e;

	// Continued fractions
	//return numeric(4) / cfpi(integer(1),integer(1),0,iter);
}
//------------------------------------------------------------------------------
} // namespace NaturalNumbers
//------------------------------------------------------------------------------
#if _MSC_VER
#undef SIZEOF_WORD
#endif
//------------------------------------------------------------------------------
#if _MSC_VER
#pragma warning(pop)
#endif
//------------------------------------------------------------------------------
#endif // NN_HPP_INCLUDED
//------------------------------------------------------------------------------

