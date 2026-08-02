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
#ifndef II_HPP_INCLUDED
#define II_HPP_INCLUDED
//------------------------------------------------------------------------------
#include "id.hpp"
#include "wk.hpp"
#include <cstdlib>
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @namespace nn (NaturalNumbers)
/// @brief Arbitrary-precision integer class
///
/// This file defines the user-facing `integer` type: an arbitrary-precision
/// signed integer built on a sign-magnitude representation backed by the
/// ref-counted `nn_integer_data` proxy pattern defined in id.hpp.
///
/// Key design decisions and references:
///   - Sign-magnitude representation using words of native width (typically
///     64-bit on x86_64, 32-bit on x86). The sign is stored in a dedicated
///     word at data_[length_] (0 = positive, ~0 = negative via two's complement
///     extension).
///   - Copy-on-write semantics: proxy_ is a shared pointer with atomic-like
///     reference counting. Mutation triggers `unique()` to create an
///     exclusive copy.
///   - Threaded multiplication: for sufficiently large operands, `operator*`
///     dispatches to a thread pool (wk.hpp) for parallel schoolbook
///     multiplication.
///
/// References:
///   - Knuth, D.E. "The Art of Computer Programming", Vol. 2:
///     Seminumerical Algorithms, 3rd ed., Addison-Wesley, 1998.
///   - Warren, H.S. "Hacker's Delight", 2nd ed., Addison-Wesley, 2013.
///   - Stein, J. "Computational problems associated with Racah algebra",
///     J. Comp. Phys., 1967 (binary GCD algorithm).
///   - Moller, N. and Granlund, T., "Improved division by invariant
///     integers", IEEE Trans. Computers, 2011.
// ---------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------
// //////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
/// @class integer
/// @brief User-facing arbitrary-precision integer type.
///
/// @details The `integer` class wraps `nn_integer_data` through a ref-counted
///   proxy pointer (`proxy_`). Construction is supported from native integer
///   types (`int`, `long`, `long long`, `unsigned*`), `double`/`long double`,
///   `std::string` (decimal and hex), and copy construction.
///
///   **Copy-on-write**: All read-only operations (comparison, bit tests)
///   operate on the shared proxy_ without copying. Any write operation first
///   calls `unique()` to duplicate the underlying data if `ref_count_ > 1`.
///
///   **Statistics**: Four static counters (`stat_iadd_`, `stat_isub_`,
///   `stat_imul_`, `stat_idiv_`) track the number of arithmetic operations
///   for profiling and benchmarking.
///
///   **Thread safety**: Individual integer objects are NOT thread-safe.
///   However, the global thread pool (`pool`) used by `tmul_*` is internally
///   synchronized.
// ---------------------------------------------------------------------------
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
		integer(const nn_integer_data * p, int) : proxy_(p->add_ref()) {}

		// at user level, must not be used
		nn_integer reset(nn_integer p){
			nn_integer r = proxy_;
			proxy_ = p->add_ref();
			return r;
		}

		integer() : integer(static_cast<long long>(0)) {}
		integer(int v) : integer(static_cast<long long>(v)) {}
		integer(unsigned int v) : integer(static_cast<unsigned long long>(v)) {}
		integer(long v) : integer(static_cast<long long>(v)) {}
		integer(unsigned long v) : integer(static_cast<unsigned long long>(v)) {}
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
			return *this = divide(v);
		}

		integer operator % (const integer & v) const {
			integer mod;
			divide(v, &mod);
			return mod;
		}

		integer & operator %= (const integer & v){
			integer mod;
			divide(v, &mod);
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

		integer & operator ++ () {
			return *this += 1;
		}

		// postfix form
		integer operator ++ (int) {
			return *this += 1;
		}

		integer & operator -- () {
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

		const std::string to_string(uintptr_t width = 0, uintptr_t base = 10) const;

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

		constexpr bool is_zero() const;
		constexpr bool is_one() const;
		constexpr bool is_ten() const;

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

		integer divide(const integer & divider, integer * p_mod = nullptr) const;

		integer nod_nok(const integer & a, integer * p_nok = nullptr) const;

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
		void sbit(uintptr_t i, const T v = 1) const {
			if( proxy_->ref_count_ > 1 || i >= proxy_->length_ * sizeof(word) * CHAR_BIT ){
				uintptr_t new_bit_size = (i + 1) + (-static_cast<intptr_t>(i + 1) & (sizeof(word) * CHAR_BIT - 1));
				uintptr_t new_size = new_bit_size / (sizeof(word) * CHAR_BIT);

				nn_integer result = nn_new(imax(new_size, proxy_->length_));

				memcpy(result->data_,proxy_->data_,sizeof(word) * proxy_->length_);

				for( uintptr_t a = proxy_->length_; a < result->length_ + 2; a++ )
					result->data_[a] = proxy_->data_[proxy_->length_];

				proxy_->release();
				proxy_ = result;
			}
			proxy_->data_[i / (sizeof(word) * CHAR_BIT)] |= static_cast<word>(v) << (i & (sizeof(word) * CHAR_BIT - 1));
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

		//static std::vector<integer> one_cache_;
		//static std::vector<integer> ten_cache_;
		//static std::vector<integer> factorial_cache_;

        template <size_t SIZE> integer tmul_a(const integer & v) const;
		integer tmul_v(const integer & v) const;
};
//------------------------------------------------------------------------------
uintptr_t integer::stat_iadd_ = 0;
uintptr_t integer::stat_isub_ = 0;
uintptr_t integer::stat_imul_ = 0;
uintptr_t integer::stat_idiv_ = 0;
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Constructs an integer from a floating-point value.
///
/// @details Decomposes the double into integer and fractional parts using
///   `modf`/`modfl`. Each decimal digit of the integer part is extracted
///   by repeated division by 10 and accumulated into `*this`. Then the
///   sign is adjusted.
///
///   Complexity: O(log10(v)) digit extractions, each requiring a decimal
///   multiplication and addition.
///
/// @param v The floating-point value to convert.
// ---------------------------------------------------------------------------
#if _MSC_VER || HAVE_LONG_DOUBLE
inline integer::integer(double v) : integer(static_cast<long double>(v)) {}
inline integer::integer(long double v) : integer(0)
#elif defined(LONG_DOUBLE)
inline integer::integer(double v) : integer(static_cast<LONG_DOUBLE>(v)) {}
inline integer::integer(LONG_DOUBLE v) : integer(0)
#else
inline integer::integer(double v) : integer(0)
#endif
{
#if _MSC_VER || HAVE_MODFL || _GLIBCXX_HAVE_MODFL
  long double ipart, m = modfl(v, &ipart);
#else
  double ipart, m = modf(v, &ipart);
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
    *this += power * static_cast<sword>(m);
    power *= integer(10);
  }

  if( v < 0 )
    *this = -*this;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Copy constructor — increments the reference count of the shared
///   proxy.
///
/// @details O(1) operation. No deep copy is performed until the first
///   write triggers copy-on-write via `unique()`.
///
/// @param v The integer to copy.
// ---------------------------------------------------------------------------
inline integer::integer(const integer & v) : proxy_(v.proxy_->add_ref())
{
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Constructs an integer from a decimal or hexadecimal string.
///
/// @details Auto-detects hex strings by the "0x" or "X" prefix. For hex:
///   each nibble (4-bit) is placed directly into the word array in O(n)
///   time. For decimal: digits are processed right-to-left using Knuth's
///   digit-by-digit accumulation in O(n^2) time since each digit requires
///   a decimal multiplication and addition across the growing word array.
///
///   Edge cases:
///   - Leading '+' or '-' sign (sign is applied after accumulation)
///   - Invalid characters throw std::invalid_argument
///   - Empty string is treated as zero
///
///   Reference: Knuth, TAOCP Vol. 2, Section 4.4 (positional notation
///   conversion).
///
/// @param s The input string.
/// @throws std::invalid_argument If s contains invalid characters.
// ---------------------------------------------------------------------------
inline integer::integer(const std::string & s) : integer(0)
{
	std::string::const_iterator i(s.cend());
	std::string::const_iterator j(s.cbegin());

	if( (s.length() > 1 && j[0] == '0' && (j[1] == 'x' || j[1] == 'X'))
		|| (s.length() > 0 && (j[0] == 'x' || j[0] == 'X')) ){

		j += j[0] == '0' ? 2 : 1;
		uintptr_t l = (i - j) * 4, p = l;
		l += -static_cast<intptr_t>(l) & (sizeof(word) * CHAR_BIT - 1);

		nn_integer result = nn_new(l / (sizeof(word) * CHAR_BIT));

		memset(result->data_, 0, result->length_ * sizeof(word));

		while( j < i ){

			word v;
			std::string::value_type a = *j;

			if( a >= '0' && a <= '9' )
				v = a - '0';
			else if( a >= 'a' && a <= 'f' )
				v = a - 'a' + 10;
			else if( a >= 'A' && a <= 'F' )
				v = a - 'A' + 10;
			else
				throw std::invalid_argument(s);

			p -= 4;
			result->data_[p / (sizeof(word) * CHAR_BIT)] |= static_cast<word>(v & 0xF) << (p & (sizeof(word) * CHAR_BIT - 1));

			j++;

		}

		result->data_[result->length_] = static_cast<sword>(result->data_[result->length_ - 1]) >> (sizeof(word) * CHAR_BIT - 1);
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

			if( std::isdigit(static_cast<unsigned char>(*i)) ){
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
// ---------------------------------------------------------------------------
/// @brief Copy-assignment operator.
///
/// @details Increments the source proxy's ref count, releases this proxy,
///   and copies the pointer. O(1).
///
/// @param v The source integer.
/// @return Reference to this.
// ---------------------------------------------------------------------------
inline integer & integer::operator = (const integer & v)
{
	v.proxy_->add_ref();
	proxy_->release();
	proxy_ = v.proxy_;
	return *this;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Extracts base-8 (octal) digits from the internal binary
///   representation.
///
/// @details Reads groups of 3 bits from the word array (LSB-first packing)
///   and pushes each as a 0–7 value into the string. The bit-width is
///   rounded up to a multiple of 3. Leading zeros are stripped unless
///   `keep_leading_zeros` is true.
///
///   The digit values are stored as raw 0–7 values (not ASCII), needing
///   later conversion via `base8string`.
///
/// @param s Output string to receive digits.
/// @param keep_leading_zeros If true, preserves all trits including zero.
// ---------------------------------------------------------------------------
inline void integer::base8digits(std::string & s,const bool keep_leading_zeros) const
{
	uintptr_t bits = sizeof(word) * CHAR_BIT * proxy_->length_;
	if( bits % 3 )
		bits += 3 - bits % 3;
	s.reserve(bits / 3);
	uint8_t * p = reinterpret_cast<uint8_t *>(proxy_->data_);

	for( intptr_t trit = bits - 3; trit >= 0; trit -= 3 ){
		uint16_t c = (*reinterpret_cast<uint16_t *>(p + (trit >> 3))) >> (trit & 3);
		s.push_back(static_cast<char>(c & 0x7));
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
// ---------------------------------------------------------------------------
/// @brief Converts base-8 digit values to ASCII characters and writes to
///   a string.
///
/// @details This is a thin wrapper around `base8digits` that adds '0' to
///   each digit value to produce ASCII characters.
///
/// @param s Output string (will be overwritten).
/// @param keep_leading_zeros If true, retains leading zero trits.
// ---------------------------------------------------------------------------
inline void integer::base8string(std::string & s,const bool keep_leading_zeros) const
{
	base8digits(s,keep_leading_zeros);

	for( intptr_t i = s.length() - 1; i >= 0; i-- )
		s[i] += '0';
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Converts the internal binary representation to a decimal string.
///
/// @details Uses Knuth's "semi-numerical" algorithm: first extracts octal
///   digits via `base8digits`, then converts the octal string to decimal
///   in-place. The core idea: for each octal position k, double the k
///   leading octal digits using decimal arithmetic, subtract them from the
///   k+1 leading digits using decimal arithmetic. This avoids full
///   arbitrary-precision division.
///
///   A stack buffer of 4096 bytes is used for the subtrahend for strings
///   up to 4096 octal digits; longer strings allocate on the heap.
///
///   Reference: Knuth, TAOCP Vol. 2, Section 4.4 (conversion between radices).
///
///   Complexity: O(n^2) in the number of octal digits, but practical for
///   the sizes encountered in this library.
///
/// @param s Output string (will be overwritten).
/// @param keep_leading_zeros If true, retains leading zeros in the output.
// ---------------------------------------------------------------------------
inline void integer::base10string(std::string & s,const bool keep_leading_zeros) const
{
    base8digits(s, keep_leading_zeros);

    //------------------------------------------------------------------
    //  Convert the base 8 string to a base 10 string using the
    //  algorithm from "Semi-Numerical Methods", by Knuth.
    //  Double the K leading octal digits using decimal arithmetic
    //  and subtract them from the K + 1 leading digits using
    //  decimal arithmetic.
    //------------------------------------------------------------------

    char * digit_array_ptr = const_cast<char *>(s.data());
    uintptr_t digit_length = s.length();

	char * subtrahend_ptr = nullptr;

	auto gen = [&] () -> void {
		for (uintptr_t k = 0; k < digit_length - 1; k++) {
			//--------------------------------------------------------------
			//  Double the K leading octal digits using base 10 arithmetic
			//  and copy these digits into the K + 1'st location in the
			//  subtrahend array.
			//--------------------------------------------------------------

			intptr_t j = 0;

			for (j = k + 1; j >= 0; j--)
				subtrahend_ptr[j] = 0;

			for (j = k; j >= 0; j--) {
				char doubled_digit = digit_array_ptr[j] << 1;

				if (doubled_digit > 9) {
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

			for (intptr_t m = k + 1; m >= 0; m--) {
				char difference = digit_array_ptr[m] - subtrahend_ptr[m];

				if (difference < 0) {
					digit_array_ptr[m] = difference + 10;

					if (m - 1 >= 0)
						digit_array_ptr[m - 1] -= 1;
				}
				else {
					digit_array_ptr[m] = difference;
				}
			}
		}
	};

	if( digit_length < 4096 ) {
		char subtrahend_buf[4096];
		subtrahend_ptr = subtrahend_buf;
		gen();
	}
	else {
		subtrahend_ptr = new char[digit_length];
		gen();
		delete [] subtrahend_ptr;
	}

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
// ---------------------------------------------------------------------------
/// @brief Converts the internal binary representation to a hexadecimal
///   string.
///
/// @details Iterates over each word of the data array (MSB to LSB) and
///   extracts nibbles (4-bit groups) from most significant to least. For
///   SIZEOF_WORD == 8, each 64-bit word produces 16 hex digits. Leading
///   zeros are stripped.
///
/// @param s Output string (will be overwritten).
/// @param uppercase If true, uses A-F; otherwise a-f.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Converts the integer to a string in base 8, 10, or 16.
///
/// @details Main string conversion entry point. For base 10, the algorithm
///   repeatedly divides by 10^N (where N = maximum decimal digits that fit
///   in a word: 39 for 64-bit, 19 for 32-bit, 9 for 16-bit) and extracts
///   digit groups from the remainder. Key properties:
///
///   - Uses `nn_maxull` (largest representable unsigned value) as the
///     divisor for the inner loop.
///   - The volatile qualifier on the __uint128_t read (line ~721) is
///     CRITICAL: the GCC optimizer generates code that produces a runtime
///     segmentation fault at high optimization levels without it, due to
///     invalid register scheduling for the 128-bit load.
///   - Complexity: O(n^2) worst case for large integers with many digit
///     groups.
///
///   For bases 8 and 16, delegates to `base8string` and `base16string`
///   respectively.
///
///   The width parameter controls minimum digit count (zero-padded on the
///   left). If the string is wider than `width`, leading zeros are stripped.
///
/// @param width  Minimum output width (zero-padded).
/// @param base   Numeric base (8, 10, or 16).
/// @return  The formatted string.
// ---------------------------------------------------------------------------
inline const std::string integer::to_string(uintptr_t width, uintptr_t base) const
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
		integer m(abs());
		if ( m.is_zero() )
			return "0";
		t.reserve(imax(((proxy_->length_ * sizeof(word)) * 2375) / 1000, width));
		integer q;
		integer ten(10);
		do {
			m = m.divide(ten, &q);
			t.push_back(static_cast<char>(q.proxy_->data_[0] % 10) + '0');
		} while( !m.is_zero() );
		while( t.size() < width )
			t.push_back('0');
		std::reverse(t.begin(), t.end());
	}

	if( is_neg() )
		t.insert(0, "-");

	return t;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Returns the absolute value of this integer.
///
/// @details If the value is negative, returns the negation; otherwise
///   returns a copy.
///
/// @return The absolute value.
// ---------------------------------------------------------------------------
inline integer integer::abs() const
{
	return sign() < 0 ? -*this : *this;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Ensures exclusive ownership of the proxy data (copy-on-write).
///
/// @details If `ref_count_ > 1`, creates a new deep copy of the data array
///   and returns it as a new `nn_integer`. The caller must manage the
///   reference count of the returned proxy.
///
/// @return A new nn_integer with ref_count_ == 1.
// ---------------------------------------------------------------------------
inline integer integer::unique() const
{
	nn_integer result = nn_new(proxy_->length_);

	memcpy(result->data_, proxy_->data_, sizeof(word) * (result->length_ + 2));

	return result;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Checks whether the integer is zero.
///
/// @details First checks pointer equality with the global `nn_izero`
///   singleton. If the data is zero but the pointer differs, canonicalizes
///   by replacing `proxy_` with `nn_izero`. This ensures that all zeros
///   share a single proxy, improving pointer-equality fast-paths.
///
/// @return true if the value is zero.
// ---------------------------------------------------------------------------
inline constexpr bool integer::is_zero() const
{
	if( proxy_ == &nn_izero )
		return true;

	if( proxy_->z() ){
		nn_izero.add_ref();
		proxy_->release();
		proxy_ = const_cast<nn_integer>(&nn_izero);
		return true;
	}

	return false;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Checks whether the integer equals one.
///
/// @details Similar canonicalization to `is_zero()`: if the value is 1,
///   the proxy is replaced with the global `nn_ione` singleton.
///
/// @return true if the value is exactly 1.
// ---------------------------------------------------------------------------
inline constexpr bool integer::is_one() const
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
	proxy_ = const_cast<nn_integer>(&nn_ione);

	return true;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Checks whether the integer equals ten.
///
/// @details Canonicalization similar to `is_zero`/`is_one`. Replaces proxy_
///   with the global `nn_iten` singleton if the value is 10.
///
/// @return true if the value is exactly 10.
// ---------------------------------------------------------------------------
inline constexpr bool integer::is_ten() const
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
	proxy_ = const_cast<nn_integer>(&nn_iten);

	return true;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Computes GCD (НОД) and optionally LCM (НОК) of this and `a`.
///
/// @details Implements a simultaneous GCD/LCM algorithm using the identity:
///
///   gcd(x,y) * lcm(x,y) = x * y
///
///   The algorithm maintains two accumulator pairs (u,v) and (g,h) that
///   track intermediate GCD and LCM contributions respectively. At each
///   step, the larger of x,y is reduced by repeated subtraction of the
///   smaller (shifted by powers of 2 for efficiency).
///
///   Reference: Knuth, TAOCP Vol. 2, Section 4.5.2.
///
///   Edge cases:
///   - Division by zero in `a`: throws std::range_error
///   - If `*this` is zero, returns 0 and sets *p_nok = 0
///
/// @param a     The second operand.
/// @param p_nok Optional pointer to receive the LCM.
/// @return The GCD.
/// @throws std::range_error if `a` is zero.
// ---------------------------------------------------------------------------
inline integer integer::nod_nok(const integer & a, integer * p_nok) const
{
	if( a.is_zero() )
		//::divide(1, 0);
		throw std::range_error("Integer divide by zero");

	integer x(abs()), y(a.abs()), u, v, q, g;

	if( is_zero() ){
		x = 0;
		if( p_nok != nullptr ) *p_nok = 0;
	}
	else {
		intptr_t c;

		u = y;
		v = x;

		for( ;; ){
			c = x.compare(y);

			if( c > 0 ){
				if( p_nok != nullptr ) g = u;
				for( q = y; x - (q << 1) >= y; q <<= 1 ) if( p_nok != nullptr ) g <<= 1;
				x -= q;
				if( p_nok != nullptr ) v += g;
			}
			else if( c < 0 ){
				if( p_nok != nullptr ) g = v;
				for( q = x; y - (q << 1) >= x; q <<= 1 ) if( p_nok != nullptr ) g <<= 1;
				y -= q;
				if( p_nok != nullptr ) u += g;
			}
			else
				break;
		}

		if( p_nok != nullptr ) *p_nok = (u + v) >> 1;
	}
	return x;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Computes this^power using exponentiation by squaring.
///
/// @details Standard binary exponentiation (right-to-left):
///
///   result = 1
///   while power > 0:
///     if power is odd: result *= base
///     base *= base
///     power /= 2
///
///   Optimized special cases:
///   - power == 0: returns 1
///   - power == 1: returns *this
///   - power == 2: single multiplication (x*x)
///   - power == 3: two multiplications (x*x*x)
///   - power == 4: two multiplications with intermediate reuse
///   - Base == 10: uses repeated (x<<3)+(x<<1) which avoids multiplication
///     entirely (faster path for decimal power-of-ten computation)
///
///   Complexity: O(log(power)) multiplications. Each multiplication is
///   O(n^2) in the bit-size of the intermediate values.
///
/// @param power The exponent (non-negative).
/// @return this raised to the given power.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Static helper: computes (base)^power.
///
/// @details Convenience wrapper around the member `pow`.
///
/// @param power The exponent.
/// @param base  The base as a native integer.
/// @return base raised to power.
// ---------------------------------------------------------------------------
inline integer integer::pow(uintptr_t power, uintptr_t base)
{
	return integer(base).pow(power);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Computes the factorial of degree (degree!).
///
/// @details Simple iterative multiplication: result = 1*2*3*...*degree.
///   Results grow super-exponentially; for degree=1000, the result exceeds
///   2500 decimal digits. The algorithm uses O(degree) multiplications,
///   each on progressively larger integers, giving overall O(degree^3)
///   time in the bit-size (using schoolbook multiplication).
///
///   Stirling's approximation: n! ~ sqrt(2*pi*n) * (n/e)^n.
///
///   The commented-out `#if 0` block shows a planned caching strategy
///   that was never activated (factorial_cache_ is commented out in the
///   class declaration).
///
/// @param degree The factorial argument (non-negative).
/// @return degree!
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Stream output operator for integer.
///
/// @details Supports hex (0x prefix with showbase), octal (0 prefix with
///   showbase), and decimal output. For decimal output, respects showpos
///   and field width (via `out.width()`).
///
/// @param out The output stream.
/// @param v   The integer value.
/// @return Reference to the output stream.
// ---------------------------------------------------------------------------
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
		out << v.abs().to_string(uintptr_t(out.width()), 10);
	}
	return out;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Historical schoolbook multiplication implementation (commented out).
///
/// @details This block contains the original serial multiplication
///   implementation, preserved for historical reference. It uses the
///   standard O(n*m) schoolbook algorithm with a 16-entry lookup table
///   for small multipliers (b up to 16) and a nested-loop accumulator
///   for larger operands. The SSE2/AVX SIMD codepath (also commented out)
///   shows an attempted vectorization using _mm_mul_epu32.
///
///   The current active implementation (`operator *` below) uses the
///   threaded `tmul_a`/`tmul_v` dispatch for large operands, falling
///   back to `nn_integer_data::imul` for small ones.
///
///   Kept in `#if 0` per project convention — do not resurrect.
// ---------------------------------------------------------------------------
//inline integer integer::operator * (const integer & v) const
//{
//	stat_imul_++;
//#if 0
//	if( is_zero() || v.is_zero() )
//		return 0;
//	if( is_one() )
//		return v;
//	if( v.is_one() )
//		return *this;
//
//	if( proxy_->length_ < v.proxy_->length_ )
//		return v * *this;
//
//	// hard way, but simple
//	integer sum(0);
//	integer a(abs());
//	integer n(v.abs());
//	uintptr_t bit_size = n.proxy_->length_ * sizeof(word) * CHAR_BIT;
//
//	if( (is_neg() ^ v.is_neg()) != 0 ){
//		for( uintptr_t i = 0; i < bit_size; i++ )
//			if( n.bit(i) != 0 )
//				sum -= a << i;
//	}
//	else {
//		for( uintptr_t i = 0; i < bit_size; i++ )
//			if( n.bit(i) != 0 )
//				sum += a << i;
//	}
//
//	return sum;
//#endif
//	integer b(v.abs());
//
//	if( b.proxy_->length_ <= 2 ){
//		dword * p = (dword *) b.proxy_->data_;
//
//		switch( *p ){
//			case  0 : return integer(0);
//			case  1 : return *this;
//			case  2 : return *this << 1;
//			case  3 : return (*this << 1) + *this;
//			case  4 : return *this << 2;
//			case  5 : return (*this << 2) + *this;
//			case  6 : return (*this << 2) + (*this << 1);
//			case  7 : return (*this << 2) + (*this << 1) + *this;
//			case  8 : return *this << 3;
//			case  9 : return (*this << 3) + *this;
//			case 10 : return (*this << 3) + (*this << 1);
//			case 11 : return (*this << 3) + (*this << 1) + *this;
//			case 12 : return (*this << 3) + (*this << 2);
//			case 13 : return (*this << 3) + (*this << 2) + *this;
//			case 14 : return (*this << 3) + (*this << 2) + (*this << 1);
//			case 15 : return (*this << 3) + (*this << 2) + (*this << 1) + *this;
//			case 16 : return *this << 4;
//		}
//	}
//
//	integer a(abs());
//	uintptr_t jj = b.proxy_->length_, ii = a.proxy_->length_;
//	nn_integer r = nn_new(ii + jj + 1);
//
//	memset(r->data_,0,(r->length_ + 2) * sizeof(word));
//
//#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
//	register __m128i aa;
//	__m128i bb, cc;
//	uint32_t * bbp = (uint32_t *) &bb, * ccp = (uint32_t *) &cc;
//#endif
//
//	word * ap = a.proxy_->data_, * bp = b.proxy_->data_;
//	uintptr_t i = 0;
//
//	while( i < ii  ){
//#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
//		ccp[2] = ccp[0] = ap[i];
//		aa = _mm_load_si128(&cc);
//#else
//		dword h = ap[i];
//#endif
//		uintptr_t j = 0;
//
//		while( j < jj ){
//#if 0 && (__SSE2__ || _M_IX86_FP >= 2 || __AVX__ || __AVX2__) && SIZEOF_WORD == 4
//			bbp[0] = bp[j];
//			bbp[2] = bp[j + 1];
//			cc = _mm_mul_epu32(aa,bb);
//
//			uintptr_t index = i + j;
//			accumulate_with_carry(r, index + 0, ccp[0]);
//			accumulate_with_carry(r, index + 1, ccp[1]);
//			accumulate_with_carry(r, index + 1, ccp[2]);
//			accumulate_with_carry(r, index + 2, ccp[3]);
//			j += 2;
//#else
//			union {
//				struct {
//					word p0;
//					word p1;
//				};
//				dword q;
//			} d;
//			d.q = h * bp[j];
//			uintptr_t index = i + j;
//			// accumulate the the low bit sum
//			accumulate_with_carry(r, index, d.p0);
//			// accumulate the the high bit sum
//			accumulate_with_carry(r, index + 1, d.p1);
//			j++;
//#endif
//		}
//
//		i++;
//	}
//
//	r->normalize();
//
//	if( (is_neg() ^ v.is_neg()) != 0 )
//		return -integer(r);
//
//	return r;
//}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Global thread-pool and hardware-concurrency constant for
///   parallel multiplication.
///
/// @details The `pool` object is a globally-accessible `ThreadPoolT`
///   instance created with `hardware_concurrency()` worker threads. The
///   `effective_thread_count()` constant is evaluated once at program
///   start and used to decide whether parallel multiplication is
///   beneficial: only when the operand word-length >= concurrency count.
///
///   Destruction order: `pool` is destroyed at program exit (static
///   destruction) after all thread tasks have completed. No explicit
///   `join` is needed because `ThreadPoolT`'s destructor waits for
///   pending tasks.
// ---------------------------------------------------------------------------
inline unsigned effective_thread_count() {
    static unsigned cached = []{
        const char* env = std::getenv("NUMERIC_THREADS");
        if (env) {
            int val = std::atoi(env);
            return val > 0 ? static_cast<unsigned>(val) : 1u;
        }
        return std::thread::hardware_concurrency();
    }();
    return cached;
}
static ThreadPool pool(effective_thread_count());
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Performs schoolbook multiplication in parallel using a
///   stack-allocated accumulator array (template parameter SIZE).
///
/// @details Dispatches one thread per non-zero word of the left operand.
///   Each thread computes a partial product (b * m) shifted by `i` words
///   and stores the result in a separate accumulator. The accumulators are
///   allocated as a `std::array<nn_integer, SIZE>` on the stack, providing
///   fast allocation for typical operand sizes. Partial products are summed
///   sequentially after all threads complete.
///
///   Dispatch model:
///   - Only non-zero words of `a` are dispatched (skip zero words)
///   - Each thread runs `nn_integer_data::imul(b, i, m)` which computes
///     the word-at-a-time product and stores it in the accumulator
///   - All threads synchronize at `pool.wait()`
///   - Accumulation (sum or difference) is single-threaded
///
///   Complexity: O(n^2 / k) where n = words in the operand and k = number
///   of threads, plus O(n*k) for the final summation.
///
/// @tparam SIZE Maximum number of accumulators (stack allocation size).
///   Must be >= the word-length of `*this`.
/// @param v The right operand.
/// @return The product.
// ---------------------------------------------------------------------------
template <size_t SIZE>
inline integer integer::tmul_a(const integer & v) const
{
	integer s(0);
	// allocate accumulators
	std::array<nn_integer, SIZE> accums;
	uintptr_t ai = 0;

	try {
		//static std::random_device rd;
		//static std::uniform_int_distribution<int> dist(0, 255);

		integer a(abs()), b(v.abs());

		for( uintptr_t i = 0; i < a.proxy_->length_; i++ ) {
			word m = a.proxy_->data_[i];

			if( m == 0 )
				continue;

			nn_integer r = nn_new(a.proxy_->length_ + b.proxy_->length_ + 1);
			accums[ai++] = r;

			pool.run_task(std::bind([&] (nn_integer r, nn_integer b, uintptr_t i, word m) -> void {

				//std::unique_ptr<std::array<uint8_t, 1 * 1024 * 1024>> buf;
				//buf.reset(new decltype(buf)::element_type);

				//for( auto & i : *buf )
				//	i = dist(rd);

				//cdc256 ctx256;
				//ctx256.init();
				//ctx256.update(*buf);
				//auto digest256 = ctx256.final();

				//cdc512 ctx512;
				//ctx512.init();
				//ctx512.update(*buf);
				//auto digest512 = ctx512.final();

				//m = m;

				memset(r->data_, 0, (r->length_ + 2) * sizeof(word));
				r->imul(b, i, m);
			}, r, b.proxy_, i, m));
		}

		pool.wait();

		if( (is_neg() ^ v.is_neg()) != 0 ) {
			while( ai > 0 )
				s -= accums[--ai];
		}
		else {
			while( ai > 0 )
				s += accums[--ai];
		}
	}
	catch(...) {
		while( ai > 0 )
			accums[--ai]->release();
	}

	return s;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Performs schoolbook multiplication in parallel using a
///   heap-allocated accumulator vector.
///
/// @details Identical dispatch strategy to `tmul_a`, but uses
///   `std::vector<integer>` for the accumulator storage. This is used when
///   the operand word-length exceeds the stack-allocated SIZE threshold
///   (4096 in the current dispatch logic).
///
///   Using `std::vector` avoids stack overflow for extremely large
///   operands but adds heap allocation overhead for each accumulator.
///
/// @param v The right operand.
/// @return The product.
// ---------------------------------------------------------------------------
inline integer integer::tmul_v(const integer & v) const
{
	std::vector<integer> accums;

	integer a(abs()), b(v.abs());

	for( uintptr_t i = 0; i < a.proxy_->length_; i++ ) {
		word m = a.proxy_->data_[i];

		if( m == 0 )
			continue;

		nn_integer r = nn_new(a.proxy_->length_ + b.proxy_->length_ + 1);

		accums.emplace_back(r);

		pool.run_task(std::bind([&] (nn_integer r, nn_integer b, uintptr_t i, word m) -> void {
			memset(r->data_, 0, (r->length_ + 2) * sizeof(word));
			r->imul(b, i, m);
		}, r, b.proxy_, i, m));
	}

	pool.wait();

	integer s(0);

	if( (is_neg() ^ v.is_neg()) != 0 ) {
		for( const auto & i : accums )
			s -= i;
	}
	else {
		for( const auto & i : accums )
			s += i;
	}

	return s;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Multiplies two integers, dispatching to parallel or serial
///   implementation based on operand size.
///
/// @details High-level dispatch:
///   1. If `effective_thread_count() > 1` and `proxy_->length_ >=
///      effective_thread_count()`, use threaded multiplication:
///      - `tmul_a<4096>` for operands up to 4096 words
///      - `tmul_v` for larger operands (heap-allocated accumulators)
///   2. Otherwise, fall back to serial `nn_integer_data::imul` through
///      the single-threaded code path.
///
///   Sign handling: the product of absolute values is computed, then the
///   sign is negated if exactly one operand is negative.
///
///   Complexity: O(n^2 / k) in parallel, O(n^2) serial.
///
/// @param v The right operand.
/// @return The product.
// ---------------------------------------------------------------------------
inline integer integer::operator * (const integer & v) const
{
	stat_imul_++;

	if( effective_thread_count() > 1 && proxy_->length_ >= effective_thread_count() ) {
        if( proxy_->length_ <= 4096 )
            return tmul_a<4096>(v);
		return tmul_v(v);
	}

	integer a(abs()), b(v.abs());

	nn_integer r = nn_new(a.proxy_->length_ + b.proxy_->length_ + 1);

	memset(r->data_, 0, (r->length_ + 2) * sizeof(word));
	r->imul(a.proxy_, b.proxy_);

	if( (is_neg() ^ v.is_neg()) != 0 )
		return -integer(r);

	return r;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Division operator — delegates to `divide()`.
///
/// @param v The divisor.
/// @return The quotient.
/// @throws std::range_error if v is zero.
// ---------------------------------------------------------------------------
inline integer integer::operator / (const integer & v) const
{
	stat_idiv_++;
	return divide(v);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Divides this integer by `divider`, returning quotient and
///   optionally the remainder.
///
/// @details Implements the restoring division algorithm (schoolbook long
///   division). For each bit of the numerator (processed from MSB to LSB):
///
///   1. Left-shift the remainder by 1
///   2. Set the LSB of the remainder to the current numerator bit
///   3. If remainder >= divisor, subtract divisor and set corresponding
///      quotient bit to 1
///
///   Reference: Knuth, TAOCP Vol. 2, Section 4.3.1, Algorithm D (Knuth's
///   long division). This implementation uses the simpler restoring variant
///   rather than the more efficient non-restoring Algorithm D, trading
///   some performance for clarity.
///
///   Complexity: O(n * m) time where n = numerator bits, m = divisor bits.
///   Space: O(max(n,m)) for the quotient.
///
///   Edge cases:
///   - Division by zero: throws std::range_error
///   - Division by one: optimized to return *this
///   - Zero numerator: returns 0
///
///   Name note: named `divide` rather than `div` to avoid shadowing
///   `std::div` from <cstdlib>.
///
/// @param divider The divisor (must not be zero)
/// @param p_mod   Optional pointer to receive the remainder
/// @return The quotient
/// @throws std::range_error if divider is zero
// ---------------------------------------------------------------------------
inline integer integer::divide(const integer & divider, integer * p_mod) const
{
	if( divider.is_zero() )
		//::divide(1, 0);
		throw std::range_error("Integer divide by zero");

	if( divider.is_one() )
		return *this;

	integer temp_r, & r = p_mod == nullptr ? temp_r : *p_mod;

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
		r.sbit(0, n.bit(i));

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
} // namespace NaturalNumbers
//------------------------------------------------------------------------------
#endif // II_HPP_INCLUDED
//------------------------------------------------------------------------------
