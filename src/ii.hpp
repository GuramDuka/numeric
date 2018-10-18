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
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
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
		integer(const nn_integer_data * p, int) : proxy_(p->add_ref()) {}

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
		void sbit(uintptr_t i, const T v = 1) const {
			if( proxy_->ref_count_ > 1 || i >= proxy_->length_ * sizeof(word) * CHAR_BIT ){
				uintptr_t new_bit_size = (i + 1) + (-intptr_t(i + 1) & (sizeof(word) * CHAR_BIT - 1));
				uintptr_t new_size = new_bit_size / (sizeof(word) * CHAR_BIT);

				nn_integer result = nn_new(imax(new_size, proxy_->length_));

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
		proxy_ = const_cast<nn_integer>(&nn_izero);
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
	proxy_ = const_cast<nn_integer>(&nn_ione);

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
	proxy_ = const_cast<nn_integer>(&nn_iten);

	return true;
}
//------------------------------------------------------------------------------
inline integer integer::nod_nok(const integer & a, integer * p_nok) const
{
	if( a.is_zero() )
		//::div(1, 0);
		throw std::range_error("Integer divide by zero");

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
		out << v.abs().to_string(uintptr_t(out.width()), 10);
	}
	return out;
}
//------------------------------------------------------------------------------
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
static const auto thread_hardware_concurrency = std::thread::hardware_concurrency();
static ThreadPool pool(thread_hardware_concurrency);
//------------------------------------------------------------------------------
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
inline integer integer::operator * (const integer & v) const
{
	stat_imul_++;

	if( thread_hardware_concurrency > 1 && proxy_->length_ >= thread_hardware_concurrency ) {
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
inline integer integer::operator / (const integer & v) const
{
	stat_idiv_++;
	return div(v);
}
//------------------------------------------------------------------------------
inline integer integer::div(const integer & divider, integer * p_mod) const
{
	if( divider.is_zero() )
		//::div(1, 0);
		throw std::range_error("Integer divide by zero");

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
