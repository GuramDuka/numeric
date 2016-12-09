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
#ifndef NN_HPP_INCLUDED
#define NN_HPP_INCLUDED
//------------------------------------------------------------------------------
#include "ii.hpp"
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
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

		numeric & operator ++ () {
			return *this += integer(1);
		}

		// postfix form
		numeric operator ++ (int) {
			return *this += integer(1);
		}

		numeric & operator -- () {
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
		//::div(1, 0);
		throw std::range_error("Numeric divide by zero");

	numeric t(numerator_.abs(), denominator_ * v.abs());

	if( (numerator_.is_neg() ^ v.is_neg()) != 0 )
		t.numerator_ = -t.numerator_;

	return t;
}
//------------------------------------------------------------------------------
inline numeric numeric::operator / (const numeric & v) const
{
	if( v.numerator_.is_zero() )
		//::div(1, 0);
		throw std::range_error("Numeric divide by zero");

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
#endif // NN_HPP_INCLUDED
//------------------------------------------------------------------------------
