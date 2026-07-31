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
// ---------------------------------------------------------------------------
/// @namespace nn (NaturalNumbers)
/// @brief Arbitrary-precision rational arithmetic with transcendental functions
///
/// The `numeric` class represents numbers as exact rational fractions:
///   value = numerator / denominator
/// where both numerator and denominator are arbitrary-precision integers
/// (class `integer`, defined in ii.hpp).
///
/// This design preserves exactness through all arithmetic operations — unlike
/// IEEE 754 floating-point, rational arithmetic introduces no representation
/// error for rational inputs. The price is potential numerator/denominator
/// growth after each operation, which is mitigated by the normalize() method.
///
/// ### Transcendental functions
/// The class provides Taylor-series implementations of sin, cos, and atan,
/// Newton's method for square root, bisection for nth root, the BBP formula
/// for π, and exponentiation-by-squaring for integer powers.
///
/// ### Fraction growth warning
/// Exact rational Newton iteration (root with power=2) doubles numerator and
/// denominator bit size every iteration. At iter ≥ 20, million-bit integers
/// appear, causing GCD normalization (normalize(3)) to grind to a halt.
/// Keep Newton iterations ≤ 12.
///
/// References:
///   - Knuth, D.E. "The Art of Computer Programming", Vol. 2:
///     Seminumerical Algorithms, 3rd ed., Addison-Wesley, 1998.
///   - Bailey, D.H., Borwein, P.B., Plouffe, S. "On the Rapid Computation
///     of Various Polylogarithmic Constants", Mathematics of Computation,
///     Vol. 66, No. 218, pp. 903–913, 1997.
///   - Newton, I. "Methodus Fluxionum et Serierum Infinitarum" (1671).
///   - Madhava of Sangamagrama (14th c.) / Gregory (1671) / Leibniz (1676).
///     The Madhava–Gregory series for arctan and the related sine/cosine
///     expansions originate in the Kerala school of astronomy and mathematics.
///   - Brent, R.P. "Fast Multiple-Precision Evaluation of Elementary
///     Functions", J. ACM, Vol. 23, No. 2, pp. 242–251, 1976.
// ---------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------
// //////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
/// @class numeric
/// @brief Arbitrary-precision rational number.
///
/// @details Represents a rational number as a pair of arbitrary-precision
///   integers: `numerator_ / denominator_`. The denominator is always
///   positive after construction; the sign of the entire number is carried
///   by the numerator.
///
///   **Normalization**: After most arithmetic operations, normalize(3) is
///   called to reduce the fraction. normalize(3) first trims common trailing
///   zero bits (fast path, bit 1), then reduces by GCD (slower canonical
///   form, bit 2). Use normalize(1) for a fast approximate reduction when
///   exact canonical form is not yet needed.
///
///   **Conversions**: Construction from `std::string` parses decimal notation
///   (e.g., "-3.14159"). Construction from `double`/`long double` converts
///   the floating-point value to a rational by extracting the fractional part
///   as a power-of-10 fraction. Output via `to_string()` or `operator<<`.
///
///   **Known pitfalls**:
///   1. `numeric("-3.14159")` historically gave `-2.85841` due to incorrect
///      sign handling of the fractional part. Fixed: the fraction is now
///      subtracted (not added) when ipart is negative.
///   2. `integer + numeric` converts numeric to integer via truncation.
///      Always wrap as `numeric(ipart) + numeric(q, p)` to stay in rational
///      domain. This broke round() historically.
///   3. normalize(3) after Newton iteration on large rationals can be
///      extremely slow — see the fraction blowup warning above.
// ---------------------------------------------------------------------------
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

		/// @brief Test whether the value is negative.
		/// @return true if numerator_ < 0 (the denominator is always positive).
		constexpr bool is_neg() const {
			return numerator_.is_neg();
		}

		/// @brief Test whether the value is zero.
		/// @return true if numerator_ == 0 (denominator is irrelevant for zero).
		constexpr bool is_zero() const {
			return numerator_.is_zero();
		}

		/// @brief Test whether the value is exactly one.
		/// @return true if numerator_ == denominator_ == 1.
		constexpr bool is_one() const {
			return numerator_.is_one() && denominator_.is_one();
		}

		/// @brief Absolute value (magnitude only, sign stripped).
		/// @return `numeric(numerator_.abs(), denominator_)` — the magnitude
		///   with the same denominator.
		numeric abs() const {
			return numeric(numerator_.abs(), denominator_);
		}

		numeric mod(integer * p_ipart = nullptr) const;
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
		integer numerator_;   ///< Numerator (carries the sign of the whole number)
		integer denominator_; ///< Denominator (always positive after construction)
};
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Convert numeric to integer via truncation toward zero.
///
/// Constructs an integer equal to `numerator_ / denominator_` using integer
/// division (truncation toward zero). This is the conversion used by
/// `integer + numeric` and similar mixed-type operations.
///
/// @warning This truncates — use mod() or modf-style extraction to get the
///   fractional part separately.
// ---------------------------------------------------------------------------
inline integer::integer(const numeric & v) : integer(v.numerator_.divide(v.denominator_))
{
}
//------------------------------------------------------------------------------
inline integer & integer::operator = (const numeric & v)
{
	return *this = v.numerator_.divide(v.denominator_);
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
// ---------------------------------------------------------------------------
/// @brief Construct a numeric from a decimal string.
///
/// Parses strings of the form `[+-]?[0-9]*(\.[0-9]+)?` (e.g., `"3.14159"`,
/// `"-2.718"`, `"42"`).
///
/// ### Algorithm
/// 1. Scan from right to left, building the fractional part first (if a
///    decimal point is found), then the integer part.
/// 2. On encountering a decimal point, switch from fraction-building mode
///    to integer-building mode, reset the place multiplier to 1.
/// 3. On encountering '-' or '+', apply the sign to the integer part and
///    stop scanning.
/// 4. After scanning, compute the combined value:
///    @code
///      numerator_ = ipart;
///      if (fractional_part_exists) {
///          denominator_ = 10^(number_of_fraction_digits);
///          numerator_ = ipart * denominator_ ± fraction;
///          // '+' if ipart >= 0, '-' if ipart < 0
///      }
///    @endcode
///
/// ### Historical bug (fixed)
/// Earlier versions always added the fraction, even for negative numbers:
///   `numeric("-3.14159")` produced `-3 + 0.14159 = -2.85841` (wrong).
/// Fixed by conditionally subtracting the fraction when `ipart.is_neg()`.
///
/// @param s Decimal string to parse.
/// @throws std::invalid_argument if the string contains unexpected characters.
// ---------------------------------------------------------------------------
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
		if( ipart.is_neg() )
			numerator_ -= fraction;
		else
			numerator_ += fraction;
	}
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Construct a numeric from a floating-point value.
///
/// Converts a `double` or `long double` to an exact rational by:
/// 1. Splitting into integer and fractional parts via `modf()`/`modfl()`.
/// 2. Multiplying the fraction by successive powers of 10 until it becomes
///    an integer (i.e., the fractional part vanishes).
/// 3. Constructing `numerator_ = ipart * 10^power + fraction` over
///    `denominator_ = 10^power`.
///
/// This produces the exact rational representation of the floating-point
/// value, not of the real number it approximates. For example,
/// `numeric(0.1)` produces `1/10` only if `0.1` is exactly representable
/// in the floating-point format — on IEEE 754 systems, this is not the case
/// and `numeric(0.1)` will reflect the binary approximation.
///
/// @param v Floating-point value to convert.
// ---------------------------------------------------------------------------
#if _MSC_VER || HAVE_LONG_DOUBLE
inline numeric::numeric(double v) : numeric(static_cast<long double>(v)) {}
inline numeric::numeric(long double v) : numeric()
#elif defined(LONG_DOUBLE)
inline numeric::numeric(double v) : numeric(static_cast<LONG_DOUBLE>(v)) {}
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
	numerator_ = integer(static_cast<double>(ipart));
#else
	numerator_ = integer(ipart);
#endif

	if( power > 0 ){
		denominator_ = integer(10).pow(power);
		numerator_ *= denominator_;
#if __GNUC__
		numerator_ += integer(static_cast<double>(fraction));
#else
		numerator_ += integer(fraction);
#endif
	}
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Construct a numeric as an exact fraction numerator/denominator.
///
/// @warning Do NOT call normalize() here. The comparison operators
///   (operator>, operator==, etc.) rely on comparing cross-multiplied values,
///   which works correctly even on non-normalized fractions. Normalizing
///   eagerly would degrade performance by doing unnecessary GCD computations.
///
/// @param numerator   The numerator (carries the sign).
/// @param denominator The denominator (should be positive — no sign flip is
///                    performed automatically).
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Addition: numeric + integer.
///
/// Converts the integer to rational form (integer / 1) with the same
/// denominator, then adds:
///   a + b = (a + b * d_a) / d_a
///
/// @note The result is NOT automatically normalized — the caller may choose
///   when to call normalize() based on subsequent operations.
///
/// @param v Integer addend.
/// @return `*this + v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator + (const integer & v) const
{
	return numeric(numerator_ + v * denominator_, denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Addition: numeric + numeric.
///
/// Cross-multiplication:
///   a/b + c/d = (ad + bc) / bd
///
/// Both numerator and denominator grow by roughly one multiplication each.
/// The denominator b*d is computed but may share common factors with the
/// numerator; normalize(3) (called during assignment or explicitly) will
/// reduce these.
///
/// @param v Rational addend.
/// @return `*this + v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator + (const numeric & v) const
{
	return numeric(numerator_ * v.denominator_ + v.numerator_ * denominator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Subtraction: numeric - integer.
///
///   a - b = (a - b * d_a) / d_a
///
/// @param v Integer subtrahend.
/// @return `*this - v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator - (const integer & v) const
{
	return numeric(numerator_ - v * denominator_, denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Subtraction: numeric - numeric.
///
/// Cross-multiplication:
///   a/b - c/d = (ad - bc) / bd
///
/// @param v Rational subtrahend.
/// @return `*this - v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator - (const numeric & v) const
{
	return numeric(numerator_ * v.denominator_ - v.numerator_ * denominator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Multiplication: numeric * integer.
///
///   (a/b) * c = (a*c) / b
///
/// Only the numerator grows — denominator stays unchanged.
///
/// @param v Integer multiplier.
/// @return `*this * v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator * (const integer & v) const
{
	return numeric(numerator_ * v, denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Multiplication: numeric * numeric.
///
///   (a/b) * (c/d) = (a*c) / (b*d)
///
/// Both numerator and denominator grow by one multiplication each. This is
/// the cheapest cross-operation in terms of temporary size: no addition is
/// needed, only two multiplications.
///
/// @param v Rational multiplier.
/// @return `*this * v` as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::operator * (const numeric & v) const
{
	return numeric(numerator_ * v.numerator_,
		denominator_ * v.denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Division: numeric / integer.
///
///   (a/b) / c = a / (b*c)
///
/// Division by zero throws `std::range_error`. Sign is handled via XOR of
/// the sign bits: negative/positive → negative, negative/negative → positive.
///
/// @param v Integer divisor.
/// @return `*this / v` as a rational.
/// @throws std::range_error if v is zero.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Division: numeric / numeric.
///
///   (a/b) / (c/d) = (a*d) / (b*c)
///
/// Division by zero throws `std::range_error`. Sign is handled via XOR of
/// the sign bits.
///
/// @param v Rational divisor.
/// @return `*this / v` as a rational.
/// @throws std::range_error if v.numerator_ is zero.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Unary negation.
///
/// @return `numeric(-numerator_, denominator_)` — only the numerator sign flips.
// ---------------------------------------------------------------------------
inline numeric numeric::operator - () const
{
	return numeric(-numerator_, denominator_);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Greater-than comparison
///
/// Uses the common-denominator approach:
///   a/b > c/d  ⇔  a*d > c*b
///
/// Two temporary integers `u.numerator_ = a*d` and `n.numerator_ = c*b` are
/// created (via intermediate `numeric` objects). This allocates large
/// temporaries — the cross-multiplication integers can be up to twice the
/// bit-size of the operands. This is by design: no GCD reduction is needed
/// for comparison, and the common-denominator method is exact.
///
/// @param v Right-hand operand.
/// @return true if *this > v.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief OLD normalize (disabled) — threshold-based variant.
///
/// This is the original `normalize(uintptr_t threshold, int how)` that
/// accepted a word-length threshold. It is preserved for reference and
/// regression testing. The active version below uses `normalize(int how)`
/// without the threshold parameter.
///
/// The threshold acted as a gate: normalization was skipped if both
/// numerator and denominator were below `threshold` words in length.
/// This was removed to make normalization behavior uniform regardless
/// of operand size — the caller controls normalization explicitly via
/// the `how` parameter.
// ---------------------------------------------------------------------------
#if 0
inline numeric & numeric::normalize(uintptr_t threshold, int how)
{
	uintptr_t nl = numerator_.proxy_->length_;
	uintptr_t dl = denominator_.proxy_->length_;

	if( (how & 1) != 0 && (threshold == 0 || nl >= threshold || dl >= threshold) ){
		// fast path
		uintptr_t i = 0;
		uintptr_t shift = 0;
		const word * pw_num = numerator_.proxy_->data_;
		const word * pw_den = denominator_.proxy_->data_;

		while( i < nl && i < dl ){
			if( *pw_num != 0 || *pw_den != 0 ) break;
			shift += sizeof(word) * CHAR_BIT;
			pw_num++;
			pw_den++;
			i++;
		}

		i *= sizeof(word);
		nl *= sizeof(word);
		dl *= sizeof(word);

		{
			const uint8_t * pb_num = reinterpret_cast<const uint8_t *>(pw_num);
			const uint8_t * pb_den = reinterpret_cast<const uint8_t *>(pw_den);

			while( i < nl && i < dl ){
				if( *pb_num != 0 || *pb_den != 0 ) break;
				shift += sizeof(uint8_t) * CHAR_BIT;
				pb_num++;
				pb_den++;
				i++;
			}
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
			integer nod(numerator_.nod_nok(denominator_, nullptr));

			if( !nod.is_zero() && !nod.is_one() ){
				numerator_ /= nod;
				denominator_ /= nod;
			}
		}
	}
#else
// ---------------------------------------------------------------------------
/// @brief Normalize the fraction to canonical form.
///
/// normalize() is the core post-arithmetic reduction routine. It operates
/// in two independent modes, selected by bits of the `how` parameter:
///
/// | how bits | Behavior |
/// |----------|----------|
/// | 1 (bit 0) | Trim common trailing zero bits from numerator and denominator (fast — no GCD needed) |
/// | 2 (bit 1) | Reduce by GCD: divide both numerator and denominator by gcd(numerator, denominator) |
/// | 3 (bits 0+1) | Both: trim trailing zeros, then GCD-reduce (standard post-arithmetic call) |
///
/// ### Bit 1 — Trailing zero trimming (fast path)
/// Scans the least-significant words and bytes of both numerator and
/// denominator to find bits that are zero in both. Shifts both right by
/// the common count. This is a cheap "poor man's GCD" that eliminates
/// factors of 2 without a full GCD computation. It is always safe to call
/// this before a full GCD reduction because it can only reduce the size
/// of the GCD inputs.
///
/// ### Bit 2 — GCD reduction (slow path)
/// Computes `d = gcd(numerator, denominator)` using the binary GCD
/// algorithm from ii.hpp. If d > 1, divides both numerator and denominator
/// by d, producing the canonical reduced fraction.
///
/// ### Performance characteristics
/// - normalize(1): O(n) in the number of words (trailing zero scan)
/// - normalize(2): O(n^2) worst-case (GCD on large integers)
/// - normalize(3): Both, but bit-1 reduction typically shrinks the GCD inputs
///
/// ### When to call
/// - normalize(3) after every arithmetic operation that creates new fractions
/// - normalize(1) between intermediate steps of a multi-step computation
///   (e.g., inside Taylor series summation) to keep sizes manageable
/// - normalize(2) when canonical form is needed (e.g., for output)
///
/// ### Known pitfall
/// normalize(3) after Newton iteration on large rationals (root with
/// power=2, iter ≥ 20) can be extremely slow because the numerator and
/// denominator may be million-bit integers. See AGENTS.md.
///
/// @param how Bitmask controlling which normalization steps to perform:
///            1 = trailing zero trim, 2 = GCD reduce, 3 = both.
/// @return Reference to *this (for chaining).
// ---------------------------------------------------------------------------
inline numeric & numeric::normalize(int how)
{
	if( how & 2 ){
		integer d(numerator_.gcd(denominator_));
		//integer nod(numerator_.nod_nok(denominator_, nullptr));

		if( !d.is_zero() ){
			numerator_ /= d;
			denominator_ /= d;
		}
	}

	if( how & 1 ){
		uintptr_t i = 0;
		uintptr_t shift = 0;
		const word * pw_num = numerator_.proxy_->data_;
		const word * pw_den = denominator_.proxy_->data_;

		uintptr_t nl = numerator_.proxy_->length_;
		uintptr_t dl = denominator_.proxy_->length_;

		// Phase 1: Scan whole words (aligned access)
		while( i < nl && i < dl ){
			if( *pw_num != 0 || *pw_den != 0 ) break;
			shift += sizeof(word) * CHAR_BIT;
			pw_num++;
			pw_den++;
			i++;
		}

		// Phase 2: Scan individual bytes within the partial word
		i *= sizeof(word);
		nl *= sizeof(word);
		dl *= sizeof(word);

		{
			const uint8_t * pb_num = reinterpret_cast<const uint8_t *>(pw_num);
			const uint8_t * pb_den = reinterpret_cast<const uint8_t *>(pw_den);

			while( i < nl && i < dl ){
				if( *pb_num != 0 || *pb_den != 0 ) break;
				shift += sizeof(uint8_t) * CHAR_BIT;
				pb_num++;
				pb_den++;
				i++;
			}
		}

		// Phase 3: Scan individual bits (fine-grained)
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
// ---------------------------------------------------------------------------
/// @brief Extract integer and fractional parts (mod / trunc toward zero).
///
/// Computes `integer_part = floor(*this)` toward zero (i.e., truncation).
/// This is equivalent to C's `div` for positive numbers, but for negative
/// numbers it truncates toward zero (not toward -infinity).
///
/// More precisely:
/// - `integer_part = numerator_ / denominator_` (integer division, truncates
///   toward zero)
/// - `fractional_part = *this - integer_part` (a rational with the original
///   denominator, but may be negative or zero)
///
/// ### Pre/post conditions
/// - If `p_ipart` is non-null, `*p_ipart` is set to the integer part.
/// - The return value is the fractional part as a rational: the remainder
///   in [0, 1) for non-negative inputs, or (-1, 0] for non-positive inputs,
///   exactly zero for integer inputs.
/// - `*this == integer_part + fractional_part` holds exactly.
///
/// @param p_ipart If non-null, receives the integer part (truncated toward zero).
/// @return The fractional part (remainder) as a `numeric`.
// ---------------------------------------------------------------------------
inline numeric numeric::mod(integer * p_ipart) const
{
	integer ipart(numerator_.divide(denominator_));

	if( p_ipart != nullptr )
		*p_ipart = ipart;

	return *this - ipart;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Integer power (exponentiation by squaring).
///
/// Computes `*this^power` for any integer exponent (including negative).
/// Uses binary exponentiation (exponentiation by squaring) which requires
/// O(log |power|) multiplications.
///
/// ### Algorithm
/// ```
/// result = 1
/// base = *this
/// n = |power|
/// while n != 0:
///     if n is odd:
///         result *= base
///         n -= 1
///     base *= base
///     n /= 2
/// if power < 0:
///     result = 1 / result
/// return result
/// ```
///
/// ### Optimizations
/// Small powers (0, 1, 2, 3, 4) have direct fast paths that avoid loop
/// overhead and temporary objects.
///
/// ### Negative powers
/// For negative `power`, the result is the reciprocal: `x^(-p) = 1 / x^p`.
/// This is exact for rational x, producing a rational result.
///
/// @param power Integer exponent (may be negative).
/// @return `*this^power` as a rational.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Round to a specified number of decimal digits.
///
/// Rounds the number to `digits` decimal places, using either floor
/// (toLess) or ceiling (toGreater) rounding.
///
/// ### Algorithm
/// 1. Split into integer part `ipart` and fractional part via mod().
/// 2. Scale the fractional part by `10^digits` and split again into
///    integer `q` and remainder.
/// 3. Adjust `q` based on the rounding mode:
///    - toLess (round toward -∞): if remainder < 0, decrement q
///    - toGreater (round toward +∞): if remainder > 0, increment q
/// 4. Recombine: `result = ipart + q / 10^digits`
///
/// ### Sign-aware rounding
/// The remainder inherits the sign of the fractional part. For negative
/// numbers with a negative remainder, toLess decrements (making it more
/// negative), correctly implementing "round toward -infinity".
///
/// @param digits Number of decimal digits to round to.
/// @param type   Rounding direction: toLess (floor) or toGreater (ceil).
/// @return The rounded value as a rational.
// ---------------------------------------------------------------------------
inline numeric numeric::round(uintptr_t digits, round_type type) const
{
	integer ipart, q, p(integer(10).pow(digits));
	numeric fraction(mod(&ipart));

	fraction *= p;

	// q = integer part of shifted value (truncated toward zero)
	// remainder = fractional part (may be negative)
	numeric remainder = fraction.mod(&q);

	// toLess  → round toward -infinity
	// toGreater → round toward +infinity
	if( !remainder.is_zero() ){
		if( type == toLess && remainder.is_neg() )
			q -= 1;   // negative remainder → go more negative
		else if( type == toGreater && !remainder.is_neg() )
			q += 1;   // positive remainder → go more positive
	}

	return numeric(ipart) + numeric(q, p);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Sine via Taylor series (Madhava–Gregory–Leibniz).
///
/// Computes sin(x) using the Taylor series expansion:
///
///   sin(x) = Σ_{n=0}^{∞} (-1)^n * x^(2n+1) / (2n+1)!
///          = x - x^3/3! + x^5/5! - x^7/7! + ...
///
/// ### Convergence
/// This is an alternating series for any finite x. The error after `iter`
/// terms is bounded by the absolute value of the first omitted term
/// (Leibniz criterion for alternating series).
///
/// ### Implementation details
/// - Factorials are computed incrementally and cached — `f(degree)` extends
///   the running product `factorial` from `cur_degree` up to `degree`.
/// - Powers of x are similarly cached — `p(power)` extends `powerv` from
///   `cur_pow` up to `power` by repeated multiplication with `*this`.
/// - Each term is normalized with `normalize(3)` to keep intermediate
///   fractions from growing unboundedly.
///
/// ### Performance
/// - O(iter) multiplications and additions on exact rationals.
/// - Each term involves one factorial (integer) and one power (rational).
/// - Intermediate fractions grow as iter increases; normalize(3) controls
///   this at O(iter) cost.
///
/// @param iter Number of terms to compute (default: 2). More terms = more
///   precision at the cost of O(iter) growth.
/// @return Approximation of sin(*this) as a rational.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Cosine via Taylor series (Madhava–Gregory–Leibniz).
///
/// Computes cos(x) using the Taylor series expansion:
///
///   cos(x) = Σ_{n=0}^{∞} (-1)^n * x^(2n) / (2n)!
///          = 1 - x^2/2! + x^4/4! - x^6/6! + ...
///
/// ### Convergence
/// Like sin(), this is an alternating series for any finite x. Error after
/// `iter` terms is bounded by the first omitted term.
///
/// ### Implementation
/// Shares the same factorial-caching and power-caching closures as sin().
/// The only difference is that cos evaluates even powers and factorials:
/// term_n = x^(2n) / (2n)! (vs. x^(2n+1) / (2n+1)! for sin).
///
/// @param iter Number of terms to compute (default: 2).
/// @return Approximation of cos(*this) as a rational.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Arctangent via Madhava–Gregory series.
///
/// Computes atan(x) using the Taylor series:
///
///   atan(x) = Σ_{n=0}^{∞} (-1)^n * x^(2n+1) / (2n+1)
///           = x - x^3/3 + x^5/5 - x^7/7 + ...
///
/// This series converges for |x| ≤ 1. For |x| > 1, convergence is slow
/// or divergent; the caller should use identities like
/// atan(x) = π/2 - atan(1/x) to reduce the argument.
///
/// ### Implementation
/// Uses an iterative Horner-like scheme:
/// - Start with `pow2np1 = |x|`, `n2np1 = 1`
/// - Precompute `pow2 = x^2` (the common multiplier for each step)
/// - For each n: add/subtract `pow2np1 / n2np1`, then update
///   `pow2np1 *= pow2` and `n2np1 += 2`
///
/// This avoids per-term re-computation of the full power and division.
///
/// ### Performance
/// O(iter) multiplications and divisions on exact rationals. The series
/// converges linearly — each term adds about one bit of precision for
/// |x| not close to 1.
///
/// ### Historical note
/// This series was discovered by Madhava of Sangamagrama (14th century),
/// later independently by Gregory (1671) and Leibniz (1676). It is a
/// special case of the more general Taylor expansion for arctan.
///
/// @param iter Number of terms to compute (default: 2).
/// @return Approximation of atan(*this) as a rational.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief nth root computation (dual method: Newton for sqrt, bisection otherwise).
///
/// Computes the nth root of *this: `⁡√(*this)`. Uses one of two methods
/// depending on the root power:
///
/// | power | Method | Convergence | Fraction growth |
/// |-------|--------|-------------|-----------------|
/// | 2     | Newton's method | Quadratic (bit count doubles each iteration) | ⚠️ Exponential! |
/// | ≠ 2   | Bisection | Linear (half interval per iteration) | None |
///
/// ### Newton's method (power = 2)
/// The iteration is:
///   x_{n+1} = 0.5 * (x_n + a/x_n)
///
/// where a = *this, and x_0 = a. This converges quadratically: the number
/// of correct bits doubles at each step.
///
/// **⚠️ CRITICAL WARNING**: Exact rational Newton causes exponential
/// numerator/denominator growth. Each iteration DOUBLES the bit size of
/// both numerator and denominator. At iter=20, million-bit integers appear
/// and GCD normalization (normalize(3)) grinds to a halt. Keep iter ≤ 12.
///
/// Convergence guards:
/// - `x * x == r` — exact result found (rare with exact rationals)
/// - `x == prev_x` — fixed point reached (no further progress)
/// - `x.to_string() == prev_x.to_string()` — string-level convergence
///   (same decimal representation)
///
/// ### Bisection method (power ≠ 2)
/// Binary search in [0, r]:
/// ```
/// a = 0, b = r (or abs(r) for odd roots)
/// for n in [0, iter):
///     m = a + (b - a) / 2
///     if m^power < r: a = m
///     elif m^power > r: b = m
///     else: break
/// result = m
/// ```
///
/// This has no fraction blowup problem — each iteration divides the interval
/// in half, and the midpoint computation does not compound fraction size.
///
/// ### Sign handling
/// For odd powers (3, 5, 7, ...), negative inputs yield negative roots.
/// For even powers (2, 4, 6, ...), only the absolute value is used.
///
/// @param power Root degree (1 = identity, 2 = square root, 3 = cube root, ...).
/// @param iter  Maximum number of iterations (default: 16).
/// @return Approximation of the nth root as a rational.
// ---------------------------------------------------------------------------
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
		// Newton's method: x_{n+1} = 0.5 * (x_n + a/x_n)
		// WARNING: exact rational Newton causes exponential fraction growth.
		// Use modest iteration counts (<12) or switch to floating-point mode.
		numeric x(r), prev_x, hone(1, 2);

		for( uintptr_t n = 0; n < iter; n++ ){
			prev_x = x;
			x = hone * (x + (r / x));
			x.normalize(3);
			// Exact rational convergence check
			if( x * x == r || x == prev_x )
				break;
			// Guard against fraction blowup: break when precision exceeds needed
			if( x.to_string() == prev_x.to_string() )
				break;
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
// ---------------------------------------------------------------------------
/// @brief Convert to decimal string.
///
/// Formats the rational as a decimal string with the specified width and
/// precision:
///   `[sign]integer_part[.fractional_part]`
///
/// ### Algorithm
/// 1. Extract integer part via abs().mod(&ipart).
/// 2. Format integer part with ipart.to_string(width).
/// 3. Compute fractional part by scaling the remainder by 10^precision
///    and extracting the integer part of the scaled fraction.
/// 4. Zero-pad the fractional part to the requested precision.
/// 5. Trim trailing zeros if precision was not explicitly specified (i.e.,
///    precision == 0 means "auto" — show enough digits to distinguish).
///
/// ### Edge cases
/// - sign: prepended to the integer part string
/// - zero padding: fractional part is left-padded with '0' to reach precision
/// - trimming: trailing '0's are removed when precision == 0
/// - whole numbers: no decimal point is appended when fraction is empty
///
/// @param width     Minimum width of the integer part (passed to integer::to_string).
/// @param precision Number of decimal digits (0 = automatic, trim trailing zeros).
/// @return Formatted decimal string.
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Stream output operator.
///
/// Supports three output modes controlled by `std::ios_base` flags:
/// - `hex` flag: Output as hexadecimal fraction: `[-]x[ numerator_abs ]/x[ denominator_ ]`
/// - `dec` flag: Output as decimal string via to_string(), respecting
///   stream width() and precision() settings.
///
/// @param out Output stream.
/// @param v   Numeric value to output.
/// @return Reference to the output stream.
// ---------------------------------------------------------------------------
inline std::ostream & operator << (std::ostream & out, const nn::numeric & v)
{
	if( out.flags() & std::ios_base::hex ){
		out << (v.numerator_.is_neg() ? "-" : "")
			<< "x" << std::hex << v.numerator_.abs()
			<< "/" << "x" << std::hex << v.denominator_;
	}
	else if( out.flags() & std::ios_base::dec ){
		out << v.to_string(uintptr_t(out.width()), uintptr_t(out.precision()));
	}

	return out;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief BBP (Bailey–Borwein–Plouffe) π computation state.
///
/// Encapsulates the running state of the BBP summation. The formula is:
///
///   π = Σ_{k=0}^{∞} (1/16^k) * [4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6)]
///
/// Each term contributes approximately log10(16) ≈ 1.204 decimal digits.
/// The terms converge quickly (geometric at 1/16^k) so no denominator
/// blowup occurs.
///
/// ### Advantages of BBP for arbitrary-precision
/// Unlike Machin's formula or Ramanujan's series, BBP can compute
/// individual hex digits of π without computing preceding digits (though
/// this implementation does not exploit the digit-extraction property).
///
/// ### References
/// - Bailey, D.H., Borwein, P.B., Plouffe, S. "On the Rapid Computation
///   of Various Polylogarithmic Constants", Mathematics of Computation,
///   Vol. 66, No. 218, pp. 903–913, 1997.
///
/// ### Alternative approaches (commented out)
/// The pi() function also contains commented-out implementations of:
/// - **Machin's formula**: π/4 = 4·arctan(1/5) - arctan(1/239)
///   (requires atan() with high precision — rational atan is expensive)
/// - **Ramanujan's series**: 1/π = (2√2/9801) Σ (4k)!(1103+26390k) / (k!^4 · 396^(4k))
///   (faster convergence but more complex terms)
/// - **Continued fractions**: cfpi() — simple but slow convergence
// ---------------------------------------------------------------------------
struct pi_cont {
	numeric e, four, five, six, one, two;
	integer sixteen, eight;
	uintptr_t i;

	pi_cont();

	pi_cont & to_iter(uintptr_t iter);
};
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Continued fraction for π (alternative method, documented only).
///
/// Computes an approximation to π using the generalized continued fraction:
///   π = 4 / cfpi(1, 1, 0, iter)
///
/// where `cfpi(q, v, i, iter)` evaluates one level of:
///   cfpi(q, v, i, iter) = q + v^2 / cfpi(q+2, v+1, i+1, iter)
///
/// This is a straightforward recursive continued fraction with terminal
/// depth `iter`. Convergence is slower than BBP — this is provided as
/// a documented alternative for educational purposes.
///
/// @param q    Current linear term (q, q+2, q+4, ...).
/// @param v    Current value (v, v+1, v+2, ...).
/// @param i    Current iteration index.
/// @param iter Maximum depth of the continued fraction.
/// @return Approximation of the continued fraction value.
// ---------------------------------------------------------------------------
static inline numeric cfpi(const integer & q, const integer & v, uintptr_t i, uintptr_t iter)
{
	if( i >= iter )
		return 1;
	return numeric(q) + numeric(v.pow(2)) / cfpi(q + 2, v + 1, i + 1, iter);
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Initialize the BBP accumulator with the first term (k = 0).
///
/// The initial term (k=0) is:
///   e = 4 - 2/4 - 1/5 - 1/6
///     = 4 - 0.5 - 0.2 - 0.166666...
///     = 4 - 0.866666...
///     = 3.133333...
///
/// Subsequent terms are added by to_iter().
// ---------------------------------------------------------------------------
inline pi_cont::pi_cont() :
	e(0), four(4), five(5), six(6), one(1), two(2), sixteen(16), eight(8), i(1)
{
	e = four - two / four - one / five - one / six;
}
//------------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/// @brief Advance the BBP summation by one or more terms.
///
/// Each call to to_iter(iter) adds terms k = i .. iter-1 of the BBP series:
///
///   π = Σ_{k=0}^{∞} (1/16^k) * [4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6)]
///
/// The `sixteen` field tracks 16^i — it is left-shifted by 4 bits each
/// iteration, effectively multiplying by 16.
///
/// @param iter Target number of terms (k = 0..iter-1). Must be ≥ current i.
/// @return Reference to this pi_cont (for chaining).
// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------
/// @brief Compute π using the BBP formula.
///
/// Computes an approximation of π using the Bailey–Borwein–Plouffe series:
///
///   π = Σ_{k=0}^{∞} (1/16^k) * [4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6)]
///
/// Each iteration contributes ~1.2 decimal digits of precision. So `iter=20`
/// yields approximately 24 correct decimal digits.
///
/// ### State persistence
/// If `p_cnt` is non-null, the BBP state is stored externally, allowing
/// incremental computation across multiple calls (e.g., to add more terms
/// later without restarting from k=0). If null, a temporary accumulator
/// is used and discarded.
///
/// ### Alternative approaches (commented out in source)
/// The function body contains commented-out implementations of Machin's
/// formula and Ramanujan's series, preserved for reference:
///
/// - **Machin** (1706): π/4 = 4·arctan(1/5) - arctan(1/239)
///   Requires atan() on rationals 1/5 and 1/239 → expensive for high precision.
///
/// - **Ramanujan** (1914): 1/π = (2√2/9801) Σ (4k)!(1103+26390k) / (k!^4 · 396^(4k))
///   Each term adds ~8 decimal digits — much faster than BBP — but the terms
///   involve large factorials and high powers.
///
/// - **Continued fraction**: cfpi() — simplest but slowest.
///
/// @param iter   Number of BBP terms to sum (k = 0..iter-1).
/// @param p_cnt  Optional persistent state pointer; if non-null, the state
///               is updated incrementally. If null, a fresh accumulator is used.
/// @return Approximation of π as a rational.
// ---------------------------------------------------------------------------
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
	if( p_cnt == nullptr ){
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
