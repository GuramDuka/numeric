#include "nn.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <string>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name, expr) do {                                   \
    if (!(expr)) {                                              \
        fprintf(stderr, "FAIL: %s (%s)\n", name, #expr);        \
        tests_failed++;                                         \
    } else {                                                    \
        tests_passed++;                                         \
    }                                                           \
} while(0)

#define TEST_CLOSE(name, a, b, eps_str) do {                    \
    auto diff = ((a) - (b)).abs();                              \
    nn::numeric eps(eps_str);                                   \
    if (!(diff < eps)) {                                        \
        fprintf(stderr, "FAIL: %s | diff=%s\n",                \
                name, diff.to_string(0, 20).c_str());          \
        tests_failed++;                                         \
    } else {                                                    \
        tests_passed++;                                         \
    }                                                           \
} while(0)

// ── Numeric Tests ─────────────────────────────────────────────────────────────

static void test_numeric_construct() {
    nn::numeric a;
    TEST("default constructor is zero", a.is_zero());

    nn::numeric b(42);
    TEST("int constructor", b.to_string() == "42");

    nn::numeric c(3.14159);
    TEST("double constructor non-empty", !c.is_zero());

    nn::numeric d(std::string("123.456"));
    TEST("string constructor", d.to_string() == "123.456");

    nn::numeric e(d);
    TEST("copy constructor", e.to_string() == "123.456");

    nn::numeric f(nn::integer(10), nn::integer(3));
    TEST("numerator/denominator constructor", f.to_string(0, 6) == "3.333333");
}

static void test_numeric_arithmetic() {
    nn::numeric a(10);
    nn::numeric b(3);

    TEST("add", (a + b).to_string() == "13");
    TEST("sub", (a - b).to_string() == "7");
    TEST("mul", (a * b).to_string() == "30");

    // Division: 10/3 = 3.333...
    nn::numeric q = a / b;
    TEST("div non-integer", q.to_string(0, 4) == "3.3333");

    // Large rational
    nn::numeric big("12345678901234567890.123456789");
    TEST("big string construct", !big.is_zero());
}

static void test_numeric_comparison() {
    nn::numeric a(5);
    nn::numeric b(10);
    nn::numeric c(5);

    TEST("less", a < b);
    TEST("greater", b > a);
    TEST("equal", a == c);
    TEST("not equal", a != b);
}

static void test_numeric_abs_neg() {
    nn::numeric a(-42.5);
    TEST("abs", a.abs().to_string() == "42.5");
    TEST("negate", (-a).to_string() == "42.5");
    TEST("double negate", (-(-a)).to_string() == "-42.5");
}

static void test_numeric_pow() {
    nn::numeric a(2);
    TEST("pow(2, 2)", a.pow(2).to_string() == "4");
    TEST("pow(2, 3)", a.pow(3).to_string() == "8");
    TEST("pow(2, 10)", a.pow(10).to_string() == "1024");
    TEST("pow(2, 0)", a.pow(0).to_string() == "1");

    nn::numeric b(5);
    nn::numeric pow5 = b.pow(-1);
    pow5.normalize(3);
    std::string sp = pow5.to_string(0, 6);
    TEST("pow(5, -1) starts with 0.2", sp.substr(0, 3) == "0.2");
}

static void test_numeric_sqrt() {
    // sqrt(4) = 2  (converges in ~5 Newton iterations)
    nn::numeric four(4);
    nn::numeric r = four.root(2, 6);
    r.normalize(3);
    TEST("sqrt(4)", r.to_string() == "2");

    // sqrt(2) ~ 1.414213562373095...
    nn::numeric two(2);
    nn::numeric sqrt2 = two.root(2, 8);
    sqrt2.normalize(3);
    nn::numeric approx("1.414213562373095");
    nn::numeric diff = (sqrt2 - approx).abs();
    TEST("sqrt(2) close", diff < nn::numeric("0.0001"));
}

static void test_numeric_trig() {
    // sin(0) = 0
    nn::numeric zero(0);
    TEST("sin(0)", zero.sin(10).to_string() == "0");

    // cos(0) = 1
    TEST("cos(0)", zero.cos(10).to_string() == "1");
}

static void test_numeric_round() {
    nn::numeric a("3.14159");

    nn::numeric r1 = a.round(2, nn::numeric::toLess);
    TEST("round toLess 2 digits", r1.to_string(0, 4) == "3.1400");

    nn::numeric r2 = a.round(2, nn::numeric::toGreater);
    TEST("round toGreater 2 digits", r2.to_string(0, 4) == "3.1500");

    // Negative round tests
    nn::numeric neg("-3.14159");
    nn::numeric rn1 = neg.round(2, nn::numeric::toLess);
    TEST("neg round toLess 2 digits", rn1.to_string(0, 4) == "-3.1500");

    nn::numeric rn2 = neg.round(2, nn::numeric::toGreater);
    TEST("neg round toGreater 2 digits", rn2.to_string(0, 4) == "-3.1400");
}

static void test_numeric_to_string() {
    nn::numeric a("123.456789");
    std::string s = a.to_string(0, 3);
    TEST("to_string precision 3", s == "123.456");
}

static void test_numeric_edge_cases() {
    nn::numeric z(0);
    nn::numeric one(1);

    TEST("zero + 1 = 1", (z + one).to_string() == "1");
    TEST("1 - 1 = 0", (one - one).is_zero());
    TEST("0 * anything = 0", (z * one).is_zero());

    // Division by zero should throw
    bool caught = false;
    try {
        nn::numeric(1) / nn::numeric(0);
    } catch (const std::range_error &) {
        caught = true;
    }
    TEST("division by zero throws", caught);
}

static void test_numeric_big_pi() {
    // Quick test: just 10 iterations of pi formula
    nn::numeric pi_val = nn::pi(10, nullptr);
    pi_val.normalize(3);
    TEST("pi(10) is positive", pi_val > nn::numeric(3));
    TEST("pi(10) is less than 4", pi_val < nn::numeric(4));
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    test_numeric_construct();
    test_numeric_arithmetic();
    test_numeric_comparison();
    test_numeric_abs_neg();
    test_numeric_pow();
    test_numeric_sqrt();
    test_numeric_trig();
    test_numeric_round();
    test_numeric_to_string();
    test_numeric_edge_cases();
    test_numeric_big_pi();

    printf("\nResults: %d passed, %d failed out of %d\n",
           tests_passed, tests_failed, tests_passed + tests_failed);

    return tests_failed > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
