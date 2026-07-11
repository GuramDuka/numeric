#include "nn.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

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

// ── Integer Tests ─────────────────────────────────────────────────────────────

static void test_integer_construct() {
    nn::integer a;
    TEST("default constructor is zero", a.is_zero());

    nn::integer b(42);
    TEST("int constructor", b.to_string() == "42");

    nn::integer c(-123);
    TEST("negative int constructor", c.to_string() == "-123");

    nn::integer d(std::string("9876543210"));
    TEST("string constructor", d.to_string() == "9876543210");

    nn::integer e(d);
    TEST("copy constructor", e.to_string() == "9876543210");
}

static void test_integer_arithmetic() {
    nn::integer a(100);
    nn::integer b(3);
    nn::integer c(0);

    TEST("add", (a + b).to_string() == "103");
    TEST("sub", (a - b).to_string() == "97");
    TEST("mul", (a * b).to_string() == "300");
    TEST("div", (a / b).to_string() == "33");

    nn::integer mod;
    a.div(b, &mod);
    TEST("mod", mod.to_string() == "1");

    TEST("negate", (-a).to_string() == "-100");
    TEST("abs of neg", (-a).abs().to_string() == "100");

    // Large numbers
    nn::integer big1("12345678901234567890");
    nn::integer big2("98765432109876543210");
    TEST("big add", (big1 + big2).to_string() == "111111111011111111100");
    TEST("big sub", (big2 - big1).to_string() == "86419753208641975320");
}

static void test_integer_comparison() {
    nn::integer a(5);
    nn::integer b(10);
    nn::integer c(5);

    TEST("less", a < b);
    TEST("greater", b > a);
    TEST("equal", a == c);
    TEST("not equal", a != b);
    TEST("less or eq", a <= c);
    TEST("greater or eq", b >= a);
}

static void test_integer_bitops() {
    nn::integer a(0b1010);  // 10
    nn::integer b(0b1100);  // 12

    TEST("bitwise and", (a & b).to_string() == "8");
    TEST("bitwise or", (a | b).to_string() == "14");
    TEST("bitwise xor", (a ^ b).to_string() == "6");

    // Shift
    nn::integer c(1);
    TEST("shift left", (c << 10).to_string() == "1024");
    TEST("shift right", (nn::integer(1024) >> 10).to_string() == "1");
}

static void test_integer_pow_factorial() {
    TEST("pow(2, 10)", nn::integer::pow(10, 2).to_string() == "1024");
    TEST("pow(2, 0)", nn::integer::pow(0, 2).to_string() == "1");
    TEST("pow(2, 1)", nn::integer::pow(1, 2).to_string() == "2");

    TEST("factorial(0)", nn::integer::factorial(0).to_string() == "1");
    TEST("factorial(1)", nn::integer::factorial(1).to_string() == "1");
    TEST("factorial(5)", nn::integer::factorial(5).to_string() == "120");
    TEST("factorial(10)", nn::integer::factorial(10).to_string() == "3628800");
}

static void test_integer_gcd_lcm() {
    nn::integer a(12);
    nn::integer b(18);

    TEST("gcd", a.gcd(b).to_string() == "6");
    TEST("lcm", a.lcm(b).to_string() == "36");

    TEST("gcd with zero", a.gcd(nn::integer(0)).to_string() == "12");
    TEST("gcd of zero and zero", nn::integer(0).gcd(nn::integer(0)).to_string() == "0");
}

static void test_integer_inc_dec() {
    nn::integer a(5);
    ++a;
    TEST("pre-increment", a.to_string() == "6");
    a--;
    TEST("post-decrement (after pre-inc)", a.to_string() == "5");

    nn::integer b(-1);
    --b;
    TEST("pre-decrement negative", b.to_string() == "-2");
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    test_integer_construct();
    test_integer_arithmetic();
    test_integer_comparison();
    test_integer_bitops();
    test_integer_pow_factorial();
    test_integer_gcd_lcm();
    test_integer_inc_dec();

    printf("\nResults: %d passed, %d failed out of %d\n",
           tests_passed, tests_failed, tests_passed + tests_failed);

    return tests_failed > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
