# Numeric

C++ arbitrary-precision integer and floating-point arithmetic library.

> A C++ class that implements infinite-length integers and floating-point numbers using integer fractions (rational arithmetic). Numbers are represented as `numerator / denominator` pairs of arbitrary-precision integers.

## Features

- **Arbitrary-precision integers** — add, subtract, multiply, divide, modulo, bitwise ops, shift, GCD, LCM, pow, factorial
- **Arbitrary-precision floating point** — represented as exact rational fractions (`num/den`), with transcendental functions
- **Transcendental functions** — `sin`, `cos`, `atan`, `root` (Newton's method), `pow`
- **π computation** — Bailey-Borwein-Plouffe (BBP) formula
- **Custom allocator** — Two-Level Segregated Fit (TLSF) memory allocator for deterministic performance
- **Thread pool** — parallel multiplication for large integers
- **Header-only** — drop in and compile (plus CMake build system)

## Build

### Prerequisites

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.14+
- pthreads (on Linux/macOS)

### Quick start

```bash
git clone <repo>
cd numeric
cmake -B build -DNUMERIC_BUILD_TESTS=ON
cmake --build build
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `NUMERIC_BUILD_TESTS` | `ON` | Build unit tests |
| `NUMERIC_BUILD_MAIN` | `ON` | Build demo/benchmark executable |

## Run

### Demo / Benchmark

```bash
./build/numeric_demo
```

Computes and validates:
- √2 (square root via Newton's method)
- ∛2 (cube root via bisection)
- sin(1), cos(1) (via Taylor series)
- π to 265+ digits (BBP formula)

### Tests

```bash
ctest --test-dir build
```

Or run directly:

```bash
./build/tests/test_integer
./build/tests/test_numeric
```

## Usage

The library lives in the `nn` (NaturalNumbers) namespace and is entirely header-based.

```cpp
#include "nn.hpp"
#include <iostream>

int main() {
    using namespace nn;

    // Integer arithmetic
    integer a("12345678901234567890");
    integer b("98765432109876543210");
    std::cout << "a + b = " << a + b << std::endl;
    std::cout << "a * b = " << a * b << std::endl;
    std::cout << "a / b = " << a / b << std::endl;
    std::cout << "gcd   = " << a.gcd(b) << std::endl;
    std::cout << "10!   = " << integer::factorial(10) << std::endl;

    // Floating point via rational fractions
    numeric pi = nn::pi(100, nullptr);
    std::cout << "π ≈ " << pi.to_string(0, 50) << std::endl;

    numeric x("1.57079632679489661923");  // π/2
    std::cout << "sin(π/2) ≈ " << x.sin(20).to_string(0, 30) << std::endl;

    numeric two(2);
    std::cout << "√2 ≈ " << two.root(2, 10).to_string(0, 30) << std::endl;
}
```

### API Overview

**`integer` class:**
- Construction: from `int`, `long`, `unsigned long long`, `std::string`, `numeric`
- Arithmetic: `+`, `-`, `*`, `/`, `%`, `+=`, `-=`, `*=`, `/=`, `%=`
- Bitwise: `&`, `|`, `^`, `~`, `<<`, `>>`, `&=`, `|=`, `^=`, `<<=`, `>>=`
- Comparison: `==`, `!=`, `<`, `>`, `<=`, `>=`
- Utilities: `abs()`, `sign()`, `pow()`, `factorial()`, `gcd()`, `lcm()`
- String: `to_string(width, base)` — base 10, 8, or 16

**`numeric` class:**
- Construction: from `int`, `double`, `std::string`, `integer`, `(numerator, denominator)`
- Arithmetic: `+`, `-`, `*`, `/`, unary `-`, `++`, `--`
- Comparison: `==`, `!=`, `<`, `>`, `<=`, `>=`
- Functions: `abs()`, `pow(n)`, `root(power, iter)`, `sin(iter)`, `cos(iter)`, `atan(iter)`
- Rounding: `round(digits, toLess | toGreater)`
- String: `to_string(width, precision)`
- Free function: `pi(iter, &state)` — BBP π formula

## Bugs Fixed

- **Negative decimal string constructor** — `-3.14159` was incorrectly parsed as `-2.85841` (fraction was added instead of subtracted for negative numbers)
- **`round()` logic** — used `integer + numeric` which triggered implicit conversion of the numeric fraction to integer (truncating to 0). Also used the wrong variable for comparison in the rounding decision
- **Newton iteration blowup** — `root()` with large iteration counts caused exponential growth of rational numerator/denominator, leading to extremely slow GCD operations. Added early-exit convergence checks
- **Missing includes** — `wk.hpp` lacked `<cstring>` include for `memcpy`/`memset`, causing LSP errors

## Limitations

- **Newton iteration fraction blowup** — Newton's method with exact rational arithmetic causes numerator/denominator sizes to double each iteration. For `root()`, keep iterations ≤ 12, or use the bisection method (for powers ≠ 2). The function includes a `to_string`-based convergence guard to mitigate this
- **Thread pool** — parallel multiplication is always enabled when `hardware_concurrency > 1`

## License

MIT — see [LICENSE](LICENSE).

Copyright (c) 2014, 2015, 2016, 2017 Guram Duka
