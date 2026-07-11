# Numeric — Agent Guide

## What This Is

Header-only C++ arbitrary-precision arithmetic library. Numbers are `numerator/denominator` rational fractions backed by arbitrary-precision integers. Provides `sqrt`, `sin`, `cos`, `atan`, `pow`, π (BBP formula), GCD, LCM, factorial.

Namespace: `nn` (NaturalNumbers). Two main classes: `integer` and `numeric`.

## Build System

CMake-based. Three targets:

```bash
cmake -B build                     # configure
cmake --build build                # build all
ctest --test-dir build --timeout 30  # run tests
```

| Target | Purpose |
|--------|---------|
| `numeric_demo` | Benchmark: sqrt, sin, cos, pi computation |
| `test_integer` | Integer unit tests |
| `test_numeric` | Numeric unit tests |

Options: `-DNUMERIC_BUILD_TESTS=OFF` / `-DNUMERIC_BUILD_MAIN=OFF`

## Code Architecture (DAG)

```
tlsf.hpp  →  id.hpp  →  ii.hpp  →  nn.hpp  →  main.cpp
(id.hpp)  →  wk.hpp  →  ii.hpp
```

| File | Responsibility |
|------|---------------|
| `tlsf.hpp` | TLSF memory allocator (custom malloc/free) |
| `id.hpp` | `nn_integer_data` — word typedefs, ref-counted big-int core, bit ops, shift, compare |
| `wk.hpp` | Thread pool (`ThreadPoolT`), cryptographic sponge (cdc256/cdc512) |
| `ii.hpp` | `integer` class — arithmetic, bitwise, string conversion, GCD, pow, factorial |
| `nn.hpp` | `numeric` class — rational fraction arithmetic, sin/cos/atan/root/pow, π, rounding |

## Known Pitfalls (Read Before Editing)

### 1. Newton fraction blowup
`numeric::root(power=2, iter)` uses Newton's method on exact rationals. Each iteration **doubles** numerator/denominator bit-size. `iter=20` creates million-bit integers → GCD in `normalize(3)` grinds to a halt. **Keep iter ≤ 12** for Newton. Has `to_string` convergence guard but don't rely on it — prevent the blowup upstream.

Bisection method (power ≠ 2) does NOT have this problem.

### 2. String constructor negative bug
Was: `numeric("-3.14159")` gave `-2.85841`. Fixed: when `ipart.is_neg()`, fraction is **subtracted**, not added. If re-touching this code, verify both signs.

### 3. `integer + numeric` ambiguity
`integer + numeric` converts numeric → integer via `integer(numeric)` which does `numerator / denominator` (truncation). Always wrap as `numeric(ipart) + numeric(q, p)` to stay in rational domain. This broke `round()`.

### 4. Header-only globals
`nn_izero`, `nn_ione`, `nn_iten`, `nn_maxull`, and `integer::stat_iadd_` etc. are defined in headers. Compiling multiple TUs that include `id.hpp`/`ii.hpp` will cause **linker errors** (multiple definition). Each test/executable must be a separate CMake target — never link two TUs that both include these headers.

### 5. `normalize(3) = normalize(how=3)`
`how=3` means bit 1 (`how & 1`) → trim common trailing zero bits, and bit 2 (`how & 2`) → GCD reduce. `normalize(1)` skips GCD. Use `normalize(3)` after arithmetic, `normalize(2)` if you only want GCD reduction.

### 6. Thread pool
`wk.hpp` creates a global `ThreadPool` with `hardware_concurrency` threads. `integer::operator*` uses it for large multiplications. This is always-on when concurrency > 1.

## Development Conventions

- **All code is in headers** — no `.cpp` files except `main.cpp` and test files
- **C++17** required (for `std::index_sequence`, `std::shared_mutex`)
- **Add tests** in `tests/test_integer.cpp` or `tests/test_numeric.cpp` using the `TEST(name, expr)` macro
- **Test each new file separately** — never link multiple test TUs due to (4)
- **Edit existing files** — this is an existing project, prefer targeted edits over rewrites
- **CMake config** in root `CMakeLists.txt` only, tests in `tests/CMakeLists.txt`
