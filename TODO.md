# TODO

Analysis findings for the `xpolynomial` project (header-only C++23 polynomial library).
Items are grouped by priority. Findings marked `[verified]` were re-checked by hand.

## High priority (correctness)

- [x] **H1: `legendre_polynomial::create(order)` returns wrong polynomials for `order >= 10`**
  - `include/legendre_polynomial.tpp:48-59`
  - Recurrence uses `n = order` instead of the loop variable `i` (line 53), and
    `P_nm2 = P_nm1` (line 59) copies the wrong value instead of shifting `P_{n-1}/P_{n-2}`.
  - Fix the loop (`n = i`, proper shift `P_nm2 = P_nm1; P_nm1 = P_n;`) and add a
    regression test for `create(9..15)`. `[verified]`

- [x] **H2: `real_polynomial_root_finder::find_roots()` throws instead of returning empty lists**
  - `include/real_polynomial_root_finder.tpp:297`
  - `yun_algorithm(p).value()` throws `std::bad_optional_access` when the coefficients are
    not all integers (`square_free_decomposition.tpp:25`), contradicting the documented
    behavior in `include/real_polynomial_root_finder.h:119-120`.
  - Return `{ {}, {} }` for non-integer input; add a test. `[verified]`

- [x] **H3: `polynomial::divide()`/`operator/` has no zero-divisor guard**
  - `include/polynomial.tpp:465-479`
  - Dividing by the zero polynomial leads to NaN/Inf. The parser guards against it
    (`polynomial_parser.h:425-428`), the public API does not.
  - Document the invariant or guard with an error/assert; add a test.

- [x] **H4: Headers are not self-contained**
  - `<optional>` used but only included transitively (`root_finder.h:13`) where directly used in
    `square_free_decomposition.h:40,45,51`, `real_polynomial_root_finder.h:34,59,74,80,102,107`,
    `polynomial_interpolation.h:36`.
  - `<tuple>` missing in `polynomial.h:325`, `euclidean_algorithm.h:45`.
  - `<type_traits>` missing in `chebyshev_polynomial.h:29`, `polynomial_interpolation.h:25`.
  - Add direct includes for every type used; add an "include single header alone" compile test.

- [x] **H5: `polynomial_parser.h` relies on a foreign `using std::views::transform`**
  - `include/polynomial_parser.h:321,535,557`
  - `transform` is only in scope via `using std::views::transform;` from
    `include/polynomial.tpp:18`. Coupled to the umbrella/include order.
  - Introduce the `using`-declaration inside `polynomial_parser.h`. `[verified]`

## Medium priority (design / maintenance / docs)

- [x] **M1: Missing `constexpr` despite C++23 header-only positioning**
  - Only static constants are `constexpr`. Consider `constexpr` for `polynomial` ctors,
    `evaluate`, `degree`, operators, `utils::nearly_*`, `interval` methods.

- [x] **M2: Missing `noexcept` in a non-throwing library**
  - `noexcept` only on move ctor/assignment (`polynomial.h:77,179`) and parser lambdas.
  - Add `noexcept` to pure accessors (`is_*`, `degree()`, `leading_coefficient()`, `coefficients()`).

- [x] **M3: Copy-paste errors in doc comments of `polynomial.h`**
  - `operator*=` documented as "Divides the polynomial with a scalar" (`:250-251`).
  - `operator/` `@param scalar ... to be multiplied with` (`:254-255`).
  - `operator-=` documented as "Subtraction assignment operator (+=)" (`:284-286`).
  - `derive()` instead of `derivative` throughout (`polynomial.h:337-339`,
    `root_finder.h:49,58`, `real_polynomial_root_finder.h:147`).

- [x] **M4: Inconsistent include-guard style**
  - `#pragma once` in `interval.h:9`, `polynomial_parser.h:9`, all `.tpp`; include guards
    everywhere else. Unify (recommended: include guards).

- [x] **M5: `interval<T>` privately inherits `std::pair<T,T>`**
  - `include/interval.h:48`; consider composition (`lower_`, `upper_` fields) instead.
  - Inconsistent naming `opened` vs `is_lower_open()`; German typo "ist empty" (`:95`);
    wrong `@param upper Now the lower boundary` (`:74-75`).

- [x] **M6: No `static_assert` guard for unsupported coefficient types**
  - Primary `polynomial_specification<T>` is empty (`polynomial.h:23-25`) so `polynomial<int>`
    fails with cryptic errors; `complex_polynomial.h:24` shows the pattern.

- [x] **M7: Non-const `at()`/`operator[]` out-of-bounds write (UB)**
  - `include/polynomial.tpp:168-170`: `return coeffs_[index];` without bounds check, e.g.
    `p[10] = 5`. Const overload returns `spec::zero` for out-of-range. Fix or document. `[verified]`

- [x] **M8: `std::function` overhead in numerics**
  - `root_finder.h:30-45,56-61`, `chebyshev_polynomial.h:65` take `std::function`; consider
    templated callables.

- [x] **M9: Test CMake setup**
  - `file(GLOB TESTS ...)` (`test/CMakeLists.txt:13`) is fragile.
  - `GTest::gtest_main` and `GTest::gmock_main` linked although `test_main.cpp` provides `main`.
  - `enable_testing()` without `add_test(...)` -> `ctest` runs nothing.

- [x] **M10: Root CMake lacks build options**
  - `add_subdirectory(examples)` + `add_subdirectory(test)` unconditional
    (`CMakeLists.txt:54-55`) -> consumers without GoogleTest cannot configure.
  - Consider `BUILD_TESTING` / examples option; `CMAKE_CXX_EXTENSIONS ON` (`:16`);
    `set(PROJECT_NAME ...)` redundant (`:18`).

- [x] **M11: Namespace-local `static` objects in parser header**
  - `polynomial_parser.h:67,109,112-130`: one copy per TU, dynamic init.
  - Duplicated `SIGN_MINUS` key in `precedence_map` (`:338-339`).

- [x] **M12: Wrong `@file` names (copy-paste)**
  - `complex_polynomial.h:1`, `polynomial.ostream.tpp:1`, `complex_polynomial.tpp:1`,
    `examples/root_finder.cpp:1`, `examples/cubic_roots.cpp:1`, `examples/parser.cpp:1`,
    `test/complex_polynomial_tests.cpp:1`.

## Low priority (cosmetic / CI / README / tests)

- [x] **N1: No CI** - `.github/workflows` missing; add GCC-14/Clang-17 build + test job.

- [x] **N2: README gaps** - no sections for `interval`, `legendre_polynomial`,
  `complex_polynomial`/`separate`, `complex_polynomial_root_finder`, `polynomial_parser`;
  no error-handling notes; no API overview/changelog.

- [x] **N3: Stale `@date`/copyright years** - `LICENSE:3` says 2024, headers 2026;
  `@date` mostly 2023/2024.

- [x] **N4: Test quality**
  - Sanity-only checks for complex root finders (no golden values / tolerance).
  - Duplicated assertions: `real_polynomial_root_finder_tests.cpp:98-107` (three times
    `roots[0]`), `polynomial_tests.cpp:445-459` (identical DivideTest).
  - Typos in test names: `ZeroPloynomialTest`, `OnePloynomialTest`,
    `SubstrationWithScalarTest`, `PolynomialSubstractionTest`,
    `QuadraticPolynomialWhitoutRealRootsTest`.

- [x] **N5: Example gaps** - no examples for `interval`, `legendre_polynomial`,
  `complex_polynomial`, parser error handling, non-integer `find_roots`.

- [x] **N6: Minor issues**
  - `euclidean_algorithm.h:38` doc example uses `EuclideanAlgorithm` not `euclidean_algorithm`.
  - `square_free_decomposition.h:2` "for for" typo; `:40,45` "optional has not a value".
  - Missing `@brief` in `chebyshev_polynomial.h:36,48,53`, `polynomial_interpolation.h:27,32`.
  - `polynomial_interpolation.tpp:24` truncates `size_t` to `uint16_t`.
  - `utils.h:93-102` `is_even`/`is_odd` only for `long`.
  - `.gitignore` still contains `/vcpkg_installed/` (dead entry).

## Round 2 findings (re-analysis, hand-verified)

A second full pass over the code base. Items below were re-checked by reading the
actual source. Findings that turned out to be already fixed (e.g. self-contained
headers, `lagrange` bounds) were dropped.

### High priority (correctness)

- [x] **H6: `chebyshev_nodes(N = 0)` falls through / divides by zero**
  - `include/chebyshev_polynomial.tpp:80-82`
  - The `N == 0` branch constructs a temporary but never `return`s it, so execution
    continues to `(2. * k - 1.) / N` (`:86`) and the node transform with `N == 0`.
  - Add `return {};` (or `return xmath::chebyshev_polynomial<T>::values_type{};`). `[verified]`

- [x] **H7: `clenshaw()` / `chebyshev_series()` read `alphas[0]` on empty input (OOB)**
  - `include/chebyshev_polynomial.tpp:107-113`, `:123-129`
  - The loop is skipped for an empty vector but the final
    `return alphas[0] + x * beta1 - beta2;` still indexes element 0.
  - Guard against `alphas.empty()` and return zero. `[verified]`

- [x] **H8: Parser pops an empty operator stack on an unmatched `)` (UB)**
  - `include/polynomial_parser.h:384-389`
  - In `process_parenthesis` the `CLOSED` case calls `operator_stack.pop()`
    unconditionally. For input like `"1 + )"` the stack is empty (or holds only
    operators) and `pop()` / `top()` on an empty `std::stack` is undefined.
  - Validate `)` against a matching `(` (emit a parse error) before popping. `[verified]`

- [x] **H9: Cardano formula uses `std::pow(negative, 1./3.)` -> NaN**
  - `include/real_polynomial_root_finder.tpp:123-124`
  - `A`/`B` are computed with `std::pow(v, 1. / 3.)`. For `v < 0` this yields NaN
    even though the real cube root is well defined (casus irreducibilis aside).
  - Use `std::cbrt(v)` (or sign-aware `std::pow(std::abs(v), 1./3.)`). `[verified]`

- [x] **H10: FetchContent consumers without GoogleTest fail to configure**
  - `CMakeLists.txt:19,59-61`; `README.md` integration section
  - `include(CTest)` makes `BUILD_TESTING` default to `ON`; `add_subdirectory(test)`
    then runs `find_package(GTest REQUIRED)`. A consumer using `FetchContent` (as the
    README advertises) will abort at configure time unless it also has GoogleTest.
  - Guard development-only targets with `if(PROJECT_IS_TOP_LEVEL)` (and keep
    `BUILD_TESTING`/examples usable when built standalone). `[verified]`

### Medium priority (design / maintenance / docs)

- [x] **M13: `operator*=(scalar)` / `operator/=(scalar)` skip `trim_coefficients()`**
  - `include/polynomial.tpp:247-252,262-267`
  - `operator*=` / `operator/=` mutate each coefficient in place but never trim, so
    multiplying by `0` leaves the degree unchanged while all coefficients are zero,
    breaking the "degree == normalized size" invariant (cf. `operator==`).
  - Call `trim_coefficients()` before returning. `[verified]`

- [x] **M14: `normalize()` / scalar division lack a zero-leading-coefficient guard**
  - `include/polynomial.tpp:255-267,451-453`
  - `normalize()` divides by `leading_coefficient()`; for the zero polynomial this is
    `0` -> NaN/Inf. Document the precondition (non-zero polynomial) or return `*this`. `[verified]`

- [x] **M15: `is_linear()` returns `degree() <= 1`**
  - `include/polynomial.tpp:118-120`
  - Includes constants (degree 0) and the zero polynomial, contradicting the doc
    ("degree 1"); sibling `is_quadratic`/`is_cubic` use `== 2` / `== 3`.
  - Use `degree() == 1`. `[verified]`

- [x] **M16: Iterative root finders lack a maximum-iteration guard**
  - `include/root_finder.tpp:28-43` (bisection), `:57-72` (regula falsi)
  - `bisection`/`regula_falsi` loop only on interval width; for numeric edge cases
    (denominator `func(b) - func(a)` near zero, stagnation) they can diverge or spin.
  - Add a `max_iterations` parameter (like `newton_raphson`) and guard the
    regula-falsi denominator. `[verified]`

- [x] **M17: Aberth-Ehrlich denominator can be zero**
  - `include/complex_polynomial_root_finder.tpp:100`
  - `p_prim(z) - p_norm(z) * S(z)` is unguarded; for multiple/clustered roots it can
    vanish -> division by zero (Inf/NaN propagated into the fixed-point iteration).
  - Guard/skip or fall back to a different step when the denominator is ~0. `[verified]`

- [ ] **M18: Parser treats `^` as left-associative**
  - `include/polynomial_parser.h:355-362`
  - `top_precedence_greater_or_equal` uses `>=` for every operator, so `"2^3^2"` parses
    as `(2^3)^2 = 64` instead of the conventional right-associative `2^(3^2) = 512`.
  - Special-case right-associativity for `POWER` (`>` instead of `>=`). `[verified]`

- [ ] **M19: CI lacks hardening (warnings-as-errors, sanitizers, coverage, MSVC)**
  - `.github/workflows/ci.yml`
  - Builds GCC 14 / Clang 18 with `-Wall -Wextra -Wpedantic` but no `-Werror`, no
    ASan/UBSan job, no coverage, and no MSVC coverage despite the README claim.
  - Add a `-Werror` build, a sanitizer job, and (optionally) MSVC. `[verified]`

### Low priority (cosmetic / docs)

- [ ] **N7: Residual N6 fixes were never applied**
  - `include/square_free_decomposition.tpp:2` still says "for for" (only the `.h` was fixed).
  - `include/root_finder.h:44,61` still say "the returned optional<> has not a value"
    (N6 only rewrote the copies in `square_free_decomposition.h`). `[verified]`

- [ ] **N8: `constexpr` (M1) only landed for helpers/`interval`**
  - `5767068` made numeric helpers and `interval` `constexpr`, but `polynomial` ctors,
    `evaluate`, `degree` and operators are still non-`constexpr` although the TODO
    item listed them. Either extend or narrow the M1 description. `[verified]`