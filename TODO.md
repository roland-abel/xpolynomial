# TODO

Analysis findings for the `xpolynomial` project (header-only C++23 polynomial library).
Items are grouped by priority. Findings marked `[verified]` were re-checked by hand.

## High priority (correctness)

- [ ] **H1: `legendre_polynomial::create(order)` returns wrong polynomials for `order >= 10`**
  - `include/legendre_polynomial.tpp:48-59`
  - Recurrence uses `n = order` instead of the loop variable `i` (line 53), and
    `P_nm2 = P_nm1` (line 59) copies the wrong value instead of shifting `P_{n-1}/P_{n-2}`.
  - Fix the loop (`n = i`, proper shift `P_nm2 = P_nm1; P_nm1 = P_n;`) and add a
    regression test for `create(9..15)`. `[verified]`

- [ ] **H2: `real_polynomial_root_finder::find_roots()` throws instead of returning empty lists**
  - `include/real_polynomial_root_finder.tpp:297`
  - `yun_algorithm(p).value()` throws `std::bad_optional_access` when the coefficients are
    not all integers (`square_free_decomposition.tpp:25`), contradicting the documented
    behavior in `include/real_polynomial_root_finder.h:119-120`.
  - Return `{ {}, {} }` for non-integer input; add a test. `[verified]`

- [ ] **H3: `polynomial::divide()`/`operator/` has no zero-divisor guard**
  - `include/polynomial.tpp:465-479`
  - Dividing by the zero polynomial leads to NaN/Inf. The parser guards against it
    (`polynomial_parser.h:425-428`), the public API does not.
  - Document the invariant or guard with an error/assert; add a test.

- [ ] **H4: Headers are not self-contained**
  - `<optional>` used but only included transitively (`root_finder.h:13`) where directly used in
    `square_free_decomposition.h:40,45,51`, `real_polynomial_root_finder.h:34,59,74,80,102,107`,
    `polynomial_interpolation.h:36`.
  - `<tuple>` missing in `polynomial.h:325`, `euclidean_algorithm.h:45`.
  - `<type_traits>` missing in `chebyshev_polynomial.h:29`, `polynomial_interpolation.h:25`.
  - Add direct includes for every type used; add an "include single header alone" compile test.

- [ ] **H5: `polynomial_parser.h` relies on a foreign `using std::views::transform`**
  - `include/polynomial_parser.h:321,535,557`
  - `transform` is only in scope via `using std::views::transform;` from
    `include/polynomial.tpp:18`. Coupled to the umbrella/include order.
  - Introduce the `using`-declaration inside `polynomial_parser.h`. `[verified]`

## Medium priority (design / maintenance / docs)

- [ ] **M1: Missing `constexpr` despite C++23 header-only positioning**
  - Only static constants are `constexpr`. Consider `constexpr` for `polynomial` ctors,
    `evaluate`, `degree`, operators, `utils::nearly_*`, `interval` methods.

- [ ] **M2: Missing `noexcept` in a non-throwing library**
  - `noexcept` only on move ctor/assignment (`polynomial.h:77,179`) and parser lambdas.
  - Add `noexcept` to pure accessors (`is_*`, `degree()`, `leading_coefficient()`, `coefficients()`).

- [ ] **M3: Copy-paste errors in doc comments of `polynomial.h`**
  - `operator*=` documented as "Divides the polynomial with a scalar" (`:250-251`).
  - `operator/` `@param scalar ... to be multiplied with` (`:254-255`).
  - `operator-=` documented as "Subtraction assignment operator (+=)" (`:284-286`).
  - `derive()` instead of `derivative` throughout (`polynomial.h:337-339`,
    `root_finder.h:49,58`, `real_polynomial_root_finder.h:147`).

- [ ] **M4: Inconsistent include-guard style**
  - `#pragma once` in `interval.h:9`, `polynomial_parser.h:9`, all `.tpp`; include guards
    everywhere else. Unify (recommended: include guards).

- [ ] **M5: `interval<T>` privately inherits `std::pair<T,T>`**
  - `include/interval.h:48`; consider composition (`lower_`, `upper_` fields) instead.
  - Inconsistent naming `opened` vs `is_lower_open()`; German typo "ist empty" (`:95`);
    wrong `@param upper Now the lower boundary` (`:74-75`).

- [ ] **M6: No `static_assert` guard for unsupported coefficient types**
  - Primary `polynomial_specification<T>` is empty (`polynomial.h:23-25`) so `polynomial<int>`
    fails with cryptic errors; `complex_polynomial.h:24` shows the pattern.

- [ ] **M7: Non-const `at()`/`operator[]` out-of-bounds write (UB)**
  - `include/polynomial.tpp:168-170`: `return coeffs_[index];` without bounds check, e.g.
    `p[10] = 5`. Const overload returns `spec::zero` for out-of-range. Fix or document. `[verified]`

- [ ] **M8: `std::function` overhead in numerics**
  - `root_finder.h:30-45,56-61`, `chebyshev_polynomial.h:65` take `std::function`; consider
    templated callables.

- [ ] **M9: Test CMake setup**
  - `file(GLOB TESTS ...)` (`test/CMakeLists.txt:13`) is fragile.
  - `GTest::gtest_main` and `GTest::gmock_main` linked although `test_main.cpp` provides `main`.
  - `enable_testing()` without `add_test(...)` -> `ctest` runs nothing.

- [ ] **M10: Root CMake lacks build options**
  - `add_subdirectory(examples)` + `add_subdirectory(test)` unconditional
    (`CMakeLists.txt:54-55`) -> consumers without GoogleTest cannot configure.
  - Consider `BUILD_TESTING` / examples option; `CMAKE_CXX_EXTENSIONS ON` (`:16`);
    `set(PROJECT_NAME ...)` redundant (`:18`).

- [ ] **M11: Namespace-local `static` objects in parser header**
  - `polynomial_parser.h:67,109,112-130`: one copy per TU, dynamic init.
  - Duplicated `SIGN_MINUS` key in `precedence_map` (`:338-339`).

- [ ] **M12: Wrong `@file` names (copy-paste)**
  - `complex_polynomial.h:1`, `polynomial.ostream.tpp:1`, `complex_polynomial.tpp:1`,
    `examples/root_finder.cpp:1`, `examples/cubic_roots.cpp:1`, `examples/parser.cpp:1`,
    `test/complex_polynomial_tests.cpp:1`.

## Low priority (cosmetic / CI / README / tests)

- [ ] **N1: No CI** - `.github/workflows` missing; add GCC-14/Clang-17 build + test job.

- [ ] **N2: README gaps** - no sections for `interval`, `legendre_polynomial`,
  `complex_polynomial`/`separate`, `complex_polynomial_root_finder`, `polynomial_parser`;
  no error-handling notes; no API overview/changelog.

- [ ] **N3: Stale `@date`/copyright years** - `LICENSE:3` says 2024, headers 2026;
  `@date` mostly 2023/2024.

- [ ] **N4: Test quality**
  - Sanity-only checks for complex root finders (no golden values / tolerance).
  - Duplicated assertions: `real_polynomial_root_finder_tests.cpp:98-107` (three times
    `roots[0]`), `polynomial_tests.cpp:445-459` (identical DivideTest).
  - Typos in test names: `ZeroPloynomialTest`, `OnePloynomialTest`,
    `SubstrationWithScalarTest`, `PolynomialSubstractionTest`,
    `QuadraticPolynomialWhitoutRealRootsTest`.

- [ ] **N5: Example gaps** - no examples for `interval`, `legendre_polynomial`,
  `complex_polynomial`, parser error handling, non-integer `find_roots`.

- [ ] **N6: Minor issues**
  - `euclidean_algorithm.h:38` doc example uses `EuclideanAlgorithm` not `euclidean_algorithm`.
  - `square_free_decomposition.h:2` "for for" typo; `:40,45` "optional has not a value".
  - Missing `@brief` in `chebyshev_polynomial.h:36,48,53`, `polynomial_interpolation.h:27,32`.
  - `polynomial_interpolation.tpp:24` truncates `size_t` to `uint16_t`.
  - `utils.h:93-102` `is_even`/`is_odd` only for `long`.
  - `.gitignore` still contains `/vcpkg_installed/` (dead entry).