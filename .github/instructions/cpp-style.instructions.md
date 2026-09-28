---
description: "C++ coding style for the audi library (headers, tests, pyaudi bindings)."
applyTo: "**/*.{hpp,cpp,h,hpp.in}"
---

# C++ coding style

## Formatting

- Follow the repository `.clang-format`: 4-space indent, no tabs, 120-column limit, Linux brace style (opening brace on its own line for namespaces, classes and functions; same line for control statements).
- `template <...>` always on its own line, above the declaration.
- Pointer/reference bind to the name: `const T &x`, `T *p`.
- Break long expressions *before* binary operators (`|| ...`, `+ ...`).
- Namespace contents are not indented; close with `} // end of namespace audi` (or `} // namespace detail`).
- Short `if` bodies still use braces in practice; always use braces for loops.

## File layout

- Every file starts with the standard audi copyright/dual GPL3-LGPL3 license header (copy it from an existing file).
- Include guards of the form `AUDI_<FILENAME>_HPP` (no `#pragma once`).
- Include order: standard library headers, then third-party (`boost`, `obake`, `mp++`, `Eigen`), then `audi/` headers, each group alphabetically sorted and separated by a blank line. Use angle brackets for all, including `<audi/...>`.
- Annotate non-obvious includes with a short trailing comment (e.g. `// for audi::abs`).
- Everything lives in `namespace audi`; implementation helpers go in `namespace audi::detail` (written as nested `namespace detail`).

## Naming

- `snake_case` for types, functions, variables, and template aliases (`gdual`, `taylor_map`, `get_order`, `constant_cf`).
- Private data members prefixed with `m_` (`m_p`, `m_order`, `m_c`).
- Template parameters: short `CamelCase` or single letters (`Cf`, `Monomial`, `T`, `M`, `U`).
- Getters use `get_` prefix; predicates use `is_` prefix (`is_zero`, `is_vectorized`).
- Type aliases inside classes end in `_type` (`cf_type`, `key_type`, `p_type`, `v_type`).
- Low-level/internal accessors are prefixed with `_` (`_poly()`, `_container()`).
- Macros are `UPPER_CASE`; `#undef` header-local macros at the end of the header.

## Templates and generic code

- Header-only: free functions are `template` + `inline`.
- Constrain overloads with SFINAE via `enable_if_t<..., int> = 0` as a defaulted template parameter; name reusable enablers as private alias templates (`generic_ctor_enabler`, `gdual_if_enabled`, `operator_enabler`).
- Use `audi::is_arithmetic`, `is_vectorized`, and `obake::is_*_v` traits rather than raw `std::is_arithmetic` when mppp/vectorized types may be involved.
- Use `static_assert` in class templates to document coefficient-type requirements.
- Use `if constexpr` for compile-time branching inside function bodies.
- Specialise for `mppp::real128` inside `#if defined(AUDI_WITH_QUADMATH)` blocks.

## Classes

- Explicitly `= default` copy/move constructors and assignment operators.
- Mark single-argument and converting constructors `explicit`.
- Binary operators are `friend` functions inside the class, forwarding to private static helpers (`add`, `sub`, `mul`, `div`).
- Compound assignment operators are implemented in terms of the binary ones (`return *this = *this + d1;`), or vice versa for small value types (`vectorized`).
- Mark non-mutating methods `const`.
- Boost.Serialization: private `friend class boost::serialization::access;` and a templated `serialize(Archive &ar, const unsigned int)`.

## Idioms

- Prefer `auto` for locals initialised from expressions; use `decltype(x.size())` / `decltype(d.get_order())` for loop counters to match container index types.
- Use unsigned literals in loops/comparisons (`0u`, `1u`, `2u`).
- Write floating literals with a trailing dot (`1.`, `0.5`, `-1.`) and construct coefficients via `T(1.)` / `Cf(value)` to avoid precision loss (e.g. double vs real128).
- Call math functions qualified as `audi::exp`, `audi::abs`, etc. so the audi overloads are selected.
- Use `std::move` when returning/forwarding freshly built polynomials into constructors.
- Use STL algorithms (`std::transform`, `std::all_of`, `std::accumulate`) with lambdas for element-wise work.
- Use `static_cast` / `boost::numeric_cast`; avoid C-style casts.

## Errors

- Validate inputs and throw `std::invalid_argument` (or `std::domain_error`) with a descriptive message; include offending values via `std::to_string` when helpful.
- In code using `audi/exceptions.hpp`, prefer the `audi_throw(ExceptionType, "message")` macro.
- Use `assert` only for internal invariants (e.g. in destructors).

## Documentation

- Doxygen comments on public API: a `///` brief line followed by a `/** ... */` block.
- Use `@param`, `@return`, `@throws` (listing each exception with its condition), `\note`, and `@code ... @endcode` examples.
- Put the math in LaTeX with `\f[ ... \f]` / `\f$ ... \f$`, typically writing the expansion around \f$T_f = f_0 + \hat f\f$.
- Group related members with `/** @name ... */ //@{ ... //@}`.
- Inline `//` comments are short and explain *why*, not *what*.
