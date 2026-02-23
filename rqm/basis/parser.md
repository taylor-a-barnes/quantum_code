# Feature: Parse QCSchema Basis Set File <!-- rq-f28bff90 -->

This feature provides functions for loading a Gaussian basis set into memory from a QCSchema JSON
file as produced by the BSE (see `rqm/basis/bse.md`). It exposes both a low-level file parser and
a high-level convenience function that fetches and parses in one call.

Exponents and contraction coefficients are converted to `f64` during parsing. SP shells (a single
shell entry with `angular_momentum: [0, 1]`) are automatically split into separate S and P shells,
each sharing the same exponents. ECP data, if present, is ignored; however if an element has no
electron shells at all, parsing fails.

Only single-element files are supported. The `elements` object in the JSON must contain exactly one
atomic-number key.

## Feature API <!-- rq-bb3e6797 -->

### Functions <!-- rq-a95c7ed4 -->

- `parse_basis(path: &Path) -> Result<BasisSet, ParseError>` <!-- rq-810be7c3 -->
  - Reads and parses a QCSchema basis set JSON file at `path`.
  - Validates that the `elements` object contains exactly one entry.
  - Derives the element symbol from the atomic-number key (1–118).
  - Splits any SP shell (multiple entries in `angular_momentum`) into one shell per angular momentum
    value, each inheriting the full set of exponents and its corresponding coefficient vector.
  - Converts all exponent and coefficient strings to `f64`; returns an error on any unparseable
    value.
  - Returns an error if the exponent count and the length of any coefficient vector differ.
  - Returns an error if `electron_shells` is absent or empty.

- `load_basis(element: &str, basis_name: &str) -> Result<BasisSet, LoadError>` <!-- rq-e279cb2c -->
  - Calls `fetch_basis(element, basis_name)` to obtain the cached file path (downloading if needed).
  - Calls `parse_basis` on the returned path.
  - Propagates `BseError` as `LoadError::Fetch` and `ParseError` as `LoadError::Parse`.

### Types <!-- rq-df930b4f -->

- `BasisSet` — the parsed representation of a single-element basis set: <!-- rq-a59faf69 -->
  - `element: String` — title-case element symbol (e.g. `"H"`, `"Au"`).
  - `atomic_number: u32` — atomic number (1–118).
  - `shells: Vec<ElectronShell>` — the contraction shells in file order, with SP shells already
    split.

- `ElectronShell` — a single contracted Gaussian shell: <!-- rq-5e83e69f -->
  - `angular_momentum: u32` — the angular momentum quantum number (0 = s, 1 = p, 2 = d, …).
  - `exponents: Vec<f64>` — primitive Gaussian exponents.
  - `coefficients: Vec<f64>` — contraction coefficients, one per exponent.

- `ParseError` — error type returned by `parse_basis`: <!-- rq-47577318 -->
  - `IoError(String)` — the file could not be read.
  - `InvalidJson(String)` — the file contents are not valid JSON.
  - `MultipleElements { found: usize }` — the `elements` object contained more than one key.
  - `NoElements` — the `elements` object is empty.
  - `InvalidAtomicNumber(String)` — the element key could not be parsed as a valid atomic number
    (1–118).
  - `NoElectronShells` — the element has no `electron_shells` key, or the array is empty.
  - `MalformedShell { index: usize, reason: String }` — a shell entry is structurally invalid
    (empty `angular_momentum`, mismatched `coefficients` count, mismatched coefficient-vector
    length, or an unparseable exponent or coefficient string).

- `LoadError` — error type returned by `load_basis`: <!-- rq-987b9f09 -->
  - `Fetch(BseError)` — the underlying `fetch_basis` call failed.
  - `Parse(ParseError)` — the file was fetched but could not be parsed.

## QCSchema Shell Format <!-- rq-9f03d983 -->

A shell entry in the `electron_shells` array has the form:

```json
{
  "function_type": "gto",
  "angular_momentum": [0],
  "exponents": ["3.4252509", "0.6239137", "0.1688554"],
  "coefficients": [["0.1543290", "0.5353281", "0.4446345"]]
}
```

- `angular_momentum` is a list of integers; one entry for a pure shell, two entries (e.g. `[0, 1]`)
  for an SP shell.
- `coefficients` is a list of lists, with one inner list per angular momentum entry. Each inner
  list has the same length as `exponents`.

---

## Gherkin Scenarios <!-- rq-b3b9e833 -->

```gherkin
Feature: Parse QCSchema basis set file

  # --- Happy paths ---

  @rq-5c07644f
  Scenario: Parse a valid single-shell file
    Given a valid QCSchema file for element "H" (Z=1) with one s-shell
    And the shell has 3 exponents and 3 coefficients
    When parse_basis is called with the file path
    Then the result is Ok(BasisSet)
    And BasisSet.element is "H"
    And BasisSet.atomic_number is 1
    And BasisSet.shells has 1 entry
    And shells[0].angular_momentum is 0
    And shells[0].exponents has 3 values
    And shells[0].coefficients has 3 values

  @rq-150a9f68
  Scenario: Parse a file with multiple shells
    Given a valid QCSchema file for element "C" (Z=6) with one s-shell and one p-shell
    When parse_basis is called with the file path
    Then BasisSet.shells has 2 entries
    And shells[0].angular_momentum is 0
    And shells[1].angular_momentum is 1

  @rq-61bece13
  Scenario: SP shell is split into separate S and P shells
    Given a valid QCSchema file for element "Li" (Z=3) with one SP shell
    And the SP shell has angular_momentum [0, 1], 3 exponents, and 2 coefficient vectors
    When parse_basis is called with the file path
    Then BasisSet.shells has 2 entries
    And shells[0].angular_momentum is 0
    And shells[0].exponents equals the original exponents
    And shells[0].coefficients equals the first coefficient vector
    And shells[1].angular_momentum is 1
    And shells[1].exponents equals the original exponents
    And shells[1].coefficients equals the second coefficient vector

  @rq-1509641c
  Scenario: ECP data is ignored when electron shells are present
    Given a valid QCSchema file for element "Cu" (Z=29) with one s-shell and an ecp_potentials key
    When parse_basis is called with the file path
    Then the result is Ok(BasisSet)
    And BasisSet.shells contains only the electron shell

  # --- File and JSON errors ---

  @rq-e18403b7
  Scenario: File does not exist
    Given no file exists at the given path
    When parse_basis is called with that path
    Then parse_basis returns Err(ParseError::IoError(_))

  @rq-bf33f678
  Scenario: File is not valid JSON
    Given a file containing "{ not valid json }"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::InvalidJson(_))

  # --- Element key errors ---

  @rq-d174c435
  Scenario: elements object contains two keys
    Given a QCSchema file whose elements object has keys "1" and "2"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MultipleElements { found: 2 })

  @rq-f141bbbf
  Scenario: elements object is empty
    Given a QCSchema file whose elements object is {}
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::NoElements)

  @rq-524efd00
  Scenario: Atomic number key is not a number
    Given a QCSchema file whose elements object has key "X"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::InvalidAtomicNumber("X"))

  @rq-f3750503
  Scenario: Atomic number key is zero
    Given a QCSchema file whose elements object has key "0"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::InvalidAtomicNumber("0"))

  @rq-bb7411a6
  Scenario: Atomic number key is out of range
    Given a QCSchema file whose elements object has key "119"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::InvalidAtomicNumber("119"))

  # --- Electron shell errors ---

  @rq-4e6861cf
  Scenario: electron_shells key is absent
    Given a QCSchema file where the element entry has no electron_shells key
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::NoElectronShells)

  @rq-5b34bd39
  Scenario: electron_shells array is empty
    Given a QCSchema file where electron_shells is []
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::NoElectronShells)

  # --- Malformed shell errors ---

  @rq-69dc5ba5
  Scenario: Shell has empty angular_momentum array
    Given a QCSchema file where shell 0 has angular_momentum []
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 0, .. })

  @rq-2aae5416
  Scenario: SP shell has wrong number of coefficient vectors
    Given a QCSchema file where shell 0 has angular_momentum [0, 1] but only 1 coefficient vector
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 0, .. })

  @rq-51537e9f
  Scenario: Coefficient vector length does not match exponent count
    Given a QCSchema file where shell 0 has 3 exponents but a coefficient vector of length 2
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 0, .. })

  @rq-3866b4e9
  Scenario: Exponent string cannot be parsed as f64
    Given a QCSchema file where shell 0 has an exponent value "abc"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 0, .. })

  @rq-45119d86
  Scenario: Coefficient string cannot be parsed as f64
    Given a QCSchema file where shell 0 has a coefficient value "abc"
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 0, .. })

  @rq-08ba2d8b
  Scenario: Error report identifies the correct shell index
    Given a QCSchema file where shell 1 (the second shell) has an unparseable exponent
    When parse_basis is called with the file path
    Then parse_basis returns Err(ParseError::MalformedShell { index: 1, .. })

  # --- load_basis ---

  @rq-ee630e3a
  Scenario: load_basis fetches and parses successfully
    Given fetch_basis("H", "sto-3g") would succeed and return a valid cached file path
    And the file at that path is a valid QCSchema basis for H
    When load_basis("H", "sto-3g") is called
    Then the result is Ok(BasisSet) with element "H"

  @rq-da61179a
  Scenario: load_basis propagates a fetch error
    Given fetch_basis("H", "unknown-basis") would return Err(BseError::UnknownBasisSet(_))
    When load_basis("H", "unknown-basis") is called
    Then load_basis returns Err(LoadError::Fetch(BseError::UnknownBasisSet(_)))

  @rq-cbf624d5
  Scenario: load_basis propagates a parse error
    Given fetch_basis("H", "sto-3g") succeeds and returns a path to a corrupted JSON file
    When load_basis("H", "sto-3g") is called
    Then load_basis returns Err(LoadError::Parse(_))
```
