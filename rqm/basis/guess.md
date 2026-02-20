# Feature: Create Initial Atomic Orbital Guess

This feature implements a function that generates an initial guess for the molecular orbital (MO)
coefficient matrix using the core Hamiltonian diagonalisation method. The returned matrix provides
the starting point for a self-consistent field (SCF) optimisation. For details about the atomic
orbital data structure, see `rqm/basis/initialization.md`.

## Background: Core Hamiltonian Guess

The core Hamiltonian H_core = T + V_ne (kinetic-energy + nuclear-electron attraction matrices) is
diagonalised in the AO basis via the symmetric generalised eigenvalue problem:

    H_core C = S C ε

where S is the AO overlap matrix, C is the MO coefficient matrix (columns = MOs), and ε is a
diagonal matrix of orbital energies. The columns of C are orthonormal in the S-metric:
C^T S C = I.

The Cholesky decomposition of S is used to reduce this to a standard symmetric eigenvalue problem:

1. Factorise S = L L^T  (L lower-triangular, Cholesky).
2. Form H' = L^{-1} H_core L^{-T}.
3. Diagonalise H' (symmetric eigendecomposition), obtaining eigenvectors U' and eigenvalues ε,
   sorted in ascending order.
4. Back-transform: C = L^{-T} U'.

The returned matrix C is n_basis × n_basis; each column is one MO, sorted left-to-right by
ascending orbital energy. The linear-algebra backend is `faer`.

---

## Feature API

### Functions

- `guess_hcore(s: &Mat<f64>, t: &Mat<f64>, v: &Mat<f64>, n_alpha: usize, n_beta: usize) -> Result<Mat<f64>, GuessError>`
  - `s` — n_basis × n_basis symmetric AO overlap matrix.
  - `t` — n_basis × n_basis symmetric kinetic-energy matrix.
  - `v` — n_basis × n_basis symmetric nuclear-attraction matrix.
  - `n_alpha`, `n_beta` — number of alpha and beta electrons; used only for input validation.
  - Forms H_core = T + V (element-wise sum), then follows the four-step algorithm above.
  - Returns the MO coefficient matrix C (n_basis × n_basis) with columns sorted by ascending
    orbital energy. The first `n_alpha` columns (resp. `n_beta`) are the occupied alpha
    (resp. beta) orbitals; the remainder are virtual.
  - All three matrices must be square and share the same dimension; returns `DimensionMismatch`
    otherwise.
  - Returns `TooManyElectrons` if `n_alpha > n_basis` or `n_beta > n_basis`.
  - Returns `SingularOverlap` if S is not positive definite (Cholesky fails).
  - `Mat<f64>` is `faer::Mat<f64>`.

### Types

- `GuessError` — error type returned by `guess_hcore`:
  - `DimensionMismatch { s_shape: (usize, usize), t_shape: (usize, usize), v_shape: (usize, usize) }`
    — the three input matrices do not all have the same n × n shape. Each field records the
    (nrows, ncols) of the corresponding matrix. Triggered by any non-square matrix or any size
    disagreement among S, T, and V.
  - `TooManyElectrons { n_alpha: usize, n_beta: usize, n_basis: usize }` — `n_alpha > n_basis`
    or `n_beta > n_basis`.
  - `SingularOverlap` — the Cholesky factorisation of S failed; S is not positive definite.

---

## Dependencies

This feature adds `faer` as a production dependency. The crate provides `faer::Mat<f64>` and the
symmetric eigendecomposition and Cholesky routines required by the algorithm.

---

## Gherkin Scenarios

```gherkin
Feature: Create initial HF/SCF guess from core Hamiltonian

  # Notation:
  #   S, T, V are always understood to be real symmetric matrices.
  #   H_core = T + V.
  #   All floating-point comparisons use a tolerance of 1e-9 unless noted.

  # --- Happy paths ---

  Scenario: 2×2 system produces orthonormal MOs (S-metric)
    Given S = [[1.0, 0.5], [0.5, 1.0]]
    And T = [[0.760, 0.236], [0.236, 0.760]]
    And V = [[-1.883, -1.190], [-1.190, -1.883]]
    And n_alpha = 1, n_beta = 1
    When guess_hcore(S, T, V, n_alpha, n_beta) is called
    Then the result is Ok(C) where C is a 2×2 matrix
    And C^T S C is approximately the 2×2 identity matrix

  Scenario: Columns of C are sorted by ascending orbital energy
    Given a 2×2 system where H_core has one eigenvalue near -1.5 and one near +0.2
    And n_alpha = 1, n_beta = 1
    When guess_hcore is called
    Then the first column of C corresponds to the lower-energy MO (≈ -1.5)
    And the second column of C corresponds to the higher-energy MO (≈ +0.2)

  Scenario: S equal to identity gives C that diagonalises H_core
    Given S is the 3×3 identity matrix
    And H_core is a known symmetric 3×3 matrix with eigenvalues [-2.0, -0.5, 1.0] (ascending)
    And n_alpha = 1, n_beta = 1
    When guess_hcore is called
    Then the result is Ok(C)
    And C^T H_core C is approximately diagonal with entries [-2.0, -0.5, 1.0]

  Scenario: 1×1 system returns a 1×1 coefficient matrix
    Given S = [[1.0]], T = [[0.5]], V = [[-1.5]]
    And n_alpha = 1, n_beta = 0
    When guess_hcore is called
    Then the result is Ok(C) where C is [[1.0]] (up to sign)

  Scenario: Zero electrons is accepted and all MOs are virtual
    Given a 3×3 well-conditioned system
    And n_alpha = 0, n_beta = 0
    When guess_hcore is called
    Then the result is Ok(C) and C is 3×3
    And C^T S C is approximately the 3×3 identity matrix

  Scenario: n_alpha = n_beta = n_basis (fully occupied) is accepted
    Given a 3×3 well-conditioned system
    And n_alpha = 3, n_beta = 3
    When guess_hcore is called
    Then the result is Ok(C) and C is 3×3

  Scenario: Unrestricted system (n_alpha ≠ n_beta) is accepted
    Given a 4×4 well-conditioned system
    And n_alpha = 3, n_beta = 2
    When guess_hcore is called
    Then the result is Ok(C) and C is 4×4
    And C^T S C is approximately the 4×4 identity matrix

  # --- Dimension mismatch errors ---

  Scenario: T has a different size from S returns DimensionMismatch
    Given S is 3×3, T is 2×2, and V is 3×3
    When guess_hcore is called
    Then the result is Err(GuessError::DimensionMismatch {
        s_shape: (3, 3), t_shape: (2, 2), v_shape: (3, 3) })

  Scenario: V has a different size from S returns DimensionMismatch
    Given S is 3×3, T is 3×3, and V is 4×4
    When guess_hcore is called
    Then the result is Err(GuessError::DimensionMismatch {
        s_shape: (3, 3), t_shape: (3, 3), v_shape: (4, 4) })

  Scenario: Non-square S returns DimensionMismatch
    Given S has 3 rows and 2 columns, T is 3×3, and V is 3×3
    When guess_hcore is called
    Then the result is Err(GuessError::DimensionMismatch {
        s_shape: (3, 2), t_shape: (3, 3), v_shape: (3, 3) })

  Scenario: Non-square T returns DimensionMismatch
    Given S is 3×3, T has 3 rows and 2 columns, and V is 3×3
    When guess_hcore is called
    Then the result is Err(GuessError::DimensionMismatch {
        s_shape: (3, 3), t_shape: (3, 2), v_shape: (3, 3) })

  # --- Electron count errors ---

  Scenario: n_alpha exceeds n_basis returns TooManyElectrons
    Given a 3×3 well-conditioned system
    And n_alpha = 4, n_beta = 1
    When guess_hcore is called
    Then the result is Err(GuessError::TooManyElectrons { n_alpha: 4, n_beta: 1, n_basis: 3 })

  Scenario: n_beta exceeds n_basis returns TooManyElectrons
    Given a 3×3 well-conditioned system
    And n_alpha = 1, n_beta = 4
    When guess_hcore is called
    Then the result is Err(GuessError::TooManyElectrons { n_alpha: 1, n_beta: 4, n_basis: 3 })

  # Dimension checks are performed before electron-count checks.
  Scenario: Dimension mismatch takes priority over TooManyElectrons
    Given S is 3×3, T is 2×2, and V is 3×3
    And n_alpha = 5, n_beta = 5
    When guess_hcore is called
    Then the result is Err(GuessError::DimensionMismatch { .. })

  # --- Numerical errors ---

  Scenario: Singular overlap matrix (rank deficient) returns SingularOverlap
    Given S = [[1.0, 1.0], [1.0, 1.0]] (rank 1, not positive definite)
    And T and V are well-formed 2×2 matrices
    And n_alpha = 0, n_beta = 0
    When guess_hcore is called
    Then the result is Err(GuessError::SingularOverlap)

  Scenario: Negative definite S (diagonal entries negative) returns SingularOverlap
    Given S = [[-1.0, 0.0], [0.0, -1.0]]
    And T and V are well-formed 2×2 matrices
    When guess_hcore is called
    Then the result is Err(GuessError::SingularOverlap)
```
