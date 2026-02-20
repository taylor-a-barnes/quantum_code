use faer::{Mat, Side};

// ── Public types ──────────────────────────────────────────────────────────────

/// Error type returned by `guess_hcore`.
#[derive(Debug, PartialEq)]
pub enum GuessError {
  /// The three input matrices do not all have the same n × n shape. Each field
  /// records the (nrows, ncols) of the corresponding matrix.
  DimensionMismatch {
    s_shape: (usize, usize),
    t_shape: (usize, usize),
    v_shape: (usize, usize),
  },
  /// `n_alpha > n_basis` or `n_beta > n_basis`.
  TooManyElectrons {
    n_alpha: usize,
    n_beta: usize,
    n_basis: usize,
  },
  /// S is not positive definite; Cholesky factorisation would fail.
  SingularOverlap,
}

// ── Public functions ──────────────────────────────────────────────────────────

/// Returns the initial MO coefficient matrix C (n_basis × n_basis) using the
/// core Hamiltonian diagonalisation method.
///
/// Algorithm (canonical orthogonalisation):
///   1. S = U_s Λ_s U_s^T    (eigendecomposition of S)
///   2. X = U_s Λ_s^{−1/2}   (orthogonaliser; equivalent to L^{−T} from Cholesky)
///   3. H' = X^T (T + V) X
///   4. H' = U' ε U'^T       (eigendecomposition, sorted ascending in ε)
///   5. C = X U'
///
/// The returned C satisfies C^T S C = I. Its columns are the MOs sorted by
/// ascending orbital energy; n_alpha and n_beta select which are occupied.
pub fn guess_hcore(
  s: &Mat<f64>,
  t: &Mat<f64>,
  v: &Mat<f64>,
  n_alpha: usize,
  n_beta: usize,
) -> Result<Mat<f64>, GuessError> {
  let s_shape = (s.nrows(), s.ncols());
  let t_shape = (t.nrows(), t.ncols());
  let v_shape = (v.nrows(), v.ncols());

  // Dimension check: all three must be square and share the same side length.
  let n = s.nrows();
  if s.nrows() != s.ncols()
    || t.nrows() != t.ncols()
    || v.nrows() != v.ncols()
    || t.nrows() != n
    || v.nrows() != n
  {
    return Err(GuessError::DimensionMismatch { s_shape, t_shape, v_shape });
  }

  // Electron count check (performed after dimension check).
  if n_alpha > n || n_beta > n {
    return Err(GuessError::TooManyElectrons { n_alpha, n_beta, n_basis: n });
  }

  // Degenerate case: no basis functions → empty coefficient matrix.
  if n == 0 {
    return Ok(Mat::zeros(0, 0));
  }

  // H_core = T + V.
  let h_core: Mat<f64> = t + v;

  // Eigendecompose S: S = U_s Λ_s U_s^T.
  let evd_s = s.selfadjoint_eigendecomposition(Side::Lower);
  let lambdas_s: Vec<f64> = (0..n)
    .map(|i| evd_s.s().column_vector().read(i))
    .collect();

  // S must be positive definite (all eigenvalues strictly positive).
  if lambdas_s.iter().any(|&l| l <= 0.0) {
    return Err(GuessError::SingularOverlap);
  }

  // Build orthogonaliser X = U_s * diag(λ_s^{−1/2}).
  // Column j of X is the j-th eigenvector of S scaled by λ_s[j]^{-1/2}.
  let u_s = evd_s.u().to_owned();
  let mut x: Mat<f64> = Mat::zeros(n, n);
  for j in 0..n {
    let scale = lambdas_s[j].powf(-0.5);
    for i in 0..n {
      x.write(i, j, u_s[(i, j)] * scale);
    }
  }

  // H' = X^T H_core X.
  let xt_hc: Mat<f64> = x.transpose() * &h_core;
  let h_prime: Mat<f64> = &xt_hc * &x;

  // Eigendecompose H' and obtain indices sorted by ascending orbital energy.
  let evd_h = h_prime.selfadjoint_eigendecomposition(Side::Lower);
  let energies: Vec<f64> = (0..n)
    .map(|i| evd_h.s().column_vector().read(i))
    .collect();

  let mut order: Vec<usize> = (0..n).collect();
  order.sort_by(|&a, &b| {
    energies[a]
      .partial_cmp(&energies[b])
      .unwrap_or(std::cmp::Ordering::Equal)
  });

  // Build U' with columns reordered from low to high energy.
  let u_prime = evd_h.u();
  let mut u_sorted: Mat<f64> = Mat::zeros(n, n);
  for (new_j, &old_j) in order.iter().enumerate() {
    for i in 0..n {
      u_sorted.write(i, new_j, u_prime[(i, old_j)]);
    }
  }

  // C = X U'_sorted.
  Ok(&x * &u_sorted)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
  use super::*;

  // Build an n×n Mat<f64> from a flat row-major array.
  fn mat2(data: [f64; 4]) -> Mat<f64> {
    Mat::from_fn(2, 2, |i, j| data[i * 2 + j])
  }
  fn mat3(data: [f64; 9]) -> Mat<f64> {
    Mat::from_fn(3, 3, |i, j| data[i * 3 + j])
  }
  fn mat4(data: [f64; 16]) -> Mat<f64> {
    Mat::from_fn(4, 4, |i, j| data[i * 4 + j])
  }

  // Verify C^T S C ≈ identity to within `tol`.
  fn assert_orthonormal(c: &Mat<f64>, s: &Mat<f64>, tol: f64) {
    let n = c.nrows();
    let ct_s: Mat<f64> = c.transpose() * s;
    let ct_s_c: Mat<f64> = &ct_s * c;
    for i in 0..n {
      for j in 0..n {
        let expected = if i == j { 1.0 } else { 0.0 };
        let diff = (ct_s_c[(i, j)] - expected).abs();
        assert!(
          diff < tol,
          "C^T S C [{i},{j}] = {:.6}, expected {expected:.1}, diff {diff:.2e}",
          ct_s_c[(i, j)]
        );
      }
    }
  }

  // H₂-like STO-3G matrices at ~1.4 Bohr.
  fn h2() -> (Mat<f64>, Mat<f64>, Mat<f64>) {
    let s = mat2([1.0, 0.5, 0.5, 1.0]);
    let t = mat2([0.760, 0.236, 0.236, 0.760]);
    let v = mat2([-1.883, -1.190, -1.190, -1.883]);
    (s, t, v)
  }

  // Simple 3×3 well-conditioned system (identity overlap).
  fn three_by_three() -> (Mat<f64>, Mat<f64>, Mat<f64>) {
    let s = Mat::identity(3, 3);
    let t = mat3([1.0, 0.2, 0.1, 0.2, 1.5, 0.3, 0.1, 0.3, 2.0]);
    let v = mat3([-3.0, -0.4, -0.2, -0.4, -2.5, -0.5, -0.2, -0.5, -2.0]);
    (s, t, v)
  }

  // ── Happy paths ────────────────────────────────────────────────────────────

  /// Scenario: 2×2 system produces orthonormal MOs (S-metric).
  #[test]
  fn two_by_two_c_is_orthonormal() {
    let (s, t, v) = h2();
    let c = guess_hcore(&s, &t, &v, 1, 1).expect("should succeed");
    assert_eq!((c.nrows(), c.ncols()), (2, 2));
    assert_orthonormal(&c, &s, 1e-6);
  }

  /// Scenario: Columns of C are sorted by ascending orbital energy.
  #[test]
  fn columns_sorted_by_ascending_energy() {
    // S = I, H_core = diag(0.2, -1.5).
    // Natural eigenvalue order from faer is unspecified, but guess_hcore must
    // sort ascending, so column 0 → -1.5, column 1 → 0.2.
    let s = Mat::<f64>::identity(2, 2);
    let t = mat2([0.2, 0.0, 0.0, 0.0]);
    let v = mat2([0.0, 0.0, 0.0, -1.5]);
    let c = guess_hcore(&s, &t, &v, 1, 1).expect("should succeed");
    let h_core = mat2([0.2, 0.0, 0.0, -1.5]);
    let temp: Mat<f64> = c.transpose() * &h_core;
    let ct_h_c: Mat<f64> = &temp * &c;
    assert!(
      ct_h_c[(0, 0)] < ct_h_c[(1, 1)],
      "first MO energy ({}) must be < second ({})",
      ct_h_c[(0, 0)],
      ct_h_c[(1, 1)]
    );
    assert!((ct_h_c[(0, 0)] - (-1.5)).abs() < 1e-6, "first MO energy ≈ -1.5");
    assert!((ct_h_c[(1, 1)] - 0.2).abs() < 1e-6, "second MO energy ≈ 0.2");
  }

  /// Scenario: S = identity gives C that diagonalises H_core.
  #[test]
  fn identity_overlap_diagonalises_h_core() {
    // H_core = diag(1.0, -2.0, -0.5); sorted eigenvalues: -2.0, -0.5, 1.0.
    let s = Mat::<f64>::identity(3, 3);
    let t = mat3([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let v = mat3([0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -0.5]);
    let c = guess_hcore(&s, &t, &v, 1, 1).expect("should succeed");
    let h_core = mat3([1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -0.5]);
    let temp: Mat<f64> = c.transpose() * &h_core;
    let ct_h_c: Mat<f64> = &temp * &c;
    for (k, &expected) in [-2.0_f64, -0.5, 1.0].iter().enumerate() {
      assert!(
        (ct_h_c[(k, k)] - expected).abs() < 1e-6,
        "diag[{k}] = {}, expected {expected}",
        ct_h_c[(k, k)]
      );
    }
  }

  /// Scenario: 1×1 system returns a 1×1 coefficient matrix.
  #[test]
  fn one_by_one_system() {
    let s = Mat::from_fn(1, 1, |_, _| 1.0_f64);
    let t = Mat::from_fn(1, 1, |_, _| 0.5_f64);
    let v = Mat::from_fn(1, 1, |_, _| -1.5_f64);
    let c = guess_hcore(&s, &t, &v, 1, 0).expect("should succeed");
    assert_eq!((c.nrows(), c.ncols()), (1, 1));
    // With S = [1], the normalised coefficient must be ±1.
    assert!((c[(0, 0)].abs() - 1.0).abs() < 1e-9, "C[0,0] must be ±1");
  }

  /// Scenario: Zero electrons is accepted; C is still a full orthonormal basis.
  #[test]
  fn zero_electrons_accepted() {
    let (s, t, v) = three_by_three();
    let c = guess_hcore(&s, &t, &v, 0, 0).expect("should succeed");
    assert_eq!((c.nrows(), c.ncols()), (3, 3));
    assert_orthonormal(&c, &s, 1e-6);
  }

  /// Scenario: n_alpha = n_beta = n_basis (fully occupied) is accepted.
  #[test]
  fn fully_occupied_accepted() {
    let (s, t, v) = three_by_three();
    let c = guess_hcore(&s, &t, &v, 3, 3).expect("should succeed");
    assert_eq!((c.nrows(), c.ncols()), (3, 3));
  }

  /// Scenario: Unrestricted (n_alpha ≠ n_beta) is accepted; C is orthonormal.
  #[test]
  fn unrestricted_spin_accepted() {
    let s = Mat::<f64>::identity(4, 4);
    let t = mat4([
      1.0, 0.1, 0.05, 0.02,
      0.1, 1.5, 0.08, 0.03,
      0.05, 0.08, 2.0, 0.1,
      0.02, 0.03, 0.1, 2.5,
    ]);
    let v = mat4([
      -3.0, -0.2, -0.1, -0.05,
      -0.2, -2.5, -0.15, -0.06,
      -0.1, -0.15, -2.0, -0.2,
      -0.05, -0.06, -0.2, -1.5,
    ]);
    let c = guess_hcore(&s, &t, &v, 3, 2).expect("should succeed");
    assert_eq!((c.nrows(), c.ncols()), (4, 4));
    assert_orthonormal(&c, &s, 1e-6);
  }

  // ── Dimension mismatch errors ──────────────────────────────────────────────

  /// Scenario: T has a different size from S → DimensionMismatch.
  #[test]
  fn t_different_size() {
    let s = Mat::<f64>::zeros(3, 3);
    let t = Mat::<f64>::zeros(2, 2);
    let v = Mat::<f64>::zeros(3, 3);
    assert_eq!(
      guess_hcore(&s, &t, &v, 0, 0),
      Err(GuessError::DimensionMismatch {
        s_shape: (3, 3),
        t_shape: (2, 2),
        v_shape: (3, 3),
      })
    );
  }

  /// Scenario: V has a different size from S → DimensionMismatch.
  #[test]
  fn v_different_size() {
    let s = Mat::<f64>::zeros(3, 3);
    let t = Mat::<f64>::zeros(3, 3);
    let v = Mat::<f64>::zeros(4, 4);
    assert_eq!(
      guess_hcore(&s, &t, &v, 0, 0),
      Err(GuessError::DimensionMismatch {
        s_shape: (3, 3),
        t_shape: (3, 3),
        v_shape: (4, 4),
      })
    );
  }

  /// Scenario: Non-square S → DimensionMismatch.
  #[test]
  fn non_square_s() {
    let s = Mat::<f64>::zeros(3, 2);
    let t = Mat::<f64>::zeros(3, 3);
    let v = Mat::<f64>::zeros(3, 3);
    assert_eq!(
      guess_hcore(&s, &t, &v, 0, 0),
      Err(GuessError::DimensionMismatch {
        s_shape: (3, 2),
        t_shape: (3, 3),
        v_shape: (3, 3),
      })
    );
  }

  /// Scenario: Non-square T → DimensionMismatch.
  #[test]
  fn non_square_t() {
    let s = Mat::<f64>::zeros(3, 3);
    let t = Mat::<f64>::zeros(3, 2);
    let v = Mat::<f64>::zeros(3, 3);
    assert_eq!(
      guess_hcore(&s, &t, &v, 0, 0),
      Err(GuessError::DimensionMismatch {
        s_shape: (3, 3),
        t_shape: (3, 2),
        v_shape: (3, 3),
      })
    );
  }

  // ── Electron count errors ──────────────────────────────────────────────────

  /// Scenario: n_alpha > n_basis → TooManyElectrons.
  #[test]
  fn n_alpha_exceeds_n_basis() {
    let (s, t, v) = three_by_three();
    assert_eq!(
      guess_hcore(&s, &t, &v, 4, 1),
      Err(GuessError::TooManyElectrons { n_alpha: 4, n_beta: 1, n_basis: 3 })
    );
  }

  /// Scenario: n_beta > n_basis → TooManyElectrons.
  #[test]
  fn n_beta_exceeds_n_basis() {
    let (s, t, v) = three_by_three();
    assert_eq!(
      guess_hcore(&s, &t, &v, 1, 4),
      Err(GuessError::TooManyElectrons { n_alpha: 1, n_beta: 4, n_basis: 3 })
    );
  }

  /// Scenario: Dimension mismatch takes priority over TooManyElectrons.
  #[test]
  fn dimension_mismatch_priority() {
    let s = Mat::<f64>::zeros(3, 3);
    let t = Mat::<f64>::zeros(2, 2);
    let v = Mat::<f64>::zeros(3, 3);
    assert!(matches!(
      guess_hcore(&s, &t, &v, 5, 5),
      Err(GuessError::DimensionMismatch { .. })
    ));
  }

  // ── Numerical errors ───────────────────────────────────────────────────────

  /// Scenario: Rank-deficient S → SingularOverlap.
  /// S = [[1,1],[1,1]] has eigenvalues 0 and 2; the zero eigenvalue triggers the error.
  #[test]
  fn singular_overlap_rank_deficient() {
    let s = mat2([1.0, 1.0, 1.0, 1.0]);
    let t = Mat::<f64>::identity(2, 2);
    let v = Mat::<f64>::zeros(2, 2);
    assert_eq!(guess_hcore(&s, &t, &v, 0, 0), Err(GuessError::SingularOverlap));
  }

  /// Scenario: Negative definite S → SingularOverlap.
  #[test]
  fn singular_overlap_negative_definite() {
    let s = mat2([-1.0, 0.0, 0.0, -1.0]);
    let t = Mat::<f64>::identity(2, 2);
    let v = Mat::<f64>::zeros(2, 2);
    assert_eq!(guess_hcore(&s, &t, &v, 0, 0), Err(GuessError::SingularOverlap));
  }
}
