use crate::basis::{BasisSet, LoadError};
use crate::input::CartesianGeometry;

// ─── Public types ─────────────────────────────────────────────────────────────

/// Structure-of-arrays representation of the contracted Cartesian AO basis.
///
/// All "per basis function" vectors have length `n_basis`.
/// All "per shell" vectors have length `n_shells`.
/// Primitive arrays (`exponents`, `coefficients`) have total length equal to
/// the sum of all `n_primitives` entries.
#[derive(Debug)]
pub struct AoBasis {
  /// Total number of Cartesian basis functions.
  pub n_basis: usize,
  /// Total number of contracted shells.
  pub n_shells: usize,

  // Per basis function
  pub center_x: Vec<f64>,
  pub center_y: Vec<f64>,
  pub center_z: Vec<f64>,
  pub lx: Vec<u32>,
  pub ly: Vec<u32>,
  pub lz: Vec<u32>,
  pub shell_index: Vec<usize>,
  pub atom_index: Vec<usize>,

  // Per shell
  pub prim_offset: Vec<usize>,
  pub n_primitives: Vec<usize>,

  // Per primitive (flat)
  pub exponents: Vec<f64>,
  pub coefficients: Vec<f64>,
}

/// Error type returned by `init_basis`.
#[derive(Debug)]
pub enum InitError {
  /// `load_basis` failed for the named element.
  BasisLoad { element: String, source: LoadError },
}

// ─── Public functions ─────────────────────────────────────────────────────────

/// Builds the contracted Cartesian AO basis from a molecular geometry and
/// a named basis set.
///
/// Calls `load_basis` once per unique element symbol and returns a flat,
/// structure-of-arrays representation ordered atom-major, shell-minor, with
/// Cartesian components ordered by descending lx, then descending ly.
pub fn init_basis(
  geometry: &CartesianGeometry,
  basis_name: &str,
) -> Result<AoBasis, InitError> {
  init_basis_impl(geometry, |element| {
    crate::basis::load_basis(element, basis_name).map_err(|e| InitError::BasisLoad {
      element: element.to_string(),
      source: e,
    })
  })
}

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Testable core: the load function is injected so tests can bypass I/O.
fn init_basis_impl<F>(
  geometry: &CartesianGeometry,
  load_fn: F,
) -> Result<AoBasis, InitError>
where
  F: Fn(&str) -> Result<BasisSet, InitError>,
{
  // Load a BasisSet for each unique element symbol, in first-occurrence order.
  let mut element_basis: std::collections::HashMap<String, BasisSet> =
    std::collections::HashMap::new();
  for symbol in &geometry.symbols {
    if !element_basis.contains_key(symbol) {
      let bs = load_fn(symbol)?;
      element_basis.insert(symbol.clone(), bs);
    }
  }

  let mut n_basis = 0usize;
  let mut n_shells = 0usize;
  let mut center_x = Vec::new();
  let mut center_y = Vec::new();
  let mut center_z = Vec::new();
  let mut lx_vec = Vec::new();
  let mut ly_vec = Vec::new();
  let mut lz_vec = Vec::new();
  let mut shell_index_vec = Vec::new();
  let mut atom_index_vec = Vec::new();
  let mut prim_offset_vec = Vec::new();
  let mut n_primitives_vec = Vec::new();
  let mut exponents_vec = Vec::new();
  let mut coefficients_vec = Vec::new();

  let mut prim_offset = 0usize;

  for (atom_idx, symbol) in geometry.symbols.iter().enumerate() {
    let cx = geometry.x[atom_idx];
    let cy = geometry.y[atom_idx];
    let cz = geometry.z[atom_idx];
    let bs = &element_basis[symbol];

    for shell in &bs.shells {
      let shell_idx = n_shells;
      n_shells += 1;

      let n_prim = shell.exponents.len();
      prim_offset_vec.push(prim_offset);
      n_primitives_vec.push(n_prim);
      exponents_vec.extend_from_slice(&shell.exponents);
      coefficients_vec.extend_from_slice(&shell.coefficients);
      prim_offset += n_prim;

      for (lx, ly, lz) in cartesian_components(shell.angular_momentum) {
        center_x.push(cx);
        center_y.push(cy);
        center_z.push(cz);
        lx_vec.push(lx);
        ly_vec.push(ly);
        lz_vec.push(lz);
        shell_index_vec.push(shell_idx);
        atom_index_vec.push(atom_idx);
        n_basis += 1;
      }
    }
  }

  Ok(AoBasis {
    n_basis,
    n_shells,
    center_x,
    center_y,
    center_z,
    lx: lx_vec,
    ly: ly_vec,
    lz: lz_vec,
    shell_index: shell_index_vec,
    atom_index: atom_index_vec,
    prim_offset: prim_offset_vec,
    n_primitives: n_primitives_vec,
    exponents: exponents_vec,
    coefficients: coefficients_vec,
  })
}

/// Returns the number of Cartesian components for angular momentum `l`.
/// n_cart(l) = (l+1)(l+2)/2.
fn n_cart(l: u32) -> usize {
  ((l + 1) * (l + 2) / 2) as usize
}

/// Returns all (lx, ly, lz) triples for angular momentum `l` in canonical
/// order: lx descending from l to 0; for each lx, ly descending from l−lx
/// to 0; lz = l − lx − ly.
fn cartesian_components(l: u32) -> Vec<(u32, u32, u32)> {
  let mut result = Vec::with_capacity(n_cart(l));
  let mut lx = l;
  loop {
    let mut ly = l - lx;
    loop {
      result.push((lx, ly, l - lx - ly));
      if ly == 0 { break; }
      ly -= 1;
    }
    if lx == 0 { break; }
    lx -= 1;
  }
  result
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
  use super::*;
  use crate::basis::{BasisSet, BseError, ElectronShell, LoadError};
  use crate::input::CartesianGeometry;
  use std::sync::Arc;
  use std::sync::atomic::{AtomicUsize, Ordering};

  // ── Test helpers ────────────────────────────────────────────────────────────

  fn shell(l: u32, exponents: Vec<f64>, coefficients: Vec<f64>) -> ElectronShell {
    ElectronShell { angular_momentum: l, exponents, coefficients }
  }

  fn uniform_shell(l: u32, n_prim: usize) -> ElectronShell {
    shell(l, vec![1.0; n_prim], vec![1.0; n_prim])
  }

  fn make_basis(element: &str, shells: Vec<ElectronShell>) -> BasisSet {
    BasisSet { element: element.to_string(), atomic_number: 1, shells }
  }

  fn geometry(symbols: Vec<&str>, xs: Vec<f64>, ys: Vec<f64>, zs: Vec<f64>) -> CartesianGeometry {
    CartesianGeometry {
      symbols: symbols.into_iter().map(|s| s.to_string()).collect(),
      x: xs,
      y: ys,
      z: zs,
    }
  }

  fn single_atom(sym: &str, x: f64, y: f64, z: f64) -> CartesianGeometry {
    geometry(vec![sym], vec![x], vec![y], vec![z])
  }

  // Load function that always returns Ok with the basis from the map.
  // For tests that only have one element.
  fn fixed_load(bs: BasisSet) -> impl Fn(&str) -> Result<BasisSet, InitError> {
    move |_| Ok(bs.clone())
  }

  // ── Edge cases ──────────────────────────────────────────────────────────────

  /// Scenario: Empty molecule returns an empty AoBasis
  #[test]
  fn empty_molecule_returns_empty_basis() {
    let geom = geometry(vec![], vec![], vec![], vec![]);
    let result = init_basis_impl(&geom, |e| panic!("load_fn called for {}", e));
    let b = result.expect("should succeed");
    assert_eq!(b.n_basis, 0);
    assert_eq!(b.n_shells, 0);
    assert!(b.center_x.is_empty());
    assert!(b.center_y.is_empty());
    assert!(b.center_z.is_empty());
    assert!(b.lx.is_empty());
    assert!(b.shell_index.is_empty());
    assert!(b.atom_index.is_empty());
    assert!(b.prim_offset.is_empty());
    assert!(b.n_primitives.is_empty());
    assert!(b.exponents.is_empty());
    assert!(b.coefficients.is_empty());
  }

  // ── Basis function count ─────────────────────────────────────────────────────

  /// Scenario: Single hydrogen atom with one s-shell gives one basis function
  #[test]
  fn single_h_one_s_shell_gives_one_basis_function() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let bs = make_basis("H", vec![uniform_shell(0, 3)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 1);
    assert_eq!(b.n_shells, 1);
  }

  /// Scenario: Single carbon atom with one s-shell and one p-shell gives four basis functions
  #[test]
  fn single_c_one_s_one_p_gives_four_basis_functions() {
    let geom = single_atom("C", 0.0, 0.0, 0.0);
    let bs = make_basis("C", vec![uniform_shell(0, 3), uniform_shell(1, 3)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 4);
    assert_eq!(b.n_shells, 2);
  }

  /// Scenario: Single d-shell gives six basis functions
  #[test]
  fn single_d_shell_gives_six_basis_functions() {
    let geom = single_atom("X", 0.0, 0.0, 0.0);
    let bs = make_basis("X", vec![uniform_shell(2, 2)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 6);
    assert_eq!(b.n_shells, 1);
  }

  /// Scenario: Two-atom molecule sums basis function counts from both atoms
  /// O has 2 s-shells and 1 p-shell; H has 1 s-shell → n_basis=6, n_shells=4
  #[test]
  fn two_atom_molecule_sums_basis_function_counts() {
    let geom = geometry(
      vec!["O", "H"],
      vec![0.0, 1.0],
      vec![0.0, 0.0],
      vec![0.0, 0.0],
    );
    let o_basis = make_basis("O", vec![uniform_shell(0, 3), uniform_shell(0, 3), uniform_shell(1, 3)]);
    let h_basis = make_basis("H", vec![uniform_shell(0, 3)]);
    let b = init_basis_impl(&geom, |e| {
      let bs = if e == "O" { o_basis.clone() } else { h_basis.clone() };
      Ok(bs)
    }).expect("should succeed");
    assert_eq!(b.n_basis, 6);
    assert_eq!(b.n_shells, 4);
  }

  // ── Atom-major ordering ──────────────────────────────────────────────────────

  /// Scenario: Basis functions for atom 0 appear before those for atom 1
  #[test]
  fn basis_functions_for_atom0_before_atom1() {
    let geom = geometry(
      vec!["H", "H"],
      vec![0.0, 1.0],
      vec![0.0, 0.0],
      vec![0.0, 0.0],
    );
    let bs = make_basis("H", vec![uniform_shell(0, 3)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 2);
    assert!((b.center_x[0] - 0.0).abs() < 1e-12);
    assert!((b.center_x[1] - 1.0).abs() < 1e-12);
  }

  /// Scenario: atom_index records the correct atom for each basis function
  #[test]
  fn atom_index_records_correct_atom() {
    // O (2 s-shells + 1 p-shell = 5 functions), H (1 s-shell = 1), H (1 s-shell = 1)
    let geom = geometry(
      vec!["O", "H", "H"],
      vec![0.0, 1.0, 2.0],
      vec![0.0, 0.0, 0.0],
      vec![0.0, 0.0, 0.0],
    );
    let o_basis = make_basis("O", vec![
      uniform_shell(0, 3),
      uniform_shell(0, 3),
      uniform_shell(1, 3),
    ]);
    let h_basis = make_basis("H", vec![uniform_shell(0, 3)]);
    let b = init_basis_impl(&geom, |e| {
      let bs = if e == "O" { o_basis.clone() } else { h_basis.clone() };
      Ok(bs)
    }).expect("should succeed");
    assert_eq!(b.n_basis, 7);
    assert_eq!(b.atom_index[0], 0);
    assert_eq!(b.atom_index[1], 0);
    assert_eq!(b.atom_index[2], 0);
    assert_eq!(b.atom_index[3], 0);
    assert_eq!(b.atom_index[4], 0);
    assert_eq!(b.atom_index[5], 1);
    assert_eq!(b.atom_index[6], 2);
  }

  // ── Shell ordering within an atom ────────────────────────────────────────────

  /// Scenario: Shells for one atom appear in basis-file order
  /// Basis [s, p, s]: shell_index = [0, 1, 1, 1, 2]
  #[test]
  fn shells_appear_in_basis_file_order() {
    let geom = single_atom("X", 0.0, 0.0, 0.0);
    let bs = make_basis("X", vec![
      uniform_shell(0, 3), // shell 0: 1 function
      uniform_shell(1, 3), // shell 1: 3 functions
      uniform_shell(0, 3), // shell 2: 1 function
    ]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 5);
    assert_eq!(b.shell_index[0], 0);
    assert_eq!(b.shell_index[1], 1);
    assert_eq!(b.shell_index[2], 1);
    assert_eq!(b.shell_index[3], 1);
    assert_eq!(b.shell_index[4], 2);
  }

  // ── Cartesian component ordering ─────────────────────────────────────────────

  /// Scenario: s-shell emits component (0, 0, 0)
  #[test]
  fn s_shell_emits_one_component() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let bs = make_basis("H", vec![uniform_shell(0, 1)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.lx[0], 0);
    assert_eq!(b.ly[0], 0);
    assert_eq!(b.lz[0], 0);
  }

  /// Scenario: p-shell emits components in order (1,0,0), (0,1,0), (0,0,1)
  #[test]
  fn p_shell_emits_components_in_correct_order() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let bs = make_basis("H", vec![uniform_shell(1, 1)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 3);
    assert_eq!((b.lx[0], b.ly[0], b.lz[0]), (1, 0, 0));
    assert_eq!((b.lx[1], b.ly[1], b.lz[1]), (0, 1, 0));
    assert_eq!((b.lx[2], b.ly[2], b.lz[2]), (0, 0, 1));
  }

  /// Scenario: d-shell emits six components in the correct order
  #[test]
  fn d_shell_emits_six_components_in_correct_order() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let bs = make_basis("H", vec![uniform_shell(2, 1)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 6);
    let expected = [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)];
    for (i, (ex, ey, ez)) in expected.iter().enumerate() {
      assert_eq!((b.lx[i], b.ly[i], b.lz[i]), (*ex, *ey, *ez),
        "mismatch at index {i}");
    }
  }

  // ── Primitive storage ────────────────────────────────────────────────────────

  /// Scenario: Exponents and coefficients are copied from the ElectronShell
  #[test]
  fn exponents_and_coefficients_are_copied() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let exps = vec![3.425, 0.624, 0.169];
    let coeffs = vec![0.154, 0.535, 0.445];
    let bs = make_basis("H", vec![shell(0, exps.clone(), coeffs.clone())]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_primitives[0], 3);
    assert_eq!(b.prim_offset[0], 0);
    for i in 0..3 {
      assert!((b.exponents[i] - exps[i]).abs() < 1e-9);
      assert!((b.coefficients[i] - coeffs[i]).abs() < 1e-9);
    }
  }

  /// Scenario: prim_offset is correct when multiple shells are present
  #[test]
  fn prim_offset_is_correct_for_multiple_shells() {
    let geom = single_atom("X", 0.0, 0.0, 0.0);
    let bs = make_basis("X", vec![
      uniform_shell(0, 3), // shell 0: 3 primitives
      uniform_shell(0, 2), // shell 1: 2 primitives
    ]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.prim_offset[0], 0);
    assert_eq!(b.prim_offset[1], 3);
  }

  /// Scenario: All Cartesian functions of a shell share the same shell_index
  #[test]
  fn all_cartesian_functions_of_shell_share_shell_index() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let bs = make_basis("H", vec![uniform_shell(1, 1)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 3);
    assert_eq!(b.shell_index[0], 0);
    assert_eq!(b.shell_index[1], 0);
    assert_eq!(b.shell_index[2], 0);
  }

  // ── Atom coordinates ─────────────────────────────────────────────────────────

  /// Scenario: center_x/y/z are taken from the CartesianGeometry coordinates
  #[test]
  fn center_coordinates_taken_from_geometry() {
    let geom = single_atom("H", 1.5, 2.5, 3.5);
    let bs = make_basis("H", vec![uniform_shell(0, 1)]);
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert!((b.center_x[0] - 1.5).abs() < 1e-12);
    assert!((b.center_y[0] - 2.5).abs() < 1e-12);
    assert!((b.center_z[0] - 3.5).abs() < 1e-12);
  }

  /// Scenario: All basis functions for the same atom share that atom's coordinates
  #[test]
  fn all_functions_for_same_atom_share_coordinates() {
    let geom = single_atom("X", 1.0, 0.0, 0.0);
    let bs = make_basis("X", vec![uniform_shell(1, 1)]); // p-shell: 3 functions
    let b = init_basis_impl(&geom, fixed_load(bs)).expect("should succeed");
    assert_eq!(b.n_basis, 3);
    for i in 0..3 {
      assert!((b.center_x[i] - 1.0).abs() < 1e-12, "center_x[{i}] wrong");
      assert!((b.center_y[i] - 0.0).abs() < 1e-12, "center_y[{i}] wrong");
    }
  }

  // ── Repeated element types ───────────────────────────────────────────────────

  /// Scenario: load_basis is called once per unique element, not once per atom
  #[test]
  fn load_basis_called_once_per_unique_element() {
    let geom = geometry(
      vec!["H", "H"],
      vec![0.0, 1.0],
      vec![0.0, 0.0],
      vec![0.0, 0.0],
    );
    let h_basis = make_basis("H", vec![
      shell(0, vec![3.425, 0.624, 0.169], vec![0.154, 0.535, 0.445]),
    ]);
    let call_count = Arc::new(AtomicUsize::new(0));
    let count_clone = call_count.clone();
    let bs_clone = h_basis.clone();
    let b = init_basis_impl(&geom, move |_| {
      count_clone.fetch_add(1, Ordering::SeqCst);
      Ok(bs_clone.clone())
    }).expect("should succeed");

    assert_eq!(call_count.load(Ordering::SeqCst), 1, "load_fn should be called once");
    // Both H atoms share the same exponents
    assert_eq!(b.n_basis, 2);
    for atom_offset in [0, 1] {
      // shell primitive data is the same for both atoms
      let prim_start = b.prim_offset[atom_offset];
      assert!((b.exponents[prim_start] - 3.425).abs() < 1e-6);
    }
  }

  // ── Error handling ───────────────────────────────────────────────────────────

  /// Scenario: InitError::BasisLoad is returned when load_basis fails for an element
  #[test]
  fn basis_load_error_returned_on_failure() {
    let geom = single_atom("H", 0.0, 0.0, 0.0);
    let result = init_basis_impl(&geom, |element| {
      Err(InitError::BasisLoad {
        element: element.to_string(),
        source: LoadError::Fetch(BseError::UnknownBasisSet("unknown-basis".to_string())),
      })
    });
    assert!(
      matches!(result, Err(InitError::BasisLoad { ref element, .. }) if element == "H"),
      "expected BasisLoad error for H"
    );
  }

  /// Scenario: The element name in BasisLoad error identifies the failing element
  #[test]
  fn basis_load_error_identifies_failing_element() {
    let geom = geometry(
      vec!["H", "C"],
      vec![0.0, 1.0],
      vec![0.0, 0.0],
      vec![0.0, 0.0],
    );
    let h_basis = make_basis("H", vec![uniform_shell(0, 3)]);
    let result = init_basis_impl(&geom, move |element| {
      if element == "H" {
        Ok(h_basis.clone())
      } else {
        Err(InitError::BasisLoad {
          element: element.to_string(),
          source: LoadError::Fetch(BseError::UnknownBasisSet("sto-3g".to_string())),
        })
      }
    });
    assert!(
      matches!(result, Err(InitError::BasisLoad { ref element, .. }) if element == "C"),
      "expected BasisLoad error identifying element C"
    );
  }
}
