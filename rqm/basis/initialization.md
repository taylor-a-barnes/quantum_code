# Feature: Initialization of the AO Basis

This feature builds the contracted Cartesian Gaussian-type orbital (CGTO) basis from a molecular
geometry and a named basis set. It consumes a `CartesianGeometry` (from `rqm/input/parser.md`) and
the basis set name, calls `load_basis` (from `rqm/basis/parser.md`) once per unique element, and
produces a flat, structure-of-arrays representation suitable for use in integral evaluation and
other high-performance routines.

Z-matrix geometry is not accepted; the caller is responsible for converting to Cartesian
coordinates before calling this feature.

## Background: Cartesian GTOs

A contracted Cartesian GTO centered on atom A with nuclear position **R**_A is:

    φ(r) = (x − Rx)^lx (y − Ry)^ly (z − Rz)^lz  Σ_i  c_i exp(−α_i |r − R_A|²)

where lx + ly + lz = l (the total angular momentum of the shell) and the sum runs over all
primitives in the contracted shell.

For a shell of angular momentum l there are n_cart(l) = (l+1)(l+2)/2 Cartesian components:

| l | Label | n_cart | Components (lx, ly, lz) |
| - | ----- | ------ | ----------------------- |
| 0 | s     | 1      | (0,0,0) |
| 1 | p     | 3      | (1,0,0) (0,1,0) (0,0,1) |
| 2 | d     | 6      | (2,0,0) (1,1,0) (1,0,1) (0,2,0) (0,1,1) (0,0,2) |
| 3 | f     | 10     | … |

The ordering within a shell is: lx descending from l to 0; for each lx, ly descending from l−lx
to 0; lz = l − lx − ly.

## Basis function ordering

Basis functions are stored in **atom-major, shell-minor** order:

1. Iterate over atoms in the order they appear in `CartesianGeometry.symbols`.
2. For each atom, iterate over the shells of that element's `BasisSet` in file order.
3. For each shell, emit all n_cart(l) Cartesian components in the order defined above.

`load_basis` is called once for each **unique** element symbol that appears in the molecule. The
same `BasisSet` is reused for every atom of that element.

---

## Feature API

### Functions

- `init_basis(geometry: &CartesianGeometry, basis_name: &str) -> Result<AoBasis, InitError>`
  - Collects the unique element symbols from `geometry.symbols`.
  - For each unique element calls `load_basis(element, basis_name)`, propagating any `LoadError`
    as `InitError::BasisLoad`.
  - Iterates over atoms in order, then shells in file order, then Cartesian components as defined
    above, building all parallel SoA vectors simultaneously.
  - Returns an `AoBasis` whose length fields (`n_basis`, `n_shells`) match the counts implied by
    the molecule and the loaded basis sets.

### Types

- `AoBasis` — structure of arrays; all "per basis function" vectors have length `n_basis`; all
  "per shell" vectors have length `n_shells`; primitive vectors have total length equal to the sum
  of all `n_primitives` entries.

  **Scalar counts**
  - `n_basis: usize` — total number of Cartesian basis functions.
  - `n_shells: usize` — total number of contracted shells (sum over all atoms of the number of
    shells for that element's basis set).

  **Per basis function** (index `i` in 0..n_basis)
  - `center_x: Vec<f64>` — x-coordinate of the atom this function is centred on (Bohr).
  - `center_y: Vec<f64>` — y-coordinate (Bohr).
  - `center_z: Vec<f64>` — z-coordinate (Bohr).
  - `lx: Vec<u32>` — Cartesian angular momentum exponent for x.
  - `ly: Vec<u32>` — Cartesian angular momentum exponent for y.
  - `lz: Vec<u32>` — Cartesian angular momentum exponent for z.
  - `shell_index: Vec<usize>` — index into the per-shell arrays for the contracted shell this
    function belongs to.
  - `atom_index: Vec<usize>` — index of the atom (in `geometry.symbols` order) this function is
    centred on.

  **Per shell** (index `s` in 0..n_shells)
  - `prim_offset: Vec<usize>` — starting index into the flat primitive arrays for shell `s`.
  - `n_primitives: Vec<usize>` — number of primitives in shell `s`.

  **Per primitive** (flat arrays, total length = Σ n_primitives[s])
  - `exponents: Vec<f64>` — primitive Gaussian exponents α.
  - `coefficients: Vec<f64>` — contraction coefficients c.

- `InitError` — error type returned by `init_basis`:
  - `BasisLoad { element: String, source: LoadError }` — `load_basis` failed for the named
    element; `source` carries the underlying `LoadError`.

---

## Gherkin Scenarios

```gherkin
Feature: Initialise the AO basis

  # Notation: a "mock BasisSet" for element X is a BasisSet value constructed directly
  # in the test (not via load_basis) to avoid network/filesystem dependencies.

  # --- Edge cases ---

  Scenario: Empty molecule returns an empty AoBasis
    Given a CartesianGeometry with no atoms
    When init_basis is called
    Then the result is Ok(AoBasis)
    And n_basis is 0
    And n_shells is 0
    And all vectors are empty
    And load_basis is never called

  # --- Basis function count ---

  Scenario: Single hydrogen atom with one s-shell gives one basis function
    Given a CartesianGeometry with one H atom at (0, 0, 0) in Bohr
    And the basis for H contains one s-shell (l=0) with 3 primitives
    When init_basis is called
    Then the result is Ok(AoBasis)
    And n_basis is 1
    And n_shells is 1

  Scenario: Single carbon atom with one s-shell and one p-shell gives four basis functions
    Given a CartesianGeometry with one C atom
    And the basis for C contains one s-shell (l=0) and one p-shell (l=1)
    When init_basis is called
    Then n_basis is 4
    And n_shells is 2

  Scenario: Single d-shell gives six basis functions
    Given a CartesianGeometry with one atom whose basis contains one d-shell (l=2)
    When init_basis is called
    Then n_basis is 6
    And n_shells is 1

  Scenario: Two-atom molecule sums basis function counts from both atoms
    Given a CartesianGeometry with atoms [O, H] where O has 2 s-shells and 1 p-shell,
      and H has 1 s-shell
    When init_basis is called
    Then n_basis is 6
    And n_shells is 4

  # --- Atom-major ordering ---

  Scenario: Basis functions for atom 0 appear before those for atom 1
    Given a CartesianGeometry with two H atoms: H0 at (0, 0, 0) and H1 at (1, 0, 0) in Bohr
    And each H has one s-shell
    When init_basis is called
    Then the basis function at index 0 has center_x = 0.0
    And the basis function at index 1 has center_x = 1.0

  Scenario: atom_index records the correct atom for each basis function
    Given a CartesianGeometry with atoms [O, H, H]
    And O has 2 s-shells and 1 p-shell, each H has 1 s-shell
    When init_basis is called
    Then atom_index[0] is 0 (O s-function)
    And atom_index[1] is 0 (O s-function from second s-shell)
    And atom_index[2] is 0
    And atom_index[3] is 0
    And atom_index[4] is 0
    And atom_index[5] is 1 (first H)
    And atom_index[6] is 2 (second H)

  # --- Shell ordering within an atom ---

  Scenario: Shells for one atom appear in basis-file order
    Given a CartesianGeometry with one atom whose basis contains shells [s, p, s] in that order
    When init_basis is called
    Then shell_index for basis functions 0, 1, 2, 3, 4 is 0, 1, 1, 1, 2 respectively

  # --- Cartesian component ordering ---

  Scenario: s-shell emits component (0, 0, 0)
    Given a CartesianGeometry with one atom whose basis has one s-shell
    When init_basis is called
    Then lx[0] = 0, ly[0] = 0, lz[0] = 0

  Scenario: p-shell emits components in order (1,0,0), (0,1,0), (0,0,1)
    Given a CartesianGeometry with one atom whose basis has one p-shell (l=1)
    When init_basis is called
    Then (lx[0], ly[0], lz[0]) = (1, 0, 0)
    And (lx[1], ly[1], lz[1]) = (0, 1, 0)
    And (lx[2], ly[2], lz[2]) = (0, 0, 1)

  Scenario: d-shell emits six components in the correct order
    Given a CartesianGeometry with one atom whose basis has one d-shell (l=2)
    When init_basis is called
    Then the six (lx, ly, lz) tuples in index order are:
      (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)

  # --- Primitive storage ---

  Scenario: Exponents and coefficients are copied from the ElectronShell
    Given a CartesianGeometry with one H atom
    And the basis for H has one s-shell with exponents [3.425, 0.624, 0.169]
      and coefficients [0.154, 0.535, 0.445]
    When init_basis is called
    Then n_primitives[0] is 3
    And prim_offset[0] is 0
    And exponents[0..3] is [3.425, 0.624, 0.169]
    And coefficients[0..3] is [0.154, 0.535, 0.445]

  Scenario: prim_offset is correct when multiple shells are present
    Given a CartesianGeometry with one atom whose basis has:
      shell 0 with 3 primitives, shell 1 with 2 primitives
    When init_basis is called
    Then prim_offset[0] is 0
    And prim_offset[1] is 3

  Scenario: All Cartesian functions of a shell share the same shell_index
    Given a CartesianGeometry with one atom whose basis has one p-shell
    When init_basis is called
    Then shell_index[0] = shell_index[1] = shell_index[2] = 0

  # --- Atom coordinates ---

  Scenario: center_x/y/z are taken from the CartesianGeometry coordinates
    Given a CartesianGeometry with one H atom at (1.5, 2.5, 3.5) in Bohr
    And the basis for H has one s-shell
    When init_basis is called
    Then center_x[0] = 1.5
    And center_y[0] = 2.5
    And center_z[0] = 3.5

  Scenario: All basis functions for the same atom share that atom's coordinates
    Given a CartesianGeometry with one atom at (1.0, 0.0, 0.0) whose basis has one p-shell
    When init_basis is called
    Then center_x[0] = center_x[1] = center_x[2] = 1.0
    And center_y[0] = center_y[1] = center_y[2] = 0.0

  # --- Repeated element types ---

  Scenario: load_basis is called once per unique element, not once per atom
    Given a CartesianGeometry with two H atoms
    And the basis for H is available
    When init_basis is called
    Then load_basis was invoked exactly once for "H"
    And both H atoms use the same shell exponents and coefficients

  # --- Error handling ---

  Scenario: InitError::BasisLoad is returned when load_basis fails for an element
    Given a CartesianGeometry with one atom of element "H"
    And load_basis("H", "unknown-basis") returns Err(LoadError::Fetch(BseError::UnknownBasisSet(_)))
    When init_basis is called with basis_name "unknown-basis"
    Then the result is Err(InitError::BasisLoad { element: "H", .. })

  Scenario: The element name in BasisLoad error identifies the failing element
    Given a CartesianGeometry with atoms [H, C]
    And load_basis succeeds for "H" but fails for "C"
    When init_basis is called
    Then the result is Err(InitError::BasisLoad { element: "C", .. })
```
