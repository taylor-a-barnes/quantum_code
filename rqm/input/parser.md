# Feature: Parse an Input File <!-- rq-203014af -->

This feature provides functions for loading and validating a simulation input file. The file is a
YAML document that specifies an initial molecular geometry and the parameters needed to run a
simulation. Supported geometry formats are Cartesian (flat coordinate array) and Z-matrix (internal
coordinates). All coordinates are stored in Bohr after parsing, being converted if necessary.

## YAML Input File Format <!-- rq-d9499ece -->

The input file is a YAML document. The four recognised top-level keys are:

| Key        | Required                   | Description                                         |
| ---------- | -------------------------- | --------------------------------------------------- |
| `driver`   | Yes                        | One of `energy`, `gradient`, `hessian`, `md`        |
| `molecule` | Yes                        | Molecular geometry block                            |
| `model`    | Yes                        | Quantum chemistry method and basis set              |
| `keywords` | Yes when `driver` is `md`  | Simulation control parameters                       |

Unknown top-level keys result in an error.

### Example — Cartesian, MD <!-- rq-87455c7f -->

```yaml
driver: md

molecule:
  symbols: [O, H, H]
  geometry: [0.0, 0.0, 0.221, 0.0, 1.431, -0.884, 0.0, -1.431, -0.884]
  units: angstrom        # optional; angstrom (default) or bohr
  charge: 0              # optional; integer, default 0
  multiplicity: 1        # optional; positive integer >= 1, default 1

model:
  method: b3lyp
  basis: sto-3g

keywords:
  timestep_fs: 0.5       # required for md; f64 > 0
  n_steps: 1000          # required for md; integer > 0
  temperature_k: 300.0   # optional; f64 >= 0, default 0.0
  thermostat: velocity_rescaling  # optional; none (default) or velocity_rescaling
```

### Example — Z-matrix, energy <!-- rq-dbdcd5b6 -->

```yaml
driver: energy

molecule:
  z_matrix:
    - symbol: O
    - symbol: H
      bond_atom: 1
      bond_length: 0.9572
    - symbol: H
      bond_atom: 1
      bond_length: 0.9572
      angle_atom: 2
      angle: 104.52
    - symbol: C
      bond_atom: 3
      bond_length: 1.5
      angle_atom: 1
      angle: 109.5
      dihedral_atom: 2
      dihedral: 120.0
  units: angstrom
  charge: 0
  multiplicity: 1

model:
  method: hf
  basis: sto-3g
```

### `molecule` Block <!-- rq-6b195950 -->

The `molecule` block must contain exactly one of:

- **Cartesian format**: both `symbols` and `geometry` keys are present.
- **Z-matrix format**: the `z_matrix` key is present.

If both sets of keys are present, the parser returns `InputError::AmbiguousGeometry`. If neither is
present, the parser returns `InputError::MissingField("molecule.geometry")`. If only one of
`symbols` / `geometry` is present without the other, the parser returns `MissingField` for the
absent key.

#### Cartesian Format

- `symbols`: YAML sequence of element symbol strings (Z = 1–118). Case-insensitive; stored in title
  case (e.g. `"o"` → `"O"`, `"FE"` → `"Fe"`).
- `geometry`: flat YAML sequence of floats in reading order `[x₁, y₁, z₁, x₂, y₂, z₂, …]`. Its
  length must equal `3 × len(symbols)`.
- `units` (optional): `angstrom` (default) or `bohr`. Coordinates are converted to Bohr on parse.
- `charge` (optional): integer, default `0`.
- `multiplicity` (optional): integer ≥ 1, default `1`.

#### Z-matrix Format

Each entry in `z_matrix` defines one atom. Atom reference indices are **1-based** and must refer to
a preceding row. The required and forbidden fields depend on the row's position:

| Row index | Required fields                                                          | Forbidden fields                                             |
| --------- | ------------------------------------------------------------------------ | ------------------------------------------------------------ |
| 0         | `symbol`                                                                 | `bond_atom`, `bond_length`, `angle_atom`, `angle`, `dihedral_atom`, `dihedral` |
| 1         | `symbol`, `bond_atom`, `bond_length`                                     | `angle_atom`, `angle`, `dihedral_atom`, `dihedral`          |
| 2         | `symbol`, `bond_atom`, `bond_length`, `angle_atom`, `angle`              | `dihedral_atom`, `dihedral`                                  |
| ≥ 3       | `symbol`, `bond_atom`, `bond_length`, `angle_atom`, `angle`, `dihedral_atom`, `dihedral` | —                                               |

Additional constraints:

- `bond_atom`, `angle_atom`, and `dihedral_atom` must be mutually distinct within the same row.
- All reference indices must be ≥ 1 and strictly less than the 1-based index of the current row
  (i.e. they must refer to a preceding row, never to the current or a future row).
- `bond_length` must be > 0.
- `angle` must satisfy 0 < angle < 180 (degrees).
- `dihedral` must satisfy −180 ≤ dihedral ≤ 180 (degrees).
- Bond lengths are converted from the molecule-level `units` to Bohr using the same factor as
  Cartesian coordinates. Angles and dihedrals are stored in degrees without conversion.
- Element symbols are validated and normalised to title case, same as in Cartesian format.
- An empty `z_matrix` sequence (zero atoms) is treated as
  `InputError::MissingField("molecule.z_matrix")`.

### `model` Block <!-- rq-6767c54c -->

Both sub-fields are required and must be non-empty strings.

- `method`: quantum chemistry method (e.g. `hf`, `b3lyp`). Stored as provided (no normalisation).
- `basis`: basis set name (e.g. `sto-3g`). Stored as provided (no normalisation).

### `keywords` Block <!-- rq-8bdee305 -->

The `keywords` block is required when `driver` is `md`. For all other drivers it is optional and its
contents are not parsed. Unknown keys within `keywords` are silently ignored.

| Field           | Required for `md` | Type   | Constraint                           | Default               |
| --------------- | ----------------- | ------ | ------------------------------------ | --------------------- |
| `timestep_fs`   | Yes               | f64    | > 0                                  | —                     |
| `n_steps`       | Yes               | usize  | > 0                                  | —                     |
| `temperature_k` | No                | f64    | ≥ 0                                  | `0.0`                 |
| `thermostat`    | No                | string | `none` or `velocity_rescaling`       | `none`                |

---

## Feature API <!-- rq-7feafa0a -->

### Functions <!-- rq-25dd2d83 -->

- `parse_input(path: &Path) -> Result<SimulationInput, InputError>` <!-- rq-379667b1 -->
  - Reads the file at `path` and delegates to `parse_input_str`.
  - Returns `InputError::IoError` if the file cannot be read.

- `parse_input_str(yaml: &str) -> Result<SimulationInput, InputError>` <!-- rq-34317111 -->
  - Parses and fully validates a YAML string.
  - Performs all validation and unit conversion described in the format section above.
  - Returns a fully validated `SimulationInput` on success.

### Types <!-- rq-3f2b44d5 -->

- `SimulationInput` <!-- rq-101b9d3d -->
  - `molecule: Molecule`
  - `model: Model`
  - `driver: Driver`
  - `keywords: Option<MdKeywords>` — `Some` when `driver` is `Md`; `None` otherwise.

- `Driver` (enum) <!-- rq-00d0ee08 -->
  - `Energy`
  - `Gradient`
  - `Hessian`
  - `Md`

- `Molecule` <!-- rq-345232c6 -->
  - `geometry: Geometry`
  - `charge: i32`
  - `multiplicity: u32`

- `Geometry` (enum) <!-- rq-40282a3b -->
  - `Cartesian(CartesianGeometry)`
  - `ZMatrix(ZMatrixGeometry)`

- `CartesianGeometry` — structure of arrays; `symbols`, `x`, `y`, `z` all have the same length <!-- rq-ba1a781f -->
  - `symbols: Vec<String>` — title-case element symbols
  - `x: Vec<f64>` — x-coordinates in Bohr
  - `y: Vec<f64>` — y-coordinates in Bohr
  - `z: Vec<f64>` — z-coordinates in Bohr

- `ZMatrixGeometry` — structure of arrays; all vectors have the same length (number of atoms) <!-- rq-6d1de1ee -->
  - `symbols: Vec<String>` — title-case element symbols
  - `bond_atoms: Vec<Option<usize>>` — 1-based indices; `None` for row 0
  - `bond_lengths_bohr: Vec<Option<f64>>` — bond lengths in Bohr; `None` for row 0
  - `angle_atoms: Vec<Option<usize>>` — `None` for rows 0–1
  - `angles_deg: Vec<Option<f64>>` — angles in degrees; `None` for rows 0–1
  - `dihedral_atoms: Vec<Option<usize>>` — `None` for rows 0–2
  - `dihedrals_deg: Vec<Option<f64>>` — dihedral angles in degrees; `None` for rows 0–2

- `Model` <!-- rq-e31d4d32 -->
  - `method: String`
  - `basis: String`

- `MdKeywords` <!-- rq-6af18856 -->
  - `timestep_fs: f64`
  - `n_steps: usize`
  - `temperature_k: f64`
  - `thermostat: Thermostat`

- `Thermostat` (enum) <!-- rq-e1ddad6c -->
  - `None`
  - `VelocityRescaling`

- `InputError` — error type returned by both functions: <!-- rq-95b0715f -->
  - `IoError(String)` — the file could not be read (`parse_input` only).
  - `InvalidYaml(String)` — the string is not valid YAML.
  - `MissingField(String)` — a required field is absent; the string identifies the dotted field path
    (e.g. `"keywords.timestep_fs"`).
  - `InvalidValue { field: String, reason: String }` — a field has an unacceptable value (wrong
    type, out-of-range, or unrecognised keyword string).
  - `AmbiguousGeometry` — both Cartesian keys (`symbols` and/or `geometry`) and `z_matrix` are
    present in the molecule block simultaneously.
  - `CoordinateMismatch { n_symbols: usize, n_coords: usize }` — the flat geometry array length
    does not equal `3 × n_symbols`.
  - `InvalidElement(String)` — an element symbol does not correspond to a known element (Z = 1–118).
  - `InvalidZMatrix { row: usize, reason: String }` — a Z-matrix row violates a structural
    constraint; `row` is 0-based.
  - `UnknownField(String)` — an unrecognised key is present at the top level of the document; the
    string names the offending key.

---

## Gherkin Scenarios <!-- rq-9164f9b1 -->

```gherkin
Feature: Parse a simulation input file

  # --- Happy paths: Cartesian ---

  @rq-a5c42421
  Scenario: Parse a minimal Cartesian energy input
    Given a YAML string with driver "energy", a molecule with symbols ["H"] and geometry [0.0, 0.0, 0.0],
      model method "hf" and basis "sto-3g", and no keywords block
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And driver is Energy
    And molecule.geometry is Cartesian with symbols ["H"]
    And molecule.charge is 0
    And molecule.multiplicity is 1
    And keywords is None

  @rq-e2b4bb49
  Scenario: Parse a Cartesian MD input with all fields
    Given a YAML string with driver "md", symbols ["O", "H", "H"],
      geometry [0.0, 0.0, 0.221, 0.0, 1.431, -0.884, 0.0, -1.431, -0.884], units "angstrom",
      charge -1, multiplicity 2, method "b3lyp", basis "sto-3g",
      timestep_fs 0.5, n_steps 1000, temperature_k 300.0, thermostat "velocity_rescaling"
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And driver is Md
    And molecule.charge is -1
    And molecule.multiplicity is 2
    And keywords is Some with timestep_fs 0.5, n_steps 1000, temperature_k 300.0,
      thermostat VelocityRescaling

  @rq-98b05632
  Scenario: Angstrom coordinates are converted to Bohr
    Given a YAML string with units "angstrom" and geometry [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      for two hydrogen atoms
    When parse_input_str is called
    Then CartesianGeometry.x[0] is approximately 1.8897259886
    And CartesianGeometry.x[1] is 0.0

  @rq-86a704d2
  Scenario: Bohr coordinates are stored unchanged
    Given a YAML string with units "bohr" and geometry [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      for two hydrogen atoms
    When parse_input_str is called
    Then CartesianGeometry.x[0] is 1.0
    And CartesianGeometry.x[1] is 0.0

  @rq-06f86585
  Scenario: Units field absent defaults to Angstrom
    Given a YAML string with no units field and geometry [1.0, 0.0, 0.0] for one hydrogen atom
    When parse_input_str is called
    Then CartesianGeometry.x[0] is approximately 1.8897259886

  @rq-e8de5a2c
  Scenario: Element symbols are normalised to title case
    Given a YAML string with symbols ["o", "H", "FE"]
    When parse_input_str is called
    Then CartesianGeometry.symbols is ["O", "H", "Fe"]

  @rq-211c6ac1
  Scenario: Driver "gradient" is accepted
    Given a YAML string with driver "gradient" and no keywords block
    When parse_input_str is called
    Then the result is Ok(SimulationInput) and driver is Gradient

  @rq-a693bf67
  Scenario: Driver "hessian" is accepted
    Given a YAML string with driver "hessian" and no keywords block
    When parse_input_str is called
    Then the result is Ok(SimulationInput) and driver is Hessian

  @rq-ce3fe4c1
  Scenario: MD keywords temperature_k and thermostat default correctly when absent
    Given a YAML string with driver "md", timestep_fs 0.5, n_steps 100, and no temperature_k
      or thermostat fields
    When parse_input_str is called
    Then keywords.temperature_k is 0.0
    And keywords.thermostat is Thermostat::None

  @rq-a0fbf0d9
  Scenario: keywords block is ignored for non-MD drivers
    Given a YAML string with driver "energy" and a keywords block containing timestep_fs 0.5
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And keywords is None

  # --- Happy paths: Z-matrix ---

  @rq-8bc6b4f9
  Scenario: Parse a Z-matrix with one atom
    Given a YAML string with z_matrix containing one entry: symbol "O"
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And geometry is ZMatrix with symbols ["O"]
    And bond_atoms[0] is None
    And bond_lengths_bohr[0] is None

  @rq-cad74df2
  Scenario: Parse a Z-matrix with two atoms
    Given a YAML string with z_matrix: [{symbol: O}, {symbol: H, bond_atom: 1, bond_length: 0.9572}]
      and units "angstrom"
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And symbols are ["O", "H"]
    And bond_atoms[1] is Some(1)
    And bond_lengths_bohr[1] is approximately 1.8897259886 * 0.9572

  @rq-8ec1d10b
  Scenario: Parse a Z-matrix with three atoms
    Given a YAML string with z_matrix containing O, H (bond to 1, length 0.96), H (bond to 1,
      length 0.96, angle_atom 2, angle 104.5)
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And angle_atoms[2] is Some(2)
    And angles_deg[2] is 104.5
    And dihedral_atoms[2] is None

  @rq-17848ed0
  Scenario: Parse a Z-matrix with four atoms including dihedral
    Given a YAML string with z_matrix containing 4 atoms where atom 3 (0-based) has
      bond_atom 1, angle_atom 2, dihedral_atom 3
    When parse_input_str is called
    Then the result is Ok(SimulationInput)
    And dihedral_atoms[3] is Some(3)
    And dihedrals_deg[3] is the specified value

  @rq-7fd5ee69
  Scenario: Z-matrix bond lengths are converted from Angstrom to Bohr
    Given a YAML string with z_matrix, units "angstrom", and a bond_length of 1.0
    When parse_input_str is called
    Then the stored bond_lengths_bohr value is approximately 1.8897259886

  # --- File I/O ---

  @rq-47e399ea
  Scenario: parse_input reads and parses a valid file
    Given a file on disk containing a valid simulation input YAML
    When parse_input is called with the file's path
    Then the result is Ok(SimulationInput)

  @rq-39f344a4
  Scenario: parse_input returns IoError for a missing file
    Given no file exists at the given path
    When parse_input is called with that path
    Then parse_input returns Err(InputError::IoError(_))

  # --- YAML errors ---

  @rq-a07bed47
  Scenario: Invalid YAML returns InvalidYaml
    Given the string "driver: md\nmolecule: :\n  bad:"
    When parse_input_str is called
    Then the result is Err(InputError::InvalidYaml(_))

  # --- Missing required top-level fields ---

  @rq-756fd335
  Scenario: Missing driver returns MissingField
    Given a valid YAML string with the driver key omitted
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("driver"))

  @rq-94e1b8d6
  Scenario: Missing molecule returns MissingField
    Given a valid YAML string with the molecule key omitted
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("molecule"))

  @rq-3473d09f
  Scenario: Missing model returns MissingField
    Given a valid YAML string with the model key omitted
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("model"))

  @rq-fb41de1a
  Scenario: Missing model.method returns MissingField
    Given a YAML string where model is present but method is absent
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("model.method"))

  @rq-d84f068c
  Scenario: Missing model.basis returns MissingField
    Given a YAML string where model is present but basis is absent
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("model.basis"))

  @rq-36b8bb6e
  Scenario: Missing keywords block when driver is md returns MissingField
    Given a YAML string with driver "md" and no keywords block
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("keywords"))

  @rq-989109f7
  Scenario: Missing keywords.timestep_fs when driver is md returns MissingField
    Given a YAML string with driver "md" and a keywords block that omits timestep_fs
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("keywords.timestep_fs"))

  @rq-bda1df68
  Scenario: Missing keywords.n_steps when driver is md returns MissingField
    Given a YAML string with driver "md" and a keywords block that omits n_steps
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("keywords.n_steps"))

  @rq-d35d5c25
  Scenario: Missing molecule.geometry when only symbols is present
    Given a YAML string where molecule contains symbols but no geometry key and no z_matrix key
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("molecule.geometry"))

  @rq-db979c67
  Scenario: Missing molecule.symbols when only geometry is present
    Given a YAML string where molecule contains geometry but no symbols key and no z_matrix key
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("molecule.symbols"))

  @rq-e21f8b10
  Scenario: Missing geometry when neither Cartesian nor z_matrix keys are present
    Given a YAML string where molecule contains only charge and multiplicity
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("molecule.geometry"))

  @rq-76843345
  Scenario: Empty z_matrix sequence returns MissingField
    Given a YAML string where z_matrix is an empty list []
    When parse_input_str is called
    Then the result is Err(InputError::MissingField("molecule.z_matrix"))

  # --- Unknown fields ---

  @rq-e2c517c7
  Scenario: Unknown top-level key returns UnknownField
    Given a YAML string that is otherwise valid but contains an extra top-level key "extra_key"
    When parse_input_str is called
    Then the result is Err(InputError::UnknownField("extra_key"))

  # --- Invalid values ---

  @rq-d0efc5b0
  Scenario: Unrecognised driver string returns InvalidValue
    Given a YAML string with driver "optimize"
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "driver", .. })

  @rq-1c40ba26
  Scenario: Unrecognised units string returns InvalidValue
    Given a YAML string with units "nanometer"
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "molecule.units", .. })

  @rq-4d09f400
  Scenario: Empty model.method string returns InvalidValue
    Given a YAML string where model.method is ""
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "model.method", .. })

  @rq-5b73711a
  Scenario: Empty model.basis string returns InvalidValue
    Given a YAML string where model.basis is ""
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "model.basis", .. })

  @rq-a07acae9
  Scenario: timestep_fs of zero returns InvalidValue
    Given a YAML string with driver "md" and keywords.timestep_fs set to 0.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "keywords.timestep_fs", .. })

  @rq-9fa5aa7a
  Scenario: Negative timestep_fs returns InvalidValue
    Given a YAML string with driver "md" and keywords.timestep_fs set to -1.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "keywords.timestep_fs", .. })

  @rq-a2e288de
  Scenario: n_steps of zero returns InvalidValue
    Given a YAML string with driver "md" and keywords.n_steps set to 0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "keywords.n_steps", .. })

  @rq-35f7fca9
  Scenario: Negative temperature_k returns InvalidValue
    Given a YAML string with driver "md" and keywords.temperature_k set to -1.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "keywords.temperature_k", .. })

  @rq-179c77ad
  Scenario: Unrecognised thermostat string returns InvalidValue
    Given a YAML string with driver "md" and keywords.thermostat set to "langevin"
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "keywords.thermostat", .. })

  @rq-916925f6
  Scenario: multiplicity of zero returns InvalidValue
    Given a YAML string where molecule.multiplicity is 0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidValue { field: "molecule.multiplicity", .. })

  # --- Geometry format errors ---

  @rq-f30173a0
  Scenario: Both Cartesian keys and z_matrix present returns AmbiguousGeometry
    Given a YAML string where molecule contains symbols, geometry, and z_matrix simultaneously
    When parse_input_str is called
    Then the result is Err(InputError::AmbiguousGeometry)

  @rq-4a594aa0
  Scenario: Geometry array length not divisible by 3 per atom returns CoordinateMismatch
    Given a YAML string with 2 symbols and a geometry array of length 5
    When parse_input_str is called
    Then the result is Err(InputError::CoordinateMismatch { n_symbols: 2, n_coords: 5 })

  @rq-3df263f0
  Scenario: Geometry array too short returns CoordinateMismatch
    Given a YAML string with 3 symbols and a geometry array of length 6
    When parse_input_str is called
    Then the result is Err(InputError::CoordinateMismatch { n_symbols: 3, n_coords: 6 })

  # --- Element validation ---

  @rq-08af8a53
  Scenario: Unknown element symbol in Cartesian symbols returns InvalidElement
    Given a YAML string with symbols ["H", "Xx"]
    When parse_input_str is called
    Then the result is Err(InputError::InvalidElement("Xx"))

  @rq-8e45c18b
  Scenario: Unknown element symbol in z_matrix returns InvalidElement
    Given a YAML string with z_matrix containing an entry with symbol "Zz"
    When parse_input_str is called
    Then the result is Err(InputError::InvalidElement("Zz"))

  # --- Z-matrix structural errors ---

  @rq-6b64811f
  Scenario: Row 0 with bond_atom present returns InvalidZMatrix for row 0
    Given a YAML string where z_matrix row 0 specifies bond_atom: 1
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 0, .. })

  @rq-1969d663
  Scenario: Row 1 missing bond_atom returns InvalidZMatrix for row 1
    Given a YAML string where z_matrix row 1 has bond_length but no bond_atom
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-75635399
  Scenario: Row 1 with angle_atom present returns InvalidZMatrix for row 1
    Given a YAML string where z_matrix row 1 specifies angle_atom and angle in addition to
      bond_atom and bond_length
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-6ea890dc
  Scenario: Row 2 missing angle returns InvalidZMatrix for row 2
    Given a YAML string where z_matrix row 2 has bond fields and angle_atom but no angle
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 2, .. })

  @rq-599edb01
  Scenario: Row 2 with dihedral_atom present returns InvalidZMatrix for row 2
    Given a YAML string where z_matrix row 2 specifies dihedral_atom and dihedral
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 2, .. })

  @rq-7f3bcecd
  Scenario: Row 3 missing dihedral_atom returns InvalidZMatrix for row 3
    Given a YAML string where z_matrix row 3 (0-based) omits the dihedral_atom field
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 3, .. })

  @rq-3246595c
  Scenario: Reference index 0 (out of 1-based range) returns InvalidZMatrix
    Given a YAML string where z_matrix row 1 has bond_atom: 0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-ed57fd1b
  Scenario: Forward reference returns InvalidZMatrix
    Given a YAML string where z_matrix row 1 has bond_atom: 2 (referring to a future row)
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-4fd55f64
  Scenario: Duplicate reference indices within a row returns InvalidZMatrix
    Given a YAML string where z_matrix row 3 has bond_atom and angle_atom both set to 1
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 3, .. })

  @rq-07ff5f51
  Scenario: bond_length of zero returns InvalidZMatrix
    Given a YAML string where z_matrix row 1 has bond_length: 0.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-adc7637f
  Scenario: Negative bond_length returns InvalidZMatrix
    Given a YAML string where z_matrix row 1 has bond_length: -1.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 1, .. })

  @rq-c03e19e3
  Scenario: angle of zero degrees returns InvalidZMatrix
    Given a YAML string where z_matrix row 2 has angle: 0.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 2, .. })

  @rq-a400431d
  Scenario: angle of 180 degrees returns InvalidZMatrix
    Given a YAML string where z_matrix row 2 has angle: 180.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 2, .. })

  @rq-80e1dfda
  Scenario: dihedral outside [-180, 180] returns InvalidZMatrix
    Given a YAML string where z_matrix row 3 has dihedral: 181.0
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 3, .. })

  @rq-e9eadc08
  Scenario: InvalidZMatrix row index identifies the correct row
    Given a YAML string with a 5-atom z_matrix where only row 4 (0-based) has an invalid
      bond_length of -1.0 and all prior rows are valid
    When parse_input_str is called
    Then the result is Err(InputError::InvalidZMatrix { row: 4, .. })
```
