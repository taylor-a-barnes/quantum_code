use std::path::Path;

const ANGSTROM_TO_BOHR: f64 = 1.8897259886;

/// All 118 known element symbols in title case, indexed by atomic number (1-based).
const ELEMENTS: &[&str] = &[
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
  "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
  "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

// ── Error type ────────────────────────────────────────────────────────────────

#[derive(Debug, PartialEq)]
pub enum InputError {
  IoError(String),
  InvalidYaml(String),
  MissingField(String),
  InvalidValue { field: String, reason: String },
  AmbiguousGeometry,
  CoordinateMismatch { n_symbols: usize, n_coords: usize },
  InvalidElement(String),
  InvalidZMatrix { row: usize, reason: String },
  UnknownField(String),
}

impl std::fmt::Display for InputError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      InputError::IoError(s) =>
        write!(f, "I/O error: {}", s),
      InputError::InvalidYaml(s) =>
        write!(f, "invalid YAML: {}", s),
      InputError::MissingField(s) =>
        write!(f, "missing required field: {}", s),
      InputError::InvalidValue { field, reason } =>
        write!(f, "invalid value for {}: {}", field, reason),
      InputError::AmbiguousGeometry =>
        write!(f, "molecule block contains both Cartesian and Z-matrix keys"),
      InputError::CoordinateMismatch { n_symbols, n_coords } =>
        write!(f, "geometry has {} coordinates but expected {} (3 × {})",
          n_coords, 3 * n_symbols, n_symbols),
      InputError::InvalidElement(s) =>
        write!(f, "unknown element symbol: {:?}", s),
      InputError::InvalidZMatrix { row, reason } =>
        write!(f, "invalid z_matrix row {}: {}", row, reason),
      InputError::UnknownField(s) =>
        write!(f, "unknown top-level field: {:?}", s),
    }
  }
}

// ── Public types ──────────────────────────────────────────────────────────────

#[derive(Debug, PartialEq)]
pub enum Driver {
  Energy,
  Gradient,
  Hessian,
  Md,
}

#[derive(Debug, PartialEq)]
pub struct Molecule {
  pub geometry: Geometry,
  pub charge: i32,
  pub multiplicity: u32,
}

#[derive(Debug, PartialEq)]
pub enum Geometry {
  Cartesian(CartesianGeometry),
  ZMatrix(ZMatrixGeometry),
}

/// Structure of arrays; all coordinate vectors have the same length as `symbols`.
#[derive(Debug, PartialEq)]
pub struct CartesianGeometry {
  pub symbols: Vec<String>,
  pub x: Vec<f64>,
  pub y: Vec<f64>,
  pub z: Vec<f64>,
}

/// Structure of arrays; all vectors have the same length (number of atoms).
#[derive(Debug, PartialEq)]
pub struct ZMatrixGeometry {
  pub symbols: Vec<String>,
  /// 1-based reference indices; `None` for row 0.
  pub bond_atoms: Vec<Option<usize>>,
  /// Bond lengths in Bohr; `None` for row 0.
  pub bond_lengths_bohr: Vec<Option<f64>>,
  /// `None` for rows 0–1.
  pub angle_atoms: Vec<Option<usize>>,
  /// Angles in degrees; `None` for rows 0–1.
  pub angles_deg: Vec<Option<f64>>,
  /// `None` for rows 0–2.
  pub dihedral_atoms: Vec<Option<usize>>,
  /// Dihedral angles in degrees; `None` for rows 0–2.
  pub dihedrals_deg: Vec<Option<f64>>,
}

#[derive(Debug, PartialEq)]
pub struct Model {
  pub method: String,
  pub basis: String,
}

#[derive(Debug, PartialEq)]
pub struct MdKeywords {
  pub timestep_fs: f64,
  pub n_steps: usize,
  pub temperature_k: f64,
  pub thermostat: Thermostat,
}

#[derive(Debug, PartialEq)]
pub enum Thermostat {
  None,
  VelocityRescaling,
}

#[derive(Debug, PartialEq)]
pub struct SimulationInput {
  pub molecule: Molecule,
  pub model: Model,
  pub driver: Driver,
  /// `Some` when `driver` is `Md`; `None` otherwise.
  pub keywords: Option<MdKeywords>,
}

// ── Public functions ──────────────────────────────────────────────────────────

/// Reads the file at `path` and delegates to [`parse_input_str`].
pub fn parse_input(path: &Path) -> Result<SimulationInput, InputError> {
  let yaml = std::fs::read_to_string(path)
    .map_err(|e| InputError::IoError(e.to_string()))?;
  parse_input_str(&yaml)
}

/// Parses and fully validates a YAML simulation input string.
pub fn parse_input_str(yaml: &str) -> Result<SimulationInput, InputError> {
  let value: serde_yaml::Value = serde_yaml::from_str(yaml)
    .map_err(|e| InputError::InvalidYaml(e.to_string()))?;

  let mapping = match &value {
    serde_yaml::Value::Mapping(m) => m,
    _ => return Err(InputError::InvalidYaml(
      "expected a mapping at top level".to_string()
    )),
  };

  // Reject unknown top-level keys.
  const KNOWN: &[&str] = &["driver", "molecule", "model", "keywords"];
  for (key, _) in mapping.iter() {
    let k = key.as_str().unwrap_or("");
    if !KNOWN.contains(&k) {
      return Err(InputError::UnknownField(k.to_string()));
    }
  }

  let driver = parse_driver(mapping)?;
  let molecule = parse_molecule(mapping)?;
  let model = parse_model(mapping)?;

  let keywords = if driver == Driver::Md {
    let kw_val = map_get(mapping, "keywords")
      .ok_or_else(|| InputError::MissingField("keywords".to_string()))?;
    Some(parse_keywords(kw_val)?)
  } else {
    Option::None
  };

  Ok(SimulationInput { molecule, model, driver, keywords })
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Look up a string key in a YAML mapping.
fn map_get<'a>(
  m: &'a serde_yaml::Mapping,
  key: &str,
) -> Option<&'a serde_yaml::Value> {
  m.iter()
   .find(|(k, _)| k.as_str() == Some(key))
   .map(|(_, v)| v)
}

/// Normalise `sym` to title case and validate against the periodic table.
/// Returns the normalised symbol on success.
fn normalize_element(sym: &str) -> Result<String, InputError> {
  if sym.is_empty() {
    return Err(InputError::InvalidElement(sym.to_string()));
  }
  let mut chars = sym.chars();
  let normalized = match chars.next() {
    Some(c) => {
      let upper: String = c.to_uppercase().collect();
      upper + &chars.as_str().to_lowercase()
    }
    Option::None => return Err(InputError::InvalidElement(sym.to_string())),
  };
  if ELEMENTS.contains(&normalized.as_str()) {
    Ok(normalized)
  } else {
    Err(InputError::InvalidElement(sym.to_string()))
  }
}

fn parse_driver(m: &serde_yaml::Mapping) -> Result<Driver, InputError> {
  let v = map_get(m, "driver")
    .ok_or_else(|| InputError::MissingField("driver".to_string()))?;
  let s = v.as_str().ok_or_else(|| InputError::InvalidValue {
    field: "driver".to_string(),
    reason: "expected a string".to_string(),
  })?;
  match s {
    "energy" => Ok(Driver::Energy),
    "gradient" => Ok(Driver::Gradient),
    "hessian" => Ok(Driver::Hessian),
    "md" => Ok(Driver::Md),
    other => Err(InputError::InvalidValue {
      field: "driver".to_string(),
      reason: format!("unrecognised driver {:?}", other),
    }),
  }
}

fn parse_molecule(m: &serde_yaml::Mapping) -> Result<Molecule, InputError> {
  let v = map_get(m, "molecule")
    .ok_or_else(|| InputError::MissingField("molecule".to_string()))?;
  let mol_map = v.as_mapping().ok_or_else(|| InputError::InvalidValue {
    field: "molecule".to_string(),
    reason: "expected a mapping".to_string(),
  })?;

  let charge = if let Some(cv) = map_get(mol_map, "charge") {
    cv.as_i64().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.charge".to_string(),
      reason: "expected an integer".to_string(),
    })? as i32
  } else {
    0
  };

  let multiplicity = if let Some(mv) = map_get(mol_map, "multiplicity") {
    let n = mv.as_i64().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.multiplicity".to_string(),
      reason: "expected an integer".to_string(),
    })?;
    if n < 1 {
      return Err(InputError::InvalidValue {
        field: "molecule.multiplicity".to_string(),
        reason: format!("must be >= 1, got {}", n),
      });
    }
    n as u32
  } else {
    1
  };

  let units_factor = parse_units(mol_map)?;

  let has_symbols = map_get(mol_map, "symbols").is_some();
  let has_geometry = map_get(mol_map, "geometry").is_some();
  let has_zmatrix = map_get(mol_map, "z_matrix").is_some();

  if has_zmatrix && (has_symbols || has_geometry) {
    return Err(InputError::AmbiguousGeometry);
  }

  let geometry = if has_zmatrix {
    Geometry::ZMatrix(parse_zmatrix(mol_map, units_factor)?)
  } else if has_symbols && has_geometry {
    Geometry::Cartesian(parse_cartesian(mol_map, units_factor)?)
  } else if has_symbols {
    return Err(InputError::MissingField("molecule.geometry".to_string()));
  } else if has_geometry {
    return Err(InputError::MissingField("molecule.symbols".to_string()));
  } else {
    return Err(InputError::MissingField("molecule.geometry".to_string()));
  };

  Ok(Molecule { geometry, charge, multiplicity })
}

/// Returns the coordinate conversion factor (raw → Bohr).
fn parse_units(mol_map: &serde_yaml::Mapping) -> Result<f64, InputError> {
  if let Some(u) = map_get(mol_map, "units") {
    let s = u.as_str().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.units".to_string(),
      reason: "expected a string".to_string(),
    })?;
    match s {
      "angstrom" => Ok(ANGSTROM_TO_BOHR),
      "bohr" => Ok(1.0),
      other => Err(InputError::InvalidValue {
        field: "molecule.units".to_string(),
        reason: format!("unrecognised units {:?}", other),
      }),
    }
  } else {
    Ok(ANGSTROM_TO_BOHR)
  }
}

fn parse_cartesian(
  mol_map: &serde_yaml::Mapping,
  factor: f64,
) -> Result<CartesianGeometry, InputError> {
  let sym_seq = map_get(mol_map, "symbols")
    .unwrap()
    .as_sequence()
    .ok_or_else(|| InputError::InvalidValue {
      field: "molecule.symbols".to_string(),
      reason: "expected a sequence".to_string(),
    })?;

  let mut symbols = Vec::with_capacity(sym_seq.len());
  for sv in sym_seq {
    let s = sv.as_str().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.symbols".to_string(),
      reason: "element symbols must be strings".to_string(),
    })?;
    symbols.push(normalize_element(s)?);
  }

  let geo_seq = map_get(mol_map, "geometry")
    .unwrap()
    .as_sequence()
    .ok_or_else(|| InputError::InvalidValue {
      field: "molecule.geometry".to_string(),
      reason: "expected a sequence".to_string(),
    })?;

  let n_symbols = symbols.len();
  let n_coords = geo_seq.len();
  if n_coords != 3 * n_symbols {
    return Err(InputError::CoordinateMismatch { n_symbols, n_coords });
  }

  let mut x = Vec::with_capacity(n_symbols);
  let mut y = Vec::with_capacity(n_symbols);
  let mut z = Vec::with_capacity(n_symbols);

  for chunk in geo_seq.chunks(3) {
    let cx = chunk[0].as_f64().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.geometry".to_string(),
      reason: "coordinates must be numbers".to_string(),
    })? * factor;
    let cy = chunk[1].as_f64().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.geometry".to_string(),
      reason: "coordinates must be numbers".to_string(),
    })? * factor;
    let cz = chunk[2].as_f64().ok_or_else(|| InputError::InvalidValue {
      field: "molecule.geometry".to_string(),
      reason: "coordinates must be numbers".to_string(),
    })? * factor;
    x.push(cx);
    y.push(cy);
    z.push(cz);
  }

  Ok(CartesianGeometry { symbols, x, y, z })
}

fn parse_zmatrix(
  mol_map: &serde_yaml::Mapping,
  factor: f64,
) -> Result<ZMatrixGeometry, InputError> {
  let rows = map_get(mol_map, "z_matrix")
    .unwrap()
    .as_sequence()
    .ok_or_else(|| InputError::InvalidValue {
      field: "molecule.z_matrix".to_string(),
      reason: "expected a sequence".to_string(),
    })?;

  if rows.is_empty() {
    return Err(InputError::MissingField("molecule.z_matrix".to_string()));
  }

  let n = rows.len();
  let mut symbols = Vec::with_capacity(n);
  let mut bond_atoms = Vec::with_capacity(n);
  let mut bond_lengths_bohr = Vec::with_capacity(n);
  let mut angle_atoms = Vec::with_capacity(n);
  let mut angles_deg = Vec::with_capacity(n);
  let mut dihedral_atoms = Vec::with_capacity(n);
  let mut dihedrals_deg = Vec::with_capacity(n);

  for (i, row_val) in rows.iter().enumerate() {
    let row_map = row_val.as_mapping().ok_or_else(|| izm(i,
      "each z_matrix entry must be a mapping"))?;

    let sym_str = map_get(row_map, "symbol")
      .ok_or_else(|| izm(i, "missing required field 'symbol'"))?
      .as_str()
      .ok_or_else(|| izm(i, "'symbol' must be a string"))?;

    let sym = normalize_element(sym_str)
      .map_err(|_| InputError::InvalidElement(sym_str.to_string()))?;
    symbols.push(sym);

    let has_bond_atom    = map_get(row_map, "bond_atom").is_some();
    let has_bond_length  = map_get(row_map, "bond_length").is_some();
    let has_angle_atom   = map_get(row_map, "angle_atom").is_some();
    let has_angle        = map_get(row_map, "angle").is_some();
    let has_dihedral_atom = map_get(row_map, "dihedral_atom").is_some();
    let has_dihedral     = map_get(row_map, "dihedral").is_some();

    match i {
      0 => {
        if has_bond_atom || has_bond_length || has_angle_atom
          || has_angle || has_dihedral_atom || has_dihedral
        {
          return Err(izm(0, "row 0 must only contain 'symbol'"));
        }
        bond_atoms.push(Option::None);
        bond_lengths_bohr.push(Option::None);
        angle_atoms.push(Option::None);
        angles_deg.push(Option::None);
        dihedral_atoms.push(Option::None);
        dihedrals_deg.push(Option::None);
      }
      1 => {
        if has_angle_atom || has_angle || has_dihedral_atom || has_dihedral {
          return Err(izm(1, "row 1 must not contain angle or dihedral fields"));
        }
        if !has_bond_atom { return Err(izm(1, "missing required field 'bond_atom'")); }
        if !has_bond_length { return Err(izm(1, "missing required field 'bond_length'")); }
        let ba = zmat_ref_idx(row_map, "bond_atom", 1)?;
        let bl = zmat_bond_length(row_map, 1, factor)?;
        bond_atoms.push(Some(ba));
        bond_lengths_bohr.push(Some(bl));
        angle_atoms.push(Option::None);
        angles_deg.push(Option::None);
        dihedral_atoms.push(Option::None);
        dihedrals_deg.push(Option::None);
      }
      2 => {
        if has_dihedral_atom || has_dihedral {
          return Err(izm(2, "row 2 must not contain dihedral fields"));
        }
        if !has_bond_atom   { return Err(izm(2, "missing required field 'bond_atom'")); }
        if !has_bond_length { return Err(izm(2, "missing required field 'bond_length'")); }
        if !has_angle_atom  { return Err(izm(2, "missing required field 'angle_atom'")); }
        if !has_angle       { return Err(izm(2, "missing required field 'angle'")); }
        let ba  = zmat_ref_idx(row_map, "bond_atom", 2)?;
        let bl  = zmat_bond_length(row_map, 2, factor)?;
        let aa  = zmat_ref_idx(row_map, "angle_atom", 2)?;
        let ang = zmat_angle(row_map, 2)?;
        check_distinct(2, ba, Some(aa), Option::None)?;
        bond_atoms.push(Some(ba));
        bond_lengths_bohr.push(Some(bl));
        angle_atoms.push(Some(aa));
        angles_deg.push(Some(ang));
        dihedral_atoms.push(Option::None);
        dihedrals_deg.push(Option::None);
      }
      _ => {
        if !has_bond_atom     { return Err(izm(i, "missing required field 'bond_atom'")); }
        if !has_bond_length   { return Err(izm(i, "missing required field 'bond_length'")); }
        if !has_angle_atom    { return Err(izm(i, "missing required field 'angle_atom'")); }
        if !has_angle         { return Err(izm(i, "missing required field 'angle'")); }
        if !has_dihedral_atom { return Err(izm(i, "missing required field 'dihedral_atom'")); }
        if !has_dihedral      { return Err(izm(i, "missing required field 'dihedral'")); }
        let ba  = zmat_ref_idx(row_map, "bond_atom", i)?;
        let bl  = zmat_bond_length(row_map, i, factor)?;
        let aa  = zmat_ref_idx(row_map, "angle_atom", i)?;
        let ang = zmat_angle(row_map, i)?;
        let da  = zmat_ref_idx(row_map, "dihedral_atom", i)?;
        let dih = zmat_dihedral(row_map, i)?;
        check_distinct(i, ba, Some(aa), Some(da))?;
        bond_atoms.push(Some(ba));
        bond_lengths_bohr.push(Some(bl));
        angle_atoms.push(Some(aa));
        angles_deg.push(Some(ang));
        dihedral_atoms.push(Some(da));
        dihedrals_deg.push(Some(dih));
      }
    }
  }

  Ok(ZMatrixGeometry {
    symbols,
    bond_atoms,
    bond_lengths_bohr,
    angle_atoms,
    angles_deg,
    dihedral_atoms,
    dihedrals_deg,
  })
}

/// Parse a 1-based Z-matrix reference index that must refer to a preceding row.
fn zmat_ref_idx(
  row_map: &serde_yaml::Mapping,
  field: &str,
  row: usize,
) -> Result<usize, InputError> {
  let v = map_get(row_map, field).unwrap();
  let idx_i = v.as_i64().ok_or_else(|| izm(row,
    &format!("'{}' must be a positive integer", field)))?;
  // Valid range: 1-based index strictly less than 1-based current row, i.e. 1..=row.
  if idx_i <= 0 || idx_i as usize > row {
    return Err(izm(row, &format!(
      "'{}' = {} is out of range; must be 1 to {}",
      field, idx_i, row
    )));
  }
  Ok(idx_i as usize)
}

fn zmat_bond_length(
  row_map: &serde_yaml::Mapping,
  row: usize,
  factor: f64,
) -> Result<f64, InputError> {
  let v = map_get(row_map, "bond_length").unwrap();
  let bl = v.as_f64()
    .ok_or_else(|| izm(row, "'bond_length' must be a number"))?;
  if bl <= 0.0 {
    return Err(izm(row, &format!("'bond_length' must be > 0, got {}", bl)));
  }
  Ok(bl * factor)
}

fn zmat_angle(row_map: &serde_yaml::Mapping, row: usize) -> Result<f64, InputError> {
  let v = map_get(row_map, "angle").unwrap();
  let a = v.as_f64()
    .ok_or_else(|| izm(row, "'angle' must be a number"))?;
  if a <= 0.0 || a >= 180.0 {
    return Err(izm(row, &format!("'angle' must satisfy 0 < angle < 180, got {}", a)));
  }
  Ok(a)
}

fn zmat_dihedral(row_map: &serde_yaml::Mapping, row: usize) -> Result<f64, InputError> {
  let v = map_get(row_map, "dihedral").unwrap();
  let d = v.as_f64()
    .ok_or_else(|| izm(row, "'dihedral' must be a number"))?;
  if d < -180.0 || d > 180.0 {
    return Err(izm(row, &format!(
      "'dihedral' must satisfy -180 <= dihedral <= 180, got {}", d
    )));
  }
  Ok(d)
}

/// Verify that bond_atom, angle_atom, and dihedral_atom are mutually distinct.
fn check_distinct(
  row: usize,
  ba: usize,
  aa: Option<usize>,
  da: Option<usize>,
) -> Result<(), InputError> {
  if let Some(a) = aa {
    if ba == a {
      return Err(izm(row, &format!(
        "bond_atom ({}) and angle_atom ({}) must be distinct", ba, a
      )));
    }
    if let Some(d) = da {
      if ba == d {
        return Err(izm(row, &format!(
          "bond_atom ({}) and dihedral_atom ({}) must be distinct", ba, d
        )));
      }
      if a == d {
        return Err(izm(row, &format!(
          "angle_atom ({}) and dihedral_atom ({}) must be distinct", a, d
        )));
      }
    }
  }
  Ok(())
}

/// Convenience constructor for `InvalidZMatrix`.
fn izm(row: usize, reason: &str) -> InputError {
  InputError::InvalidZMatrix { row, reason: reason.to_string() }
}

fn parse_model(m: &serde_yaml::Mapping) -> Result<Model, InputError> {
  let v = map_get(m, "model")
    .ok_or_else(|| InputError::MissingField("model".to_string()))?;
  let model_map = v.as_mapping().ok_or_else(|| InputError::InvalidValue {
    field: "model".to_string(),
    reason: "expected a mapping".to_string(),
  })?;

  let method = map_get(model_map, "method")
    .ok_or_else(|| InputError::MissingField("model.method".to_string()))?
    .as_str()
    .ok_or_else(|| InputError::InvalidValue {
      field: "model.method".to_string(),
      reason: "expected a string".to_string(),
    })?
    .to_string();
  if method.is_empty() {
    return Err(InputError::InvalidValue {
      field: "model.method".to_string(),
      reason: "must not be empty".to_string(),
    });
  }

  let basis = map_get(model_map, "basis")
    .ok_or_else(|| InputError::MissingField("model.basis".to_string()))?
    .as_str()
    .ok_or_else(|| InputError::InvalidValue {
      field: "model.basis".to_string(),
      reason: "expected a string".to_string(),
    })?
    .to_string();
  if basis.is_empty() {
    return Err(InputError::InvalidValue {
      field: "model.basis".to_string(),
      reason: "must not be empty".to_string(),
    });
  }

  Ok(Model { method, basis })
}

fn parse_keywords(v: &serde_yaml::Value) -> Result<MdKeywords, InputError> {
  let kw_map = v.as_mapping().ok_or_else(|| InputError::InvalidValue {
    field: "keywords".to_string(),
    reason: "expected a mapping".to_string(),
  })?;

  let timestep_fs = {
    let tv = map_get(kw_map, "timestep_fs")
      .ok_or_else(|| InputError::MissingField("keywords.timestep_fs".to_string()))?;
    let t = tv.as_f64().ok_or_else(|| InputError::InvalidValue {
      field: "keywords.timestep_fs".to_string(),
      reason: "expected a number".to_string(),
    })?;
    if t <= 0.0 {
      return Err(InputError::InvalidValue {
        field: "keywords.timestep_fs".to_string(),
        reason: format!("must be > 0, got {}", t),
      });
    }
    t
  };

  let n_steps = {
    let nv = map_get(kw_map, "n_steps")
      .ok_or_else(|| InputError::MissingField("keywords.n_steps".to_string()))?;
    let n = nv.as_i64().ok_or_else(|| InputError::InvalidValue {
      field: "keywords.n_steps".to_string(),
      reason: "expected an integer".to_string(),
    })?;
    if n <= 0 {
      return Err(InputError::InvalidValue {
        field: "keywords.n_steps".to_string(),
        reason: format!("must be > 0, got {}", n),
      });
    }
    n as usize
  };

  let temperature_k = if let Some(tv) = map_get(kw_map, "temperature_k") {
    let t = tv.as_f64().ok_or_else(|| InputError::InvalidValue {
      field: "keywords.temperature_k".to_string(),
      reason: "expected a number".to_string(),
    })?;
    if t < 0.0 {
      return Err(InputError::InvalidValue {
        field: "keywords.temperature_k".to_string(),
        reason: format!("must be >= 0, got {}", t),
      });
    }
    t
  } else {
    0.0
  };

  let thermostat = if let Some(tv) = map_get(kw_map, "thermostat") {
    let s = tv.as_str().ok_or_else(|| InputError::InvalidValue {
      field: "keywords.thermostat".to_string(),
      reason: "expected a string".to_string(),
    })?;
    match s {
      "none" => Thermostat::None,
      "velocity_rescaling" => Thermostat::VelocityRescaling,
      other => return Err(InputError::InvalidValue {
        field: "keywords.thermostat".to_string(),
        reason: format!("unrecognised thermostat {:?}", other),
      }),
    }
  } else {
    Thermostat::None
  };

  Ok(MdKeywords { timestep_fs, n_steps, temperature_k, thermostat })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
  use super::*;
  use std::io::Write;

  const BOHR: f64 = ANGSTROM_TO_BOHR;

  fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
  }

  // ── YAML helpers ────────────────────────────────────────────────────────────

  fn energy_yaml() -> &'static str {
    "driver: energy\n\
     molecule:\n\
     \x20 symbols: [H]\n\
     \x20 geometry: [0.0, 0.0, 0.0]\n\
     model:\n\
     \x20 method: hf\n\
     \x20 basis: sto-3g\n"
  }

  fn md_yaml() -> &'static str {
    "driver: md\n\
     molecule:\n\
     \x20 symbols: [H]\n\
     \x20 geometry: [0.0, 0.0, 0.0]\n\
     model:\n\
     \x20 method: hf\n\
     \x20 basis: sto-3g\n\
     keywords:\n\
     \x20 timestep_fs: 0.5\n\
     \x20 n_steps: 100\n"
  }

  fn zmat_energy_yaml() -> &'static str {
    "driver: energy\n\
     molecule:\n\
     \x20 z_matrix:\n\
     \x20   - symbol: O\n\
     \x20   - symbol: H\n\
     \x20     bond_atom: 1\n\
     \x20     bond_length: 0.96\n\
     \x20   - symbol: H\n\
     \x20     bond_atom: 1\n\
     \x20     bond_length: 0.96\n\
     \x20     angle_atom: 2\n\
     \x20     angle: 104.5\n\
     \x20   - symbol: C\n\
     \x20     bond_atom: 1\n\
     \x20     bond_length: 1.5\n\
     \x20     angle_atom: 2\n\
     \x20     angle: 109.5\n\
     \x20     dihedral_atom: 3\n\
     \x20     dihedral: 120.0\n\
     model:\n\
     \x20 method: hf\n\
     \x20 basis: sto-3g\n"
  }

  // ── Happy paths: Cartesian ──────────────────────────────────────────────────

  #[test]
  fn test_minimal_cartesian_energy() {
    let r = parse_input_str(energy_yaml()).unwrap();
    assert_eq!(r.driver, Driver::Energy);
    assert_eq!(r.keywords, Option::None);
    assert_eq!(r.molecule.charge, 0);
    assert_eq!(r.molecule.multiplicity, 1);
    match &r.molecule.geometry {
      Geometry::Cartesian(c) => assert_eq!(c.symbols, vec!["H"]),
      _ => panic!("expected Cartesian"),
    }
  }

  #[test]
  fn test_full_cartesian_md() {
    let yaml = "driver: md\n\
      molecule:\n\
      \x20 symbols: [O, H, H]\n\
      \x20 geometry: [0.0, 0.0, 0.221, 0.0, 1.431, -0.884, 0.0, -1.431, -0.884]\n\
      \x20 units: angstrom\n\
      \x20 charge: -1\n\
      \x20 multiplicity: 2\n\
      model:\n\
      \x20 method: b3lyp\n\
      \x20 basis: sto-3g\n\
      keywords:\n\
      \x20 timestep_fs: 0.5\n\
      \x20 n_steps: 1000\n\
      \x20 temperature_k: 300.0\n\
      \x20 thermostat: velocity_rescaling\n";
    let r = parse_input_str(yaml).unwrap();
    assert_eq!(r.driver, Driver::Md);
    assert_eq!(r.molecule.charge, -1);
    assert_eq!(r.molecule.multiplicity, 2);
    let kw = r.keywords.unwrap();
    assert!(approx(kw.timestep_fs, 0.5));
    assert_eq!(kw.n_steps, 1000);
    assert!(approx(kw.temperature_k, 300.0));
    assert_eq!(kw.thermostat, Thermostat::VelocityRescaling);
  }

  #[test]
  fn test_angstrom_converted_to_bohr() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H, H]\n\
      \x20 geometry: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\
      \x20 units: angstrom\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::Cartesian(c) => {
        assert!(approx(c.x[0], BOHR));
        assert!(approx(c.x[1], 0.0));
      }
      _ => panic!("expected Cartesian"),
    }
  }

  #[test]
  fn test_bohr_stored_unchanged() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H, H]\n\
      \x20 geometry: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\
      \x20 units: bohr\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::Cartesian(c) => {
        assert!(approx(c.x[0], 1.0));
        assert!(approx(c.x[1], 0.0));
      }
      _ => panic!("expected Cartesian"),
    }
  }

  #[test]
  fn test_missing_units_defaults_to_angstrom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [1.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::Cartesian(c) => assert!(approx(c.x[0], BOHR)),
      _ => panic!("expected Cartesian"),
    }
  }

  #[test]
  fn test_element_symbols_normalised_to_title_case() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [o, H, FE]\n\
      \x20 geometry: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::Cartesian(c) => {
        assert_eq!(c.symbols, vec!["O", "H", "Fe"]);
      }
      _ => panic!("expected Cartesian"),
    }
  }

  #[test]
  fn test_driver_gradient_accepted() {
    let yaml = energy_yaml().replace("energy", "gradient");
    let r = parse_input_str(&yaml).unwrap();
    assert_eq!(r.driver, Driver::Gradient);
  }

  #[test]
  fn test_driver_hessian_accepted() {
    let yaml = energy_yaml().replace("energy", "hessian");
    let r = parse_input_str(&yaml).unwrap();
    assert_eq!(r.driver, Driver::Hessian);
  }

  #[test]
  fn test_md_keywords_optional_fields_default() {
    let r = parse_input_str(md_yaml()).unwrap();
    let kw = r.keywords.unwrap();
    assert!(approx(kw.temperature_k, 0.0));
    assert_eq!(kw.thermostat, Thermostat::None);
  }

  #[test]
  fn test_keywords_ignored_for_non_md_driver() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n\
      keywords:\n\
      \x20 timestep_fs: 0.5\n\
      \x20 n_steps: 100\n";
    let r = parse_input_str(yaml).unwrap();
    assert_eq!(r.driver, Driver::Energy);
    assert_eq!(r.keywords, Option::None);
  }

  // ── Happy paths: Z-matrix ───────────────────────────────────────────────────

  #[test]
  fn test_zmatrix_one_atom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::ZMatrix(z) => {
        assert_eq!(z.symbols, vec!["O"]);
        assert_eq!(z.bond_atoms[0], Option::None);
        assert_eq!(z.bond_lengths_bohr[0], Option::None);
      }
      _ => panic!("expected ZMatrix"),
    }
  }

  #[test]
  fn test_zmatrix_two_atoms() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.9572\n\
      \x20 units: angstrom\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::ZMatrix(z) => {
        assert_eq!(z.symbols, vec!["O", "H"]);
        assert_eq!(z.bond_atoms[1], Some(1));
        assert!(approx(z.bond_lengths_bohr[1].unwrap(), 0.9572 * BOHR));
      }
      _ => panic!("expected ZMatrix"),
    }
  }

  #[test]
  fn test_zmatrix_three_atoms() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 104.5\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::ZMatrix(z) => {
        assert_eq!(z.angle_atoms[2], Some(2));
        assert!(approx(z.angles_deg[2].unwrap(), 104.5));
        assert_eq!(z.dihedral_atoms[2], Option::None);
      }
      _ => panic!("expected ZMatrix"),
    }
  }

  #[test]
  fn test_zmatrix_four_atoms_with_dihedral() {
    let r = parse_input_str(zmat_energy_yaml()).unwrap();
    match &r.molecule.geometry {
      Geometry::ZMatrix(z) => {
        assert_eq!(z.dihedral_atoms[3], Some(3));
        assert!(approx(z.dihedrals_deg[3].unwrap(), 120.0));
      }
      _ => panic!("expected ZMatrix"),
    }
  }

  #[test]
  fn test_zmatrix_bond_length_converted_to_bohr() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 1.0\n\
      \x20 units: angstrom\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    let r = parse_input_str(yaml).unwrap();
    match &r.molecule.geometry {
      Geometry::ZMatrix(z) => {
        assert!(approx(z.bond_lengths_bohr[1].unwrap(), BOHR));
      }
      _ => panic!("expected ZMatrix"),
    }
  }

  // ── File I/O ────────────────────────────────────────────────────────────────

  #[test]
  fn test_parse_input_valid_file() {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(energy_yaml().as_bytes()).unwrap();
    assert!(parse_input(f.path()).is_ok());
  }

  #[test]
  fn test_parse_input_missing_file() {
    let path = std::path::Path::new("/tmp/nonexistent_input_file_xyz.yaml");
    assert!(matches!(parse_input(path), Err(InputError::IoError(_))));
  }

  // ── YAML errors ─────────────────────────────────────────────────────────────

  #[test]
  fn test_invalid_yaml() {
    let yaml = "driver: md\nmolecule: :\n  bad:";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidYaml(_))
    ));
  }

  // ── Missing required top-level fields ───────────────────────────────────────

  #[test]
  fn test_missing_driver() {
    let yaml = "molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("driver".to_string()))
    );
  }

  #[test]
  fn test_missing_molecule() {
    let yaml = "driver: energy\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("molecule".to_string()))
    );
  }

  #[test]
  fn test_missing_model() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("model".to_string()))
    );
  }

  #[test]
  fn test_missing_model_method() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("model.method".to_string()))
    );
  }

  #[test]
  fn test_missing_model_basis() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("model.basis".to_string()))
    );
  }

  #[test]
  fn test_missing_keywords_block_for_md() {
    let yaml = "driver: md\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("keywords".to_string()))
    );
  }

  #[test]
  fn test_missing_keywords_timestep_fs() {
    let yaml = "driver: md\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n\
      keywords:\n\
      \x20 n_steps: 100\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("keywords.timestep_fs".to_string()))
    );
  }

  #[test]
  fn test_missing_keywords_n_steps() {
    let yaml = "driver: md\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n\
      keywords:\n\
      \x20 timestep_fs: 0.5\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("keywords.n_steps".to_string()))
    );
  }

  #[test]
  fn test_missing_geometry_only_symbols_present() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("molecule.geometry".to_string()))
    );
  }

  #[test]
  fn test_missing_symbols_only_geometry_present() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("molecule.symbols".to_string()))
    );
  }

  #[test]
  fn test_missing_geometry_neither_cartesian_nor_zmatrix() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 charge: 0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("molecule.geometry".to_string()))
    );
  }

  #[test]
  fn test_empty_z_matrix_sequence() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix: []\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::MissingField("molecule.z_matrix".to_string()))
    );
  }

  // ── Unknown fields ───────────────────────────────────────────────────────────

  #[test]
  fn test_unknown_top_level_key() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n\
      extra_key: oops\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::UnknownField("extra_key".to_string()))
    );
  }

  // ── Invalid values ───────────────────────────────────────────────────────────

  #[test]
  fn test_unrecognised_driver() {
    let yaml = energy_yaml().replace("energy", "optimize");
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "driver"
    ));
  }

  #[test]
  fn test_unrecognised_units() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      \x20 units: nanometer\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "molecule.units"
    ));
  }

  #[test]
  fn test_empty_method() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: \"\"\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "model.method"
    ));
  }

  #[test]
  fn test_empty_basis() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: \"\"\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "model.basis"
    ));
  }

  #[test]
  fn test_timestep_fs_zero() {
    let yaml = md_yaml().replace("timestep_fs: 0.5", "timestep_fs: 0.0");
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "keywords.timestep_fs"
    ));
  }

  #[test]
  fn test_timestep_fs_negative() {
    let yaml = md_yaml().replace("timestep_fs: 0.5", "timestep_fs: -1.0");
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "keywords.timestep_fs"
    ));
  }

  #[test]
  fn test_n_steps_zero() {
    let yaml = md_yaml().replace("n_steps: 100", "n_steps: 0");
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "keywords.n_steps"
    ));
  }

  #[test]
  fn test_temperature_k_negative() {
    let yaml = format!("{}  temperature_k: -1.0\n", md_yaml());
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "keywords.temperature_k"
    ));
  }

  #[test]
  fn test_unrecognised_thermostat() {
    let yaml = format!("{}  thermostat: langevin\n", md_yaml());
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "keywords.thermostat"
    ));
  }

  #[test]
  fn test_multiplicity_zero() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      \x20 multiplicity: 0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidValue { field, .. }) if field == "molecule.multiplicity"
    ));
  }

  // ── Geometry format errors ───────────────────────────────────────────────────

  #[test]
  fn test_ambiguous_geometry() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H]\n\
      \x20 geometry: [0.0, 0.0, 0.0]\n\
      \x20 z_matrix:\n\
      \x20   - symbol: H\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(parse_input_str(yaml), Err(InputError::AmbiguousGeometry));
  }

  #[test]
  fn test_coordinate_mismatch_not_multiple_of_three() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H, H]\n\
      \x20 geometry: [0.0, 0.0, 0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::CoordinateMismatch { n_symbols: 2, n_coords: 5 })
    );
  }

  #[test]
  fn test_coordinate_mismatch_too_short() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H, H, H]\n\
      \x20 geometry: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::CoordinateMismatch { n_symbols: 3, n_coords: 6 })
    );
  }

  // ── Element validation ───────────────────────────────────────────────────────

  #[test]
  fn test_unknown_element_in_cartesian() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 symbols: [H, Xx]\n\
      \x20 geometry: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::InvalidElement("Xx".to_string()))
    );
  }

  #[test]
  fn test_unknown_element_in_zmatrix() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: Zz\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert_eq!(
      parse_input_str(yaml),
      Err(InputError::InvalidElement("Zz".to_string()))
    );
  }

  // ── Z-matrix structural errors ───────────────────────────────────────────────

  #[test]
  fn test_zmat_row0_with_bond_atom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20     bond_atom: 1\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 0, .. })
    ));
  }

  #[test]
  fn test_zmat_row1_missing_bond_atom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_length: 0.96\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_row1_with_angle_atom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 1\n\
      \x20     angle: 90.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_row2_missing_angle() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 2, .. })
    ));
  }

  #[test]
  fn test_zmat_row2_with_dihedral_atom() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 104.5\n\
      \x20     dihedral_atom: 1\n\
      \x20     dihedral: 0.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 2, .. })
    ));
  }

  #[test]
  fn test_zmat_row3_missing_dihedral_atom() {
    // Use the 4-atom template and remove dihedral_atom from row 3.
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 104.5\n\
      \x20   - symbol: C\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 1.5\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 109.5\n\
      \x20     dihedral: 120.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 3, .. })
    ));
  }

  #[test]
  fn test_zmat_ref_index_zero() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 0\n\
      \x20     bond_length: 0.96\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_forward_reference() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 2\n\
      \x20     bond_length: 0.96\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_duplicate_reference_indices() {
    // bond_atom == angle_atom at row 3 (0-based).
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 104.5\n\
      \x20   - symbol: C\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 1.5\n\
      \x20     angle_atom: 1\n\
      \x20     angle: 109.5\n\
      \x20     dihedral_atom: 3\n\
      \x20     dihedral: 120.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 3, .. })
    ));
  }

  #[test]
  fn test_zmat_bond_length_zero() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_bond_length_negative() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: -1.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 1, .. })
    ));
  }

  #[test]
  fn test_zmat_angle_zero() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 0.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 2, .. })
    ));
  }

  #[test]
  fn test_zmat_angle_180() {
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 180.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 2, .. })
    ));
  }

  #[test]
  fn test_zmat_dihedral_out_of_range() {
    let yaml = zmat_energy_yaml().replace("dihedral: 120.0", "dihedral: 181.0");
    assert!(matches!(
      parse_input_str(&yaml),
      Err(InputError::InvalidZMatrix { row: 3, .. })
    ));
  }

  #[test]
  fn test_zmat_error_identifies_correct_row() {
    // 5-atom z-matrix where only row 4 has an invalid bond_length.
    let yaml = "driver: energy\n\
      molecule:\n\
      \x20 z_matrix:\n\
      \x20   - symbol: O\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20   - symbol: H\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 0.96\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 104.5\n\
      \x20   - symbol: C\n\
      \x20     bond_atom: 1\n\
      \x20     bond_length: 1.5\n\
      \x20     angle_atom: 2\n\
      \x20     angle: 109.5\n\
      \x20     dihedral_atom: 3\n\
      \x20     dihedral: 120.0\n\
      \x20   - symbol: N\n\
      \x20     bond_atom: 4\n\
      \x20     bond_length: -1.0\n\
      \x20     angle_atom: 1\n\
      \x20     angle: 90.0\n\
      \x20     dihedral_atom: 2\n\
      \x20     dihedral: 0.0\n\
      model:\n\
      \x20 method: hf\n\
      \x20 basis: sto-3g\n";
    assert!(matches!(
      parse_input_str(yaml),
      Err(InputError::InvalidZMatrix { row: 4, .. })
    ));
  }
}
