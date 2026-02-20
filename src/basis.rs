use std::path::{Path, PathBuf};

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

const BSE_BASE_URL: &str = "https://www.basissetexchange.org";
const DEFAULT_CACHE_ROOT: &str = "data/basis";

#[derive(Debug, PartialEq)]
pub enum BseError {
  InvalidElement(String),
  InvalidBasisSetName(String),
  ElementNotInBasisSet { element: String, basis_name: String },
  UnknownBasisSet(String),
  NetworkError(String),
  IoError(String),
  InvalidResponse(String),
}

impl std::fmt::Display for BseError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      BseError::InvalidElement(s) =>
        write!(f, "invalid element: {:?}", s),
      BseError::InvalidBasisSetName(s) =>
        write!(f, "invalid basis set name: {:?}", s),
      BseError::ElementNotInBasisSet { element, basis_name } =>
        write!(f, "element {} not found in basis set {}", element, basis_name),
      BseError::UnknownBasisSet(s) =>
        write!(f, "unknown basis set: {:?}", s),
      BseError::NetworkError(s) =>
        write!(f, "network error: {}", s),
      BseError::IoError(s) =>
        write!(f, "I/O error: {}", s),
      BseError::InvalidResponse(s) =>
        write!(f, "invalid response: {}", s),
    }
  }
}

/// Downloads (if needed) and returns the path to the cached QCSchema JSON file
/// for `element` in `basis_name`, using `data/basis` as the cache root.
pub fn fetch_basis(element: &str, basis_name: &str) -> Result<PathBuf, BseError> {
  fetch_basis_impl(element, basis_name, BSE_BASE_URL, Path::new(DEFAULT_CACHE_ROOT))
}

fn normalize_element(element: &str) -> Result<String, BseError> {
  let mut chars = element.chars();
  let normalized = match chars.next() {
    None => String::new(),
    Some(c) => {
      let upper: String = c.to_uppercase().collect();
      upper + &chars.as_str().to_lowercase()
    }
  };
  if ELEMENTS.contains(&normalized.as_str()) {
    Ok(normalized)
  } else {
    Err(BseError::InvalidElement(element.to_string()))
  }
}

fn is_valid_cache(path: &Path) -> bool {
  match std::fs::read_to_string(path) {
    Err(_) => false,
    Ok(content) => {
      !content.is_empty()
        && serde_json::from_str::<serde_json::Value>(&content).is_ok()
    }
  }
}

fn elements_field_is_empty(json: &serde_json::Value) -> bool {
  match json.get("elements") {
    None => true,
    Some(serde_json::Value::Object(m)) => m.is_empty(),
    Some(serde_json::Value::Array(a)) => a.is_empty(),
    Some(_) => false,
  }
}

fn fetch_basis_impl(
  element: &str,
  basis_name: &str,
  base_url: &str,
  cache_root: &Path,
) -> Result<PathBuf, BseError> {
  if basis_name.is_empty() {
    return Err(BseError::InvalidBasisSetName(basis_name.to_string()));
  }

  let element_norm = normalize_element(element)?;
  let basis_norm = basis_name.to_lowercase();

  let path = cache_root
    .join(&basis_norm)
    .join(format!("{}.json", element_norm));

  if is_valid_cache(&path) {
    return Ok(path);
  }

  let url = format!(
    "{}/api/basis/{}/format/qcschema?elements={}",
    base_url, basis_norm, element_norm
  );

  let response = reqwest::blocking::get(&url)
    .map_err(|e| BseError::NetworkError(e.to_string()))?;

  match response.status().as_u16() {
    404 => return Err(BseError::UnknownBasisSet(basis_norm)),
    200 => {}
    code => return Err(BseError::NetworkError(format!("unexpected HTTP status {}", code))),
  }

  let body = response.text()
    .map_err(|e| BseError::NetworkError(e.to_string()))?;

  let json: serde_json::Value = serde_json::from_str(&body)
    .map_err(|e| BseError::InvalidResponse(e.to_string()))?;

  if elements_field_is_empty(&json) {
    return Err(BseError::ElementNotInBasisSet {
      element: element_norm,
      basis_name: basis_norm,
    });
  }

  let dir = path.parent().expect("cache path always has a parent");
  std::fs::create_dir_all(dir)
    .map_err(|e| BseError::IoError(e.to_string()))?;
  std::fs::write(&path, &body)
    .map_err(|e| BseError::IoError(e.to_string()))?;

  Ok(path)
}

// ============================================================================
// Types for parsed basis sets
// ============================================================================

#[derive(Debug, PartialEq)]
pub struct ElectronShell {
  pub angular_momentum: u32,
  pub exponents: Vec<f64>,
  pub coefficients: Vec<f64>,
}

#[derive(Debug, PartialEq)]
pub struct BasisSet {
  pub element: String,
  pub atomic_number: u32,
  pub shells: Vec<ElectronShell>,
}

#[derive(Debug, PartialEq)]
pub enum ParseError {
  IoError(String),
  InvalidJson(String),
  MultipleElements { found: usize },
  NoElements,
  InvalidAtomicNumber(String),
  NoElectronShells,
  MalformedShell { index: usize, reason: String },
}

impl std::fmt::Display for ParseError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      ParseError::IoError(s) => write!(f, "I/O error: {}", s),
      ParseError::InvalidJson(s) => write!(f, "invalid JSON: {}", s),
      ParseError::MultipleElements { found } =>
        write!(f, "expected 1 element, found {}", found),
      ParseError::NoElements => write!(f, "elements object is absent or empty"),
      ParseError::InvalidAtomicNumber(s) =>
        write!(f, "invalid atomic number: {:?}", s),
      ParseError::NoElectronShells => write!(f, "no electron shells found"),
      ParseError::MalformedShell { index, reason } =>
        write!(f, "shell {}: {}", index, reason),
    }
  }
}

#[derive(Debug)]
pub enum LoadError {
  Fetch(BseError),
  Parse(ParseError),
}

impl std::fmt::Display for LoadError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      LoadError::Fetch(e) => write!(f, "fetch error: {}", e),
      LoadError::Parse(e) => write!(f, "parse error: {}", e),
    }
  }
}

// ============================================================================
// parse_basis
// ============================================================================

/// Reads and parses a QCSchema basis set JSON file.
///
/// The file must contain exactly one element. SP shells (multiple entries in
/// `angular_momentum`) are split into separate shells, one per angular momentum
/// value, each sharing the original exponents.
pub fn parse_basis(path: &Path) -> Result<BasisSet, ParseError> {
  let content = std::fs::read_to_string(path)
    .map_err(|e| ParseError::IoError(e.to_string()))?;

  let json: serde_json::Value = serde_json::from_str(&content)
    .map_err(|e| ParseError::InvalidJson(e.to_string()))?;

  let elements_obj = json
    .get("elements")
    .and_then(|v| v.as_object())
    .ok_or(ParseError::NoElements)?;

  let z_str: &str = match elements_obj.len() {
    0 => return Err(ParseError::NoElements),
    1 => elements_obj.keys().next().unwrap(),
    n => return Err(ParseError::MultipleElements { found: n }),
  };

  let z: u32 = match z_str.parse::<u32>() {
    Ok(n) if (1..=118).contains(&n) => n,
    _ => return Err(ParseError::InvalidAtomicNumber(z_str.to_string())),
  };

  let symbol = ELEMENTS[(z - 1) as usize].to_string();
  let element_data = &elements_obj[z_str];

  let shells_raw = element_data
    .get("electron_shells")
    .and_then(|v| v.as_array())
    .ok_or(ParseError::NoElectronShells)?;

  if shells_raw.is_empty() {
    return Err(ParseError::NoElectronShells);
  }

  let mut shells: Vec<ElectronShell> = Vec::new();

  for (idx, shell_val) in shells_raw.iter().enumerate() {
    let am_arr = shell_val
      .get("angular_momentum")
      .and_then(|v| v.as_array())
      .ok_or_else(|| ParseError::MalformedShell {
        index: idx,
        reason: "missing or invalid angular_momentum".to_string(),
      })?;

    if am_arr.is_empty() {
      return Err(ParseError::MalformedShell {
        index: idx,
        reason: "angular_momentum is empty".to_string(),
      });
    }

    let angular_momenta: Vec<u32> = am_arr
      .iter()
      .map(|v| {
        v.as_u64().map(|n| n as u32).ok_or_else(|| ParseError::MalformedShell {
          index: idx,
          reason: format!("angular_momentum entry {:?} is not a non-negative integer", v),
        })
      })
      .collect::<Result<_, _>>()?;

    let exp_arr = shell_val
      .get("exponents")
      .and_then(|v| v.as_array())
      .ok_or_else(|| ParseError::MalformedShell {
        index: idx,
        reason: "missing or invalid exponents".to_string(),
      })?;

    let exponents: Vec<f64> = exp_arr
      .iter()
      .map(|v| {
        let s = v.as_str().ok_or_else(|| ParseError::MalformedShell {
          index: idx,
          reason: format!("exponent {:?} is not a string", v),
        })?;
        s.parse::<f64>().map_err(|_| ParseError::MalformedShell {
          index: idx,
          reason: format!("cannot parse exponent {:?} as f64", s),
        })
      })
      .collect::<Result<_, _>>()?;

    let coeff_arr = shell_val
      .get("coefficients")
      .and_then(|v| v.as_array())
      .ok_or_else(|| ParseError::MalformedShell {
        index: idx,
        reason: "missing or invalid coefficients".to_string(),
      })?;

    if coeff_arr.len() != angular_momenta.len() {
      return Err(ParseError::MalformedShell {
        index: idx,
        reason: format!(
          "expected {} coefficient vector(s) to match angular_momentum, found {}",
          angular_momenta.len(),
          coeff_arr.len()
        ),
      });
    }

    for (am, coeff_vec_val) in angular_momenta.iter().zip(coeff_arr.iter()) {
      let coeff_vec =
        coeff_vec_val.as_array().ok_or_else(|| ParseError::MalformedShell {
          index: idx,
          reason: "coefficient set is not an array".to_string(),
        })?;

      if coeff_vec.len() != exponents.len() {
        return Err(ParseError::MalformedShell {
          index: idx,
          reason: format!(
            "coefficient vector has {} values but there are {} exponents",
            coeff_vec.len(),
            exponents.len()
          ),
        });
      }

      let coefficients: Vec<f64> = coeff_vec
        .iter()
        .map(|v| {
          let s = v.as_str().ok_or_else(|| ParseError::MalformedShell {
            index: idx,
            reason: format!("coefficient {:?} is not a string", v),
          })?;
          s.parse::<f64>().map_err(|_| ParseError::MalformedShell {
            index: idx,
            reason: format!("cannot parse coefficient {:?} as f64", s),
          })
        })
        .collect::<Result<_, _>>()?;

      shells.push(ElectronShell {
        angular_momentum: *am,
        exponents: exponents.clone(),
        coefficients,
      });
    }
  }

  Ok(BasisSet { element: symbol, atomic_number: z, shells })
}

// ============================================================================
// load_basis
// ============================================================================

/// Fetches (if needed) and parses the QCSchema basis set for `element` in
/// `basis_name`, using the live BSE API and `data/basis` as the cache root.
pub fn load_basis(element: &str, basis_name: &str) -> Result<BasisSet, LoadError> {
  load_basis_impl(element, basis_name, BSE_BASE_URL, Path::new(DEFAULT_CACHE_ROOT))
}

fn load_basis_impl(
  element: &str,
  basis_name: &str,
  base_url: &str,
  cache_root: &Path,
) -> Result<BasisSet, LoadError> {
  let path = fetch_basis_impl(element, basis_name, base_url, cache_root)
    .map_err(LoadError::Fetch)?;
  parse_basis(&path).map_err(LoadError::Parse)
}

#[cfg(test)]
mod tests {
  use super::*;
  use tempfile::TempDir;

  const VALID_RESPONSE: &str = r#"{"elements":{"1":{"electron_shells":[]}}}"#;

  fn temp_dir() -> TempDir {
    TempDir::new().expect("create temp dir")
  }

  // -------------------------------------------------------------------------
  // Happy paths
  // -------------------------------------------------------------------------

  // Scenario: Download a basis set that is not cached
  #[test]
  fn download_basis_not_cached() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("H", "sto-3g", &server.url(), dir.path());
    let path = result.expect("should succeed");

    assert_eq!(path, dir.path().join("sto-3g").join("H.json"));
    let content = std::fs::read_to_string(&path).expect("file written");
    assert_eq!(content, VALID_RESPONSE);
  }

  // Scenario: Return cached file when a valid cache exists
  // (also verifies no HTTP request is made)
  #[test]
  fn return_cached_when_valid() {
    let dir = temp_dir();
    let cache_path = dir.path().join("sto-3g").join("H.json");
    std::fs::create_dir_all(cache_path.parent().unwrap()).unwrap();
    std::fs::write(&cache_path, VALID_RESPONSE).unwrap();

    // "http://localhost:1" is unreachable; any HTTP request would produce NetworkError
    let result = fetch_basis_impl("H", "sto-3g", "http://localhost:1", dir.path());
    assert_eq!(result.expect("should succeed"), cache_path);
  }

  // Scenario: Create data/basis directory if it does not exist
  #[test]
  fn creates_missing_data_basis_directory() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    let cache_root = dir.path().join("data").join("basis");
    // cache_root does not exist yet
    fetch_basis_impl("H", "sto-3g", &server.url(), &cache_root)
      .expect("should succeed");

    assert!(cache_root.join("sto-3g").join("H.json").exists());
  }

  // Scenario: Create basis-set subdirectory if it does not exist
  #[test]
  fn creates_missing_basis_subdirectory() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    // cache root exists; subdirectory does not
    std::fs::create_dir_all(dir.path()).unwrap();

    fetch_basis_impl("H", "sto-3g", &server.url(), dir.path())
      .expect("should succeed");

    assert!(dir.path().join("sto-3g").exists());
    assert!(dir.path().join("sto-3g").join("H.json").exists());
  }

  // Scenario: Basis set name is normalized to lowercase in the file path
  #[test]
  fn basis_name_normalized_to_lowercase() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("H", "STO-3G", &server.url(), dir.path());
    assert_eq!(result.expect("should succeed"), dir.path().join("sto-3g").join("H.json"));
  }

  // Scenario: Element symbol is normalized to title case in the file path
  #[test]
  fn element_normalized_to_title_case() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("h", "sto-3g", &server.url(), dir.path());
    assert_eq!(result.expect("should succeed"), dir.path().join("sto-3g").join("H.json"));
  }

  // -------------------------------------------------------------------------
  // Input validation
  // -------------------------------------------------------------------------

  // Scenario: Reject an unrecognised element symbol
  #[test]
  fn reject_unrecognised_element() {
    let dir = temp_dir();
    let result = fetch_basis_impl("Xx", "sto-3g", "http://localhost:1", dir.path());
    assert!(matches!(result, Err(BseError::InvalidElement(s)) if s == "Xx"));
  }

  // Scenario: Reject an empty element symbol
  #[test]
  fn reject_empty_element() {
    let dir = temp_dir();
    let result = fetch_basis_impl("", "sto-3g", "http://localhost:1", dir.path());
    assert!(matches!(result, Err(BseError::InvalidElement(s)) if s.is_empty()));
  }

  // Scenario: Reject an empty basis set name
  #[test]
  fn reject_empty_basis_name() {
    let dir = temp_dir();
    let result = fetch_basis_impl("H", "", "http://localhost:1", dir.path());
    assert!(matches!(result, Err(BseError::InvalidBasisSetName(s)) if s.is_empty()));
  }

  // -------------------------------------------------------------------------
  // BSE API error responses
  // -------------------------------------------------------------------------

  // Scenario: Basis set name is not known to the BSE (HTTP 404)
  #[test]
  fn unknown_basis_set_on_404() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/unknown-basis/format/qcschema?elements=H")
      .with_status(404)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("H", "unknown-basis", &server.url(), dir.path());
    assert!(matches!(result, Err(BseError::UnknownBasisSet(s)) if s == "unknown-basis"));
    assert!(!dir.path().join("unknown-basis").join("H.json").exists());
  }

  // Scenario: Element is not included in the requested basis set
  #[test]
  fn element_not_in_basis_set_on_empty_elements() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=Au")
      .with_status(200)
      .with_body(r#"{"elements":{}}"#)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("Au", "sto-3g", &server.url(), dir.path());
    assert!(
      matches!(result, Err(BseError::ElementNotInBasisSet { ref element, ref basis_name })
        if element == "Au" && basis_name == "sto-3g"),
      "expected ElementNotInBasisSet, got {:?}",
      result
    );
    assert!(!dir.path().join("sto-3g").join("Au.json").exists());
  }

  // Scenario: BSE API returns an unexpected HTTP status code
  #[test]
  fn network_error_on_unexpected_status() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(500)
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("H", "sto-3g", &server.url(), dir.path());
    assert!(matches!(result, Err(BseError::NetworkError(_))));
    assert!(!dir.path().join("sto-3g").join("H.json").exists());
  }

  // Scenario: BSE API is unreachable
  #[test]
  fn network_error_when_unreachable() {
    let dir = temp_dir();
    // Port 59999 is highly unlikely to be bound
    let result = fetch_basis_impl("H", "sto-3g", "http://127.0.0.1:59999", dir.path());
    assert!(matches!(result, Err(BseError::NetworkError(_))));
  }

  // Scenario: BSE API returns a response that is not valid JSON
  #[test]
  fn invalid_response_on_non_json_body() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body("this is not json")
      .create();

    let dir = temp_dir();
    let result = fetch_basis_impl("H", "sto-3g", &server.url(), dir.path());
    assert!(matches!(result, Err(BseError::InvalidResponse(_))));
    assert!(!dir.path().join("sto-3g").join("H.json").exists());
  }

  // -------------------------------------------------------------------------
  // Cache validation
  // -------------------------------------------------------------------------

  // Scenario: Re-download when cached file is empty
  #[test]
  fn redownload_when_cache_is_empty() {
    let mut server = mockito::Server::new();
    let mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .expect(1)
      .create();

    let dir = temp_dir();
    let cache_path = dir.path().join("sto-3g").join("H.json");
    std::fs::create_dir_all(cache_path.parent().unwrap()).unwrap();
    std::fs::write(&cache_path, "").unwrap();

    fetch_basis_impl("H", "sto-3g", &server.url(), dir.path())
      .expect("should succeed");

    mock.assert();
    assert_eq!(std::fs::read_to_string(&cache_path).unwrap(), VALID_RESPONSE);
  }

  // Scenario: Re-download when cached file contains invalid JSON
  #[test]
  fn redownload_when_cache_contains_invalid_json() {
    let mut server = mockito::Server::new();
    let mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .expect(1)
      .create();

    let dir = temp_dir();
    let cache_path = dir.path().join("sto-3g").join("H.json");
    std::fs::create_dir_all(cache_path.parent().unwrap()).unwrap();
    std::fs::write(&cache_path, "{ invalid json }").unwrap();

    fetch_basis_impl("H", "sto-3g", &server.url(), dir.path())
      .expect("should succeed");

    mock.assert();
    assert_eq!(std::fs::read_to_string(&cache_path).unwrap(), VALID_RESPONSE);
  }

  // -------------------------------------------------------------------------
  // Filesystem errors
  // -------------------------------------------------------------------------

  // Scenario: Filesystem error when creating the directory
  // Simulate by placing a regular file at the path where the directory should be.
  #[test]
  fn io_error_when_dir_cannot_be_created() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    // Place a file at the path where the "sto-3g" directory should be created
    std::fs::write(dir.path().join("sto-3g"), "").unwrap();

    let result = fetch_basis_impl("H", "sto-3g", &server.url(), dir.path());
    assert!(matches!(result, Err(BseError::IoError(_))));
  }

  // Scenario: Filesystem error when writing the downloaded file
  // Simulate by placing a directory at the path where the JSON file should be written.
  #[test]
  fn io_error_when_file_cannot_be_written() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(VALID_RESPONSE)
      .create();

    let dir = temp_dir();
    // Place a directory at the path where H.json should be written
    std::fs::create_dir_all(dir.path().join("sto-3g").join("H.json")).unwrap();

    let result = fetch_basis_impl("H", "sto-3g", &server.url(), dir.path());
    assert!(matches!(result, Err(BseError::IoError(_))));
  }

  // ==========================================================================
  // parse_basis tests
  // ==========================================================================

  // JSON fixtures used across parse tests.
  // A real STO-3G hydrogen s-shell.
  const H_1S: &str = r#"{"elements":{"1":{"electron_shells":[
    {"function_type":"gto","angular_momentum":[0],
     "exponents":["3.4252509","0.6239137","0.1688554"],
     "coefficients":[["0.1543290","0.5353281","0.4446345"]]}
  ]}}}"#;

  // STO-3G carbon: one s-shell (Z=6) and one p-shell.
  const C_1S_1P: &str = r#"{"elements":{"6":{"electron_shells":[
    {"function_type":"gto","angular_momentum":[0],
     "exponents":["71.6168370","13.0450963","3.5305122"],
     "coefficients":[["0.1543290","0.5353281","0.4446345"]]},
    {"function_type":"gto","angular_momentum":[1],
     "exponents":["2.9412494","0.6834831","0.2222899"],
     "coefficients":[["0.2364600","0.8768660","0.2364600"]]}
  ]}}}"#;

  // Lithium with a single SP shell (angular_momentum [0,1], 2 coeff vectors).
  const LI_SP: &str = r#"{"elements":{"3":{"electron_shells":[
    {"function_type":"gto","angular_momentum":[0,1],
     "exponents":["16.1195750","2.9362007","0.7946505"],
     "coefficients":[
       ["0.1543290","0.5353281","0.4446345"],
       ["0.2494820","0.8657560","0.0000000"]
     ]}
  ]}}}"#;

  // Copper with one s-shell and a spurious ecp_potentials key.
  const CU_WITH_ECP: &str = r#"{"elements":{"29":{"electron_shells":[
    {"function_type":"gto","angular_momentum":[0],
     "exponents":["1.2","0.5"],
     "coefficients":[["0.8","0.4"]]}
  ],"ecp_potentials":[{"ecp_type":"scalar"}]}}}"#;

  fn write_json(dir: &TempDir, content: &str) -> std::path::PathBuf {
    let path = dir.path().join("test.json");
    std::fs::write(&path, content).unwrap();
    path
  }

  // Scenario: Parse a valid single-shell file
  #[test]
  fn parse_valid_single_shell() {
    let dir = temp_dir();
    let path = write_json(&dir, H_1S);
    let bs = parse_basis(&path).expect("should succeed");
    assert_eq!(bs.element, "H");
    assert_eq!(bs.atomic_number, 1);
    assert_eq!(bs.shells.len(), 1);
    assert_eq!(bs.shells[0].angular_momentum, 0);
    assert_eq!(bs.shells[0].exponents.len(), 3);
    assert_eq!(bs.shells[0].coefficients.len(), 3);
  }

  // Scenario: Parse a file with multiple shells
  #[test]
  fn parse_multiple_shells() {
    let dir = temp_dir();
    let path = write_json(&dir, C_1S_1P);
    let bs = parse_basis(&path).expect("should succeed");
    assert_eq!(bs.shells.len(), 2);
    assert_eq!(bs.shells[0].angular_momentum, 0);
    assert_eq!(bs.shells[1].angular_momentum, 1);
  }

  // Scenario: SP shell is split into separate S and P shells
  #[test]
  fn sp_shell_is_split() {
    let dir = temp_dir();
    let path = write_json(&dir, LI_SP);
    let bs = parse_basis(&path).expect("should succeed");
    assert_eq!(bs.shells.len(), 2);
    assert_eq!(bs.shells[0].angular_momentum, 0);
    assert_eq!(bs.shells[1].angular_momentum, 1);
    // Both shells share the same exponents
    let exp = &[16.1195750_f64, 2.9362007, 0.7946505];
    for (actual, expected) in bs.shells[0].exponents.iter().zip(exp.iter()) {
      assert!((actual - expected).abs() < 1e-6);
    }
    for (actual, expected) in bs.shells[1].exponents.iter().zip(exp.iter()) {
      assert!((actual - expected).abs() < 1e-6);
    }
    // S shell gets first coefficient vector
    assert!((bs.shells[0].coefficients[0] - 0.1543290).abs() < 1e-6);
    // P shell gets second coefficient vector
    assert!((bs.shells[1].coefficients[0] - 0.2494820).abs() < 1e-6);
  }

  // Scenario: ECP data is ignored when electron shells are present
  #[test]
  fn ecp_data_is_ignored() {
    let dir = temp_dir();
    let path = write_json(&dir, CU_WITH_ECP);
    let bs = parse_basis(&path).expect("should succeed");
    assert_eq!(bs.element, "Cu");
    assert_eq!(bs.shells.len(), 1);
  }

  // Scenario: File does not exist
  #[test]
  fn file_does_not_exist() {
    let dir = temp_dir();
    let result = parse_basis(&dir.path().join("nonexistent.json"));
    assert!(matches!(result, Err(ParseError::IoError(_))));
  }

  // Scenario: File is not valid JSON
  #[test]
  fn file_is_not_valid_json() {
    let dir = temp_dir();
    let path = write_json(&dir, "{ not valid json }");
    assert!(matches!(parse_basis(&path), Err(ParseError::InvalidJson(_))));
  }

  // Scenario: elements object contains two keys
  #[test]
  fn multiple_elements_error() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[]},"2":{"electron_shells":[]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MultipleElements { found: 2 }))
    );
  }

  // Scenario: elements object is empty
  #[test]
  fn empty_elements_object() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{}}"#);
    assert!(matches!(parse_basis(&path), Err(ParseError::NoElements)));
  }

  // Scenario: Atomic number key is not a number
  #[test]
  fn atomic_number_not_a_number() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{"X":{"electron_shells":[]}}}"#);
    assert!(
      matches!(parse_basis(&path), Err(ParseError::InvalidAtomicNumber(s)) if s == "X")
    );
  }

  // Scenario: Atomic number key is zero
  #[test]
  fn atomic_number_zero() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{"0":{"electron_shells":[]}}}"#);
    assert!(
      matches!(parse_basis(&path), Err(ParseError::InvalidAtomicNumber(s)) if s == "0")
    );
  }

  // Scenario: Atomic number key is out of range
  #[test]
  fn atomic_number_out_of_range() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{"119":{"electron_shells":[]}}}"#);
    assert!(
      matches!(parse_basis(&path), Err(ParseError::InvalidAtomicNumber(s)) if s == "119")
    );
  }

  // Scenario: electron_shells key is absent
  #[test]
  fn electron_shells_absent() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{"1":{}}}"#);
    assert!(matches!(parse_basis(&path), Err(ParseError::NoElectronShells)));
  }

  // Scenario: electron_shells array is empty
  #[test]
  fn electron_shells_empty() {
    let dir = temp_dir();
    let path = write_json(&dir, r#"{"elements":{"1":{"electron_shells":[]}}}"#);
    assert!(matches!(parse_basis(&path), Err(ParseError::NoElectronShells)));
  }

  // Scenario: Shell has empty angular_momentum array
  #[test]
  fn empty_angular_momentum() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[],"exponents":["1.0"],"coefficients":[["1.0"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 0, .. }))
    );
  }

  // Scenario: SP shell has wrong number of coefficient vectors
  #[test]
  fn sp_shell_wrong_coefficient_count() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[0,1],"exponents":["1.0","2.0"],
         "coefficients":[["0.5","0.5"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 0, .. }))
    );
  }

  // Scenario: Coefficient vector length does not match exponent count
  #[test]
  fn coefficient_length_mismatch() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[0],"exponents":["1.0","2.0","3.0"],
         "coefficients":[["0.5","0.5"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 0, .. }))
    );
  }

  // Scenario: Exponent string cannot be parsed as f64
  #[test]
  fn unparseable_exponent() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[0],"exponents":["abc"],"coefficients":[["1.0"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 0, .. }))
    );
  }

  // Scenario: Coefficient string cannot be parsed as f64
  #[test]
  fn unparseable_coefficient() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[0],"exponents":["1.0"],"coefficients":[["abc"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 0, .. }))
    );
  }

  // Scenario: Error report identifies the correct shell index
  #[test]
  fn error_identifies_correct_shell_index() {
    let dir = temp_dir();
    let path = write_json(
      &dir,
      r#"{"elements":{"1":{"electron_shells":[
        {"angular_momentum":[0],"exponents":["1.0"],"coefficients":[["0.5"]]},
        {"angular_momentum":[0],"exponents":["bad"],"coefficients":[["0.5"]]}
      ]}}}"#,
    );
    assert!(
      matches!(parse_basis(&path), Err(ParseError::MalformedShell { index: 1, .. }))
    );
  }

  // ==========================================================================
  // load_basis tests
  // ==========================================================================

  // A full valid QCSchema response usable by both fetch_basis and parse_basis.
  const H_STO3G_FULL: &str = r#"{"elements":{"1":{"electron_shells":[
    {"function_type":"gto","angular_momentum":[0],
     "exponents":["3.4252509","0.6239137","0.1688554"],
     "coefficients":[["0.1543290","0.5353281","0.4446345"]]}
  ]}}}"#;

  // Scenario: load_basis fetches and parses successfully
  #[test]
  fn load_basis_success() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/sto-3g/format/qcschema?elements=H")
      .with_status(200)
      .with_body(H_STO3G_FULL)
      .create();

    let dir = temp_dir();
    let result = load_basis_impl("H", "sto-3g", &server.url(), dir.path());
    let bs = result.expect("should succeed");
    assert_eq!(bs.element, "H");
    assert_eq!(bs.shells.len(), 1);
  }

  // Scenario: load_basis propagates a fetch error
  #[test]
  fn load_basis_propagates_fetch_error() {
    let mut server = mockito::Server::new();
    let _mock = server
      .mock("GET", "/api/basis/unknown-basis/format/qcschema?elements=H")
      .with_status(404)
      .create();

    let dir = temp_dir();
    let result = load_basis_impl("H", "unknown-basis", &server.url(), dir.path());
    assert!(matches!(result, Err(LoadError::Fetch(BseError::UnknownBasisSet(_)))));
  }

  // Scenario: load_basis propagates a parse error
  // Pre-seed the cache with valid JSON that fails semantic parse_basis validation
  // (empty elements object passes is_valid_cache but fails parse_basis).
  #[test]
  fn load_basis_propagates_parse_error() {
    let dir = temp_dir();
    let cache_path = dir.path().join("sto-3g").join("H.json");
    std::fs::create_dir_all(cache_path.parent().unwrap()).unwrap();
    // Valid JSON but semantically invalid for parse_basis (no electron_shells key)
    std::fs::write(&cache_path, r#"{"elements":{"1":{}}}"#).unwrap();

    let result = load_basis_impl("H", "sto-3g", "http://localhost:1", dir.path());
    assert!(matches!(result, Err(LoadError::Parse(_))));
  }
}
