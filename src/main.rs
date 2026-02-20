mod basis;
mod input;
mod orbital;

use std::path::Path;
use std::process;

use input::{parse_input, Driver, Geometry};

/// Core CLI logic. Takes the arguments (excluding argv[0]) and returns either
/// the success line to print on stdout, or the error message to print on stderr
/// (without the "error: " prefix — that is added by `main`).
fn run(args: &[String]) -> Result<String, String> {
  if args.len() != 1 {
    return Err("usage: electron <input-file>".to_string());
  }

  let path = Path::new(&args[0]);

  match parse_input(path) {
    Ok(sim) => {
      let driver = match sim.driver {
        Driver::Energy   => "energy",
        Driver::Gradient => "gradient",
        Driver::Hessian  => "hessian",
        Driver::Md       => "md",
      };
      let atoms = match &sim.molecule.geometry {
        Geometry::Cartesian(c) => c.symbols.len(),
        Geometry::ZMatrix(z)   => z.symbols.len(),
      };
      Ok(format!(
        "Parsed: driver={}, method={}, basis={}, atoms={}",
        driver, sim.model.method, sim.model.basis, atoms
      ))
    }
    Err(e) => Err(e.to_string()),
  }
}

fn main() {
  let args: Vec<String> = std::env::args().skip(1).collect();
  match run(&args) {
    Ok(msg) => println!("{}", msg),
    Err(e) => {
      eprintln!("error: {}", e);
      process::exit(1);
    }
  }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
  use super::*;
  use std::io::Write;

  fn arg(s: &str) -> Vec<String> {
    vec![s.to_string()]
  }

  fn temp_file(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
  }

  const ENERGY_H2: &str = "\
driver: energy
molecule:
  symbols: [H, H]
  geometry: [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]
  units: bohr
model:
  method: hf
  basis: sto-3g
";

  const MD_WATER: &str = "\
driver: md
molecule:
  symbols: [O, H, H]
  geometry: [0.0, 0.0, 0.117, 0.0, 0.757, -0.469, 0.0, -0.757, -0.469]
  units: angstrom
model:
  method: b3lyp
  basis: sto-3g
keywords:
  timestep_fs: 0.5
  n_steps: 100
";

  const ZMAT_4: &str = "\
driver: energy
molecule:
  z_matrix:
    - symbol: O
    - symbol: H
      bond_atom: 1
      bond_length: 0.96
    - symbol: H
      bond_atom: 1
      bond_length: 0.96
      angle_atom: 2
      angle: 104.5
    - symbol: C
      bond_atom: 1
      bond_length: 1.5
      angle_atom: 2
      angle: 109.5
      dihedral_atom: 3
      dihedral: 120.0
model:
  method: hf
  basis: sto-3g
";

  // ── Happy paths ─────────────────────────────────────────────────────────────

  /// Scenario: Parse a valid Cartesian energy input file
  #[test]
  fn test_cartesian_energy_file() {
    let f = temp_file(ENERGY_H2);
    let result = run(&arg(f.path().to_str().unwrap()));
    assert_eq!(
      result.unwrap(),
      "Parsed: driver=energy, method=hf, basis=sto-3g, atoms=2"
    );
  }

  /// Scenario: Parse a valid MD input file
  #[test]
  fn test_md_file() {
    let f = temp_file(MD_WATER);
    let result = run(&arg(f.path().to_str().unwrap()));
    assert_eq!(
      result.unwrap(),
      "Parsed: driver=md, method=b3lyp, basis=sto-3g, atoms=3"
    );
  }

  /// Scenario: Relative path is resolved from the working directory
  #[test]
  fn test_relative_path() {
    // Write a file with a unique name into the current working directory,
    // then pass just the filename as a relative path.
    let filename = format!("__cli_test_{}.yaml", std::process::id());
    std::fs::write(&filename, ENERGY_H2).unwrap();
    let result = run(&arg(&filename));
    std::fs::remove_file(&filename).ok();
    assert!(result.is_ok(), "expected Ok, got: {:?}", result);
  }

  /// Scenario: Absolute path is accepted
  #[test]
  fn test_absolute_path() {
    let f = temp_file(ENERGY_H2);
    let abs = f.path().canonicalize().unwrap();
    let result = run(&arg(abs.to_str().unwrap()));
    assert!(result.is_ok(), "expected Ok, got: {:?}", result);
  }

  /// Scenario: Atom count for Z-matrix geometry is the number of rows
  #[test]
  fn test_zmatrix_atom_count() {
    let f = temp_file(ZMAT_4);
    let msg = run(&arg(f.path().to_str().unwrap())).unwrap();
    assert!(msg.contains("atoms=4"), "output was: {msg}");
  }

  // ── Argument errors ─────────────────────────────────────────────────────────

  /// Scenario: No argument given
  #[test]
  fn test_no_arguments() {
    let result = run(&[]);
    assert_eq!(result.unwrap_err(), "usage: electron <input-file>");
  }

  /// Scenario: More than one argument given
  #[test]
  fn test_too_many_arguments() {
    let result = run(&["a.yaml".to_string(), "b.yaml".to_string()]);
    assert_eq!(result.unwrap_err(), "usage: electron <input-file>");
  }

  // ── File and parse errors ───────────────────────────────────────────────────

  /// Scenario: Input file does not exist
  #[test]
  fn test_missing_file() {
    let result = run(&arg("/tmp/nonexistent_electron_cli_test.yaml"));
    assert!(result.is_err());
  }

  /// Scenario: Input file contains invalid YAML
  #[test]
  fn test_invalid_yaml() {
    let f = temp_file("driver: md\nmolecule: :\n  bad:");
    let result = run(&arg(f.path().to_str().unwrap()));
    assert!(result.is_err());
  }

  /// Scenario: Input file is missing a required field
  #[test]
  fn test_missing_required_field() {
    let f = temp_file(
      "molecule:\n  symbols: [H]\n  geometry: [0.0, 0.0, 0.0]\nmodel:\n  method: hf\n  basis: sto-3g\n",
    );
    let result = run(&arg(f.path().to_str().unwrap()));
    assert!(result.is_err());
  }

  /// Scenario: Input file contains an invalid value
  #[test]
  fn test_invalid_driver_value() {
    let f = temp_file(
      "driver: optimize\nmolecule:\n  symbols: [H]\n  geometry: [0.0, 0.0, 0.0]\nmodel:\n  method: hf\n  basis: sto-3g\n",
    );
    let result = run(&arg(f.path().to_str().unwrap()));
    assert!(result.is_err());
  }
}
