# Feature: Command Line Interface

The program (`electron`) is invoked from the command line with a single positional argument: the
path to a simulation input file. Relative paths are resolved against the working directory from
which the command was run. Absolute paths are used as-is.

On launch the program parses the input file (see `rqm/input/parser.md`). If parsing succeeds the
program writes a one-line summary to standard output and exits with code 0. If any error occurs
(missing argument, too many arguments, or any `InputError` from the parser) the program writes a
human-readable error message to standard error and exits with code 1. Nothing is written to
standard output on failure.

## Invocation

```
electron <input-file>
```

`<input-file>` is the only accepted argument. Exactly one positional argument is required; zero
arguments or more than one argument are both errors.

## Exit Codes

| Code | Meaning                                                     |
| ---- | ----------------------------------------------------------- |
| 0    | Input file was parsed successfully.                         |
| 1    | Any error: wrong number of arguments, I/O failure, or parse error. |

## Standard Output (success)

When parsing succeeds, the program writes a single line to standard output of the form:

```
Parsed: driver=<driver>, method=<method>, basis=<basis>, atoms=<n>
```

Where:
- `<driver>` is one of `energy`, `gradient`, `hessian`, `md`.
- `<method>` and `<basis>` are the strings from the `model` block, exactly as stored (no
  normalisation is applied here).
- `<n>` is the total number of atoms in the molecule (i.e. the length of the `symbols` vector for
  Cartesian geometry, or the length of the `symbols` vector for Z-matrix geometry).

No other output is written to standard output on success.

## Standard Error (failure)

All error output is written to standard error. The format is:

```
error: <message>
```

Error conditions and their messages:

| Condition                        | Message                                                              |
| -------------------------------- | -------------------------------------------------------------------- |
| No argument given                | `usage: electron <input-file>`                                       |
| More than one argument given     | `usage: electron <input-file>`                                       |
| Any `InputError` from the parser | The `Display` representation of the `InputError` (see parser spec)   |

Nothing is written to standard output when an error occurs.

---

## Gherkin Scenarios

```gherkin
Feature: Command line interface

  # --- Happy paths ---

  Scenario: Parse a valid Cartesian energy input file
    Given a valid simulation input file at "input.yaml" with driver "energy",
      2 atoms, method "hf", and basis "sto-3g"
    When the program is invoked with the argument "input.yaml"
    Then the program exits with code 0
    And standard output contains exactly "Parsed: driver=energy, method=hf, basis=sto-3g, atoms=2"
    And standard error is empty

  Scenario: Parse a valid MD input file
    Given a valid simulation input file at "water_md.yaml" with driver "md",
      3 atoms, method "b3lyp", and basis "sto-3g"
    When the program is invoked with the argument "water_md.yaml"
    Then the program exits with code 0
    And standard output contains exactly "Parsed: driver=md, method=b3lyp, basis=sto-3g, atoms=3"

  Scenario: Relative path is resolved from the working directory
    Given a valid simulation input file exists at "<cwd>/subdir/input.yaml"
    And the working directory is "<cwd>"
    When the program is invoked with the argument "subdir/input.yaml"
    Then the program exits with code 0

  Scenario: Absolute path is accepted
    Given a valid simulation input file exists at "/tmp/input.yaml"
    When the program is invoked with the argument "/tmp/input.yaml"
    Then the program exits with code 0

  Scenario: Atom count for Z-matrix geometry is the number of rows
    Given a valid simulation input file with a z_matrix containing 4 rows,
      method "hf", and basis "sto-3g"
    When the program is invoked with that file
    Then standard output contains "atoms=4"

  # --- Argument errors ---

  Scenario: No argument given
    When the program is invoked with no arguments
    Then the program exits with code 1
    And standard error contains "error: usage: electron <input-file>"
    And standard output is empty

  Scenario: More than one argument given
    When the program is invoked with two arguments "a.yaml" and "b.yaml"
    Then the program exits with code 1
    And standard error contains "error: usage: electron <input-file>"
    And standard output is empty

  # --- File and parse errors ---

  Scenario: Input file does not exist
    Given no file exists at "nonexistent.yaml"
    When the program is invoked with the argument "nonexistent.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  Scenario: Input file contains invalid YAML
    Given a file at "bad.yaml" containing invalid YAML
    When the program is invoked with the argument "bad.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  Scenario: Input file is missing a required field
    Given a file at "incomplete.yaml" that omits the required "driver" key
    When the program is invoked with the argument "incomplete.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  Scenario: Input file contains an invalid value
    Given a file at "bad_driver.yaml" with driver set to "optimize"
    When the program is invoked with the argument "bad_driver.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty
```
