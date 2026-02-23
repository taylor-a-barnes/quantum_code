# Feature: Command Line Interface <!-- rq-bce67599 -->

The program (`electron`) is invoked from the command line with a single positional argument: the
path to a simulation input file. Relative paths are resolved against the working directory from
which the command was run. Absolute paths are used as-is.

On launch the program parses the input file (see `rqm/input/parser.md`). If parsing succeeds the
program writes a one-line summary to standard output and exits with code 0. If any error occurs
(missing argument, too many arguments, or any `InputError` from the parser) the program writes a
human-readable error message to standard error and exits with code 1. Nothing is written to
standard output on failure.

## Feature API <!-- rq-fbe1b430 -->

### Functions <!-- rq-6f4eafd6 -->

- `run(args: &[String]) -> Result<String, String>` <!-- rq-a0a7c2fa -->
  - Accepts the command-line arguments with argv[0] (the binary name) already stripped.
  - Returns `Err("usage: electron <input-file>".to_string())` if `args` does not contain
    exactly one element.
  - Otherwise resolves `args[0]` as a path (relative paths are resolved against the process
    working directory) and delegates to `parse_input`.
  - On success returns `Ok` containing the formatted summary line (without a trailing newline).
  - On failure returns `Err` containing the `Display` representation of the `InputError`
    (without the `"error: "` prefix).

- `main()` <!-- rq-72795519 -->
  - Collects `std::env::args()`, strips argv[0], and passes the remainder to `run`.
  - Prints the `Ok` string to standard output followed by a newline.
  - Prints `"error: {msg}"` to standard error and exits with code 1 on `Err`.

---

## Invocation <!-- rq-d6f2c78a -->

```
electron <input-file>
```

`<input-file>` is the only accepted argument. Exactly one positional argument is required; zero
arguments or more than one argument are both errors.

## Exit Codes <!-- rq-1f01b9f3 -->

| Code | Meaning                                                     |
| ---- | ----------------------------------------------------------- |
| 0    | Input file was parsed successfully.                         |
| 1    | Any error: wrong number of arguments, I/O failure, or parse error. |

## Standard Output (success) <!-- rq-e1c899d8 -->

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

## Standard Error (failure) <!-- rq-c56f7119 -->

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

## Gherkin Scenarios <!-- rq-b6b0de3d -->

```gherkin
Feature: Command line interface

  # --- Happy paths ---

  @rq-abf35b0a
  Scenario: Parse a valid Cartesian energy input file
    Given a valid simulation input file at "input.yaml" with driver "energy",
      2 atoms, method "hf", and basis "sto-3g"
    When the program is invoked with the argument "input.yaml"
    Then the program exits with code 0
    And standard output contains exactly "Parsed: driver=energy, method=hf, basis=sto-3g, atoms=2"
    And standard error is empty

  @rq-d53adf01
  Scenario: Parse a valid MD input file
    Given a valid simulation input file at "water_md.yaml" with driver "md",
      3 atoms, method "b3lyp", and basis "sto-3g"
    When the program is invoked with the argument "water_md.yaml"
    Then the program exits with code 0
    And standard output contains exactly "Parsed: driver=md, method=b3lyp, basis=sto-3g, atoms=3"

  @rq-8c9e0b32
  Scenario: Relative path is resolved from the working directory
    Given a valid simulation input file exists at "<cwd>/subdir/input.yaml"
    And the working directory is "<cwd>"
    When the program is invoked with the argument "subdir/input.yaml"
    Then the program exits with code 0

  @rq-e1160c83
  Scenario: Absolute path is accepted
    Given a valid simulation input file exists at "/tmp/input.yaml"
    When the program is invoked with the argument "/tmp/input.yaml"
    Then the program exits with code 0

  @rq-2509fb27
  Scenario: Atom count for Z-matrix geometry is the number of rows
    Given a valid simulation input file with a z_matrix containing 4 rows,
      method "hf", and basis "sto-3g"
    When the program is invoked with that file
    Then standard output contains "atoms=4"

  # --- Argument errors ---

  @rq-679028f8
  Scenario: No argument given
    When the program is invoked with no arguments
    Then the program exits with code 1
    And standard error contains "error: usage: electron <input-file>"
    And standard output is empty

  @rq-eb02af74
  Scenario: More than one argument given
    When the program is invoked with two arguments "a.yaml" and "b.yaml"
    Then the program exits with code 1
    And standard error contains "error: usage: electron <input-file>"
    And standard output is empty

  # --- File and parse errors ---

  @rq-4686ba54
  Scenario: Input file does not exist
    Given no file exists at "nonexistent.yaml"
    When the program is invoked with the argument "nonexistent.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  @rq-36dfa760
  Scenario: Input file contains invalid YAML
    Given a file at "bad.yaml" containing invalid YAML
    When the program is invoked with the argument "bad.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  @rq-2001bb71
  Scenario: Input file is missing a required field
    Given a file at "incomplete.yaml" that omits the required "driver" key
    When the program is invoked with the argument "incomplete.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty

  @rq-a7a57f4e
  Scenario: Input file contains an invalid value
    Given a file at "bad_driver.yaml" with driver set to "optimize"
    When the program is invoked with the argument "bad_driver.yaml"
    Then the program exits with code 1
    And standard error contains "error: "
    And standard output is empty
```
