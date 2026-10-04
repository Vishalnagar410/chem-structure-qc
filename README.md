# Chemical Structure Validation, Canonicalization & QC

## 1. Overview

`chem-structure-qc` is a small Python/RDKit package and command-line tool for assessing chemical structures in CSV datasets. It produces canonical SMILES, common molecular descriptors, duplicate-group information, and a machine-readable QC summary.

## 2. Why this matters in chemical data workflows

Structure errors and duplicate records can distort compound counts, descriptor analyses, and downstream model datasets. A repeatable QC pass makes those issues visible before analysis or model preparation; it does not make scientific decisions about which record should be retained.

## 3. Features

- Validate required input columns and parse SMILES with RDKit.
- Keep one result row per input and report invalid/missing SMILES with a reason.
- Generate RDKit canonical SMILES for valid molecules.
- Identify duplicate rows and groups using exact canonical-SMILES matches.
- Calculate molecular weight, Crippen logP, TPSA, HBD, HBA, rotatable bonds, heavy atoms, and ring count.
- Write CSV outputs and a JSON summary, including validity and duplicate rates.

Validation, canonicalization, duplicate detection, and descriptor calculation are separate steps in the workflow.

## 4. Workflow

```text
CSV
→ SMILES validation
→ canonicalization
→ duplicate detection
→ descriptor calculation
→ QC outputs
```

## 5. Example input

The small educational dataset in `data/example_compounds.csv` contains common example structures, one intentionally invalid SMILES, and equivalent ethanol notations (`CCO` and `OCC`) to demonstrate duplicate detection. It is not experimental or proprietary data.

```csv
compound_id,smiles
CMP001,CCO
CMP002,OCC
CMP003,c1ccccc1
```

The required columns are `compound_id` and `smiles`. A custom SMILES column can be selected with `--smiles-column`.

## 6. Example QC result

Running the included data reports:

```text
Input rows          : 6
Valid structures    : 5
Invalid structures  : 1
Unique structures   : 4
Duplicate rows      : 2
Duplicate groups    : 1
```

## 7. Installation

Python 3.10 or newer is required. Install the package and test dependency in a virtual environment:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
```

RDKit wheels are available for many supported Python/platform combinations. If pip cannot resolve an RDKit wheel for your platform, use a conda-forge environment with `rdkit` and `pandas`, then install this project with `python -m pip install -e ".[test]"`.

## Interactive Dashboard

The optional Streamlit dashboard is a presentation layer over the existing QC package; it does not duplicate or extend the scientific processing. Install the UI extra and launch it from the repository root:

```powershell
python -m pip install -e ".[ui]"
streamlit run app/app.py
```

The dashboard accepts the included example dataset or an uploaded CSV, lets you choose the SMILES column, explores valid and invalid structures, and provides downloads of the same CSV and JSON QC outputs as the core package. Uploaded data is processed in memory only; the app has no login, database, external API calls, or persistent user storage. Scientific limitations described below also apply to the dashboard.

## 8. Usage

```powershell
chem-structure-qc --input data/example_compounds.csv --output outputs
```

Use `chem-structure-qc --help` to see the available arguments. The output directory is created automatically. Missing input files, malformed CSVs, and missing required columns produce a concise error and nonzero exit status.

## 9. Output files

- `qc_results.csv` — every input row, validity, canonical SMILES, descriptors, and duplicate flag.
- `invalid_structures.csv` — invalid or missing structures with `error_reason`.
- `duplicate_groups.csv` — canonical SMILES, group size, and member compound IDs encoded as a JSON array.
- `qc_report.json` — row counts, unique structure count, duplicate group/row counts, and rates.

Representative outputs from the included dataset are checked in under `examples/example_output/`.

## 10. Scientific notes / limitations

Canonicalization is **not** equivalent to complete chemical standardization. Duplicate detection uses exact RDKit canonical-SMILES strings under the installed RDKit behavior. The package does not attempt production-grade chemical standardization, tautomer normalization, comprehensive salt stripping, charge normalization, isotope normalization, or a policy for stereochemistry. Descriptor failures are represented as missing values rather than aborting the dataset run.

## 11. Testing

```powershell
python -m pytest -q
```

The deterministic tests use only local example/synthetic data and do not require internet access.

## 12. Project structure

```text
chem-structure-qc/
├── .streamlit/config.toml
├── CHANGELOG.md
├── LICENSE
├── README.md
├── app/
│   ├── app.py
│   └── ui_helpers.py
├── data/example_compounds.csv
├── examples/example_output/
├── pyproject.toml
├── src/chem_structure_qc/
│   ├── __init__.py
│   ├── cli.py
│   └── qc.py
└── tests/
    ├── test_cli.py
    ├── test_qc.py
    └── test_ui_helpers.py
```

## 13. Future enhancements

Potential extensions include configurable salt/fragment handling, charge normalization, tautomer handling, stereochemistry policies, and richer validation rules. These require explicit scientific policies and dedicated tests before implementation.

## 14. License

This project is distributed under the MIT License. See [LICENSE](LICENSE).
